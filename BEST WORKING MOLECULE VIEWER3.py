import sys
import logging
import sqlite3
import random
import requests
import webbrowser
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFileDialog, QLineEdit, QMessageBox, QTableWidget, QTableWidgetItem, QHeaderView, QAction, QDialog, QCheckBox, QStyleFactory, QTreeWidget, QTreeWidgetItem, QTextEdit
from PyQt5.QtGui import QPixmap, QImage, QPalette, QColor
from PyQt5.QtCore import Qt, QTimer
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from PIL import Image
import io
import pubchempy as pcp

# Set up logging
logging.basicConfig(filename='molecule_viewer.log', level=logging.INFO, format='%(asctime)s - %(message)s')

# Settings Dialog for Dark Mode

class SettingsDialog(QDialog):
    def __init__(self, parent=None):
        super(SettingsDialog, self).__init__(parent)
        self.setWindowTitle("Settings")
        self.setGeometry(300, 300, 200, 200)
        self.settings_manager = SettingsManager()

        layout = QVBoxLayout()

        self.dark_mode_checkbox = QCheckBox("Enable Dark Mode", self)
        layout.addWidget(self.dark_mode_checkbox)

        # Load existing setting for dark mode
        self.dark_mode_checkbox.setChecked(parent.is_dark_mode)

        # Checkboxes for showing/hiding columns
        self.show_xlogp_checkbox = QCheckBox("Show XLogP", self)
        self.show_h_bond_donors_checkbox = QCheckBox("Show H-Bond Donors", self)
        self.show_h_bond_acceptors_checkbox = QCheckBox("Show H-Bond Acceptors", self)
        layout.addWidget(self.show_xlogp_checkbox)
        layout.addWidget(self.show_h_bond_donors_checkbox)
        layout.addWidget(self.show_h_bond_acceptors_checkbox)

        # Load existing settings for showing/hiding columns
        self.show_xlogp_checkbox.setChecked(self.settings_manager.get_setting('show_xlogp') == 'True')
        self.show_h_bond_donors_checkbox.setChecked(self.settings_manager.get_setting('show_h_bond_donors') == 'True')
        self.show_h_bond_acceptors_checkbox.setChecked(self.settings_manager.get_setting('show_h_bond_acceptors') == 'True')

        # Apply button
        apply_button = QPushButton("Apply", self)
        apply_button.clicked.connect(self.apply_settings)
        layout.addWidget(apply_button)

        self.setLayout(layout)

    def apply_settings(self):
        # Toggle dark mode based on checkbox
        self.parent().toggle_dark_mode(self.dark_mode_checkbox.isChecked())

        # Save column visibility settings
        self.settings_manager.save_setting('show_xlogp', str(self.show_xlogp_checkbox.isChecked()))
        self.settings_manager.save_setting('show_h_bond_donors', str(self.show_h_bond_donors_checkbox.isChecked()))
        self.settings_manager.save_setting('show_h_bond_acceptors', str(self.show_h_bond_acceptors_checkbox.isChecked()))

        self.accept()


# Molecule Analysis Class
class MoleculeAnalysis:
    def __init__(self, mol):
        self.mol = mol

    def get_atom_counts(self):
        atom_counts = {}
        for atom in self.mol.GetAtoms():
            element = atom.GetSymbol()
            atom_counts[element] = atom_counts.get(element, 0) + 1
        return atom_counts

    def calculate_bond_info(self):
        bond_info = []
        electronegativity = {
            'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
            'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'I': 2.66
        }
        for bond in self.mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            elem1 = atom1.GetSymbol()
            elem2 = atom2.GetSymbol()
            if elem1 in electronegativity and elem2 in electronegativity:
                en_diff = abs(electronegativity[elem1] - electronegativity[elem2])
                bond_length = Chem.rdMolTransforms.GetBondLength(self.mol.GetConformer(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                bond_info.append((f"{elem1}-{elem2}", f"{en_diff:.2f}", f"{bond_length:.2f} Å"))
        return bond_info

# Main Application Window
class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Molecule Viewer")
        self.setGeometry(100, 100, 1400, 800)  # Larger window size

        # Initialize the Settings Manager
        self.settings_manager = SettingsManager()

        # Load dark mode setting
        dark_mode_setting = self.settings_manager.get_setting('dark_mode')
        self.is_dark_mode = dark_mode_setting == 'True'
        self.toggle_dark_mode(self.is_dark_mode)  # Apply dark mode if enabled

        # Initialize the Database Manager
        self.db_manager = DatabaseManager()

        # Central Widget
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)

        # Layouts
        self.main_layout = QVBoxLayout(self.central_widget)
        self.top_layout = QHBoxLayout()
        self.main_layout.addLayout(self.top_layout)
        self.bottom_layout = QVBoxLayout()
        self.main_layout.addLayout(self.bottom_layout)

        # Input field
        self.input_field = QLineEdit(self)
        self.input_field.setPlaceholderText("Enter SMILES, InChI, or common name...")
        self.input_field.returnPressed.connect(self.draw_molecule)  # Trigger drawing on Enter key
        self.top_layout.addWidget(self.input_field)
        
##        # Caps Lock indicator check using a QTimer
##        self.caps_lock_timer = QTimer(self)
##        self.caps_lock_timer.timeout.connect(self.check_caps_lock_state)
##        self.caps_lock_timer.start(500)  # Check every 500ms

        # Draw button
        self.draw_button = QPushButton("Draw Molecule", self)
        self.draw_button.clicked.connect(self.draw_molecule)
        self.top_layout.addWidget(self.draw_button)

        # Save button
        self.save_button = QPushButton("Save Image", self)
        self.save_button.clicked.connect(self.save_image)
        self.top_layout.addWidget(self.save_button)

        # Update all database button
        self.update_button = QPushButton("Update All Database", self)
        self.update_button.clicked.connect(self.update_all_database)
        self.top_layout.addWidget(self.update_button)
        
        # Add the new button to open the PubChem page
        self.open_pubchem_button = QPushButton("Open PubChem Page", self)
        self.open_pubchem_button.clicked.connect(self.open_pubchem_page)
        self.top_layout.addWidget(self.open_pubchem_button)
        
        # Left layout for Atom Count Table
        self.left_layout = QVBoxLayout()
        self.top_layout.addLayout(self.left_layout)

        # Atom Count Table
        self.atom_count_tree = QTreeWidget(self)
        self.atom_count_tree.setColumnCount(2)
        self.atom_count_tree.setHeaderLabels(["Element", "Count"])
        self.left_layout.addWidget(self.atom_count_tree)

        # Center layout for Molecule Image
        self.center_layout = QVBoxLayout()
        self.top_layout.addLayout(self.center_layout)

        # Display area for molecule image
        self.image_label = QLabel(self)
        self.center_layout.addWidget(self.image_label)

        # Right layout for Bond Information Table
        self.right_layout = QVBoxLayout()
        self.top_layout.addLayout(self.right_layout)

        # Bond Information Table
        self.bond_info_tree = QTreeWidget(self)
        self.bond_info_tree.setColumnCount(4)
        self.bond_info_tree.setHeaderLabels(["Bond", "Electronegativity Difference", "Length", "Polarity"])
        self.bottom_layout.addWidget(self.bond_info_tree)

        # Database view (bottom)
        self.db_table = QTableWidget(self)
        self.db_table.setColumnCount(4)
        self.db_table.setHorizontalHeaderLabels(["SMILES", "Atomic Mass", "IUPAC Name", "Formula"])
        self.db_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.db_table.cellClicked.connect(self.on_db_select)  # Automatically draw when selecting from database
        self.bottom_layout.addWidget(self.db_table)
        self.populate_database_view()
        #self.db_table.setSortingEnabled(True)  # Enable sorting on columns
        
        # Molecule Info
        self.info_label = QLabel(self)
        self.top_layout.addWidget(self.info_label)

        # Current image (for saving)
        self.current_image = None

        # Dark Mode State
        self.is_dark_mode = False

        # Create a menu for settings
        self.create_menu()

        # Set modern theme
        QApplication.setStyle(QStyleFactory.create("Fusion"))

        # Install event filter for Caps Lock indicator
        # self.input_field.installEventFilter(self)

    def open_pubchem_page(self):
        """Open the PubChem page for the selected molecule using the SMILES query."""
        selected_smiles = self.input_field.text().strip()
        if selected_smiles:
            url = f"https://pubchem.ncbi.nlm.nih.gov/#query={selected_smiles}&input_type=smiles"
            webbrowser.open(url)
        else:
            QMessageBox.warning(self, "No SMILES Input", "Please enter a valid SMILES string or select a molecule from the database.")

    def create_menu(self):
        """Create a menu for accessing settings and viewing logs."""
        menubar = self.menuBar()

        # Settings action
        settings_action = QAction("Settings", self)
        settings_action.triggered.connect(self.open_settings_dialog)
        settings_menu = menubar.addMenu("Settings")
        settings_menu.addAction(settings_action)

        # Log viewer action
        log_viewer_action = QAction("View Logs", self)
        log_viewer_action.triggered.connect(self.open_log_viewer)
        menubar.addAction(log_viewer_action)
        
    def open_settings_dialog(self):
        """Open the settings dialog."""
        dialog = SettingsDialog(self)
        dialog.exec_()

    def open_log_viewer(self):
        """Open the log viewer."""
        log_viewer = QDialog(self)
        log_viewer.setWindowTitle("Log Viewer")
        log_viewer.setGeometry(100, 100, 600, 400)

        layout = QVBoxLayout()

        log_text_edit = QTextEdit(log_viewer)
        log_text_edit.setReadOnly(True)
        try:
            with open('molecule_viewer.log', 'r') as log_file:
                log_text_edit.setPlainText(log_file.read())
        except FileNotFoundError:
            log_text_edit.setPlainText("Log file not found.")

        layout.addWidget(log_text_edit)
        log_viewer.setLayout(layout)
        log_viewer.exec_()

            
    def toggle_dark_mode(self, enable):
        """Toggle dark mode on or off."""
        if enable:
            palette = QPalette()
            palette.setColor(QPalette.Window, QColor(53, 53, 53))
            palette.setColor(QPalette.WindowText, Qt.white)
            palette.setColor(QPalette.Base, QColor(25, 25, 25))
            palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
            palette.setColor(QPalette.ToolTipBase, Qt.white)
            palette.setColor(QPalette.ToolTipText, Qt.white)
            palette.setColor(QPalette.Text, Qt.white)
            palette.setColor(QPalette.Button, QColor(53, 53, 53))
            palette.setColor(QPalette.ButtonText, Qt.white)
            palette.setColor(QPalette.BrightText, Qt.red)
            palette.setColor(QPalette.Link, QColor(42, 130, 218))
            palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
            palette.setColor(QPalette.HighlightedText, Qt.black)
            QApplication.setPalette(palette)
        else:
            QApplication.setPalette(QApplication.style().standardPalette())

        self.is_dark_mode = enable
        self.settings_manager.save_setting('dark_mode', str(enable))  # Save dark mode state to settings

    def draw_molecule(self):
        user_input = self.input_field.text().strip()
        mol = None

        # Try to interpret input as SMILES
        try:
            mol = Chem.MolFromSmiles(user_input)
            if mol:
                Chem.SanitizeMol(mol)
                logging.info(f"SMILES input recognized: {user_input}")
        except Exception as e:
            logging.error(f"Error interpreting SMILES: {e}")

        # Process molecule
        if mol:
            try:
                AllChem.Compute2DCoords(mol)
                img = Draw.MolToImage(mol)
                img_bytes = io.BytesIO()
                img.save(img_bytes, format='PNG')  # Save as PNG
                img_bytes.seek(0)  # Move cursor to the beginning of the file-like object

                # Convert image to QImage and display
                qimage = QImage.fromData(img_bytes.getvalue())  # Convert BytesIO to QImage
                pixmap = QPixmap.fromImage(qimage)

                # Display the image on the QLabel
                self.image_label.setPixmap(pixmap)
                self.current_image = Image.open(io.BytesIO(img_bytes.getvalue()))  # Reopen BytesIO for saving

                # Fetch additional data from PubChem
                pubchem_data = self.fetch_pubchem_data(user_input)
                if pubchem_data:
                    # **Corrected variables**
                    molecular_formula = pubchem_data.molecular_formula or "N/A"
                    atomic_mass = pubchem_data.molecular_weight or "N/A"
                    iupac_name = pubchem_data.iupac_name or "N/A"

                    # **Ensure all arguments, including formula, are passed to add_molecule**
                    self.db_manager.add_molecule(user_input, atomic_mass, iupac_name, molecular_formula)
                    self.populate_database_view()

                    # **Display molecule information in the right layout**
                    self.info_label.setText(f"IUPAC: {iupac_name}\nFormula: {molecular_formula}\nWeight: {atomic_mass} g/mol")

                    # Populate atom and bond info
                    self.populate_atom_bond_info(mol)

                    # Detect functional groups using FunctionalGroupDetector
                    detector = FunctionalGroupDetector()  # Create an instance of the detector
                    detected_groups = detector.detect(mol)
                    if detected_groups:
                        functional_groups_text = "\n".join(detected_groups)
                        self.info_label.setText(self.info_label.text() + f"\nFunctional Groups:\n{functional_groups_text}")
                    else:
                        self.info_label.setText(self.info_label.text() + "\nFunctional Groups: None")

            except Exception as e:
                logging.error(f"Error processing molecule: {e}")
                QMessageBox.critical(self, "Error", f"Could not process the molecule: {e}")
        else:
            QMessageBox.warning(self, "Invalid Input", "Could not interpret input as a valid molecule.")

    def fetch_pubchem_data(self, smiles):
        """Fetch data from PubChem using pubchempy."""
        try:
            compound = pcp.get_compounds(smiles, namespace='smiles')
            if compound:
                return compound[0]
            else:
                logging.warning(f"No data found for SMILES {smiles}.")
                return None
        except Exception as e:
            logging.error(f"Error fetching data from PubChem: {e}")
            return None

    def populate_atom_bond_info(self, mol):
        """Populate atom counts and bond information in the tree widgets."""
        self.atom_count_tree.clear()  # Clear previous entries
        self.bond_info_tree.clear()  # Clear previous entries

        # Analyze the molecule
        analysis = MoleculeAnalysis(mol)
        atom_counts = analysis.get_atom_counts()
        bond_info = analysis.calculate_bond_info()

        # Color map for elements
        element_colors = {
            'H': 'white', 'C': 'gray', 'N': 'blue', 'O': 'red', 'F': 'green', 
            'Cl': 'green', 'Br': 'brown', 'I': 'purple', 'S': 'yellow', 'P': 'orange'
        }

        # Populate atom counts
        for atom, count in atom_counts.items():
            atom_item = QTreeWidgetItem([atom, str(count)])
            color = element_colors.get(atom, 'black')
            atom_item.setForeground(0, QColor(color))
            self.atom_count_tree.addTopLevelItem(atom_item)

        # Populate bond information
        for bond in bond_info:
            bond_item = QTreeWidgetItem(bond)
            self.bond_info_tree.addTopLevelItem(bond_item)

    def save_image(self):
        """Save the currently displayed molecule image."""
        if self.current_image:
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Image", "", "PNG files (*.png)")
            if file_path:
                self.current_image.save(file_path)
                QMessageBox.information(self, "Image Saved", "Molecule image has been saved.")
                logging.info(f"Image saved: {file_path}")
        else:
            QMessageBox.warning(self, "No Image", "No image available to save.")
        #caps
            ''''''
    def eventFilter(self, obj, event):
        """Handle Caps Lock indicator."""
        if obj == self.input_field and event.type() == QEvent.KeyPress:
            self.update_caps_lock_indicator()
        return super().eventFilter(obj, event)

    def update_caps_lock_indicator(self):
        """Update Caps Lock indicator based on keyboard state."""
        caps_on = QApplication.keyboardModifiers() & Qt.CapsLockModifier
        self.input_field.setPlaceholderText("Enter SMILES... [CAPS ON]" if caps_on else "Enter SMILES...")

    def update_all_database(self):
        """Update all entries in the database with missing information."""
        self.db_manager.update_all_molecules()
        self.populate_database_view()
        QMessageBox.information(self, "Update Complete", "All database entries have been updated.")

    def populate_database_view(self):
        # Load settings for which columns to show
        show_xlogp = self.settings_manager.get_setting('show_xlogp') == 'True'
        show_h_bond_donors = self.settings_manager.get_setting('show_h_bond_donors') == 'True'
        show_h_bond_acceptors = self.settings_manager.get_setting('show_h_bond_acceptors') == 'True'

        # Set headers and column count based on the desired order
        columns = ["SMILES", "IUPAC Name", "Formula", "Atomic Mass"]
        if show_xlogp:
            columns.append("XLogP")
        if show_h_bond_donors:
            columns.append("H-Bond Donors")
        if show_h_bond_acceptors:
            columns.append("H-Bond Acceptors")
        
        self.db_table.setColumnCount(len(columns))
        self.db_table.setHorizontalHeaderLabels(columns)

        # Fetch data from the database
        molecules = self.db_manager.get_all_molecules()
        self.db_table.setRowCount(len(molecules))

        for row_idx, molecule in enumerate(molecules):
            # Ensure data matches header order: SMILES, IUPAC Name, Formula, Atomic Mass
            self.db_table.setItem(row_idx, 0, QTableWidgetItem(str(molecule[0])))  # SMILES
            self.db_table.setItem(row_idx, 1, QTableWidgetItem(str(molecule[2])))  # IUPAC Name
            self.db_table.setItem(row_idx, 2, QTableWidgetItem(str(molecule[3])))  # Formula
            self.db_table.setItem(row_idx, 3, QTableWidgetItem(str(molecule[1])))  # Atomic Mass
            
            if show_xlogp and len(molecule) > 4:
                self.db_table.setItem(row_idx, 4, QTableWidgetItem(str(molecule[4])))  # XLogP
            
            if show_h_bond_donors and len(molecule) > 5:
                self.db_table.setItem(row_idx, 5, QTableWidgetItem(str(molecule[5])))  # H-Bond Donors
            
            if show_h_bond_acceptors and len(molecule) > 6:
                self.db_table.setItem(row_idx, 6, QTableWidgetItem(str(molecule[6])))  # H-Bond Acceptors

    def on_db_select(self, row, column):
        """Handle selection in the database table to automatically draw the molecule and update missing data."""
        selected_smiles = self.db_table.item(row, 0).text()  # Get the SMILES string from the first column
        
        # Fetch the corresponding data from the database
        self.db_manager.cursor.execute("SELECT atomic_mass, iupac_name, molecular_formula FROM molecules WHERE smiles = ?", (selected_smiles,))
        result = self.db_manager.cursor.fetchone()

        if result:
            # **Corrected variables**
            atomic_mass, iupac_name, molecular_formula = result

            # **Display the fetched information in the labels**
            self.info_label.setText(f"IUPAC: {iupac_name or 'N/A'}\nFormula: {molecular_formula or 'N/A'}\nWeight: {atomic_mass or 'N/A'} g/mol")
            
        self.input_field.setText(selected_smiles)  # Set the input field with the selected SMILES
        self.draw_molecule()  # Draw the molecule based on the selected SMILES
            
    ##    def check_caps_lock_state(self):
    ##        """Check the state of Caps Lock and update the input field's placeholder text."""
    ##        caps_on = QApplication.keyboardModifiers() & Qt.CapsLockModifier
    ##        placeholder_text = "Enter SMILES... [CAPS ON]" if caps_on else "Enter SMILES, InChI, or common name..."
    ##        self.input_field.setPlaceholderText(placeholder_text)
##        
class FunctionalGroupDetector:
    def __init__(self):
        # Initialize the SMARTS patterns for common functional groups
        self.functional_groups = {
            'Alkane (single bond)': '[CX4]',
            'Cycloalkane': '[R][CX4][R]',  # Detects rings with single-bonded carbon atoms
            'Alkene (double bond)': '[CX3]=[CX3]',
            'Alkyne (triple bond)': '[CX2]#C',
            'Aromatic (benzene ring)': 'c1ccccc1',  # Aromatic benzene ring
            'Non-aromatic six-membered ring': 'C1=CC=CC=C1',  # Non-aromatic version
            'Alcohol (hydroxyl group)': '[OX2H]',
            'Ether (R-O-R)': '[OD2]([#6])[#6]',
            'Aldehyde (formyl group)': '[CX3H1](=O)[#6]',
            'Ketone (carbonyl group)': '[CX3](=O)[#6]',
            'Carboxylic Acid (carboxyl group)': '[CX3](=O)[OX2H1]',
            'Ester (COOR group)': '[CX3](=O)[OX2H0][#6]',
            'Amine (primary)': '[NX3;H2,H1][#6]',
            'Amine (secondary)': '[NX3;H1][#6][#6]',
            'Amine (tertiary)': '[NX3]([#6])[#6][#6]',
            'Amide (CONH2)': '[NX3][CX3](=[OX1])[#6]',
            'Nitro (NO2 group)': '[NX3](=O)[O-]',
            'Thiol (mercaptan)': '[SX2H]',
            'Halide (F, Cl, Br, I)': '[F,Cl,Br,I]',
            'Nitrile (C≡N)': '[CX2]#N',
            'Imine (C=N)': '[CX2]=[NX3]',
            'Carboxylate (COO-)': '[CX3](=O)[O-]',
            'Phenol (aromatic OH)': 'c1ccc(O)cc1',
            'Sulfonic Acid (SO3H)': 'S(=O)(=O)[OH]'
        }
    
    def detect(self, mol):
        """
        Detect functional groups in a given molecule.
        
        Parameters:
        mol (rdkit.Chem.Mol): The molecule to analyze.
        
        Returns:
        list: A list of detected functional group names.
        """
        detected_groups = []
        for group_name, smarts in self.functional_groups.items():
            pattern = Chem.MolFromSmarts(smarts)
            if mol.HasSubstructMatch(pattern):
                detected_groups.append(group_name)
        return detected_groups

class SettingsManager:
    def __init__(self, db_name='settings.db'):
        """Initialize the settings database."""
        self.db_name = db_name
        self.connection = sqlite3.connect(self.db_name)
        self.cursor = self.connection.cursor()
        self.create_table()

    def create_table(self):
        """Create settings table if it doesn't exist."""
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS settings (
                id INTEGER PRIMARY KEY,
                setting_name TEXT UNIQUE,
                setting_value TEXT
            )
        ''')
        self.connection.commit()

    def save_setting(self, setting_name, setting_value):
        """Save or update a setting in the database."""
        self.cursor.execute('''
            INSERT OR REPLACE INTO settings (setting_name, setting_value)
            VALUES (?, ?)
        ''', (setting_name, setting_value))
        self.connection.commit()

    def get_setting(self, setting_name):
        """Retrieve a setting from the database."""
        self.cursor.execute('SELECT setting_value FROM settings WHERE setting_name = ?', (setting_name,))
        result = self.cursor.fetchone()
        return result[0] if result else None

    def close(self):
        """Close the database connection."""
        self.connection.close()


# Database Manager Class
class DatabaseManager:
    def __init__(self, db_name='molecules.db'):
        print("Initializing database manager")
        self.db_name = db_name
        self.connection = sqlite3.connect(self.db_name)
        self.cursor = self.connection.cursor()
        self.create_table()
        self.update_existing_table()  # Ensure the schema is updated
        
    def create_table(self):
        print("Creating molecules table")
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS molecules (
                id INTEGER PRIMARY KEY,
                smiles TEXT UNIQUE,
                iupac_name TEXT,
                molecular_formula TEXT,
                atomic_mass REAL
            )
        ''')
        self.connection.commit()

    def add_molecule(self, smiles, atomic_mass, iupac_name, molecular_formula):
        print(f"Adding molecule: {smiles}")
        # Convert "N/A" to None for consistency
        iupac_name = None if iupac_name == "N/A" else iupac_name
        molecular_formula = None if molecular_formula == "N/A" else molecular_formula
        
        try:
            self.cursor.execute(
                "INSERT INTO molecules (smiles, iupac_name, molecular_formula, atomic_mass) VALUES (?, ?, ?, ?)",
                (smiles, iupac_name, molecular_formula, atomic_mass)
            )
            self.connection.commit()
        except sqlite3.IntegrityError:
            logging.warning(f"Molecule with SMILES {smiles} already exists in the database.")

            
    def update_existing_table(self):
 
        # Check and add missing columns if they don't exist
 
        self.cursor.execute("PRAGMA table_info(molecules)")
 
        columns = [column[1] for column in self.cursor.fetchall()]
 
 
 
        if 'molecular_formula' not in columns:
 
            self.cursor.execute("ALTER TABLE molecules ADD COLUMN molecular_formula TEXT")
 
            self.connection.commit()
 
            logging.info("Added 'molecular_formula' column to 'molecules' table.")
 
        else:
 
            logging.info("'molecular_formula' column already exists.")
    def update_molecule(self, smiles, atomic_mass, iupac_name, molecular_formula):
        print(f"Updating molecule: {smiles}")
        iupac_name = None if iupac_name == "N/A" else iupac_name
        molecular_formula = None if molecular_formula == "N/A" else molecular_formula
        
        self.cursor.execute(
            "UPDATE molecules SET atomic_mass = ?, iupac_name = ?, molecular_formula = ? WHERE smiles = ?",
            (atomic_mass, iupac_name, molecular_formula, smiles)
        )
        self.connection.commit()

    def get_all_molecules(self):
        print("Fetching all molecules from the database")
        self.cursor.execute("SELECT smiles, atomic_mass, iupac_name, molecular_formula FROM molecules")
        return self.cursor.fetchall()

    def update_all_molecules(self):
        molecules = self.get_all_molecules()
        for molecule in molecules:
            smiles, atomic_mass, iupac_name, molecular_formula = molecule
            updated = False

            # Debugging: Log the molecule being processed
            logging.info(f"Processing molecule with SMILES: {smiles}")

            if iupac_name is None or iupac_name == "N/A" or molecular_formula is None or molecular_formula == "N/A" or atomic_mass is None:
                pubchem_data = self.fetch_pubchem_data(smiles)
                if pubchem_data:
                    # Log fetched data
                    logging.info(f"Fetched data for {smiles}: IUPAC Name: {pubchem_data.iupac_name}, "
                                 f"Molecular Formula: {pubchem_data.molecular_formula}, "
                                 f"Molecular Weight: {pubchem_data.molecular_weight}")

                    if (iupac_name is None or iupac_name == "N/A") and pubchem_data.iupac_name:
                        iupac_name = pubchem_data.iupac_name
                        updated = True
                    if (molecular_formula is None or molecular_formula == "N/A") and pubchem_data.molecular_formula:
                        molecular_formula = pubchem_data.molecular_formula
                        updated = True
                    if atomic_mass is None and pubchem_data.molecular_weight:
                        atomic_mass = pubchem_data.molecular_weight
                        updated = True

                if updated:
                    self.update_molecule(smiles, atomic_mass, iupac_name, molecular_formula)
                    logging.info(f"Updated molecule in database with SMILES: {smiles}")
                else:
                    logging.info(f"No updates made for molecule with SMILES: {smiles}")

            else:
                logging.info(f"No missing data for molecule with SMILES: {smiles}")



    def fetch_pubchem_data(self, smiles):
        try:
            compound = pcp.get_compounds(smiles, namespace='smiles')
            if compound:
                return compound[0]
            else:
                logging.warning(f"No data found for SMILES {smiles}.")
                return None
        except Exception as e:
            logging.error(f"Error fetching data from PubChem: {e}")
            return None

    def close(self):
        self.connection.close()
        
# Main execution
if __name__ == '__main__':
    app = QApplication(sys.argv)
    viewer = MoleculeViewer()
    viewer.show()
    sys.exit(app.exec_())

