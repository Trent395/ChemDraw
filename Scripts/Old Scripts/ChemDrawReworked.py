import sys
import logging
import sqlite3
import requests
import webbrowser
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFileDialog, QLineEdit, QMessageBox, QTableWidget, QTableWidgetItem, QHeaderView, QCheckBox, QStyleFactory, QTextEdit
from PyQt5.QtGui import QPixmap, QImage, QPalette, QColor
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem
import pubchempy as pcp
import io
from PIL import Image

# Set up logging
logging.basicConfig(filename='molecule_viewer.log', level=logging.INFO, format='%(asctime)s - %(message)s')


class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecule Viewer")
        self.setGeometry(100, 100, 1400, 800)

        self.connection = sqlite3.connect('molecules.db')
        self.cursor = self.connection.cursor()
        self.create_table()

        self.is_dark_mode = False
        self.initialize_ui()

    def initialize_ui(self):
        # Central Widget
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)

        # Layouts
        main_layout = QVBoxLayout(self.central_widget)
        top_layout = QHBoxLayout()
        main_layout.addLayout(top_layout)
        bottom_layout = QVBoxLayout()
        main_layout.addLayout(bottom_layout)

        # Input field
        self.input_field = QLineEdit(self)
        self.input_field.setPlaceholderText("Enter SMILES, InChI, or common name...")
        self.input_field.returnPressed.connect(self.draw_molecule)
        top_layout.addWidget(self.input_field)

        # Draw button
        draw_button = QPushButton("Draw Molecule", self)
        draw_button.clicked.connect(self.draw_molecule)
        top_layout.addWidget(draw_button)

        # Save button
        save_button = QPushButton("Save Image", self)
        save_button.clicked.connect(self.save_image)
        top_layout.addWidget(save_button)

        # Update all database button
        update_button = QPushButton("Update All Database", self)
        update_button.clicked.connect(self.update_all_database)
        top_layout.addWidget(update_button)

        # Open PubChem page button
        open_pubchem_button = QPushButton("Open PubChem Page", self)
        open_pubchem_button.clicked.connect(self.open_pubchem_page)
        top_layout.addWidget(open_pubchem_button)

        # ** Add this QLabel for displaying the molecule image **
        self.image_label = QLabel(self)
        self.image_label.setAlignment(Qt.AlignCenter)  # Center the image
        main_layout.addWidget(self.image_label)

        # Atom Count Table
        self.atom_count_table = QTableWidget(self)
        self.atom_count_table.setColumnCount(2)
        self.atom_count_table.setHorizontalHeaderLabels(["Element", "Count"])
        bottom_layout.addWidget(self.atom_count_table)

        # Bond Information Table
        self.bond_info_table = QTableWidget(self)
        self.bond_info_table.setColumnCount(4)
        self.bond_info_table.setHorizontalHeaderLabels(["Bond", "Electronegativity Difference", "Length", "Polarity"])
        bottom_layout.addWidget(self.bond_info_table)

        # Database view
        self.db_table = QTableWidget(self)
        self.db_table.setColumnCount(4)
        self.db_table.setHorizontalHeaderLabels(["SMILES", "IUPAC Name", "Formula", "Atomic Mass"])
        self.db_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.db_table.cellClicked.connect(self.on_db_select)
        bottom_layout.addWidget(self.db_table)
        self.populate_database_view()

        # Molecule Info
        self.info_label = QLabel(self)
        top_layout.addWidget(self.info_label)

        # Current image (for saving)
        self.current_image = None

        # Menu
        self.create_menu()

        # Set modern theme
        QApplication.setStyle(QStyleFactory.create("Fusion"))
    def create_menu(self):
        menubar = self.menuBar()
        settings_action = menubar.addAction("Toggle Dark Mode")
        settings_action.triggered.connect(self.toggle_dark_mode)

        log_viewer_action = menubar.addAction("View Logs")
        log_viewer_action.triggered.connect(self.open_log_viewer)

    def create_table(self):
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

    def fetch_pubchem_data(self, smiles):
        data = {'molecular_formula': 'N/A', 'atomic_mass': 'N/A', 'iupac_name': 'N/A', 'xlogp': 'N/A', 'h_bond_donors': 'N/A', 'h_bond_acceptors': 'N/A'}
        try:
            compound = pcp.get_compounds(smiles, namespace='smiles')
            if compound:
                compound = compound[0]
                data['molecular_formula'] = compound.molecular_formula or 'N/A'
                data['atomic_mass'] = compound.molecular_weight or 'N/A'
                data['iupac_name'] = compound.iupac_name or 'N/A'
                data['xlogp'] = compound.xlogp or 'N/A'

            mol = Chem.MolFromSmiles(smiles)
            if mol:
                data['h_bond_donors'] = Descriptors.NumHDonors(mol)
                data['h_bond_acceptors'] = Descriptors.NumHAcceptors(mol)

        except Exception as e:
            logging.error(f"Error fetching data for SMILES {smiles}: {e}")
        return data

    def upsert_molecule_in_database(self, smiles, data):
        try:
            self.cursor.execute("SELECT atomic_mass, iupac_name, molecular_formula FROM molecules WHERE smiles = ?", (smiles,))
            existing_data = self.cursor.fetchone()

            if existing_data:
                atomic_mass, iupac_name, molecular_formula = existing_data
                atomic_mass = data['atomic_mass'] if data['atomic_mass'] != 'N/A' and data['atomic_mass'] is not None else atomic_mass
                iupac_name = data['iupac_name'] if data['iupac_name'] != 'N/A' and data['iupac_name'] is not None else iupac_name
                molecular_formula = data['molecular_formula'] if data['molecular_formula'] != 'N/A' and data['molecular_formula'] is not None else molecular_formula

                self.cursor.execute(
                    "UPDATE molecules SET atomic_mass = ?, iupac_name = ?, molecular_formula = ? WHERE smiles = ?",
                    (atomic_mass, iupac_name, molecular_formula, smiles)
                )
            else:
                self.cursor.execute(
                    "INSERT INTO molecules (smiles, atomic_mass, iupac_name, molecular_formula) VALUES (?, ?, ?, ?)",
                    (smiles, data['atomic_mass'], data['iupac_name'], data['molecular_formula'])
                )
            self.connection.commit()

        except sqlite3.Error as e:
            logging.error(f"Database error: {e}")

    def draw_molecule(self):
        user_input = self.input_field.text().strip()
        mol = None

        try:
            mol = Chem.MolFromSmiles(user_input)
            if mol:
                Chem.SanitizeMol(mol)
                logging.info(f"SMILES input recognized: {user_input}")

            AllChem.Compute2DCoords(mol)
            img = Draw.MolToImage(mol)
            img_bytes = io.BytesIO()
            img.save(img_bytes, format='PNG')
            img_bytes.seek(0)

            qimage = QImage.fromData(img_bytes.getvalue())
            pixmap = QPixmap.fromImage(qimage)
            self.image_label.setPixmap(pixmap)
            self.current_image = Image.open(io.BytesIO(img_bytes.getvalue()))

            molecule_data = self.fetch_pubchem_data(user_input)
            self.upsert_molecule_in_database(user_input, molecule_data)

            info_text = f"IUPAC: {molecule_data['iupac_name']}\nFormula: {molecule_data['molecular_formula']}\nWeight: {molecule_data['atomic_mass']} g/mol"
            info_text += f"\nXLogP: {molecule_data['xlogp']}\nH-Bond Donors: {molecule_data['h_bond_donors']}\nH-Bond Acceptors: {molecule_data['h_bond_acceptors']}"
            self.info_label.setText(info_text)

            self.populate_atom_bond_info(mol)

        except Exception as e:
            logging.error(f"Error processing molecule: {e}")
            QMessageBox.critical(self, "Error", f"Could not process the molecule: {e}")

    def populate_atom_bond_info(self, mol):
        atom_counts = {}
        for atom in mol.GetAtoms():
            element = atom.GetSymbol()
            atom_counts[element] = atom_counts.get(element, 0) + 1

        bond_info = []
        electronegativity = {'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44}
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            elem1 = atom1.GetSymbol()
            elem2 = atom2.GetSymbol()
            en_diff = abs(electronegativity.get(elem1, 0) - electronegativity.get(elem2, 0))
            bond_length = Chem.rdMolTransforms.GetBondLength(mol.GetConformer(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            bond_info.append((f"{elem1}-{elem2}", f"{en_diff:.2f}", f"{bond_length:.2f} Å"))

        self.atom_count_table.setRowCount(len(atom_counts))
        for row_idx, (element, count) in enumerate(atom_counts.items()):
            self.atom_count_table.setItem(row_idx, 0, QTableWidgetItem(element))
            self.atom_count_table.setItem(row_idx, 1, QTableWidgetItem(str(count)))

        self.bond_info_table.setRowCount(len(bond_info))
        for row_idx, bond in enumerate(bond_info):
            for col_idx, value in enumerate(bond):
                self.bond_info_table.setItem(row_idx, col_idx, QTableWidgetItem(value))

    def save_image(self):
        if self.current_image:
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Image", "", "PNG files (*.png)")
            if file_path:
                self.current_image.save(file_path)
                QMessageBox.information(self, "Image Saved", "Molecule image has been saved.")
                logging.info(f"Image saved: {file_path}")
        else:
            QMessageBox.warning(self, "No Image", "No image available to save.")

    def open_pubchem_page(self):
        selected_smiles = self.input_field.text().strip()
        if selected_smiles:
            url = f"https://pubchem.ncbi.nlm.nih.gov/#query={selected_smiles}&input_type=smiles"
            webbrowser.open(url)
        else:
            QMessageBox.warning(self, "No SMILES Input", "Please enter a valid SMILES string or select a molecule from the database.")

    def toggle_dark_mode(self):
        self.is_dark_mode = not self.is_dark_mode
        palette = self.get_dark_mode_palette() if self.is_dark_mode else QApplication.style().standardPalette()
        QApplication.setPalette(palette)

    def get_dark_mode_palette(self):
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
        return palette

    def open_log_viewer(self):
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

    def populate_database_view(self):
        self.db_table.setRowCount(0)
        self.cursor.execute("SELECT smiles, iupac_name, molecular_formula, atomic_mass FROM molecules")
        molecules = self.cursor.fetchall()

        for row_idx, molecule in enumerate(molecules):
            self.db_table.insertRow(row_idx)
            for col_idx, value in enumerate(molecule):
                self.db_table.setItem(row_idx, col_idx, QTableWidgetItem(str(value)))

    def on_db_select(self, row, column):
        selected_smiles = self.db_table.item(row, 0).text()
        self.input_field.setText(selected_smiles)
        self.draw_molecule()

    def update_all_database(self):
        self.cursor.execute("SELECT smiles FROM molecules")
        smiles_list = [row[0] for row in self.cursor.fetchall()]

        for smiles in smiles_list:
            molecule_data = self.fetch_pubchem_data(smiles)
            self.upsert_molecule_in_database(smiles, molecule_data)

        self.populate_database_view()
        QMessageBox.information(self, "Update Complete", "All database entries have been updated.")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    viewer = MoleculeViewer()
    viewer.show()
    sys.exit(app.exec_())
