import os
import logging
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QLineEdit, QTableWidget, QTableWidgetItem, QHeaderView, QAction, QDialog, QTreeWidget, QTreeWidgetItem, QMessageBox, QFileDialog, QStyleFactory, QTextEdit, QApplication
from PyQt5.QtGui import QPixmap, QImage, QPalette, QColor
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from PIL import Image
import io
import pubchempy as pcp
from settings_manager import SettingsManager
from database_manager import DatabaseManager
from functional_group_detector import FunctionalGroupDetector

class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecule Viewer")
        self.setGeometry(100, 100, 1400, 800)

        # Initialize managers
        self.settings_manager = SettingsManager()
        self.db_manager = DatabaseManager()

        # Load dark mode setting
        dark_mode_setting = self.settings_manager.get_setting('dark_mode')
        self.is_dark_mode = dark_mode_setting == 'True'
        self.toggle_dark_mode(self.is_dark_mode)

        # Central widget
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
        self.input_field.returnPressed.connect(self.draw_molecule)
        self.top_layout.addWidget(self.input_field)

        # Buttons
        self.draw_button = QPushButton("Draw Molecule", self)
        self.draw_button.clicked.connect(self.draw_molecule)
        self.top_layout.addWidget(self.draw_button)

        self.save_button = QPushButton("Save Image", self)
        self.save_button.clicked.connect(self.save_image)
        self.top_layout.addWidget(self.save_button)

        self.update_button = QPushButton("Update All Database", self)
        self.update_button.clicked.connect(self.update_all_database)
        self.top_layout.addWidget(self.update_button)

        self.open_pubchem_button = QPushButton("Open PubChem Page", self)
        self.open_pubchem_button.clicked.connect(self.open_pubchem_page)
        self.top_layout.addWidget(self.open_pubchem_button)

        # Tree widgets for atom counts and bond information
        self.atom_count_tree = QTreeWidget(self)
        self.atom_count_tree.setColumnCount(2)
        self.atom_count_tree.setHeaderLabels(["Element", "Count"])
        self.top_layout.addWidget(self.atom_count_tree)

        self.bond_info_tree = QTreeWidget(self)
        self.bond_info_tree.setColumnCount(3)
        self.bond_info_tree.setHeaderLabels(["Bond", "Electronegativity Difference", "Length"])
        self.bottom_layout.addWidget(self.bond_info_tree)

        # Database view
        self.db_table = QTableWidget(self)
        self.db_table.setColumnCount(4)
        self.db_table.setHorizontalHeaderLabels(["SMILES", "IUPAC Name", "Formula", "Atomic Mass"])
        self.db_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.db_table.cellClicked.connect(self.on_db_select)
        self.bottom_layout.addWidget(self.db_table)
        self.populate_database_view()

        # Molecule info label
        self.info_label = QLabel(self)
        self.top_layout.addWidget(self.info_label)

        # Current image
        self.current_image = None

        # Create menu
        self.create_menu()

        # Set modern theme
        QApplication.setStyle(QStyleFactory.create("Fusion"))

    def create_menu(self):
        menubar = self.menuBar()

        settings_action = QAction("Settings", self)
        settings_action.triggered.connect(self.open_settings_dialog)
        menubar.addAction(settings_action)

        log_viewer_action = QAction("View Logs", self)
        log_viewer_action.triggered.connect(self.open_log_viewer)
        menubar.addAction(log_viewer_action)

    def toggle_dark_mode(self, enable):
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
        self.settings_manager.save_setting('dark_mode', str(enable))


    def draw_molecule(self):
        user_input = self.input_field.text().strip()
        mol = Chem.MolFromSmiles(user_input)
        if mol:
            try:
                AllChem.Compute2DCoords(mol)
                img = Draw.MolToImage(mol)
                img_bytes = io.BytesIO()
                img.save(img_bytes, format='PNG')
                img_bytes.seek(0)

                qimage = QImage.fromData(img_bytes.getvalue())
                pixmap = QPixmap.fromImage(qimage)
                self.image_label.setPixmap(pixmap)
                self.current_image = Image.open(io.BytesIO(img_bytes.getvalue()))

                pubchem_data = self.fetch_pubchem_data(user_input)
                if pubchem_data:
                    molecular_formula = pubchem_data.molecular_formula or "N/A"
                    atomic_mass = pubchem_data.molecular_weight or "N/A"
                    iupac_name = pubchem_data.iupac_name or "N/A"
                    self.db_manager.add_molecule(user_input, atomic_mass, iupac_name, molecular_formula)
                    self.populate_database_view()
                    self.info_label.setText(f"IUPAC: {iupac_name}\nFormula: {molecular_formula}\nWeight: {atomic_mass} g/mol")

                    self.populate_atom_bond_info(mol)
                    detector = FunctionalGroupDetector()
                    detected_groups = detector.detect(mol)
                    functional_groups_text = "\n".join(detected_groups)
                    self.info_label.setText(self.info_label.text() + f"\nFunctional Groups:\n{functional_groups_text}")
            except Exception as e:
                logging.error(f"Error processing molecule: {e}")
                QMessageBox.critical(self, "Error", f"Could not process the molecule: {e}")
        else:
            QMessageBox.warning(self, "Invalid Input", "Could not interpret input as a valid molecule.")

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

    def populate_atom_bond_info(self, mol):
        self.atom_count_tree.clear()
        self.bond_info_tree.clear()

        atom_counts = {}
        for atom in mol.GetAtoms():
            element = atom.GetSymbol()
            atom_counts[element] = atom_counts.get(element, 0) + 1

        for atom, count in atom_counts.items():
            atom_item = QTreeWidgetItem([atom, str(count)])
            self.atom_count_tree.addTopLevelItem(atom_item)

        bond_info = []
        electronegativity = {
            'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
            'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'I': 2.66
        }
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            elem1 = atom1.GetSymbol()
            elem2 = atom2.GetSymbol()
            en_diff = abs(electronegativity[elem1] - electronegativity[elem2])
            bond_length = Chem.rdMolTransforms.GetBondLength(mol.GetConformer(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            bond_info.append((f"{elem1}-{elem2}", f"{en_diff:.2f}", f"{bond_length:.2f} Å"))

        for bond in bond_info:
            bond_item = QTreeWidgetItem(bond)
            self.bond_info_tree.addTopLevelItem(bond_item)

    def save_image(self):
        if self.current_image:
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Image", "", "PNG files (*.png)")
            if file_path:
                self.current_image.save(file_path)
                QMessageBox.information(self, "Image Saved", "Molecule image has been saved.")
                logging.info(f"Image saved: {file_path}")
        else:
            QMessageBox.warning(self, "No Image", "No image available to save.")

    def update_all_database(self):
        self.db_manager.update_all_molecules()
        self.populate_database_view()
        QMessageBox.information(self, "Update Complete", "All database entries have been updated.")

    def populate_database_view(self):
        self.db_table.setRowCount(0)
        molecules = self.db_manager.get_all_molecules()
        for row_idx, molecule in enumerate(molecules):
            self.db_table.insertRow(row_idx)
            self.db_table.setItem(row_idx, 0, QTableWidgetItem(molecule[0]))
            self.db_table.setItem(row_idx, 1, QTableWidgetItem(molecule[2]))
            self.db_table.setItem(row_idx, 2, QTableWidgetItem(molecule[3]))
            self.db_table.setItem(row_idx, 3, QTableWidgetItem(molecule[1]))

    def on_db_select(self, row, column):
        selected_smiles = self.db_table.item(row, 0).text()
        self.input_field.setText(selected_smiles)
        self.draw_molecule()

    def open_settings_dialog(self):
        dialog = SettingsDialog(self)
        dialog.exec_()

    def open_log_viewer(self):
        log_viewer = QDialog(self)
        log_viewer.setWindowTitle("Log Viewer")
        log_viewer.setGeometry(100, 100, 600, 400)
        layout = QVBoxLayout()

        log_text_edit = QTextEdit(log_viewer)
        log_text_edit.setReadOnly(True)
        try:
            with open(LOG_FILE_PATH, 'r') as log_file:
                log_text_edit.setPlainText(log_file.read())
        except FileNotFoundError:
            log_text_edit.setPlainText("Log file not found.")

        layout.addWidget(log_text_edit)
        log_viewer.setLayout(layout)
        log_viewer.exec_()

    def open_pubchem_page(self):
        selected_smiles = self.input_field.text().strip()
        if selected_smiles:
            url = f"https://pubchem.ncbi.nlm.nih.gov/#query={selected_smiles}&input_type=smiles"
            webbrowser.open(url)
        else:
            QMessageBox.warning(self, "No SMILES Input", "Please enter a valid SMILES string or select a molecule from the database.")
