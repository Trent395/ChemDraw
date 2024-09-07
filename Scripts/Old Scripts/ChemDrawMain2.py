import sys
import logging
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QHBoxLayout, QLabel, QPushButton, QFileDialog, QLineEdit, QMessageBox, QTreeWidget, QTreeWidgetItem
from PyQt5.QtGui import QPixmap, QImage, QColor
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from PIL import Image
import io

# Import from other scripts
from molecule_db import MoleculeDatabase
from molecule_functions import generate_random_smiles, get_molecule_info, get_longest_carbon_chain, fetch_pubchem_data
from database_manager import DatabaseManager
from molecule_viewer_gui import MoleculeViewer

# Set up logging
logging.basicConfig(filename='molecule_viewer.log', level=logging.INFO, format='%(asctime)s - %(message)s')

class ChemDrawMain(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Molecule Viewer")
        self.setGeometry(100, 100, 1400, 800)

        # Database Manager
        self.db_manager = MoleculeDatabase()

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
        self.input_field.setPlaceholderText("Enter SMILES...")
        self.input_field.returnPressed.connect(self.draw_molecule)
        self.top_layout.addWidget(self.input_field)

        # Draw button
        self.draw_button = QPushButton("Draw Molecule", self)
        self.draw_button.clicked.connect(self.draw_molecule)
        self.top_layout.addWidget(self.draw_button)

        # Image label for displaying the molecule
        self.image_label = QLabel(self)
        self.top_layout.addWidget(self.image_label)  # Add the image label to the layout


        # Save button
        self.save_button = QPushButton("Save Image", self)
        self.save_button.clicked.connect(self.save_image)
        self.top_layout.addWidget(self.save_button)

        # Database view (bottom)
        self.db_table = QTreeWidget(self)
        self.db_table.setColumnCount(4)
        self.db_table.setHeaderLabels(["SMILES", "Atomic Mass", "IUPAC Name", "Formula"])
        self.db_table.itemClicked.connect(self.on_db_select)
        self.bottom_layout.addWidget(self.db_table)

        # Molecule Info
        self.info_label = QLabel(self)
        self.main_layout.addWidget(self.info_label)

        # Current image (for saving)
        self.current_image = None

        # Populate database view
        self.populate_database_view()

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
                img.save(img_bytes, format='PNG')
                img_bytes.seek(0)

                qimage = QImage.fromData(img_bytes.getvalue())
                pixmap = QPixmap.fromImage(qimage)

                # Display the image on the QLabel
                self.image_label.setPixmap(pixmap)
                self.current_image = Image.open(io.BytesIO(img_bytes.getvalue()))

                # Fetch additional data from PubChem and display molecule info
                pubchem_data = fetch_pubchem_data(user_input)
                if pubchem_data:
                    molecular_formula = pubchem_data.get('molecular_formula', "N/A")
                    atomic_mass = pubchem_data.get('atomic_mass', "N/A")
                    iupac_name = pubchem_data.get('iupac_name', "N/A")

                    self.db_manager.add_molecule(user_input, atomic_mass, iupac_name, molecular_formula)
                    self.populate_database_view()

                    self.info_label.setText(f"IUPAC: {iupac_name}\nFormula: {molecular_formula}\nWeight: {atomic_mass} g/mol")

                    # Populate atom and bond info (can be modularized further if necessary)
                    self.populate_atom_bond_info(mol)

            except Exception as e:
                logging.error(f"Error processing molecule: {e}")
                QMessageBox.critical(self, "Error", f"Could not process the molecule: {e}")
        else:
            QMessageBox.warning(self, "Invalid Input", "Could not interpret input as a valid molecule.")

    def populate_atom_bond_info(self, mol):
        """Populate atom counts and bond information in the tree widgets."""
        self.atom_count_tree.clear()
        self.bond_info_tree.clear()

        # Analyze the molecule
        molecule_info = get_molecule_info(mol)
        atom_counts = molecule_info['atom_counts']
        bond_info = molecule_info['bond_info']

        # Populate atom counts
        for atom, count in atom_counts.items():
            atom_item = QTreeWidgetItem([atom, str(count)])
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

    def populate_database_view(self):
        """Populate the Treeview with SMILES from the database."""
        self.db_table.clear()

        molecules = self.db_manager.get_all_molecules()

        for molecule in molecules:
            item = QTreeWidgetItem(molecule)
            self.db_table.addTopLevelItem(item)

    def on_db_select(self, item, column):
        """Handle selection in the database view."""
        selected_smiles = item.text(0)
        self.input_field.setText(selected_smiles)
        self.draw_molecule()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_window = ChemDrawMain()
    main_window.show()
    sys.exit(app.exec_())
