from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, 
    QLineEdit, QCheckBox, QTreeWidget, QTreeWidgetItem, QGroupBox, QFormLayout, 
    QFileDialog, QMessageBox, QSpinBox
)
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem
import io
import logging
import molecule_functions as mf  # Importing the modular functions

class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Organic Chemistry Molecule Viewer")
        self.setGeometry(100, 100, 1000, 800)

        self.initUI()

    def initUI(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)

        # Input section
        input_layout = QHBoxLayout()
        layout.addLayout(input_layout)

        self.input_label = QLabel("Enter SMILES, InChI, or common name:")
        input_layout.addWidget(self.input_label)

        self.input_entry = QLineEdit()
        input_layout.addWidget(self.input_entry)

        self.draw_button = QPushButton("Draw Molecule")
        self.draw_button.clicked.connect(self.draw_molecule)
        input_layout.addWidget(self.draw_button)

        # Element selection
        element_group = QGroupBox("Include Elements:")
        element_layout = QVBoxLayout(element_group)
        self.element_checkboxes = {}
        elements = ['C', 'O', 'N', 'F', 'Cl', 'Br']
        for elem in elements:
            checkbox = QCheckBox(elem)
            checkbox.setChecked(True)
            self.element_checkboxes[elem] = checkbox
            element_layout.addWidget(checkbox)

        layout.addWidget(element_group)

        # Bond type selection
        bond_group = QGroupBox("Include Bonds:")
        bond_layout = QVBoxLayout(bond_group)
        self.bond_checkboxes = {}
        bond_types = ['Single', 'Double', 'Triple']
        for bond_type in bond_types:
            checkbox = QCheckBox(bond_type)
            checkbox.setChecked(True)
            self.bond_checkboxes[bond_type] = checkbox
            bond_layout.addWidget(checkbox)

        layout.addWidget(bond_group)

        # Molecule generation section
        generate_layout = QHBoxLayout()
        layout.addLayout(generate_layout)

        self.num_carbons_label = QLabel("Number of Carbons:")
        generate_layout.addWidget(self.num_carbons_label)
        self.num_carbons_input = QSpinBox()
        self.num_carbons_input.setMinimum(1)
        generate_layout.addWidget(self.num_carbons_input)

        self.num_hydrogens_label = QLabel("Number of Hydrogens:")
        generate_layout.addWidget(self.num_hydrogens_label)
        self.num_hydrogens_input = QSpinBox()
        self.num_hydrogens_input.setMinimum(0)
        generate_layout.addWidget(self.num_hydrogens_input)

        self.generate_button = QPushButton("Generate Random Molecule")
        self.generate_button.clicked.connect(self.generate_molecule)
        generate_layout.addWidget(self.generate_button)

        # Molecule display section
        self.canvas_label = QLabel()
        self.canvas_label.setFixedSize(400, 400)
        layout.addWidget(self.canvas_label)

        # Molecule information display
        self.info_label = QLabel()
        layout.addWidget(self.info_label)

        self.api_data_label = QLabel()
        layout.addWidget(self.api_data_label)

        self.longest_chain_label = QLabel()
        layout.addWidget(self.longest_chain_label)

        # Database view section
        self.db_tree = QTreeWidget()
        self.db_tree.setHeaderLabels(["SMILES"])
        layout.addWidget(self.db_tree)
        self.db_tree.itemClicked.connect(self.on_db_select)

    def generate_random_smiles(self):
        """Generate a random SMILES string based on selected elements, bonds, and user input."""
        num_carbons = self.num_carbons_input.value()
        num_hydrogens = self.num_hydrogens_input.value()

        try:
            smiles = mf.generate_random_smiles(num_carbons, num_hydrogens)
            logging.info(f"Generated SMILES: {smiles}")
            return smiles
        except Exception as e:
            logging.error(f"Error generating SMILES: {e}")
            QMessageBox.critical(self, "Generation Error", f"An error occurred while generating SMILES: {e}")
            return None

    def generate_molecule(self):
        smiles = self.generate_random_smiles()
        if not smiles:
            return

        self.input_entry.setText(smiles)
        self.draw_molecule()

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

    def display_molecule_info(self, mol):
        try:
            mol_with_hs = Chem.AddHs(mol)
            mol_weight = Descriptors.MolWt(mol_with_hs)
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol_with_hs)
            total_atoms = sum(1 for atom in mol_with_hs.GetAtoms())

            info_text = f"Molecular Formula: {formula}\nMolecular Weight: {mol_weight:.2f} g/mol\nTotal Atoms: {total_atoms}"
            self.info_label.setText(info_text)
            logging.info(f"Molecule Info - Formula: {formula}, Weight: {mol_weight}, Total Atoms: {total_atoms}")

            longest_chain_name = mf.get_longest_carbon_chain(mol)
            self.longest_chain_label.setText(f"Longest Carbon Chain: {longest_chain_name}")
        except Exception as e:
            logging.error(f"Error displaying molecule info: {e}")
            QMessageBox.critical(self, "Error", f"Could not display molecule information: {e}")

    def on_db_select(self, item):
        """Handle selection in the database Treeview."""
        selected_item = self.db_tree.currentItem()
        if selected_item:
            selected_smiles = selected_item.text(0)
            self.input_entry.setText(selected_smiles)
            self.draw_molecule()

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    app = QApplication([])
    viewer = MoleculeViewer()
    viewer.show()
    app.exec_()
