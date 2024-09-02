import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QGridLayout, QComboBox, QPushButton, QRadioButton, QButtonGroup, QMessageBox
from PyQt5.QtGui import QPixmap, QImage, QDesktopServices
from PyQt5.QtCore import Qt, QUrl
import sqlite3
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors

class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.conn = sqlite3.connect("molecules.db")
        self.cursor = self.conn.cursor()
        self.initialize_database()

    def initUI(self):
        self.setWindowTitle("Molecule Viewer")
        self.setGeometry(100, 100, 1200, 800)

        # Main layout
        layout = QVBoxLayout()

        # Dropdown for selecting carbon count
        self.carbon_count_label = QLabel("Select Carbon Count (x):")
        layout.addWidget(self.carbon_count_label)
        self.carbon_count_dropdown = QComboBox()
        self.carbon_count_dropdown.addItems([str(i) for i in range(1, 21)])  # For x = 1 to 20
        layout.addWidget(self.carbon_count_dropdown)

        # Radio buttons for molecule type
        self.molecule_type_label = QLabel("Select Molecule Type:")
        layout.addWidget(self.molecule_type_label)
        self.molecule_type_group = QButtonGroup(self)
        self.linear_radio = QRadioButton("Linear")
        self.cyclic_radio = QRadioButton("Cyclic")
        self.molecule_type_group.addButton(self.linear_radio)
        self.molecule_type_group.addButton(self.cyclic_radio)
        layout.addWidget(self.linear_radio)
        layout.addWidget(self.cyclic_radio)
        self.linear_radio.setChecked(True)

        # Radio buttons for branching type
        self.branching_type_label = QLabel("Select Branching Type:")
        layout.addWidget(self.branching_type_label)
        self.branching_type_group = QButtonGroup(self)
        self.unbranched_radio = QRadioButton("Unbranched")
        self.branched_radio = QRadioButton("Branched")
        self.branching_type_group.addButton(self.unbranched_radio)
        self.branching_type_group.addButton(self.branched_radio)
        layout.addWidget(self.unbranched_radio)
        layout.addWidget(self.branched_radio)
        self.unbranched_radio.setChecked(True)

        # Button to generate and display molecules
        self.generate_button = QPushButton("Generate Molecules")
        self.generate_button.clicked.connect(self.generate_and_display_molecules)
        layout.addWidget(self.generate_button)

        # Grid layout for molecule images
        self.grid_layout = QGridLayout()
        layout.addLayout(self.grid_layout)

        # Set the main layout to the central widget
        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def initialize_database(self):
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS molecules (
                id INTEGER PRIMARY KEY,
                smiles TEXT UNIQUE,
                formula TEXT,
                iupac_name TEXT
            )
        """)
        self.conn.commit()

    def generate_and_display_molecules(self):
        # Clear previous grid
        for i in reversed(range(self.grid_layout.count())):
            widget = self.grid_layout.itemAt(i).widget()
            if widget is not None:
                widget.setParent(None)

        x = int(self.carbon_count_dropdown.currentText())
        is_cyclic = self.cyclic_radio.isChecked()
        is_branched = self.branched_radio.isChecked()

        # Generate structures
        structures = self.generate_structures(x, is_cyclic, is_branched)

        row, col = 0, 0
        for smiles, formula, iupac_name in structures:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))
                img_qt = self.pil_image_to_qimage(img)
                pixmap = QPixmap.fromImage(img_qt)
                img_label = QLabel(self)
                img_label.setPixmap(pixmap)
                img_label.setAlignment(Qt.AlignCenter)

                # Make image clickable
                img_label.mousePressEvent = lambda event, url=f"https://pubchem.ncbi.nlm.nih.gov/#query={smiles}": self.open_pubchem_url(url)

                # Add image and details to grid
                self.grid_layout.addWidget(img_label, row, col)

                # Add formula and IUPAC name below the image
                details_label = QLabel(f"Formula: {formula}\nIUPAC: {iupac_name}")
                details_label.setAlignment(Qt.AlignCenter)
                self.grid_layout.addWidget(details_label, row + 1, col)

                col += 1
                if col > 3:  # Assuming 4 columns per row
                    col = 0
                    row += 2

    def pil_image_to_qimage(self, pil_image):
        """Convert PIL Image to QImage."""
        buffer = pil_image.tobytes("raw", "RGB")
        qimage = QImage(buffer, pil_image.size[0], pil_image.size[1], QImage.Format_RGB888)
        return qimage

    def generate_structures(self, x, is_cyclic, is_branched):
        """Generate a list of tuples containing SMILES, formula, and IUPAC names for x carbon atoms."""
        structures = []

        if is_cyclic:
            # Generate cyclic structures
            if x >= 3:
                mol = Chem.MolFromSmiles(f"C1{'C' * (x-1)}1")
                if mol:
                    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                    iupac_name = Chem.MolToName(mol)
                    smiles = Chem.MolToSmiles(mol)
                    structures.append((smiles, formula, iupac_name))
                    self.store_molecule(smiles, formula, iupac_name)
        else:
            if is_branched:
                # Generate branched structures using RDKit's enumeration tools
                mol = Chem.MolFromSmiles(f"C{'C' * (x - 1)}")
                if mol:
                    # Add branching logic (example for simplicity)
                    if x == 5:
                        # Example branch for pentane
                        branched_smiles = "CC(C)CC"
                        branched_mol = Chem.MolFromSmiles(branched_smiles)
                        if branched_mol:
                            formula = Chem.rdMolDescriptors.CalcMolFormula(branched_mol)
                            iupac_name = Chem.MolToName(branched_mol)
                            smiles = Chem.MolToSmiles(branched_mol)
                            structures.append((smiles, formula, iupac_name))
                            self.store_molecule(smiles, formula, iupac_name)
            else:
                # Generate linear unbranched structure
                mol = Chem.MolFromSmiles(f"C{'C' * (x - 1)}")
                if mol:
                    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                    iupac_name = Chem.MolToName(mol)
                    smiles = Chem.MolToSmiles(mol)
                    structures.append((smiles, formula, iupac_name))
                    self.store_molecule(smiles, formula, iupac_name)

        return structures

    def store_molecule(self, smiles, formula, iupac_name):
        """Store molecule data in the database."""
        try:
            self.cursor.execute("INSERT INTO molecules (smiles, formula, iupac_name) VALUES (?, ?, ?)",
                                (smiles, formula, iupac_name))
            self.conn.commit()
        except sqlite3.IntegrityError:
            pass  # Ignore duplicates

    def open_pubchem_url(self, url):
        """Open the PubChem page for the given SMILES."""
        QDesktopServices.openUrl(QUrl(url))

    def closeEvent(self, event):
        self.conn.close()


def main():
    app = QApplication(sys.argv)
    window = MoleculeViewer()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
