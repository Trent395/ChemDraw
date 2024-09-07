# molecule_analyzer.py

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, QPushButton, QVBoxLayout, QWidget, QMessageBox
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class MoleculeHydrogenAnalyzer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecule Hydrogen Analyzer")
        self.setGeometry(100, 100, 500, 400)
        
        # Set up the main layout
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        layout = QVBoxLayout(self.central_widget)
        
        # Input for SMILES
        self.smiles_input = QLineEdit(self)
        self.smiles_input.setPlaceholderText("Enter SMILES string")
        layout.addWidget(self.smiles_input)
        
        # Calculate Button
        self.calculate_button = QPushButton("Analyze Molecule", self)
        self.calculate_button.clicked.connect(self.analyze_molecule)
        layout.addWidget(self.calculate_button)
        
        # Label for displaying the hydrogen count
        self.hydrogen_count_label = QLabel(self)
        layout.addWidget(self.hydrogen_count_label)
        
        # Label for showing molecule image
        self.molecule_image_label = QLabel(self)
        layout.addWidget(self.molecule_image_label)
        
    # Function to calculate hydrogen atoms based on your hypothesis
    def calculate_hydrogens(self, c, n, o, x, b):
        return 4 * c + 3 * n + 2 * o + x - 2 * b

    # Function to analyze molecule structure
    def analyze_molecule(self):
        smiles = self.smiles_input.text()
        if not smiles:
            QMessageBox.warning(self, "Input Error", "Please enter a valid SMILES string")
            return
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Getting the counts of atoms and bonds
                atom_counts = {'C': 0, 'N': 0, 'O': 0, 'X': 0}  # X for halogen atoms (generalized)
                bond_count = 0

                for atom in mol.GetAtoms():
                    symbol = atom.GetSymbol()
                    if symbol == 'C':
                        atom_counts['C'] += 1
                    elif symbol == 'N':
                        atom_counts['N'] += 1
                    elif symbol == 'O':
                        atom_counts['O'] += 1
                    elif symbol in ['F', 'Cl', 'Br', 'I']:  # Halogens
                        atom_counts['X'] += 1

                bond_count = len(mol.GetBonds())

                # Calculate expected hydrogen atoms
                expected_hydrogens = self.calculate_hydrogens(
                    atom_counts['C'], atom_counts['N'], atom_counts['O'], atom_counts['X'], bond_count
                )
                
                # Update hydrogen count label
                self.hydrogen_count_label.setText(f"Expected Hydrogen Atoms: {expected_hydrogens}")
                
                # Render the molecule
                AllChem.Compute2DCoords(mol)
                fig, ax = plt.subplots(figsize=(3, 3))
                ax.imshow(Draw.MolToImage(mol))
                ax.axis('off')
                
                # Create a canvas to show the plot
                canvas = FigureCanvas(fig)
                self.molecule_image_label.setPixmap(canvas.grab())  # Display image as pixmap
                plt.close(fig)
            else:
                QMessageBox.warning(self, "SMILES Error", "Invalid SMILES string")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MoleculeHydrogenAnalyzer()
    window.show()
    sys.exit(app.exec_())
