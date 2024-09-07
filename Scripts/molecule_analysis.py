# molecule_analysis.py

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem
from PyQt5.QtGui import QPixmap, QImage
import io
import logging

class MoleculeAnalysis:
    def __init__(self, smiles):
        self.smiles = smiles
        self.mol = self.get_molecule()

    def get_molecule(self):
        if self.smiles and isinstance(self.smiles, str):
            mol = Chem.MolFromSmiles(self.smiles)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    AllChem.Compute2DCoords(mol)
                    return mol
                except Exception as e:
                    logging.error(f"Error sanitizing molecule: {e}")
        else:
            logging.error("Invalid SMILES input or empty string.")
        return None

    def get_pixmap(self):
        if self.mol:
            img = Draw.MolToImage(self.mol)
            img_bytes = io.BytesIO()
            img.save(img_bytes, format='PNG')
            img_bytes.seek(0)

            qimage = QImage.fromData(img_bytes.getvalue())
            pixmap = QPixmap.fromImage(qimage)
            return pixmap
        else:
            # Return an empty pixmap or a placeholder image if the molecule is None
            return QPixmap()

    def get_atom_bond_info(self):
        if self.mol:
            atom_counts = {}
            for atom in self.mol.GetAtoms():
                element = atom.GetSymbol()
                atom_counts[element] = atom_counts.get(element, 0) + 1

            bond_info = []
            electronegativity = {'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44} #Make this reference my elements script
            for bond in self.mol.GetBonds():
                atom1 = bond.GetBeginAtom()
                atom2 = bond.GetEndAtom()
                elem1 = atom1.GetSymbol()
                elem2 = atom2.GetSymbol()
                en_diff = abs(electronegativity.get(elem1, 0) - electronegativity.get(elem2, 0))
                bond_length = Chem.rdMolTransforms.GetBondLength(self.mol.GetConformer(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                bond_info.append((f"{elem1}-{elem2}", f"{en_diff:.2f}", f"{bond_length:.2f} Å"))

            return atom_counts, bond_info
        else:
            return {}, []

    def get_info_text(self, molecule_data):
        if molecule_data:
            return f"IUPAC: {molecule_data.get('iupac_name', 'N/A')}\nFormula: {molecule_data.get('molecular_formula', 'N/A')}\nWeight: {molecule_data.get('atomic_mass', 'N/A')} g/mol\nXLogP: {molecule_data.get('xlogp', 'N/A')}\nH-Bond Donors: {molecule_data.get('h_bond_donors', 'N/A')}\nH-Bond Acceptors: {molecule_data.get('h_bond_acceptors', 'N/A')}"
        else:
            return "No molecule data available."
