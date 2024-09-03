# molecule_analysis.py

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem
import io
from PIL import Image
from PyQt5.QtGui import QPixmap, QImage

class MoleculeAnalysis:
    def __init__(self, smiles):
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)
        Chem.SanitizeMol(self.mol)
        AllChem.Compute2DCoords(self.mol)

    def get_molecule_pixmap(self):
        img = Draw.MolToImage(self.mol)
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        img_bytes.seek(0)

        qimage = QImage.fromData(img_bytes.getvalue())
        pixmap = QPixmap.fromImage(qimage)
        return pixmap

    def get_atom_bond_info(self):
        atom_counts = {}
        for atom in self.mol.GetAtoms():
            element = atom.GetSymbol()
            atom_counts[element] = atom_counts.get(element, 0) + 1

        bond_info = []
        electronegativity = {'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44}
        for bond in self.mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            elem1 = atom1.GetSymbol()
            elem2 = atom2.GetSymbol()
            en_diff = abs(electronegativity.get(elem1, 0) - electronegativity.get(elem2, 0))
            bond_length = Chem.rdMolTransforms.GetBondLength(self.mol.GetConformer(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            bond_info.append((f"{elem1}-{elem2}", f"{en_diff:.2f}", f"{bond_length:.2f} Å"))

        return atom_counts, bond_info

    def get_info_text(self, molecule_data):
        return f"IUPAC: {molecule_data['iupac_name']}\nFormula: {molecule_data['molecular_formula']}\nWeight: {molecule_data['atomic_mass']} g/mol\nXLogP: {molecule_data['xlogp']}\nH-Bond Donors: {molecule_data['h_bond_donors']}\nH-Bond Acceptors: {molecule_data['h_bond_acceptors']}"
