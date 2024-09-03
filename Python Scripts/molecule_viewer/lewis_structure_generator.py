# lewis_structure_generator.py

from rdkit import Chem
from PIL import Image, ImageDraw, ImageFont
import io
from PyQt5.QtGui import QPixmap, QImage

class LewisStructureGenerator:
    def __init__(self, smiles):
        self.mol = Chem.MolFromSmiles(smiles)
        Chem.SanitizeMol(self.mol)

    def draw_lewis_structure(self):
        # Create a blank canvas for the Lewis structure
        img = Image.new("RGBA", (400, 400), (255, 255, 255, 0))
        draw = ImageDraw.Draw(img)
        font = ImageFont.load_default()

        # Simple geometric placement based on bonds (this is a naive approach)
        atom_positions = self.get_atom_positions()
        bond_lines = self.get_bond_lines(atom_positions)

        # Draw bonds first
        for line in bond_lines:
            draw.line(line, fill="black", width=3)

        # Draw atoms
        for atom, pos in atom_positions.items():
            draw.text(pos, atom, fill="black", font=font)

        # Draw lone pairs (simplified)
        lone_pairs = self.calculate_lone_pairs(atom_positions)
        for pos in lone_pairs:
            draw.text(pos, "•", fill="black", font=font)

        return img

    def get_atom_positions(self):
        # Placeholder: a simplistic approach to position atoms
        atom_positions = {}
        for i, atom in enumerate(self.mol.GetAtoms()):
            x = 200 + 50 * (i % 4)
            y = 200 + 50 * (i // 4)
            atom_positions[atom.GetSymbol()] = (x, y)
        return atom_positions

    def get_bond_lines(self, atom_positions):
        # Placeholder: simple bond lines between atoms
        bond_lines = []
        for bond in self.mol.GetBonds():
            start_atom = self.mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
            end_atom = self.mol.GetAtomWithIdx(bond.GetEndAtomIdx())
            start_pos = atom_positions[start_atom.GetSymbol()]
            end_pos = atom_positions[end_atom.GetSymbol()]
            bond_lines.append([start_pos, end_pos])
        return bond_lines

    def calculate_lone_pairs(self, atom_positions):
        # Placeholder: a simplistic approach to calculate lone pairs
        lone_pairs = []
        for atom in self.mol.GetAtoms():
            valence = atom.GetTotalValence() + atom.GetTotalNumHs()
            # Simplified logic: assuming lone pairs need to fill to 8 valence electrons
            needed_lone_pairs = (8 - valence) // 2
            for i in range(needed_lone_pairs):
                offset = (10 * i, -10 * i)
                pos = (atom_positions[atom.GetSymbol()][0] + offset[0], atom_positions[atom.GetSymbol()][1] + offset[1])
                lone_pairs.append(pos)
        return lone_pairs

    def get_lewis_structure_pixmap(self):
        img = self.draw_lewis_structure()
        img_bytes = io.BytesIO()
        img.save(img_bytes, format='PNG')
        img_bytes.seek(0)

        qimage = QImage.fromData(img_bytes.getvalue())
        pixmap = QPixmap.fromImage(qimage)
        return pixmap
