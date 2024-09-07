# hydrogen_calculator.py
import logging
from rdkit import Chem

# Fix to work for O=O
class HydrogenCalculator:
    def __init__(self):
        self.atom_counts = {'C': 0, 'N': 0, 'O': 0, 'F': 0, 'Cl': 0, 'Br': 0, 'I': 0}
        self.bond_count = 0
        self.hydrogen_count = 0

    def calculate_hydrogen_count(self, smiles):
        """Calculate hydrogen count using the 4C + 3N + 2O + X - 2B formula."""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError("Invalid SMILES string")

        self.atom_counts = {'C': 0, 'N': 0, 'O': 0, 'F': 0, 'Cl': 0, 'Br': 0, 'I': 0}
        self.bond_count = len(mol.GetBonds())  # Bonds between non-hydrogen atoms

        # Debug: Print atom and bond counts
        logging.info(f"SMILES: {smiles}")
        logging.info(f"Bond count: {self.bond_count}")

        # Count atoms in the molecule
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in self.atom_counts:
                self.atom_counts[symbol] += 1

        logging.info(f"Atom counts: {self.atom_counts}")

        # Extract variables
        c = self.atom_counts['C']
        n = self.atom_counts['N']
        o = self.atom_counts['O']
        x = self.atom_counts['F'] + self.atom_counts['Cl'] + self.atom_counts['Br'] + self.atom_counts['I']  # Halogens
        b = self.bond_count

        # Debug: Log the variables before hydrogen calculation
        logging.info(f"C: {c}, N: {n}, O: {o}, X: {x}, B: {b}")

        # Calculate hydrogens
        self.hydrogen_count = 4 * c + 3 * n + 2 * o + x - 2 * b
        logging.info(f"Calculated Hydrogens: {self.hydrogen_count}")

        return self.hydrogen_count


    def get_atom_counts(self):
        """Return the atom counts and bond count."""
        return self.atom_counts

    def get_bond_count(self):
        """Return the number of bonds between non-hydrogen atoms."""
        return self.bond_count

    def get_hydrogen_count(self):
        """Return the calculated hydrogen count."""
        return self.hydrogen_count
