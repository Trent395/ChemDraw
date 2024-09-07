# hydrogen_calculator.py

import logging
from rdkit import Chem

class HydrogenCalculator:
    """
    A class for calculating hydrogen count and determining the number of bonds (b) between non-hydrogen atoms in a molecule.
    """

    def __init__(self):
        """
        Initialize the HydrogenCalculator, set default atom counts, and bond/hydrogen counts.
        """
        # Initialize atom counts for common elements
        self.atom_counts = {'C': 0, 'N': 0, 'O': 0, 'F': 0, 'Cl': 0, 'Br': 0, 'I': 0}
        self.bond_count = 0
        self.hydrogen_count = 0

    def calculate_hydrogen_count(self, smiles):
        """
        Calculate hydrogen count using the 4C + 3N + 2O + X - 2B formula, where:
        C = number of carbon atoms, N = number of nitrogen atoms, O = number of oxygen atoms,
        X = number of halogen atoms, and B = number of bonds between non-hydrogen atoms.

        Args:
            smiles (str): The SMILES code for the molecule.

        Returns:
            int: The calculated number of hydrogen atoms.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError("Invalid SMILES string")

        # Reset atom counts and bond count
        self.atom_counts = {'C': 0, 'N': 0, 'O': 0, 'F': 0, 'Cl': 0, 'Br': 0, 'I': 0}
        self.bond_count = 0  # Bonds between non-hydrogen atoms

        # Count atoms in the molecule
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in self.atom_counts:
                self.atom_counts[symbol] += 1

        # Count bonds between non-hydrogen atoms
        self.bond_count = self.calculate_bonds_between_non_hydrogens(mol)

        # Debug: Log atom and bond counts
        logging.info(f"SMILES: {smiles}")
        logging.info(f"Bond count: {self.bond_count}")
        logging.info(f"Atom counts: {self.atom_counts}")

        # Extract variables for the formula
        c = self.atom_counts['C']
        n = self.atom_counts['N']
        o = self.atom_counts['O']
        x = self.atom_counts['F'] + self.atom_counts['Cl'] + self.atom_counts['Br'] + self.atom_counts['I']  # Halogens
        b = self.bond_count

        # Debug: Log the variables before hydrogen calculation
        logging.info(f"C: {c}, N: {n}, O: {o}, X: {x}, B: {b}")

        # Calculate hydrogens using the formula 4C + 3N + 2O + X - 2B
        self.hydrogen_count = 4 * c + 3 * n + 2 * o + x - 2 * b
        logging.info(f"Calculated Hydrogens: {self.hydrogen_count}")

        return self.hydrogen_count

    def calculate_bonds_between_non_hydrogens(self, mol):
        """
        Calculate the number of bonds between non-hydrogen atoms.
        
        Args:
            mol (rdkit.Chem.Mol): The RDKit molecule object.

        Returns:
            int: The number of bonds between non-hydrogen atoms.
        """
        bond_count = 0

        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()

            # Ignore hydrogen atoms
            if atom1.GetAtomicNum() != 1 and atom2.GetAtomicNum() != 1:
                bond_type = bond.GetBondType()

                # Add bond count based on bond type
                if bond_type == Chem.BondType.SINGLE:
                    bond_count += 1
                elif bond_type == Chem.BondType.DOUBLE:
                    bond_count += 2
                elif bond_type == Chem.BondType.TRIPLE:
                    bond_count += 3
                elif bond_type == Chem.BondType.AROMATIC:
                    bond_count += 1  # Aromatic bonds count as single bonds for simplicity

        return bond_count

    def get_atom_counts(self):
        """
        Return the atom counts for C, N, O, F, Cl, Br, I.

        Returns:
            dict: The atom counts.
        """
        return self.atom_counts

    def get_bond_count(self):
        """
        Return the number of bonds between non-hydrogen atoms.

        Returns:
            int: The number of bonds (b).
        """
        return self.bond_count

    def get_hydrogen_count(self):
        """
        Return the calculated number of hydrogen atoms.

        Returns:
            int: The calculated hydrogen count.
        """
        return self.hydrogen_count


# Example usage
if __name__ == "__main__":
    # Example SMILES strings
    smiles_input = 'C#C'  # Acetylene (ethyne)
    
    calculator = HydrogenCalculator()
    hydrogen_count = calculator.calculate_hydrogen_count(smiles_input)
    bond_count = calculator.get_bond_count()

    print(f"SMILES: {smiles_input}")
    print(f"Number of bonds between non-hydrogen atoms (b): {bond_count}")
    print(f"Calculated number of hydrogen atoms: {hydrogen_count}")
