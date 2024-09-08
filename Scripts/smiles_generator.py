from rdkit import Chem
from rdkit.Chem import AllChem
import itertools
import logging
from collections import Counter
from PIL import ImageQt, Image  # Added ImageQt

class SMILESGenerator:
    def __init__(self, elements, counts, molecule_name="Unnamed Molecule"):
        """
        Initialize the SMILESGenerator with elements, their counts, and an optional molecule name.
        
        :param elements: List of element symbols (e.g., ['C', 'H', 'O'])
        :param counts: List of counts corresponding to the elements (e.g., [2, 6, 1])
        :param molecule_name: Optional name for the molecule (e.g., 'Ethanol')
        """
        self.elements = elements
        self.counts = counts
        self.molecule_name = molecule_name
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        # Create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)

        # Create formatter and add it to the handlers
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)

        # Add the handlers to the logger
        self.logger.addHandler(ch)
    
    def generate_combinations(self):
        """Generate unique permutations of atoms based on element counts."""
        atoms = [atom for atom, count in zip(self.elements, self.counts) for _ in range(count)]
        unique_combinations = set(itertools.permutations(atoms, len(atoms)))
        return unique_combinations

    def is_valid_molecule(self, mol):
        """Check if the molecule is valid by examining valences and connectivity."""
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
            return True
        except Exception as e:
            self.logger.debug(f"Invalid molecule: {e}")
            return False

    def build_molecule(self, atoms):
        """Attempt to build a molecule from a sequence of atoms."""
        mol = Chem.RWMol()
        atom_indices = []
        atom_valences = {"H": 1, "C": 4, "O": 2, "N": 3}  # Valences of common elements, add more as needed
        current_valence = {atom: 0 for atom in atoms}  # Track the current valence for each atom

        for atom in atoms:
            atom_idx = mol.AddAtom(Chem.Atom(atom))
            atom_indices.append(atom_idx)

        # Attempt to connect atoms without exceeding their valence
        for i in range(len(atom_indices) - 1):
            if current_valence[atoms[i]] < atom_valences[atoms[i]] and current_valence[atoms[i + 1]] < atom_valences[atoms[i + 1]]:
                mol.AddBond(atom_indices[i], atom_indices[i + 1], Chem.BondType.SINGLE)
                current_valence[atoms[i]] += 1
                current_valence[atoms[i + 1]] += 1

        if self.is_valid_molecule(mol):
            AllChem.Compute2DCoords(mol)
            return mol
        else:
            return None


    def generate_valid_smiles(self):
        """Generate all valid SMILES codes for molecules with specific atom counts."""
        valid_smiles = set()
        combinations = self.generate_combinations()
        
        for combo in combinations:
            mol = self.build_molecule(combo)
            if mol:
                try:
                    smiles = Chem.MolToSmiles(mol)
                    if '.' not in smiles:  # Simple filter to avoid fragments, but further checks are needed
                        valid_smiles.add(smiles)
                        self.logger.info(f"Generated SMILES: {smiles}")
                except Exception as e:
                    self.logger.error(f"Failed to generate SMILES: {e}")
        
        return list(valid_smiles)
    
    def process_smiles_list(self, valid_smiles_list):
        """Process each SMILES string, generate an image, and convert it to ImageQt."""
        images = []
        for smiles in valid_smiles_list:
            try:
                img = self.generate_image_from_smiles(smiles)
                if img:
                    qim = ImageQt(img)
                    images.append(qim)
            except Exception as e:
                logging.error(f"Failed to process SMILES {smiles}: {e}")
                continue
        return images

    
    def save_smiles_to_file(self, filename="generated_smiles.txt"):
        """Save the generated SMILES codes to a text file."""
        smiles_list = self.generate_valid_smiles()
        with open(filename, 'w') as file:
            file.write(f"SMILES codes for {self.molecule_name}:\n")
            for smiles in smiles_list:
                file.write(f"{smiles}\n")
        self.logger.info(f"SMILES codes saved to {filename}")

# Example usage
if __name__ == "__main__":
    # Example usage: Generate SMILES codes for molecules with 2 carbons and 6 hydrogens (ethane-like structures)
    elements = ['C', 'H']
    counts = [2, 6]
    generator = SMILESGenerator(elements, counts, molecule_name="Ethane")
    generator.save_smiles_to_file("ethane_smiles.txt")
