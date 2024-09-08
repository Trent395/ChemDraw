from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem
from PyQt5.QtGui import QPixmap, QImage
import io
import logging
from elements import Elements  # Import the Elements class for atomic properties

class MoleculeAnalysis:
    def __init__(self, smiles):
        self.smiles = smiles
        self.mol = self.get_molecule()

    # SMARTS patterns for a wide range of functional groups
    # Updated SMARTS patterns for a wide range of functional groups
        self.functional_groups = {
            "Alcohol": "[OX2H]",  # -OH (hydroxyl)
            "Aldehyde": "[CX3H1](=O)[#6]",  # Updated -CHO (specific for aldehyde)
            "Ketone": "[CX3](=O)[C][C]",  # Updated -C=O with two separate carbon atoms
            "Carboxylic Acid": "[CX3](=O)[OX2H1]",  # -COOH
            "Ester": "[CX3](=O)[OX2][#6]",  # -COOR
            "Amine": "[NX3;H2,H1;!$(NC=O)]",  # -NH2
            "Amide": "[NX3][CX3](=O)[#6]",  # -CONH2 or similar (R-C=O-N)
            "Ether": "[OD2]([#6])[#6]",  # -O-
            "Nitrile": "[NX1]#[CX2]",  # -C≡N
            "Alkyne": "[CX2]#C",  # Updated -C≡C- (corrected pattern for alkyne)
            "Alkene": "[CX3]=[CX3]",  # C=C
            "Haloalkane": "[CX4][F,Cl,Br,I]",  # C-X (halogen)
            "Thiol": "[#16X2H]",  # -SH (sulfhydryl)
            "Sulfide": "[#16X2]",  # -S-
            "Disulfide": "[#16X2][#16X2]",  # -S-S-
            "Sulfoxide": "[#16X2][#16X1](=O)",  # -S(=O)-
            "Sulfonyl": "[#16X4](=O)(=O)[#6]",  # -SO2-
            "Imine": "[CX2]=[NX2]",  # -C=N
            "Isocyanate": "[NX2]=[CX2]=[OX1]",  # -N=C=O
            "Isothiocyanate": "[NX2]=[CX2]=[SX1]",  # -N=C=S
            "Hydrazine": "[NX3][NX3]",  # -NH-NH-
            "Azide": "[NX1]=[NX2]=[NX1]",  # -N=N=N
            "Diazo": "[NX2]=[NX2+]=[CX3-]",  # -N=N-
            "Anhydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",  # -CO-O-CO-
            "Peroxide": "[OX2][OX2]",  # -O-O-
            "Acetal": "[CX4]([OX2R])[OX2R]",  # R-O-R
            "Enol": "[OX2H][CX3]=[CX3]",  # -OH attached to a C=C
            "Carbamate": "[NX3][CX3](=[OX1])[OX2]",  # -NHCOO-
            "Phosphate": "[OX1]=[PX4](O)(O)O",  # -PO4
            "Phosphonate": "[OX1]=[PX3](O)(O)[CX4]",  # -PO3
            "Nitro": "[NX3](=O)[OX1]",  # -NO2
            "Nitroso": "[NX2]=[OX1]",  # -NO
            "Sulfonamide": "[SX4](=O)(=O)[NX3]",  # -SO2NH2
            "Sulfonate": "[SX4](=O)(=O)[OX2][#6]",  # -SO3-
            "Carbodiimide": "[NX2]=[CX2]=[NX2]",  # -N=C=N-
            "Azo": "[NX2]=[NX2]",  # -N=N-
            "Cyanate": "[NX1]=[CX2]=[OX1]",  # -OCN
            "Thiocyanate": "[NX1]=[CX2]=[SX1]",  # -SCN
        }

    def detect_functional_groups(self):
        """
        Detects functional groups in the molecule using predefined SMARTS patterns.
        Returns a list of detected functional groups.
        """
        if not self.mol:
            return []

        detected_groups = []
        for group_name, smarts_pattern in self.functional_groups.items():
            try:
                pattern = Chem.MolFromSmarts(smarts_pattern)
                if pattern is None:
                    logging.error(f"Invalid SMARTS pattern for {group_name}: {smarts_pattern}")
                    continue

                if self.mol.HasSubstructMatch(pattern):
                    detected_groups.append(group_name)
            except Exception as e:
                logging.error(f"Error matching {group_name} with SMARTS {smarts_pattern}: {e}")
                continue

        # Log detected functional groups
        logging.info(f"Detected functional groups in {self.smiles}: {detected_groups}")
        return detected_groups

    def get_bond_polarity(self, en_diff):
        """Classify the bond type based on electronegativity difference (ΔEN)."""
        if en_diff < 0.5:
            return "Non-Polar"
        elif 0.5 <= en_diff < 1.7:
            return "Polar"
        else:
            return "Ionic"

    def get_molecule(self):
        """Generate the molecule object from SMILES."""
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
        """Generate a QPixmap image for display."""
        if self.mol:
            img = Draw.MolToImage(self.mol)
            img_bytes = io.BytesIO()
            img.save(img_bytes, format='PNG')
            img_bytes.seek(0)

            qimage = QImage.fromData(img_bytes.getvalue())
            pixmap = QPixmap.fromImage(qimage)
            return pixmap
        else:
            return QPixmap()  # Return an empty pixmap if no molecule
    def get_atom_bond_info(self):
        """
        Get atom counts and bond information.
        Return atom counts including atomic mass and electronegativity, 
        and bond information including bond length and polarity.
        """
        if self.mol:
            atom_counts = {}
            bond_info = []

            try:
                elements = Elements()  # Create an instance of the Elements class
                
                # Loop through atoms to group them by element and count occurrences, including hydrogen
                for atom in self.mol.GetAtoms():
                    element = atom.GetSymbol()
                    if element in atom_counts:
                        atom_counts[element]['count'] += 1
                    else:
                        # Fetch atomic mass and electronegativity from Elements class
                        atomic_mass = elements.get_atomic_mass(element)
                        electronegativity = elements.get_electronegativity(element)
                        atom_counts[element] = {
                            'count': 1,
                            'atomic_mass': f"{atomic_mass:.2f}",
                            'electronegativity': electronegativity
                        }
                
                # Loop through bonds to get bond information
                for bond in self.mol.GetBonds():
                    atom1 = bond.GetBeginAtom().GetSymbol()
                    atom2 = bond.GetEndAtom().GetSymbol()
                    
                    # Get bond length
                    bond_length = Chem.rdMolTransforms.GetBondLength(
                        self.mol.GetConformer(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                    )
                    
                    # Calculate electronegativity difference and bond polarity
                    en_diff = abs(elements.get_electronegativity(atom1) - elements.get_electronegativity(atom2))
                    bond_type = self.get_bond_polarity(en_diff)
                    
                    # Ensure that bond_info has exactly 4 values per bond tuple
                    bond_info.append((f"{atom1}-{atom2}", f"{en_diff:.2f}", f"{bond_length:.2f} Å", bond_type))
                    
                return atom_counts, bond_info
            
            except Exception as e:
                logging.error(f"Error calculating atom/bond info: {e}")
                return {}, []
        else:
            logging.error("No valid molecule loaded for atom/bond info.")
            return {}, []


    def get_info_text(self, molecule_data):
        """
        Returns a string with all molecule data from the database, including functional groups.
        """   
        if molecule_data:
            # Detect functional groups
            functional_groups = self.detect_functional_groups()
            functional_group_text = ", ".join(functional_groups) if functional_groups else "None"

            # Build the info text dynamically from the database fields
            info_lines = []
            for key, value in molecule_data.items():
                # If the key is "functional_groups", append our detected groups
                if key == 'functional_groups':
                    continue  # Avoid appending it here, as we'll add it separately below
                # Skip any keys that are None or 'N/A'
                if value is not None and value != 'N/A':
                    # Format the key to a readable form and append to the list
                    info_lines.append(f"{key.replace('_', ' ').title()}: {value}")
            
            # Add functional groups to the info
            info_lines.append(f"Functional Groups: {functional_group_text}")

            # Join the info lines into a single string
            return "\n".join(info_lines)

        else:
            return "No molecule data available."
