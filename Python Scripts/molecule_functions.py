import logging
import sqlite3
import random
import requests
from rdkit import Chem
from rdkit.Chem import rdmolops, Descriptors
from PIL import Image

# Set up logging
logging.basicConfig(filename='molecule_viewer.log', level=logging.INFO, format='%(asctime)s - %(message)s')

bond_symbols = {
    'Single': '',     # Single bond: No symbol needed
    'Double': '=',    # Double bond: '=' symbol
    'Triple': '#'     # Triple bond: '#' symbol
}

color_table = {
    'H': 'White',
    'C': 'Gray',
    'N': 'Blue',
    'O': 'Red',
    'F': 'Green',
    'Cl': 'Green',
    'Br': 'Brown',
    'Default': 'Gray'
}

carbon_chain_names = {
    1: "Methane",
    2: "Ethane",
    3: "Propane",
    4: "Butane",
    5: "Pentane",
    6: "Hexane",
    7: "Heptane",
    8: "Octane",
    9: "Nonane",
    10: "Decane",
    11: "Undecane",
    12: "Dodecane"
}

def generate_random_smiles(num_carbons, num_hydrogens=None):
    if num_carbons < 1:
        return ""

    # Initialize an empty molecule with carbon atoms
    mol = Chem.RWMol()
    carbon_indices = []

    # Add the specified number of carbon atoms
    for _ in range(num_carbons):
        carbon_idx = mol.AddAtom(Chem.Atom('C'))
        carbon_indices.append(carbon_idx)

    # Add bonds between carbons, ensuring no more than 3 bonds between any two carbons
    for i in range(num_carbons - 1):
        bond_type = Chem.rdchem.BondType.SINGLE
        mol.AddBond(carbon_indices[i], carbon_indices[i + 1], bond_type)

    # Ensure maximum 4 bonds on any one carbon, fill with hydrogens as needed
    for carbon_idx in carbon_indices:
        attached_atoms = mol.GetAtomWithIdx(carbon_idx).GetTotalNumHs()
        num_bonds = mol.GetAtomWithIdx(carbon_idx).GetDegree()

        # Calculate remaining bonds available for this carbon
        remaining_bonds = 4 - num_bonds - attached_atoms

        if num_hydrogens is not None:
            # If the number of hydrogens is specified, distribute them
            hydrogens_to_add = min(remaining_bonds, num_hydrogens)
            num_hydrogens -= hydrogens_to_add
        else:
            # If not specified, add hydrogens to saturate the carbon
            hydrogens_to_add = remaining_bonds

        for _ in range(hydrogens_to_add):
            mol.AddAtom(Chem.Atom('H'))
            mol.AddBond(carbon_idx, mol.GetNumAtoms() - 1, Chem.rdchem.BondType.SINGLE)

    # Update property cache to ensure correct valence and hydrogen information
    mol.UpdatePropertyCache(strict=False)

    # Check if any hydrogens are still unallocated
    if num_hydrogens is not None and num_hydrogens > 0:
        for carbon_idx in carbon_indices:
            if num_hydrogens <= 0:
                break
            attached_atoms = mol.GetAtomWithIdx(carbon_idx).GetTotalNumHs()
            num_bonds = mol.GetAtomWithIdx(carbon_idx).GetDegree()
            remaining_bonds = 4 - num_bonds - attached_atoms

            if remaining_bonds > 0:
                mol.AddAtom(Chem.Atom('H'))
                mol.AddBond(carbon_idx, mol.GetNumAtoms() - 1, Chem.rdchem.BondType.SINGLE)
                num_hydrogens -= 1

    # Remove hydrogens and return the SMILES string of the molecule
    mol = Chem.RemoveHs(mol)
    smiles = Chem.MolToSmiles(mol)
    return smiles


def get_element_color(element):
    return color_table.get(element, color_table['Default'])

def get_molecule_info(mol):
    mol_with_hs = Chem.AddHs(mol)
    mol_weight = Descriptors.MolWt(mol_with_hs)
    atom_counts = {}
    atom_masses = {}
    for atom in mol_with_hs.GetAtoms():
        elem = atom.GetSymbol()
        atom_counts[elem] = atom_counts.get(elem, 0) + 1
        atom_masses[elem] = atom_masses.get(elem, 0) + atom.GetMass()
    total_atoms = sum(atom_counts.values())
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol_with_hs)
    info_text = f"Molecular Formula: {formula}\nMolecular Weight: {mol_weight:.2f} g/mol\nTotal Atoms: {total_atoms}"
    log_text = f"Formula: {formula}, Weight: {mol_weight}, Total Atoms: {total_atoms}"
    
    return {'info_text': info_text, 'log_text': log_text}

def get_longest_carbon_chain(mol):
    def dfs(atom, visited, length):
        visited.add(atom.GetIdx())
        max_len = length
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited and neighbor.GetSymbol() == 'C':
                max_len = max(max_len, dfs(neighbor, visited, length + 1))
        visited.remove(atom.GetIdx())
        return max_len

    max_length = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            max_length = max(max_length, dfs(atom, set(), 1))
    
    return max_length

def fetch_pubchem_data(smiles):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        return data
    except requests.exceptions.RequestException as e:
        logging.error(f"Error fetching data from PubChem: {e}")
        return None

def display_molecule_info(self, mol):
        # Add explicit hydrogens to the molecule
        mol_with_hs = Chem.AddHs(mol)

        # Calculate molecular weight and number of atoms
        mol_weight = Descriptors.MolWt(mol_with_hs)
        atom_counts = {}
        atom_masses = {}
        electronegativity = {
            'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
            'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'I': 2.66
        }
        for atom in mol_with_hs.GetAtoms():
            elem = atom.GetSymbol()
            atom_counts[elem] = atom_counts.get(elem, 0) + 1
            atom_masses[elem] = atom_masses.get(elem, 0) + atom.GetMass()
        total_atoms = sum(atom_counts.values())

        # Molecular formula
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol_with_hs)

        # Calculate molar mass contributions and percentages
        molar_mass_contributions = {elem: mass for elem, mass in atom_masses.items()}
        molar_mass_percentages = {
            elem: (mass / mol_weight) * 100 for elem, mass in molar_mass_contributions.items()
        }

        # Detect functional groups
        functional_groups = []
        smarts_patterns = {
            'Alcohol': '[CX4][OH]',
            'Amine': '[NX3][CX4]',
            'Amide': '[CX3](=[OX1])[NX3]',
            'Ester': '[CX3](=[OX1])[OX2H0]',
            'Carboxylic Acid': '[CX3](=O)[OX2H1]',
            'Aromatic': 'c1ccccc1',
            'Alkene': 'C=C',
            'Alkyne': 'C#C'
        }
        for group, smarts in smarts_patterns.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                functional_groups.append(group)

        # Populate Atom Count Table
        for i in self.atom_count_tree.get_children():
            self.atom_count_tree.delete(i)  # Clear previous entries

        for elem, count in atom_counts.items():
            color = get_element_color(elem)
            self.atom_count_tree.insert("", "end", values=(elem, count, f"{molar_mass_contributions[elem]:.2f}", f"{molar_mass_percentages[elem]:.2f}%"), tags=(elem,))

        # Configure tags for Treeview based on color table
        for elem, color_name in color_table.items():
            self.atom_count_tree.tag_configure(elem, foreground=color_name.lower())

        # Populate Functional Groups Table
        for i in self.functional_groups_tree.get_children():
            self.functional_groups_tree.delete(i)  # Clear previous entries

        for group in functional_groups:
            self.functional_groups_tree.insert("", "end", values=(group,))

        # Populate Bond Information Table
        bond_info = []
        for i in self.bond_info_tree.get_children():
            self.bond_info_tree.delete(i)  # Clear previous entries

        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            elem1 = atom1.GetSymbol()
            elem2 = atom2.GetSymbol()
            if elem1 in electronegativity and elem2 in electronegativity:
                en_diff = abs(electronegativity[elem1] - electronegativity[elem2])
                bond_length = Chem.rdMolTransforms.GetBondLength(mol_with_hs.GetConformer(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                bond_info.append((f"{elem1}-{elem2}", f"{en_diff:.2f}", f"{bond_length:.2f} Å"))
        
        # Insert bond information with color tags
        for bond in bond_info:
            bond_tag = bond[0].split('-')[0]  # Get the first element of the bond
            self.bond_info_tree.insert("", "end", values=bond, tags=(bond_tag,))

        # Configure tags for Treeview based on color table
        for elem, color_name in color_table.items():
            self.bond_info_tree.tag_configure(elem, foreground=color_name)

        # Display additional molecule information
        info_text = f"Molecular Formula: {formula}\nMolecular Weight: {mol_weight:.2f} g/mol\nTotal Atoms: {total_atoms}\n"
        self.info_label.config(text=info_text)
        logging.info(f"Molecule Info - Formula: {formula}, Weight: {mol_weight}, Total Atoms: {total_atoms}, Atom Count: {atom_counts}, Groups: {functional_groups}")
