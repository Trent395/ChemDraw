import logging
import sqlite3
import random
import requests
from rdkit import Chem
from rdkit.Chem import rdmolops, inchi, Descriptors, AllChem
import io
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


def get_element_color(element):
    return color_table.get(element, color_table['Default'])

class DatabaseManager:
    def __init__(self, db_name='molecules.db'):
        self.db_name = db_name
        self.connection = sqlite3.connect(self.db_name)
        self.cursor = self.connection.cursor()
        self.create_table()

    def create_table(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS molecules (
                id INTEGER PRIMARY KEY,
                smiles TEXT UNIQUE,
                atomic_mass REAL,
                iupac_name TEXT
            )
        ''')
        self.connection.commit()

    def update_existing_table(self):
        try:
            self.cursor.execute("ALTER TABLE molecules ADD COLUMN atomic_mass REAL")
            logging.info("Added atomic_mass column to molecules table.")
        except sqlite3.OperationalError:
            logging.info("Column atomic_mass already exists.")

        try:
            self.cursor.execute("ALTER TABLE molecules ADD COLUMN iupac_name TEXT")
            logging.info("Added iupac_name column to molecules table.")
        except sqlite3.OperationalError:
            logging.info("Column iupac_name already exists.")

        self.connection.commit()

    def add_molecule(self, smiles, atomic_mass, iupac_name):
        try:
            self.cursor.execute("INSERT INTO molecules (smiles, atomic_mass, iupac_name) VALUES (?, ?, ?)", 
                                (smiles, atomic_mass, iupac_name))
            self.connection.commit()
        except sqlite3.IntegrityError:
            logging.warning(f"Molecule with SMILES {smiles} already exists in the database.")

    def add_smiles(self, smiles):
        try:
            self.cursor.execute("INSERT INTO molecules (smiles) VALUES (?)", (smiles,))
            self.connection.commit()
        except sqlite3.IntegrityError:
            logging.warning(f"SMILES {smiles} already exists in the database.")

    def get_all_smiles(self):
        self.cursor.execute("SELECT smiles FROM molecules")
        return [row[0] for row in self.cursor.fetchall()]

    def get_all_molecules(self):
        self.cursor.execute("SELECT smiles, atomic_mass, iupac_name FROM molecules")
        return self.cursor.fetchall()

    def close(self):
        self.connection.close()

def get_longest_carbon_chain(mol):
    max_length = 0
    longest_path = []

    for atom1 in mol.GetAtoms():
        if atom1.GetSymbol() == 'C': 
            for atom2 in mol.GetAtoms():
                if atom2.GetSymbol() == 'C' and atom1.GetIdx() != atom2.GetIdx():
                    path = rdmolops.GetShortestPath(mol, atom1.GetIdx(), atom2.GetIdx())
                    if len(path) > max_length:
                        max_length = len(path)
                        longest_path = path

    if max_length > 0 and max_length <= len(carbon_chain_names):
        chain_name = carbon_chain_names[max_length]
    else:
        chain_name = f"Chain of {max_length} carbons"

    return chain_name

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

def generate_random_smiles(selected_elements, selected_bond_types):
    bond_choices = [bond_symbols[bond] for bond in selected_bond_types if selected_bond_types[bond]]

    if not selected_elements or not bond_choices:
        return ""

    smiles = ''
    molecule_length = random.randint(3, 10)
    for _ in range(molecule_length):
        atom = random.choice(selected_elements)
        bond = random.choice(bond_choices)
        smiles += bond + atom

    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return smiles
    else:
        return ""

