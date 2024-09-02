import os
import sqlite3
import logging
from rdkit import Chem
from rdkit.Chem import Descriptors
import pubchempy as pcp

CHEMDRAW_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ChemDraw")
DATABASE_DIR = os.path.join(CHEMDRAW_DIR, "Databases")
os.makedirs(DATABASE_DIR, exist_ok=True)

class DatabaseManager:
    def __init__(self, db_name='molecules.db'):
        self.db_name = os.path.join(DATABASE_DIR, db_name)
        self.connection = sqlite3.connect(self.db_name)
        self.cursor = self.connection.cursor()
        self.create_table()
        self.update_existing_table()

    def create_table(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS molecules (
                id INTEGER PRIMARY KEY,
                smiles TEXT UNIQUE,
                iupac_name TEXT,
                molecular_formula TEXT,
                atomic_mass REAL
            )
        ''')
        self.connection.commit()

    def add_molecule(self, smiles, atomic_mass, iupac_name, molecular_formula):
        try:
            self.cursor.execute(
                "INSERT INTO molecules (smiles, iupac_name, molecular_formula, atomic_mass) VALUES (?, ?, ?, ?)",
                (smiles, iupac_name, molecular_formula, atomic_mass)
            )
            self.connection.commit()
        except sqlite3.IntegrityError:
            logging.warning(f"Molecule with SMILES {smiles} already exists in the database.")

    def update_existing_table(self):
        self.cursor.execute("PRAGMA table_info(molecules)")
        columns = [column[1] for column in self.cursor.fetchall()]

        if 'molecular_formula' not in columns:
            self.cursor.execute("ALTER TABLE molecules ADD COLUMN molecular_formula TEXT")
            self.connection.commit()
            logging.info("Added 'molecular_formula' column to 'molecules' table.")
        else:
            logging.info("'molecular_formula' column already exists.")

    def update_molecule(self, smiles, atomic_mass, iupac_name, molecular_formula):
        self.cursor.execute(
            "UPDATE molecules SET atomic_mass = ?, iupac_name = ?, molecular_formula = ? WHERE smiles = ?",
            (atomic_mass, iupac_name, molecular_formula, smiles)
        )
        self.connection.commit()

    def get_all_molecules(self):
        self.cursor.execute("SELECT smiles, atomic_mass, iupac_name, molecular_formula FROM molecules")
        return self.cursor.fetchall()

    def update_all_molecules(self):
        molecules = self.get_all_molecules()
        for molecule in molecules:
            smiles, atomic_mass, iupac_name, molecular_formula = molecule
            updated = False

            if iupac_name is None or molecular_formula is None or atomic_mass is None:
                pubchem_data = self.fetch_pubchem_data(smiles)
                if pubchem_data:
                    if iupac_name is None and pubchem_data['iupac_name'] != "N/A":
                        iupac_name = pubchem_data['iupac_name']
                        updated = True
                    if molecular_formula is None and pubchem_data['molecular_formula'] != "N/A":
                        molecular_formula = pubchem_data['molecular_formula']
                        updated = True
                    if atomic_mass is None and pubchem_data['atomic_mass'] != "N/A":
                        atomic_mass = pubchem_data['atomic_mass']
                        updated = True

                if updated:
                    self.update_molecule(smiles, atomic_mass, iupac_name, molecular_formula)

    def fetch_pubchem_data(self, smiles):
        try:
            compound = pcp.get_compounds(smiles, namespace='smiles')
            if compound:
                return {
                    'molecular_formula': compound[0].molecular_formula or 'N/A',
                    'atomic_mass': compound[0].molecular_weight or 'N/A',
                    'iupac_name': compound[0].iupac_name or 'N/A',
                }
            else:
                logging.warning(f"No data found for SMILES {smiles}.")
                return {'molecular_formula': 'N/A', 'atomic_mass': 'N/A', 'iupac_name': 'N/A'}
        except Exception as e:
            logging.error(f"Error fetching data for SMILES {smiles}: {e}")
            return {'molecular_formula': 'N/A', 'atomic_mass': 'N/A', 'iupac_name': 'N/A'}

    def close(self):
        self.connection.close()
