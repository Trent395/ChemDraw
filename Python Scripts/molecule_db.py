import os
from database_manager import DatabaseManager
import pubchempy as pcp

# Directory setup
MOLECULES_DB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "MoleculeDB")
os.makedirs(MOLECULES_DB_DIR, exist_ok=True)

# Schema specific for molecules
MOLECULES_SCHEMA = {
    'smiles': 'TEXT UNIQUE',
    'iupac_name': 'TEXT',
    'molecular_formula': 'TEXT',
    'atomic_mass': 'REAL'
}

class MoleculeDatabase:
    def __init__(self, db_name='molecules.db'):
        # Initialize the DatabaseManager with the molecules-specific schema
        self.db_manager = DatabaseManager(
            db_name=db_name,
            table_name='molecules',
            schema=MOLECULES_SCHEMA,
            directory=MOLECULES_DB_DIR
        )

    def add_molecule(self, smiles, iupac_name=None, molecular_formula=None, atomic_mass=None):
        record = {
            'smiles': smiles,
            'iupac_name': iupac_name,
            'molecular_formula': molecular_formula,
            'atomic_mass': atomic_mass
        }
        self.db_manager.add_record(record)

    def fetch_molecule(self, smiles):
        return self.db_manager.fetch_record(where_clause={'smiles': smiles})

    def update_molecule(self, smiles, iupac_name=None, molecular_formula=None, atomic_mass=None):
        record = {
            'iupac_name': iupac_name,
            'molecular_formula': molecular_formula,
            'atomic_mass': atomic_mass
        }
        self.db_manager.update_record(record, where_clause={'smiles': smiles})

    def fetch_pubchem_data(self, smiles):
        data = {'molecular_formula': 'N/A', 'atomic_mass': 'N/A', 'iupac_name': 'N/A', 'xlogp': 'N/A', 'h_bond_donors': 'N/A', 'h_bond_acceptors': 'N/A'}
        try:
            compound = pcp.get_compounds(smiles, namespace='smiles')
            if compound:
                compound = compound[0]
                data['molecular_formula'] = compound.molecular_formula or 'N/A'
                data['atomic_mass'] = compound.molecular_weight or 'N/A'
                data['iupac_name'] = compound.iupac_name or 'N/A'
                data['xlogp'] = compound.xlogp or 'N/A'

            mol = Chem.MolFromSmiles(smiles)
            if mol:
                data['h_bond_donors'] = Descriptors.NumHDonors(mol)
                data['h_bond_acceptors'] = Descriptors.NumHAcceptors(mol)

        except Exception as e:
            logging.error(f"Error fetching data for SMILES {smiles}: {e}")
        return data
    def update_all_molecules(self):
        molecules = self.db_manager.fetch_records()
        for molecule in molecules:
            smiles, iupac_name, molecular_formula, atomic_mass = molecule
            pubchem_data = self.fetch_pubchem_data(smiles)
            self.update_molecule(smiles, **pubchem_data)

    def get_all_molecules(self):
        """Fetch all molecules from the database."""
        return self.db_manager.fetch_records(columns='smiles, atomic_mass, iupac_name, molecular_formula')

    def close(self):
        self.db_manager.close()
