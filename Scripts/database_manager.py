# database_manager.py

import os
import sqlite3
import logging
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem

class DatabaseManager:
    def __init__(self, db_name='molecules.db'):
        """
        Initialize the DatabaseManager class and ensure the database is saved in the Databases folder.
        """
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # Ensure the Databases folder exists
        db_folder = os.path.join(root_dir, 'Databases')
        if not os.path.exists(db_folder):
            os.makedirs(db_folder)

        # Set the database file path in the Databases folder
        self.db_name = os.path.join(db_folder, db_name)

        # Connect to the database
        self.connection = sqlite3.connect(self.db_name)
        self.cursor = self.connection.cursor()
        self.create_table()

    def create_table(self):
        """
        Create the molecules table if it doesn't exist. Adds melting_point, boiling_point, and nickname fields.
        """
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS molecules (
                id INTEGER PRIMARY KEY,
                smiles TEXT UNIQUE,
                common_name TEXT,
                nickname TEXT,  -- New column for user-defined nickname
                iupac_name TEXT,
                molecular_formula TEXT,
                atomic_mass REAL,
                xlogp REAL,
                tpsa REAL,
                rotatable_bonds INTEGER,
                molecular_volume REAL,
                exact_mass REAL,
                complexity REAL,
                h_bond_donors INTEGER,
                h_bond_acceptors INTEGER,
                heavy_atom_count INTEGER,
                canonical_smiles TEXT,
                inchi TEXT,
                inchikey TEXT,
                isomeric_smiles TEXT,
                formal_charge INTEGER,
                charge INTEGER,
                num_rings INTEGER,
                num_aromatic_rings INTEGER,
                fraction_csp3 REAL,
                num_aliphatic_carbons INTEGER,
                num_aliphatic_heteroatoms INTEGER,
                num_aliphatic_rings INTEGER,
                num_stereocenters INTEGER,
                melting_point REAL,  -- New column for melting point
                boiling_point REAL,  -- New column for boiling point
                starred INTEGER DEFAULT 0
            )
        ''')
        self.connection.commit()

    def upsert_molecule(self, smiles, data):
        """
        Insert or update molecule data in the database, including melting_point, boiling_point, and nickname.
        """
        try:
            self.cursor.execute("SELECT id FROM molecules WHERE smiles = ?", (smiles,))
            existing_data = self.cursor.fetchone()

            if existing_data:
                self.cursor.execute('''
                    UPDATE molecules SET
                    common_name = ?, iupac_name = ?, molecular_formula = ?, atomic_mass = ?,
                    xlogp = ?, tpsa = ?, rotatable_bonds = ?, molecular_volume = ?,
                    exact_mass = ?, complexity = ?, h_bond_donors = ?, h_bond_acceptors = ?,
                    heavy_atom_count = ?, canonical_smiles = ?, inchi = ?, inchikey = ?,
                    isomeric_smiles = ?, formal_charge = ?, charge = ?, num_rings = ?,
                    num_aromatic_rings = ?, fraction_csp3 = ?, num_aliphatic_carbons = ?,
                    num_aliphatic_heteroatoms = ?, num_aliphatic_rings = ?, num_stereocenters = ?,
                    melting_point = ?, boiling_point = ?  -- Include new fields in update
                    WHERE smiles = ?
                ''', (
                    data['common_name'], data['iupac_name'], data['molecular_formula'], data['atomic_mass'],
                    data['xlogp'], data['tpsa'], data['rotatable_bonds'], data['molecular_volume'],
                    data['exact_mass'], data['complexity'], data['h_bond_donors'], data['h_bond_acceptors'],
                    data['heavy_atom_count'], data['canonical_smiles'], data['inchi'], data['inchikey'],
                    data['isomeric_smiles'], data['formal_charge'], data['charge'], data['num_rings'],
                    data['num_aromatic_rings'], data['fraction_csp3'], data['num_aliphatic_carbons'],
                    data['num_aliphatic_heteroatoms'], data['num_aliphatic_rings'], data['num_stereocenters'],
                    data['melting_point'], data['boiling_point'], smiles
                ))
            else:
                self.cursor.execute('''
                    INSERT INTO molecules (smiles, common_name, iupac_name, molecular_formula, atomic_mass,
                    xlogp, tpsa, rotatable_bonds, molecular_volume, exact_mass, complexity, h_bond_donors,
                    h_bond_acceptors, heavy_atom_count, canonical_smiles, inchi, inchikey, isomeric_smiles,
                    formal_charge, charge, num_rings, num_aromatic_rings, fraction_csp3, num_aliphatic_carbons,
                    num_aliphatic_heteroatoms, num_aliphatic_rings, num_stereocenters, melting_point, boiling_point)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    smiles, data['common_name'], data['iupac_name'], data['molecular_formula'], data['atomic_mass'],
                    data['xlogp'], data['tpsa'], data['rotatable_bonds'], data['molecular_volume'], data['exact_mass'],
                    data['complexity'], data['h_bond_donors'], data['h_bond_acceptors'], data['heavy_atom_count'],
                    data['canonical_smiles'], data['inchi'], data['inchikey'], data['isomeric_smiles'], data['formal_charge'],
                    data['charge'], data['num_rings'], data['num_aromatic_rings'], data['fraction_csp3'],
                    data['num_aliphatic_carbons'], data['num_aliphatic_heteroatoms'], data['num_aliphatic_rings'],
                    data['num_stereocenters'], data['melting_point'], data['boiling_point']
                ))

            self.connection.commit()
        except sqlite3.Error as e:
            logging.error(f"Database error: {e}")


    def fetch_pubchem_data(self, smiles):
        data = {
            'common_name': 'N/A',
            'molecular_formula': 'N/A',
            'atomic_mass': 'N/A',
            'iupac_name': 'N/A',
            'xlogp': 'N/A',
            'tpsa': 'N/A',
            'rotatable_bonds': 'N/A',
            'molecular_volume': 'N/A',
            'exact_mass': 'N/A',
            'complexity': 'N/A',
            'h_bond_donors': 'N/A',
            'h_bond_acceptors': 'N/A',
            'heavy_atom_count': 'N/A',
            'canonical_smiles': 'N/A',
            'inchi': 'N/A',
            'inchikey': 'N/A',
            'isomeric_smiles': 'N/A',
            'formal_charge': 'N/A',
            'charge': 'N/A',
            'num_rings': 'N/A',
            'num_aromatic_rings': 'N/A',
            'fraction_csp3': 'N/A',
            'num_aliphatic_carbons': 'N/A',
            'num_aliphatic_heteroatoms': 'N/A',
            'num_aliphatic_rings': 'N/A',
            'num_stereocenters': 'N/A'
        }

        try:
            compound = pcp.get_compounds(smiles, namespace='smiles')
            if not compound:
                logging.warning(f"No valid data found for SMILES: {smiles}")
                return data

            compound = compound[0]
            data['common_name'] = compound.synonyms[0] if compound.synonyms else 'N/A'
            data['molecular_formula'] = compound.molecular_formula or 'N/A'
            data['atomic_mass'] = compound.molecular_weight or 'N/A'
            data['iupac_name'] = compound.iupac_name or 'N/A'
            data['xlogp'] = compound.xlogp or 'N/A'
            data['tpsa'] = compound.tpsa or 'N/A'
            data['rotatable_bonds'] = compound.rotatable_bond_count or 'N/A'
            data['molecular_volume'] = compound.volume_3d or 'N/A'
            data['exact_mass'] = compound.exact_mass or 'N/A'
            data['complexity'] = compound.complexity or 'N/A'
            data['h_bond_donors'] = compound.h_bond_donor_count or 'N/A'
            data['h_bond_acceptors'] = compound.h_bond_acceptor_count or 'N/A'
            data['heavy_atom_count'] = compound.heavy_atom_count or 'N/A'
            data['canonical_smiles'] = compound.canonical_smiles or 'N/A'
            data['inchi'] = compound.inchi or 'N/A'
            data['inchikey'] = compound.inchikey or 'N/A'
            data['isomeric_smiles'] = compound.isomeric_smiles or 'N/A'
            data['formal_charge'] = compound.charge or 'N/A'
            data['charge'] = compound.charge or 'N/A'

            mol = Chem.MolFromSmiles(smiles)
            if mol:
                data['num_rings'] = Descriptors.RingCount(mol)
                data['num_aromatic_rings'] = Descriptors.NumAromaticRings(mol)
                data['fraction_csp3'] = Descriptors.FractionCSP3(mol)
                data['num_aliphatic_carbons'] = Descriptors.NumAliphaticCarbocycles(mol)
                data['num_aliphatic_heteroatoms'] = Descriptors.NumAliphaticHeterocycles(mol)
                data['num_aliphatic_rings'] = Descriptors.NumAliphaticRings(mol)
                data['num_stereocenters'] = self.count_stereocenters(mol)

        except Exception as e:
            logging.error(f"Error fetching data for SMILES {smiles}: {e}")
        
        return data

    def count_stereocenters(self, mol):
        chiral_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        return chiral_centers

    def upsert_molecule(self, smiles, data):
        try:
            self.cursor.execute("SELECT atomic_mass, common_name, iupac_name, molecular_formula, starred FROM molecules WHERE smiles = ?", (smiles,))
            existing_data = self.cursor.fetchone()

            if existing_data:
                atomic_mass, common_name, iupac_name, molecular_formula, starred = existing_data
                atomic_mass = data['atomic_mass'] if data['atomic_mass'] != 'N/A' and data['atomic_mass'] is not None else atomic_mass
                common_name = data['common_name'] if data['common_name'] != 'N/A' and data['common_name'] is not None else common_name
                iupac_name = data['iupac_name'] if data['iupac_name'] != 'N/A' and data['iupac_name'] is not None else iupac_name
                molecular_formula = data['molecular_formula'] if data['molecular_formula'] != 'N/A' and data['molecular_formula'] is not None else molecular_formula

                self.cursor.execute(
                    "UPDATE molecules SET atomic_mass = ?, common_name = ?, iupac_name = ?, molecular_formula = ?, starred = ? WHERE smiles = ?",
                    (atomic_mass, common_name, iupac_name, molecular_formula, starred, smiles)
                )
            else:
                self.cursor.execute(
                    "INSERT INTO molecules (smiles, common_name, iupac_name, molecular_formula, atomic_mass, starred) VALUES (?, ?, ?, ?, ?, ?)",
                    (smiles, data['common_name'], data['iupac_name'], data['molecular_formula'], data['atomic_mass'], 0)
                )
            self.connection.commit()

        except sqlite3.Error as e:
            logging.error(f"Database error: {e}")

    def get_all_molecules(self):
        self.cursor.execute("SELECT smiles, common_name, iupac_name, molecular_formula, atomic_mass, starred FROM molecules ORDER BY starred DESC, common_name ASC")
        return self.cursor.fetchall()
    
    def get_all_smiles(self):
        self.cursor.execute("SELECT smiles FROM molecules")
        return [row[0] for row in self.cursor.fetchall()]
    
    def upsert_important_smiles(self, smiles, common_name, starred=0):
        try:
            self.cursor.execute("SELECT id FROM important_smiles WHERE smiles = ?", (smiles,))
            existing_data = self.cursor.fetchone()

            if existing_data:
                self.cursor.execute(
                    "UPDATE important_smiles SET common_name = ?, starred = ? WHERE smiles = ?",
                    (common_name, starred, smiles)
                )
            else:
                self.cursor.execute(
                    "INSERT INTO important_smiles (smiles, common_name, starred) VALUES (?, ?, ?)",
                    (smiles, common_name, starred)
                )
            self.connection.commit()
        except sqlite3.Error as e:
            logging.error(f"Database error: {e}")

    def get_important_smiles(self):
        self.cursor.execute("SELECT smiles, common_name, starred FROM important_smiles ORDER BY starred DESC, common_name ASC")
        return self.cursor.fetchall()

    def update_all_database(self):
        smiles_list = self.get_all_smiles()
        for smiles in smiles_list:
            molecule_data = self.fetch_pubchem_data(smiles)
            self.upsert_molecule(smiles, molecule_data)
        logging.info("All database entries have been updated.")

    def clean_up_database(self):
        try:
            self.cursor.execute("DELETE FROM molecules WHERE smiles IS NULL OR smiles = 'None'")
            self.connection.commit()
            logging.info("Cleaned up invalid SMILES entries from the database.")
        except sqlite3.Error as e:
            logging.error(f"Error cleaning up the database: {e}")

    def close(self):
        self.connection.close()
