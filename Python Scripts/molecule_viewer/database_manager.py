# database_manager.py

import sqlite3
import logging

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
                common_name TEXT,
                iupac_name TEXT,
                molecular_formula TEXT,
                atomic_mass REAL,
                starred INTEGER DEFAULT 0  -- Add a 'starred' column
            )
        ''')
        self.connection.commit()

    def upsert_molecule(self, smiles, data):
        try:
            self.cursor.execute("SELECT atomic_mass, iupac_name, molecular_formula, common_name, starred FROM molecules WHERE smiles = ?", (smiles,))
            existing_data = self.cursor.fetchone()

            if existing_data:
                atomic_mass, iupac_name, molecular_formula, common_name, starred = existing_data
                atomic_mass = data['atomic_mass'] if data['atomic_mass'] != 'N/A' and data['atomic_mass'] is not None else atomic_mass
                iupac_name = data['iupac_name'] if data['iupac_name'] != 'N/A' and data['iupac_name'] is not None else iupac_name
                molecular_formula = data['molecular_formula'] if data['molecular_formula'] != 'N/A' and data['molecular_formula'] is not None else molecular_formula
                common_name = data['common_name'] if data['common_name'] != 'N/A' and data['common_name'] is not None else common_name

                self.cursor.execute(
                    "UPDATE molecules SET atomic_mass = ?, iupac_name = ?, molecular_formula = ?, common_name = ?, starred = ? WHERE smiles = ?",
                    (atomic_mass, iupac_name, molecular_formula, common_name, starred, smiles)
                )
            else:
                self.cursor.execute(
                    "INSERT INTO molecules (smiles, atomic_mass, iupac_name, molecular_formula, common_name, starred) VALUES (?, ?, ?, ?, ?, ?)",
                    (smiles, data['atomic_mass'], data['iupac_name'], data['molecular_formula'], data['common_name'], 0)
                )
            self.connection.commit()

        except sqlite3.Error as e:
            logging.error(f"Database error: {e}")

    def star_molecule(self, smiles):
        """Toggle the star status of a molecule."""
        self.cursor.execute("UPDATE molecules SET starred = NOT starred WHERE smiles = ?", (smiles,))
        self.connection.commit()

    def get_all_molecules(self):
        self.cursor.execute("SELECT smiles, common_name, iupac_name, molecular_formula, atomic_mass, starred FROM molecules")
        return self.cursor.fetchall()

    def get_all_smiles(self):
        self.cursor.execute("SELECT smiles FROM molecules")
        return [row[0] for row in self.cursor.fetchall()]

    def close(self):
        self.connection.close()
