import sqlite3

class SQLiteMerger:
    def __init__(self, new_db_name):
        """Initialize the SQLiteMerger with the name of the new database."""
        self.new_db_name = new_db_name
        self.conn = sqlite3.connect(new_db_name)
        self.cursor = self.conn.cursor()

    def attach_database(self, db_file, alias):
        """Attach an external SQLite database to the current connection."""
        self.cursor.execute(f"ATTACH DATABASE '{db_file}' AS {alias}")

    def detach_database(self, alias):
        """Detach an attached SQLite database from the current connection."""
        self.cursor.execute(f"DETACH DATABASE '{alias}'")

    def get_table_row_count(self, alias, table_name):
        """Get the number of entries in a specific table of the attached database."""
        self.cursor.execute(f"SELECT COUNT(*) FROM {alias}.{table_name}")
        return self.cursor.fetchone()[0]

    def get_overlap_count(self, table_name, primary_key='id'):
        """Calculate the number of overlapping entries between the new database and the attached database."""
        self.cursor.execute(f"""
        SELECT COUNT(*)
        FROM {table_name} AS new
        INNER JOIN to_merge.{table_name} AS old
        ON new.{primary_key} = old.{primary_key}
        """)
        return self.cursor.fetchone()[0]

    def merge_table(self, alias, table_name):
        """Merge a table from the attached database into the new database."""
        # Create the table in the new database if it doesn't exist
        self.cursor.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} AS 
        SELECT * FROM {alias}.{table_name} WHERE 0
        """)
        
        # Insert data from the attached table to the new table
        self.cursor.execute(f"INSERT INTO {table_name} SELECT * FROM {alias}.{table_name}")

    def merge_database(self, db_file, alias='to_merge'):
        """Attach, merge all tables, and then detach a database."""
        self.attach_database(db_file, alias)
        
        # Get the list of tables in the attached database
        self.cursor.execute(f"SELECT name FROM {alias}.sqlite_master WHERE type='table'")
        tables = self.cursor.fetchall()
        
        for table_name in tables:
            row_count = self.get_table_row_count(alias, table_name[0])
            overlap_count = self.get_overlap_count(table_name[0])
            print(f"Table: {table_name[0]}, Rows: {row_count}, Overlaps: {overlap_count}")
            self.merge_table(alias, table_name[0])
        
        self.detach_database(alias)

    def merge_multiple_databases(self, db_files):
        """Merge multiple databases into the new database."""
        for db_file in db_files:
            self.merge_database(db_file)
    
    def commit_and_close(self):
        """Commit changes and close the database connection."""
        self.conn.commit()
        self.conn.close()
        print(f"Merge complete and database '{self.new_db_name}' saved.")

    def get_merge_info(self, db_files):
        """Get the merge info (number of rows and overlaps) for the databases."""
        info = ""
        for db_file in db_files:
            self.attach_database(db_file, 'to_merge')
            
            self.cursor.execute(f"SELECT name FROM to_merge.sqlite_master WHERE type='table'")
            tables = self.cursor.fetchall()
            
            for table_name in tables:
                row_count = self.get_table_row_count('to_merge', table_name[0])
                overlap_count = self.get_overlap_count(table_name[0])
                info += f"Database: {db_file} | Table: {table_name[0]} | Rows: {row_count} | Overlaps: {overlap_count}\n"
            
            self.detach_database('to_merge')
        return info

# Standalone usage
if __name__ == '__main__':
    db_files = ['molecules-DESKTOP-TRENT.db', 'molecules-TRENT-FRAMEWORK.db', 'molecules.db']
    merger = SQLiteMerger('merged_molecules.db')
    print(merger.get_merge_info(db_files))
    merger.merge_multiple_databases(db_files)
    merger.commit_and_close()
