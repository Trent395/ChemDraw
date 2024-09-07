import os
import sqlite3
import logging

class DatabaseManager:
    def __init__(self, db_name='database.db', table_name='data', schema=None, directory=None):
        """
        Initialize the database manager.

        :param db_name: Name of the database file.
        :param table_name: Name of the table to operate on.
        :param schema: Dictionary representing the table schema. Format: {'column_name': 'column_type'}
        :param directory: Directory where the database file will be stored. If None, the current directory is used.
        """
        if directory:
            os.makedirs(directory, exist_ok=True)
            self.db_name = os.path.join(directory, db_name)
        else:
            self.db_name = db_name

        self.table_name = table_name
        self.connection = sqlite3.connect(self.db_name)
        self.cursor = self.connection.cursor()

        if schema:
            self.create_table(schema)

    def create_table(self, schema):
        """
        Create a table with the given schema.

        :param schema: Dictionary representing the table schema. Format: {'column_name': 'column_type'}
        """
        columns_definition = ', '.join(f"{col_name} {col_type}" for col_name, col_type in schema.items())
        self.cursor.execute(f'''
            CREATE TABLE IF NOT EXISTS {self.table_name} (
                id INTEGER PRIMARY KEY,
                {columns_definition}
            )
        ''')
        self.connection.commit()

    def add_record(self, record):
        """
        Add a new record to the table.

        :param record: Dictionary representing the record to add. Format: {'column_name': value}
        """
        placeholders = ', '.join('?' * len(record))
        columns = ', '.join(record.keys())
        values = tuple(record.values())
        
        try:
            self.cursor.execute(
                f"INSERT INTO {self.table_name} ({columns}) VALUES ({placeholders})", values
            )
            self.connection.commit()
        except sqlite3.IntegrityError as e:
            logging.warning(f"Failed to add record {record}: {e}")

    def update_record(self, record, where_clause):
        """
        Update an existing record in the table.

        :param record: Dictionary representing the updated values. Format: {'column_name': value}
        :param where_clause: Dictionary representing the condition to identify the record(s) to update. Format: {'column_name': value}
        """
        columns = ', '.join(f"{key} = ?" for key in record.keys())
        values = tuple(record.values())
        where_clause_str = ' AND '.join(f"{key} = ?" for key in where_clause.keys())
        where_values = tuple(where_clause.values())
        
        try:
            self.cursor.execute(
                f"UPDATE {self.table_name} SET {columns} WHERE {where_clause_str}", values + where_values
            )
            self.connection.commit()
        except sqlite3.Error as e:
            logging.warning(f"Failed to update record {record}: {e}")

    def fetch_records(self, columns='*', where_clause=None):
        """
        Fetch records from the table.

        :param columns: Columns to fetch, separated by commas, or '*' for all columns.
        :param where_clause: Dictionary representing the condition to filter the records. Format: {'column_name': value}
        :return: List of tuples containing the fetched records.
        """
        where_clause_str = ''
        where_values = ()
        if where_clause:
            where_clause_str = 'WHERE ' + ' AND '.join(f"{key} = ?" for key in where_clause.keys())
            where_values = tuple(where_clause.values())
        
        self.cursor.execute(f"SELECT {columns} FROM {self.table_name} {where_clause_str}", where_values)
        return self.cursor.fetchall()

    def fetch_record(self, columns='*', where_clause=None):
        """
        Fetch a single record from the table.

        :param columns: Columns to fetch, separated by commas, or '*' for all columns.
        :param where_clause: Dictionary representing the condition to filter the record. Format: {'column_name': value}
        :return: Tuple containing the fetched record, or None if no record is found.
        """
        records = self.fetch_records(columns=columns, where_clause=where_clause)
        return records[0] if records else None

    def update_schema(self, schema_updates):
        """
        Dynamically update the table schema by adding new columns if they don't already exist.

        :param schema_updates: Dictionary representing the new columns to add. Format: {'column_name': 'column_type'}
        """
        existing_columns = self.get_columns()
        for column_name, column_type in schema_updates.items():
            if column_name not in existing_columns:
                self.cursor.execute(f"ALTER TABLE {self.table_name} ADD COLUMN {column_name} {column_type}")
                self.connection.commit()
                logging.info(f"Added '{column_name}' column to '{self.table_name}' table.")

    def get_columns(self):
        """
        Get the list of columns in the table.

        :return: List of column names.
        """
        self.cursor.execute(f"PRAGMA table_info({self.table_name})")
        return [column[1] for column in self.cursor.fetchall()]

    def close(self):
        """
        Close the database connection.
        """
        self.connection.close()
