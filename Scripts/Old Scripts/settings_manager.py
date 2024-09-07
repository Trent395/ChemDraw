import os
import sqlite3

CHEMDRAW_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ChemDraw")
DATABASE_DIR = os.path.join(CHEMDRAW_DIR, "Databases")
os.makedirs(DATABASE_DIR, exist_ok=True)

class SettingsManager:
    def __init__(self, db_name='settings.db'):
        self.db_name = os.path.join(DATABASE_DIR, db_name)
        self.connection = sqlite3.connect(self.db_name)
        self.cursor = self.connection.cursor()
        self.create_table()

    def create_table(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS settings (
                id INTEGER PRIMARY KEY,
                setting_name TEXT UNIQUE,
                setting_value TEXT
            )
        ''')
        self.connection.commit()

    def save_setting(self, setting_name, setting_value):
        self.cursor.execute('''
            INSERT OR REPLACE INTO settings (setting_name, setting_value)
            VALUES (?, ?)
        ''', (setting_name, setting_value))
        self.connection.commit()

    def get_setting(self, setting_name):
        self.cursor.execute('SELECT setting_value FROM settings WHERE setting_name = ?', (setting_name,))
        result = self.cursor.fetchone()
        return result[0] if result else None

    def close(self):
        self.connection.close()
