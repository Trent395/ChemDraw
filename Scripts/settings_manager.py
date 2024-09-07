# settings_manager.py

import sqlite3
import os
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QCheckBox, QLabel

class SettingsManager(QDialog):
    """
    Manages saving and loading application settings using SQLite in the 'settings.db' database.
    Provides a GUI for managing settings like dark mode.
    """

    def __init__(self, db_name='settings.db', parent=None):
        """
        Initialize the SettingsManager, ensure the settings table exists in the SQLite database,
        and create a settings GUI.

        Args:
        db_name (str): The name of the database (default is 'settings.db').
        parent: The parent window to inherit the current theme.
        """
        super().__init__(parent)

        # Ensure the Databases folder exists
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        db_folder = os.path.join(root_dir, 'Databases')
        if not os.path.exists(db_folder):
            os.makedirs(db_folder)

        # Set the database file path in the Databases folder
        self.db_path = os.path.join(db_folder, db_name)

        # Connect to the database
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

        # Create the settings table if it doesn't exist
        self.create_settings_table()

        # Load settings from the database
        self.settings = self.load_settings()

        # Create GUI elements
        self.setWindowTitle("Application Settings")
        self.setGeometry(300, 300, 400, 200)

        # Main layout
        self.layout = QVBoxLayout()

        # Dark Mode Toggle
        self.dark_mode_checkbox = QCheckBox("Enable Dark Mode", self)
        self.dark_mode_checkbox.setChecked(self.get_setting('dark_mode', 'True') == 'True')
        self.layout.addWidget(self.dark_mode_checkbox)

        # Save and Cancel buttons
        button_layout = QHBoxLayout()

        save_button = QPushButton("Save", self)
        save_button.clicked.connect(self.save_and_close)
        button_layout.addWidget(save_button)

        cancel_button = QPushButton("Cancel", self)
        cancel_button.clicked.connect(self.close)
        button_layout.addWidget(cancel_button)

        self.layout.addLayout(button_layout)
        self.setLayout(self.layout)

        # Inherit the parent's theme (dark mode or light mode)
        if parent:
            self.setPalette(parent.palette())

    def create_settings_table(self):
        """
        Create the settings table in the database if it doesn't exist.
        """
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS settings (
                key TEXT PRIMARY KEY,
                value TEXT
            )
        ''')
        self.connection.commit()

    def load_settings(self):
        """
        Load all settings from the database.
        Returns the settings as a dictionary.
        """
        self.cursor.execute('SELECT key, value FROM settings')
        settings = {key: value for key, value in self.cursor.fetchall()}
        return settings

    def save_setting(self, key, value):
        """
        Save a setting in the database (insert or update).
        Args:
            key (str): The setting key.
            value (str): The setting value.
        """
        self.cursor.execute('''
            INSERT OR REPLACE INTO settings (key, value) VALUES (?, ?)
        ''', (key, value))
        self.connection.commit()

    def get_setting(self, key, default=None):
        """
        Get a specific setting value by key from the database.
        Args:
            key (str): The setting key.
            default: The default value to return if the key isn't found.
        Returns:
            The setting value if it exists, otherwise the default value.
        """
        return self.settings.get(key, default)

    def set_setting(self, key, value):
        """
        Set a specific setting value in the database.
        Args:
            key (str): The setting key.
            value: The new value for the setting.
        """
        self.save_setting(key, value)

    def save_and_close(self):
        """
        Save the settings and close the settings window.
        """
        # Save the dark mode setting to the database
        self.set_setting('dark_mode', 'True' if self.dark_mode_checkbox.isChecked() else 'False')

        self.connection.commit()
        self.close()

    def closeEvent(self, event):
        """
        Override the close event to ensure the database connection is closed properly.
        """
        self.connection.close()
        event.accept()
