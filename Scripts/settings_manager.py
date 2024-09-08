# settings_manager.py

import sqlite3
import os
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QCheckBox
from log_viewer import LogViewer

class SettingsManager(QDialog):
    """
    Manages saving and loading application settings using SQLite in the 'settings.db' database.
    Provides a GUI for managing settings like dark mode and viewing logs.
    """

    def __init__(self, db_name='settings.db', parent=None):
        """
        Initialize the SettingsManager, ensure the settings table exists in the SQLite database,
        and create a settings GUI with a button to view logs.
        """
        super().__init__(parent)

        # Database path
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        db_folder = os.path.join(root_dir, 'Databases')
        if not os.path.exists(db_folder):
            os.makedirs(db_folder)

        self.db_path = os.path.join(db_folder, db_name)
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

        # Create settings table
        self.create_settings_table()

        # Load settings
        self.settings = self.load_settings()

        # Set up the GUI
        self.setWindowTitle("Application Settings")
        self.setGeometry(300, 300, 400, 250)
        self.layout = QVBoxLayout()

        # Dark Mode Checkbox
        self.dark_mode_checkbox = QCheckBox("Enable Dark Mode", self)
        self.dark_mode_checkbox.setChecked(self.get_setting('dark_mode', 'False') == 'True')
        self.layout.addWidget(self.dark_mode_checkbox)

        # Buttons for saving and viewing logs
        button_layout = QHBoxLayout()

        save_button = QPushButton("Save", self)
        save_button.clicked.connect(self.save_and_close)
        button_layout.addWidget(save_button)

        logs_button = QPushButton("View Logs", self)
        logs_button.clicked.connect(self.open_log_viewer)
        button_layout.addWidget(logs_button)

        cancel_button = QPushButton("Cancel", self)
        cancel_button.clicked.connect(self.close)
        button_layout.addWidget(cancel_button)

        self.layout.addLayout(button_layout)
        self.setLayout(self.layout)

    def create_settings_table(self):
        """Create the settings table in the database."""
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS settings (
                key TEXT PRIMARY KEY,
                value TEXT
            )
        ''')
        self.connection.commit()

    def load_settings(self):
        """Load all settings from the database."""
        self.cursor.execute('SELECT key, value FROM settings')
        settings = {key: value for key, value in self.cursor.fetchall()}
        return settings

    def save_setting(self, key, value):
        """Save a specific setting in the database."""
        self.cursor.execute('''
            INSERT OR REPLACE INTO settings (key, value) VALUES (?, ?)
        ''', (key, value))
        self.connection.commit()

    def get_setting(self, key, default=None):
        """Get a specific setting by key."""
        return self.settings.get(key, default)

    def save_and_close(self):
        """Save the dark mode setting and close the window."""
        self.save_setting('dark_mode', 'True' if self.dark_mode_checkbox.isChecked() else 'False')
        self.connection.commit()
        self.close()

    def open_log_viewer(self):
        """Open the log viewer dialog."""
        log_viewer = LogViewer(parent=self)
        log_viewer.exec_()

    def closeEvent(self, event):
        """Close database connection on window close."""
        self.connection.close()
        event.accept()
