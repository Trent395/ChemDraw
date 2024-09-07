import os
import sqlite3
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QCheckBox, QPushButton, QApplication
from dark_mode import DarkMode  # Import DarkMode class

class Settings(QDialog):
    def __init__(self, parent=None, base_dir=None, db_name='settings.db'):
        super(Settings, self).__init__(parent)
        
        # Set up the window
        self.setWindowTitle("Settings")
        self.setGeometry(300, 300, 200, 200)
        
        # Set up base directory, defaulting to a general "Settings" folder in the user's home directory
        if base_dir is None:
            base_dir = os.path.join(os.path.expanduser("~"), "Settings")
        os.makedirs(base_dir, exist_ok=True)
        
        # Set up database directory and connection
        self.db_name = os.path.join(base_dir, db_name)
        self.connection = sqlite3.connect(self.db_name)
        self.cursor = self.connection.cursor()
        self.create_table()

        # UI Setup
        layout = QVBoxLayout()

        # Dark Mode Checkbox
        self.dark_mode_checkbox = QCheckBox("Enable Dark Mode", self)
        self.dark_mode_checkbox.setChecked(self.get_setting('dark_mode') == 'True')
        layout.addWidget(self.dark_mode_checkbox)

        # Apply Button
        apply_button = QPushButton("Apply", self)
        apply_button.clicked.connect(self.apply_settings)
        layout.addWidget(apply_button)

        self.setLayout(layout)

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

    def apply_settings(self):
        # Save dark mode setting
        dark_mode_enabled = self.dark_mode_checkbox.isChecked()
        self.save_setting('dark_mode', str(dark_mode_enabled))
        
        # Apply dark mode using DarkMode class
        if dark_mode_enabled:
            DarkMode.apply(self.parent())
        else:
            self.parent().setStyleSheet("")  # Reset to default style

        self.accept()

    def closeEvent(self, event):
        self.connection.close()
        event.accept()

# Example usage within a PyQt5 application
if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QMainWindow, QPushButton

    class MainApp(QMainWindow):
        def __init__(self):
            super(MainApp, self).__init__()
            self.initUI()

        def initUI(self):
            self.setWindowTitle("Main App")
            self.setGeometry(100, 100, 400, 300)

            settings_button = QPushButton("Settings", self)
            settings_button.setGeometry(50, 50, 100, 30)
            settings_button.clicked.connect(self.open_settings)

            # Apply initial dark mode based on saved settings
            self.apply_initial_settings()

        def apply_initial_settings(self):
            # Example usage: Pass a custom base directory if needed
            settings_manager = Settings(self, base_dir=os.path.join(os.path.expanduser("~"), "MyAppSettings"))
            dark_mode_enabled = settings_manager.get_setting('dark_mode') == 'True'
            if dark_mode_enabled:
                DarkMode.apply(self)

        def open_settings(self):
            # Example usage: Pass a custom base directory if needed
            settings_dialog = Settings(self, base_dir=os.path.join(os.path.expanduser("~"), "MyAppSettings"))
            settings_dialog.exec_()

    app = QApplication(sys.argv)
    main_app = MainApp()
    main_app.show()
    sys.exit(app.exec_())
