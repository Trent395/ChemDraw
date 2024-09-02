import sys
import os
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QFileDialog, QMessageBox, QHBoxLayout
from PyQt5.QtCore import Qt
from theme_manager import ThemeManager
import sqlite3

class MainApp(QWidget):
    def __init__(self):
        super().__init__()
        self.theme_manager = ThemeManager()
        self.initUI()
        self.load_settings()

    def initUI(self):
        self.setWindowTitle('Script Launcher')
        self.setGeometry(100, 100, 400, 300)

        self.layout = QVBoxLayout(self)
        self.layout.setAlignment(Qt.AlignTop)

        self.theme_manager.apply_theme(self)

        # Layout for buttons
        self.script_layout = QVBoxLayout()
        self.layout.addLayout(self.script_layout)

        # Layout for add/remove buttons
        self.bottom_layout = QHBoxLayout()

        # Add Script Button (Plus Sign)
        self.add_button = QPushButton('+', self)
        self.add_button.setFixedSize(50, 50)
        self.add_button.setStyleSheet("""
            QPushButton {
                font-size: 24px;
                text-align: center;
                color: black;
            }
        """ + self.theme_manager.get_button_style())
        self.add_button.clicked.connect(self.add_script)
        self.bottom_layout.addWidget(self.add_button)

        # Remove Script Button (Minus Sign)
        self.remove_button = QPushButton('-', self)
        self.remove_button.setFixedSize(50, 50)
        self.remove_button.setStyleSheet("""
            QPushButton {
                font-size: 24px;
                text-align: center;
                color: black;
            }
        """ + self.theme_manager.get_button_style())
        self.remove_button.clicked.connect(self.remove_script)
        self.bottom_layout.addWidget(self.remove_button)

        self.layout.addLayout(self.bottom_layout)

        # Theme Settings Button (small box in the corner)
        self.theme_button = QPushButton('', self)
        self.theme_button.setFixedSize(20, 20)
        self.theme_button.setStyleSheet(self.theme_manager.get_button_style())
        self.theme_button.clicked.connect(lambda: self.theme_manager.open_theme_settings(self))
        self.theme_button.move(370, 270)  # Place it at the bottom right corner
        self.layout.addWidget(self.theme_button)
        
        self.load_settings()

    def add_script(self):
        options = QFileDialog.Options()
        script, _ = QFileDialog.getOpenFileName(self, "Select Python Script", "", "Python Files (*.py);;All Files (*)", options=options)
        if script:
            script_name = os.path.basename(script).replace('.py', '')
            button = QPushButton(script_name, self)
            button.setStyleSheet(self.theme_manager.get_button_style())
            button.clicked.connect(lambda: self.run_script(script))
            self.script_layout.addWidget(button)
            self.save_settings()

    def remove_script(self):
        # Remove selected script
        if self.script_layout.count() > 0:
            script_button = self.script_layout.itemAt(self.script_layout.count() - 1).widget()
            if script_button:
                script_button.setParent(None)
            self.save_settings()

    def run_script(self, script_path):
        try:
            os.system(f'python "{script_path}"')
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to run script: {str(e)}")

    def load_settings(self):
        try:
            conn = sqlite3.connect('settings.db')
            c = conn.cursor()
            c.execute('''CREATE TABLE IF NOT EXISTS scripts (name TEXT, path TEXT)''')
            c.execute('SELECT name, path FROM scripts')
            rows = c.fetchall()
            for name, path in rows:
                button = QPushButton(name, self)
                button.setStyleSheet(self.theme_manager.get_button_style())
                button.clicked.connect(lambda _, p=path: self.run_script(p))
                self.script_layout.addWidget(button)
            conn.close()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load settings: {str(e)}")

    def save_settings(self):
        try:
            conn = sqlite3.connect('settings.db')
            c = conn.cursor()
            c.execute('DELETE FROM scripts')  # Clear previous entries
            for i in range(self.script_layout.count()):
                script_button = self.script_layout.itemAt(i).widget()
                name = script_button.text()
                path = script_button.toolTip()
                c.execute('INSERT INTO scripts (name, path) VALUES (?, ?)', (name, path))
            conn.commit()
            conn.close()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save settings: {str(e)}")

def main():
    app = QApplication(sys.argv)
    ex = MainApp()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
