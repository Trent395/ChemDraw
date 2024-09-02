import sqlite3
from PyQt5.QtWidgets import QColorDialog
from PyQt5.QtCore import Qt

class ThemeManager:
    def __init__(self, db_path='scripts.db'):
        self.db_path = db_path
        self.create_theme_table()
        self.theme = self.load_theme()

    def create_theme_table(self):
        """Create the theme table if it doesn't exist."""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('''CREATE TABLE IF NOT EXISTS theme
                     (name TEXT PRIMARY KEY, value TEXT)''')
        conn.commit()
        conn.close()

    def open_theme_settings(self, parent=None):
        """Open color picker dialogs for background, text, and button color."""
        bg_color = QColorDialog.getColor(Qt.white, parent, "Select Background Color")
        text_color = QColorDialog.getColor(Qt.black, parent, "Select Text Color")
        button_color = QColorDialog.getColor(Qt.blue, parent, "Select Button Color")

        if bg_color.isValid() and text_color.isValid() and button_color.isValid():
            self.save_theme_to_db(bg_color.name(), text_color.name(), button_color.name())
            self.theme = self.load_theme()

    def save_theme_to_db(self, bg_color, text_color, button_color):
        """Save the selected theme colors to the database."""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('''INSERT OR REPLACE INTO theme (name, value) VALUES 
                     ('bg_color', ?), ('text_color', ?), ('button_color', ?)''', 
                     (bg_color, text_color, button_color))
        conn.commit()
        conn.close()

    def load_theme(self):
        """Load the theme from the database."""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('SELECT name, value FROM theme')
        theme = {name: value for name, value in c.fetchall()}
        conn.close()

        return {
            'bg_color': theme.get('bg_color', '#2c3e50'),
            'text_color': theme.get('text_color', '#ecf0f1'),
            'button_color': theme.get('button_color', '#3498db')
        }

    def apply_theme(self, widget):
        """Apply the loaded theme to the given widget."""
        widget.setStyleSheet(f"""
            background-color: {self.theme['bg_color']};
            color: {self.theme['text_color']};
        """)

    def get_button_style(self):
        """Return the button style based on the current theme."""
        return f"""
            QPushButton {{
                background-color: {self.theme['button_color']};
                border: none;
                padding: 10px;
                font-size: 18px;
                border-radius: 5px;
                text-align: left;
                color: {self.theme['text_color']};
            }}
            QPushButton:hover {{
                background-color: {self.theme['button_color']};
            }}
        """

