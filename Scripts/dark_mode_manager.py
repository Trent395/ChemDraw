# dark_mode_manager.py
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QApplication, QHBoxLayout, QLabel, QPushButton, QLineEdit, QFileDialog, QMessageBox, QInputDialog, QTableWidget, QTableWidgetItem, QHeaderView, QStyleFactory, QDialog, QGridLayout, QTextEdit
from PyQt5.QtGui import QPixmap, QPalette, QColor
from PyQt5.QtCore import Qt

class DarkModeManager:
    def __init__(self, is_dark_mode=True):
        """
        Initialize the DarkModeManager with a default mode (dark mode by default).
        """
        self.is_dark_mode = is_dark_mode

    def toggle_dark_mode(self):
        """
        Toggle between dark mode and light mode.
        """
        self.is_dark_mode = not self.is_dark_mode
        palette = self.get_dark_mode_palette() if self.is_dark_mode else QApplication.style().standardPalette()
        QApplication.setPalette(palette)

    def get_dark_mode_palette(self):
        """
        Return the color palette for dark mode.
        """
        palette = QPalette()
        palette.setColor(QPalette.Window, QColor(53, 53, 53))
        palette.setColor(QPalette.WindowText, Qt.white)
        palette.setColor(QPalette.Base, QColor(25, 25, 25))
        palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
        palette.setColor(QPalette.ToolTipBase, Qt.white)
        palette.setColor(QPalette.ToolTipText, Qt.white)
        palette.setColor(QPalette.Text, Qt.white)
        palette.setColor(QPalette.Button, QColor(53, 53, 53))
        palette.setColor(QPalette.ButtonText, Qt.white)
        palette.setColor(QPalette.BrightText, Qt.red)
        palette.setColor(QPalette.Link, QColor(42, 130, 218))
        palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        palette.setColor(QPalette.HighlightedText, Qt.black)
        return palette

    def apply_palette(self):
        """
        Apply the appropriate palette based on the current mode.
        """
        palette = self.get_dark_mode_palette() if self.is_dark_mode else QApplication.style().standardPalette()
        QApplication.setPalette(palette)
