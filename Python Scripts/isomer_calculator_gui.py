import sys
import os
from PyQt5 import QtWidgets
from isomer_calculator import IsomerCalculator
from dark_mode import DarkMode
from settings import Settings

class IsomerCalculatorApp(QtWidgets.QMainWindow):
    def __init__(self):
        super(IsomerCalculatorApp, self).__init__()

        # Set window properties
        self.setWindowTitle("Isomer Calculator")
        self.setGeometry(100, 100, 500, 250)

        # Apply initial settings, such as dark mode
        self.apply_initial_settings()

        # Initialize UI components
        self.initUI()

        # Settings page
        self.settings_page = None

    def initUI(self):
        # Label for number of carbon atoms
        self.label = QtWidgets.QLabel("Enter number of carbon atoms (n):", self)
        self.label.setGeometry(50, 50, 300, 20)

        # Input field for number of carbon atoms
        self.n_input = QtWidgets.QLineEdit(self)
        self.n_input.setGeometry(50, 80, 100, 30)
        self.n_input.setPlaceholderText("n")

        # Button to calculate isomers
        self.calculate_button = QtWidgets.QPushButton("Calculate Isomers", self)
        self.calculate_button.setGeometry(170, 80, 150, 30)
        self.calculate_button.clicked.connect(self.calculate_isomers)

        # Output area for displaying the result
        self.output_area = QtWidgets.QLabel(self)
        self.output_area.setGeometry(50, 120, 300, 50)

        # Button to open settings page
        self.settings_button = QtWidgets.QPushButton("Settings", self)
        self.settings_button.setGeometry(350, 10, 80, 30)
        self.settings_button.clicked.connect(self.open_settings)

    def calculate_isomers(self):
        try:
            n = int(self.n_input.text())
            calculator = IsomerCalculator()
            consider_hydrogens = False  # This could come from the settings in the future
            num_isomers = calculator.get_alkane_isomers(n, consider_hydrogens)
            self.output_area.setText(f"Number of isomers for C{n}H{2*n+2}: {num_isomers}")
        except ValueError:
            self.output_area.setText("Please enter a valid integer.")

    def apply_initial_settings(self):
        settings_manager = Settings(self, base_dir=os.path.join(os.path.expanduser("~"), "IsomerCalculatorSettings"))
        dark_mode_enabled = settings_manager.get_setting('dark_mode') == 'True'
        if dark_mode_enabled:
            DarkMode.apply(self)

    def open_settings(self):
        if self.settings_page is None:
            self.settings_page = Settings(self, base_dir=os.path.join(os.path.expanduser("~"), "IsomerCalculatorSettings"))
            self.settings_page.exec_()
            self.settings_page = None  # Reset to allow re-opening the settings dialog

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = IsomerCalculatorApp()
    window.show()
    sys.exit(app.exec_())
