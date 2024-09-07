# viewer.py

# Modules
import os
import sys
import logging
import webbrowser
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QApplication, QHBoxLayout, QLabel, QPushButton, QLineEdit, QFileDialog, QMessageBox, QInputDialog, QTableWidget, QTableWidgetItem, QHeaderView, QStyleFactory, QDialog, QGridLayout, QTextEdit
from PyQt5.QtGui import QPixmap, QPalette, QColor
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import ImageQt, Image
# My Scripts
from hydrogen_calculator import HydrogenCalculator  # Organic Chem Math
from database_manager import DatabaseManager
from molecule_analysis import MoleculeAnalysis  # Script that makes molecule image
from smiles_generator import SMILESGenerator
from log_viewer import LogViewer
from settings_manager import SettingsManager
from button_actions import ButtonActions
from table_populator import TablePopulator
from dark_mode_manager import DarkModeManager

# Logging setup (placed at the top)
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # Project root directory
logs_folder = os.path.join(root_dir, 'Logs')  # Logs folder path

# Ensure Logs folder exists
if not os.path.exists(logs_folder):
    os.makedirs(logs_folder)

# Set up logging configuration
log_file = os.path.join(logs_folder, 'molecule_viewer.log')
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,  # You can set this to DEBUG for more detailed logs
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logging.info("Logging setup complete")

class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()

        # Initialize the Database Manager and Hydrogen Calculator and table populator ad darkmode manager
        self.db_manager = DatabaseManager()
        self.h_calculator = HydrogenCalculator()
        self.table_populator = TablePopulator(self)
        self.dark_mode_manager = DarkModeManager(self)
        self.setWindowTitle("Molecule Viewer")

        # Initialize the settings manager
        self.settings_manager = SettingsManager()
        self.settings_manager.load_settings()  # Load settings at startup

        # Initialize the button actions module
        self.button_actions = ButtonActions(self)  # Pass the viewer instance to ButtonActions

        # Set the initial dark mode based on settings
        self.is_dark_mode = self.settings_manager.get_setting('dark_mode', True)
  
        # Apply dark mode if enabled
        if self.is_dark_mode:
            QApplication.setPalette(self.dark_mode_manager.get_dark_mode_palette())

        # Set the window size from saved settings
        window_size = self.settings_manager.get_setting('window_size', '1600,900')
        width, height = map(int, window_size.split(','))  # Parse width and height as integers
        self.resize(width, height)

        # Initialize UI components
        self.initialize_ui()
        self.create_menu()  # Ensure the menu is created

        # Initialize the SMILESGenerator with example molecule "Ethane"
        elements = ['C', 'H']
        counts = [2, 6]
        self.smiles_generator = SMILESGenerator(elements, counts, "Ethane")
            
    def on_db_select(self, row, column):
        selected_smiles = self.db_table.item(row, 0).text()
        if selected_smiles and selected_smiles.lower() != 'none':
            self.input_field.setText(selected_smiles)
            
            # Call draw molecule action from ButtonActions
            self.button_actions.draw_molecule(selected_smiles)

            # Now, call populate_atom_bond_info using the table_populator instance (instead of inside ButtonActions)
            mol_analysis = MoleculeAnalysis(selected_smiles)
            self.table_populator.populate_atom_bond_info(mol_analysis)

        else:
            logging.error("Invalid SMILES input selected from the database.")
            QMessageBox.warning(self, "Invalid SMILES", "The selected entry does not contain a valid SMILES code.")

    def initialize_ui(self):
        # Central Widget
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)

        # Main Layout
        main_layout = QVBoxLayout(self.central_widget)
        top_layout = QHBoxLayout()
        info_layout = QVBoxLayout()
        
        main_layout.addLayout(top_layout)
        main_layout.addLayout(info_layout)

        # Input field for SMILES
        self.input_field = QLineEdit(self)
        self.input_field.setPlaceholderText("Enter SMILES code...")
        self.input_field.returnPressed.connect(self.button_actions.draw_molecule)
        top_layout.addWidget(self.input_field)

        # Buttons
        draw_button = QPushButton("Draw Molecule", self)
        draw_button.clicked.connect(self.button_actions.draw_molecule)
        top_layout.addWidget(draw_button)

        save_button = QPushButton("Save Image", self)
        save_button.clicked.connect(self.button_actions.save_image)
        top_layout.addWidget(save_button)

        update_button = QPushButton("Update All Database", self)
        update_button.clicked.connect(self.button_actions.update_all_database)
        top_layout.addWidget(update_button)

        open_pubchem_button = QPushButton("Open PubChem Page", self)
        open_pubchem_button.clicked.connect(self.button_actions.open_pubchem_page)
        top_layout.addWidget(open_pubchem_button)

        add_nickname_button = QPushButton("Add Nickname", self)
        add_nickname_button.clicked.connect(self.button_actions.add_nickname)
        top_layout.addWidget(add_nickname_button)

        generate_smiles_button = QPushButton("Generate SMILES Grid", self)
        generate_smiles_button.clicked.connect(self.button_actions.generate_smiles_grid)
        top_layout.addWidget(generate_smiles_button)

        # Molecule Info Label
        self.info_label = QLabel(self)
        self.info_label.setAlignment(Qt.AlignTop)
        info_layout.addWidget(self.info_label)

        # Molecule Image Display
        self.image_label = QLabel(self)
        self.image_label.setAlignment(Qt.AlignCenter)
        info_layout.addWidget(self.image_label)

        # Hydrogen Count Display
        self.hydrogen_count_label = QLabel("Calculated Hydrogen Count: ", self)
        info_layout.addWidget(self.hydrogen_count_label)

        # Atom Count Table
        self.atom_count_table = QTableWidget(self)
        self.atom_count_table.setColumnCount(2)
        self.atom_count_table.setHorizontalHeaderLabels(["Element", "Count"])
        info_layout.addWidget(self.atom_count_table)

        # Bond Information Table
        self.bond_info_table = QTableWidget(self)
        self.bond_info_table.setColumnCount(4)
        self.bond_info_table.setHorizontalHeaderLabels(["Bond", "Δ EN", "Length", "Type"])
        info_layout.addWidget(self.bond_info_table)

        # Database View
        self.db_table = QTableWidget(self)
        self.db_table.setColumnCount(6)
        self.db_table.setHorizontalHeaderLabels(["SMILES", "Common Name", "IUPAC Name", "Formula", "Atomic Mass", "Starred"])
        self.db_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.db_table.cellClicked.connect(self.on_db_select)
        main_layout.addWidget(self.db_table)

        self.table_populator.populate_database_view()

        QApplication.setStyle(QStyleFactory.create("Fusion"))

    def create_menu(self):
        """Create the menu bar and add options for dark mode and log viewing."""
        menubar = self.menuBar()
        settings_action = menubar.addAction("Settings")
        settings_action.triggered.connect(self.open_settings_page) 

    def open_settings_page(self):
        """Open the settings page to modify application settings."""
        settings_page = SettingsManager(db_name='settings.db', parent=self)
        settings_page.exec_()

    def closeEvent(self, event):
        """Override close event to save window size and dark mode setting."""
        self.settings_manager.save_setting('window_size', f'{self.size().width()},{self.size().height()}')
        self.settings_manager.save_setting('dark_mode', 'True' if self.is_dark_mode else 'False')
        self.settings_manager.close()
        event.accept()
# Inside viewer.py, make sure the table_populator object is used correctly.

if __name__ == '__main__':
    app = QApplication(sys.argv)
    viewer = MoleculeViewer()
    viewer.show()
    sys.exit(app.exec_())
