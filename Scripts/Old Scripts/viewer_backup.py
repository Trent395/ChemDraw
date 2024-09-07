# viewer.py

import sys
import logging
import webbrowser
from PyQt5.QtWidgets import QBoxLayout, QMainWindow, QWidget, QVBoxLayout, QApplication, QHBoxLayout, QLabel, QPushButton, QLineEdit, QFileDialog, QMessageBox, QInputDialog, QTableWidget, QTableWidgetItem, QHeaderView, QStyleFactory, QTextEdit, QDialog, QGridLayout  # Added QGridLayout
from PyQt5.QtGui import QPixmap, QPalette, QColor
from PyQt5.QtCore import Qt
from database_manager import DatabaseManager
from molecule_analysis import MoleculeAnalysis
from smiles_generator import SMILESGenerator  # Import your SMILESGenerator class
from rdkit import Chem  # Added Chem
from rdkit.Chem import Draw, AllChem  # Added Draw
from PIL import ImageQt, Image  # Added ImageQt
from hydrogen_calculator import HydrogenCalculator

class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.db_manager = DatabaseManager()
        self.setWindowTitle("Molecule Viewer")
        self.setGeometry(100, 100, 1600, 1000)

        self.is_dark_mode = True  # Set to True by default
        self.initialize_ui()

        self.create_menu()  # Ensure the menu is created

        # Initialize the SMILESGenerator
        elements = ['C', 'H']  # Example elements, modify as needed
        counts = [2, 6]  # Example counts, modify as needed
        self.smiles_generator = SMILESGenerator(elements, counts, "Ethane")
        
        # Set the application to start in dark mode by default
        if self.is_dark_mode:
            QApplication.setPalette(self.get_dark_mode_palette())

    def initialize_ui(self):
        # Central Widget
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)

        # Main Layout
        main_layout = QVBoxLayout(self.central_widget)
        top_layout = QHBoxLayout()
        bottom_layout = QHBoxLayout()
        info_layout = QVBoxLayout()
        
        main_layout.addLayout(top_layout)
        main_layout.addLayout(bottom_layout)
        bottom_layout.addLayout(info_layout)

        # Input field
        self.input_field = QLineEdit(self)
        self.input_field.setPlaceholderText("Enter SMILES, InChI, or common name...")
        self.input_field.returnPressed.connect(self.draw_molecule)
        top_layout.addWidget(self.input_field)

        # Draw button
        draw_button = QPushButton("Draw Molecule", self)
        draw_button.clicked.connect(self.draw_molecule)
        top_layout.addWidget(draw_button)

        # Save button
        save_button = QPushButton("Save Image", self)
        save_button.clicked.connect(self.save_image)
        top_layout.addWidget(save_button)

        # Update all database button
        update_button = QPushButton("Update All Database", self)
        update_button.clicked.connect(self.update_all_database)
        top_layout.addWidget(update_button)

        # Open PubChem page button
        open_pubchem_button = QPushButton("Open PubChem Page", self)
        open_pubchem_button.clicked.connect(self.open_pubchem_page)
        top_layout.addWidget(open_pubchem_button)

        # Add Common Name button
        add_common_name_button = QPushButton("Add Common Name", self)
        add_common_name_button.clicked.connect(self.add_common_name)
        top_layout.addWidget(add_common_name_button)

        # In the `initialize_ui` method, add the following under your existing buttons:
        generate_smiles_button = QPushButton("Generate SMILES Grid", self)
        generate_smiles_button.clicked.connect(self.generate_smiles_grid)
        top_layout.addWidget(generate_smiles_button)

        # Molecule Info Label
        self.info_label = QLabel(self)
        self.info_label.setAlignment(Qt.AlignTop)
        info_layout.addWidget(self.info_label)

        # Molecule Image Display
        self.image_label = QLabel(self)
        self.image_label.setAlignment(Qt.AlignCenter)
        info_layout.addWidget(self.image_label)

        # Atom Count Table
        self.atom_count_table = QTableWidget(self)
        self.atom_count_table.setColumnCount(2)
        self.atom_count_table.setHorizontalHeaderLabels(["Element", "Count"])
        # info_layout.addWidget(self.atom_count_table)

        # Bond Information Table
        self.bond_info_table = QTableWidget(self)
        self.bond_info_table.setColumnCount(4)
        self.bond_info_table.setHorizontalHeaderLabels(["Bond", "Electronegativity Difference", "Length", "Polarity"])
        # info_layout.addWidget(self.bond_info_table)

        # Database View
        self.db_table = QTableWidget(self)
        self.db_table.setColumnCount(6)
        self.db_table.setHorizontalHeaderLabels(["SMILES", "Common Name", "IUPAC Name", "Formula", "Atomic Mass", "Starred"])
        self.db_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.db_table.cellClicked.connect(self.on_db_select)
        main_layout.addWidget(self.db_table)

        # Initialize and populate the database view
        self.populate_database_view()

        # Set modern theme
        QApplication.setStyle(QStyleFactory.create("Fusion"))

    def generate_image_from_smiles(self, smiles):
        """
        Generate an image of the molecule from a SMILES string.
        
        :param smiles: A SMILES string representing the molecule.
        :return: An image object (PIL Image) of the molecule.
        """
        try:
            # Create the RDKit molecule object from the SMILES string
            mol = Chem.MolFromSmiles(smiles)
            
            if mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
            
            # Generate the 2D coordinates for the molecule (needed for drawing)
            Chem.AllChem.Compute2DCoords(mol)
            
            # Generate the image of the molecule
            img = Draw.MolToImage(mol)
            
            # Debug: Save the image to check if it's generated correctly
            img.save("debug_image.png")
            logging.info(f"Image generated for SMILES {smiles} and saved as debug_image.png")
            
            return img
        
        except Exception as e:
            logging.error(f"Error generating image for SMILES {smiles}: {e}")
            return None
        
    def generate_smiles_grid(self):
        valid_smiles_list = self.smiles_generator.generate_valid_smiles()
        
        # Assuming you have a method that generates images from SMILES
        for smiles in valid_smiles_list:
            img = self.generate_image_from_smiles(smiles)
            qim = ImageQt(img)

        # Create a new window or widget to display the grid
        grid_window = QDialog(self)
        grid_window.setWindowTitle("SMILES Grid")
        grid_window.setGeometry(100, 100, 1000, 800)

        grid_layout = QGridLayout(grid_window)

        for idx, smiles in enumerate(valid_smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(mol)
            
            # Convert the PIL image to a format that Qt can display
            qim = ImageQt(img)
            pix = QPixmap.fromImage(qim)
            label = QLabel(self)
            label.setPixmap(pix)

            # Add the QLabel to the grid
            row, col = divmod(idx, 5)  # Display 5 molecules per row
            grid_layout.addWidget(label, row, col)

        grid_window.setLayout(grid_layout)
        grid_window.exec_()

    def create_menu(self):
        menubar = self.menuBar()
        settings_action = menubar.addAction("Toggle Dark Mode")
        settings_action.triggered.connect(self.toggle_dark_mode)

        log_viewer_action = menubar.addAction("View Logs")
        log_viewer_action.triggered.connect(self.open_log_viewer)

    def draw_molecule(self):
        user_input = self.input_field.text().strip()
        
        # Check if the input is valid
        if not user_input or user_input.lower() == 'none':
            logging.error("Invalid SMILES input: None or empty string received.")
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid SMILES string.")
            return

        try:
            mol_analysis = MoleculeAnalysis(user_input)
            molecule_data = self.db_manager.fetch_pubchem_data(user_input)  # Use self.db_manager to call the method

            # Check if PubChem data retrieval was successful
            if molecule_data['iupac_name'] == 'N/A':
                logging.warning(f"No valid data found for SMILES: {user_input}")
                QMessageBox.warning(self, "Data Not Found", "No valid data was found for the entered SMILES string.")
                return

            self.db_manager.upsert_molecule(user_input, molecule_data)

            self.info_label.setText(mol_analysis.get_info_text(molecule_data))
            self.image_label.setPixmap(mol_analysis.get_pixmap())

            self.populate_atom_bond_info(mol_analysis)

            # Update the database after each successful entry
            self.populate_database_view()

        except Exception as e:
            logging.error(f"Error processing molecule: {e}")
            QMessageBox.critical(self, "Error", f"Could not process the molecule: {e}")

    def add_common_name(self):
        selected_smiles = self.input_field.text().strip()

        if not selected_smiles or selected_smiles.lower() == 'none':
            QMessageBox.warning(self, "No SMILES Selected", "Please select or enter a valid SMILES string.")
            return

        # Prompt the user to enter a common name
        common_name, ok = QInputDialog.getText(self, "Add Common Name", "Enter a common name for this molecule:")
        if ok and common_name:
            # Check if the user wants to star this molecule
            star_check = QMessageBox.question(self, "Star", "Would you like to star this molecule?", QMessageBox.Yes | QMessageBox.No)
            starred = 1 if star_check == QMessageBox.Yes else 0

            # Update the important SMILES database
            self.db_manager.upsert_important_smiles(selected_smiles, common_name, starred)

            # Refresh the database view
            self.populate_database_view()

            QMessageBox.information(self, "Common Name Added", "The common name has been added successfully.")
        else:
            QMessageBox.warning(self, "Operation Cancelled", "No common name was added.")

    def populate_info_table(self, molecule_data):
        """Populate the info table with the fetched data."""
        properties = ["Common Name", "IUPAC Name", "Formula", "Atomic Mass", "XLogP", "H-Bond Donors", "H-Bond Acceptors"]
        for i, prop in enumerate(properties):
            self.info_table.setItem(i, 0, QTableWidgetItem(prop))
            self.info_table.setItem(i, 1, QTableWidgetItem(str(molecule_data[prop.lower().replace(" ", "_")])))

    def populate_atom_bond_info(self, mol_analysis):
        atom_counts, bond_info = mol_analysis.get_atom_bond_info()

        self.atom_count_table.setRowCount(len(atom_counts))
        for row_idx, (element, count) in enumerate(atom_counts.items()):
            self.atom_count_table.setItem(row_idx, 0, QTableWidgetItem(element))
            self.atom_count_table.setItem(row_idx, 1, QTableWidgetItem(str(count)))

        self.bond_info_table.setRowCount(len(bond_info))
        for row_idx, bond in enumerate(bond_info):
            for col_idx, value in enumerate(bond):
                self.bond_info_table.setItem(row_idx, col_idx, QTableWidgetItem(value))

    def save_image(self):
        if self.image_label.pixmap():
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Molecule Image", "", "PNG files (*.png)")
            if file_path:
                self.image_label.pixmap().save(file_path)
                QMessageBox.information(self, "Image Saved", "Molecule image has been saved.")
                logging.info(f"Molecule image saved: {file_path}")
        else:
            QMessageBox.warning(self, "No Image", "No image available to save.")

    def open_pubchem_page(self):
        selected_smiles = self.input_field.text().strip()
        if selected_smiles:
            url = f"https://pubchem.ncbi.nlm.nih.gov/#query={selected_smiles}&input_type=smiles"
            webbrowser.open(url)
        else:
            QMessageBox.warning(self, "No SMILES Input", "Please enter a valid SMILES string or select a molecule from the database.")

    def toggle_dark_mode(self):
        self.is_dark_mode = not self.is_dark_mode
        palette = self.get_dark_mode_palette() if self.is_dark_mode else QApplication.style().standardPalette()
        QApplication.setPalette(palette)

    def get_dark_mode_palette(self):
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

    def open_log_viewer(self):
        log_viewer = QDialog(self)
        log_viewer.setWindowTitle("Log Viewer")
        log_viewer.setGeometry(100, 100, 600, 400)

        layout = QVBoxLayout()

        log_text_edit = QTextEdit(log_viewer)
        log_text_edit.setReadOnly(True)
        try:
            with open('molecule_viewer.log', 'r') as log_file:
                log_text_edit.setPlainText(log_file.read())
        except FileNotFoundError:
            log_text_edit.setPlainText("Log file not found.")

        layout.addWidget(log_text_edit)
        log_viewer.setLayout(layout)
        log_viewer.exec_()

    def populate_database_view(self):
        self.db_table.setRowCount(0)
        molecules = self.db_manager.get_all_molecules()

        for row_idx, molecule in enumerate(molecules):
            self.db_table.insertRow(row_idx)
            for col_idx, value in enumerate(molecule):
                self.db_table.setItem(row_idx, col_idx, QTableWidgetItem(str(value)))


    def on_db_select(self, row, column):
        # Adjust the index to match the SMILES column in your database view.
        selected_smiles = self.db_table.item(row, 0).text()  # Assuming the SMILES column is at index 0

        if selected_smiles and selected_smiles.lower() != 'none':
            self.input_field.setText(selected_smiles)
            self.draw_molecule()
        else:
            logging.error("Invalid SMILES input selected from the database.")
            QMessageBox.warning(self, "Invalid SMILES", "The selected entry does not contain a valid SMILES code.")

    def update_all_database(self):
        self.db_manager.update_all_database()
        self.populate_database_view()
        QMessageBox.information(self, "Update Complete", "All database entries have been updated.")

        logging.info("All database entries have been updated.")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    viewer = MoleculeViewer()
    viewer.show()
    sys.exit(app.exec_())
