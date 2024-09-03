# viewer.py

import sys
import logging
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QLineEdit, QFileDialog, QMessageBox, QTableWidget, QTableWidgetItem, QHeaderView, QStyleFactory, QTextEdit, QDialog, QCheckBox
from PyQt5.QtGui import QPixmap, QPalette, QColor
from PyQt5.QtCore import Qt
from database_manager import DatabaseManager
from molecule_analysis import MoleculeAnalysis
from pubchem_fetcher import fetch_pubchem_data

class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecule Viewer")
        self.setGeometry(100, 100, 1400, 800)

        self.db_manager = DatabaseManager()
        self.is_dark_mode = False
        self.initialize_ui()

    def initialize_ui(self):
        # Central Widget
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)

        # Layouts
        main_layout = QVBoxLayout(self.central_widget)
        top_layout = QHBoxLayout()
        main_layout.addLayout(top_layout)

        # Input field and buttons
        self.input_field = QLineEdit(self)
        self.input_field.setPlaceholderText("Enter SMILES, InChI, or common name...")
        self.input_field.returnPressed.connect(self.draw_molecule)
        top_layout.addWidget(self.input_field)

        draw_button = QPushButton("Draw Molecule", self)
        draw_button.clicked.connect(self.draw_molecule)
        top_layout.addWidget(draw_button)

        save_button = QPushButton("Save Image", self)
        save_button.clicked.connect(self.save_image)
        top_layout.addWidget(save_button)

        update_button = QPushButton("Update All Database", self)
        update_button.clicked.connect(self.update_all_database)
        top_layout.addWidget(update_button)

        open_pubchem_button = QPushButton("Open PubChem Page", self)
        open_pubchem_button.clicked.connect(self.open_pubchem_page)
        top_layout.addWidget(open_pubchem_button)

        # Molecule image display
        image_info_layout = QHBoxLayout()
        main_layout.addLayout(image_info_layout)

        self.image_label_molecule = QLabel(self)
        self.image_label_molecule.setAlignment(Qt.AlignCenter)
        image_info_layout.addWidget(self.image_label_molecule, stretch=3)

        self.info_table = QTableWidget(6, 2)  # Display table for compound info
        self.info_table.setHorizontalHeaderLabels(["Property", "Value"])
        self.info_table.verticalHeader().setVisible(False)
        self.info_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        image_info_layout.addWidget(self.info_table, stretch=1)

        # Database view
        self.db_table = QTableWidget(self)
        self.db_table.setColumnCount(6)
        self.db_table.setHorizontalHeaderLabels(["Star", "SMILES", "Common Name", "IUPAC Name", "Formula", "Atomic Mass"])
        self.db_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.db_table.setSortingEnabled(True)  # Enable sorting by column
        self.db_table.cellClicked.connect(self.on_db_select)
        main_layout.addWidget(self.db_table)
        self.populate_database_view()

        # Menu
        self.create_menu()

        # Set modern theme
        QApplication.setStyle(QStyleFactory.create("Fusion"))

    def create_menu(self):
        menubar = self.menuBar()
        settings_action = menubar.addAction("Toggle Dark Mode")
        settings_action.triggered.connect(self.toggle_dark_mode)

        log_viewer_action = menubar.addAction("View Logs")
        log_viewer_action.triggered.connect(self.open_log_viewer)

    def draw_molecule(self):
        user_input = self.input_field.text().strip()

        try:
            mol_analysis = MoleculeAnalysis(user_input)
            molecule_data = fetch_pubchem_data(user_input)
            self.db_manager.upsert_molecule(user_input, molecule_data)

            self.image_label_molecule.setPixmap(mol_analysis.get_molecule_pixmap())
            self.populate_atom_bond_info(mol_analysis)
            self.populate_info_table(molecule_data)

        except Exception as e:
            logging.error(f"Error processing molecule: {e}")
            QMessageBox.critical(self, "Error", f"Could not process the molecule: {e}")

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
        if self.image_label_molecule.pixmap():
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Molecule Image", "", "PNG files (*.png)")
            if file_path:
                self.image_label_molecule.pixmap().save(file_path)
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
            with open('molecule_viewer.log', 'r') as log_file):
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
            starred = molecule[5]
            star_checkbox = QCheckBox()
            star_checkbox.setChecked(bool(starred))
            star_checkbox.stateChanged.connect(lambda _, s=molecule[0]: self.db_manager.star_molecule(s))  # Connect to star/unstar logic
            self.db_table.setCellWidget(row_idx, 0, star_checkbox)
            for col_idx, value in enumerate(molecule[1:]):
                self.db_table.setItem(row_idx, col_idx + 1, QTableWidgetItem(str(value)))

    def on_db_select(self, row, column):
        selected_smiles = self.db_table.item(row, 1).text()  # Adjusted index due to added star column
        self.input_field.setText(selected_smiles)
        self.draw_molecule()

    def update_all_database(self):
        smiles_list = self.db_manager.get_all_smiles()

        for smiles in smiles_list:
            molecule_data = fetch_pubchem_data(smiles)
            self.db_manager.upsert_molecule(smiles, molecule_data)

        self.populate_database_view()
        QMessageBox.information(self, "Update Complete", "All database entries have been updated.")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    viewer = MoleculeViewer()
    viewer.show()
    sys.exit(app.exec_())
