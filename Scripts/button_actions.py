# button_actions.py

import os
import logging
import webbrowser
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QInputDialog, QDialog, QGridLayout, QLabel
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt
from database_manager import DatabaseManager
from molecule_analysis import MoleculeAnalysis
from hydrogen_calculator import HydrogenCalculator
from log_viewer import LogViewer
from rdkit import Chem
from PIL import ImageQt


class ButtonActions:
    def __init__(self, viewer_instance):
        """
        Initialize ButtonActions with a reference to the MoleculeViewer instance.
        """
        self.viewer = viewer_instance
        self.db_manager = DatabaseManager()
        self.h_calculator = HydrogenCalculator()

    def draw_molecule(self, smiles_input):
        """
        Draw the molecule, display hydrogen count, and update the database for the entered SMILES string.
        """
        if not smiles_input:
            logging.error("Invalid SMILES input: None or empty string received.")
            QMessageBox.warning(self.viewer, "Invalid Input", "Please enter a valid SMILES string.")
            return

        try:
            mol_analysis = MoleculeAnalysis(smiles_input)
            molecule_data = self.db_manager.fetch_pubchem_data(smiles_input)

            if molecule_data['iupac_name'] == 'N/A':
                logging.warning(f"No valid data found for SMILES: {smiles_input}")
                QMessageBox.warning(self.viewer, "Data Not Found", "No valid data was found for the entered SMILES string.")
                return

            # Update UI elements
            self.viewer.info_label.setText(mol_analysis.get_info_text(molecule_data))
            self.viewer.image_label.setPixmap(mol_analysis.get_pixmap())

            # Call the MoleculeViewer's method to populate atom and bond information (viewer is responsible for this)
            self.viewer.table_populator.populate_atom_bond_info(mol_analysis)

            # Hydrogen calculation
            calculated_hydrogen_count = self.h_calculator.calculate_hydrogen_count(smiles_input)
            self.viewer.hydrogen_count_label.setText(f"Calculated Hydrogen Count: {calculated_hydrogen_count}")

            # Update database and display
            self.db_manager.upsert_molecule(smiles_input, molecule_data)
            self.viewer.table_populator.populate_database_view()  # Make sure to call this after updating

        except Exception as e:
            logging.error(f"Error processing molecule: {e}")
            QMessageBox.critical(self.viewer, "Error", f"Could not process the molecule: {e}")

    def save_image(self):
        """Save the displayed molecule image to a file."""
        if self.viewer.image_label.pixmap():
            file_path, _ = QFileDialog.getSaveFileName(self.viewer, "Save Molecule Image", "", "PNG files (*.png)")
            if file_path:
                self.viewer.image_label.pixmap().save(file_path)
                QMessageBox.information(self.viewer, "Image Saved", "Molecule image has been saved.")
                logging.info(f"Molecule image saved: {file_path}")
        else:
            QMessageBox.warning(self.viewer, "No Image", "No image available to save.")

    def add_nickname(self):
        """Allow the user to add a nickname for the molecule."""
        selected_smiles = self.viewer.input_field.text().strip()
        if not selected_smiles or selected_smiles.lower() == 'none':
            QMessageBox.warning(self.viewer, "No SMILES Selected", "Please select or enter a valid SMILES string.")
            return

        nickname, ok = QInputDialog.getText(self.viewer, "Add Nickname", "Enter a nickname for this molecule:")
        if ok and nickname:
            molecule_data = self.db_manager.fetch_pubchem_data(selected_smiles)
            molecule_data['nickname'] = nickname

            self.db_manager.upsert_molecule(selected_smiles, molecule_data)
            self.viewer.table_populator.populate_database_view()

            QMessageBox.information(self.viewer, "Nickname Added", "The nickname has been added successfully.")
        else:
            QMessageBox.warning(self.viewer, "Operation Cancelled", "No nickname was added.")

    def open_pubchem_page(self):
        """Open PubChem page for the selected SMILES string."""
        selected_smiles = self.viewer.input_field.text().strip()
        if selected_smiles:
            url = f"https://pubchem.ncbi.nlm.nih.gov/#query={selected_smiles}&input_type=smiles"
            webbrowser.open(url)
        else:
            QMessageBox.warning(self.viewer, "No SMILES Input", "Please enter a valid SMILES string or select a molecule from the database.")

    def update_all_database(self):
        """Update the database with all entries and refresh the view."""
        self.db_manager.update_all_database()
        self.viewer.table_populator.populate_database_view()  # Refresh after updating the database
        QMessageBox.information(self.viewer, "Update Complete", "All database entries have been updated.")
        logging.info("All database entries have been updated.")

    def open_log_viewer(self):
        """Open the log viewer for displaying log entries."""
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        logs_folder = os.path.join(root_dir, 'Logs')
        log_file_name = 'molecule_viewer.log'
        log_viewer = LogViewer(logs_folder, log_file_name)
        log_viewer.exec_()

    def generate_smiles_grid(self):
        """Generate a grid of SMILES images for valid SMILES strings."""
        valid_smiles_list = self.viewer.smiles_generator.generate_valid_smiles()
        grid_window = QDialog(self.viewer)
        grid_window.setWindowTitle("SMILES Grid")
        grid_window.setGeometry(100, 100, 1000, 800)
        grid_layout = QGridLayout(grid_window)

        for idx, smiles in enumerate(valid_smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(mol)  # Generate molecule image
            qim = ImageQt(img)
            pix = QPixmap.fromImage(qim)
            label = QLabel(self.viewer)
            label.setPixmap(pix)
            row, col = divmod(idx, 5)
            grid_layout.addWidget(label, row, col)

        grid_window.setLayout(grid_layout)
        grid_window.exec_()
