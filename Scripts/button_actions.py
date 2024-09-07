# button_actions.py

import logging
import webbrowser
from PyQt5.QtWidgets import QMessageBox, QFileDialog, QInputDialog
from database_manager import DatabaseManager
from molecule_analysis import MoleculeAnalysis  # Script that makes molecule image
from hydrogen_calculator import HydrogenCalculator
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

            self.populate_atom_bond_info(mol_analysis)

            # Hydrogen calculation
            calculated_hydrogen_count = self.h_calculator.calculate_hydrogen_count(smiles_input)
            self.viewer.hydrogen_count_label.setText(f"Calculated Hydrogen Count: {calculated_hydrogen_count}")

            # Update database and display
            self.db_manager.upsert_molecule(smiles_input, molecule_data)
            self.viewer.populate_database_view()

        except Exception as e:
            logging.error(f"Error processing molecule: {e}")
            QMessageBox.critical(self.viewer, "Error", f"Could not process the molecule: {e}")

    def populate_atom_bond_info(self, mol_analysis):
        """
        Populate atom and bond information tables based on the MoleculeAnalysis results.
        """
        atom_counts, bond_info = mol_analysis.get_atom_bond_info()
        self.viewer.atom_count_table.setRowCount(len(atom_counts))

        # Populate atom count table
        for row_idx, (element, count) in enumerate(atom_counts.items()):
            self.viewer.atom_count_table.setItem(row_idx, 0, QTableWidgetItem(element))
            self.viewer.atom_count_table.setItem(row_idx, 1, QTableWidgetItem(str(count)))

        # Populate bond info table
        self.viewer.bond_info_table.setRowCount(len(bond_info))
        for row_idx, bond in enumerate(bond_info):
            bond_type = self.get_bond_polarity(float(bond[1]))  # Classify bond type based on ΔEN
            for col_idx, value in enumerate(bond):
                self.viewer.bond_info_table.setItem(row_idx, col_idx, QTableWidgetItem(value))
            self.viewer.bond_info_table.setItem(row_idx, 3, QTableWidgetItem(bond_type))  # Add bond type

    def get_bond_polarity(self, en_diff):
        """Classify the bond type based on electronegativity difference (ΔEN)."""
        if en_diff < 0.5:
            return "Non-Polar"
        elif 0.5 <= en_diff < 1.7:
            return "Polar"
        else:
            return "Ionic"

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
            self.viewer.populate_database_view()

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

