import logging
from PyQt5.QtWidgets import QTableWidgetItem
from molecule_analysis import MoleculeAnalysis  # Import MoleculeAnalysis to use its get_bond_polarity method

class TablePopulator:
    def __init__(self, viewer):
        """
        Initialize the TablePopulator with a reference to the MoleculeViewer instance.
        """
        self.viewer = viewer

    def populate_atom_bond_info(self, mol_analysis):
        """
        Populate the atom and bond information tables with data from the molecule analysis.
        Use the MoleculeAnalysis object to reference bond polarity and other molecular data like molar mass, electronegativity, and bond length.
        """
        atom_counts, bond_info = mol_analysis.get_atom_bond_info()

        # Log atom counts
        logging.info(f"Populating atom counts: {atom_counts}")

        # Populate atom count table with hydrogen, molar mass, and electronegativity
        self.viewer.atom_count_table.setRowCount(len(atom_counts))
        row_idx = 0
        for element, data in atom_counts.items():
            self.viewer.atom_count_table.setItem(row_idx, 0, QTableWidgetItem(element))
            self.viewer.atom_count_table.setItem(row_idx, 1, QTableWidgetItem(str(data['count'])))
            self.viewer.atom_count_table.setItem(row_idx, 2, QTableWidgetItem(str(data['atomic_mass'])))
            self.viewer.atom_count_table.setItem(row_idx, 3, QTableWidgetItem(str(data['electronegativity'])))
            row_idx += 1
            #add logging
            
        # Populate bond info table and classify bonds using MoleculeAnalysis.get_bond_polarity()
        self.viewer.bond_info_table.setRowCount(len(bond_info))
        for row_idx, bond in enumerate(bond_info):
            bond_type = mol_analysis.get_bond_polarity(float(bond[1]))  # Get bond type from molecule analysis
            logging.info(f"Bond: {bond[0]}, EN Difference: {bond[1]}, Bond Length: {bond[2]}, Bond Type: {bond_type}")
            for col_idx, value in enumerate(bond):
                self.viewer.bond_info_table.setItem(row_idx, col_idx, QTableWidgetItem(value))
            self.viewer.bond_info_table.setItem(row_idx, 3, QTableWidgetItem(bond_type))  # Add bond type to the last column

    def populate_database_view(self):
        """
        Populate the database view with all molecule entries from the database.
        """
        self.viewer.db_table.setRowCount(0)
        molecules = self.viewer.db_manager.get_all_molecules()

        logging.info(f"Populating database with {len(molecules)} molecules.")
        
        for row_idx, molecule in enumerate(molecules):
            self.viewer.db_table.insertRow(row_idx)
            for col_idx, value in enumerate(molecule):
                self.viewer.db_table.setItem(row_idx, col_idx, QTableWidgetItem(str(value)))
