import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from rdkit import Chem
from rdkit.Chem import Draw, inchi, Descriptors, AllChem, rdmolops, rdMolDescriptors
from PIL import Image, ImageTk
import io
import sqlite3
import logging
import random
import requests
import numpy as np


#from pymol import cmd, finish_launching

# Initialize PyMOL
#finish_launching()  # Start the PyMOL instance


# Set up logging
logging.basicConfig(filename='molecule_viewer.log', level=logging.INFO, format='%(asctime)s - %(message)s')

####
####Vars

color_table = {
    'H': 'White',
    'C': 'Gray',
    'N': 'Blue',
    'O': 'Red',
    'F': 'Green',
    'Cl': 'Green',
    'Br': 'Brown',
    'Default': 'Gray'
}

# Define a dictionary to map the length of the carbon chain to its prefix
carbon_chain_names = {
    1: "Methane",
    2: "Ethane",
    3: "Propane",
    4: "Butane",
    5: "Pentane",
    6: "Hexane",
    7: "Heptane",
    8: "Octane",
    9: "Nonane",
    10: "Decane",
    11: "Undecane",
    12: "Dodecane"
}

####
#### End variable definitions
def get_element_color(element):
    return color_table.get(element, color_table['Default'])

####    
#### Database manager class
class DatabaseManager:
    def __init__(self, db_name='molecules.db'):
        self.db_name = db_name
        self.connection = sqlite3.connect(self.db_name)
        self.cursor = self.connection.cursor()
        self.create_table()

    def create_table(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS molecules (
                id INTEGER PRIMARY KEY,
                smiles TEXT UNIQUE,
                atomic_mass REAL,
                iupac_name TEXT
            )
        ''')
        self.connection.commit()
        
    def update_existing_table(self):
        # Check and add missing columns if they don't exist
        try:
            self.cursor.execute("ALTER TABLE molecules ADD COLUMN atomic_mass REAL")
            logging.info("Added atomic_mass column to molecules table.")
        except sqlite3.OperationalError:
            logging.info("Column atomic_mass already exists.")

        try:
            self.cursor.execute("ALTER TABLE molecules ADD COLUMN iupac_name TEXT")
            logging.info("Added iupac_name column to molecules table.")
        except sqlite3.OperationalError:
            logging.info("Column iupac_name already exists.")

        self.connection.commit()

  
    def add_molecule(self, smiles, atomic_mass, iupac_name):
        try:
            self.cursor.execute("INSERT INTO molecules (smiles, atomic_mass, iupac_name) VALUES (?, ?, ?)", 
                                (smiles, atomic_mass, iupac_name))
            self.connection.commit()
        except sqlite3.IntegrityError:
            logging.warning(f"Molecule with SMILES {smiles} already exists in the database.")
            
    def add_smiles(self, smiles):
        try:
            self.cursor.execute("INSERT INTO molecules (smiles) VALUES (?)", (smiles,))
            self.connection.commit()
        except sqlite3.IntegrityError:
            logging.warning(f"SMILES {smiles} already exists in the database.")

    def get_all_smiles(self):
        self.cursor.execute("SELECT smiles FROM molecules")
        return [row[0] for row in self.cursor.fetchall()]

    def close(self):
        self.connection.close()

#### End database code
####
        
# Define the main application class
class MoleculeViewer:
    def get_longest_carbon_chain(self, mol):
        """
        Detect the longest carbon chain in a molecule and return its name.
        :param mol: RDKit molecule object
        :return: Name of the longest carbon chain (IUPAC naming)
        """
        max_length = 0
        longest_path = []

        # Iterate through all atoms to find paths between carbon atoms
        for atom1 in mol.GetAtoms():
            if atom1.GetSymbol() == 'C':  # Only consider carbon atoms
                for atom2 in mol.GetAtoms():
                    if atom2.GetSymbol() == 'C' and atom1.GetIdx() != atom2.GetIdx():
                        path = rdmolops.GetShortestPath(mol, atom1.GetIdx(), atom2.GetIdx())
                        if len(path) > max_length:
                            max_length = len(path)
                            longest_path = path

        # Determine the name of the longest carbon chain based on its length
        if max_length > 0 and max_length <= len(carbon_chain_names):
            chain_name = carbon_chain_names[max_length]
        else:
            chain_name = f"Chain of {max_length} carbons"

        return chain_name

# After calculating molecular info
        longest_chain_name = self.get_longest_carbon_chain(mol)
        
# Add longest chain information to the info text
        info_text = f"Molecular Formula: {formula}\nMolecular Weight: {mol_weight:.2f} g/mol\n" \
                    f"Total Atoms: {total_atoms}\nLongest Carbon Chain: {longest_chain_name}"
        
        self.info_label.config(text=info_text)
        logging.info(f"Molecule Info - Formula: {formula}, Weight: {mol_weight}, Total Atoms: {total_atoms}, " \
                     f"Longest Carbon Chain: {longest_chain_name}")


    def __init__(self, master):
        self.master = master
        master.title("Organic Chemistry Molecule Viewer")
        
        # Database manager
        self.db_manager = DatabaseManager()
        
        # Update the existing table to ensure the schema is correct
        self.db_manager.update_existing_table()
        
        # Variables for selected elements and bonds
        self.selected_elements = {}
        self.selected_bond_types = {}

        # Checkboxes for elements
        self.elements_frame = tk.Frame(master)
        self.elements_frame.grid(row=0, column=0, padx=10, pady=5, sticky='nw')
        self.elements = {'C': tk.BooleanVar(value=True),
                         'O': tk.BooleanVar(value=True),
                         'N': tk.BooleanVar(value=True),
                         'F': tk.BooleanVar(value=True),
                         'Cl': tk.BooleanVar(value=True),
                         'Br': tk.BooleanVar(value=True)}
        tk.Label(self.elements_frame, text="Include Elements:").pack(anchor='w')
        for elem in self.elements:
            tk.Checkbutton(self.elements_frame, text=elem, variable=self.elements[elem]).pack(anchor='w')

        # Checkboxes for bond types
        self.bonds_frame = tk.Frame(master)
        self.bonds_frame.grid(row=0, column=1, padx=10, pady=5, sticky='nw')
        self.bond_types = {'Single': tk.BooleanVar(value=True),
                           'Double': tk.BooleanVar(value=True),
                           'Triple': tk.BooleanVar(value=True)}
        tk.Label(self.bonds_frame, text="Include Bonds:").pack(anchor='w')
        for bond_type in self.bond_types:
            tk.Checkbutton(self.bonds_frame, text=bond_type, variable=self.bond_types[bond_type]).pack(anchor='w')

        
        # Bond coordinates
        self.bond_coordinates = []  

        

        # Label for input
        self.label = tk.Label(master, text="Enter SMILES, InChI, or common name:", font=('Arial', 10))
        self.label.grid(row=1, column=0, columnspan=3, padx=10, pady=5, sticky='w')

        # Entry widget for user input
        self.input_entry = tk.Entry(master, width=40)
        self.input_entry.grid(row=1, column=0, columnspan=3, padx=10, pady=5)

        # Buttons for drawing, clearing, and saving
        self.draw_button = tk.Button(master, text="Draw Molecule", command=self.draw_molecule)
        self.draw_button.grid(row=2, column=0, padx=5, pady=5)

        self.clear_button = tk.Button(master, text="Clear", command=self.clear_canvas)
        self.clear_button.grid(row=2, column=1, padx=5, pady=5)

        self.save_button = tk.Button(master, text="Save Image", command=self.save_image)
        self.save_button.grid(row=2, column=2, padx=5, pady=5)

        # Add molecule generator button
        self.generate_button = tk.Button(master, text="Generate Random Molecule", command=self.generate_molecule)
        self.generate_button.grid(row=3, column=2, padx=10, pady=5)

        # Treeview for Database View
        self.db_view_frame = ttk.Frame(master)
        self.db_view_frame.grid(row=4, column=0, rowspan=4, padx=10, pady=5, sticky='n')
        ttk.Label(self.db_view_frame, text="Database", font=('Arial', 12, 'bold')).grid(row=0, column=0, padx=10, pady=5)
        self.db_tree = ttk.Treeview(self.db_view_frame, columns=("SMILES",), show='headings', height=20)
        self.db_tree.heading("SMILES", text="SMILES")
        self.db_tree.grid(row=1, column=0, padx=10, pady=5)

        db_vsb = ttk.Scrollbar(self.db_view_frame, orient="vertical", command=self.db_tree.yview)
        db_vsb.grid(row=1, column=1, sticky="ns")
        self.db_tree.configure(yscrollcommand=db_vsb.set)

        self.populate_database_view()

        # Canvas to display the molecule image
        self.canvas = tk.Canvas(master, width=400, height=400, bg='white')
        self.canvas.grid(row=4, column=1, rowspan=4, columnspan=2, padx=10, pady=5)

        # Label for displaying molecule information
        self.info_label = tk.Label(master, text="", justify=tk.LEFT, font=('Arial', 10))
        self.info_label.grid(row=8, column=1, columnspan=2, padx=10, pady=5, sticky='w')

        # Frame for displaying atom counts
        self.atom_count_frame = ttk.Frame(master)
        self.atom_count_frame.grid(row=9, column=0, columnspan=3, padx=10, pady=5)

        # Treeview for Atom Count Table
        self.atom_count_tree = ttk.Treeview(self.atom_count_frame, columns=("Element", "Count", "Molar Mass", "Mass %"), show='headings', height=8)
        self.atom_count_tree.heading("Element", text="Element")
        self.atom_count_tree.heading("Count", text="Count")
        self.atom_count_tree.heading("Molar Mass", text="Molar Mass (g/mol)")
        self.atom_count_tree.heading("Mass %", text="Mass %")
        self.atom_count_tree.grid(row=5, column=0, padx=10, pady=5)

        # Treeview for Functional Groups
        self.functional_groups_frame = ttk.Frame(master)
        self.functional_groups_frame.grid(row=10, column=0, columnspan=1, padx=10, pady=5)
        ttk.Label(self.functional_groups_frame, text="Functional Groups", font=('Arial', 12, 'bold')).grid(row=0, column=0, padx=10, pady=5)
        self.functional_groups_tree = ttk.Treeview(self.functional_groups_frame, columns=("Functional Group",), show='headings', height=5)
        self.functional_groups_tree.heading("Functional Group", text="Functional Group")
        self.functional_groups_tree.grid(row=1, column=0, padx=10, pady=5)

        # Treeview for Bond Information
        self.bond_info_frame = ttk.Frame(self.master)
        self.bond_info_frame.grid(row=9, column=0, columnspan=3, padx=10, pady=5, sticky='nsew')
        ttk.Label(self.bond_info_frame, text="Bond Information", font=('Arial', 12, 'bold')).grid(row=0, column=0, padx=10, pady=5)
        self.bond_info_tree = ttk.Treeview(self.bond_info_frame, columns=("Bond", "EN Diff", "Length"), show='headings', height=8)
        self.bond_info_tree.heading("Bond", text="Bond")
        self.bond_info_tree.heading("EN Diff", text="EN Diff")
        self.bond_info_tree.heading("Length", text="Length (Å)")
        self.bond_info_tree.column("Bond", width=100)  # Set appropriate width
        self.bond_info_tree.column("EN Diff", width=100)  # Set appropriate width
        self.bond_info_tree.column("Length", width=100)  # Set appropriate width
        self.bond_info_tree.grid(row=1, column=0, padx=10, pady=5, sticky='nsew')

        # Configure scrollbar for Treeview
        vsb = ttk.Scrollbar(self.bond_info_frame, orient="vertical", command=self.bond_info_tree.yview)
        vsb.grid(row=1, column=1, sticky="ns")
        self.bond_info_tree.configure(yscrollcommand=vsb.set)


        # To store the current image for saving
        self.current_image = None

        # Bind selection event to bond info Treeview
        self.bond_info_tree.bind("<<TreeviewSelect>>", self.on_bond_select)
        
        # Bind selection event to database Treeview
        self.db_tree.bind("<<TreeviewSelect>>", self.on_db_select)

        # Additional UI elements to display fetched data
        self.api_data_label = tk.Label(master, text="", justify=tk.LEFT, font=('Arial', 10))
        self.api_data_label.grid(row=9, column=1, columnspan=2, padx=10, pady=5, sticky='w')

        # Add a new label to display the longest carbon chain
        self.longest_chain_label = tk.Label(master, text="", justify=tk.LEFT, font=('Arial', 10))
        self.longest_chain_label.grid(row=9, column=1, columnspan=2, padx=10, pady=5, sticky='w')


    def fetch_pubchem_data(self, smiles):
        """Fetch additional data from PubChem based on SMILES."""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            return data
        except requests.exceptions.RequestException as e:
            logging.error(f"Error fetching data from PubChem: {e}")
            messagebox.showerror("API Error", "Failed to retrieve data from PubChem.")
            return None
    
    def update_selection(self):
            """Update selected elements and bond types based on checkboxes."""
            self.selected_elements = {elem: var.get() for elem, var in self.elements.items()}
            self.selected_bond_types = {bond: var.get() for bond, var in self.bond_types.items()}

    def display_pubchem_data(self, smiles):
        """Fetch and display additional data from PubChem."""
        data = self.fetch_pubchem_data(smiles)
        if data:
            properties = data.get('PropertyTable', {}).get('Properties', [{}])[0]
            formula = properties.get('MolecularFormula', 'N/A')
            weight = properties.get('MolecularWeight', 'N/A')
            iupac_name = properties.get('IUPACName', 'N/A')
            api_info = f"PubChem Data:\nFormula: {formula}\nMolecular Weight: {weight}\nIUPAC Name: {iupac_name}"
            self.api_data_label.config(text=api_info)
            logging.info(f"Fetched PubChem data - {api_info}")
        else:
            self.api_data_label.config(text="No additional data available.")
            logging.info("No additional data fetched from PubChem.")

    def generate_random_smiles(self):
        """Generate a random SMILES string based on selected elements, bonds, and functional groups."""
        self.update_selection()  # Update selected elements and bonds from checkboxes
        
        selected_elements = [elem for elem, selected in self.selected_elements.items() if selected]
        selected_bond_types = [bond for bond, selected in self.selected_bond_types.items() if selected]
        
        if not selected_elements or not selected_bond_types:
            messagebox.showwarning("Selection Error", "No elements or bond types selected for random generation.")
            return ""
        
        bond_symbols = {'Single': '', 'Double (=)': '=', 'Triple (#)': '#'}
        bond_choices = [bond_symbols[bond] for bond in selected_bond_types]

        # Generate base molecule with random atoms and bonds
        smiles = ''
        molecule_length = random.randint(3, 10)
        for _ in range(molecule_length):
            atom = random.choice(selected_elements)
            bond = random.choice(bond_choices)
            smiles += bond + atom

        # Add functional groups to the generated molecule
        if random.random() > 0.5:  # 50% chance to add a functional group
            group_name, group_info = random.choice(list(functional_groups.items()))
            group_smarts = group_info['smarts']
            
            # Convert functional group SMILES into RDKit Mol and add it to the existing molecule
            func_group_mol = Chem.MolFromSmarts(group_smarts)
            if func_group_mol:
                base_mol = Chem.MolFromSmiles(smiles)
                combined_mol = rdmolops.CombineMols(base_mol, func_group_mol)
                smiles = Chem.MolToSmiles(combined_mol)
        
        # Validate the generated SMILES string
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return smiles
        else:
            return ""


    def on_bond_select(self, event):
        """Handle bond selection in the Treeview."""
        selected_item = self.bond_info_tree.selection()
        if selected_item:
            bond_info = self.bond_info_tree.item(selected_item)["values"]
            print(f"Selected Bond Info: {bond_info}")

    def on_db_select(self, event):
        """Handle selection in the database Treeview."""
        selected_item = self.db_tree.selection()
        if selected_item:
            selected_smiles = self.db_tree.item(selected_item, "values")[0]
            self.input_entry.delete(0, tk.END)
            self.input_entry.insert(0, selected_smiles)
            self.draw_molecule()

    def populate_database_view(self):
        """Populate the Treeview with SMILES from the database."""
        for smiles in self.db_manager.get_all_smiles():
            self.db_tree.insert("", "end", values=(smiles,))

    def generate_molecule(self):
        """Generate a random molecule and display it."""
        smiles = self.generate_random_smiles()
        if not smiles:
            return

        self.input_entry.delete(0, tk.END)
        self.input_entry.insert(0, smiles)
        self.draw_molecule()
        
    def draw_molecule(self):
        user_input = self.input_entry.get().strip()
        mol = None
        
        # Try to interpret input as SMILES
        try:
            mol = Chem.MolFromSmiles(user_input)
            if mol:
                Chem.SanitizeMol(mol)
                logging.info(f"SMILES input recognized: {user_input}")
        except Exception as e:
            logging.error(f"Error interpreting SMILES: {e}")

        # Try to interpret input as InChI if SMILES fails
        if not mol:
            try:
                mol = inchi.MolFromInchi(user_input)
                if mol:
                    logging.info(f"InChI input recognized: {user_input}")
            except Exception as e:
                logging.error(f"Error interpreting InChI: {e}")

        # Try to interpret input as a common name if InChI fails
        if not mol:
            try:
                mol = Chem.MolFromMolFile(Chem.MolToMolBlock(Chem.MolFromName(user_input)))
                if mol:
                    logging.info(f"Common name input recognized: {user_input}")
            except Exception as e:
                logging.error(f"Error interpreting common name: {e}")

        # If the molecule is valid, process it
        if mol:
            try:
                # Compute 2D coordinates for molecule drawing
                AllChem.Compute2DCoords(mol)

                # Save to database if valid
                self.db_manager.add_smiles(user_input)

                # Display molecule information
                self.display_molecule_info(mol)

                # Fetch and display additional data from PubChem
                pubchem_data = self.fetch_pubchem_data(user_input)
                if pubchem_data:
                    properties = pubchem_data.get("PropertyTable", {}).get("Properties", [{}])[0]
                    mol_formula = properties.get("MolecularFormula", "N/A")
                    mol_weight = properties.get("MolecularWeight", "N/A")
                    iupac_name = properties.get("IUPACName", "N/A")

                    api_info = f"Molecular Formula: {mol_formula}\n" \
                               f"Molecular Weight: {mol_weight} g/mol\n" \
                               f"IUPAC Name: {iupac_name}"
                    self.api_data_label.config(text=api_info)

                    # Save molecule data to the database
                    self.db_manager.add_molecule(user_input, mol_weight, iupac_name)

                # Calculate and display the longest carbon chain
                longest_chain_name = self.get_longest_carbon_chain(mol)
                self.longest_chain_label.config(text=f"Longest Carbon Chain: {longest_chain_name}")

                # Convert molecule to image and display it in the GUI
                img = Draw.MolToImage(mol)
                img_bytes = io.BytesIO()
                img.save(img_bytes, format='PNG')
                img_bytes.seek(0)

                img_tk = ImageTk.PhotoImage(Image.open(img_bytes).resize((400, 400), Image.LANCZOS))
                self.canvas.create_image(200, 200, image=img_tk)
                self.canvas.image = img_tk  # Keep a reference to avoid garbage collection
                self.current_image = Image.open(img_bytes)

            except Exception as e:
                logging.error(f"Error processing molecule: {e}")
                messagebox.showerror("Error", f"Could not process the molecule: {e}")

        else:
            # If no valid molecule was found
            messagebox.showerror("Error", "Could not interpret input as a valid molecule. Please enter a valid SMILES, InChI, or common name.")
            logging.warning(f"Invalid input: {user_input}")



    def display_molecule_info(self, mol):
        # Add explicit hydrogens to the molecule
        mol_with_hs = Chem.AddHs(mol)

        # Calculate molecular weight and number of atoms
        mol_weight = Descriptors.MolWt(mol_with_hs)
        atom_counts = {}
        atom_masses = {}
        electronegativity = {
            'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
            'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'I': 2.66
        }
        for atom in mol_with_hs.GetAtoms():
            elem = atom.GetSymbol()
            atom_counts[elem] = atom_counts.get(elem, 0) + 1
            atom_masses[elem] = atom_masses.get(elem, 0) + atom.GetMass()
        total_atoms = sum(atom_counts.values())

        # Molecular formula
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol_with_hs)

        # Calculate molar mass contributions and percentages
        molar_mass_contributions = {elem: mass for elem, mass in atom_masses.items()}
        molar_mass_percentages = {
            elem: (mass / mol_weight) * 100 for elem, mass in molar_mass_contributions.items()
        }

        # Detect functional groups
        functional_groups = []
        smarts_patterns = {
            'Alcohol': '[CX4][OH]',
            'Amine': '[NX3][CX4]',
            'Amide': '[CX3](=[OX1])[NX3]',
            'Ester': '[CX3](=[OX1])[OX2H0]',
            'Carboxylic Acid': '[CX3](=O)[OX2H1]',
            'Aromatic': 'c1ccccc1',
            'Alkene': 'C=C',
            'Alkyne': 'C#C'
        }
        for group, smarts in smarts_patterns.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                functional_groups.append(group)

        # Populate Atom Count Table
        for i in self.atom_count_tree.get_children():
            self.atom_count_tree.delete(i)  # Clear previous entries

        for elem, count in atom_counts.items():
            color = get_element_color(elem)
            self.atom_count_tree.insert("", "end", values=(elem, count, f"{molar_mass_contributions[elem]:.2f}", f"{molar_mass_percentages[elem]:.2f}%"), tags=(elem,))

        # Configure tags for Treeview based on color table
        for elem, color_name in color_table.items():
            self.atom_count_tree.tag_configure(elem, foreground=color_name.lower())

        # Populate Functional Groups Table
        for i in self.functional_groups_tree.get_children():
            self.functional_groups_tree.delete(i)  # Clear previous entries

        for group in functional_groups:
            self.functional_groups_tree.insert("", "end", values=(group,))

        # Populate Bond Information Table
        bond_info = []
        for i in self.bond_info_tree.get_children():
            self.bond_info_tree.delete(i)  # Clear previous entries

        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            elem1 = atom1.GetSymbol()
            elem2 = atom2.GetSymbol()
            if elem1 in electronegativity and elem2 in electronegativity:
                en_diff = abs(electronegativity[elem1] - electronegativity[elem2])
                bond_length = Chem.rdMolTransforms.GetBondLength(mol_with_hs.GetConformer(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                bond_info.append((f"{elem1}-{elem2}", f"{en_diff:.2f}", f"{bond_length:.2f} Å"))
        
        # Insert bond information with color tags
        for bond in bond_info:
            bond_tag = bond[0].split('-')[0]  # Get the first element of the bond
            self.bond_info_tree.insert("", "end", values=bond, tags=(bond_tag,))

        # Configure tags for Treeview based on color table
        for elem, color_name in color_table.items():
            self.bond_info_tree.tag_configure(elem, foreground=color_name)

        # Display additional molecule information
        info_text = f"Molecular Formula: {formula}\nMolecular Weight: {mol_weight:.2f} g/mol\nTotal Atoms: {total_atoms}\n"
        self.info_label.config(text=info_text)
        logging.info(f"Molecule Info - Formula: {formula}, Weight: {mol_weight}, Total Atoms: {total_atoms}, Atom Count: {atom_counts}, Groups: {functional_groups}")

    def clear_canvas(self):
        """Clear the canvas and reset displayed information."""
        self.canvas.delete("all")
        self.info_label.config(text="")
        self.current_image = None
        for widget in [self.atom_count_tree, self.functional_groups_tree, self.bond_info_tree]:
            for i in widget.get_children():
                widget.delete(i)  # Clear treeview entries

    def save_image(self):
        """Save the currently displayed molecule image."""
        if self.current_image:
            file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png")])
            if file_path:
                self.current_image.save(file_path)
                messagebox.showinfo("Image Saved", "Molecule image has been saved.")
                logging.info(f"Image saved: {file_path}")
        else:
            messagebox.showwarning("No Image", "No image available to save.")


    def view_database(self):
        # Create a new window for the database viewer
        db_window = tk.Toplevel(self.master)
        db_window.title("Database Viewer")
        db_window.geometry("600x600")  # Set the window size

        # Define the columns, including the new ones
        columns = ("SMILES", "Atomic Mass", "IUPAC Name")
        
        # Create the Treeview widget and set the columns
        db_tree = ttk.Treeview(db_window, columns=columns, show='headings', height=15)
        
        # Define the headings for each column
        db_tree.heading("SMILES", text="SMILES")
        db_tree.heading("Atomic Mass", text="Atomic Mass")
        db_tree.heading("IUPAC Name", text="IUPAC Name")
        
        # Set the column widths (you can adjust these)
        db_tree.column("SMILES", width=200)
        db_tree.column("Atomic Mass", width=100)
        db_tree.column("IUPAC Name", width=200)
        
        # Pack the Treeview into the window
        db_tree.pack(fill=tk.BOTH, expand=True)
        
        # Fetch data from the database
        for row in self.db_manager.get_all_molecules():
            db_tree.insert("", "end", values=row)
        
        # Add a scrollbar
        db_vsb = ttk.Scrollbar(db_window, orient="vertical", command=db_tree.yview)
        db_vsb.pack(side=tk.RIGHT, fill=tk.Y)
        db_tree.configure(yscrollcommand=db_vsb.set)

        # Function to handle selection events
        def on_db_select(event):
            selected_item = db_tree.selection()
            if selected_item:
                selected_values = db_tree.item(selected_item, "values")
                selected_smiles = selected_values[0]  # Get the SMILES string
                
                # Insert the SMILES string into the input field
                self.input_entry.delete(0, tk.END)
                self.input_entry.insert(0, selected_smiles)
                
                # Optionally, use atomic_mass or iupac_name if needed
                # selected_atomic_mass = selected_values[1]
                # selected_iupac_name = selected_values[2]
                
                self.draw_molecule()  # Draw the molecule based on the selected SMILES
        
        # Bind the selection event to the Treeview
        db_tree.bind("<<TreeviewSelect>>", on_db_select)
        
        # Add a close button
        close_button = tk.Button(db_window, text="Close", command=db_window.destroy)
        close_button.pack(pady=10)

# Create the Tkinter window and run the application
root = tk.Tk()
app = MoleculeViewer(root)
root.mainloop()
self.db_manager = DatabaseManager()  # This will automatically update the existing table if needed

                                     
