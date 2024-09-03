from rdkit import Chem
from rdkit.Chem import Draw, inchi, Descriptors, AllChem, rdmolops, rdMolDescriptors
from PIL import Image, ImageTk
import io
import sqlite3
import logging
import random
import requests
import numpy as np
from database_manager import DatabaseManager
#from pymol import cmd, finish_launching

# Initialize PyMOL
#finish_launching()  # Start the PyMOL instance

# Set up logging
logging.basicConfig(filename='molecule_viewer.log', level=logging.INFO, format='%(asctime)s - %(message)s')


# Define the main application class
class MoleculeViewer:
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


    def on_bond_select(self, event):
        """Handle bond selection in the Treeview."""
        selected_item = self.bond_info_tree.selection()
        if selected_item:
            bond_info = self.bond_info_tree.item(selected_item)["values"]
            #print(f"Selected Bond Info: {bond_info}")

    def on_db_select(self, event):
        """Handle selection in the database Treeview."""
        selected_item = self.db_tree.selection()
        if selected_item:
            selected_smiles = self.db_tree.item(selected_item, "values")[0]
            self.input_entry.delete(0, tk.END)
            self.input_entry.insert(0, selected_smiles)
            self.draw_molecule()
            print(f"Selected molecule: {selected_smiles}")

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
                
                # Calculate and display the longest carbon chain
                
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

                                     
