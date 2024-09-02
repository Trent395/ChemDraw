import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from PIL import Image, ImageTk
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
import io
import logging  # Add this import
import molecule_functions as mf  # Import functionality from the other file

class MoleculeViewer:
    def __init__(self, master):
        self.master = master
        master.title("Organic Chemistry Molecule Viewer")
        
        self.db_manager = mf.DatabaseManager()
        self.db_manager.update_existing_table()
        
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

        self.label = tk.Label(master, text="Enter SMILES, InChI, or common name:", font=('Arial', 10))
        self.label.grid(row=1, column=0, columnspan=3, padx=10, pady=5, sticky='w')

        self.input_entry = tk.Entry(master, width=40)
        self.input_entry.grid(row=1, column=0, columnspan=3, padx=10, pady=5)

        self.draw_button = tk.Button(master, text="Draw Molecule", command=self.draw_molecule)
        self.draw_button.grid(row=2, column=0, padx=5, pady=5)

        self.clear_button = tk.Button(master, text="Clear", command=self.clear_canvas)
        self.clear_button.grid(row=2, column=1, padx=5, pady=5)

        self.save_button = tk.Button(master, text="Save Image", command=self.save_image)
        self.save_button.grid(row=2, column=2, padx=5, pady=5)

        self.generate_button = tk.Button(master, text="Generate Random Molecule", command=self.generate_molecule)
        self.generate_button.grid(row=3, column=2, padx=10, pady=5)

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

        self.canvas = tk.Canvas(master, width=400, height=400, bg='white')
        self.canvas.grid(row=4, column=1, rowspan=4, columnspan=2, padx=10, pady=5)

        self.info_label = tk.Label(master, text="", justify=tk.LEFT, font=('Arial', 10))
        self.info_label.grid(row=8, column=1, columnspan=2, padx=10, pady=5, sticky='w')

        self.longest_chain_label = tk.Label(master, text="", justify=tk.LEFT, font=('Arial', 10))
        self.longest_chain_label.grid(row=9, column=1, columnspan=2, padx=10, pady=5, sticky='w')

    def display_molecule_info(self, mol):
        # This method was missing, now added
        mol_with_hs = Chem.AddHs(mol)
        mol_weight = Descriptors.MolWt(mol_with_hs)
        atom_counts = {}
        atom_masses = {}
        for atom in mol_with_hs.GetAtoms():
            elem = atom.GetSymbol()
            atom_counts[elem] = atom_counts.get(elem, 0) + 1
            atom_masses[elem] = atom_masses.get(elem, 0) + atom.GetMass()
        total_atoms = sum(atom_counts.values())
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol_with_hs)
        info_text = f"Molecular Formula: {formula}\nMolecular Weight: {mol_weight:.2f} g/mol\nTotal Atoms: {total_atoms}"
        self.info_label.config(text=info_text)
        logging.info(f"Molecule Info - Formula: {formula}, Weight: {mol_weight}, Total Atoms: {total_atoms}")

    def update_selection(self):
        self.selected_elements = {elem: var.get() for elem, var in self.elements.items()}
        self.selected_bond_types = {bond: var.get() for bond, var in self.bond_types.items()}

    def populate_database_view(self):
        for smiles in self.db_manager.get_all_smiles():
            self.db_tree.insert("", "end", values=(smiles,))

    def generate_molecule(self):
        self.update_selection()
        smiles = mf.generate_random_smiles(self.selected_elements, self.selected_bond_types)
        if not smiles:
            return

        self.input_entry.delete(0, tk.END)
        self.input_entry.insert(0, smiles)
        self.draw_molecule()

    def draw_molecule(self):
        user_input = self.input_entry.get().strip()
        mol = None
        
        try:
            mol = Chem.MolFromSmiles(user_input)
            if mol:
                Chem.SanitizeMol(mol)
                self.db_manager.add_smiles(user_input)
                self.display_molecule_info(mol)

                pubchem_data = mf.fetch_pubchem_data(user_input)
                if pubchem_data:
                    properties = pubchem_data.get("PropertyTable", {}).get("Properties", [{}])[0]
                    mol_formula = properties.get("MolecularFormula", "N/A")
                    mol_weight = properties.get("MolecularWeight", "N/A")
                    iupac_name = properties.get("IUPACName", "N/A")

                    api_info = f"Molecular Formula: {mol_formula}\nMolecular Weight: {mol_weight} g/mol\nIUPAC Name: {iupac_name}"
                    self.info_label.config(text=api_info)

                    self.db_manager.add_molecule(user_input, mol_weight, iupac_name)

                longest_chain_name = mf.get_longest_carbon_chain(mol)
                self.longest_chain_label.config(text=f"Longest Carbon Chain: {longest_chain_name}")

                img = Draw.MolToImage(mol)
                img_bytes = io.BytesIO()
                img.save(img_bytes, format='PNG')
                img_bytes.seek(0)

                img_tk = ImageTk.PhotoImage(Image.open(img_bytes).resize((400, 400), Image.LANCZOS))
                self.canvas.create_image(200, 200, image=img_tk)
                self.canvas.image = img_tk
                self.current_image = Image.open(img_bytes)

        except Exception as e:
            logging.error(f"Error processing molecule: {e}")
            messagebox.showerror("Error", f"Could not process the molecule: {e}")

    def clear_canvas(self):
        self.canvas.delete("all")
        self.info_label.config(text="")
        self.longest_chain_label.config(text="")
        self.current_image = None

    def save_image(self):
        if self.current_image:
            file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png")])
            if file_path:
                self.current_image.save(file_path)
                messagebox.showinfo("Image Saved", "Molecule image has been saved.")
                logging.info(f"Image saved: {file_path}")
        else:
            messagebox.showwarning("No Image", "No image available to save.")

# Create the Tkinter window and run the application
root = tk.Tk()
app = MoleculeViewer(root)
root.mainloop()
