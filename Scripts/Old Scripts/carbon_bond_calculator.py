import tkinter as tk
from tkinter import ttk, messagebox
import sqlite3

class CarbonBondCalculator(tk.Tk):
    def __init__(self):
        super().__init__()

        # Set up the main window
        self.title("Carbon Bond Calculator")
        self.geometry("1000x600")

        # Initialize database
        self.initialize_database()

        # Label and Entry for number of carbon atoms (x)
        self.label_x = tk.Label(self, text="Enter the number of carbon atoms (x):")
        self.label_x.pack(pady=10)
        self.entry_x = tk.Entry(self)
        self.entry_x.pack(pady=10)

        # Button to calculate and save to database
        self.calculate_button = tk.Button(self, text="Calculate and Save to DB", command=self.calculate_and_save)
        self.calculate_button.pack(pady=10)

        # Button to view results in a sortable table
        self.view_button = tk.Button(self, text="View Results", command=self.view_results)
        self.view_button.pack(pady=10)

    def initialize_database(self):
        self.conn = sqlite3.connect("carbon_bonds.db")
        self.cursor = self.conn.cursor()
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS bonds (
                x INTEGER,
                prefix TEXT,
                bc_max_linear TEXT,
                bc_min_linear TEXT,
                delta_bc_linear TEXT,
                y_max_linear TEXT,
                y_min_linear TEXT,
                bc_max_cyclic TEXT,
                bc_min_cyclic TEXT,
                delta_bc_cyclic TEXT,
                y_max_cyclic TEXT,
                y_min_cyclic TEXT
            )
        """)
        self.conn.commit()

    def calculate_and_save(self):
        try:
            # Get the number of carbon atoms from the entry
            x = int(self.entry_x.get())

            # Generate prefix based on number of carbon atoms
            prefix = self.get_prefix(x)

            # Calculate for linear structures
            bc_min_linear = f"{x} - 1"
            if x % 2 == 0:
                bc_max_linear = f"2*{x} - 1"
                delta_bc_linear = f"{x}"
            else:
                bc_max_linear = f"2*{x} - 2"
                delta_bc_linear = f"{x} - 1"

            y_max_linear = f"4*{x} - 2*({bc_min_linear})"
            y_min_linear = f"4*{x} - 2*({bc_max_linear})"

            # Calculate for cyclic structures
            bc_min_cyclic = f"{x}"
            bc_max_cyclic = f"2*{x} - 2"
            delta_bc_cyclic = f"{bc_max_cyclic} - {bc_min_cyclic}"

            y_max_cyclic = f"4*{x} - 2*({bc_min_cyclic})"
            y_min_cyclic = f"4*{x} - 2*({bc_max_cyclic})"

            # Insert into the database
            self.cursor.execute("""
                INSERT INTO bonds (x, prefix, bc_max_linear, bc_min_linear, delta_bc_linear, y_max_linear, y_min_linear,
                bc_max_cyclic, bc_min_cyclic, delta_bc_cyclic, y_max_cyclic, y_min_cyclic)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (x, prefix, bc_max_linear, bc_min_linear, delta_bc_linear, y_max_linear, y_min_linear,
                  bc_max_cyclic, bc_min_cyclic, delta_bc_cyclic, y_max_cyclic, y_min_cyclic))
            self.conn.commit()

            messagebox.showinfo("Success", f"Calculations for x = {x} saved to the database!")

        except ValueError:
            messagebox.showerror("Input Error", "Please enter a valid integer for the number of carbon atoms.")

    def view_results(self):
        # Create a new window for viewing results
        view_window = tk.Toplevel(self)
        view_window.title("Results")
        view_window.geometry("1000x400")

        # Set up the treeview (table)
        columns = ("x", "prefix", "bc_max_linear", "bc_min_linear", "delta_bc_linear", "y_max_linear", "y_min_linear",
                   "bc_max_cyclic", "bc_min_cyclic", "delta_bc_cyclic", "y_max_cyclic", "y_min_cyclic")
        tree = ttk.Treeview(view_window, columns=columns, show="headings")

        for col in columns:
            tree.heading(col, text=col, command=lambda _col=col: self.sort_column(tree, _col, False))
            tree.column(col, width=100)

        # Fetch data from the database and insert into the treeview
        self.cursor.execute("SELECT * FROM bonds")
        rows = self.cursor.fetchall()
        for row in rows:
            tree.insert("", tk.END, values=row)

        tree.pack(fill=tk.BOTH, expand=True)

    def sort_column(self, tree, col, reverse):
        l = [(tree.set(k, col), k) for k in tree.get_children("")]
        l.sort(reverse=reverse)

        for index, (val, k) in enumerate(l):
            tree.move(k, "", index)

        tree.heading(col, command=lambda: self.sort_column(tree, col, not reverse))

    def get_prefix(self, x):
        prefixes = [
            "meth-", "eth-", "prop-", "but-", "pent-", "hex-", "hept-", "oct-", "non-", "dec-",
            "undec-", "dodec-", "tridec-", "tetradec-", "pentadec-", "hexadec-", "heptadec-",
            "octadec-", "nonadec-", "icos-"
        ]
        if 1 <= x <= 20:
            return prefixes[x - 1]
        else:
            return ""

    def on_closing(self):
        self.conn.close()
        self.destroy()

if __name__ == "__main__":
    app = CarbonBondCalculator()
    app.protocol("WM_DELETE_WINDOW", app.on_closing)
    app.mainloop()
