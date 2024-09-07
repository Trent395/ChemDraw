from Elements import Elements  # Import the Elements class from Elements.py

class BondProperties:
    def __init__(self):
        self.elements = Elements()  # Create an instance of the Elements class
        self.molecules = {}

    def add_molecule(self, name, nodes, bonds):
        """Add a molecule with nodes (elements) and bonds (connections between nodes)."""
        self.molecules[name] = {
            "nodes": nodes,  # Dictionary of nodes with element symbols
            "bonds": bonds   # List of tuples representing bonds between nodes
        }

    def calculate_electronegativity_difference(self, atom1, atom2):
        """Calculate the electronegativity difference between two atoms."""
        en1 = self.elements.get_electronegativity(atom1)
        en2 = self.elements.get_electronegativity(atom2)
        if en1 is None or en2 is None:
            return None
        return abs(en1 - en2)

    def get_electronegativity_difference(self, molecule_name):
        """Calculate and return electronegativity differences for all bonds in a molecule."""
        molecule = self.molecules.get(molecule_name)
        if not molecule:
            return "Molecule not found"

        en_diffs = []
        for bond in molecule['bonds']:
            atom1, atom2 = bond
            en_diff = self.calculate_electronegativity_difference(molecule['nodes'][atom1], molecule['nodes'][atom2])
            en_diffs.append((atom1, atom2, en_diff))

        return en_diffs

    def calculate_angle_between_bonds(self, molecule_name, node1, node2, node3):
        """Calculate the angle between two bonds: node1-node2 and node2-node3."""
        molecule = self.molecules.get(molecule_name)
        if not molecule:
            return "Molecule not found"

        # Simplified example, assuming each bond is of unit length
        bond1_vector = (1, 0)  # Assume bond 1 lies along x-axis
        bond2_vector = (math.cos(math.radians(120)), math.sin(math.radians(120)))  # Example with 120 degree angle

        dot_product = bond1_vector[0] * bond2_vector[0] + bond1_vector[1] * bond2_vector[1]
        magnitude_bond1 = math.sqrt(bond1_vector[0] ** 2 + bond1_vector[1] ** 2)
        magnitude_bond2 = math.sqrt(bond2_vector[0] ** 2 + bond2_vector[1] ** 2)

        angle = math.acos(dot_product / (magnitude_bond1 * magnitude_bond2))
        return math.degrees(angle)

# Example usage
if __name__ == "__main__":
    bond_props = BondProperties()

    # Define water molecule (H2O)
    nodes = {
        "O": "O",  # Central Oxygen atom
        "H1": "H", # Hydrogen 1
        "H2": "H"  # Hydrogen 2
    }
    bonds = [("O", "H1"), ("O", "H2")]

    bond_props.add_molecule("H2O", nodes, bonds)

    # Calculate electronegativity differences for all bonds in H2O
    en_diffs = bond_props.get_electronegativity_difference("H2O")
    print("Electronegativity Differences in H2O:", en_diffs)

    # Calculate angle between bonds in H2O
    # Assuming H1-O-H2 bond
    angle = bond_props.calculate_angle_between_bonds("H2O", "H1", "O", "H2")
    print("Bond Angle in H2O:", angle, "degrees")
