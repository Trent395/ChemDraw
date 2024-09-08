class Elements:
    def __init__(self):
        # Element Names
        self.element_names = {
            "H": "Hydrogen", "He": "Helium", "Li": "Lithium", "Be": "Beryllium", "B": "Boron",
            "C": "Carbon", "N": "Nitrogen", "O": "Oxygen", "F": "Fluorine", "Ne": "Neon",
            "Na": "Sodium", "Mg": "Magnesium", "Al": "Aluminum", "Si": "Silicon", "P": "Phosphorus",
            "S": "Sulfur", "Cl": "Chlorine", "Ar": "Argon", "K": "Potassium", "Ca": "Calcium",
            "Sc": "Scandium", "Ti": "Titanium", "V": "Vanadium", "Cr": "Chromium", "Mn": "Manganese",
            "Fe": "Iron", "Co": "Cobalt", "Ni": "Nickel", "Cu": "Copper", "Zn": "Zinc",
            "Ga": "Gallium", "Ge": "Germanium", "As": "Arsenic", "Se": "Selenium", "Br": "Bromine",
            "Kr": "Krypton", "Rb": "Rubidium", "Sr": "Strontium", "Y": "Yttrium", "Zr": "Zirconium",
            "Nb": "Niobium", "Mo": "Molybdenum", "Tc": "Technetium", "Ru": "Ruthenium", "Rh": "Rhodium",
            "Pd": "Palladium", "Ag": "Silver", "Cd": "Cadmium", "In": "Indium", "Sn": "Tin",
            "Sb": "Antimony", "Te": "Tellurium", "I": "Iodine", "Xe": "Xenon", "Cs": "Cesium",
            "Ba": "Barium", "La": "Lanthanum", "Ce": "Cerium", "Pr": "Praseodymium", "Nd": "Neodymium",
            "Pm": "Promethium", "Sm": "Samarium", "Eu": "Europium", "Gd": "Gadolinium", "Tb": "Terbium",
            "Dy": "Dysprosium", "Ho": "Holmium", "Er": "Erbium", "Tm": "Thulium", "Yb": "Ytterbium",
            "Lu": "Lutetium", "Hf": "Hafnium", "Ta": "Tantalum", "W": "Tungsten", "Re": "Rhenium",
            "Os": "Osmium", "Ir": "Iridium", "Pt": "Platinum", "Au": "Gold", "Hg": "Mercury",
            "Tl": "Thallium", "Pb": "Lead", "Bi": "Bismuth", "Po": "Polonium", "At": "Astatine",
            "Rn": "Radon", "Fr": "Francium", "Ra": "Radium", "Ac": "Actinium", "Th": "Thorium",
            "Pa": "Protactinium", "U": "Uranium", "Np": "Neptunium", "Pu": "Plutonium", "Am": "Americium",
            "Cm": "Curium", "Bk": "Berkelium", "Cf": "Californium", "Es": "Einsteinium", "Fm": "Fermium",
            "Md": "Mendelevium", "No": "Nobelium", "Lr": "Lawrencium", "Rf": "Rutherfordium", "Db": "Dubnium",
            "Sg": "Seaborgium", "Bh": "Bohrium", "Hs": "Hassium", "Mt": "Meitnerium", "Ds": "Darmstadtium",
            "Rg": "Roentgenium", "Cn": "Copernicium", "Nh": "Nihonium", "Fl": "Flerovium", "Mc": "Moscovium",
            "Lv": "Livermorium", "Ts": "Tennessine", "Og": "Oganesson"
        }
        
        #Atomic Mumbers
        self.atomic_numbers = {
            "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
            "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
            "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
            "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
            "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
            "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
            "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
            "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
            "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
            "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
            "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109,
            "Ds": 110, "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
        }


        # Element positions in the periodic table
        self.elements = {
            "H": (1, 1), "He": (1, 18),
            "Li": (2, 1), "Be": (2, 2), "B": (2, 13), "C": (2, 14), "N": (2, 15), "O": (2, 16), "F": (2, 17), "Ne": (2, 18),
            "Na": (3, 1), "Mg": (3, 2), "Al": (3, 13), "Si": (3, 14), "P": (3, 15), "S": (3, 16), "Cl": (3, 17), "Ar": (3, 18),
            "K": (4, 1), "Ca": (4, 2), "Sc": (4, 3), "Ti": (4, 4), "V": (4, 5), "Cr": (4, 6), "Mn": (4, 7), "Fe": (4, 8),
            "Co": (4, 9), "Ni": (4, 10), "Cu": (4, 11), "Zn": (4, 12), "Ga": (4, 13), "Ge": (4, 14), "As": (4, 15),
            "Se": (4, 16), "Br": (4, 17), "Kr": (4, 18),
            "Rb": (5, 1), "Sr": (5, 2), "Y": (5, 3), "Zr": (5, 4), "Nb": (5, 5), "Mo": (5, 6), "Tc": (5, 7), "Ru": (5, 8),
            "Rh": (5, 9), "Pd": (5, 10), "Ag": (5, 11), "Cd": (5, 12), "In": (5, 13), "Sn": (5, 14), "Sb": (5, 15),
            "Te": (5, 16), "I": (5, 17), "Xe": (5, 18),
            "Cs": (6, 1), "Ba": (6, 2), "La": (9, 4), "Ce": (9, 5), "Pr": (9, 6), "Nd": (9, 7), "Pm": (9, 8),
            "Sm": (9, 9), "Eu": (9, 10), "Gd": (9, 11), "Tb": (9, 12), "Dy": (9, 13),
            "Ho": (9, 14), "Er": (9, 15), "Tm": (9, 16), "Yb": (9, 17), "Lu": (9, 18),
            "Hf": (6, 4), "Ta": (6, 5), "W": (6, 6), "Re": (6, 7), "Os": (6, 8),
            "Ir": (6, 9), "Pt": (6, 10), "Au": (6, 11), "Hg": (6, 12), "Tl": (6, 13), "Pb": (6, 14), "Bi": (6, 15),
            "Po": (6, 16), "At": (6, 17), "Rn": (6, 18),
            "Fr": (7, 1), "Ra": (7, 2), "Ac": (10, 4), "Th": (10, 5), "Pa": (10, 6), "U": (10, 7), "Np": (10, 8),
            "Pu": (10, 9), "Am": (10, 10), "Cm": (10, 11), "Bk": (10, 12), "Cf": (10, 13), "Es": (10, 14), "Fm": (10, 15),
            "Md": (10, 16), "No": (10, 17), "Lr": (10, 18),
            "Rf": (7, 4), "Db": (7, 5), "Sg": (7, 6), "Bh": (7, 7), "Hs": (7, 8),
            "Mt": (7, 9), "Ds": (7, 10), "Rg": (7, 11), "Cn": (7, 12), "Nh": (7, 13),
            "Fl": (7, 14), "Mc": (7, 15), "Lv": (7, 16), "Ts": (7, 17), "Og": (7, 18)
        }

        # Colors for different categories of elements
        self.element_colors = {
            "nonmetal": "#FFCCCB",
            "noble_gas": "#A0C4FF",
            "alkali_metal": "#FFD700",
            "alkaline_earth_metal": "#FFB6C1",
            "metalloid": "#F0E68C",
            "halogen": "#FFC0CB",
            "metal": "#D3D3D3",
            "transition_metal": "#B0E0E6",
            "lanthanide": "#FFDEAD",
            "actinide": "#CD5C5C"
        }

        # Categories for each element
        self.element_categories = {
            "H": "nonmetal", "He": "noble_gas", "Li": "alkali_metal", "Be": "alkaline_earth_metal", 
            "B": "metalloid", "C": "nonmetal", "N": "nonmetal", "O": "nonmetal", "F": "halogen", "Ne": "noble_gas",
            "Na": "alkali_metal", "Mg": "alkaline_earth_metal", "Al": "metal", "Si": "metalloid", "P": "nonmetal",
            "S": "nonmetal", "Cl": "halogen", "Ar": "noble_gas", "K": "alkali_metal", "Ca": "alkaline_earth_metal",
            "Sc": "transition_metal", "Ti": "transition_metal", "V": "transition_metal", "Cr": "transition_metal",
            "Mn": "transition_metal", "Fe": "transition_metal", "Co": "transition_metal", "Ni": "transition_metal",
            "Cu": "transition_metal", "Zn": "transition_metal", "Ga": "metal", "Ge": "metalloid", "As": "metalloid",
            "Se": "nonmetal", "Br": "halogen", "Kr": "noble_gas", "Rb": "alkali_metal", "Sr": "alkaline_earth_metal",
            "Y": "transition_metal", "Zr": "transition_metal", "Nb": "transition_metal", "Mo": "transition_metal",
            "Tc": "transition_metal", "Ru": "transition_metal", "Rh": "transition_metal", "Pd": "transition_metal",
            "Ag": "transition_metal", "Cd": "transition_metal", "In": "metal", "Sn": "metal", "Sb": "metalloid",
            "Te": "metalloid", "I": "halogen", "Xe": "noble_gas", "Cs": "alkali_metal", "Ba": "alkaline_earth_metal",
            "La": "lanthanide", "Ce": "lanthanide", "Pr": "lanthanide", "Nd": "lanthanide", "Pm": "lanthanide",
            "Sm": "lanthanide", "Eu": "lanthanide", "Gd": "lanthanide", "Tb": "lanthanide", "Dy": "lanthanide",
            "Ho": "lanthanide", "Er": "lanthanide", "Tm": "lanthanide", "Yb": "lanthanide", "Lu": "lanthanide",
            "Hf": "transition_metal", "Ta": "transition_metal", "W": "transition_metal", "Re": "transition_metal",
            "Os": "transition_metal", "Ir": "transition_metal", "Pt": "transition_metal", "Au": "transition_metal",
            "Hg": "transition_metal", "Tl": "metal", "Pb": "metal", "Bi": "metal", "Po": "metalloid", "At": "halogen",
            "Rn": "noble_gas", "Fr": "alkali_metal", "Ra": "alkaline_earth_metal", "Ac": "actinide", "Th": "actinide",
            "Pa": "actinide", "U": "actinide", "Np": "actinide", "Pu": "actinide", 
            "Am": "actinide", "Cm": "actinide", "Bk": "actinide", "Cf": "actinide", 
            "Es": "actinide", "Fm": "actinide", "Md": "actinide", "No": "actinide", 
            "Lr": "actinide", "Rf": "transition_metal", "Db": "transition_metal", 
            "Sg": "transition_metal", "Bh": "transition_metal", "Hs": "transition_metal", 
            "Mt": "transition_metal", "Ds": "transition_metal", "Rg": "transition_metal", 
            "Cn": "transition_metal", "Nh": "metal", "Fl": "metal", "Mc": "metal", 
            "Lv": "metal", "Ts": "halogen", "Og": "noble_gas"
        }
    # Atomic masses for each element
        self.atomic_masses = {
            "H": 1.008, "He": 4.0026,
            "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
            "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948,
            "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845,
            "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.63, "As": 74.922, "Se": 78.96,
            "Br": 79.904, "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95,
            "Tc": 98, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
            "Sb": 121.76, "Te": 127.6, "I": 126.9, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12,
            "Pr": 140.91, "Nd": 144.24, "Pm": 145, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.5,
            "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.04, "Lu": 174.97, "Hf": 178.49, "Ta": 180.95, "W": 183.84,
            "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2,
            "Bi": 208.98, "Th": 232.04, "Pa": 231.04, "U": 238.03, "Np": 237, "Pu": 244, "Am": 243, "Cm": 247, "Bk": 247,
            "Cf": 251, "Es": 252, "Fm": 257, "Md": 258, "No": 259, "Lr": 262, "Rf": 267, "Db": 270, "Sg": 271, "Bh": 270,
            "Hs": 277, "Mt": 278, "Ds": 281, "Rg": 282, "Cn": 285, "Nh": 286, "Fl": 289, "Mc": 290, "Lv": 293, "Ts": 294, "Og": 294
        }
        
        # Electronegativities
        self.electronegativities = {
            "H": 2.20, "He": None, "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98, "Ne": None,
            "Na": 0.93, "Mg": 1.31, "Al": 1.61, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16, "Ar": None, "K": 0.82, "Ca": 1.00,
            "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66, "Mn": 1.55, "Fe": 1.83, "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65,
            "Ga": 1.81, "Ge": 2.01, "As": 2.18, "Se": 2.55, "Br": 2.96, "Kr": 3.00, "Rb": 0.82, "Sr": 0.95, "Y": 1.22, "Zr": 1.33,
            "Nb": 1.60, "Mo": 2.16, "Tc": 1.90, "Ru": 2.20, "Rh": 2.28, "Pd": 2.20, "Ag": 1.93, "Cd": 1.69, "In": 1.78, "Sn": 1.96,
            "Sb": 2.05, "Te": 2.10, "I": 2.66, "Xe": 2.60, "Cs": 0.79, "Ba": 0.89, "La": 1.10, "Ce": 1.12, "Pr": 1.13, "Nd": 1.14,
            "Pm": 1.13, "Sm": 1.17, "Eu": None, "Gd": 1.20, "Tb": 1.22, "Dy": 1.23, "Ho": 1.24, "Er": 1.24, "Tm": 1.25, "Yb": None,
            "Lu": 1.27, "Hf": 1.30, "Ta": 1.50, "W": 2.36, "Re": 1.90, "Os": 2.20, "Ir": 2.20, "Pt": 2.28, "Au": 2.54, "Hg": 2.00,
            "Tl": 1.62, "Pb": 2.33, "Bi": 2.02, "Po": 2.00, "At": 2.20, "Rn": None, "Fr": 0.70, "Ra": 0.89, "Ac": 1.10, "Th": 1.30,
            "Pa": 1.50, "U": 1.38, "Np": 1.36, "Pu": 1.28, "Am": 1.13, "Cm": 1.28, "Bk": 1.30, "Cf": 1.30, "Es": 1.30, "Fm": 1.30,
            "Md": 1.30, "No": 1.30, "Lr": None, "Rf": None, "Db": None, "Sg": None, "Bh": None, "Hs": None, "Mt": None, "Ds": None,
            "Rg": None, "Cn": None, "Nh": None, "Fl": None, "Mc": None, "Lv": None, "Ts": None, "Og": None
        }
    
        #Oxidation States
        self.oxidation_states = {
            "H": [-1, 1], "He": [0], "Li": [1], "Be": [2], "B": [3], "C": [-4, -3, -2, -1, 0, 1, 2, 3, 4], 
            "N": [-3, -2, -1, 0, 1, 2, 3, 4, 5], "O": [-2, -1, 0, 1, 2], "F": [-1], "Ne": [0],
            "Na": [1], "Mg": [2], "Al": [3], "Si": [-4, -3, -2, -1, 0, 1, 2, 3, 4], 
            "P": [-3, -2, -1, 0, 1, 2, 3, 4, 5], "S": [-2, -1, 0, 1, 2, 4, 6], "Cl": [-1, 0, 1, 3, 4, 5, 7], "Ar": [0],
            "K": [1], "Ca": [2], "Sc": [3], "Ti": [-1, 0, 1, 2, 3, 4], "V": [-3, -2, -1, 0, 1, 2, 3, 4, 5], 
            "Cr": [-4, -2, -1, 0, 1, 2, 3, 4, 5, 6], "Mn": [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7], 
            "Fe": [-4, -2, -1, 0, 1, 2, 3, 4, 5, 6],
            "Co": [-3, -1, 0, 1, 2, 3, 4, 5], "Ni": [-2, -1, 0, 1, 2, 3, 4], "Cu": [-2, 0, 1, 2, 3, 4], "Zn": [0, 1, 2], 
            "Ga": [-5, -4, -2, 1, 2, 3], "Ge": [-4, -2, 1, 2, 3, 4], "As": [-3, -2, 1, 2, 3, 5], 
            "Se": [-2, 1, 2, 3, 4, 5, 6], "Br": [-1, 0, 1, 3, 4, 5], "Kr": [0, 2, 4, 6],
            "Rb": [1], "Sr": [2], "Y": [1, 2, 3], "Zr": [0, 1, 2, 3, 4], "Nb": [-3, -1, 1, 2, 3, 4, 5], 
            "Mo": [-4, -2, -1, 0, 1, 2, 3, 4, 5, 6], "Tc": [-3, -1, 0, 1, 2, 3, 4, 5, 6, 7], 
            "Ru": [-2, 0, 1, 2, 3, 4, 5, 6, 7, 8], "Rh": [-1, 0, 1, 2, 3, 4, 5, 6],
            "Pd": [-2, 0, 1, 2, 3, 4], "Ag": [-2, -1, 0, 1, 2, 3], "Cd": [-2, 0, 1, 2], "In": [-5, -2, 1, 2, 3], 
            "Sn": [-4, -2, 0, 2, 4], "Sb": [-3, 0, 1, 2, 3, 5], "Te": [-2, 0, 2, 4, 5, 6], "I": [-1, 0, 1, 3, 5, 7], "Xe": [0, 2, 4, 6, 8],
            "Cs": [1], "Ba": [2], "La": [1, 2, 3], "Ce": [2, 3, 4], "Pr": [3, 4, 5], 
            "Nd": [2, 3], "Pm": [3], "Sm": [2, 3], "Eu": [2, 3], "Gd": [1, 2, 3], "Tb": [1, 3, 4], 
            "Dy": [2, 3], "Ho": [2, 3], "Er": [2, 3], "Tm": [2, 3], "Yb": [2, 3], "Lu": [1, 2, 3],
            "Hf": [1, 2, 3, 4], "Ta": [-1, 0, 1, 2, 3, 4, 5], "W": [-2, 0, 1, 2, 3, 4, 5, 6], 
            "Re": [-3, -1, 0, 1, 2, 3, 4, 5, 6, 7], "Os": [-2, 0, 1, 2, 3, 4, 5, 6, 7, 8], 
            "Ir": [-3, 0, 1, 2, 3, 4, 5, 6], "Pt": [-3, 0, 1, 2, 3, 4, 5, 6], "Au": [-1, 0, 1, 2, 3, 5],
            "Hg": [0, 1, 2, 4], "Tl": [-5, -2, 1, 3], "Pb": [-4, -2, 0, 2, 4], "Bi": [-3, 0, 1, 2, 3, 5], 
            "Po": [-2, 0, 2, 4, 6], "At": [-1, 0, 1, 3, 5], "Rn": [0, 2, 6],
            "Fr": [1], "Ra": [2], "Ac": [0, 3], "Th": [0, 2, 3, 4], "Pa": [2, 3, 4, 5], 
            "U": [0, 2, 3, 4, 5, 6], "Np": [2, 3, 4, 5, 6, 7], "Pu": [3, 4, 5, 6, 7], 
            "Am": [2, 3, 4, 5, 6], "Cm": [3, 4, 5, 6], "Bk": [3, 4], "Cf": [3], 
            "Es": [2, 3], "Fm": [2, 3], "Md": [2, 3], "No": [2, 3], "Lr": [3],
            "Rf": [4], "Db": [5], "Sg": [6], "Bh": [7], "Hs": [8], "Mt": [9], "Ds": [10], 
            "Rg": [11], "Cn": [12], "Nh": [13], "Fl": [14], "Mc": [15], "Lv": [16], 
            "Ts": [17], "Og": [18]
        }
        # Full electron configuration sequence
        self.electron_sequence = [
            "1s2", "2s2", "2p6", "3s2", "3p6", "4s2", "3d10", "4p6", "5s2", 
            "4d10", "5p6", "6s2", "4f14", "5d10", "6p6", "7s2", "5f14", "6d10", "7p6"
        ]

        # Stopping points for each element
        self.stop_points = {
            "H": 1, "He": 1,
            "Li": 2, "Be": 2, "B": 3, "C": 3, "N": 3, "O": 3, "F": 3, "Ne": 3,
            "Na": 4, "Mg": 4, "Al": 5, "Si": 5, "P": 5, "S": 5, "Cl": 5, "Ar": 5,
            "K": 6, "Ca": 6, "Sc": 7, "Ti": 7, "V": 7, "Cr": 7, "Mn": 7, "Fe": 7,
            "Co": 7, "Ni": 7, "Cu": 7, "Zn": 7, "Ga": 8, "Ge": 8, "As": 8, "Se": 8,
            "Br": 8, "Kr": 8, "Rb": 9, "Sr": 9, "Y": 10, "Zr": 10, "Nb": 10, "Mo": 10,
            "Tc": 10, "Ru": 10, "Rh": 10, "Pd": 10, "Ag": 10, "Cd": 10, "In": 11,
            "Sn": 11, "Sb": 11, "Te": 11, "I": 11, "Xe": 11, "Cs": 12, "Ba": 12,
            "La": 13, "Ce": 13, "Pr": 13, "Nd": 13, "Pm": 13, "Sm": 13, "Eu": 13,
            "Gd": 13, "Tb": 13, "Dy": 13, "Ho": 13, "Er": 13, "Tm": 13, "Yb": 13,
            "Lu": 13, "Hf": 14, "Ta": 14, "W": 14, "Re": 14, "Os": 14, "Ir": 14,
            "Pt": 14, "Au": 14, "Hg": 14, "Tl": 15, "Pb": 15, "Bi": 15, "Po": 15,
            "At": 15, "Rn": 15, "Fr": 16, "Ra": 16, "Ac": 17, "Th": 17, "Pa": 17,
            "U": 17, "Np": 17, "Pu": 17, "Am": 17, "Cm": 17, "Bk": 17, "Cf": 17,
            "Es": 17, "Fm": 17, "Md": 17, "No": 17, "Lr": 17, "Rf": 18, "Db": 18,
            "Sg": 18, "Bh": 18, "Hs": 18, "Mt": 18, "Ds": 18, "Rg": 18, "Cn": 18,
            "Nh": 19, "Fl": 19, "Mc": 19, "Lv": 19, "Ts": 19, "Og": 19
        }
        
        # Shorthand electron configurations
        self.shorthand_electron_configurations = {
            "H": "1s1", "He": "1s2",
            "Li": "[He] 2s1", "Be": "[He] 2s2", "B": "[He] 2s2 2p1", "C": "[He] 2s2 2p2",
            "N": "[He] 2s2 2p3", "O": "[He] 2s2 2p4", "F": "[He] 2s2 2p5", "Ne": "[He] 2s2 2p6",
            "Na": "[Ne] 3s1", "Mg": "[Ne] 3s2", "Al": "[Ne] 3s2 3p1", "Si": "[Ne] 3s2 3p2",
            "P": "[Ne] 3s2 3p3", "S": "[Ne] 3s2 3p4", "Cl": "[Ne] 3s2 3p5", "Ar": "[Ne] 3s2 3p6",
            "K": "[Ar] 4s1", "Ca": "[Ar] 4s2", "Sc": "[Ar] 3d1 4s2", "Ti": "[Ar] 3d2 4s2",
            "V": "[Ar] 3d3 4s2", "Cr": "[Ar] 3d5 4s1", "Mn": "[Ar] 3d5 4s2", "Fe": "[Ar] 3d6 4s2",
            "Co": "[Ar] 3d7 4s2", "Ni": "[Ar] 3d8 4s2", "Cu": "[Ar] 3d10 4s1", "Zn": "[Ar] 3d10 4s2",
            "Ga": "[Ar] 3d10 4s2 4p1", "Ge": "[Ar] 3d10 4s2 4p2", "As": "[Ar] 3d10 4s2 4p3",
            "Se": "[Ar] 3d10 4s2 4p4", "Br": "[Ar] 3d10 4s2 4p5", "Kr": "[Ar] 3d10 4s2 4p6",
            "Rb": "[Kr] 5s1", "Sr": "[Kr] 5s2", "Y": "[Kr] 4d1 5s2", "Zr": "[Kr] 4d2 5s2",
            "Nb": "[Kr] 4d4 5s1", "Mo": "[Kr] 4d5 5s1", "Tc": "[Kr] 4d5 5s2", "Ru": "[Kr] 4d7 5s1",
            "Rh": "[Kr] 4d8 5s1", "Pd": "[Kr] 4d10", "Ag": "[Kr] 4d10 5s1", "Cd": "[Kr] 4d10 5s2",
            "In": "[Kr] 4d10 5s2 5p1", "Sn": "[Kr] 4d10 5s2 5p2", "Sb": "[Kr] 4d10 5s2 5p3",
            "Te": "[Kr] 4d10 5s2 5p4", "I": "[Kr] 4d10 5s2 5p5", "Xe": "[Kr] 4d10 5s2 5p6",
            "Cs": "[Xe] 6s1", "Ba": "[Xe] 6s2", "La": "[Xe] 5d1 6s2", "Ce": "[Xe] 4f1 5d1 6s2",
            "Pr": "[Xe] 4f3 6s2", "Nd": "[Xe] 4f4 6s2", "Pm": "[Xe] 4f5 6s2", "Sm": "[Xe] 4f6 6s2",
            "Eu": "[Xe] 4f7 6s2", "Gd": "[Xe] 4f7 5d1 6s2", "Tb": "[Xe] 4f9 6s2", "Dy": "[Xe] 4f10 6s2",
            "Ho": "[Xe] 4f11 6s2", "Er": "[Xe] 4f12 6s2", "Tm": "[Xe] 4f13 6s2", "Yb": "[Xe] 4f14 6s2",
            "Lu": "[Xe] 4f14 5d1 6s2", "Hf": "[Xe] 4f14 5d2 6s2", "Ta": "[Xe] 4f14 5d3 6s2",
            "W": "[Xe] 4f14 5d4 6s2", "Re": "[Xe] 4f14 5d5 6s2", "Os": "[Xe] 4f14 5d6 6s2",
            "Ir": "[Xe] 4f14 5d7 6s2", "Pt": "[Xe] 4f14 5d9 6s1", "Au": "[Xe] 4f14 5d10 6s1",
            "Hg": "[Xe] 4f14 5d10 6s2", "Tl": "[Xe] 4f14 5d10 6s2 6p1", "Pb": "[Xe] 4f14 5d10 6s2 6p2",
            "Bi": "[Xe] 4f14 5d10 6s2 6p3", "Po": "[Xe] 4f14 5d10 6s2 6p4", "At": "[Xe] 4f14 5d10 6s2 6p5",
            "Rn": "[Xe] 4f14 5d10 6s2 6p6", "Fr": "[Rn] 7s1", "Ra": "[Rn] 7s2", "Ac": "[Rn] 6d1 7s2",
            "Th": "[Rn] 6d2 7s2", "Pa": "[Rn] 5f2 6d1 7s2", "U": "[Rn] 5f3 6d1 7s2",
            "Np": "[Rn] 5f4 6d1 7s2", "Pu": "[Rn] 5f6 7s2", "Am": "[Rn] 5f7 7s2", "Cm": "[Rn] 5f7 6d1 7s2",
            "Bk": "[Rn] 5f9 7s2", "Cf": "[Rn] 5f10 7s2", "Es": "[Rn] 5f11 7s2", "Fm": "[Rn] 5f12 7s2",
            "Md": "[Rn] 5f13 7s2", "No": "[Rn] 5f14 7s2", "Lr": "[Rn] 5f14 7s2 7p1",
            "Rf": "[Rn] 5f14 6d2 7s2", "Db": "[Rn] 5f14 6d3 7s2", "Sg": "[Rn] 5f14 6d4 7s2",
            "Bh": "[Rn] 5f14 6d5 7s2", "Hs": "[Rn] 5f14 6d6 7s2", "Mt": "[Rn] 5f14 6d7 7s2",
            "Ds": "[Rn] 5f14 6d9 7s1", "Rg": "[Rn] 5f14 6d10 7s1", "Cn": "[Rn] 5f14 6d10 7s2",
            "Nh": "[Rn] 5f14 6d10 7s2 7p1", "Fl": "[Rn] 5f14 6d10 7s2 7p2", "Mc": "[Rn] 5f14 6d10 7s2 7p3",
            "Lv": "[Rn] 5f14 6d10 7s2 7p4", "Ts": "[Rn] 5f14 6d10 7s2 7p5", "Og": "[Rn] 5f14 6d10 7s2 7p6"
        }

        # Noble gas full configurations
        self.noble_gas_full_configs = {
            "He": "1s2",
            "Ne": "1s2 2s2 2p6",
            "Ar": "1s2 2s2 2p6 3s2 3p6",
            "Kr": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",
            "Xe": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5p6",
            "Rn": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5p6 6s2 4f14 5d10 6p6"
        }
    # Densities of elements at STP (in g/cm³)
        self.densities = {
            "H": 0.00008988, "He": 0.0001785, "Li": 0.534, "Be": 1.85, "B": 2.34,
            "C": 2.267, "N": 0.0012506, "O": 0.001429, "F": 0.001696, "Ne": 0.0008999,
            "Na": 0.971, "Mg": 1.738, "Al": 2.70, "Si": 2.329, "P": 1.82,
            "S": 2.067, "Cl": 0.003214, "Ar": 0.0017837, "K": 0.862, "Ca": 1.54,
            "Sc": 2.989, "Ti": 4.54, "V": 6.11, "Cr": 7.19, "Mn": 7.3,
            "Fe": 7.874, "Co": 8.9, "Ni": 8.908, "Cu": 8.96, "Zn": 7.14,
            "Ga": 5.91, "Ge": 5.323, "As": 5.776, "Se": 4.809, "Br": 3.122,
            "Kr": 0.003733, "Rb": 1.532, "Sr": 2.64, "Y": 4.469, "Zr": 6.52,
            "Nb": 8.57, "Mo": 10.28, "Tc": 11, "Ru": 12.37, "Rh": 12.41,
            "Pd": 12.02, "Ag": 10.49, "Cd": 8.65, "In": 7.31, "Sn": 7.287,
            "Sb": 6.685, "Te": 6.232, "I": 4.933, "Xe": 0.005887, "Cs": 1.873,
            "Ba": 3.62, "La": 6.145, "Ce": 6.77, "Pr": 6.773, "Nd": 7.007,
            "Pm": 7.26, "Sm": 7.52, "Eu": 5.243, "Gd": 7.895, "Tb": 8.23,
            "Dy": 8.55, "Ho": 8.795, "Er": 9.066, "Tm": 9.321, "Yb": 6.965,
            "Lu": 9.84, "Hf": 13.31, "Ta": 16.69, "W": 19.25, "Re": 21.02,
            "Os": 22.59, "Ir": 22.56, "Pt": 21.45, "Au": 19.32, "Hg": 13.534,
            "Tl": 11.85, "Pb": 11.34, "Bi": 9.78, "Po": 9.32, "At": None,
            "Rn": 0.00973, "Fr": None, "Ra": 5.5, "Ac": 10.07, "Th": 11.72,
            "Pa": 15.37, "U": 18.95, "Np": 20.25, "Pu": 19.84, "Am": 13.67,
            "Cm": 13.51, "Bk": 14, "Cf": 15.1, "Es": None, "Fm": None,
            "Md": None, "No": None, "Lr": None, "Rf": None, "Db": None,
            "Sg": None, "Bh": None, "Hs": None, "Mt": None, "Ds": None,
            "Rg": None, "Cn": None, "Nh": None, "Fl": None, "Mc": None,
            "Lv": None, "Ts": None, "Og": None
        }

    def get_density(self, symbol):
        """Get the density of an element."""
        return self.densities.get(symbol, "Unknown")
    
    # Method to get atomic mass
    def get_atomic_mass(self, symbol):
        """Get the atomic mass of an element."""
        return self.atomic_masses.get(symbol, "Unknown")
    
    def get_atomic_number(self, symbol):
        """Get the atomic number of an element."""
        return self.atomic_numbers.get(symbol, "Unknown")
    
    def get_oxidation_states(self, symbol):
        """Get the oxidation states of an element."""
        return self.oxidation_states.get(symbol, "Unknown")

    def get_electronegativity(self, symbol):
        """Get the electronegativity of an element."""
        return self.electronegativities.get(symbol, "Unknown")

        
    def get_element_position(self, symbol):
        """Get the position of an element in the periodic table."""
        return self.elements.get(symbol, None)

    def get_element_color(self, symbol):
        """Get the color associated with an element's category."""
        category = self.element_categories.get(symbol, "metal")
        return self.element_colors.get(category, "#D3D3D3")  # Default to gray for metals

    def get_element_category(self, symbol):
        """Get the category of an element."""
        return self.element_categories.get(symbol, "metal")
    
    def get_electron_configuration(self, symbol):
        """Get the electron configuration of an element."""
        return self.get_full_electron_configuration(symbol)
    # Test

    def print_properties(self, symbol):
        """Print properties of a given element by symbol."""
        print(f"Name: {self.element_names.get(symbol, 'Unknown')}")
        print(f"Atomic Number: {self.atomic_numbers.get(symbol, 'Unknown')}")
        print(f"Atomic Mass: {self.get_atomic_mass(symbol)} amu")
        print(f"Electronegativity: {self.get_electronegativity(symbol)}")
        print(f"Oxidation States: {self.get_oxidation_states(symbol)}")
        print(f"Electron Configuration: {self.get_full_electron_configuration(symbol)}")
        print(f"Category: {self.get_element_category(symbol)}")
        print(f"Density: {self.get_density(symbol)} g/cm³")  # Assuming density data is available
        #print(f"Position in Periodic Table: {self.get_element_position(symbol)}")


    # Method to get the full electron configuration
    def get_full_electron_configuration(self, element_symbol):
        shorthand_config = self.shorthand_electron_configurations.get(element_symbol)
        if shorthand_config is None:
            return "Unknown element"
        
        # Extract the noble gas symbol from the shorthand config
        if shorthand_config.startswith("["):
            noble_gas = shorthand_config[1:shorthand_config.index("]")]
            return self.noble_gas_full_configs[noble_gas] + shorthand_config[shorthand_config.index("]") + 1:]
        else:
            return shorthand_config
# Create an instance of Elements
elements = Elements()

# Print properties of Gold
#elements.print_properties("Au")
