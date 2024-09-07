import os
import logging
import networkx as nx
from networkx.algorithms import isomorphism
from itertools import permutations

#### DOESNT WORK


# Ensure the directories exist
log_dir = 'ChemDraw/Logs'
os.makedirs(log_dir, exist_ok=True)

# Configure logging
logging.basicConfig(filename=os.path.join(log_dir, 'isomer_calculator.log'), level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

class IsomerCalculator:
    def __init__(self):
        logging.info("Initialized IsomerCalculator.")

    def get_alkane_isomers(self, n, consider_hydrogens=False):
        """Returns the number of constitutional isomers for alkanes with n carbon atoms."""
        logging.info(f"Calculating isomers for C{n}H{2*n+2}, consider hydrogens: {consider_hydrogens}")
        
        if n < 1:
            logging.error("Number of carbon atoms must be at least 1.")
            return 0

        # Generate possible structures using graph theory
        carbon_atoms = ['C'] * n
        possible_structures = []

        for perm in permutations(carbon_atoms):
            G = nx.Graph()
            G.add_nodes_from(range(n))

            # Add edges based on permutations
            for i in range(1, n):
                G.add_edge(perm[i-1], perm[i])

            # Check if this structure is isomorphic to any existing structure
            is_unique = True
            for existing_graph in possible_structures:
                if nx.is_isomorphic(G, existing_graph):
                    is_unique = False
                    break

            if is_unique:
                possible_structures.append(G)

        num_isomers = len(possible_structures)
        logging.info(f"Number of isomers for C{n}H{2*n+2}: {num_isomers}")
        return num_isomers
