# pubchem_fetcher.py

import pubchempy as pcp
import logging
from rdkit import Chem
from rdkit.Chem import Descriptors

def fetch_pubchem_data(smiles):
    data = {'common_name': 'N/A', 'molecular_formula': 'N/A', 'atomic_mass': 'N/A', 'iupac_name': 'N/A', 'xlogp': 'N/A', 'h_bond_donors': 'N/A', 'h_bond_acceptors': 'N/A'}
    try:
        compound = pcp.get_compounds(smiles, namespace='smiles')
        if compound:
            compound = compound[0]
            data['common_name'] = compound.synonyms[0] if compound.synonyms else 'N/A'  # Fetch common name
            data['molecular_formula'] = compound.molecular_formula or 'N/A'
            data['atomic_mass'] = compound.molecular_weight or 'N/A'
            data['iupac_name'] = compound.iupac_name or 'N/A'
            data['xlogp'] = compound.xlogp or 'N/A'

        mol = Chem.MolFromSmiles(smiles)
        if mol:
            data['h_bond_donors'] = Descriptors.NumHDonors(mol)
            data['h_bond_acceptors'] = Descriptors.NumHAcceptors(mol)

    except Exception as e:
        logging.error(f"Error fetching data for SMILES {smiles}: {e}")
    return data
