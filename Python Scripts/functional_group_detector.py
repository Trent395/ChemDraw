from rdkit import Chem

class FunctionalGroupDetector:
    def __init__(self):
        self.functional_groups = {
            'Alkane (single bond)': '[CX4]',
            'Alkene (double bond)': '[CX3]=[CX3]',
            'Alkyne (triple bond)': '[CX2]#C',
            'Aromatic (benzene ring)': 'c1ccccc1',
            'Alcohol (hydroxyl group)': '[OX2H]',
            'Ether (R-O-R)': '[OD2]([#6])[#6]',
            'Aldehyde (formyl group)': '[CX3H1](=O)[#6]',
            'Ketone (carbonyl group)': '[CX3](=O)[#6]',
            'Carboxylic Acid (carboxyl group)': '[CX3](=O)[OX2H1]',
            'Ester (COOR group)': '[CX3](=O)[OX2H0][#6]',
            'Amine (primary)': '[NX3;H2,H1][#6]',
            'Amine (secondary)': '[NX3;H1][#6][#6]',
            'Amine (tertiary)': '[NX3]([#6])[#6][#6]',
            'Amide (CONH2)': '[NX3][CX3](=[OX1])[#6]',
            'Thiol (mercaptan)': '[SX2H]',
            'Halide (F, Cl, Br, I)': '[F,Cl,Br,I]',
        }
    
    def detect(self, mol):
        detected_groups = []
        for group_name, smarts in self.functional_groups.items():
            pattern = Chem.MolFromSmarts(smarts)
            if mol.HasSubstructMatch(pattern):
                detected_groups.append(group_name)
        return detected_groups
