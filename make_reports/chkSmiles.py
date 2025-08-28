from rdkit import Chem
import sys

smi=Chem.MolFromSmiles(sys.argv[1])

if smi==None:
    print('None')
else:
    print('true')
