##   PYTHON FILE HEADER #
##
##   File:         [frontend.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2024 LEEDMOL. All rights reserved.']
##   Date:         ['14.11.2024']
##   Version:      ['1.0.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['Backend of the application.']

### IMPORTS ###

from numpy import array, concatenate
from pandas import read_csv
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.Fingerprints.FingerprintMols import FingerprintMol
from rdkit.DataStructs import FingerprintSimilarity
import pickle
from warnings import filterwarnings
filterwarnings('ignore')

### CODE ###

class Backend:
    def __init__(self):
        self._loadData()
        self._loadModels()
        self.dataBaseMol = self._generateMol(self.data)

    def _loadData(self):
        self.data = read_csv('Data/no_missing_data.csv')
        self.morganBase = read_csv('Data/morgan.csv')

    def _generateMol(self, df):
        mols = []
        for i in range(len(df)):
            mol = Chem.MolFromSmiles(df['smiles'][i])
            mols.append(mol)
        return mols
    
    def _PrepareInput(self, fp, etn):
        input = concatenate((fp, [[etn]]), axis=1)

        return input

    def _tanimoto(self, queryMol, threshold=0.6):
        fps = [FingerprintMol(m) for m in self.databaseMol]
        query = FingerprintMol(Chem.MolFromSmiles(queryMol))

        similarities = [(idx, FingerprintSimilarity(query, f)) for idx, f in enumerate(fps)]
        similarities.sort(key=lambda x: x[1], reverse=True)

        max_similarity = similarities[0][1] if similarities else 0
        trustable = max_similarity >= threshold

        return similarities, max_similarity * 100, trustable
        

    def _generateMF(self, smiles, nbits=2048, raio=4):
        mol = Chem.MolFromSmiles(smiles)
        bit = {}
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, nBits=nbits, radius=raio, bitInfo=bit)
        fp = array([x for x in fp])
        fp = fp.reshape(1, -1)
        return fp, bit
    
    def _solventToEtn(self, solvent):

        translate = {'Dimethylsulfoxide (0.444)': 0.444, 'Water (1.0)': 1.0, 'Dichloromethane (0.309)': 0.309, 'Acetonitrile (0.46)': 0.46, 'Tolune (0.099)': 0.099, 'Ethanol (0.654)': 0.654, 'Chloroform (0.259)': 0.259, 'Methanol (0.762)': 0.762, 'Dimethylformamide (0.386)': 0.386, 'Cyclochexane (0.006)': 0.006, 'Hexane (0.009)': 0.009, 'Tetrahidrofurane (0.605)': 0.605, 'Methyl Cianide (0.46)': 0.46, 'Acetate (0.355)': 0.355, 'Methyl Phenil (0.099)': 0.099, 'Isopropanol (0.546)': 0.546, 'Dioxane (0.164)': 0.164, 'mXylene (0.074)': 0.074, 'Chlorobenzene (0.333)': 0.333, 'Ethyl Acetate (0.228)': 0.228, 'DiethylEter (0.117)': 0.117, 'Octanol (0.537)': 0.537, 'Nitrobenzila (0.333)': 0.333, 'Benzene (0.111)': 0.111, 'Dichloroethane (0.194)': 0.194}

        return translate[solvent]
    
    def _loadModels(self):

        path_lgbm_maxabs = 'Models/lgbm_kf_max_abs.pkl'
        path_rf_maxabs = 'Models/rf_kf_max_abs.pkl'
        path_xgb_maxabs = 'Models/xgb_kf_max_abs.pkl'

        path_lgbm_maxem = 'Models/lgbm_kf_max_em.pkl'
        path_rf_maxem = 'Models/rf_kf_max_em.pkl'
        path_xgb_maxem = 'Models/xgb_kf_max_em.pkl'

        self.lgbm_maxabs = pickle.load(open(path_lgbm_maxabs, 'rb'))
        self.rf_maxabs = pickle.load(open(path_rf_maxabs, 'rb'))
        self.xgb_maxabs = pickle.load(open(path_xgb_maxabs, 'rb'))

        self.lgbm_maxem = pickle.load(open(path_lgbm_maxem, 'rb'))
        self.rf_maxem = pickle.load(open(path_rf_maxem, 'rb'))
        self.xgb_maxem = pickle.load(open(path_xgb_maxem, 'rb'))

