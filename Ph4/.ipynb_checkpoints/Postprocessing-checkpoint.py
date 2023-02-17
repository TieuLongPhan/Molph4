import numpy as np
import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
import os
from sklearn.preprocessing import MinMaxScaler, StandardScaler, RobustScaler
class pharmacophore_postprocess:
    
    """
    Post processing for pharmacophore searching database, using to support MOE application.
    

    Parameters
    ----------
    data : pandas.DataFrame
        Data searched Pharmacophore with "ID" and "scores" column.
    model : str
        Name of pharmacophore model.
    ref : pandas.DataFrame
        Data before Pharmacophore searching with "ID" column.
    ID : str
        Identification of Molecules.
    scores: float
        Docking score or RMSD columns in pharamacophore searching.
    rescore: str
        Rescore method: minmax, standard, quantile.

    Returns
    -------
    df: pandas.DataFrame
        Data after post processing with Active, predict and rescore columns.
    """
    
    def __init__(self,data, model, ref, ID, scores, rescore = None):
        self.data = data
        self.model = model
        self.ref = ref
        self.ID = ID
        self.scores = scores # score in model (docking score in docking and rmsd in pharmacophore)
        self.rescore = rescore
        
        
    def post_process(self):
        # Add column for predict and rmsd
        self.ref["predict"] = np.zeros(len(self.ref))
        self.ref["scores"] = self.data[self.scores].max()

        idx = []
        for key, values in enumerate(self.ref[self.ID]):
            if values in self.data[self.ID].values:
                idx.append(key)
                
        # convert "predict" and "rmsd" by index
        self.ref.loc[idx, "predict"] = 1.0
        self.ref.loc[idx, "scores"] = self.data[self.scores].values
        
    def Rescore(self):
        if self.rescore == None:
            self.ref[f'{self.model}_rescore'] = -self.ref['scores']
        elif self.rescore == 'minmax':            
            scl = MinMaxScaler().fit(-self.ref['scores'].values.reshape(-1, 1))
            self.ref[f'{self.model}_rescore'] = scl.transform(-self.ref['scores'].values.reshape(-1, 1))
        elif self.rescore == 'standard':            
            scl = StandardScaler().fit(-self.ref['scores'].values.reshape(-1, 1))
            self.ref[f'{self.model}_rescore'] = scl.transform(-self.ref['scores'].values.reshape(-1, 1))
        elif self.rescore == 'quantile':            
            scl = RobustScaler().fit(-self.data['scores'].values.reshape(-1, 1))
            self.ref[f'{self.model}_rescore'] = scl.transform(-self.ref['scores'].values.reshape(-1, 1))
        self.df = self.ref
            
    def fit(self):
        self.post_process()
        self.Rescore()
