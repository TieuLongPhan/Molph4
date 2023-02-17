import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
from Postprocessing import pharmacophore_postprocess
from Validation import pharmacophore_validation

class autoph4result:
    
    """
    Automated validation pharmacophore.

    Parameters
    ----------
    path : str
        Path to folder contains separate pharmacophore database folder (cover 5,6,7).
    ref : pandas.DataFrame
        Data before Pharmacophore searching with "ID" column.
    active : str
        Name of "Active" column (binary).
    scores: float
        Docking score, RMSD in pharamacophore searching columns.
    ID : str
        Identification of Molecules.
    rescore: str
        Rescore method: minmax, standard, quantile. 
    auc_thresh: float
        Threshold for ROC display.

    Returns
    -------
    Raw_result: folder
        Folder contains pharmacophore results.
    Raw_result_png: folder
        Folder contains ROC plot
        
    """
    
    def __init__(self, path, ref, scores, active, ID, rescore = None, auc_thresh = 0.5):
        self.path = path # path to folder contain result folders (cover5, cover6,...)
        self.ref = ref # active + decoy
        self.active = active
        self.scores = scores
        self.ID = ID
        self.rescore = rescore
        self.auc_thresh = auc_thresh
    
    def compare_model(self, path, data):
        os.chdir(path)
        data_dir = str(os.getcwd())
        data_name = []
        for i in sorted(glob.glob(f"{path}/*.txt")):
            data_name.append(i[len(data_dir)+1:-4])
        #print(data_name)


        fig = plt.figure(figsize = (16,8))
        table = pd.DataFrame()
        for i in data_name:
            data = pd.read_csv(f"{i}.txt")
            post = pharmacophore_postprocess(data = data, model =i, ID = 'ID', ref=self.ref, 
                                 scores = 'rmsd', rescore = 'minmax')
            post.fit()
            ph4 = pharmacophore_validation(data = post.df, active = "Active", predict = "predict",
                               scores = f"{i}_rescore", model =i, auc_thresh = self.auc_thresh)
            
            ph4.validation()
            if ph4.table['Sensitivity'].astype('float').values != 0:
                table = pd.concat([table, ph4.table], axis =0).reset_index(drop = True)
            #display(post.data.head(5))
        return table, fig
    
    def compare_model_multi(self):
        os.chdir(self.path)
        data_dir = str(os.getcwd())
        folder_name = []
        for i in sorted(glob.glob(f"{self.path}/*")):
            folder_name.append(i[len(data_dir)+1:])

        directory = 'Raw_result'
        self.path_raw = os.path.join(self.path, directory)

        if os.path.isdir(self.path_raw) == False:
            os.mkdir(self.path_raw)
        else:
            pass
        
        directory = 'Raw_result_png'
        self.path_png = os.path.join(self.path, directory)

        if os.path.isdir(self.path_png) == False:
            os.mkdir(self.path_png)
        else:
            pass

        for i in folder_name:
            print(f"{i} result")
            new_path = self.path+i
            table, fig = self.compare_model(new_path, self.ref)
            table.to_csv(self.path_raw+f'/{i}.csv') 
            fig.savefig(self.path_png+f'/{i}.png', dpi = 600, transparent = True)
            
        self.concat()
            
    def concat(self):
        os.chdir(self.path_raw)
        self.data_all = pd.DataFrame()
        data_dir = str(os.getcwd())

        file_name = []
        for i in sorted(glob.glob(f"{self.path_raw}/*")):
            file_name.append(i[len(data_dir)+1:])

        for i in file_name:
            df = pd.read_csv(self.path_raw + f"/{i}").drop(['Unnamed: 0'], axis = 1)
            self.data_all = pd.concat([self.data_all, df], axis = 0)
            
        self.data_all.reset_index(drop = True, inplace = True)
        self.data_all.to_csv(self.path_raw +f'/Ph4_results.csv')
