import bionetgen
import numpy as np
import os
import shutil
import pandas as pd

def run_bionetgen(n, p, p1, dir_name):
    """
    Run bionetgen simulations and process the output data.

    Parameters:
    n (int): Number of simulations to run.
    p (int): Index of the observable to extract from the output.
    dir_name (str): Directory name for the simulation files.

    Returns:
    tuple: Time and averaged observable data.
    """
    
    print(f'Running bionetgen simulations in directory: {dir_name}') 


    for i in range(n):
        result=bionetgen.run(f"{dir_name}/model.bngl", out=str(dir_name)+'/bng_output_'+str(i), suppress=True)

    pzap=[]
    psyk=[]
    time=[]
    for i in range(n):
        X=pd.read_csv(f"{dir_name}/bng_output_{i}/model.gdat", sep='\t', header=None,skiprows=1)
        X=np.asarray(X)
        time.append(X[:,0])
        pzap.append(X[:,p])
        psyk.append(X[:,p1])

    time=np.array(time).T
    pzap=np.array(pzap).T  
    psyk=np.array(psyk).T
    # Average the data
    time = np.mean(time, axis=1)
    pzap_avg = np.mean(pzap, axis=1)
    psyk_avg = np.mean(psyk, axis=1)

    # if os.path.exists(dir_name):
    #     shutil.rmtree(dir_name)
    # else:
    #     print(f"Directory {dir_name} does not exist, skipping deletion.")    

    return time, pzap_avg, psyk_avg 

