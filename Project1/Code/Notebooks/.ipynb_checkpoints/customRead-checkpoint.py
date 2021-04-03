import pandas as pd
import numpy as np
from blocking import block

def bigFile(filename):
    """
    Takes a filename and returns the mean and standard error of the first column,
    in addition to the onebody density arrays and time the calculation took
    """
    big = pd.read_csv(filename, delim_whitespace = True, engine="python")
    t = big.iloc[-1][0]
    accepted = big.iloc[-2][0]
    big = big[:-2]
    
    energy = np.array(big["Energy"])
    end = 2**int(np.log2(len(energy)))
    energy = energy[:end]
   
    mean, var = block(energy)
    std = np.sqrt(var)

    return mean, std, accepted, t

def smallFile(filename):
    """
    Takes a filename for a small simulation and returns the columns as dataframes together with the execution time
    """
    small = pd.read_csv(filename, delim_whitespace = True, skiprows = 8, engine="python")
    t = small.iloc[-1][2]
    
    small = pd.read_csv(filename, delim_whitespace = True, skiprows = 8, skipfooter = 2, engine="python")
    alpha = small["Param_1"]
    E = small["E"]
    E2 = small["E^2"]
    VAR = small["VAR"]
    accepted = small["%"]
    
    return alpha, E, E2, VAR, accepted, t

def oneBodyFile(filename):
    """
    Takes a filename for a oneBody density file and returns the radius bins with their frequency
    """
    oneBody = pd.read_csv(filename, delim_whitespace = True, engine="python")
    
    r = oneBody["r"]
    density = oneBody["Density"]
    
    return r, density
