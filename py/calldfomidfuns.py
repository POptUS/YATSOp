import numpy as np
import time
import sys
import pickle 
from calfun_sample import calfun_sample
from ecnoise import ecnoise
from dfoxsnew import dfoxsnew

# Define the variable matrix (Var)
Var = np.array([
    [1, 100, 200, 1e-7],
    [2, 100, 200, 5e-8],
    [3, 100, 200, 5e-8],
    [201, 100, 100, 4e-9],
    [202, 100, 198, 1e-8],
    [19, 100, 192, 9e-9],
    [128, 100, 100, 1e-6],
    [203, 100, 100, 4e-10],
    [204, 100, 100, 5e-10],
    [205, 125, 125, 1e-9],
    [16, 100, 100, 4e-8],
    [130, 100, 100, 4e-9],
    [131, 100, 100, 4e-9],
    [206, 98, 98, 4e-10],
    [207, 100, 100, 6e-9],
    [15, 100, 100, 4e-10],
    [216, 100, 100, 5e-9],
    [120, 100, 100, 4e-9],
    [219, 110, 110, 2e-8],
    [208, 110, 110, 6e-9],
    [220, 110, 110, 1e-8],
    [4, 100, 100, 1e-8],
    [217, 100, 198, 2e-8],
    [129, 100, 100, 1e-9],
    [21, 100, 100, 2e-7],
    [228, 100, 100, 1e-11],
    [210, 100, 100, 2e-8],
    [211, 100, 100, 2e-8],
    [112, 100, 100, 9e-9],
    [124, 100, 200, 5e-9],
    [123, 100, 101, 8e-7],
    [218, 100, 100, 3e-8],
    [113, 100, 100, 2e-8],
    [5, 100, 198, 1e-8],
    [213, 100, 164, 3e-8],
    [121, 100, 100, 4e-9],
    [125, 100, 102, 1e-8],
    [126, 100, 100, 3e-7],
    [214, 99, 99, 9e-8],
    [215, 99, 99, 9e-8],
    [212, 100, 100, 0],
    [126, 100, 100, 0],
    [209, 100, 198, 0],
    [228, 100, 100, 0]
])

noiseflag = 0

nrows = Var.shape[0]
probtype = 'smooth'
probspecs = {
    'trunc': 10**16
}



namestr = []

if noiseflag:
    sys.path.append('~/repos/randprojections21/src/testfuncs/')
    print('num     prob     n     m            f0     noise     hopt    time')
else:
    print('Problem     n    m            f0     h')

file_name = "pymidfuns100.txt"
for i in range(40):      #26
    probspecs['nprob'] = int(Var[i, 0])
    probspecs['n'] = int(Var[i, 1])
    probspecs['m'] = int(Var[i, 2])
    
    factor = 10**0  # revisit!
    X0, prob = dfoxsnew(probspecs['m'], probspecs['n'], probspecs['nprob'])
    # X0 = X0 + 100 * np.ones(X0.shape)
    namestr.append(prob['name'])
    
    X0 = factor * X0
    
    start_time = time.time()
    y, fvec = calfun_sample(X0, probspecs, probtype, file_name)
    ti = time.time() - start_time
    
    if noiseflag:
        h = 1e-11
        nf = 23
        np.random.seed(1)  # To replicate the random state

        p = np.random.rand(probspecs['n'])  # Direction to compute derivative
        fval = np.zeros(nf)
        for j in range(nf):
            fval[j] = calfun_sample(X0 + h * (j - 7) * p, probspecs, probtype)
        
        # Compute noise estimate
        fnoise, level, inform = ecnoise(nf, fval)
        fder2, s2n = calfun_sample(X0, probspecs, probtype)
        hopt = 1.68 * np.sqrt(fnoise / np.abs(fder2))
        
        p = np.random.rand(probspecs['n'])  # Compute derivative again
        fval = np.zeros(nf)
        for j in range(nf):
            fval[j] = calfun_sample(X0 + h * (j - 7) * p, probspecs, probtype)
        
        # Compute noise estimate
        fnoise2, level, inform = ecnoise(nf, fval)
        fder2, s2n = calfun_sample(X0, probspecs, probtype)
        hopt2 = 1.68 * np.sqrt(fnoise2 / np.abs(fder2))
        
        p = np.random.rand(probspecs['n'])  # Compute derivative again
        fval = np.zeros(nf)
        for j in range(nf):
            fval[j] = calfun_sample(X0 + h * (j - 7) * p, probspecs, probtype)
        
        # Compute noise estimate
        fnoise3, level, inform = ecnoise(nf, fval)
        fder2, s2n = calfun_sample(X0, probspecs, probtype)
        hopt3 = 1.68 * np.sqrt(fnoise3 / np.abs(fder2))
        
        fnoise0 = (fnoise + fnoise2 + fnoise3) / 3
        hopt0 = (hopt + hopt2 + hopt3) / 3
        
        print(f"{namestr[i]} {probspecs['n']} {probspecs['m']} {y} {fnoise0[i]} {hopt0[i]}")
    else:
        print(f"{i} {namestr[i]} {probspecs['n']} {probspecs['m']} {y} {prob['h']}")
if noiseflag:
    ratio_hopt = np.min(np.min(np.array([hopt, hopt2, hopt3]).T, axis=1) / np.max(np.array([hopt, hopt2, hopt3]).T, axis=1))
    ratio_fnoise = np.min(np.min(np.array([fnoise, fnoise2, fnoise3]).T, axis=1) / np.max(np.array([fnoise, fnoise2, fnoise3]).T, axis=1))
    hopt = np.round(hopt, decimals=-int(np.floor(np.log10(abs(hopt))) - 1))
    with open('hvals_mid.pkl', 'wb') as file:
        pickle.dump((hopt, namestr), file)
