import numpy as np
import time
from calfun_sample import calfun_sample
from ecnoise import ecnoise
from dfoxsnew import dfoxsnew

# Define the array of variables, equivalent to 'Var' in MATLAB
Var = np.array([
    [1, 2000, 4000, 0],   # ARGLALE
    [2, 2000, 4000, 0],   # ARGLBLE
    [201, 5000, 5000, 0], # ARTIF
    [202, 5000, 9998, 0], # ARWHDNE
    [128, 1000, 1000, 0], # BDVALUES + Off in 5th digit but calling OK
    [203, 4900, 4900, 0], # BRATU2D
    [204, 4900, 4900, 0], # BRATU2DT
    [205, 3375, 3375, 0], # BRATU3D
    [16, 1000, 1000, 0],  # BROWNALE
    [130, 1000, 1000, 0], # BROYDN3D
    [131, 5000, 5000, 0], # BROYDNBD + OK disagrees with table but agrees with OPM
    [206, 2888, 2888, 0], # CBRATU2D
    [207, 1000, 1000, 0], # CHANDHEQ
    [208, 2550, 2550, 0], # EIGENB
    [217, 5000, 9998, 0], # FREUROTH + (replaces FREURONE)
    [129, 1000, 1000, 0], # INTEGREQ
    [228, 1000, 1000, 0], # MOREBV   + (replaces MOREBVNE)
    [210, 4900, 4900, 0], # MSQRTA   + Off in 5th digit but calling OK
    [211, 1024, 1024, 0], # MSQRTB   ! Off in 3rd digit
    [112, 1000, 1000, 0], # OSCIGRNE - matches, but should confirm another value with LR
    [123, 1000, 1001, 0], # PENLT1NE
    [218, 1000, 1000, 0], # POWELLSE
    [213, 1000, 1664, 0], # SPMSQRT
    [125, 1000, 1002, 0], # VARDIMNE
    [214, 2600, 2600, 0], # YATP1SQ
    [215, 2600, 2600, 0], # YATP2SQ
    [126, 1000, 1000, 0], # VarTrig  (replaces ARGTRIG)
    [216, 1000, 1000, 0], # ConnBand (new)
    [113, 2000, 2000, 0], # POWELLSG (new)
    [15, 1000, 1000, 0],  # CHEBYQAD (new)
    [124, 1000, 2000, 0], # Penalty2 (new)
    [3, 2000, 4000, 0],   # ARGLCLE  (new)
    [4, 1000, 1000, 0],   # EXTROSNB (new)
    [5, 1000, 1998, 0],   # ROSENBR  (new)
    [121, 2000, 2000, 0], # SROSENBR (new)
    [19, 2000, 3992, 0],  # BDQRTIC  (new)
    [120, 2000, 2000, 0], # CUBE     (new)
    [21, 1000, 1000, 0],  # MANCINO  (new)
    [219, 3660, 3660, 0], # EIGENA   (new)
    [220, 2550, 2550, 0], # EIGENC   (new)
    [212, 1000, 1000, 0], # SEMICN2U * More' - missing
    [126, 1000, 1000, 0], # ARGTRIG  X Agrees with OPM, replaced by VarTrig
    [209, 5000, 9998, 0], # FREURONE X Disagrees, replaced by FREUROTH
    [228, 1000, 1000, 0]  # MOREBVNE X Not yet in agreement, replaced with MOREBV
])

# Set flags and parameters
noiseflag = 0
nrows = Var.shape[0]
probtype = "smooth"
probspecs = {'trunc': 10**16}


# Print headers
if noiseflag:
    print('num     prob     n     m            f0     noise     hopt    time')
else:
    print('Problem     n    m            f0     h')


namestr = []
file_name = "pydfof.txt"
# Loop over rows in Var (similar to the for-loop in MATLAB)
for i in range(40):  # Python uses 0-based indexing
    probspecs['nprob'] = int(Var[i, 0])
    probspecs['n'] = int(Var[i, 1])
    probspecs['m'] = int(Var[i, 2])
    
    factor = 10 ** (Var[i,3])  # Placeholder value

    # Get starting point and problem information
    [X0, prob] = dfoxsnew(probspecs['m'], probspecs['n'], probspecs['nprob'])
    # X0 = X0 + 100 * np.ones(X0.shape)
    namestr.append(prob['name'])  # Replace with actual name retrieval logic
    X0 = factor * X0

    # Start timing
    start_time = time.time()

    # Compute the objective function (calculation)
    y, fvec = calfun_sample(X0, probspecs, probtype, file_name)

    # Stop timing
    elapsed_time = time.time() - start_time

    # If noise flag is on, perform noise-related calculations
    if noiseflag:
        h = 1e-11
        nf = 23
        np.random.seed(1)  # Ensure reproducibility

        # Sample noise values and compute estimates
        p = np.random.rand(probspecs['n'])
        fval = [calfun_sample(X0 + h * (j - 7) * p, probspecs, probtype, file_name) for j in range(nf)]
        fnoise, _, _ = ecnoise(nf, fval)
        fder2, s2n = calfun_sample(X0, probspecs, probtype, file_name)
        hopt = 1.68 * np.sqrt(fnoise / abs(fder2))

        fnoise2, _, _ = ecnoise(nf, fval)
        fder2, s2n = calfun_sample(X0, probspecs, probtype, file_name)
        hopt2 = 1.68 * np.sqrt(fnoise2 / abs(fder2))

        fnoise3, _, _ = ecnoise(nf, fval)
        fder2, s2n = calfun_sample(X0, probspecs, probtype, file_name)
        hopt3 = 1.68 * np.sqrt(fnoise3 / abs(fder2))

        # Average the results
        fnoise0 = (fnoise + fnoise2 + fnoise3) / 3
        hopt0 = (hopt + hopt2 + hopt3) / 3

        # Print the results
        print(f'{probspecs["nprob"]} {namestr[0]} {probspecs["n"]} {probspecs["m"]} {y} {fnoise0[i]} {hopt0[i]} {elapsed_time}')
    else:
        # Print the results without noise calculations
        print(f'{i} {probspecs["nprob"]} {namestr[-1]} {probspecs["n"]} {probspecs["m"]} {y} {prob["h"]} {elapsed_time}')

if noiseflag:
    # Compute the min of the ratio of hopt, hopt2, hopt3 and fnoise, fnoise2, fnoise3
    hopt_min_ratio = np.min(np.min(np.vstack([hopt, hopt2, hopt3]), axis=1) / np.max(np.vstack([hopt, hopt2, hopt3]), axis=1))
    fnoise_min_ratio = np.min(np.min(np.vstack([fnoise, fnoise2, fnoise3]), axis=1) / np.max(np.vstack([fnoise, fnoise2, fnoise3]), axis=1))

    # Round hopt to 1 significant digit
    hopt = np.round(hopt, 1)

    # Save the values of hopt and namestr to a file (similar to the MATLAB save function)
    np.savez('hvals_mid.npz', hopt=hopt, namestr=namestr)
