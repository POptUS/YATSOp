import numpy as np
from mghvec import mghvec
def calfun_sample(x, probspecs, probtype):
    # Ensure input is a column vector
    if x.shape[0] == 1:
        x = x[:, np.newaxis]

    nin = x.shape[0]  # Problem dimension

    if nin != probspecs['n']:
        raise ValueError("Input x is not of size n by 1.")

    # Generate the vector
    fvec = mghvec(probspecs['m'], probspecs['n'], x, probspecs['nprob'])

    # Optional truncation
    if "trunc" in probspecs.keys() and np.max(np.abs(fvec)) > probspecs['trunc']:
        fvec = np.sign(fvec) * np.minimum(np.abs(fvec), probspecs['trunc'])

    # Optional for noisy problems
    if "sigma" not in probspecs:  # Default for noisy problems
        probspecs['sigma'] = 10**-3
    
    if "seed" in probspecs:  # If seed is specified
        np.random.seed(probspecs['seed'])

    # Calculate the function value
    if probtype == "absnormal":
        z = probspecs['sigma'] * np.random.randn(probspecs['m'], 1)
        fvec += z
        y = np.sum(fvec**2)
    
    elif probtype == "absuniform":
        z = (probspecs['sigma'] * np.sqrt(3)) * (2 * np.random.rand(probspecs['m'], 1) - 1)
        fvec += z
        y = np.sum(fvec**2)
    
    elif probtype == "abswild":
        z = 0.9 * np.sin(100 * np.linalg.norm(x, 1)) * np.cos(100 * np.linalg.norm(x, np.inf)) + 0.1 * np.cos(np.linalg.norm(x, 2))
        z = z * (4 * z**2 - 3)
        y = np.sum(fvec**2) + z

    elif probtype == "nondiff":
        y = np.sum(np.abs(fvec))

    elif probtype == "relnormal":
        z = probspecs['sigma'] * np.random.randn(probspecs['m'], 1)
        fvec *= (1 + z)
        y = np.sum(fvec**2)
    
    elif probtype == "reluniform":
        z = (probspecs['sigma'] * np.sqrt(3)) * (2 * np.random.rand(probspecs['m'], 1) - 1)
        fvec *= (1 + z)
        y = np.sum(fvec**2)

    elif probtype == "absnormal2":
        z = probspecs['sigma'] * np.random.randn()
        y = np.sum(fvec**2) + z

    elif probtype == "absuniform2":
        z = (probspecs['sigma'] * np.sqrt(3)) * (2 * np.random.rand() - 1)
        y = np.sum(fvec**2) + z

    elif probtype == "reluniform2":
        z = (probspecs['sigma'] * np.sqrt(3)) * (2 * np.random.rand() - 1)
        y = np.sum(fvec**2) * (1 + z)

    elif probtype == "relnormal2":
        z = probspecs['sigma'] * np.random.randn()
        y = np.sum(fvec**2) * (1 + z)

    elif probtype == "relwild":
        z = 0.9 * np.sin(100 * np.linalg.norm(x, 1)) * np.cos(100 * np.linalg.norm(x, np.inf)) + 0.1 * np.cos(np.linalg.norm(x, 2))
        z = z * (4 * z**2 - 3)
        y = (1 + probspecs['sigma'] * z) * np.sum(fvec**2)

    elif probtype == "smooth":
        y = np.sum(np.array(fvec)**2)

    return y, fvec
