import numpy as np


def ecnoise(nf, fval):
    """
    Generate noise similar to MATLAB's ECnoise function.

    Parameters:
    - nf: Number of frequency points.
    - fval: Frequency values (numpy array).

    Returns:
    - fnoise: Generated noise array.
    - level: Standard deviation (level) of the noise.
    - inform: Information array (initialized to zero).
    """
    # Ensure fval is a numpy array
    fval = np.array(fval)

    # Check if the length of fval matches nf
    if len(fval) != nf:
        raise ValueError("Length of fval must match nf")

    # Generate random noise
    fnoise = np.random.randn(nf)

    # Calculate the standard deviation (level) of the noise
    level = np.std(fnoise)

    # Initialize information array
    inform = np.zeros(nf)

    return fnoise, level, inform
