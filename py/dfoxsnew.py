import numpy as np
from math import sqrt, sin, cos, log
from scipy.sparse import eye
def dfoxsnew(m, n, nprob):
    x = np.zeros(n)  # Initialize x as a zero array of length n
    prob = {}
    # Error messages
    eid = 'Input:dimensionIncompatible'    
    err = [
        'Input m must be equal to n.',
        'Input m must be at least n.',
        'Input n must be even.',
        'Input sqrt(n) must be an integer.',
        'Input nthroot(n,3) must be an integer.',
        'Input sqrt(n/2) must be an integer.',
        'Input m must be equal to 2*n-2.',
        'Input m must be equal to 2*(n-4).',
        'Input n must be at least 5.',
        'Input n must be at least 2.',
        'Input (sqrt(1+4*n)-1)/2 must be an integer.',
        'Input m must be equal to n+1.',
        'Input m must be equal to 2*n.',
        'Input n/4 must be an integer.',
        'Input n must be at least 13.',
        'Input (n+2)/3 must be an integer.',
        'Input m must be equal to (5*n-8)/3.',
        'Input m must be equal to n+2.',
        'Input sqrt(n+1) must be an integer.'
    ]
    
    # Switch statement equivalent in Python
    if nprob == 1:  # Case for MGH 32 linear function - full rank or rank 1
        if m < n:
            raise ValueError(eid, err[1])  # Error if m < n
        x = np.ones(n)  # Set x to a vector of ones of length n
        prob['name'] = 'ARGLALE'  # Name of the problem
        prob['h'] = 2e-7  # Some parameter
        prob['hm'] = 4000  # Another parameter
        prob['hn'] = 2000  # Another parameter
        prob['fbest'] = m - n  # Compute fbest as m - n
    elif nprob == 2:  # Case for MGH 33 linear function - full rank or rank 1
        if m < n:
            raise ValueError(eid, err[1])  # Error if m < n
        x = np.ones(n)  # Set x to a vector of ones of length n
        prob['name'] = 'ARGLBLE'  # Name of the problem
        prob['h'] = 1e-7  # Some parameter
        prob['hm'] = 4000  # Another parameter
        prob['hn'] = 2000  # Another parameter
        prob['fbest'] = m * (n - 1) / (2 * (2 * m + 1))  # Compute fbest

    elif nprob == 3:  # MGH 34 
        if m < n:
            raise ValueError(eid, err[1])  
        x = np.ones(n)  
        prob['name'] = 'ARGLCLE'  
        prob['h'] = 1e-7  
        prob['hm'] = 4000  
        prob['hn'] = 2000  
        prob['fbest'] = (m**2 + 3 * m - 6) / (2 * (2 * m - 3))  

    elif nprob == 4:  # MGH 1 Rosenbrock
        if m != n:
            raise ValueError(eid, err[0])  
        x = -np.ones(n)  
        prob['name'] = 'EXTROSNB' 
        prob['h'] = 2e-8  
        prob['hm'] = 1000  
        prob['hn'] = 1000  
        prob['fbest'] = 0  
        prob['xbest'] = np.zeros(n)  

    elif nprob == 5:  
        if m != 2 * n - 2:
            raise ValueError(eid, err[6])  
        x = -np.ones(n)  
        prob['name'] = 'ROSENBR'  
        prob['h'] = 3e-8  
        prob['hm'] = 1998  
        prob['hn'] = 1000  
        prob['fbest'] = 0  
        prob['xbest'] = np.ones(n)  

    elif nprob == 121:  
        if m != n:
            raise ValueError(eid, err[0])  
        if n % 2 != 0:
            raise ValueError(eid, err[2])  
        x[0] = -1.2
        x[1] = 1
        x = np.tile(x[:2], (n // 2, 1)).flatten()  
        prob['name'] = 'SROSENBR'  
        prob['h'] = 6e-9  
        prob['hm'] = 2000  
        prob['hn'] = 2000  
        prob['fbest'] = 0  
        prob['xbest'] = np.ones(n)  

    elif nprob == 15:  
        if m < n:
            raise ValueError(eid, err[1])  
        for k in range(n):
            x[k] = (k + 1) / (n + 1)  
        prob['name'] = 'CHEBYQAD'  
        prob['h'] = 5e-11  
        prob['hm'] = 1000  
        prob['hn'] = 1000  

    elif nprob == 16:  # MGH 27 Brown 
        if m != n:
            raise ValueError(eid, err[0])  
        x = 0.5 * np.ones(n)  
        prob['name'] = 'BROWNALE' 
        prob['h'] = 7e-8  
        prob['hm'] = 1000  
        prob['hn'] = 1000  
        prob['fbest'] = 0  
        prob['xbest'] = np.ones(n)  

    elif nprob == 19:  # BDQRTIC
        if n < 5:
            raise ValueError(eid, err[8])  
        if m != 2 * (n - 4):
            raise ValueError(eid, err[7])  
        x = np.ones(n)  
        prob['name'] = 'BDQRTIC'  
        prob['h'] = 2e-8 
        prob['hm'] = 3992  
        prob['hn'] = 2000  

    elif nprob == 120:  # CUBE, cubic version of Rosenbrock
        if m != n:
            raise ValueError(eid, err[0])  
        x = np.concatenate([[-1.2], np.ones(n - 1)]) 
        prob['name'] = 'CUBE'  
        prob['h'] = 3e-9  
        prob['hm'] = 2000 
        prob['hn'] = 2000  
        prob['fbest'] = 0 
        prob['xbest'] = np.ones(n)  

    elif nprob == 21:  # MANCINO
        if m != n:
            raise ValueError(eid, err[0])  # If m is not equal to n, raise an error
        if n < 2:
            raise ValueError(eid, err[9])  # If n is less than 2, raise an error
        for i in range(n):
            ss = 0
            for j in range(1,n+1):
                ss += np.sqrt((i+1) / j) * (np.sin(np.log(np.sqrt((i+1) / j)))**5 + np.cos(np.log(np.sqrt((i+1) / j)))**5)
            x[i] = -8.710996e-4 * ((i + 1 - 50)**3 + ss)
        prob['name'] = 'MANCINO'  # Problem name
        prob['h'] = 2e-3  # Parameter h
        prob['hm'] = 1000  # Parameter hm
        prob['hn'] = 1000  # Parameter hn

    elif nprob == 217:  # FREUROTH
        if m != 2 * n - 2:
            raise ValueError(eid, err[6])  # If m is not equal to 2 * (n - 4), raise an error
        if n % 2 != 0:
            raise ValueError(eid, err[2])  # If n is not even, raise an error
        x[::2] = 0.5  # Set every second element to 0.5 (odd indices)
        x[1::2] = -2  # Set every second element to -2 (even indices)
        prob['name'] = 'FREUROTH'  # Problem name
        prob['h'] = 6e-8  # Parameter h
        prob['hm'] = 9998  # Parameter hm
        prob['hn'] = 5000  # Parameter hn

    elif nprob == 126:  # VarTrig
        if m != n:
            raise ValueError(eid, err[0])  # If m is not equal to n, raise an error
        x = np.ones(n) / n  # Set x as a vector of 1/n
        prob['name'] = 'VarTrig'  # Problem name
        prob['h'] = 2e-7  # Parameter h
        prob['hm'] = 1000  # Parameter hm
        prob['hn'] = 1000  # Parameter hn
        prob['fbest'] = 0  # Set fbest
        prob['xbest'] = np.zeros(n)  # Set xbest as a zero vector


    elif nprob == 201:  # ARTIFICIAL TURNING POINT
        if m != n:
            raise ValueError(eid, err[0])  # m must be equal to n
        x = np.ones(n)  # Set x to be a vector of ones
        prob['name'] = 'ARTIF'  # Problem name
        prob['h'] = 7e-9  # Parameter h
        prob['hm'] = 5000  # Parameter hm
        prob['hn'] = 5000  # Parameter hn
        prob['fbest'] = 0  # Best function value
        prob['xbest'] = np.zeros(n)  # Best solution (zeros)

    elif nprob == 202:  # ARWHDNE
        if m != 2 * n - 2:
            raise ValueError(eid, err[7])  # m must be equal to 2*n-2
        x = np.ones(n)  # Set x to be a vector of ones
        prob['name'] = 'ARWHDNE'  # Problem name
        prob['h'] = 2e-8  # Parameter h
        prob['hm'] = 9998  # Parameter hm
        prob['hn'] = 5000  # Parameter hn

    elif nprob == 128:  # BDVALUES
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        h = 1 / (n + 1)
        t = np.arange(1, n + 1) * h
        x = 1000 * (t * (t - 1))
        prob['name'] = 'BDVALUES'
        prob['h'] = 1e-7
        prob['hm'] = 1000
        prob['hn'] = 1000

    elif nprob == 130:  # BROYDN3D
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        x = -np.ones(n)
        prob['name'] = 'BROYDN3D'
        prob['h'] = 7e-9
        prob['hm'] = 1000
        prob['hn'] = 1000

    elif nprob == 131:  # BROYDNBD
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        x = -np.ones(n)
        prob['name'] = 'BROYDNBD'
        prob['h'] = 1e-8
        prob['hm'] = 1000
        prob['hn'] = 1000

    elif nprob == 203:  # BRATU2D
        
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        if np.sqrt(n) % 1 != 0:
            raise ValueError(eid, err[3])  # Using err array for error message
        x = np.zeros(n)
        prob['name'] = 'BRATU2D'
        prob['h'] = 7e-11
        prob['hm'] = 4900
        prob['hn'] = 4900

    elif nprob == 204:  # BRATU2DT
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        if np.sqrt(n) % 1 != 0:
            raise ValueError(eid, err[3])  # Using err array for error message
        x = np.zeros(n)
        prob['name'] = 'BRATU2DT'
        prob['h'] = 7e-11
        prob['hm'] = 4900
        prob['hn'] = 4900

    elif nprob == 205:  # BRATU3D
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        if np.cbrt(n) % 1 != 0:
            raise ValueError(eid, err[4])  # Using err array for error message
        x = np.zeros(n)
        prob['name'] = 'BRATU3D'
        prob['h'] = 7e-10
        prob['hm'] = 3375
        prob['hn'] = 3375

    elif nprob == 206:  # Complex 2D version (CBRATU2D)
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        if np.sqrt(n / 2) % 1 != 0:
            raise ValueError(eid, err[5])  # Using err array for error message
        x = np.zeros(n)
        prob['name'] = 'CBRATU2D'
        prob['h'] = 8e-11
        prob['hm'] = 2888
        prob['hn'] = 2888

    elif nprob == 207:  # CHANDHEQ
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        x = np.ones(n)
        prob['name'] = 'CHANDHEQ'
        prob['h'] = 1e-8
        prob['hm'] = 1000
        prob['hn'] = 1000

    elif nprob == 208:  # EIGENB
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        if (np.sqrt(1 + 4 * n) - 1) / 2 % 1 != 0:
            raise ValueError(eid, err[10])  # Using err array for error message
        p = 0.5 * (np.sqrt(1 + 4 * n) - 1)
        Q = eye(int(p))
        x = np.concatenate([Q.toarray().flatten(), np.ones(int(p))])
        prob['name'] = 'EIGENB'
        prob['h'] = 9e-9
        prob['hm'] = 2550
        prob['hn'] = 2550

        # Create matrix A and compute eigenvalues and eigenvectors
        A = np.diag(2 * np.ones(int(p))) - np.diag(np.ones(int(p - 1)), -1) - np.diag(np.ones(int(p - 1)), 1)
        V, D = np.linalg.eig(A)
        prob['xbest'] = np.concatenate([D.flatten(), V])
        prob['fbest'] = 0  # (1e-28 for n=2550 in Matlab)

    elif nprob == 219:  # EIGENA
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        if (np.sqrt(1 + 4 * n) - 1) / 2 % 1 != 0:
            raise ValueError(eid, err[10])  # Using err array for error message
        p = 0.5 * (np.sqrt(1 + 4 * n) - 1)  # This is what must be integer
        Q = eye(p, format='dense')
        x = np.concatenate([np.array(Q.flatten()).reshape([-1]).flatten(), np.ones(int(p))])
        prob['name'] = 'EIGENA'
        prob['h'] = 1e-8
        prob['hm'] = 3660
        prob['hn'] = 3660
        prob['xbest'] = np.concatenate([np.eye(int(p)).reshape(int(p)*int(p)), np.arange(1, int(p) + 1)])
        prob['fbest'] = 0

    elif nprob == 220:  # EIGENC
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        if (np.sqrt(1 + 4 * n) - 1) / 2 % 1 != 0:
            raise ValueError(eid, err[10])  # Using err array for error message
        p = 0.5 * (np.sqrt(1 + 4 * n) - 1)  # This is what must be integer
        Q = eye(int(p), format='dense')
        x = np.concatenate([np.array(Q.flatten()).reshape([-1]).flatten(), np.ones(int(p))])
        prob['name'] = 'EIGENC'
        prob['h'] = 1e-8
        prob['hm'] = 2550
        prob['hn'] = 2550
        A = np.diag(np.arange(p, 0, -1)) + np.diag(np.ones(int(p) - 1), -1) + np.diag(np.ones(int(p) - 1), 1)
        V, D = np.linalg.eig(A)
        prob['xbest'] = np.concatenate([D.flatten(), V])
        prob['fbest'] = 0

    elif nprob == 129:  # INTEGREQ
        if m != n:
            raise ValueError(eid, err[0])  # Using err array for error message
        h = 1 / (n + 1)
        x = h * (np.arange(1, n+1) * (h * np.arange(1, n+1) - 1))
        prob['name'] = 'INTEGREQ'
        prob['h'] = 2e-9
        prob['hm'] = 2550
        prob['hn'] = 2550

    elif nprob == 228:
        if m != n:
            raise ValueError(f"{eid}: {err[0]}")
        h = 1 / (n + 1)
        t = np.arange(1, n+1).reshape(-1, 1) * h  # Create t as an (n, 1) array
        x = t * (t - 1)
        prob = {
            'name': 'MOREBV',
            'h': 5e-13,
            'hm': 1000,
            'hn': 1000,
        }

    elif nprob == 210:  # Case for MSQRTA
        if m != n:
            raise ValueError(eid, err[0])  # Error if m != n
        if int(sqrt(n)) != sqrt(n):  # Check if sqrt(n) is an integer
            raise ValueError(eid, err[3])  # Error if sqrt(n) is not an integer

        Npar = int(sqrt(n))  # Calculate Npar
        xx = np.sin(np.arange(1, n+1)**2)  # sin([1:n].^2)
        xx[Npar * 2] = 0  # Set the value at the specific index to 0 (corresponds to MATLAB indexing)
        x = (xx - 0.8 * np.sin(np.arange(1, n+1)**2))  # Compute x as per the formula
        prob['name'] = 'MSQRTA'  # Name of the problem
        prob['h'] = 7e-8  # Some parameter
        prob['hm'] = 4900  # Another parameter
        prob['hn'] = 4900  # Another parameter
        prob['xbest'] = np.sin(np.arange(1, n+1)**2)  # Initialize prob.xbest
        prob['xbest'][2 * Npar] = 0  # Set specific index to 0
        prob['fbest'] = 0  # Set fbest to 0

    elif nprob == 211:  # Case for MSQRTB
        if m != n:
            raise ValueError(eid, err[0])  # Error if m != n
        if int(sqrt(n)) != sqrt(n):  # Check if sqrt(n) is an integer
            raise ValueError(eid, err[3])  # Error if sqrt(n) is not an integer

        xx = np.sin(np.arange(1, n+1)**2)  # sin([1:n].^2)
        x = (xx - 0.8 * np.sin(np.arange(1, n+1)**2))  # Compute x as per the formula
        prob['name'] = 'MSQRTB'  # Name of the problem
        prob['h'] = 9e-8  # Some parameter
        prob['hm'] = 1024  # Another parameter
        prob['hn'] = 1024  # Another parameter
        prob['xbest'] = np.sin(np.arange(1, n+1)**2)  # Initialize prob.xbest
        prob['fbest'] = 0  # Set fbest to 0

    elif nprob == 112:  # Case for OSCIGRNE
        if m != n:
            raise ValueError(eid, err[0])  # Error if m != n

        x = np.zeros(n)  # Initialize x as a zero array of length n
        x[0] = -2  # Set the first element to -2, the rest stay 1
        x[1:] = 1  # Set the remaining elements to 1

        prob['name'] = 'OSCIGRNE'  # Name of the problem
        prob['h'] = 4e-9  # Some parameter
        prob['hm'] = 1000  # Another parameter
        prob['hn'] = 1000  # Another parameter
        prob['xbest'] = np.ones(n)  # Initialize prob.xbest to a vector of ones
        prob['fbest'] = 0  # Set fbest to 0

    elif nprob == 123:  # Case for PENLT1NE (MGH 23)
        if m != n + 1:
            raise ValueError(eid, err[11])  # Error if m != n + 1

        x = np.arange(1, n+1)  # Create a vector [1, 2, ..., n]

        prob['name'] = 'PENLT1NE'  # Name of the problem
        prob['h'] = 1e-5  # Some parameter
        prob['hm'] = 1001  # Another parameter
        prob['hn'] = 1000  # Another parameter

    elif nprob == 124:  # Case for Penalty2 (MGH 24)
        if m != 2 * n:
            raise ValueError(eid, err[12])  # Error if m != 2 * n

        x = np.full(n, 0.5)  # Create a vector of length n filled with 0.5

        prob['name'] = 'Penalty2'  # Name of the problem
        prob['h'] = 5e-9  # Some parameter
        prob['hm'] = 2000  # Another parameter
        prob['hn'] = 1000  # Another parameter

    elif nprob == 113:  # Case for POWELLSG (MGH 13)
        if m != n:
            raise ValueError(eid, err[0])  

        if n % 4 != 0:
            raise ValueError(eid, err[13])  

        
        x = np.array([3, -1, 0, 1])
        x = np.tile(x, (n // 4, 1)).flatten()  

        prob['name'] = 'POWELLSG'  
        prob['h'] = 4e-8  
        prob['hm'] = 2000  
        prob['hn'] = 2000 

    elif nprob == 218:  # Case for POWELLSE (MGH 13)
        if m != n:
            raise ValueError(eid, err[0])  

        if n % 4 != 0:
            raise ValueError(eid, err[13])  

        x = np.array([3, -1, 0, 1])
        x = np.tile(x, (n // 4, 1)).flatten()  

        prob['name'] = 'POWELLSE'  
        prob['h'] = 5e-8  
        prob['hm'] = 1000  
        prob['hn'] = 1000  

    elif nprob == 213:  # Case for SPMSQRT
        if n < 13:
            raise ValueError(eid, err[14])  # Error if n < 13

        if (n + 2) % 3 != 0:
            raise ValueError(eid, err[15])  # Error if (n+2)/3 is not an integer

        if m != (5 * n - 8) / 3:
            raise ValueError(eid, err[16])  # Error if m != (5*n - 8)/3

        # Initialize x
        x = 0.2 * np.sin(np.arange(1, n + 1) ** 2)

        # Set problem parameters
        prob['name'] = 'SPMSQRT'
        prob['h'] = 4e-8
        prob['hm'] = 1664
        prob['hn'] = 1000

    elif nprob == 125:  # Case for VARDIMNE (MGH 25)
        if m != n + 2:
            raise ValueError(eid, err[17])  # Error if m != n + 2

        # Initialize x
        x = 1 - np.arange(1, n + 1) / n  # Create a vector [1/n, 2/n, ..., n/n] and subtract from 1

        # Set problem parameters
        prob['name'] = 'VARDIMNE'
        prob['h'] = 2e-8
        prob['hm'] = 1002
        prob['hn'] = 1000
        prob['xbest'] = np.ones(n)  # Vector of ones of length n
        prob['fbest'] = 0

    elif nprob == 214:  # Case for YATP1SQ (MGH variant of YATP1LS)
        if m != n:
            raise ValueError(eid, err[0])  # Error if m != n

        if int(np.sqrt(n + 1)) ** 2 != n + 1:
            raise ValueError(eid, err[18])  # Error if sqrt(n + 1) is not an integer

        Npar = int(np.sqrt(n + 1) - 1)
        x = np.concatenate([np.zeros(2 * Npar), 6 * np.ones(Npar**2)])

        prob['name'] = 'YATP1SQ'
        prob['h'] = 2e-7
        prob['hm'] = 2600
        prob['hn'] = 2600

    elif nprob == 216:  # Case for ConnBand (m=n)
        if m != n:
            raise ValueError(eid, err[0])  # Error if m != n
        
        x = np.zeros(n)  # Initialize x as a zero array of length n
        x[0::2] = 1  # Set the odd-indexed elements to 1
        x[1::2] = -1  # Set the even-indexed elements to -1

        prob['name'] = 'ConnBand'
        prob['h'] = 3e-9
        prob['hm'] = 1000
        prob['hn'] = 1000
        prob['xbest'] = np.zeros(n)  # Initialize xbest as a zero array of length n
        prob['fbest'] = 0

    elif nprob == 209:  # Case for FREURONE (freudenstein and roth function)
        if m != 2 * n - 2:
            raise ValueError(eid, err[6])  # Error if m != 2 * n - 2
        if n % 2 != 0:
            raise ValueError(eid, err[2])  # Error if n is not even
        x[::2] = 0.5  # Set every other element starting from index 0 to 0.5
        x[1::2] = -2   # Set every other element starting from index 1 to -2
        prob['name'] = 'FREURONE'  # Name of the problem

    elif nprob == 212:
        x = np.zeros(n)

    elif nprob == 215:       # YATP2SQ (Think this is YATP2LS)
        if m != n:
            raise ValueError(err[0].format(eid))

           
        if not np.isclose(np.sqrt(n + 1) % 1, 0):
            raise ValueError(err[18]) 

         
        Npar = int(np.sqrt(n + 1) - 1)

       
        x = np.concatenate((np.zeros(2 * Npar), 10 * np.ones(Npar ** 2)))

     
        prob = {
            'name': 'YATP2SQ',  # CUTEST
            'h': 1e-7,
            'hm': 2600,
            'hn': 2600
        }


    else:
        print('Unknown problem.')
        del x




    return x, prob  # Return the x array and error messages

# Example of calling the function
