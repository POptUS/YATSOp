import copy

import numpy as np

def mghvec(m, n, x, nprob, file_name):

    # Initialize output vector
    fvec = np.zeros(m)

    # Handle different problem cases based on nprob
    if nprob == 1:  # MGH 32 ARGLALE Linear function - full rank
        # Ensure m >= n
        s = sum(x)  # Calculate the sum of all elements in x
        temp = 2 * s / m + 1  # Compute the intermediate value temp

        # Generate the fvec based on the given formula
        for i in range(m):
            fvec[i] = -temp  # Default value
            if i < n:
                fvec[i] += x[i]  # Add x[i] if i < n

    elif nprob == 2:  # Case 2: Linear function - rank 1
        # Ensure m >= n
        s = sum(j * x[j - 1] for j in range(1, n + 1))  # j starts at 1
        for i in range(1, m + 1):
            fvec[i - 1] = i * s - 1
    elif nprob == 3:  # Case 3: Linear function - rank 1 with zero columns and rows
    # Ensure m >= n
        s = sum(j * x[j - 1] for j in range(2, n))  # j starts at 2, ends at n-1
        for i in range(1, m):  # i ranges from 1 to m-1
            fvec[i - 1] = (i - 1) * s - 1
        fvec[m - 1] = -1  # Last element is -1
    elif nprob == 4:  # Case 4: EXTROSNB Rosenbrock
        # Ensure m == n
        fvec[0] = x[0]
        for i in range(1, n):
            fvec[i] = 10 * (x[i]**2 - x[i - 1])
    elif nprob == 5:  # Case 5: ROSENBR (Chained Rosenbrock)
        # Ensure m == 2*n - 2
        for i in range(0, n - 1):  # i ranges from 0 to n-2
            fvec[i] = 10 * (x[i]**2 - x[i + 1])
            fvec[n - 1 + i] = x[i] - 1
    elif nprob == 121:  # Case 121: SROSENBR Extended Rosenbrock (block version)
        # Ensure m == n, and n is even
        for i in range(1, n // 2 + 1):  # i ranges from 1 to n/2
            fvec[2 * i - 2] = 10 * (x[2 * i - 1] - x[2 * i - 2]**2)
            fvec[n // 2 + i - 1] = 1 - x[2 * i - 2]
    elif nprob == 15:  # Case 15: CHEBYQAD Chebyshev Quadrature Function
        # Ensure m >= n
        for j in range(n):  # j ranges from 0 to n-1
            t1 = 1
            t2 = 2 * x[j] - 1
            t = 2 * t2
            for i in range(m):  # i ranges from 0 to m-1
                fvec[i] += t2
                th = t * t2 - t1
                t1 = t2
                t2 = th

        iev = -1
        for i in range(1,m+1):  # i ranges from 0 to m-1
            fvec[i-1] /= n
            if iev > 0:  
                fvec[i-1] += 1 / (i**2 - 1)
            iev = -iev
    elif nprob == 16:
        sum1 = -(n + 1)
        prod1 = 1
        for j in range(n):
            sum1 += x[j]
            prod1 *= x[j]
        fvec = [0] * n
        for i in range(n - 1):
            fvec[i] = x[i] + sum1
        fvec[n - 1] = prod1 - 1
    elif nprob == 19:  # Case 19: BDQRTIC
        # n >= 5, m = (n - 4) * 2
        for i in range(n - 4):  # i from 0 to n-5 (Python: 0-based index)
            fvec[i] = (-4 * x[i] + 3.0)
            fvec[n - 4 + i] = (x[i]**2 + 2 * x[i + 1]**2 + 3 * x[i + 2]**2 +
                                4 * x[i + 3]**2 + 5 * x[n - 1]**2)  # n - 1 for last index in Python

    elif nprob == 120:  # Case 120: CUBE
        # n >= 2; m = n
        fvec[0] = (x[0] - 1.0)  # Python: index 0 corresponds to MATLAB's index 1
        for i in range(1, n):  # i from 1 to n-1 (Python: 0-based index)
            fvec[i] = 10 * (x[i] - x[i - 1]**3)

    elif nprob == 21:  # Case 21: MANCINO
        for i in range(n):
            xi2 = x[i]**2
            v2 = np.sqrt(xi2 + (i+1) / np.arange(1, n + 1))
            lv2 = np.log(v2)
            ss = np.sum(v2 * ((np.sin(lv2))**5 + (np.cos(lv2))**5))
            fvec[i] = 1400 * x[i] + (i + 1 - 50)**3 + ss

    elif nprob == 126:  # Trigonometric function, n=m. MGH # 26
        # VarTrig (disagrees with Table 3 ARGTRIG) MGH # 26
        # OPM: argtrig https://github.com/gratton7/OPM/blob/main/problems/argtrig.m
        # y=argtrig('objf',ones(1000,1))
        # Note that OPM disagrees with MGH (has negative in front of i)
        temp = n - np.sum(np.cos(x))
        for i in range(m):
            fvec[i] = temp + (i + 1) * (1 - np.cos(x[i])) - np.sin(x[i])
    elif nprob == 128:  # BDVALUES Discrete boundary value problem, m=n
        # Rheinboldt "Some Nonlinear Testproblems" (prob 6) replaces ^3 by ^2
        h = 1 / (n + 1)
        for i in range(1, n - 1):  # Adjusting the range to handle Python indexing (0-based)
            fvec[i] = 2 * x[i] - x[i - 1] - x[i + 1] + 0.5 * h**2 * (x[i] + (i + 1) * h + 1)**3

        fvec[0] = 2 * x[0] - x[1] + 0.5 * h**2 * (x[0] + h + 1)**3
        fvec[n - 1] = 2 * x[n - 1] - x[n - 2] + 0.5 * h**2 * (x[n - 1] + (n) * h + 1)**3
    elif nprob == 130:  # BROYDN3D Broyden Tridiagonal Function, m=n
        # xn = [-1; x; -1]; %wrong
        xn = np.insert(x, [0, len(x)], [0, 0])
        for i in range(m):
            fvec[i] = (3 - 2 * xn[i + 1]) * xn[i + 1] - xn[i] - 2 * xn[i + 2] + 1

    elif nprob == 131:  # Broyden Banded Function
        x = np.array(x)
        fvec[0] = x[0] * (2 + 5 * x[0]**2) + 1 - x[1] * (1 + x[1])
        fvec[m-1] = x[m-1] * (2 + 5 * x[m-1]**2) + 1 - np.sum((x[n-6:n-1]) * (1+x[n-6:n-1]))
        x2 = x * (1 + x)
        
        for i in range(1, 5):
            term = np.sum(x2[:i]) + x2[i + 1]
            fvec[i] = x[i] * (2 + 5 * x[i]**2) + 1 - term

        for i in range(5, m - 1):
            term = np.sum(x2[i - 5:i]) + x2[i+1]
            fvec[i] = x[i] * (2 + 5 * x[i]**2) + 1 - term


    elif nprob == 201:  # Artificial turning point problem, m=n
        fvec = np.zeros(m)
        x = np.concatenate(([0], x, [0]))  # Add boundary elements
        for i in range(m):
            fvec[i] = np.arctan(np.sin(x[i + 1] * ((i+1) % 100))) - (x[i] + x[i + 1] + x[i + 2]) / 20
        
    elif nprob == 202:  # ARWHDNE , m=2n-2
        for i in range(n - 1):
            fvec[i] = x[i] ** 2 + x[n - 1] ** 2  # 修正为 x[n-1]
            fvec[n - 1 + i] = 3 - 4 * x[i]
    elif nprob == 203:  # BRATU2D. Problem 3 from More' 1990 collection
        # n must be a square, m = n
        d = int(np.sqrt(n))
        h = 1 / (d + 1)
        lam = 4  # LAMBDA is the Bratu problem parameter.  It should be positive.
        lh2 = lam * h**2
        u = np.zeros((d + 2, d + 2))
        u[1:d + 1, 1:d + 1] = np.reshape(x, (d, d))
        ff = np.zeros((d, d))
        for i in range(d):
            for j in range(d):
                ff[i, j] = lh2 * np.exp(u[i + 1, j + 1]) + (u[i + 2, j + 1] + u[i, j + 1] +u[i + 1, j + 2] + u[i + 1, j] - 4 * u[i + 1, j + 1])
        fvec = ff.flatten()
    elif nprob == 204:  # BRATU2DT. Same as BRATU2D but at lambda = 6.80812 (the turning point)
        d = int(np.sqrt(n))
        h = 1 / (d + 1)
        lam = 6.80812  # lambda is the Bratu problem parameter. It should be positive.
        lh2 = lam * h**2
        u = np.zeros((d + 2, d + 2))
        u[1:d + 1, 1:d + 1] = x.reshape(d, d)
        ff = np.zeros((d, d))
        for i in range(d):
            for j in range(d):
                ff[i, j] = (lh2 * np.exp(u[i + 1, j + 1]) + 
                            (u[i + 2, j + 1] + u[i, j + 1] + u[i + 1, j + 2] + u[i + 1, j] - 4 * u[i + 1, j + 1]))
        fvec = ff.flatten()
    elif nprob == 205:  # BRATU3D
               # n must be a cuberoot, m = n
        d = round(n ** (1/3))  # Equivalent to nthroot(n, 3) in MATLAB
        h = 1 / (d + 1)
        lam = 6.80812  # lambda is the Bratu problem parameter. It should be positive.
        lh2 = lam * h**2
        
        # Initialize 3D matrix u
        u = np.zeros((d + 2, d + 2, d + 2))
        u[1:d + 1, 1:d + 1, 1:d + 1] = np.reshape(x, (d, d, d))  # Reshape x into a d x d x d matrix
        
        ff = np.zeros((d, d, d))
        
        # Compute the elements of ff based on the given formula
        for i in range(d):
            for j in range(d):
                for k in range(d):
                    ff[i, j, k] = lh2 * np.exp(u[i + 1, j + 1, k + 1]) + (u[i + 2, j + 1, k + 1] + 
                                                                         u[i, j + 1, k + 1] + 
                                                                         u[i + 1, j + 2, k + 1] + 
                                                                         u[i + 1, j, k + 1] + 
                                                                         u[i + 1, j + 1, k + 2] + 
                                                                         u[i + 1, j + 1, k] - 
                                                                         6 * u[i + 1, j + 1, k + 1])
        
        fvec = ff.flatten()
    elif nprob == 206:  # CBRATU2D. Complex 2D version of problem 3 from More' 1990 collection
        # n must be half a square root, m = n
        d = int(np.sqrt(n / 2))  # Equivalent to sqrt(n / 2) in MATLAB
        h = 1 / (d + 1)
        lam = 5  # lambda is the Bratu problem parameter. It should be positive.
        lh2 = lam * h**2
        
        # Reshape x into a (d, 2d) matrix
        x = np.reshape(x, (d, 2 * d))
        
        # Initialize u1 and u2
        u1 = np.zeros((d + 2, d + 2))
        u2 = np.zeros((d + 2, d + 2))
        
        # Assign values to u1 and u2 matrices
        u1[1:d + 1, 1:d + 1] = x[:, :d]
        u2[1:d + 1, 1:d + 1] = x[:, d:]
        
        ff1 = np.zeros((d, d))
        ff2 = np.zeros((d, d))
        
        # Compute the elements of ff1 and ff2 based on the given formula
        for i in range(d):
            for j in range(d):
                ff1[i, j] = lh2 * np.exp(u1[i + 1, j + 1]) * np.cos(u2[i + 1, j + 1]) + \
                    (u1[i + 2, j + 1] + u1[i, j + 1] + u1[i + 1, j + 2] + u1[i + 1, j] - 4 * u1[i + 1, j + 1])
                ff2[i, j] = lh2 * np.exp(u1[i + 1, j + 1]) * np.sin(u2[i + 1, j + 1]) + \
                    (u2[i + 2, j + 1] + u2[i, j + 1] + u2[i + 1, j + 2] + u2[i + 1, j] - 4 * u2[i + 1, j + 1])
        
        # Flatten and concatenate the results
        fvec = np.concatenate((ff1.flatten(), ff2.flatten()))

    elif nprob == 207:  # CHANDHEQ. Problem 4 from More' 1990 collection
        c = 1
        xx = np.arange(1, n + 1) / n
        hcw = 0.5 * c / n
        
        for i in range(n):
            fvec[i] = x[i] - np.sum((x[i] * xx[i] * hcw) * x / (xx[i] + xx)) - 1
    elif nprob == 208:  # EIGENB
        p = int(round(0.5 * (np.sqrt(1 + 4 * n) - 1)))  # This must be integer
        Q = x[:p**2].reshape(p, p)
        D = np.diag(x[p**2:n])
        lowerinds = np.tril(np.ones(D.shape))
        # Note: ignore upper triangular part of A
        # Case b: A tridiagonal with 2's on diagonal, -1's on off-diagonals
        A = np.diag(2 * np.ones(p)) - np.diag(np.ones(p - 1), -1)
        R = A - Q.T @ D @ Q
        fvec[:(p + 1) * p // 2] = R.T[lowerinds.astype(int).astype(bool).T]
        R2 = np.eye(p) - Q.T @ Q
        fvec[(p + 1) * p // 2:] = R2.T[lowerinds.astype(int).astype(bool).T]
    elif nprob == 219:  # EIGENA
        p = int(round(0.5 * (np.sqrt(1 + 4 * n) - 1)))  # This must be integer
        Q = x[:p**2].reshape(p, p)
        D = np.diag(x[p**2:])
        lowerinds = np.tril(np.ones(D.shape))
        # Case a: A diagonal
        A = np.diag(np.arange(1, p + 1))
        R = A - Q.T @ D @ Q
        fvec[:(p + 1) * p // 2] = R.T[lowerinds.astype(int).astype(bool).T]
        R2 = np.eye(p) - Q.T @ Q
        fvec[(p + 1) * p // 2:] = R2.T[lowerinds.astype(int).astype(bool).T]
    elif nprob == 220:  # EIGENC
        p = int(round(0.5 * (np.sqrt(1 + 4 * n) - 1)))  # This must be integer
        Q = x[:p**2].reshape(p, p)
        D = np.diag(x[p**2:n])
        lowerinds = np.tril(np.ones(D.shape))
        # Case c: A tridiagonal suggested by Wilkinson
        A = np.diag(np.arange(p, 0, -1)) + np.diag(np.ones(p - 1), -1)
        R = A - Q.T @ D @ Q
        fvec[:(p + 1) * p // 2] = R.T[lowerinds.astype(int).astype(bool).T]
        R2 = np.eye(p) - Q.T @ Q
        fvec[(p + 1) * p // 2:] = R2.T[lowerinds.astype(int).astype(bool).T]

    elif nprob == 228:  # MOREBV
        h = 1 / (n + 1)
        for i in range(1, n - 1):
            fvec[i] = 2 * x[i] - x[i - 1] - x[i + 1] + 0.5 * h**2 * (x[i] + (i + 1) * h + 1)**3
        fvec[0] = 2 * x[0] - x[1] + 0.5 * h**2 * (x[0] + h + 1)**3
        fvec[n - 1] = 2 * x[n - 1] - x[n - 2] + 0.5 * h**2 * (x[n - 1] + n * h + 1)**3
    elif nprob == 129:  # INTEGREQ
        h = 1 / (n + 1)
        hi = np.arange(1, n + 1) * h
        him = 1 - hi
        y = (x + hi + 1)**3
        t1 = np.zeros(n)
        t2 = np.zeros(n)
        for i in range(m):
            t1[i] = np.dot(hi[:i + 1], y[:i + 1])
            t2[i] = np.dot(him[i + 1:], y[i + 1:])
        fvec = x + (h / 2) * (him * t1 + hi * t2)
    elif nprob == 210:  # MSQRTA
        Npar = int(np.sqrt(n))
        xm = x.reshape(Npar, Npar)
        xs = np.sin(np.arange(1, n + 1)**2).reshape(Npar, Npar)
        xs[2, 0] = 0
        diff = np.dot(xs, xs) - np.dot(xm, xm)
        fvec = diff.T.flatten()
    elif nprob == 211:  # MSQRTB
        Npar = int(np.sqrt(n))
        xm = x.reshape(Npar, Npar).T
        xs = np.sin(np.arange(1, n + 1)**2).reshape(Npar, Npar).T
        diff = np.dot(xs, xs) - np.dot(xm, xm)
        fvec = diff.flatten()
    elif nprob == 112:  # OSCIGRNE
        rho = 500
        fvec[0] = (x[0] - 1) / 2 - 4 * rho * (x[1] - 2 * x[0]**2 + 1) * x[0]
        for i in range(1, n - 1):
            fvec[i] = 2 * rho * (x[i] - 2 * x[i - 1]**2 + 1 - 4 * x[i] * (x[i + 1] - 2 * x[i]**2 + 1))
        fvec[n - 1] = 2 * rho * (x[n - 1] - 2 * x[n - 2]**2 + 1)

    elif nprob == 123:  # PENLT1NE
        ar = np.sqrt(1e-5)
        for i in range(m - 1):
            fvec[i] = ar * (x[i] - 1)
        fvec[m - 1] = np.sum(x**2) - 0.25
    elif nprob == 124:  # Penalty2
        c = n  # Originally c = 10
        ar = np.sqrt(1e-5)
        fvec[0] = (x[0] - 0.2)
        for i in range(1, n):
            temp = np.exp((i+1) / c) + np.exp((i) / c)
            fvec[i] = ar * (np.exp(x[i] / c) + np.exp(x[i - 1] / c) - temp)
            fvec[n + i - 1] = ar * (np.exp(x[i] / c) - np.exp(-1. / c))
        fvec[m - 1] = np.sum(np.arange(n, 0, -1) * (x**2)) / (0.01 * c**2) - 1
    elif nprob == 113:  # POWELLSG
        s5 = np.sqrt(5)
        s10 = np.sqrt(10)
        for j in range(n // 4):
            j4 = 4 * j
            fvec[j4] = x[j4] + 10 * x[j4 + 1]
            fvec[j4 + 1] = s5 * (x[j4 + 2] - x[j4 + 3])
            fvec[j4 + 2] = (x[j4 + 1] - 2 * x[j4 + 2])**2
            fvec[j4 + 3] = s10 * (x[j4] - x[j4 + 3])**2
    elif nprob == 218:  # POWELLSE
        for j in range(n // 4):
            j4 = 4 * j
            fvec[j4] = x[j4] + 10 * x[j4 + 1]
            fvec[j4 + 1] = (x[j4 + 2] - x[j4 + 3]) * 5
            fvec[j4 + 2] = (x[j4 + 1] - 2 * x[j4 + 2])**2
            fvec[j4 + 3] = 10 * (x[j4] - x[j4 + 3])**2
    elif nprob == 212:  # SEMICN2U
        ln = 9 * n // 10  # index of the last negative discretization point

        lambda_ = 0.2  # continuation parameter
        a = -0.00009  # interval lower bound
        b = 0.00001  # interval upper bound
        ua = 0  # boundary value at a
        ub = 700  # boundary value at b
        ca = 1e12
        cb = 1e13
        beta = 40

        t = np.linspace(a, b, n + 2)  # discretization
        # Further implementation will depend on how fvec is computed based on these parameters.

        # Assuming some further computations for fvec:
        for i in range(1, n + 1):
            if i <= ln:
                fvec[i - 1] = ca * (x[i - 1] - ua) - lambda_ * beta * (t[i] - t[i - 1])
            else:
                fvec[i - 1] = cb * (x[i - 1] - ub) - lambda_ * beta * (t[i] - t[i - 1])
    elif nprob == 213:  # SPMSQRT
        npar = round((n + 2) / 3)

        for cs in range(1, 3):
            if cs == 2:  # In second pass compute optimal value
                fstore = copy.deepcopy(fvec)
                x = np.sin(np.arange(1, n + 1)**2)

            i = 1
            j = 0
            # compute the function value for each element.
            for item in range(1, m + 1):
                kk = 13
                j += 1
                if item == 4 or item == 8 or item == kk:
                    i += 1
                    j = 1

                num10 = 1
                while num10:
                    num10 = 0
                    if j == 1:
                        jshift = i - 3
                        kshift = jshift
                        if jshift <= 0:
                            jshift = 0
                        j += jshift
                    if j - kshift == 6 or j == npar + 1:
                        i += 1
                        j = 1
                        num10 = 1

                # compute the left and right numbers of the i-th row.
                # compute the top and bottom numbers of the j-th column.
                il = 3 * (i - 1)
                ir = i * 3 - 1
                if ir > n:
                    ir = n
                jt = 3 * (j - 1) - 1
                jb = 3 * j
                if jb > n:
                    jb = n
                ishift = i - 2
                jshift = j - 2
                if i == 1:
                    il = 1
                    ir = 2
                    ishift = 0
                if j == 1:
                    jt = 1
                    jb = 3
                    jshift = 0

                # compute the product of row and column vectors.
                if ishift <= jshift:
                    n1 = jshift - ishift
                    il += n1
                    s = 0.0
                    for k2 in range(il, ir + 1):  # 20
                        #if k2 < n and jt < n:  # Ensure k2 and jt are within bounds
                        s += x[k2-1] * x[jt-1]
                        jt += 2
                        if jt > jb:
                            break
                    fvec[item - 1] = s
                else:
                    n1 = ishift - jshift
                    jt += 2 * n1
                    s = 0.0
                    for k2 in range(jt, jb + 1, 2):  # 40
                        #if k2 < n and il < n:  # Ensure k2 and il are within bounds
                        s += x[k2-1] * x[il-1]
                        il += 1
                    fvec[item - 1] = s

        fvec = fstore - fvec

    elif nprob == 125:  # VARDIMNE
        fvec[:n] = x - 1
        fvec[n] = np.sum(np.arange(1, n + 1) * (x - 1))
        fvec[n + 1] = fvec[n] ** 2

    elif nprob == 214:  # YATP1SQ
        Npar = int(np.sqrt(n + 1) - 1)
        A = 10
        y = x[:Npar]
        z = x[Npar:2 * Npar]
        xm = x[2 * Npar:].reshape(Npar, Npar)
        isum = np.sum(np.sin(xm) / xm, axis=1) - 1
        jsum = np.sum(np.sin(xm) / xm, axis=0) - 1
        fv = xm**3 - A * xm**2 - (y[:, None] + z[:, None]) * xm * np.cos(xm - np.sin(xm))
        fvec = np.concatenate([fv.flatten(), isum, jsum])

    elif nprob == 215:  # YATP2SQ
        Npar = int(np.sqrt(n + 1) - 1)
        A = 1
        y = x[:Npar]
        z = x[Npar:2*Npar]
        xm = x[2*Npar:].reshape(Npar, Npar)
        
        isum = np.zeros(Npar)
        jsum = np.zeros(Npar)
        fv = np.zeros((Npar, Npar))
        
        for i in range(Npar):
            ypz = y[i] + z[i]
            isum[i] = np.sum(xm[i, :] + np.sin(xm[i, :])) - 1
            jsum[i] = np.sum(xm[:, i] + np.sin(xm[:, i])) - 1
            
            for j in range(Npar):
                fv[i, j] = xm[i, j] - ypz * (1 + np.cos(xm[i, j])) - A
        
        fvec = np.concatenate([fv.flatten(), isum, jsum])

    elif nprob == 216:  # ConnBand
        for i in range(n - 2):
            fvec[i] = (x[i] + x[i + 1] + x[n - 1])**2
        fvec[n - 2] = x[0] - x[1]
        fvec[n - 1] = x[n - 2] - x[n - 1]

    elif nprob == 217:  # FREUROTH
        for i in range(1, n):
            fvec[i - 1] = -13 + x[i - 1] + ((5 - x[i]) * x[i] - 2) * x[i]
            fvec[n - 1 + i - 1] = -29 + x[i - 1] + ((1 + x[i]) * x[i] - 14) * x[i]

    elif nprob == 209:  # FREURONE
        for i in range(1, n):
            fvec[i - 1] = -13 + x[i - 1] + ((5 - x[i]) * x[i] - 2) * x[i]
            fvec[n - 1 + i - 1] = -29 + x[i - 1] + ((1 + x[i]) * x[i] - 14) * x[i]
 
    if len(fvec) != m:
        print('Wrong m')
        print(n)
        print(m)
        print(len(fvec))


    fvec = np.array(fvec)
    list = fvec.tolist()
    with open(file_name, 'a+') as file:
        file.write(str(list))
        file.write('\n')


    return fvec




