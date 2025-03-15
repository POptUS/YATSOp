import numpy as np
import time
from calfun_sample import calfun_sample
from dfoxsnew import dfoxsnew

# This assumes that you have first run regressiontest.m in ../m/

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
nrows = Var.shape[0]
probtype = "smooth"
probspecs = {'trunc': 10**16}

# Print headers
print('Problem     n    m            f0     h')


namestr = []
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

    # Compute the objective function (calculation)
    y, fvec = calfun_sample(X0, probspecs, probtype)

    fvec = np.array(fvec)
    list = fvec.tolist()

    with open("pydfof.txt", 'a+') as file:
        file.write(str(list))
        file.write('\n')

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

# Set flags and parameters
nrows = Var.shape[0]
probtype = "smooth"
probspecs = {'trunc': 10**16}

# Print headers
print('Problem     n    m            f0     h')


namestr = []
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

    # Compute the objective function (calculation)
    y, fvec = calfun_sample(X0, probspecs, probtype)

    fvec = np.array(fvec)
    list = fvec.tolist()

    with open("pydfomidf.txt", 'a+') as file:
        file.write(str(list))
        file.write('\n')

   

def read_pydata(file_path):
    data = []
    lines = []
    with open(file_path, 'r') as file:
        for line in file:
            lines.append(line)
            
            line = line.strip().strip('[]')
            
            parsed_line = [float(item) for item in line.split(',')]
            data.append(parsed_line)

    # array = np.array(data)
    # return array
    return data

def read_matdata(file_path):
    data = []
    lines = []
    with open(file_path, 'r') as file:
        for line in file:
            lines.append(line)
 
            line = line.strip().strip('[]')

            parsed_line = [float(item) for item in line.split(';')]
            data.append(parsed_line)


    # array = np.array(data)
    # return array

    return data

#################################################

file_path = 'pydfof.txt'  
data1 = read_pydata(file_path)
file_path2 = '../m/dfof.txt'
data2 = read_matdata(file_path2)


print("calldfofuns file")
#print(array)

diff = np.zeros(40)
for i in range(40):
    d1 = np.array(data1[i])
    d2 = np.array(data2[i])
    diff[i] = np.sqrt(np.sum((d1-d2)**2))  /  np.sqrt(np.sum((d1)**2))
    #print(i, diff[i])
    # if np.isnan(diff[i]):
    #     print(i, np.sum((d1-d2)**2))
    #     print(i, np.sum((d1)**2))

print("max=",np.max(diff)," median=",np.median(diff))

file_path = 'pydfomidf.txt' 
data1 = read_pydata(file_path)
file_path2 = '../m/dfomidf.txt'
data2 = read_matdata(file_path2)


print("calldmidfuns file")
#print(array)


diff = np.zeros(40)
for i in range(40):
    d1 = np.array(data1[i])
    d2 = np.array(data2[i])
    diff[i] = np.sqrt(np.sum((d1-d2)**2))  /  np.sqrt(np.sum((d1)**2))
    #print(i, diff[i])
    # if np.isnan(diff[i]):
    #     print(i, np.sum((d1-d2)**2))
    #     print(i, np.sum((d1)**2))

print("max=",np.max(diff)," median=",np.median(diff))
