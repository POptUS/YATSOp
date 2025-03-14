import numpy as np




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
file_path2 = 'dfof.txt'
data2 = read_matdata(file_path2)


print("calldfofuns file")
#print(array)


for i in range(40):
    d1 = np.array(data1[i])
    d2 = np.array(data2[i])
    diff = np.sum((d1-d2)**2)  /  np.sum((d1)**2)
    print(i, diff)
    # if np.isnan(diff):
    #     print(i, np.sum((d1-d2)**2))
    #     print(i, np.sum((d1)**2))


file_path = 'pydfomidf.txt' 
data1 = read_pydata(file_path)
file_path2 = 'dfomidf.txt'
data2 = read_matdata(file_path2)


print("calldmidfuns file")
#print(array)


for i in range(40):
    d1 = np.array(data1[i])
    d2 = np.array(data2[i])
    diff = np.sum((d1-d2)**2)  /  np.sum((d1)**2)
    print(i, diff)
    # if np.isnan(diff):
    #     print(i, np.sum((d1-d2)**2))
    #     print(i, np.sum((d1)**2))
