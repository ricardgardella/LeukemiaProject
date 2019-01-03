import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from sklearn.cross_decomposition import PLSRegression


def plsr(data, response, ncomp, validation):
    nrows = data.shape[0]
    error_comp = []
    components = list(range(1, ncomp + 1))
    for c in components:
        pls1 = PLSRegression(n_components = c)
        errors = []
        for row in range(nrows):
            rows = [i for i in range(0, nrows) if i != row]
            pls1.fit(data.iloc[rows], response[:row] + response[row+1:])
            value = pls1.predict(data.iloc[[row]])
            error = (value[0][0] - response[row])**2
            #threshold = 0.5
            #prediction = 1 if value[0][0] > threshold else 0
            #error = 0 if prediction == response[row] else 1
            errors.append(error)
        error_comp.append(np.sqrt(np.mean(errors)))
    print(error_comp)
    plt.plot(components, error_comp)
    plt.show()

train = pd.read_csv("data/data_set_ALL_AML_train.csv",sep=";")
test = pd.read_csv("data/data_set_ALL_AML_independent.csv",sep=";")

# pre-process train
train = train.T
rows = [str(i) for i in range(1, 39)]
train = train[train.index.isin(rows)]
response = {1: 'ALL', 2: 'ALL', 3: 'ALL', 4: 'ALL', 5: 'ALL', 6: 'ALL', 7: 'ALL', 8: 'ALL', 9: 'ALL', 
    10: 'ALL', 11: 'ALL', 12: 'ALL', 13: 'ALL', 14: 'ALL', 15: 'ALL', 16: 'ALL', 17: 'ALL', 18: 'ALL',
    19: 'ALL', 20: 'ALL', 21: 'ALL', 22: 'ALL', 23: 'ALL', 24: 'ALL', 25: 'ALL', 26: 'ALL', 27: 'ALL',
    28: 'AML', 29: 'AML', 30: 'AML', 31: 'AML', 32: 'AML', 33: 'AML', 34: 'AML', 35: 'AML', 36: 'AML',
    37: 'AML', 38: 'AML'}

response_ordered = [1 if response[int(indx)] == 'ALL' else 0 for indx in list(train.index)]
train = train.assign(response = response_ordered)

plsr(train, response_ordered, 10, 'loo')

# pre-process test
test = test.T
rows_test = [str(i) for i in range(39, 73)]
test = test[test.index.isin(rows_test)]
response_test = {39: 'ALL', 40: 'ALL', 41: 'ALL', 42: 'ALL', 43: 'ALL', 44: 'ALL', 45: 'ALL', 46: 'ALL',
    47: 'ALL', 48: 'ALL', 49: 'ALL', 50: 'AML', 51: 'AML', 52: 'AML', 53: 'AML', 54: 'AML', 55: 'ALL',  
    56: 'ALL', 57: 'AML', 58: 'AML', 59: 'ALL', 60: 'AML', 61: 'AML', 62: 'AML', 63: 'AML', 64: 'AML',
    65: 'AML', 66: 'AML', 67: 'ALL', 68: 'ALL', 69: 'ALL', 70: 'ALL', 71: 'ALL', 72: 'ALL'}
response_test_ordered = [1 if response_test[int(indx)] == 'ALL' else 0 for indx in list(test.index)]
test = test.assign(response = response_test_ordered)

