import numpy as np
import matplotlib.pylab as pl
import ot
import csv
from numpy import genfromtxt

N = 2
d = 2

directory = '/home/hinouzhe/pCloudDrive/R/Capitulo_Libro_Cytometry/'

x1 = np.genfromtxt(directory + 'data1.csv', delimiter=',', skip_header = 1)
x2 = np.genfromtxt(directory + 'data2.csv', delimiter=',', skip_header = 1)

measures_locations = [x1, x2]
measures_weights = [ot.unif(x1.shape[0]), ot.unif(x2.shape[0])]

k = 1500  # number of Diracs of the barycenter
X_init = x1[np.random.random_integers(0, x1.shape[0] - 1, k)]
b = np.ones((k,)) / k  # weights of the barycenter (it will not be optimized, only the locations are optimized)

X = ot.lp.free_support_barycenter(measures_locations, measures_weights, X_init, b, numItermax = 150, verbose = True)

np.savetxt(directory + 'wasser_bar.csv', X, delimiter=",")

x1 = np.genfromtxt(directory + 'data1C1.csv', delimiter=',', skip_header = 1)
x2 = np.genfromtxt(directory + 'data2C1.csv', delimiter=',', skip_header = 1)

measures_locations = [x1, x2]
measures_weights = [ot.unif(x1.shape[0]), ot.unif(x2.shape[0])]

k = 931  # number of Diracs of the barycenter
X_init = x1[np.random.random_integers(0, x1.shape[0] - 1, k)]
b = np.ones((k,)) / k  # weights of the barycenter (it will not be optimized, only the locations are optimized)

X = ot.lp.free_support_barycenter(measures_locations, measures_weights, X_init, b, numItermax = 150, verbose = True)

np.savetxt(directory + 'wasser_bar_C1.csv2', X, delimiter=",")

x1 = np.genfromtxt(directory + 'data1C2.csv', delimiter=',', skip_header = 1)
x2 = np.genfromtxt(directory + 'data2C2.csv', delimiter=',', skip_header = 1)

measures_locations = [x1, x2]
measures_weights = [ot.unif(x1.shape[0]), ot.unif(x2.shape[0])]

k = 482  # number of Diracs of the barycenter
X_init = x1[np.random.random_integers(0, x1.shape[0] - 1, k)]
b = np.ones((k,)) / k  # weights of the barycenter (it will not be optimized, only the locations are optimized)

X = ot.lp.free_support_barycenter(measures_locations, measures_weights, X_init, b, numItermax = 150, verbose = True)

np.savetxt(directory + 'wasser_bar_C2.csv', X, delimiter=",")

x1 = np.genfromtxt(directory + 'data1C3.csv', delimiter=',', skip_header = 1)
x2 = np.genfromtxt(directory + 'data2C3.csv', delimiter=',', skip_header = 1)

measures_locations = [x1, x2]
measures_weights = [ot.unif(x1.shape[0]), ot.unif(x2.shape[0])]

k = 80  # number of Diracs of the barycenter
X_init = x1[np.random.random_integers(0, x1.shape[0] - 1, k)]
b = np.ones((k,)) / k  # weights of the barycenter (it will not be optimized, only the locations are optimized)

X = ot.lp.free_support_barycenter(measures_locations, measures_weights, X_init, b, numItermax = 150, verbose = True)

np.savetxt(directory + 'wasser_bar_C3.csv', X, delimiter=",")

x1 = np.genfromtxt(directory + 'data1C4.csv', delimiter=',', skip_header = 1)
x2 = np.genfromtxt(directory + 'data2C4.csv', delimiter=',', skip_header = 1)

measures_locations = [x1, x2]
measures_weights = [ot.unif(x1.shape[0]), ot.unif(x2.shape[0])]

k = 7  # number of Diracs of the barycenter
X_init = x1[np.random.random_integers(0, x1.shape[0] - 1, k)]
b = np.ones((k,)) / k  # weights of the barycenter (it will not be optimized, only the locations are optimized)

X = ot.lp.free_support_barycenter(measures_locations, measures_weights, X_init, b, numItermax = 150, verbose = True)

np.savetxt(directory + 'wasser_bar_C4.csv', X, delimiter=",")


