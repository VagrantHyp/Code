#This files contains methods that acts on the parsed data
#The purpose is to avoid parsing the same data over and over (which wastes a lot of time)
#It also creates intermediate data that allows more flexible changes

import numpy as np
import matplotlib.pylab as plt

#Global Ree parameter
Nlpf = 10009
M = 100; N = 100
Ndata = int(M*N)
Ncols = 9