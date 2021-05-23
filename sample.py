import pandas as pd
import numpy as np
import scipy.stats
import glob
import collections


# in order to reorder methylation samples
sample = pd.read_csv("sample.csv", sep="\t", header=None)


def remove(s):
    return float(s.replace('GR', '').replace('_V', ''))


a = list(map(remove, sample.values[:, 1]))
index = np.argsort(a)

# methylation
x = pd.read_csv("GSE161020_series_matrix.txt.gz", sep="\t", comment="!")
y = scipy.stats.zscore(x.values[:, 1:][:, index].astype(np.float64))
np.save("y_methyl",y)
X = y.T@y
np.save("X",X)

# gene expression
path = './*.count.txt.gz'
files = glob.glob(path)
files = sorted(files)
for i in range(len(files)):
    print(i)
    x = pd.read_csv(files[i], sep="\t", header=None)
    if i == 0:
        x_all = x.values[:, 1]
    else:
        x_all = np.vstack((x_all, x.values[:, 1]))

x_all = x_all.T
y = scipy.stats.zscore(x_all.astype(np.float64))
np.save("y_exp",y)
X0 = y.T@y
np.save("X0",X0)

# proteome
x = pd.read_csv("GR01,04,09,10,11,13,15,17,18,19.txt", sep="\t")
x1 = pd.read_csv("GR02,03,05,06,07.txt", sep="\t")
x.values[0:3, 1] = ["a", "b", "c"]
x1.values[0:3, 1] = ["a", "b", "c"]
x0 = pd.merge(x, x1, on='Unnamed: 1', how='outer')
cells = np.array([l[0] for l in x0.columns.str.split(".")])
GR = x0.values[0, ]
visit = x0.values[1, ]
np.save("x0",x0)

def remove2(s):
    return s.replace('_x', '').replace('_y', '')


cells = np.array(list(map(remove2, cells)))
labels_cells = [label for label in collections.Counter(cells).keys()]
labels_cells = sorted(np.array(labels_cells))
labels_cells = [labels_cells[i] for i in [0, 3]]
labels_GR = [label for label in collections.Counter(GR).keys()]
labels_GR = sorted(np.array(labels_GR))
labels_GR = labels_GR[:-2]
labels_visit = [label for label in collections.Counter(visit).keys()]
labels_visit = sorted(np.array(labels_visit))
labels_visit = labels_visit[:-2]
# use np.isin
visit[np.bitwise_not(np.isin(visit, labels_visit))] = "X"
cells[np.bitwise_not(np.isin(cells, labels_cells))] = "X"
GR[np.bitwise_not(np.isin(GR, labels_GR))] = "X"

#visit[[i not in labels_visit for i in visit]] = "X"
#cells[[i not in labels_cells for i in cells]] = "X"
#GR[[i not in labels_GR for i in GR]] = "X"

id_WBC = None
id_Plasma = None
for j in range(len(labels_GR)):
    for i in range(len(labels_visit)):        
        index = [x == labels_visit[i] and y == labels_GR[j] for (x, y) in zip(visit, GR)]
        index_WBC = [x and y == "WBC" for (x, y) in zip(index, cells)]
        index_Plasma = [x and y == "Plasma" for (x, y) in zip(index, cells)]
        id_WBC = np.hstack((id_WBC, np.arange(183)[index_WBC][0]))
        id_Plasma = np.hstack((id_Plasma, np.arange(183)[index_Plasma][0]))
       # print(i,j,np.arange(183)[index_WBC],np.arange(183)[index_Plasma])

id_WBC = id_WBC[1:]
id_Plasma = id_Plasma[1:]


x0_WBC = x0.values[:, list(map(int, id_WBC))]
x0_Plasma = x0.values[:, list(map(int, id_Plasma))]

x0_WBC = x0_WBC[4:, :].astype(np.float64)
x0_Plasma = x0_Plasma[4:, :].astype(np.float64)

x0_WBC = np.nan_to_num(x0_WBC)
x0_Plasma = np.nan_to_num(x0_Plasma)

y = scipy.stats.zscore(x0_WBC)
np.save("y_WBC",y)
X1_WBC = y.T@y

y = scipy.stats.zscore(x0_Plasma)
np.save("y_Plasma",y)
X1_Plasma = y.T@y

np.save("X1_WBC",X1_WBC)
np.save("X1_Plasma",X1_Plasma)

#---- 
X = np.load("X.npy")
X0 = np.load("X0.npy")
X1_WBC = np.load("X1_WBC.npy")
X1_Plasma = np.load("X1_Plasma.npy")
X = np.reshape(X,[15,5,15,5])
X = X.T
X = X/np.mean(X)
X0 = np.reshape(X0,[15,5,15,5])
X0 = X0.T
X0 = X0/np.mean(X0)
X1_WBC = np.reshape(X1_WBC,[15,5,15,5])
X1_WBC = X1_WBC.T
X1_WBC = X1_WBC/np.mean(X1_WBC)
X1_Plasma = np.reshape(X1_Plasma,[15,5,15,5])
X1_Plasma = X1_Plasma.T
X1_Plasma = X1_Plasma/np.mean(X1_Plasma)

import tensorly as tl

Z = tl.tensor(np.zeros(75*75*4).reshape(5,15,5,15,4))
Z[:,:,:,:,0]=X
Z[:,:,:,:,1]=X0
Z[:,:,:,:,2]=X1_WBC
Z[:,:,:,:,3]=X1_Plasma
from tensorly.decomposition import tucker
core, factors = tucker(Z, rank=[5,15,5,15,4],n_iter_max=1,init='svd')

np.save("core",core)
np.save("factors",factors)

#-----
import statsmodels.stats.multitest

core = np.load("core.npy")
factors = np.load("factors.npy",allow_pickle=True)
u= np.outer(factors[0][:,1],factors[1][:,0])
u.reshape(75)


#methylation
y = np.load("y_methyl.npy")
P = scipy.stats.chi2.sf(scipy.stats.zscore(scipy.stats.zscore(y)@u.T.reshape(75))**2,1)
P0 = np.zeros(len(P))
P0.fill(0.01)
collections.Counter(statsmodels.stats.multitest.multipletests(P,0.01,method="fdr_bh")[0])

#Counter({False: 685505, True: 2077})

y = pd.read_csv("GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz",sep=",",skiprows=7)
gene = y.values[:,14][np.isin(y.values[:,0],x.values[statsmodels.stats.multitest.multipletests(P,0.01,method="fdr_bh")[0],0])]

#gene expression
y=np.load("y_exp.npy")
P = scipy.stats.chi2.sf(scipy.stats.zscore(scipy.stats.zscore(y)@u.T.reshape(75))**2,1)
P0 = np.zeros(len(P))
P0.fill(0.01)
collections.Counter(statsmodels.stats.multitest.multipletests(P,0.01,method="fdr_bh")[0])

#Counter({False: 58296, True: 11})
x = pd.read_csv(files[0], sep="\t", header=None)
gene = x.values[statsmodels.stats.multitest.multipletests(P,0.01,method="fdr_bh")[0],0]

#Proteome 
#WBC
y=np.load("y_WBC.npy")
P = scipy.stats.chi2.sf(scipy.stats.zscore(scipy.stats.zscore(y)@u.T.reshape(75))**2,1)
P0 = np.zeros(len(P))
P0.fill(0.05)
collections.Counter(statsmodels.stats.multitest.multipletests(P,0.01,method="fdr_bh")[0])

#Counter({False: 1579, True: 5})

x0=np.load("x0.npy",allow_pickle=True)
gene = x0[4:,0][statsmodels.stats.multitest.multipletests(P,0.01,method="fdr_bh")[0]]

#Plasma
y=np.load("y_Plasma.npy")
P = scipy.stats.chi2.sf(scipy.stats.zscore(scipy.stats.zscore(y)@u.T.reshape(75))**2,1)
P0 = np.zeros(len(P))
P0.fill(0.05)
collections.Counter(statsmodels.stats.multitest.multipletests(P,0.01,method="fdr_bh")[0])

#Counter({False: 1579, True: 5})

gene = x0[4:,0][statsmodels.stats.multitest.multipletests(P,0.01,method="fdr_bh")[0]]
