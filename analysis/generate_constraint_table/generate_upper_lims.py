
# this file generates upper limit tables for the zero isocurvature model

import matplotlib

import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table, Column
from loadMontePython import load as loadMCMC
import glob

from scipy.optimize import brentq
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde

basedir = '../chains/'

def read_chain(folder):
    files = glob.glob(folder + "/*__*.txt")
    params = glob.glob(folder + "/*_.paramnames")

    datalist = []
    for f in files:
        datalist.append(  loadMCMC(f, params[0]) )

    # print(datalist)
    data = astropy.table.vstack( datalist )


    return data

def plotVariable95UpperLim(ax, data_input, variable, inputLabel, fill=True):
    N = len(data_input)
    Z = data_input[variable]
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)

    interpFunc = interp1d(X2, F2 - 0.95)
    upperLim95 = brentq(interpFunc,np.min(Z), np.max(Z) )
    # plt.plot(X2, F2) # plt.vlines(upperLim95, ymin=0.0, ymax=1.0)
    # plt.show()

    density = gaussian_kde(Z)
    xs = np.linspace(np.min(Z), np.max(Z),200)
    density.covariance_factor = lambda : .25
    density._compute_covariance()

    filled_xs = xs[xs < upperLim95]
    if fill:
        ax.fill_between(filled_xs,density(filled_xs) / np.max(density(filled_xs)),\
        alpha=0.2, label=inputLabel)#,color='b'
    else:
        ax.fill_between( [], [],\
        alpha=0.2, label=inputLabel)#,color='b'
    ax.plot(xs,density(xs) /
        np.max(density(xs)),markersize=0.0)
    ax.set_ylim(0.0,1.0)
    ax.set_xlim(0.0,np.max(Z))
    ax.set_xlabel( r'$' + variable + r'$')

    return upperLim95, inputLabel


def plotVariable95Symmetric(ax, data_input, variable, inputLabel, fill=True):
    N = len(data_input)
    Z = data_input[variable]
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)

    interpFunc = interp1d(X2, F2 - 0.95)
    mean = np.mean(Z)
    std = np.std(Z)
    upperLim95 = mean + 2 * std
    lowerLim95 = mean - 2 * std
    # plt.plot(X2, F2) # plt.vlines(upperLim95, ymin=0.0, ymax=1.0)
    # plt.show()

    density = gaussian_kde(Z)
    xs = np.linspace(np.min(Z), np.max(Z),200)
    density.covariance_factor = lambda : .25
    density._compute_covariance()

    filled_xs =  xs[ np.logical_and(xs > lowerLim95,xs < upperLim95) ]
    if fill:
        ax.fill_between(filled_xs,density(filled_xs) / np.max(density(filled_xs)),\
        alpha=0.2, label=inputLabel)#,color='b'
    else:
        ax.fill_between( [], [],\
        alpha=0.2, label=inputLabel)#,color='b'
    ax.plot(xs,density(xs) /
        np.max(density(xs)),markersize=0.0)
    # ax.set_ylim(0.0,1.0)
    # ax.set_xlim(0.0,np.max(Z))
    ax.set_xlabel( r'$' + variable + r'$')

    return 2*std, inputLabel

def plotVariable95(ax, data_input, variable, inputLabel, fill=True):
    if variable == 'P_{RI}^1':
        return plotVariable95Symmetric(ax, data_input, variable, inputLabel, fill)
    else:
        return plotVariable95UpperLim(ax, data_input, variable, inputLabel, fill)


plt.clf()
f, axarr = plt.subplots(3)

expLabelList = ['planck TEB lowl + planck TT', \
    'Planck TEB', 'Planck lowl + S4','PIXIE + Planck highl','PIXIE + S4']
expFolderList = ['planckdata', \
     'planck_zero_iso',\
     'planck_s4_zero_iso',\
     'pixie_planck_zero_iso',\
     'pixie_s4_zero_iso']

print("reading in data...")
dataFiles = [read_chain(basedir + folder) for folder in expFolderList]

print("plotting...")
construction_dict = {'P_{II}^1':[],'P_{II}^2':[],'P_{RI}^1':[] }
for axNum, var in enumerate(['P_{II}^1','P_{II}^2','P_{RI}^1']):

    for expLabel, dataFile in zip(expLabelList, dataFiles):
        upper, label = plotVariable95(axarr[axNum], \
            dataFile,\
            var, expLabel)
        construction_dict[var].append(upper)



# %% now generate the figure

# print(construction_dict)
t = Table()
t['experiments'] = expLabelList
t['P_{II}^1'] = construction_dict['P_{II}^1']
t['P_{II}^2'] = construction_dict['P_{II}^2']
plt.tight_layout()
plt.legend()
# plt.savefig('debug_get_95_percentiles.pdf')
plt.show()

print(t)

# now we need to write something for P_{RI}^1
