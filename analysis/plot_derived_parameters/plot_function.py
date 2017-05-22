import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from scipy.stats import gaussian_kde
from pprint import pprint

import sys
import os
from astropy.io import ascii
from astropy.table import vstack

# THIS FILE: UTILITY FUNCTIONS FOR PLOTTING!

def loadChainFolder(chainfolder):

    for filename in os.listdir(chainfolder):
        if '.paramnames' in filename:
            paramfile = os.path.join(chainfolder, filename)

    # print(paramfile)


    params = np.array(ascii.read(paramfile,delimiter="\t", format="no_header"))['col1']
    # print(params)

    data_all = None

    # print(chainfolder)
    for filename in os.listdir(chainfolder):
        if filename.endswith(".txt"):

            chainfile = os.path.join(chainfolder, filename)
            # print(chainfile)
            data = (ascii.read(chainfile, delimiter="\s"))[100:]

            # set up column names (read in from param file)
            data['col1'].name = 'acceptance'
            data['col2'].name = 'likelihood'
            for i in range(3,len(params)+3):
                data['col' + str(i)].name = params[i-3]


            if data_all == None:
                data_all = data
            else:
                data_all =  vstack( [data_all, data] )
                # print(len(data), len(data_all))

    return data_all


def repeatRows(data, acceptance):
    newData = data[:]

    # for each row in data, add rows based on acceptance number
    for row, acc in zip(data, acceptance):
        for i in range(acc-1):
            newData.append( row )

    return newData

def denplot( list_data, ax, acc, name="data", \
    lower=0.0, upper=0.25, nbins=20, extend=False, \
    extent=0.1, cov=0.2, fmt="k-", mylabel="label" ):

    # print("repeating")
    # list_data = np.array(list_data).tolist()
    # list_data = repeatRows(list_data, acc)
    # list_data = np.array(list_data)

    x = np.linspace(lower, upper, 300)

    # new_weights = data['acceptance']
    if extend:
        new_list_data = np.hstack(  (list_data,-list_data) )
        density = gaussian_kde(new_list_data)
    else:
        density = gaussian_kde( list_data )
    density.covariance_factor = lambda : cov
    density._compute_covariance()
    ax.plot( x, density(x) / np.max(density(x)), fmt, label=mylabel )

    # counts, bins = np.histogram( list_data, bins=x, weights=new_weights, density=True )
    #ax.plot( x[:-1], counts, "r." )
    ax.get_yaxis().set_ticks([])
   # ax.set_ylim( 0.0, counts.max() )
    ax.set_xlim( lower, upper )
    ax.set_xlabel( name )


def plotRow(data, ax1, ax2, ax3, ax4, c, mylabel):

    prr1 = data['P_{RR}^1']; pii1 = data['P_{II}^1'];  pri1 = data['P_{RI}^1'];
    prr2 = data['P_{RR}^2']; pii2 = data['P_{II}^2'];  pri2 = pri1 * np.sqrt(pii2 * prr2 / (pii1 * prr1))
    # make density plot
    # sc(x,y)
    beta_iso1 = pii1 / (prr1 + pii1)
    beta_iso2 = pii2 / (prr2 + pii2)
    alpha = pri1 / np.sqrt( pii1 * prr1 )
    # \frac{\log( P_{AB}^2 / P_{AB}^1 )}{\log ( k_2 / k_1 )
    k1 = 0.002 # Mpc^{-1}
    k2 = 0.1 # Mpc^{-1}
    nRR = np.log(prr2/prr1) / np.log(k2/k1)
    nRI = np.log(pri2/pri1) / np.log(k2/k1)
    nII = np.log(pii2/pii1) / np.log(k2/k1)

    denplot( beta_iso1, ax1, data['acceptance'], r"$\beta_{iso}(k_{low})$", 0.0, 0.1, extend=True, fmt=c )
    denplot( beta_iso2, ax2, data['acceptance'], r"$\beta_{iso}(k_{high})$", 0.0, 0.8, extend=True, fmt=c)
    denplot( alpha,     ax3, data['acceptance'], r"$\cos \Delta$", -0.5, 0.5, fmt=c)
    denplot( nII,     ax4, data['acceptance'], r"$n_{II}$", -1.0, 2.8, fmt=c, mylabel=mylabel  )
    ax4.legend()
