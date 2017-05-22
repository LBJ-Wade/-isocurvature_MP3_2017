import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from scipy.stats import gaussian_kde
from pprint import pprint

import sys
import os
from astropy.io import ascii
from astropy.table import vstack

# chainfile = "/Users/zequnl/Installs/montepython_public/chains/example/2016-10-18_10000__1.txt"

# CONFIGURATION -------------
# chainfile = "chains/CDI_2/2016-11-02_1000000__1.txt"
# paramfile = "chains/CDI_2/2016-11-02_1000000_.paramnames"

chainfile = "chains/CDI_2/2016-11-02_1000000__1.txt"
paramfile = "chains/CDI_2/2016-11-02_1000000_.paramnames"

xname = 'P_{II}^1'
yname = 'P_{RI}^1'
options = ""

chainfolder = "chains/CDI_2/"

if len(sys.argv) >= 3:
    chainfile = sys.argv[1]
    paramfile = sys.argv[2]
    print(paramfile)

if len(sys.argv) >= 5:
    xname = sys.argv[3]
    yname = sys.argv[4]
    options = sys.argv[5:]
elif len(sys.argv) == 2:
    if sys.argv[1] == "info":
        params = np.array(ascii.read('params', delimiter="\t", format="no_header")['col1']).tolist()
        print(params)
        sys.exit(0)

# ---------------------------



params = np.array(ascii.read(paramfile, delimiter="\t", format="no_header")['col1'])


data_all = None

for filename in os.listdir(chainfolder):
    if filename.startswith("201") and filename.endswith(".txt"):

        chainfile = os.path.join(chainfolder, filename)
        print(chainfile)
        data = (ascii.read(chainfile, delimiter="\s"))[300:]

        # set up column names (read in from param file)
        data['col1'].name = 'acceptance'
        data['col2'].name = 'likelihood'
        for i in range(3,len(params)+3):
            data['col' + str(i)].name = params[i-3]


        if data_all == None:
            data_all = data
        else:
            data_all =  vstack( [data_all, data] )
            print(len(data), len(data_all))

data = data_all

print(len(data), "rows")



x = data[xname]
y = data[yname]
t = np.array( range(len(x)))

# we look for the -s option, and then find the number afterwards.
# that's where we start
if "-s" in options:
    s = int(sys.argv[sys.argv.index("-s")+1])
    x = x[s:]
    y = y[s:]


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

def denplot( list_data, ax, name="data", lower=0.0, upper=0.25, \
    nbins=20, extend=False, extent=0.1, cov=0.2 ):
    """
    plot a smoothed histogram
    """
    x = np.linspace(lower, upper, 150)

    if extend:

        bools = list_data < extent
        new_list_data = np.hstack(  (list_data,-list_data) )
        new_weights = np.hstack(  (data['acceptance'], (data['acceptance']) ) )

        density = gaussian_kde(new_list_data)

    else:
        density = gaussian_kde( list_data )
    density.covariance_factor = lambda : cov
    density._compute_covariance()
    ax.plot( x, density(x), "k--" )

    counts, bins = np.histogram( list_data, bins=x, weights=data['acceptance'], density=True )
    #ax.plot( x[:-1], counts, "r." )
    ax.get_yaxis().set_ticks([])
   # ax.set_ylim( 0.0, counts.max() )
    ax.set_xlim( lower, upper )
    ax.set_xlabel( name )




fig = plt.figure(figsize=(12,3))
ax1 = fig.add_subplot(141)
ax2 = fig.add_subplot(142)
ax3 = fig.add_subplot(143)
ax4 = fig.add_subplot(144)

denplot( beta_iso1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.25, extend=True )
denplot( beta_iso2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, extend=True)
denplot( alpha,     ax3, r"$\cos \Delta$", -0.5, 0.5 )
denplot( nII,     ax4, r"$n_{II}$", -1.0, 2.8 )

plt.tight_layout()
plt.savefig("../../figures/beta_planck.pdf")
plt.show()

## TESTING
