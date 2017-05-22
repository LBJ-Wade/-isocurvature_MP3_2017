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

colorlist = ['b','r','g','c','m', 'k','y']
colorInd = 0


basedir = '../chains/'
# CONFIGURATION -------------
chainfile = "../chains/planckdata/r1.txt"
paramfile = "../chains/planckdata/param"
xname = 'P_{II}^1'
yname = 'P_{RI}^1'
options = ""


chainfolder = "chains/CDI_2/"
# chainfolder = "/u/zequnl/Installs/MontePython/chains/pa/"
for filename in os.listdir(chainfolder):
    if '.paramnames' in filename:
        paramfile = os.path.join(chainfolder, filename)


# ---------------------------
params = np.array(ascii.read(paramfile, delimiter="\t", format="no_header")['col1'])

data_all = None

for filename in os.listdir(chainfolder):
    if filename.endswith(".txt"):

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

def denplot( list_data, ax, name="data", \
    lower=0.0, upper=0.25, acc=data['acceptance'], nbins=20, extend=False, \
    extent=0.1, cov=0.2, fmt="k--", mylabel="label" ):
    x = np.linspace(lower, upper, 300)


    if extend:

        bools = list_data < extent
        new_list_data = np.hstack(  (list_data,-list_data) )
        new_weights = np.hstack(  (data['acceptance'], (data['acceptance']) ) )

        density = gaussian_kde(new_list_data)

    else:
        density = gaussian_kde( list_data )
    density.covariance_factor = lambda : cov
    density._compute_covariance()
    ax.plot( x, density(x) / np.max(density(x)), fmt, label=mylabel )

    counts, bins = np.histogram( list_data, bins=x, weights=acc, density=True )
    #ax.plot( x[:-1], counts, "r." )
    ax.get_yaxis().set_ticks([])
   # ax.set_ylim( 0.0, counts.max() )
    ax.set_xlim( lower, upper )
    ax.set_xlabel( name )




bp1 = beta_iso1[:]
bp2 = beta_iso2[:]
ap1 = alpha[:]
npi = nII[:]
pacc = data['acceptance'][:]

fig = plt.figure(figsize=(12,18))
# ax1 = fig.add_subplot(141)
# ax2 = fig.add_subplot(142)
# ax3 = fig.add_subplot(143)
# ax4 = fig.add_subplot(144)

# c = 'k-'
# denplot( bp1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=pacc, extend=True, fmt=c )
# denplot( bp2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=pacc, extend=True, fmt=c)
# denplot( ap1,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=pacc, fmt=c)
# denplot( npi,     ax4, r"$n_{II}$", -1.0, 2.8, fmt=c, acc=pacc, mylabel="Planck"  )


## NOW ACT PART

basedir = '../chains/'
# CONFIGURATION ------------------------------------------------
chainfile = ""
paramfile = ""
xname = 'P_{II}^1'
yname = 'P_{RI}^1'
options = ""


# chainfolder = "chains/fF/"
# mylabel="PIXIE lowl, S4"

chainfolder = basedir + 'planck_zero_iso'
mylabel="planck lowl, planck forecast for highl"

# chainfolder = "chains/fC/"
# mylabel="Planck lowl, S4"

# chainfolder = "chains/fE/"
# mylabel="PIXIE lowl, Planck highl"


outname = mylabel

for filename in os.listdir(chainfolder):
    if '.paramnames' in filename:
        paramfile = os.path.join(chainfolder, filename)

print(paramfile)


# ---------------------------, delimiter="\t", format="no_header"
params = np.array(ascii.read(paramfile,delimiter="\t", format="no_header"))['col1']
print(params)

data_all = None

print(chainfolder)
for filename in os.listdir(chainfolder):
    if filename.startswith("201") and filename.endswith(".txt"):

        chainfile = os.path.join(chainfolder, filename)
        print(chainfile)
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
            print(len(data), len(data_all))

data = data_all


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



ax1 = fig.add_subplot(541)
ax2 = fig.add_subplot(542)
ax3 = fig.add_subplot(543)
ax4 = fig.add_subplot(544)

c = 'k-'
denplot( bp1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=pacc, extend=True, fmt=c )
denplot( bp2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=pacc, extend=True, fmt=c)
denplot( ap1,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=pacc, fmt=c)
denplot( npi,     ax4, r"$n_{II}$", -1.0, 2.8, fmt=c, acc=pacc, mylabel="Planck"  )


c = colorlist[colorInd] + '-'
colorInd += 1
denplot( beta_iso1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=data['acceptance'], extend=True, fmt=c )
denplot( beta_iso2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=data['acceptance'], extend=True, fmt=c)
denplot( alpha,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=data['acceptance'], fmt=c)
denplot( nII,     ax4, r"$n_{II}$", -1.0, 2.8, acc=data['acceptance'], fmt=c, mylabel=mylabel  )

plt.legend()


## NOW ACT PART
# CONFIGURATION ------------------------------------------------
chainfile = ""
paramfile = ""
xname = 'P_{II}^1'
yname = 'P_{RI}^1'
options = ""


# chainfolder = "chains/fF/"
# mylabel="PIXIE lowl, S4"

chainfolder = basedir + 'fB'
mylabel="planck lowl, planck+pol forecast for highl"

# chainfolder = "chains/fC/"
# mylabel="Planck lowl, S4"

# chainfolder = "chains/fE/"
# mylabel="PIXIE lowl, Planck highl"


outname = mylabel

for filename in os.listdir(chainfolder):
    if '.paramnames' in filename:
        paramfile = os.path.join(chainfolder, filename)

print(paramfile)


# ---------------------------, delimiter="\t", format="no_header"
params = np.array(ascii.read(paramfile,delimiter="\t", format="no_header"))['col1']
print(params)

data_all = None

for filename in os.listdir(chainfolder):
    if filename.startswith("201") and filename.endswith(".txt"):

        chainfile = os.path.join(chainfolder, filename)
        print(chainfile)
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
            print(len(data), len(data_all))

data = data_all


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


ax1 = fig.add_subplot(545)
ax2 = fig.add_subplot(546)
ax3 = fig.add_subplot(547)
ax4 = fig.add_subplot(548)

c = 'k-'
denplot( bp1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=pacc, extend=True, fmt=c )
denplot( bp2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=pacc, extend=True, fmt=c)
denplot( ap1,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=pacc, fmt=c)
denplot( npi,     ax4, r"$n_{II}$", -1.0, 2.8, fmt=c, acc=pacc, mylabel="Planck"  )

c = colorlist[colorInd] + '-'
colorInd += 1
denplot( beta_iso1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=data['acceptance'], extend=True, fmt=c )
denplot( beta_iso2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=data['acceptance'], extend=True, fmt=c)
denplot( alpha,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=data['acceptance'], fmt=c)
denplot( nII,     ax4, r"$n_{II}$", -1.0, 2.8, acc=data['acceptance'], fmt=c, mylabel=mylabel  )

plt.legend()


## NOW ACT PART
# CONFIGURATION ------------------------------------------------
chainfile = ""
paramfile = ""
xname = 'P_{II}^1'
yname = 'P_{RI}^1'
options = ""


# chainfolder = "chains/fF/"
# mylabel="PIXIE lowl, S4"

# chainfolder = "chains/oldF/fA/"
# mylabel="planck lowl, planck forecast for highl"

chainfolder = basedir + 'fC'
mylabel="Planck lowl, S4"

# chainfolder = "chains/fE/"
# mylabel="PIXIE lowl, Planck highl"


outname = mylabel

for filename in os.listdir(chainfolder):
    if '.paramnames' in filename:
        paramfile = os.path.join(chainfolder, filename)

print(paramfile)


# ---------------------------, delimiter="\t", format="no_header"
params = np.array(ascii.read(paramfile,delimiter="\t", format="no_header"))['col1']
print(params)

data_all = None

for filename in os.listdir(chainfolder):
    if filename.startswith("201") and filename.endswith(".txt"):

        chainfile = os.path.join(chainfolder, filename)
        print(chainfile)
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
            print(len(data), len(data_all))

data = data_all


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


ax1 = fig.add_subplot(549)
ax2 = fig.add_subplot(5,4,10)
ax3 = fig.add_subplot(5,4,11)
ax4 = fig.add_subplot(5,4,12)

c = 'k-'
denplot( bp1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=pacc, extend=True, fmt=c )
denplot( bp2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=pacc, extend=True, fmt=c)
denplot( ap1,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=pacc, fmt=c)
denplot( npi,     ax4, r"$n_{II}$", -1.0, 2.8, fmt=c, acc=pacc, mylabel="Planck"  )

c = colorlist[colorInd] + '-'
colorInd += 1
denplot( beta_iso1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=data['acceptance'], extend=True, fmt=c )
denplot( beta_iso2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=data['acceptance'], extend=True, fmt=c)
denplot( alpha,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=data['acceptance'], fmt=c)
denplot( nII,     ax4, r"$n_{II}$", -1.0, 2.8, acc=data['acceptance'], fmt=c, mylabel=mylabel  )

plt.legend()


## NOW ACT PART
# CONFIGURATION ------------------------------------------------
chainfile = ""
paramfile = ""
xname = 'P_{II}^1'
yname = 'P_{RI}^1'
options = ""


# chainfolder = "chains/fF/"
# mylabel="PIXIE lowl, S4"

# chainfolder = "chains/oldF/fA/"
# mylabel="planck lowl, planck forecast for highl"

# chainfolder = "chains/fC/"
# mylabel="Planck lowl, S4"

chainfolder = basedir + 'fE'
mylabel="PIXIE lowl, Planck highl"


outname = mylabel

for filename in os.listdir(chainfolder):
    if '.paramnames' in filename:
        paramfile = os.path.join(chainfolder, filename)

print(paramfile)


# ---------------------------, delimiter="\t", format="no_header"
params = np.array(ascii.read(paramfile,delimiter="\t", format="no_header"))['col1']
print(params)

data_all = None

for filename in os.listdir(chainfolder):
    if filename.startswith("201") and filename.endswith(".txt"):

        chainfile = os.path.join(chainfolder, filename)
        print(chainfile)
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
            print(len(data), len(data_all))

data = data_all


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

ax1 = fig.add_subplot(5,4,13)
ax2 = fig.add_subplot(5,4,14)
ax3 = fig.add_subplot(5,4,15)
ax4 = fig.add_subplot(5,4,16)

c = 'k-'
denplot( bp1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=pacc, extend=True, fmt=c )
denplot( bp2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=pacc, extend=True, fmt=c)
denplot( ap1,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=pacc, fmt=c)
denplot( npi,     ax4, r"$n_{II}$", -1.0, 2.8, fmt=c, acc=pacc, mylabel="Planck"  )


c = colorlist[colorInd] + '-'
colorInd += 1
denplot( beta_iso1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=data['acceptance'], extend=True, fmt=c )
denplot( beta_iso2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=data['acceptance'], extend=True, fmt=c)
denplot( alpha,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=data['acceptance'], fmt=c)
denplot( nII,     ax4, r"$n_{II}$", -1.0, 2.8, acc=data['acceptance'], fmt=c, mylabel=mylabel  )

plt.legend()



## NOW ACT PART
# CONFIGURATION ------------------------------------------------
chainfile = ""
paramfile = ""
xname = 'P_{II}^1'
yname = 'P_{RI}^1'
options = ""


chainfolder = basedir + 'fF'
mylabel="PIXIE lowl, S4"

# chainfolder = "chains/oldF/fA/"
# mylabel="planck lowl, planck forecast for highl"

# chainfolder = "chains/fC/"
# mylabel="Planck lowl, S4"

# chainfolder = "chains/fE/"
# mylabel="PIXIE lowl, Planck highl"


outname = mylabel

for filename in os.listdir(chainfolder):
    if '.paramnames' in filename:
        paramfile = os.path.join(chainfolder, filename)

print(paramfile)


# ---------------------------, delimiter="\t", format="no_header"
params = np.array(ascii.read(paramfile,delimiter="\t", format="no_header"))['col1']
print(params)

data_all = None

for filename in os.listdir(chainfolder):
    if filename.startswith("201") and filename.endswith(".txt"):

        chainfile = os.path.join(chainfolder, filename)
        print(chainfile)
        data = ascii.read(chainfile, delimiter="\s")

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

ax1 = fig.add_subplot(5,4,17)
ax2 = fig.add_subplot(5,4,18)
ax3 = fig.add_subplot(5,4,19)
ax4 = fig.add_subplot(5,4,20
    )

c = 'k-'
denplot( bp1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=pacc, extend=True, fmt=c )
denplot( bp2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=pacc, extend=True, fmt=c)
denplot( ap1,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=pacc, fmt=c)
denplot( npi,     ax4, r"$n_{II}$", -1.0, 2.8, fmt=c, acc=pacc, mylabel="Planck"  )

c = colorlist[colorInd] + '-'
colorInd += 1
denplot( beta_iso1, ax1, r"$\beta_{iso}(k_{low})$", 0.0, 0.1, acc=data['acceptance'], extend=True, fmt=c )
denplot( beta_iso2, ax2, r"$\beta_{iso}(k_{high})$", 0.0, 0.8, acc=data['acceptance'], extend=True, fmt=c)
denplot( alpha,     ax3, r"$\cos \Delta$", -0.5, 0.5, acc=data['acceptance'], fmt=c)
denplot( nII,     ax4, r"$n_{II}$", -1.0, 2.8, acc=data['acceptance'], fmt=c, mylabel=mylabel  )

plt.legend()


plt.tight_layout()
plt.savefig("../../figures/all_derived_forecast.pdf")
plt.show()
