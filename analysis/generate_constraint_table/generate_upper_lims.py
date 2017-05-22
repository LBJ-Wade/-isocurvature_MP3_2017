


import matplotlib

import numpy as np import matplotlib.pyplot as plt import astropy from
loadMontePython import load as loadMCMC import glob

from scipy.optimize import brentq from scipy.interpolate import interp1d from
scipy.stats import gaussian_kde

basedir = '../chains/zero_isocurvature_model/'
#"/home/zequnl/Projects/isocurvature_2017/analysis/plot_triangle/nonzero/"

burnin = 1000 data1 = loadMCMC('../chains/planckdata/r1.txt',
'../chains/planckdata/param') data2 = loadMCMC('../chains/planckdata/r2.txt',
'../chains/planckdata/param') data = astropy.table.vstack( [data1[burnin:],
data2[burnin:]]) data_planck = data[:] weights_planck = data['acceptance'][:]

## A folder = basedir + 'fA/' files = glob.glob(folder + "*__*.txt") params =
## glob.glob(folder + "*_.paramnames")

datalist = [] for f in files: datalist.append(  loadMCMC(f, params[0]) )

data = astropy.table.vstack( datalist ) data_A = data[:] weights_act_A =
data['acceptance'][:]

## B folder = basedir + 'fB/' files = glob.glob(folder + "*__*.txt") params =
## glob.glob(folder + "*_.paramnames")

datalist = [] for f in files: datalist.append(  loadMCMC(f, params[0]) )

data = astropy.table.vstack( datalist ) data_B = data[:] weights_act_B =
data['acceptance'][:]

## C folder = basedir + 'fC/' files = glob.glob(folder + "*__*.txt") params =
## glob.glob(folder + "*_.paramnames")

datalist = [] for f in files: datalist.append(  loadMCMC(f, params[0]) )

data = astropy.table.vstack( datalist ) data_C = data[:] weights_act_C =
data['acceptance'][:]

## E folder = basedir + 'fE/' files = glob.glob(folder + "*__*.txt") params =
## glob.glob(folder + "*_.paramnames") datalist = [] for f in files:
## datalist.append(  loadMCMC(f, params[0]) ) data = astropy.table.vstack(
## datalist ) data_E = data[:] weights_act_E = data['acceptance'][:]

## F folder = basedir + 'fF/' files = glob.glob(folder + "*__*.txt") params =
## glob.glob(folder + "*_.paramnames") datalist = [] for f in files:
## datalist.append(  loadMCMC(f, params[0]) ) data = astropy.table.vstack(
## datalist ) data_F = data[:] weights_act_F = data['acceptance'][:]

data_A.colnames

def plotVariable95(ax, data_input, variable, inputLabel, fill=True): N =
len(data_input) Z = data_input[variable] X2 = np.sort(Z) F2 =
np.array(range(N))/float(N)

	interpFunc = interp1d(X2, F2 - 0.95) upperLim95 = brentq(interpFunc,
	np.min(Z), np.max(Z) ) # plt.plot(X2, F2) # plt.vlines(upperLim95, ymin=0.0,
	ymax=1.0) # plt.show()

	density = gaussian_kde(Z) xs = np.linspace(np.min(Z), np.max(Z),200)
	density.covariance_factor = lambda : .25 density._compute_covariance()

	filled_xs = xs[xs < upperLim95] if fill:
	ax.fill_between(filled_xs,density(filled_xs) / np.max(density(filled_xs)),\
	alpha=0.2, label=inputLabel)#,color='b' else: ax.fill_between( [], [],\
	alpha=0.2, label=inputLabel)#,color='b' ax.plot(xs,density(xs) /
	np.max(density(xs)),markersize=0.0) ax.set_ylim(0.0,1.0)
	ax.set_xlim(0.0,np.max(Z)) ax.set_xlabel( r'$' + variable + r'$')

plt.clf() f, axarr = plt.subplots(2)
for axNum, var in enumerate(['P_{II}^1','P_{II}^2']):

    plotVariable95(axarr[axNum],data_A,var, 'planck TEB lowl + planck TT')
    plotVariable95(axarr[axNum],data_B,var, 'planck TEB')
    plotVariable95(axarr[axNum],data_C,var, 'planck lowl + S4')
    plotVariable95(axarr[axNum],data_E,var, 'PIXIE lowl + planck TEB highl')
    plotVariable95(axarr[axNum],data_F,var, 'PIXIE lowl + S4' )
    plt.tight_layout() plt.legend() plt.savefig('debug_get_95_percentiles.pdf')
    plt.show()
