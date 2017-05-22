

import matplotlib
from getdist import plots, MCSamples

import getdist

import numpy as np
import matplotlib.pyplot as plt
import astropy
from loadMontePython import load as loadMCMC
import glob

basedir = '../chains/nonzero_model/'
#"/home/zequnl/Projects/isocurvature_2017/analysis/plot_triangle/nonzero/"

## F
folder = basedir + 'fF/'
files = glob.glob(folder + "*__*.txt")
params = glob.glob(folder + "*_.paramnames")
datalist = []
for f in files:
	datalist.append(  loadMCMC(f, params[0]) )
data = astropy.table.vstack( datalist )
data_sim = data[:]
weights_act = data['acceptance'][:]
for col in ['likelihood', 'acceptance','omega_b','omega_cdm','100theta_s','tau_reio']:
	data.remove_column(col)
nparr_act = np.array(data.as_array().tolist()[:])
pixie_s4 = MCSamples(samples=nparr_act,names = data.colnames, labels = data.colnames, name_tag='PIXIE low_l + S4')

## G
folder = basedir + 'fG/'
files = glob.glob(folder + "*__*.txt")
params = glob.glob(folder + "*_.paramnames")
datalist = []
for f in files:
	datalist.append(  loadMCMC(f, params[0]) )
data = astropy.table.vstack( datalist )
data_sim = data[:]
weights_act = data['acceptance'][:]
for col in ['likelihood', 'acceptance','omega_b','omega_cdm','100theta_s','tau_reio']:
	data.remove_column(col)
nparr_act = np.array(data.as_array().tolist()[:])
pixie_s4_nonoise = MCSamples(samples=nparr_act,names = data.colnames, labels = data.colnames, name_tag='PIXIE low_l + S4')


#Triangle plot
g = plots.getSubplotPlotter()
g.triangle_plot([ pixie_s4, pixie_s4_nonoise], filled=True)

# now we add some boundaries

# P_II^1
for ax in g.subplots[:,2]:
	if ax != None:
		ax.set_xlim(0,ax.get_xlim()[1])

for ax in g.subplots[2,:]:
	if ax != None:
		ax.set_ylim(0,ax.get_ylim()[1])


# P_II^2
for ax in g.subplots[:,3]:
	if ax != None:
		ax.set_xlim(0,ax.get_xlim()[1])


for ax in g.subplots[3,:]:
	if ax != None:
		ax.set_ylim(0,ax.get_ylim()[1])

plt.savefig('../../figures/compare_noise.pdf')
plt.show()
