# this file generates the overplots!

import matplotlib
from getdist import plots, MCSamples
import getdist

import numpy as np
import matplotlib.pyplot as plt
import astropy
from loadMontePython import load as loadMCMC
import glob

basedir = '../chains/'
#"/home/zequnl/Projects/isocurvature_2017/analysis/plot_triangle/nonzero/"

def plot_chain(folder, name_tag):
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
	if 'A_planck' in data.colnames:
		data.remove_column('A_planck')
	nparr_act = np.array(data.as_array().tolist()[:])
	return MCSamples(samples=nparr_act,names = data.colnames, labels = data.colnames, name_tag=name_tag)

def fix_boundaries(g):

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

planck = plot_chain(basedir + 'planckdata/', 'Planck 2015')
planckpol = plot_chain(basedir + 'planck_zero_iso/', 'Forecasted Planck Pol')
planck_s4 = plot_chain(basedir + 'planck_s4_zero_iso/', 'Planck + S4')
pixie_planck = plot_chain(basedir + 'pixie_planck_zero_iso/', 'PIXIE + Planck')
pixie_s4 = plot_chain(basedir + 'pixie_s4_zero_iso/', 'PIXIE + S4')

#Triangle plot -- zero iso model
g = plots.getSubplotPlotter()
g.triangle_plot([ planckpol, planck_s4, pixie_planck, pixie_s4], filled=True)
fix_boundaries(g)

plt.savefig('../../figures/overplot_zero_iso.pdf')
plt.clf()

# Triangle plot -- zero iso model with real planck data overplotted
g = plots.getSubplotPlotter()
g.triangle_plot([ planck, planckpol, planck_s4, pixie_planck, pixie_s4], filled=True)
fix_boundaries(g)

plt.savefig('../../figures/overplot_zero_iso_with_planck2015.pdf')
plt.clf()



#Triangle plot with nonzero iso model
planckpol = plot_chain(basedir + 'planck_nonzero_iso/', 'Forecasted Planck Pol')
planck_s4 = plot_chain(basedir + 'planck_s4_nonzero_iso/', 'Planck + S4')
pixie_planck = plot_chain(basedir + 'pixie_planck_nonzero_iso/', 'PIXIE + Planck')
pixie_s4 = plot_chain(basedir + 'pixie_s4_nonzero_iso/', 'PIXIE + S4')

g = plots.getSubplotPlotter()
g.triangle_plot([ planckpol, planck_s4, pixie_planck, pixie_s4], filled=True)
fix_boundaries(g)


plt.savefig('../../figures/overplot_nonzero_iso.pdf')
