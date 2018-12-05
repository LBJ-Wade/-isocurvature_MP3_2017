import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from scipy.stats import gaussian_kde
from pprint import pprint

import sys
import os
from astropy.io import ascii
from astropy.table import vstack

from plot_function import loadChainFolder, denplot, plotRow

fig = plt.figure(figsize=(8.5,11))
chainLocation = '/Users/zequnl/Desktop/dunkley2016/isocurvature_2017/analysis/chains/'


planck2015Data = loadChainFolder( 'chains/CDI_2/' )
def plotPlanck(ax1,ax2,ax3,ax4):
    plotRow(planck2015Data, ax1, ax2, ax3, ax4, 'k--', "Planck 2015")

figNum = 0

chainData = loadChainFolder( chainLocation + 'planck_nonzero_iso' )
ax1, ax2, ax3, ax4 = [fig.add_subplot(4,4,figNum+1),fig.add_subplot(4,4,figNum+2),\
    fig.add_subplot(4,4,figNum+3),fig.add_subplot(4,4,figNum+4)]
plotRow(chainData, ax1, ax2, ax3, ax4, 'b-', "Forecasted Planck Pol")
plotPlanck(ax1,ax2,ax3,ax4)
figNum += 4
top_axis = ax1

chainData = loadChainFolder( chainLocation + 'planck_s4_nonzero_iso' )
ax1, ax2, ax3, ax4 = [fig.add_subplot(4,4,figNum+1),fig.add_subplot(4,4,figNum+2),\
    fig.add_subplot(4,4,figNum+3),fig.add_subplot(4,4,figNum+4)]
plotRow(chainData, ax1, ax2, ax3, ax4, 'r-', "Forecasted Planck low\_l + S4")
plotPlanck(ax1,ax2,ax3,ax4)
figNum += 4


chainData = loadChainFolder( chainLocation + 'planck_s4_nonzero_iso' )
ax1, ax2, ax3, ax4 = [fig.add_subplot(4,4,figNum+1),fig.add_subplot(4,4,figNum+2),\
    fig.add_subplot(4,4,figNum+3),fig.add_subplot(4,4,figNum+4)]
plotRow(chainData, ax1, ax2, ax3, ax4, 'g-', "Forecasted PIXIE + Planck Pol high\_l")
plotPlanck(ax1,ax2,ax3,ax4)
figNum += 4


chainData = loadChainFolder( chainLocation + 'pixie_s4_nonzero_iso' )
ax1, ax2, ax3, ax4 = [fig.add_subplot(4,4,figNum+1),fig.add_subplot(4,4,figNum+2),\
    fig.add_subplot(4,4,figNum+3),fig.add_subplot(4,4,figNum+4)]
plotRow(chainData, ax1, ax2, ax3, ax4, 'c-', "Forecasted PIXIE + S4")
plotPlanck(ax1,ax2,ax3,ax4)
figNum += 4


plt.tight_layout()
# plt.show()
plt.savefig("../../figures/all_derived_forecast_nonzero.pdf")




# now generate an overplot ------------------------------

plt.clf()
fig = plt.figure(figsize=(8.5,3))

figNum = 0
chainData = loadChainFolder( chainLocation + 'planck_nonzero_iso' )
ax1, ax2, ax3, ax4 = [fig.add_subplot(1,4,figNum+1),fig.add_subplot(1,4,figNum+2),\
    fig.add_subplot(1,4,figNum+3),fig.add_subplot(1,4,figNum+4)]
plotRow(chainData, ax1, ax2, ax3, ax4, 'b-', "Forecasted Planck Pol")

chainData = loadChainFolder( chainLocation + 'planck_s4_nonzero_iso' )
ax1, ax2, ax3, ax4 = [fig.add_subplot(1,4,figNum+1),fig.add_subplot(1,4,figNum+2),\
    fig.add_subplot(1,4,figNum+3),fig.add_subplot(1,4,figNum+4)]
plotRow(chainData, ax1, ax2, ax3, ax4, 'r-', "Forecasted Planck low\_l + S4")

chainData = loadChainFolder( chainLocation + 'planck_s4_nonzero_iso' )
ax1, ax2, ax3, ax4 = [fig.add_subplot(1,4,figNum+1),fig.add_subplot(1,4,figNum+2),\
    fig.add_subplot(1,4,figNum+3),fig.add_subplot(1,4,figNum+4)]
plotRow(chainData, ax1, ax2, ax3, ax4, 'g-', "Forecasted PIXIE + Planck Pol high\_l")

chainData = loadChainFolder( chainLocation + 'pixie_s4_nonzero_iso' )
ax1, ax2, ax3, ax4 = [fig.add_subplot(1,4,figNum+1),fig.add_subplot(1,4,figNum+2),\
    fig.add_subplot(1,4,figNum+3),fig.add_subplot(1,4,figNum+4)]
plotRow(chainData, ax1, ax2, ax3, ax4, 'c-', "Forecasted PIXIE + S4")

plt.tight_layout()
# plt.show()
plt.savefig("../../figures/overplotted_derived_nonzero.pdf")
