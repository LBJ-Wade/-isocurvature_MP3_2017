# must use python 2

from classy import Class
import matplotlib.pyplot as plt 
import numpy as np 
import math

max_l = 5000

max_scalars = '5000'

ell = np.array( range(1, max_l+1) ) 

def getDl( pii1=0.5e-10, pii2=1e-9, pri1=1e-13 ):
    # Define your cosmology (what is not specified will be set to CLASS default parameters)
    params = {
        'output': 'tCl lCl pCl', 
        'modes': 's', # scalar perturbations
        'lensing': 'yes',
        'ic': 'ad&cdi',
        'l_max_scalars':max_scalars,
        'P_k_ini type': 'two_scales',
        'k1': 0.002,
        'k2': 0.1,
        'P_{RR}^1': 2.34e-9,
        'P_{RR}^2': 2.115e-9,
        'P_{II}^1' : pii1,
        'P_{II}^2' : pii2,
        'P_{RI}^1' : pri1}

    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()

    # print(dir(cosmo)) # use this command to see what is in the cosmo

    # It is a dictionary that contains the fields: tt, te, ee, bb, pp, tp
    cls = cosmo.raw_cl(max_l)  # Access the  cl until l=1000

    yy = np.array( cls['ee'][1:] )
    zz = np.array( cls['tt'][1:] ) 
    yz = np.array( cls['te'][1:] ) 

    ee = ((ell)*(ell+1) * yy / (2 * math.pi))
    tt =  ((ell)*(ell+1) * zz / (2 * math.pi))
    te = ((ell)*(ell+1) * yz / (2 * math.pi))

    
    cosmo.struct_cleanup()
    return tt, te, ee

# Print on screen to see the output
# print len(cls['tt'])


pii1 = 0.5e-10
pii2 = 1e-9
pri1 = 1e-13
dpii1 = pii1 / 10000.0
dpii2 = pii2 / 10000.0
dpri1 = pri1 / 10000.0

pii1_tt1, pii1_te1, pii1_ee1 = getDl( pii1 = pii1 - dpii1 )
pii1_tt2, pii1_te2, pii1_ee2 = getDl( pii1 = pii1 + dpii1 )

pii2_tt1, pii2_te1, pii2_ee1 = getDl( pii2 = pii2 - dpii2 )
pii2_tt2, pii2_te2, pii2_ee2 = getDl( pii2 = pii2 + dpii2 )

pri1_tt1, pri1_te1, pri1_ee1 = getDl( pri1 = pri1 - dpri1 )
pri1_tt2, pri1_te2, pri1_ee2 = getDl( pri1 = pri1 + dpri1 )


# plot something with matplotlib...


plt.plot( (pii1_tt2 - pii1_tt1)/(2 * dpii1), label='$P_{II}^1$', markersize=0 )
plt.plot( (pii2_tt2 - pii2_tt1)/(2 * dpii2), label='$P_{II}^2$', markersize=0 )
# plt.plot( (pri1_tt2 - pri1_tt1)/(2 * dpri1), label='$P_{RI}^1$', markersize=0 )
plt.title('TT Derivatives')
plt.ylabel(r'$d \mathcal{D}_l / d P_{II}^j$')
plt.xlabel(r'$l$')
plt.legend()
plt.savefig('tt.pdf')
plt.clf()


plt.plot( (pii1_te2 - pii1_te1)/(2 * dpii1), label='$P_{II}^1$', markersize=0 )
plt.plot( (pii2_te2 - pii2_te1)/(2 * dpii2), label='$P_{II}^2$', markersize=0 )
# plt.plot( (pri1_te2 - pri1_te1)/(2 * dpri1), label='$P_{RI}^1$', markersize=0 )
plt.title('TE Derivatives')
plt.ylabel(r'$d \mathcal{D}_l / d P_{II}^j$')
plt.xlabel(r'$l$')
plt.legend()
plt.savefig('te.pdf')
plt.clf()


plt.plot( (pii1_ee2 - pii1_ee1)/(2 * dpii1), label='$P_{II}^1$', markersize=0 )
plt.plot( (pii2_ee2 - pii2_ee1)/(2 * dpii2), label='$P_{II}^2$', markersize=0 )
# plt.plot( (pri1_ee2 - pri1_ee1)/(2 * dpri1), label='$P_{RI}^1$', markersize=0 )
plt.title('EE Derivatives')
plt.ylabel(r'$d \mathcal{D}_l / d P_{II}^j$')
plt.xlabel(r'$l$')
plt.legend()
plt.savefig('ee.pdf')
plt.clf()

plt.plot( (pii1_ee2 - pii1_ee1)/(2 * dpii1), label='$P_{II}^1$', markersize=0 )
plt.plot( (pii2_ee2 - pii2_ee1)/(2 * dpii2), label='$P_{II}^2$', markersize=0 )
# plt.plot( (pri1_ee2 - pri1_ee1)/(2 * dpri1), label='$P_{RI}^1$', markersize=0 )
plt.title('EE Derivatives')
plt.ylabel(r'$d \mathcal{D}_l / d P_{II}^j$')
plt.xlabel(r'$l$')
plt.yscale('log')
plt.legend()
plt.savefig('logee.pdf')
plt.clf()


plt.plot( (pii1_tt2 - pii1_tt1)/(2 * dpii1), label='$P_{II}^1$', markersize=0 )
plt.plot( (pii2_tt2 - pii2_tt1)/(2 * dpii2), label='$P_{II}^2$', markersize=0 )
# plt.plot( (pri1_tt2 - pri1_tt1)/(2 * dpri1), label='$P_{RI}^1$', markersize=0 )
plt.title('TT Derivatives')
plt.ylabel(r'$d \mathcal{D}_l / d P_{II}^j$')
plt.xlabel(r'$l$')
plt.yscale('log')
plt.legend()
plt.savefig('logtt.pdf')
plt.clf()

