

import matplotlib.pyplot as plt 
from astropy.io import ascii

outdir = 'output/'
namelist =  ( 'l', 'TT', 'EE', 'TE', 'BB', 'phiphi',  'TPhi', 'Ephi' )



cl_file = outdir + 'act100_cl.dat'
ad_file = outdir + 'act100_cls_ad.dat'
cdi_file = outdir + 'act100_cls_cdi.dat'
ad_cdi_file = outdir + 'act100_cls_ad_cdi.dat'

total_cl1 = ascii.read(cl_file, delimiter="\s", comment='#', names=namelist )
ad1 = ascii.read(ad_file, delimiter="\s", comment='#', names=namelist )
cdi1 = ascii.read(cdi_file, delimiter="\s", comment='#', names=namelist )
ad_cdi1 = ascii.read(ad_cdi_file, delimiter="\s", comment='#', names=namelist )

cl_file = outdir + 'act200_cl.dat'
ad_file = outdir + 'act200_cls_ad.dat'
cdi_file = outdir + 'act200_cls_cdi.dat'
ad_cdi_file = outdir + 'act200_cls_ad_cdi.dat'

total_cl2 = ascii.read(cl_file, delimiter="\s", comment='#', names=namelist )
ad2 = ascii.read(ad_file, delimiter="\s", comment='#', names=namelist )
cdi2 = ascii.read(cdi_file, delimiter="\s", comment='#', names=namelist )
ad_cdi2 = ascii.read(ad_cdi_file, delimiter="\s", comment='#', names=namelist )


l = total_cl1['l']


# plt.plot( l, cl['TT'], label="TT" )
# plt.plot( cl['l'], , 'r-' )


# plt.show()


fig, axarr = plt.subplots(3, sharex=True, figsize = (6,12))

plt.xlabel('$l$')


for ax, cor_type in zip( axarr, ('TT', 'TE', 'EE') ):
	# ax.plot(l, total_cl1[cor_type], "k--", label='C_l')

	ax.plot(l, ad1[cor_type]/100, 'r-', label='AD/100')
	ax.plot(l, cdi1[cor_type], 'g-', label='CDI')
	ax.plot(l, ad_cdi1[cor_type], 'b-', label='AD\_CDI')


	# ax.plot(l, ad2[cor_type]/100, 'r--')
	# ax.plot(l, cdi2[cor_type], 'g--')
	# ax.plot(l, ad_cdi2[cor_type], 'b--')


	ax.set_ylabel('$(l(l+1)/2\pi) C_l $')
	ax.set_title(cor_type)


	if cor_type=='TT':
		ax.legend()

plt.tight_layout()
plt.savefig('scaled_ad.pdf')


plt.clf()


fig, axarr = plt.subplots(3, sharex=True, figsize = (6,12))

plt.xlabel('$l$')


for ax, cor_type in zip( axarr, ('TT', 'TE', 'EE') ):
	# ax.plot(l, total_cl1[cor_type], "k--", label='C_l')

	ax.plot(l, ad1[cor_type]/100 - ad2[cor_type]/100, 'r-', label='$\Delta$AD/100')
	ax.plot(l, cdi1[cor_type] - cdi2[cor_type], 'g-', label='$\Delta$CDI')
	ax.plot(l, ad_cdi1[cor_type] - ad_cdi2[cor_type], 'b-', label='$\Delta$AD\_CDI')


	ax.set_ylabel('$(l(l+1)/2\pi) C_l $')
	ax.set_title(cor_type)


	if cor_type=='TT':
		ax.legend()

plt.tight_layout()
plt.savefig('delta_plotted.pdf')

# plt.show()