import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec as gd
import os
import pickle as pc
import utils1 as utl
from astropy.io import fits
from astroquery.mast import Observations as obs
from scipy import interpolate as inp
import matplotlib.cm as cm
import matplotlib.colors as cls
from scipy.optimize import curve_fit as cft
from pylab import *
import seaborn as sns


#path1 = input('Enter the path of this folder: ')
path1 = '/home/jayshil/Documents/Dissertation'

#---------------------------------------------------------------------------------------------
#--------------------------Taking Data from the data file-------------------------------------
#---------------------------------------------------------------------------------------------

name = np.loadtxt('data2.dat', dtype = str, usecols = 0, unpack = True)
teff, lg, mh, vturb, p, pperr, pnerr, tc, aste, asteperr, astenerr, ecc, ome, rprst, rprstperr, rprstnerr, tce = np.loadtxt('data2.dat', usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), unpack = True)

#--------------------------------------------------------------------------------------------------
#-------------------------Now its time to plot some amazing results--------------------------------
#--------------------------------------------------------------------------------------------------

sns.set_context("talk")
sns.set_style("ticks")

# Fonts:
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size':12})
plt.rc('legend', **{'fontsize':12})

# Ticks to the outside:
rcParams['axes.linewidth'] = 1.2 
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'


print('-------------------------------------------------------------------------')
print('---------------------Plotting Stellar Parameters-------------------------')
print('-------------------------------------------------------------------------')

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.0175

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a rectangular Figure
stellar_prop = plt.figure(figsize=(9, 9))

ax_scatter = plt.axes(rect_scatter)
ax_scatter.tick_params(direction='in', top=True, right=True, labelrotation=67.5)
ax_histx = plt.axes(rect_histx)
ax_histx.tick_params(direction='in', labelbottom=False)
ax_histy = plt.axes(rect_histy)
ax_histy.tick_params(direction='in', labelleft=False)

# the scatter plot:
ax_scatter.grid()
ax_scatter.scatter(teff, lg, marker='o', c='royalblue')
ax_scatter.set_xlabel('Effective Temperature of stellar host (in K)')
ax_scatter.set_ylabel('Surface Gravity of stellar host (log (g))')

ax_scatter.set_xlim((3000, 7000))
ax_scatter.set_ylim((3.6, 5.2))

#bins = np.arange(-lim, lim + binwidth, binwidth)
ax_histx.hist(teff, bins=utl.freedman_diaconis(data=teff, returnas="bins"), color='orangered')
ax_histx.set_ylabel('Number of stars')
ax_histy.hist(lg, bins=utl.freedman_diaconis(data=lg, returnas="bins"), orientation='horizontal', color='orangered')
ax_histy.set_xlabel('Number of stars')
ax_histx.set_xlim(ax_scatter.get_xlim())
ax_histy.set_ylim(ax_scatter.get_ylim())

plt.savefig(path1 + '/Results/stellar_prop/1.pdf')
plt.close(stellar_prop)


print('-------------------------------------------------------------------------')
print('----------------Plots to check how good was the fit----------------------')
print('-------------------------------------------------------------------------')

amax_t = np.max(aste)
pmax_t = np.max(p)
rmax_t = np.max(rprst)
amin_t = np.min(aste)
pmin_t = np.min(p)
rmin_t = np.min(rprst)

#------------------------------------------------
print('#-----------------For a/R*-----------------------')
#------------------------------------------------

a_j, a_jp, a_jn = np.loadtxt(path1 + '/Results/comp_a_r_p/aste.dat', usecols = (1,2,3), unpack = True)

amax_j = np.max(a_j)
amin_j = np.min(a_j)

xla_j = np.minimum(amin_t, amin_j)
xua_j = np.maximum(amax_t, amax_j)

x1a = y1a = np.linspace(xla_j, xua_j, 100)
y11a = np.zeros(len(x1a))

diff_a = np.array([])
diff_ae = np.array([])

for i in range(len(a_j)):
	at1 = np.random.normal(aste[i],asteperr[i],10000)
	ac1 = np.random.normal(a_j[i], a_jp[i], 10000)
	diff1 = at1 - ac1
	am1 = np.median(diff1)
	ae1 = np.std(diff1)
	diff_a = np.hstack((diff_a,am1))
	diff_ae = np.hstack((diff_ae,ae1))

fig_a = plt.figure(figsize = (8,10))
gs_a = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_a = plt.subplot(gs_a[0])

ax_a.errorbar(aste, a_j, xerr = [astenerr, asteperr], yerr = [a_jn, a_jp], fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xla_j, xua_j])
plt.ylim([xla_j, xua_j])
ax_a.plot(x1a, y1a, 'k--')
ax_a.grid()
plt.ylabel(r'$a/R_*$ (Observed)')
#plt.title('Comparison between literature values and calculated values of a/R*')

ax1_a = plt.subplot(gs_a[1], sharex = ax_a)

ax1_a.errorbar(aste, diff_a, xerr = [astenerr, asteperr], yerr = diff_ae, fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xla_j, xua_j])
plt.ylim([-2,2])
ax1_a.plot(x1a, y11a, 'k--')
ax1_a.grid()
plt.ylabel('Residuals')
plt.xlabel(r'Values of $a/R_*$ taken the from the literature')

plt.subplots_adjust(hspace = 0.2)
plt.savefig(path1 + '/Results/comp_a_r_p/a.pdf')

diff_a1 = np.abs(diff_a)
erra1 = 3*(asteperr + a_jp)

fa = open(path1 + '/Results/comp_a_r_p/off_a.dat', 'w')
fa.write('#These are the systems which have residulas larger than 3-sigma\n')

for i in range(len(diff_a)):
	if diff_a1[i] > erra1[i]:
		fa.write(name[i] + '\t' + str(a_j[i]) + '+-' + str(a_jp[i]) + '\t' + str(aste[i]) + '+-' + str(asteperr[i]) + '\n')

fa.close()
plt.close(fig_a)

#------------------------------------------------
print('#-----------------For Period---------------------')
#------------------------------------------------

p_j, p_jp, p_jn = np.loadtxt(path1 + '/Results/comp_a_r_p/period.dat', usecols = (1,2,3), unpack = True)

pmax_j = np.max(p_j)
pmin_j = np.min(p_j)

xlp_j = np.minimum(pmin_t, pmin_j)
xup_j = np.maximum(pmax_t, pmax_j)

x1p = y1p = np.linspace(xlp_j, xup_j, 100)
y11p = np.zeros(len(x1p))

diff_p = np.array([])
diff_pe = np.array([])

for i in range(len(p_j)):
	pt1 = np.random.normal(p[i],pperr[i],10000)
	pc1 = np.random.normal(p_j[i], p_jp[i], 10000)
	diff1 = (pt1 - pc1)*86400
	pm1 = np.median(diff1)
	pe1 = np.std(diff1)
	diff_p = np.hstack((diff_p,pm1))
	diff_pe = np.hstack((diff_pe,pe1))

fig_p = plt.figure(figsize = (8,10))
gs_p = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_p = plt.subplot(gs_p[0])

ax_p.errorbar(p, p_j, xerr = [pnerr, pperr], yerr = [p_jn, p_jp], fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xlp_j, xup_j])
plt.ylim([xlp_j, xup_j])
ax_p.plot(x1p, y1p, 'k--')
ax_p.grid()

plt.ylabel('Period (Observed - in days)')
#plt.title('Comparison between literature values and calculated values of period')

ax1_p = plt.subplot(gs_p[1], sharex = ax_p)

ax1_p.errorbar(p, diff_p, xerr = [pnerr, pperr], yerr = diff_pe, fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xlp_j, xup_j])
plt.ylim([-3,3])
ax1_p.plot(x1p, y11p, 'k--')
ax1_p.grid()
plt.ylabel('Residuals (in sec)')
plt.xlabel('Values of period taken from the literature')

plt.subplots_adjust(hspace = 0.2)
plt.savefig(path1 + '/Results/comp_a_r_p/period.pdf')

diff_p1 = np.abs(diff_p)
errp1 = 3*(pperr + p_jp)*86400

fp = open(path1 + '/Results/comp_a_r_p/off_p.dat', 'w')
fp.write('These are the systems which have residulas larger than 3-sigma\n')

for i in range(len(diff_p)):
	if diff_p1[i] > errp1[i]:
		fp.write(name[i] + '\t' + str(p_j[i]) + '+-' + str(p_jp[i]) + '\t' + str(p[i]) + '+-' + str(pperr[i]) + '\n')

fp.close()
plt.close(fig_p)

#------------------------------------------------
print('#-----------------For Rp/R*----------------------')
#------------------------------------------------

r_j, r_jp, r_jn = np.loadtxt(path1 + '/Results/comp_a_r_p/rprste.dat', usecols = (1,2,3), unpack = True)

rmax_j = np.max(r_j)
rmin_j = np.min(r_j)

xlr_j = np.minimum(rmin_t, rmin_j)
xur_j = np.maximum(rmax_t, rmax_j)

x1r = y1r = np.linspace(xlr_j, xur_j, 100)
y11r = np.zeros(len(x1r))


diff_r = np.array([])
diff_re = np.array([])

for i in range(len(r_j)):
	rt1 = np.random.normal(rprst[i],rprstperr[i],10000)
	rc1 = np.random.normal(r_j[i], r_jp[i], 10000)
	diff1 = rt1 - rc1
	rm1 = np.median(diff1)
	re1 = np.std(diff1)
	diff_r = np.hstack((diff_r,rm1))
	diff_re = np.hstack((diff_re,re1))

fig_r = plt.figure(figsize = (8,10))
gs_r = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_r = plt.subplot(gs_r[0])

ax_r.errorbar(rprst, r_j, xerr = [rprstnerr, rprstperr], yerr = [r_jn, r_jp], fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xlr_j, xur_j])
plt.ylim([xlr_j, xur_j])
ax_r.plot(x1r, y1r, 'k--')
ax_r.grid()

plt.ylabel(r'$R_p/R_*$ (Observed)')
#plt.title('Comparison between literature values and calculated values of Rp/R*')

ax1_r = plt.subplot(gs_r[1], sharex = ax_r)

ax1_r.errorbar(rprst, diff_r, xerr = [rprstnerr, rprstperr], yerr = diff_re, fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xlr_j, xur_j])
plt.ylim([-0.03,0.03])
ax1_r.plot(x1r, y11r, 'k--')
ax1_r.grid()
plt.ylabel('Residuals')
plt.xlabel(r'Values of $R_p/R_*$ taken from the literature')

plt.subplots_adjust(hspace = 0.2)
plt.savefig(path1 + '/Results/comp_a_r_p/r.pdf')

diff_r1 = np.abs(diff_r)
errr1 = 3*(rprstperr + r_jp)

fr = open(path1 + '/Results/comp_a_r_p/off_r.dat', 'w')
fr.write('These are the systems which have residulas larger than 3-sigma\n')

for i in range(len(diff_r)):
	if diff_r1[i] > errr1[i]:
		fr.write(name[i] + '\t' + str(r_j[i]) + '+-' + str(r_jp[i]) + '\t' + str(rprst[i]) + '+-' + str(rprstperr[i]) + '\n')

fr.close()
plt.close(fig_r)

#------------------------------------------------
print('#----------------For Tc--------------------------')
#------------------------------------------------
tc_t1, tc_te1 = np.loadtxt(path1 + '/Results/comp_a_r_p/to_the.dat', usecols = (1,2), unpack = True)

tc_t = np.array([])
tc_te = np.array([])

for i in range(len(tc_t1)):
	t1 = np.random.normal(tc_t1[i], tc_te1[i], 10000)
	t2 = t1 - 2457000
	t3 = np.median(t2)
	t4 = np.std(t2)
	tc_t = np.hstack((tc_t,t3))
	tc_te = np.hstack((tc_te,t4))

tc_j1, tc_jp1, tc_jn1 = np.loadtxt(path1 + '/Results/comp_a_r_p/tc.dat', usecols = (1,2,3), unpack = True)

tc_j = np.array([])
tc_jp = np.array([])
tc_jn = tc_jn1

for i in range(len(tc_j1)):
	t5 = np.random.normal(tc_j1[i], tc_jp1[i], 10000)
	t6 = t5 - 2457000
	t7 = np.median(t6)
	t8 = np.std(t6)
	tc_j = np.hstack((tc_j,t7))
	tc_jp = np.hstack((tc_jp,t8))

tmax_t = np.max(tc_t)
tmin_t = np.min(tc_t)

tcmax_j = np.max(tc_j)
tcmin_j = np.min(tc_j)

xlt_j = np.minimum(tmin_t, tcmin_j)
xut_j = np.maximum(tmax_t, tcmax_j)

x1t = y1t = np.linspace(xlt_j, xut_j, 100)
y11t = np.zeros(len(x1t))

diff_t = np.array([])
diff_te = np.array([])

for i in range(len(tc_j)):
	tt1 = np.random.normal(tc_t[i],tc_te[i],10000)
	tc1 = np.random.normal(tc_j[i], tc_jp[i], 10000)
	diff1 = (tt1 - tc1)*1440
	tm1 = np.median(diff1)
	te1 = np.std(diff1)
	diff_t = np.hstack((diff_t,tm1))
	diff_te = np.hstack((diff_te,te1))

fig_t = plt.figure(figsize = (8,10))
gs_t = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_t = plt.subplot(gs_t[0])

ax_t.errorbar(tc_t, tc_j, xerr = tc_te, yerr = [tc_jn, tc_jp], fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xlt_j, xut_j])
plt.ylim([xlt_j, xut_j])
ax_t.plot(x1t, y1t, 'k--')
ax_t.grid()

plt.ylabel(r'$T_c$ (Observed - in TJD)')
#plt.title('Comparison between literature values and calculated values of Tc')

ax1_t = plt.subplot(gs_t[1], sharex = ax_t)

ax1_t.errorbar(tc_t, diff_t, xerr = tc_te, yerr = diff_te, fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xlt_j, xut_j])
plt.ylim([-15,15])
ax1_t.plot(x1t, y11t, 'k--')
ax1_t.grid()
plt.ylabel('Residuals (min)')
plt.xlabel(r'Values of $T_c$ taken from the literature (in TJD)')

plt.subplots_adjust(hspace = 0.2)
plt.savefig(path1 + '/Results/comp_a_r_p/tc.pdf')
plt.close(fig_t)

diff_t1 = np.abs(diff_t)
errt1 = 3*(tc_te + tc_jp)*1440

ft = open(path1 + '/Results/comp_a_r_p/off_t.dat', 'w')
ft.write('These are the systems which have residulas larger than 3-sigma\n')

for i in range(len(diff_t)):
	if diff_t1[i] > errt1[i]:
		ft.write(name[i] + '\t' + str(tc_j[i]) + '+-' + str(tc_jp[i]) + '\t' + str(tc_t[i]) + '+-' + str(tc_te[i]) + '\n')

ft.close()


print('-------------------------------------------------------------------------')
print('-------------------------Plots to compare LDCs---------------------------')
print('-------------------------------------------------------------------------')

u1_j, u1_jp, u1_jn, u2_j, u2_jp, u2_jn = np.loadtxt(path1 + '/Results/cal_us_and_evidance/cal_u1_u2.dat', usecols = (1,2,3,4,5,6), unpack = True)

u1max_j = np.max(u1_j)
u1min_j = np.min(u1_j)
u2max_j = np.max(u2_j)
u2min_j = np.min(u2_j)

#---------------------------------------------------------
print('--How good are the Claret(2017) (With Phoenix-s) LDCs?---')
#---------------------------------------------------------

u1_c_p, u2_c_p = np.loadtxt(path1 + '/Phoenix/claret_limiting_LDC_pho.dat', usecols = (1,2), unpack = True)

u1max_c_p = np.max(u1_c_p)
u1min_c_p = np.min(u1_c_p)
u2max_c_p = np.max(u2_c_p)
u2min_c_p = np.min(u2_c_p)

u1_c_a, u2_c_a = np.loadtxt(path1 + '/Atlas/claret_limiting_LDC_ata.dat', usecols = (1,2), unpack = True)

u1max_c_a = np.max(u1_c_a)
u1min_c_a = np.min(u1_c_a)
u2max_c_a = np.max(u2_c_a)
u2min_c_a = np.min(u2_c_a)

#-----------------------
print('#--------u1-------------')
#-----------------------

xlu1_c_p = np.minimum(u1min_j, u1min_c_p)
xuu1_c_p = np.maximum(u1max_j, u1max_c_p)

xlu1_c_a = np.minimum(u1min_j, u1min_c_a)
xuu1_c_a = np.maximum(u1max_j, u1max_c_a)

xlo = np.minimum(xlu1_c_p, xlu1_c_a)
xup = np.maximum(xuu1_c_p, xuu1_c_a)

x1u1_c_p = y1u1_c_p = np.linspace(xlo, xup, 100)
y11u1_c_p = np.zeros(len(x1u1_c_p))

diff_u1_c_p = np.array([])#--------------------------------------------------------------------------------------------
diff_u1_c_pe = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u1_j)):
	u11_c_p = np.random.normal(u1_c_p[i], 0, 10000)
	u11_j = np.random.normal(u1_j[i], u1_jp[i], 10000)
	diff1 = u11_c_p - u11_j
	u11_m = np.median(diff1)
	u11_e = np.std(diff1)
	diff_u1_c_p = np.hstack((diff_u1_c_p, u11_m))
	diff_u1_c_pe = np.hstack((diff_u1_c_pe, u11_e))

diff_u1_c_a = np.array([])#--------------------------------------------------------------------------------------------
diff_u1_c_ae = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u1_j)):
	u11_c_a = np.random.normal(u1_c_a[i], 0, 10000)
	u11_j = np.random.normal(u1_j[i], u1_jp[i], 10000)
	diff1 = u11_c_a - u11_j
	u11_m = np.median(diff1)
	u11_e = np.std(diff1)
	diff_u1_c_a = np.hstack((diff_u1_c_a, u11_m))
	diff_u1_c_ae = np.hstack((diff_u1_c_ae, u11_e))


fig_u1_c_p = plt.figure(figsize=(8,10))
gs_u1_c_p = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_u1_c_p = plt.subplot(gs_u1_c_p[0])

ax_u1_c_p.errorbar(u1_j, u1_c_p, xerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
ax_u1_c_p.errorbar(u1_j, u1_c_a, xerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

ax_u1_c_p.plot(x1u1_c_p, y1u1_c_p, 'k--')
ax_u1_c_p.grid()

ax_u1_c_p.set_xlim([xlo, xup])
ax_u1_c_p.set_ylim([xlo, xup])

plt.legend(loc='best')
plt.ylabel(r'$u_1$ (Theoretical)')
plt.xlabel(r'$u_1$ (Empirical)')
plt.title('Values from Claret(2017)')

ax1_u1_c_p = plt.subplot(gs_u1_c_p[1])#, sharex = ax_u1_c_p)

ax1_u1_c_p.hist(diff_u1_c_p, bins=utl.freedman_diaconis(data=diff_u1_c_p, returnas="bins"), alpha=0.7, color='orangered', zorder=5)
ax1_u1_c_p.hist(diff_u1_c_a, bins=utl.freedman_diaconis(data=diff_u1_c_a, returnas="bins"), alpha=0.7, color='cornflowerblue', zorder=5)

plt.ylabel('Count')
plt.xlabel('Residuals')

plt.subplots_adjust(hspace = 0.3)
plt.savefig(path1 + '/Results/cal_us_and_evidance/u1_cla.pdf')
plt.close(fig_u1_c_p)

#-----------------------
print('#--------u2-------------')
#-----------------------

xlu2_c_p = np.minimum(u2min_j, u2min_c_p)
xuu2_c_p = np.maximum(u2max_j, u2max_c_p)

xlu2_c_a = np.minimum(u2min_j, u2min_c_a)
xuu2_c_a = np.maximum(u2max_j, u2max_c_a)

xlo = np.minimum(xlu2_c_p, xlu2_c_a)
xup = np.maximum(xuu2_c_p, xuu2_c_a)

x1u2_c_p = y1u2_c_p = np.linspace(xlo, xup, 100)
y11u2_c_p = np.zeros(len(x1u2_c_p))

diff_u2_c_p = np.array([])#--------------------------------------------------------------------------------------------
diff_u2_c_pe = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u2_j)):
	u22_c_p = np.random.normal(u2_c_p[i], 0, 10000)
	u22_j = np.random.normal(u2_j[i], u2_jp[i], 10000)
	diff2 = u22_c_p - u22_j
	u22_m = np.median(diff2)
	u22_e = np.std(diff2)
	diff_u2_c_p = np.hstack((diff_u2_c_p, u22_m))
	diff_u2_c_pe = np.hstack((diff_u2_c_pe, u22_e))

diff_u2_c_a = np.array([])#--------------------------------------------------------------------------------------------
diff_u2_c_ae = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u2_j)):
	u22_c_a = np.random.normal(u2_c_a[i], 0, 10000)
	u22_j = np.random.normal(u2_j[i], u2_jp[i], 10000)
	diff2 = u22_c_a - u22_j
	u22_m = np.median(diff2)
	u22_e = np.std(diff2)
	diff_u2_c_a = np.hstack((diff_u2_c_a, u22_m))
	diff_u2_c_ae = np.hstack((diff_u2_c_ae, u22_e))


fig_u2_c_p = plt.figure(figsize=(8,10))
gs_u2_c_p = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_u2_c_p = plt.subplot(gs_u2_c_p[0])

ax_u2_c_p.errorbar(u2_j, u2_c_p, xerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
ax_u2_c_p.errorbar(u2_j, u2_c_a, xerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

ax_u2_c_p.plot(x1u2_c_p, y1u2_c_p, 'k--')
ax_u2_c_p.grid()

ax_u2_c_p.set_xlim([xlo, xup])
ax_u2_c_p.set_ylim([xlo, xup])

plt.legend(loc='best')
plt.ylabel(r'$u_2$ (Theoretical)')
plt.xlabel(r'$u_2$ (Empirical)')
plt.title('Values from Claret(2017)')

ax1_u2_c_p = plt.subplot(gs_u2_c_p[1])#, sharex = ax_u2_c_p)

ax1_u2_c_p.hist(diff_u2_c_p, bins=utl.freedman_diaconis(data=diff_u2_c_p, returnas="bins"), alpha=0.7, color='orangered', zorder=0.5)
ax1_u2_c_p.hist(diff_u2_c_a, bins=utl.freedman_diaconis(data=diff_u2_c_a, returnas="bins"), alpha=0.7, color='cornflowerblue', zorder=0.5)

plt.ylabel('Count')
plt.xlabel('Residuals')

plt.subplots_adjust(hspace = 0.3)
plt.savefig(path1 + '/Results/cal_us_and_evidance/u2_cla.pdf')
plt.close(fig_u2_c_p)

#---------------------------------------------------------
print('--How good are the Claret(2017) (With Phoenix-r) LDCs?---')
#---------------------------------------------------------

u1_c_pr, u2_c_pr = np.loadtxt(path1 + '/Phoenix/claret_limiting_LDC_pho.dat', usecols = (1,2), unpack = True)

u1max_c_pr = np.max(u1_c_pr)
u1min_c_pr = np.min(u1_c_pr)
u2max_c_pr = np.max(u2_c_pr)
u2min_c_pr = np.min(u2_c_pr)

u1_c_ar, u2_c_ar = np.loadtxt(path1 + '/Phoenix/claret_limiting_LDC_pho_r.dat', usecols = (1,2), unpack = True)

u1max_c_ar = np.max(u1_c_ar)
u1min_c_ar = np.min(u1_c_ar)
u2max_c_ar = np.max(u2_c_ar)
u2min_c_ar = np.min(u2_c_ar)

#-----------------------
print('#--------u1-------------')
#-----------------------

xlu1_c_pr = np.minimum(u1min_j, u1min_c_pr)
xuu1_c_pr = np.maximum(u1max_j, u1max_c_pr)

xlu1_c_ar = np.minimum(u1min_j, u1min_c_ar)
xuu1_c_ar = np.maximum(u1max_j, u1max_c_ar)

xlo = np.minimum(xlu1_c_pr, xlu1_c_ar)
xup = np.maximum(xuu1_c_pr, xuu1_c_ar)

x1u1_c_pr = y1u1_c_pr = np.linspace(xlo, xup, 100)
y11u1_c_pr = np.zeros(len(x1u1_c_pr))

diff_u1_c_pr = np.array([])#--------------------------------------------------------------------------------------------
diff_u1_c_per = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u1_j)):
	u11_c_pr = np.random.normal(u1_c_pr[i], 0, 10000)
	u11_jr = np.random.normal(u1_j[i], u1_jp[i], 10000)
	diff1r = u11_c_pr - u11_jr
	u11_m = np.median(diff1r)
	u11_e = np.std(diff1r)
	diff_u1_c_pr = np.hstack((diff_u1_c_pr, u11_m))
	diff_u1_c_per = np.hstack((diff_u1_c_per, u11_e))

diff_u1_c_ar = np.array([])#--------------------------------------------------------------------------------------------
diff_u1_c_aer = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u1_j)):
	u11_c_ar = np.random.normal(u1_c_ar[i], 0, 10000)
	u11_j = np.random.normal(u1_j[i], u1_jp[i], 10000)
	diff1 = u11_c_a - u11_j
	u11_m = np.median(diff1)
	u11_e = np.std(diff1)
	diff_u1_c_ar = np.hstack((diff_u1_c_ar, u11_m))
	diff_u1_c_aer = np.hstack((diff_u1_c_aer, u11_e))


fig_u1_c_pr = plt.figure(figsize=(8,10))
gs_u1_c_pr = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_u1_c_pr = plt.subplot(gs_u1_c_pr[0])

ax_u1_c_pr.errorbar(u1_j, u1_c_pr, xerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='deepskyblue', zorder=5, label = 'PHOENIX LDCs - s method')
ax_u1_c_pr.errorbar(u1_j, u1_c_ar, xerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='darkred',zorder=5, label = 'PHOENIX LDCs - r method')

ax_u1_c_pr.plot(x1u1_c_pr, y1u1_c_pr, 'k--')
ax_u1_c_pr.grid()

ax_u1_c_pr.set_xlim([xlo, xup])
ax_u1_c_pr.set_ylim([xlo, xup])

plt.legend(loc='best')
plt.ylabel(r'$u_1$ (Theoretical)')
plt.xlabel(r'$u_1$ (Empirical)')
plt.title('Values from Claret(2017)')

ax1_u1_c_pr = plt.subplot(gs_u1_c_pr[1])#, sharex = ax_u1_c_pr)

ax1_u1_c_pr.hist(diff_u1_c_pr, bins=utl.freedman_diaconis(data=diff_u1_c_pr, returnas="bins"), alpha=0.7, color='deepskyblue', zorder=5)
ax1_u1_c_pr.hist(diff_u1_c_ar, bins=utl.freedman_diaconis(data=diff_u1_c_ar, returnas="bins"), alpha=0.7, color='darkred',zorder=5)

plt.xlabel('Residuals')
plt.ylabel('Count')

plt.subplots_adjust(hspace = 0.3)
plt.savefig(path1 + '/Results/cal_us_and_evidance/u1_cla_r.pdf')
plt.close(fig_u1_c_pr)

#-----------------------
print('#--------u2-------------')
#-----------------------

xlu2_c_pr = np.minimum(u2min_j, u2min_c_pr)
xuu2_c_pr = np.maximum(u2max_j, u2max_c_pr)

xlu2_c_ar = np.minimum(u2min_j, u2min_c_ar)
xuu2_c_ar = np.maximum(u2max_j, u2max_c_ar)

xlo = np.minimum(xlu2_c_pr, xlu2_c_ar)
xup = np.maximum(xuu2_c_pr, xuu2_c_ar)

x1u2_c_pr = y1u2_c_pr = np.linspace(xlo, xup, 100)
y11u2_c_pr = np.zeros(len(x1u2_c_pr))

diff_u2_c_pr = np.array([])#--------------------------------------------------------------------------------------------
diff_u2_c_per = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u2_j)):
	u22_c_pr = np.random.normal(u2_c_pr[i], 0, 10000)
	u22_j = np.random.normal(u2_j[i], u2_jp[i], 10000)
	diff2 = u22_c_pr - u22_j
	u22_m = np.median(diff2)
	u22_e = np.std(diff2)
	diff_u2_c_pr = np.hstack((diff_u2_c_pr, u22_m))
	diff_u2_c_per = np.hstack((diff_u2_c_per, u22_e))

diff_u2_c_ar = np.array([])#--------------------------------------------------------------------------------------------
diff_u2_c_aer = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u2_j)):
	u22_c_ar = np.random.normal(u2_c_ar[i], 0, 10000)
	u22_j = np.random.normal(u2_j[i], u2_jp[i], 10000)
	diff2 = u22_c_ar - u22_j
	u22_m = np.median(diff2)
	u22_e = np.std(diff2)
	diff_u2_c_ar = np.hstack((diff_u2_c_ar, u22_m))
	diff_u2_c_aer = np.hstack((diff_u2_c_aer, u22_e))


fig_u2_c_pr = plt.figure(figsize=(8,10))
gs_u2_c_pr = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_u2_c_pr = plt.subplot(gs_u2_c_pr[0])

ax_u2_c_pr.errorbar(u2_j, u2_c_pr, xerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='deepskyblue', zorder=5, label = 'PHOENIX LDCs - s method')
ax_u2_c_pr.errorbar(u2_j, u2_c_ar, xerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='darkred',zorder=5, label = 'PHOENIX LDCs - r method')

ax_u2_c_pr.plot(x1u2_c_pr, y1u2_c_pr, 'k--')
ax_u2_c_pr.grid()

ax_u2_c_pr.set_xlim([xlo, xup])
ax_u2_c_pr.set_ylim([xlo, xup])

plt.legend(loc='best')
plt.ylabel(r'$u_2$ (Theoretical)')
plt.xlabel(r'$u_2$ (Empirical)')
plt.title('Values from Claret(2017)')

ax1_u2_c_pr = plt.subplot(gs_u2_c_pr[1])#, sharex = ax_u2_c_pr)

ax1_u2_c_pr.hist(diff_u2_c_pr, bins=utl.freedman_diaconis(data=diff_u2_c_pr, returnas="bins"), alpha=0.7, color='deepskyblue', zorder=5)
ax1_u2_c_pr.hist(diff_u2_c_ar, bins=utl.freedman_diaconis(data=diff_u2_c_ar, returnas="bins"), alpha=0.7, color='darkred',zorder=5)

plt.xlabel('Residuals')
plt.ylabel('Count')

plt.subplots_adjust(hspace = 0.3)
plt.savefig(path1 + '/Results/cal_us_and_evidance/u2_cla_r.pdf')
plt.close(fig_u2_c_pr)


#------------------------------------------------
print('--------How good are the Code LDCs?-------------')
#------------------------------------------------

u1_co_p, u2_co_p = np.loadtxt(path1 + '/Phoenix/code_limiting_LDC_pho.dat', usecols = (1,2), unpack = True)

u1max_co_p = np.max(u1_co_p)
u1min_co_p = np.min(u1_co_p)
u2max_co_p = np.max(u2_co_p)
u2min_co_p = np.min(u2_co_p)

u1_co_a, u2_co_a = np.loadtxt(path1 + '/Atlas/code_limiting_LDC_ata.dat', usecols = (1,2), unpack = True)

u1max_co_a = np.max(u1_co_a)
u1min_co_a = np.min(u1_co_a)
u2max_co_a = np.max(u2_co_a)
u2min_co_a = np.min(u2_co_a)

#-----------------------
print('#--------u1-------------')
#-----------------------

xlu1_co_p = np.minimum(u1min_j, u1min_co_p)
xuu1_co_p = np.maximum(u1max_j, u1max_co_p)

xlu1_co_a = np.minimum(u1min_j, u1min_co_a)
xuu1_co_a = np.maximum(u1max_j, u1max_co_a)

xlo = np.minimum(xlu1_co_p, xlu1_co_a)
xup = np.maximum(xuu1_co_p, xuu1_co_a)

x1u1_co_p = y1u1_co_p = np.linspace(xlo, xup, 100)
y11u1_co_p = np.zeros(len(x1u1_co_p))

diff_u1_co_p = np.array([])#--------------------------------------------------------------------------------------------
diff_u1_co_pe = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u1_j)):
	u11_co_p = np.random.normal(u1_co_p[i], 0, 10000)
	u11_j = np.random.normal(u1_j[i], u1_jp[i], 10000)
	diff1 = u11_co_p - u11_j
	u11_m = np.median(diff1)
	u11_e = np.std(diff1)
	diff_u1_co_p = np.hstack((diff_u1_co_p, u11_m))
	diff_u1_co_pe = np.hstack((diff_u1_co_pe, u11_e))

diff_u1_co_a = np.array([])#--------------------------------------------------------------------------------------------
diff_u1_co_ae = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u1_j)):
	u11_co_a = np.random.normal(u1_co_a[i], 0, 10000)
	u11_j = np.random.normal(u1_j[i], u1_jp[i], 10000)
	diff1 = u11_co_a - u11_j
	u11_m = np.median(diff1)
	u11_e = np.std(diff1)
	diff_u1_co_a = np.hstack((diff_u1_co_a, u11_m))
	diff_u1_co_ae = np.hstack((diff_u1_co_ae, u11_e))


fig_u1_co_p = plt.figure(figsize=(8,10))
gs_u1_co_p = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_u1_co_p = plt.subplot(gs_u1_co_p[0])

ax_u1_co_p.errorbar(u1_j, u1_co_p, xerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
ax_u1_co_p.errorbar(u1_j, u1_co_a, xerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

ax_u1_co_p.plot(x1u1_co_p, y1u1_co_p, 'k--')
ax_u1_co_p.grid()

ax_u1_co_p.set_xlim([xlo, xup])
ax_u1_co_p.set_ylim([xlo, xup])

plt.legend(loc='best')
plt.ylabel(r'$u_1$ (Theoretical)')
plt.xlabel(r'$u_1$ (Empirical)')
plt.title('Values from Espinoza \& Jordan(2015)')

ax1_u1_co_p = plt.subplot(gs_u1_co_p[1])#, sharex = ax_u1_co_p)

ax1_u1_co_p.hist(diff_u1_co_p, bins=utl.freedman_diaconis(data=diff_u1_co_p, returnas="bins"), alpha=0.7, color='orangered', zorder=5)
ax1_u1_co_p.hist(diff_u1_co_a, bins=utl.freedman_diaconis(data=diff_u1_co_a, returnas="bins"), alpha=0.7, color='cornflowerblue',zorder=5)

plt.xlabel('Residuals')
plt.ylabel('Count')

plt.subplots_adjust(hspace = 0.3)
plt.savefig(path1 + '/Results/cal_us_and_evidance/u1_code.pdf')
plt.close(fig_u1_co_p)

#-----------------------
print('#--------u2-------------')
#-----------------------

xlu2_co_p = np.minimum(u2min_j, u2min_co_p)
xuu2_co_p = np.maximum(u2max_j, u2max_co_p)

xlu2_co_a = np.minimum(u2min_j, u2min_co_a)
xuu2_co_a = np.maximum(u2max_j, u2max_co_a)

xlo = np.minimum(xlu2_co_p, xlu2_co_a)
xup = np.maximum(xuu2_co_p, xuu2_co_a)

x1u2_co_p = y1u2_co_p = np.linspace(xlo, xup, 100)
y11u2_co_p = np.zeros(len(x1u2_co_p))

diff_u2_co_p = np.array([])#--------------------------------------------------------------------------------------------
diff_u2_co_pe = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u2_j)):
	u22_co_p = np.random.normal(u2_co_p[i], 0, 10000)
	u22_j = np.random.normal(u2_j[i], u2_jp[i], 10000)
	diff2 = u22_co_p - u22_j
	u22_m = np.median(diff2)
	u22_e = np.std(diff2)
	diff_u2_co_p = np.hstack((diff_u2_co_p, u22_m))
	diff_u2_co_pe = np.hstack((diff_u2_co_pe, u22_e))

diff_u2_co_a = np.array([])#--------------------------------------------------------------------------------------------
diff_u2_co_ae = np.array([])#-------------------------------------------------------------------------------------------

for i in range(len(u2_j)):
	u22_co_a = np.random.normal(u2_co_a[i], 0, 10000)
	u22_j = np.random.normal(u2_j[i], u2_jp[i], 10000)
	diff2 = u22_co_a - u22_j
	u22_m = np.median(diff2)
	u22_e = np.std(diff2)
	diff_u2_co_a = np.hstack((diff_u2_co_a, u22_m))
	diff_u2_co_ae = np.hstack((diff_u2_co_ae, u22_e))


fig_u2_co_p = plt.figure(figsize=(8,10))
gs_u2_co_p = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_u2_co_p = plt.subplot(gs_u2_co_p[0])

ax_u2_co_p.errorbar(u2_j, u2_co_p, xerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
ax_u2_co_p.errorbar(u2_j, u2_co_a, xerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

ax_u2_co_p.plot(x1u2_co_p, y1u2_co_p, 'k--')
ax_u2_co_p.grid()

ax_u2_co_p.set_xlim([xlo, xup])
ax_u2_co_p.set_ylim([xlo, xup])

plt.legend(loc='best')
plt.ylabel(r'$u_2$ (Theoretical)')
plt.xlabel(r'$u_2$ (Empirical)')
plt.title('Values from Espinoza \& Jordan(2015)')

ax1_u2_co_p = plt.subplot(gs_u2_co_p[1])#, sharex = ax_u2_co_p)

ax1_u2_co_p.hist(diff_u2_co_p, bins=utl.freedman_diaconis(data=diff_u2_co_p, returnas="bins"), color='orangered', alpha=0.7, zorder=5)
ax1_u2_co_p.hist(diff_u2_co_a, bins=utl.freedman_diaconis(data=diff_u2_co_a, returnas="bins"), color='cornflowerblue', alpha=0.7, zorder=5)

plt.xlabel('Residuals')
plt.ylabel('Count')

plt.subplots_adjust(hspace = 0.3)
plt.savefig(path1 + '/Results/cal_us_and_evidance/u2_code.pdf')
plt.close(fig_u2_co_p)

print('--------------------------------------------------------------------------------------------------')
print('----------------------------Making Table of mean offset in LDCs-----------------------------------')
print('--------------------------------------------------------------------------------------------------')

f101 = open(path1 + '/Results/mean_off.dat', 'w')

print('-----------diff_u1_c_p')

diff = np.random.normal(0,0,10000)

for i in range(len(diff_u1_c_p)):
	diff1 = np.random.normal(diff_u1_c_p[i], diff_u1_c_pe[i], 10000)
	diff = np.vstack((diff, diff1))

diffi = diff[1:]
mean1 = np.mean(diffi, axis=0)

med1 = np.median(mean1)
std1 = np.std(mean1)

print('-----------diff_u2_c_p')

diff = np.random.normal(0,0,10000)

for i in range(len(diff_u2_c_p)):
	diff1 = np.random.normal(diff_u2_c_p[i], diff_u2_c_pe[i], 10000)
	diff = np.vstack((diff, diff1))

diffi = diff[1:]
mean2 = np.mean(diffi, axis=0)

med2 = np.median(mean2)
std2 = np.std(mean2)

print('-----------diff_u1_c_a')

diff = np.random.normal(0,0,10000)

for i in range(len(diff_u1_c_a)):
	diff1 = np.random.normal(diff_u1_c_a[i], diff_u1_c_ae[i], 10000)
	diff = np.vstack((diff, diff1))

diffi = diff[1:]
mean3 = np.mean(diffi, axis=0)

med3 = np.median(mean3)
std3 = np.std(mean3)

print('-----------diff_u2_c_a')

diff = np.random.normal(0,0,10000)

for i in range(len(diff_u2_c_a)):
	diff1 = np.random.normal(diff_u2_c_a[i], diff_u2_c_ae[i], 10000)
	diff = np.vstack((diff, diff1))

diffi = diff[1:]
mean4 = np.mean(diffi, axis=0)

med4 = np.median(mean4)
std4 = np.std(mean4)

print('-----------diff_u1_c_ar')

diff = np.random.normal(0,0,10000)

for i in range(len(diff_u1_c_ar)):
	diff1 = np.random.normal(diff_u1_c_ar[i], diff_u1_c_aer[i], 10000)
	diff = np.vstack((diff, diff1))

diffi = diff[1:]
mean3r = np.mean(diffi, axis=0)

med3r = np.median(mean3r)
std3r = np.std(mean3r)

print('-----------diff_u2_c_ar')

diff = np.random.normal(0,0,10000)

for i in range(len(diff_u2_c_ar)):
	diff1 = np.random.normal(diff_u2_c_ar[i], diff_u2_c_aer[i], 10000)
	diff = np.vstack((diff, diff1))

diffi = diff[1:]
mean4r = np.mean(diffi, axis=0)

med4r = np.median(mean4r)
std4r = np.std(mean4r)

###############################

print('-----------diff_u1_co_p')

diff = np.random.normal(0,0,10000)

for i in range(len(diff_u1_co_p)):
	diff1 = np.random.normal(diff_u1_co_p[i], diff_u1_co_pe[i], 10000)
	diff = np.vstack((diff, diff1))

diffi = diff[1:]
mean5 = np.mean(diffi, axis=0)

med5 = np.median(mean5)
std5 = np.std(mean5)

print('-----------diff_u2_co_p')

diff = np.random.normal(0,0,10000)

for i in range(len(diff_u2_co_p)):
	diff1 = np.random.normal(diff_u2_co_p[i], diff_u2_co_pe[i], 10000)
	diff = np.vstack((diff, diff1))

diffi = diff[1:]
mean6 = np.mean(diffi, axis=0)

med6 = np.median(mean6)
std6 = np.std(mean6)

print('-----------diff_u1_co_a')

diff = np.random.normal(0,0,10000)

for i in range(len(diff_u1_co_a)):
	diff1 = np.random.normal(diff_u1_co_a[i], diff_u1_co_ae[i], 10000)
	diff = np.vstack((diff, diff1))

diffi = diff[1:]
mean7 = np.mean(diffi, axis=0)

med7 = np.median(mean7)
std7 = np.std(mean7)

print('-----------diff_u2_co_a')

diff = np.random.normal(0,0,10000)

for i in range(len(diff_u2_co_a)):
	diff1 = np.random.normal(diff_u2_co_a[i], diff_u2_co_ae[i], 10000)
	diff = np.vstack((diff, diff1))

diffi = diff[1:]
mean8 = np.mean(diffi, axis=0)

med8 = np.median(mean8)
std8 = np.std(mean8)

f101.write('\t\t\t\tu1\t\t\t\tu2\n')
f101.write('Claret(2017), PHOENIX - s method\t\t' + str(med1) + ' +/- ' + str(std1) + '\t\t' + str(med2) + ' +/- ' + str(std2) + '\n')
f101.write('Claret(2017), PHOENIX - r method\t\t' + str(med3r) + ' +/- ' + str(std3r) + '\t\t' + str(med4r) + ' +/- ' + str(std4r) + '\n')
f101.write('Claret(2017), ATLAS\t\t' + str(med3) + ' +/- ' + str(std3) + '\t\t' + str(med4) + ' +/- ' + str(std4) + '\n')
f101.write('EJ(2015), PHOENIX\t\t' + str(med5) + ' +/- ' + str(std5) + '\t\t' + str(med6) + ' +/- ' + str(std6) + '\n')
f101.write('EJ(2015), ATLAS\t\t\t' + str(med7) + ' +/- ' + str(std7) + '\t\t' + str(med8) + ' +/- ' + str(std8) + '\n')
f101.close()

print('--------------------------------------------------------------------------------------------------')
print('-------------------------------Plots of LDCs with temperature variation---------------------------')
print('--------------------------------------------------------------------------------------------------')

tmin = np.min(teff)
tmax = np.max(teff)

t11 = np.linspace(tmin, tmax, 1000)

x = np.linspace(tmin, tmax, 100)
y = np.zeros(len(x))

def line(x,m,c):
	function = m*x + c
	return function

def constant(x, c):
	function = c + x*0
	return function

def quadratic(x, a, b, c):
	function = a*x*x + b*x + c
	return function

model = ['constant', 'linear', 'quadratic']

#------------------------------------------------
print('---------u1 - Claret(2017) ---------------------')
#------------------------------------------------

fig1 = plt.figure(figsize=(12,10))

plt.errorbar(teff, diff_u1_c_p, yerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
plt.errorbar(teff, diff_u1_c_a, yerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

#------------------------------------------------------------#
#abc1 = np.vstack((teff, diff_u1_c_p))                       #
#abc = np.transpose(np.vstack((abc1, u1_jp)))                #
#np.savetxt(path1 + '/1.dat', abc)                           #
#------------------------------------------------------------#

f111 = open(path1 + '/Results/variation_with_temp/Off/off_u1_cla_p.dat', 'w')
for i in range(len(teff)):
	ddd = 5*(u1_jp[i]) - np.abs(diff_u1_c_p[i])
	if ddd<0:
		f111.write(name[i] + '\t' + str(u1_j[i]) + '+-' + str(u1_jp[i]) + '\t' + str(np.abs(diff_u1_c_p[i])) + '\n')

f111.close()

f222 = open(path1 + '/Results/variation_with_temp/Off/off_u1_cla_a.dat', 'w')
for i in range(len(teff)):
	ddd = 5*(u1_jp[i]) - np.abs(diff_u1_c_a[i])
	if ddd<0:
		f222.write(name[i] + '\t' + str(u1_j[i]) + '+-' + str(u1_jp[i]) + '\t' + str(np.abs(diff_u1_c_a[i])) + '\n')

f222.close()

#plt.plot(teff, np.linspace(med1, med1, len(teff)), 'red', teff, np.linspace(med3, med3, len(teff)), 'blue')
#plt.fill_between(teff, med1 - std1, med1 + std1, color = 'red', alpha = 0.5)
#plt.fill_between(teff, med3 + std3, med3 - std3, color = 'blue', alpha = 0.5)

#----Constant fitting
popt_c, pcov_c = cft(constant, teff, diff_u1_c_p)
popt1_c, pcov1_c = cft(constant, teff, diff_u1_c_a)

rss_u1_c_p_c = rss_u1_c_a_c = 0
for i in range(len(teff)):
	r111 = (diff_u1_c_p[i] - constant(teff[i], *popt_c))**2
	rss_u1_c_p_c = rss_u1_c_p_c + r111
	r222 = (diff_u1_c_a[i] - constant(teff[i], *popt1_c))**2
	rss_u1_c_a_c = rss_u1_c_a_c + r222

bic_u1_c_p_c = len(diff_u1_c_p)*np.log((rss_u1_c_p_c)/(len(diff_u1_c_p))) + np.log(len(diff_u1_c_p))
bic_u1_c_a_c = len(diff_u1_c_a)*np.log((rss_u1_c_a_c)/(len(diff_u1_c_a))) + np.log(len(diff_u1_c_a))

#----Linear fitting
popt_l, pcov_l = cft(line, teff, diff_u1_c_p)
popt1_l, pcov1_l = cft(line, teff, diff_u1_c_a)

rss_u1_c_p_l = rss_u1_c_a_l = 0
for i in range(len(teff)):
	r111 = (diff_u1_c_p[i] - line(teff[i], *popt_l))**2
	rss_u1_c_p_l = rss_u1_c_p_l + r111
	r222 = (diff_u1_c_a[i] - line(teff[i], *popt1_l))**2
	rss_u1_c_a_l = rss_u1_c_a_l + r222

bic_u1_c_p_l = len(diff_u1_c_p)*np.log((rss_u1_c_p_l)/(len(diff_u1_c_p))) + 2*np.log(len(diff_u1_c_p))
bic_u1_c_a_l = len(diff_u1_c_a)*np.log((rss_u1_c_a_l)/(len(diff_u1_c_a))) + 2*np.log(len(diff_u1_c_a))

#----Quadratic fitting
popt_q, pcov_q = cft(quadratic, teff, diff_u1_c_p)
popt1_q, pcov1_q = cft(quadratic, teff, diff_u1_c_a)

rss_u1_c_p_q = rss_u1_c_a_q = 0
for i in range(len(teff)):
	r111 = (diff_u1_c_p[i] - quadratic(teff[i], *popt_q))**2
	rss_u1_c_p_q = rss_u1_c_p_q + r111
	r222 = (diff_u1_c_a[i] - quadratic(teff[i], *popt1_q))**2
	rss_u1_c_a_q = rss_u1_c_a_q + r222

bic_u1_c_p_q = len(diff_u1_c_p)*np.log((rss_u1_c_p_q)/(len(diff_u1_c_p))) + 3*np.log(len(diff_u1_c_p))
bic_u1_c_a_q = len(diff_u1_c_a)*np.log((rss_u1_c_a_q)/(len(diff_u1_c_a))) + 3*np.log(len(diff_u1_c_a))

#----Plotting models
bic_u1_c_p = [bic_u1_c_p_c, bic_u1_c_p_l, bic_u1_c_p_q]
bic_u1_c_a = [bic_u1_c_a_c, bic_u1_c_a_l, bic_u1_c_a_q]

dict_bic_u1_c_p = {}
dict_bic_u1_c_a = {}
for i in range(3):
	dict_bic_u1_c_p[model[i]] = bic_u1_c_p[i]
	dict_bic_u1_c_a[model[i]] = bic_u1_c_a[i]

best_bic_u1_c_p = utl.lowest_bic(dict_bic_u1_c_p)
if 'constant' in best_bic_u1_c_p:
	plt.plot(t11, constant(t11, *popt_c), color = 'orangered', ls = '-.')
elif 'linear' in best_bic_u1_c_p:
	plt.plot(t11, line(t11, *popt_l), color='orangered', ls='-.')
elif 'quadratic' in best_bic_u1_c_p:
	plt.plot(t11, quadratic(t11, *popt_q), color = 'orangered', ls='-.')

best_bic_u1_c_a = utl.lowest_bic(dict_bic_u1_c_a)
if 'constant' in best_bic_u1_c_a:
	plt.plot(t11, constant(t11, *popt1_c), color = 'cornflowerblue', ls = '-.')
elif 'linear' in best_bic_u1_c_a:
	plt.plot(t11, line(t11, *popt1_l), color='cornflowerblue', ls='-.')
elif 'quadratic' in best_bic_u1_c_a:
	plt.plot(t11, quadratic(t11, *popt1_q), color = 'cornflowerblue', ls='-.')

#---------------------
plt.plot(x, y, 'k--')
plt.grid()
plt.legend(loc='best')
plt.ylabel(r'Residuals in $u_1$ Claret(2017) LDCs')
plt.xlabel('Effective Temperature')
plt.savefig(path1 + '/Results/variation_with_temp/u1_cla_te.pdf')
plt.close(fig1)

#------------------------------------------------
print('---------u2 - Claret(2017) ---------------------')
#------------------------------------------------

fig2 = plt.figure(figsize=(12,10))

plt.errorbar(teff, diff_u2_c_p, yerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
plt.errorbar(teff, diff_u2_c_a, yerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

#plt.plot(teff, med2, 'red', teff, med4, 'blue')
#plt.fill_between(teff, med2 + std2, med2 - std2, color = 'red', alpha = 0.5)
#plt.fill_between(teff, med4 + std4, med4 - std4, color = 'blue', alpha = 0.5)

f333 = open(path1 + '/Results/variation_with_temp/Off/off_u2_cla_p.dat', 'w')
for i in range(len(teff)):
	ddd = 5*(u2_jp[i]) - np.abs(diff_u2_c_p[i])
	if ddd<0:
		f333.write(name[i] + '\t' + str(u2_j[i]) + '+-' + str(u2_jp[i]) + '\t' + str(np.abs(diff_u2_c_p[i])) + '\n')

f333.close()

f444 = open(path1 + '/Results/variation_with_temp/Off/off_u2_cla_a.dat', 'w')
for i in range(len(teff)):
	ddd = 5*(u2_jp[i]) - np.abs(diff_u2_c_a[i])
	if ddd<0:
		f444.write(name[i] + '\t' + str(u2_j[i]) + '+-' + str(u2_jp[i]) + '\t' + str(np.abs(diff_u2_c_a[i])) + '\n')

f444.close()

#----Constant fitting
popt2_c, pcov2_c = cft(constant, teff, diff_u2_c_p)
popt3_c, pcov3_c = cft(constant, teff, diff_u2_c_a)

rss_u2_c_p_c = rss_u2_c_a_c = 0
for i in range(len(teff)):
	r111 = (diff_u2_c_p[i] - constant(teff[i], *popt2_c))**2
	rss_u2_c_p_c = rss_u2_c_p_c + r111
	r222 = (diff_u2_c_a[i] - constant(teff[i], *popt3_c))**2
	rss_u2_c_a_c = rss_u2_c_a_c + r222

bic_u2_c_p_c = len(diff_u2_c_p)*np.log((rss_u2_c_p_c)/(len(diff_u2_c_p))) + np.log(len(diff_u2_c_p))
bic_u2_c_a_c = len(diff_u2_c_a)*np.log((rss_u2_c_a_c)/(len(diff_u2_c_a))) + np.log(len(diff_u2_c_a))

#----Linear fitting
popt2_l, pcov2_l = cft(line, teff, diff_u2_c_p)
popt3_l, pcov3_l = cft(line, teff, diff_u2_c_a)

rss_u2_c_p_l = rss_u2_c_a_l = 0
for i in range(len(teff)):
	r111 = (diff_u2_c_p[i] - line(teff[i], *popt2_l))**2
	rss_u2_c_p_l = rss_u2_c_p_l + r111
	r222 = (diff_u2_c_a[i] - line(teff[i], *popt3_l))**2
	rss_u2_c_a_l = rss_u2_c_a_l + r222

bic_u2_c_p_l = len(diff_u2_c_p)*np.log((rss_u2_c_p_l)/(len(diff_u2_c_p))) + 2*np.log(len(diff_u2_c_p))
bic_u2_c_a_l = len(diff_u2_c_a)*np.log((rss_u2_c_a_l)/(len(diff_u2_c_a))) + 2*np.log(len(diff_u2_c_a))

#----Quadratic fitting
popt2_q, pcov2_q = cft(quadratic, teff, diff_u2_c_p)
popt3_q, pcov3_q = cft(quadratic, teff, diff_u2_c_a)

rss_u2_c_p_q = rss_u2_c_a_q = 0
for i in range(len(teff)):
	r111 = (diff_u2_c_p[i] - quadratic(teff[i], *popt2_q))**2
	rss_u2_c_p_q = rss_u2_c_p_q + r111
	r222 = (diff_u2_c_a[i] - quadratic(teff[i], *popt3_q))**2
	rss_u2_c_a_q = rss_u2_c_a_q + r222

bic_u2_c_p_q = len(diff_u2_c_p)*np.log((rss_u2_c_p_q)/(len(diff_u2_c_p))) + 3*np.log(len(diff_u2_c_p))
bic_u2_c_a_q = len(diff_u2_c_a)*np.log((rss_u2_c_a_q)/(len(diff_u2_c_a))) + 3*np.log(len(diff_u2_c_a))

#----Plotting models
bic_u2_c_p = [bic_u2_c_p_c, bic_u2_c_p_l, bic_u2_c_p_q]
bic_u2_c_a = [bic_u2_c_a_c, bic_u2_c_a_l, bic_u2_c_a_q]

dict_bic_u2_c_p = {}
dict_bic_u2_c_a = {}
for i in range(3):
	dict_bic_u2_c_p[model[i]] = bic_u2_c_p[i]
	dict_bic_u2_c_a[model[i]] = bic_u2_c_a[i]

best_bic_u2_c_p = utl.lowest_bic(dict_bic_u2_c_p)
if 'constant' in best_bic_u2_c_p:
	plt.plot(t11, constant(t11, *popt2_c), color='orangered', ls='-.')
elif 'linear' in best_bic_u2_c_p:
	plt.plot(t11, line(t11, *popt2_l), color='orangered',ls='-.')
elif 'quadratic' in best_bic_u2_c_p:
	plt.plot(t11, quadratic(t11, *popt2_q), color='orangered', ls='-.')

best_bic_u2_c_a = utl.lowest_bic(dict_bic_u2_c_a)
if 'constant' in best_bic_u2_c_a:
	plt.plot(t11, constant(t11, *popt3_c), color='cornflowerblue', ls='-.')
elif 'linear' in best_bic_u2_c_a:
	plt.plot(t11, line(t11, *popt3_l), color='cornflowerblue',ls='-.')
elif 'quadratic' in best_bic_u2_c_a:
	plt.plot(t11, quadratic(t11, *popt3_q), color='cornflowerblue', ls='-.')

#---------------------
plt.plot(x, y, 'k--')
plt.grid()
plt.legend(loc='best')
plt.ylabel(r'Residuals in $u_2$ Claret(2017) LDCs')
plt.xlabel('Effective Temperature')
plt.savefig(path1 + '/Results/variation_with_temp/u2_cla_te.pdf')
plt.close(fig2)

#------------------------------------------------
print('-----------u1 - Code ---------------------------')
#------------------------------------------------

fig5 = plt.figure(figsize=(12,10))

plt.errorbar(teff, diff_u1_co_p, yerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
plt.errorbar(teff, diff_u1_co_a, yerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

#plt.plot(teff, med5, 'red', teff, med7, 'blue')
#plt.fill_between(teff, med5 + std5, med5 - std5, color = 'red', alpha = 0.5)
#plt.fill_between(teff, med7 + std7, med7 - std7, color = 'blue', alpha = 0.5)

f555 = open(path1 + '/Results/variation_with_temp/Off/off_u1_co_p.dat', 'w')
for i in range(len(teff)):
	ddd = 5*(u1_jp[i]) - np.abs(diff_u1_co_p[i])
	if ddd<0:
		f555.write(name[i] + '\t' + str(u1_j[i]) + '+-' + str(u1_jp[i]) + '\t' + str(np.abs(diff_u1_co_p[i])) + '\n')

f555.close()

f666 = open(path1 + '/Results/variation_with_temp/Off/off_u1_co_a.dat', 'w')
for i in range(len(teff)):
	ddd = 5*(u1_jp[i]) - np.abs(diff_u1_co_a[i])
	if ddd<0:
		f666.write(name[i] + '\t' + str(u1_j[i]) + '+-' + str(u1_jp[i]) + '\t' + str(np.abs(diff_u1_co_a[i])) + '\n')

f666.close()

#----Constant fitting
popt4_c, pcov4_c = cft(constant, teff, diff_u1_co_p)
popt5_c, pcov5_c = cft(constant, teff, diff_u1_co_a)

rss_u1_co_p_c = rss_u1_co_a_c = 0
for i in range(len(teff)):
	r111 = (diff_u1_co_p[i] - constant(teff[i], *popt4_c))**2
	rss_u1_co_p_c = rss_u1_co_p_c + r111
	r222 = (diff_u1_co_a[i] - constant(teff[i], *popt5_c))**2
	rss_u1_co_a_c = rss_u1_co_a_c + r222

bic_u1_co_p_c = len(diff_u1_co_p)*np.log((rss_u1_co_p_c)/(len(diff_u1_co_p))) + np.log(len(diff_u1_co_p))
bic_u1_co_a_c = len(diff_u1_co_a)*np.log((rss_u1_co_a_c)/(len(diff_u1_co_a))) + np.log(len(diff_u1_co_a))

#----Linear fitting
popt4_l, pcov4_l = cft(line, teff, diff_u1_co_p)
popt5_l, pcov5_l = cft(line, teff, diff_u1_co_a)

rss_u1_co_p_l = rss_u1_co_a_l = 0
for i in range(len(teff)):
	r111 = (diff_u1_co_p[i] - line(teff[i], *popt4_l))**2
	rss_u1_co_p_l = rss_u1_co_p_l + r111
	r222 = (diff_u1_co_a[i] - line(teff[i], *popt5_l))**2
	rss_u1_co_a_l = rss_u1_co_a_l + r222

bic_u1_co_p_l = len(diff_u1_co_p)*np.log((rss_u1_co_p_l)/(len(diff_u1_co_p))) + 2*np.log(len(diff_u1_co_p))
bic_u1_co_a_l = len(diff_u1_co_a)*np.log((rss_u1_co_a_l)/(len(diff_u1_co_a))) + 2*np.log(len(diff_u1_co_a))

#----Quadratic fitting
popt4_q, pcov4_q = cft(quadratic, teff, diff_u1_co_p)
popt5_q, pcov5_q = cft(quadratic, teff, diff_u1_co_a)

rss_u1_co_p_q = rss_u1_co_a_q = 0
for i in range(len(teff)):
	r111 = (diff_u1_co_p[i] - quadratic(teff[i], *popt4_q))**2
	rss_u1_co_p_q = rss_u1_co_p_q + r111
	r222 = (diff_u1_co_a[i] - quadratic(teff[i], *popt5_q))**2
	rss_u1_co_a_q = rss_u1_co_a_q + r222

bic_u1_co_p_q = len(diff_u1_co_p)*np.log((rss_u1_co_p_q)/(len(diff_u1_co_p))) + 3*np.log(len(diff_u1_co_p))
bic_u1_co_a_q = len(diff_u1_co_a)*np.log((rss_u1_co_a_q)/(len(diff_u1_co_a))) + 3*np.log(len(diff_u1_co_a))

#----Plotting models
bic_u1_co_p = [bic_u1_co_p_c, bic_u1_co_p_l, bic_u1_co_p_q]
bic_u1_co_a = [bic_u1_co_a_c, bic_u1_co_a_l, bic_u1_co_a_q]

dict_bic_u1_co_p = {}
dict_bic_u1_co_a = {}
for i in range(3):
	dict_bic_u1_co_p[model[i]] = bic_u1_co_p[i]
	dict_bic_u1_co_a[model[i]] = bic_u1_co_a[i]

best_bic_u1_co_p = utl.lowest_bic(dict_bic_u1_co_p)
if 'constant' in best_bic_u1_co_p:
	plt.plot(t11, constant(t11, *popt4_c), color='orangered', ls='-.')
if 'linear' in best_bic_u1_co_p:
	plt.plot(t11, line(t11, *popt4_l), color='orangered',ls='-.')
if 'quadratic' in best_bic_u1_co_p:
	plt.plot(t11, quadratic(t11, *popt4_q), color='orangered', ls='-.')

best_bic_u1_co_a = utl.lowest_bic(dict_bic_u1_co_a)
if 'constant' in best_bic_u1_co_a:
	plt.plot(t11, constant(t11, *popt5_c), color='cornflowerblue', ls='-.')
if 'linear' in best_bic_u1_co_a:
	plt.plot(t11, line(t11, *popt5_l), color='cornflowerblue',ls='-.')
if 'quadratic' in best_bic_u1_co_a:
	plt.plot(t11, quadratic(t11, *popt5_q), color='cornflowerblue', ls='-.')

#---------------------
plt.plot(x, y, 'k--')
plt.grid()
plt.legend(loc='best')
plt.ylabel(r'Residuals in $u_1$ - Espinoza \& Jordan LDCs')
plt.xlabel('Effective Temperature')
plt.savefig(path1 + '/Results/variation_with_temp/u1_code_te.pdf')
plt.close(fig5)

#------------------------------------------------
print('------------u2 - Code ------------------')
#------------------------------------------------

fig6 = plt.figure(figsize=(12,10))

plt.errorbar(teff, diff_u2_co_p, yerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
plt.errorbar(teff, diff_u2_co_a, yerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

#plt.plot(teff, med6, 'red', teff, med8, 'blue')
#plt.fill_between(teff, med8 + std8, med8 - std8, color = 'blue', alpha = 0.5)
#plt.fill_between(teff, med6 + std6, med6 - std6, color = 'red', alpha = 0.5)

f777 = open(path1 + '/Results/variation_with_temp/Off/off_u2_co_p.dat', 'w')
for i in range(len(teff)):
	ddd = 5*(u2_jp[i]) - np.abs(diff_u2_co_p[i])
	if ddd<0:
		f777.write(name[i] + '\t' + str(u2_j[i]) + '+-' + str(u2_jp[i]) + '\t' + str(np.abs(diff_u2_co_p[i])) + '\n')

f777.close()

f888 = open(path1 + '/Results/variation_with_temp/Off/off_u2_co_a.dat', 'w')
for i in range(len(teff)):
	ddd = 5*(u2_jp[i]) - np.abs(diff_u2_co_a[i])
	if ddd<0:
		f888.write(name[i] + '\t' + str(u2_j[i]) + '+-' + str(u2_jp[i]) + '\t' + str(np.abs(diff_u2_co_a[i])) + '\n')

f888.close()

#----Constant fitting
popt6_c, pcov6_c = cft(constant, teff, diff_u2_co_p)
popt7_c, pcov7_c = cft(constant, teff, diff_u2_co_a)

rss_u2_co_p_c = rss_u2_co_a_c = 0
for i in range(len(teff)):
	r111 = (diff_u2_co_p[i] - constant(teff[i], *popt6_c))**2
	rss_u2_co_p_c = rss_u2_co_p_c + r111
	r222 = (diff_u2_co_a[i] - constant(teff[i], *popt7_c))**2
	rss_u2_co_a_c = rss_u2_co_a_c + r222

bic_u2_co_p_c = len(diff_u2_co_p)*np.log((rss_u2_co_p_c)/(len(diff_u2_co_p))) + np.log(len(diff_u2_co_p))
bic_u2_co_a_c = len(diff_u2_co_a)*np.log((rss_u2_co_a_c)/(len(diff_u2_co_a))) + np.log(len(diff_u2_co_a))

#----Linear fitting
popt6_l, pcov6_l = cft(line, teff, diff_u2_co_p)
popt7_l, pcov7_l = cft(line, teff, diff_u2_co_a)

rss_u2_co_p_l = rss_u2_co_a_l = 0
for i in range(len(teff)):
	r111 = (diff_u2_co_p[i] - line(teff[i], *popt6_l))**2
	rss_u2_co_p_l = rss_u2_co_p_l + r111
	r222 = (diff_u2_co_a[i] - line(teff[i], *popt7_l))**2
	rss_u2_co_a_l = rss_u2_co_a_l + r222

bic_u2_co_p_l = len(diff_u2_co_p)*np.log((rss_u2_co_p_l)/(len(diff_u2_co_p))) + 2*np.log(len(diff_u2_co_p))
bic_u2_co_a_l = len(diff_u2_co_a)*np.log((rss_u2_co_a_l)/(len(diff_u2_co_a))) + 2*np.log(len(diff_u2_co_a))

#----Quadratic fitting
popt6_q, pcov6_q = cft(quadratic, teff, diff_u2_co_p)
popt7_q, pcov7_q = cft(quadratic, teff, diff_u2_co_a)

rss_u2_co_p_q = rss_u2_co_a_q = 0
for i in range(len(teff)):
	r111 = (diff_u2_co_p[i] - quadratic(teff[i], *popt6_q))**2
	rss_u2_co_p_q = rss_u2_co_p_q + r111
	r222 = (diff_u2_co_a[i] - quadratic(teff[i], *popt7_q))**2
	rss_u2_co_a_q = rss_u2_co_a_q + r222

bic_u2_co_p_q = len(diff_u2_co_p)*np.log((rss_u2_co_p_q)/(len(diff_u2_co_p))) + 3*np.log(len(diff_u2_co_p))
bic_u2_co_a_q = len(diff_u2_co_a)*np.log((rss_u2_co_a_q)/(len(diff_u2_co_a))) + 3*np.log(len(diff_u2_co_a))

#----Plotting Models
bic_u2_co_p = [bic_u2_co_p_c, bic_u2_co_p_l, bic_u2_co_p_q]
bic_u2_co_a = [bic_u2_co_a_c, bic_u2_co_a_l, bic_u2_co_a_q]

dict_bic_u2_co_p = {}
dict_bic_u2_co_a = {}
for i in range(3):
	dict_bic_u2_co_p[model[i]] = bic_u2_co_p[i]
	dict_bic_u2_co_a[model[i]] = bic_u2_co_a[i]

best_bic_u2_co_p = utl.lowest_bic(dict_bic_u2_co_p)
if 'constant' in best_bic_u2_co_p:
	plt.plot(t11, constant(t11, *popt6_c), color='orangered', ls='-.')
elif 'linear' in best_bic_u2_co_p:
	plt.plot(t11, line(t11, *popt6_l), color='orangered', ls='-.')
elif 'quadratic' in best_bic_u2_co_p:
	plt.plot(t11, quadratic(t11, *popt6_q), color='orangered', ls='-.')

best_bic_u2_co_a = utl.lowest_bic(dict_bic_u2_co_a)
if 'constant' in best_bic_u2_co_a:
	plt.plot(t11, constant(t11, *popt7_c), color='cornflowerblue', ls='-.')
elif 'linear' in best_bic_u2_co_a:
	plt.plot(t11, line(t11, *popt7_l), color='cornflowerblue', ls='-.')
elif 'quadratic' in best_bic_u2_co_a:
	plt.plot(t11, quadratic(t11, *popt7_q), color='cornflowerblue', ls='-.')

#--------------------
plt.plot(x, y, 'k--')
plt.grid()
plt.legend(loc='best')
plt.ylabel(r'Residuals in $u_2$ Espinoza \& Jordan (2015) LDCs')
plt.xlabel('Effective Temperature')
plt.savefig(path1 + '/Results/variation_with_temp/u2_code_te.pdf')
plt.close(fig6)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
f102 = open(path1 + '/Results/variation_with_temp/bic.dat', 'w')
f102.write('\t\t\t\tu1\t\t\t\tu2\n')
f102.write('-------------------------------------------------------------------------------------------\n')
f102.write('                                                                                           \n')
f102.write('-------------------------------------------------------------------------------------------\n')
f102.write('Claret(2017), PHOENIX-s\t\t' + 'Constant fit  \t' + str(bic_u1_c_p_c) + '\t' + str(bic_u2_c_p_c) + '\n')
f102.write('                       \t\t' + 'Linear fit    \t' + str(bic_u1_c_p_l) + '\t' + str(bic_u2_c_p_l) + '\n')
f102.write('                       \t\t' + 'Quadratic fit \t' + str(bic_u1_c_p_q) + '\t' + str(bic_u2_c_p_q) + '\n')
f102.write('-------------------------------------------------------------------------------------------\n')
f102.write('-------------------------------------------------------------------------------------------\n')
f102.write('Claret(2017), ATLAS    \t\t' + 'Constant fit  \t' + str(bic_u1_c_a_c) + '\t' + str(bic_u2_c_a_c) + '\n')
f102.write('                       \t\t' + 'Linear fit    \t' + str(bic_u1_c_a_l) + '\t' + str(bic_u2_c_a_l) + '\n')
f102.write('                       \t\t' + 'Quadratic fit \t' + str(bic_u1_c_a_q) + '\t' + str(bic_u2_c_a_q) + '\n')
f102.write('-------------------------------------------------------------------------------------------\n')
f102.write('-------------------------------------------------------------------------------------------\n')
f102.write('EJ(2015), PHOENIX      \t\t' + 'Constant fit  \t' + str(bic_u1_co_p_c) + '\t' + str(bic_u2_co_p_c) + '\n')
f102.write('                       \t\t' + 'Linear fit    \t' + str(bic_u1_co_p_l) + '\t' + str(bic_u2_co_p_l) + '\n')
f102.write('                       \t\t' + 'Quadratic fit \t' + str(bic_u1_co_p_q) + '\t' + str(bic_u2_co_p_q) + '\n')
f102.write('-------------------------------------------------------------------------------------------\n')
f102.write('-------------------------------------------------------------------------------------------\n')
f102.write('EJ(2015), ATLAS        \t\t' + 'Constant fit  \t' + str(bic_u1_co_a_c) + '\t' + str(bic_u2_co_a_c) + '\n')
f102.write('                       \t\t' + 'Linear fit    \t' + str(bic_u1_co_a_l) + '\t' + str(bic_u2_co_a_l) + '\n')
f102.write('                       \t\t' + 'Quadratic fit \t' + str(bic_u1_co_a_q) + '\t' + str(bic_u2_co_a_q) + '\n')
f102.close()

print("--------------------------------------------------------------------------------")
print("---------------------Your ploting task is complete!!----------------------------")
print("--------------------------------------------------------------------------------")
