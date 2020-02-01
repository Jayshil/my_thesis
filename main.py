import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec as gd
import os
import pickle as pc
import utils1 as utl
from exotoolbox import utils
from astropy.io import fits
from astroquery.mast import Observations as obs
from scipy import interpolate as inp
import matplotlib.cm as cm
import matplotlib.colors as cls
from scipy.optimize import curve_fit as cft
from pylab import *
import seaborn as sns


path1 = os.getcwd()#input('Enter the path of this folder: ')
#path1 = '/home/jayshil/Documents/Dissertation'

os.system('mkdir ' + path1 + '/Light-curve')

os.system('mkdir ' + path1 + '/Results')
os.system('mkdir ' + path1 + '/Results/cal_us_and_evidance')
os.system('mkdir ' + path1 + '/Results/comp_a_r_p')
os.system('mkdir ' + path1 + '/Results/stellar_prop')
os.system('mkdir ' + path1 + '/Results/variation_with_temp')
os.system('mkdir ' + path1 + '/limb-darkening-master/results')

#--------------------------------------------------------------------------------------------------
#--------------------------------Claret (2017) PHOENIX LDCs----------------------------------------
#--------------------------------------------------------------------------------------------------
f1 = open(path1 + '/Phoenix/claret_us_nl_pho.dat','w')#---------------------Non-linear-------------
f1.write('#Name\t\tc1\t\t\tc2\t\t\tc3\t\t\tc4\n')

f2 = open(path1 + '/Phoenix/claret_limiting_LDC_pho.dat','w')#-------------Limiting----------------
f2.write('#Name\t\tu1\t\tu2\n')

#--------------------------------------------------------------------------------------------------
#--------------------------------Claret (2017) ATLAS LDCs------------------------------------------
#--------------------------------------------------------------------------------------------------
f3 = open(path1 + '/Atlas/claret_us_nl_ata.dat','w')#---------------------Non-linear---------------
f3.write('#Name\t\tc1\t\t\tc2\t\t\tc3\t\t\tc4\n')

f4 = open(path1 + '/Atlas/claret_limiting_LDC_ata.dat','w')#-------------Limiting------------------
f4.write('#Name\t\tu1\t\tu2\n')

#--------------------------------------------------------------------------------------------------
#----------------------------------Code ATLAS LDCs-------------------------------------------------
#--------------------------------------------------------------------------------------------------
f33 = open(path1 + '/Atlas/code_us_nl_ata.dat','w')#-------------------Non-linear------------------
f33.write('#Name\t\tc1\t\t\tc2\t\t\tc3\t\t\tc4\n')

f44 = open(path1 + '/Atlas/code_limiting_LDC_ata.dat','w')#-----------------Limiting---------------
f44.write('#Name\t\tu1\t\t\t\tu2\n')

#--------------------------------------------------------------------------------------------------
#-----------------------------------Code PHOENIX LDCs----------------------------------------------
#--------------------------------------------------------------------------------------------------
f11 = open(path1 + '/Phoenix/code_us_nl_pho.dat','w')#----------------------Non-linear-------------
f11.write('#Name\t\tc1\t\t\tc2\t\t\tc3\t\t\tc4\n')

f22 = open(path1 + '/Phoenix/code_limiting_LDC_pho.dat','w')#--------------Limiting----------------
f22.write('#Name\t\tu1\t\tu2\n')

#--------------------------------------------------------------------------------------------------
#--------------------------------For a/R*, period and Rp/R*----------------------------------------
#--------------------------------------------------------------------------------------------------
f7 = open(path1 + '/Results/comp_a_r_p/aste.dat','w')
f7.write('#Name\ta_jul\t+err_jul\t-err_jul\n')

f8 = open(path1 + '/Results/comp_a_r_p/period.dat','w')
f8.write('#Name\tP_jul\t+err_jul\t-err_jul\n')

f9 = open(path1 + '/Results/comp_a_r_p/rprste.dat','w')
f9.write('#Name\tp_jul\t+err_jul\t-err_jul\n')

f10 = open(path1 + '/Results/comp_a_r_p/tc.dat','w')
f10.write('#Name\ttc_jul\t+err_jul\t-err_jul\n')

f19 = open(path1 + '/Results/cal_us_and_evidance/cal_u1_u2.dat','w')
f19.write('#Name\t\tu1\t\t\t+err\t\t\t-err\t\t\tu2\t\t\t+err\t\t\t-err\n')

f1_to = open(path1 + '/Results/comp_a_r_p/to_the.dat', 'w')

#---------------------------------------------------------------------------------------------
#--------------------------Taking Data from the data file-------------------------------------
#---------------------------------------------------------------------------------------------

name = np.loadtxt('data.dat', dtype = str, usecols = 0, unpack = True)
teff, lg, mh, vturb, p, pperr, pnerr, tc, aste, asteperr, astenerr, ecc, ome, rprst, rprstperr, rprstnerr, tce = np.loadtxt('data.dat', usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), unpack = True)
ra, dec = np.loadtxt('data.dat', dtype = str, usecols = (18,19), unpack = True)

f_data = open('data2.dat','w')
f_data.write('#Name\tTeff\tlog(g)\t[M/H]\tVturb\tPeriod\tP+err\tP-err\tTc\ta/R*\ta+err\ta-err\tEccentricity\tOmega\tRp/R*\tr+err\tr-err\tTc-err\n')

print('------------------------------------------------------------------------------------------------------------------')
print('------------------------ Starting First Iteration: To Download data and analyse them------------------------------')
print('------------------------------------------------------------------------------------------------------------------')

for i in range(len(teff)):
	#--------------------------------------
	#--------Downloading Data--------------
	#--------------------------------------
	print('----------------------------------------------------------------')
	print('---------------Working on '+ name[i]+' -------')
	print('----------------------------------------------------------------')
	ra1 = utils.RA_to_deg(ra[i])
	dec1 = utils.DEC_to_deg(dec[i])
	coord = str(ra1) + ' ' + str(dec1)
	obt = obs.query_region(coord,radius=0.001)
	b = np.array([])
	for j in range(len(obt['intentType'])):
		if obt['obs_collection'][j] == 'TESS' and obt['dataproduct_type'][j] == 'timeseries':
			b = np.hstack((b,j))
	if len(b) == 0:
		continue
	obsid1 = np.array([])
	for j1 in range(len(b)):
		ob = int(b[j1])
		obsid1 = np.hstack((obsid1,obt['obsid'][ob]))
	for i11 in range(len(obsid1)):
		print('-------------------------------------------------')
		print('--------Working on '+ name[i] + ' - Sector ' + str(i11) + ' -------')
		print('-------------------------------------------------')
		try:
			tab = obs.download_products(str(obsid1[i11]),extension='fits')
		except:
			continue
		for j2 in range(len(tab['Local Path'])):
			b1 = tab['Local Path'][j2]
			if b1[-7:] == 'lc.fits':
				c1 = j2
		try:
			d = tab['Local Path'][int(c1)]
		except:
			continue
		os.system('mv ' + path1 + d[1:] + ' ' + path1 + '/Light-curve/' + name[i] + '_sector' + str(i11) + '.fits')
		os.system('rm -r mastDownload')
		#--------------------------------------
		#-------Data Downloaded----------------
		#-----And stored at proper folder------
		#--------------------------------------
		#--------------------------------------
		#----Getting Data from fits file-------
		#--------------------------------------
		hdul = fits.open(path1 + '/Light-curve/' + name[i] + '_sector' + str(i11) + '.fits')
		h = hdul[1].data
		#--------------------------------------
		#------Saving Data to Numpy Array------
		#--------------------------------------
		time_np = np.array([])
		flux_np = np.array([])
		fluxerr_np = np.array([])
		for j3 in h['TIME']:
			time_np = np.hstack((time_np,j3))
		for j4 in h['PDCSAP_FLUX']:
			flux_np = np.hstack((flux_np,j4))
		for j5 in h['PDCSAP_FLUX_ERR']:
			fluxerr_np = np.hstack((fluxerr_np,j5))
		time_bjd = time_np + 2457000
		#--------------------------------------
		#--------Removing NaN values-----------
		#--------------------------------------
		nan_val_in_flux = np.isnan(flux_np)
		for j6 in range(len(flux_np)):
			j61 = np.bool(nan_val_in_flux[j6])
			if j61 is True:
				time_bjd[j6] = flux_np[j6]
				fluxerr_np[j6] = flux_np[j6]
			else:
				continue
		time_bjd = time_bjd[np.isfinite(time_bjd)]
		flux_np = flux_np[np.isfinite(flux_np)]
		fluxerr_np = fluxerr_np[np.isfinite(fluxerr_np)]
		#--------------------------------------
		#--------Calculating Relative flux-----
		#--And generating light-curve from it--
		#--------------------------------------
		median_flux = np.median(flux_np)
		rel_flux = flux_np/median_flux
		rel_flux_err = fluxerr_np/median_flux
		fig_lightcurve = plt.figure()
		plt.errorbar(time_bjd, rel_flux, yerr=rel_flux_err, fmt='.')
		plt.title('Light-curve for ' + name[i] + '_sector' + str(i11) + ' system')
		plt.xlabel('Time (in BJD)')
		plt.ylabel('Relative Flux')
		plt.grid()
		plt.savefig('Fig.png')
		plt.close(fig_lightcurve)
		utl.move_file(in_path = path1, fi_name='/Fig.png', out_path = path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/', new_name='Fig.png')
		#--------------------------------------
		#-------Creating Data file and---------
		#-------External Parameter file--------
		#--------------------------------------
		j7 = np.vstack((time_bjd, rel_flux))
		j77 = np.vstack((j7, rel_flux_err))
		data_to_be_dumped = np.transpose(j77)
		np.savetxt(name[i] + '_sector' + str(i11) +'_lc.dat', data_to_be_dumped, newline='\tTESS\n', delimiter='\t')
		utl.move_file(in_path = path1 + '/', fi_name = name[i] + '_sector' + str(i11) + '_lc.dat', out_path = path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/data/', new_name = name[i] + '_sector' + str(i11) + '_lc.dat')
		np.savetxt(name[i] + '_sector' + str(i11) + '_lceparam.dat', time_bjd, newline = '\tTESS\n', delimiter = '\t')
		utl.move_file(in_path = path1 + '/', fi_name = name[i] + '_sector' + str(i11) + '_lceparam.dat', out_path = path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) +'/data/', new_name = name[i] + '_sector' + str(i11) + '_lceparam.dat')
		#--------------------------------------
		#------Creating Prior File-------------
		#--------------------------------------
		phase = (time_bjd - tc[i])/p[i]
		phase1 = int(phase[0])+1
		for j8 in range(len(time_bjd)):
			phase2 = phase[j8]-phase1
			if phase2<0.0001:
				j88 = j8
		t0 = time_bjd[j88]
		f1_p = open(name[i] + '_sector' + str(i11) +'_priors_exm.dat','w')#For Exponential-Matern Kernel
		f2_p = open(name[i] + '_sector' + str(i11) +'_priors_qp.dat','w')#For Quasi-Periodic Kernel
		f3_p = open(name[i] + '_sector' + str(i11) +'_priors_exm1.dat','w')#For Ex-M Kernel with 'a' uniformally distributed
		f4_p = open(name[i] + '_sector' + str(i11) +'_priors_qp1.dat','w')#For QP Kernel with 'a' uniformally distributed
		f1_p.write('#Physical parameters of the transiting system:\n')
		f1_p.write('P_p1\t\t\tNormal\t\t' + str(p[i]) + ',' + str(pperr[i]) + '\n')
		f1_p.write('r1_p1\t\t\tUniform\t\t0.0,1.0\n')
		f1_p.write('r2_p1\t\t\tUniform\t\t0.0,1.0\n')
		f1_p.write('a_p1\t\t\tNormal\t\t' + str(aste[i]) + ',' + str(asteperr[i]) + '\n')
		f1_p.write('t0_p1\t\t\tNormal\t\t' + str(t0) + ',0.1\n')
		f1_p.write('ecc_p1\t\t\tFIXED\t\t' + str(ecc[i]) + '\n')
		f1_p.write('omega_p1\t\tFIXED\t\t' + str(ome[i]) + '\n')
		f1_p.write('#Photometric priors for TESS photometry:\n')
		f1_p.write('q1_TESS\t\t\tUniform\t\t0.0,1.0\n')
		f1_p.write('q2_TESS\t\t\tUniform\t\t0.0,1.0\n')
		f1_p.write('mdilution_TESS\t\tFIXED\t\t1.0\n')
		f1_p.write('mflux_TESS\t\tNormal\t\t0.0,1.0\n')
		f1_p.write('sigma_w_TESS\t\tJeffreys\t0.1,10000\n')
		f1_p.write('GP_sigma_TESS\t\tJeffreys\t0.1,10000\n')
		f1_p.write('GP_rho_TESS\t\tUniform\t\t0.001,1e4\n')
		f1_p.write('GP_timescale_TESS\tUniform\t\t0.001,1e4')
		f1_p.close()
		f2_p.write('#Physical parameters of the transiting system:\n')
		f2_p.write('P_p1\t\t\tNormal\t\t' + str(p[i]) + ',' + str(pperr[i]) + '\n')
		f2_p.write('r1_p1\t\t\tUniform\t\t0.0,1.0\n')
		f2_p.write('r2_p1\t\t\tUniform\t\t0.0,1.0\n')
		f2_p.write('a_p1\t\t\tNormal\t\t' + str(aste[i]) + ',' + str(asteperr[i]) + '\n')
		f2_p.write('t0_p1\t\t\tNormal\t\t' + str(t0) + ',0.1\n')
		f2_p.write('ecc_p1\t\t\tFIXED\t\t' + str(ecc[i]) + '\n')
		f2_p.write('omega_p1\t\tFIXED\t\t' + str(ome[i]) + '\n')
		f2_p.write('#Photometric priors for TESS photometry:\n')
		f2_p.write('q1_TESS\t\t\tUniform\t\t0.0,1.0\n')
		f2_p.write('q2_TESS\t\t\tUniform\t\t0.0,1.0\n')
		f2_p.write('mdilution_TESS\t\tFIXED\t\t1.0\n')
		f2_p.write('mflux_TESS\t\tNormal\t\t0.0,1.0\n')
		f2_p.write('sigma_w_TESS\t\tJeffreys\t0.1,10000\n')
		f2_p.write('GP_B_TESS\t\tJeffreys\t1e-5,1000\n')
		f2_p.write('GP_C_TESS\t\tJeffreys\t1e-5,1000\n')
		f2_p.write('GP_L_TESS\t\tJeffreys\t1e-5,1000\n')
		f2_p.write('GP_Prot_TESS\t\tJeffreys\t0.1,30\n')
		f2_p.close()
		f3_p.write('#Physical parameters of the transiting system:\n')
		f3_p.write('P_p1\t\t\tNormal\t\t' + str(p[i]) + ',' + str(pperr[i]) + '\n')
		f3_p.write('r1_p1\t\t\tUniform\t\t0.0,1.0\n')
		f3_p.write('r2_p1\t\t\tUniform\t\t0.0,1.0\n')
		f3_p.write('a_p1\t\t\tUniform\t\t1,100\n')
		f3_p.write('t0_p1\t\t\tNormal\t\t' + str(t0) + ',0.1\n')
		f3_p.write('ecc_p1\t\t\tFIXED\t\t' + str(ecc[i]) + '\n')
		f3_p.write('omega_p1\t\tFIXED\t\t' + str(ome[i]) + '\n')
		f3_p.write('#Photometric priors for TESS photometry:\n')
		f3_p.write('q1_TESS\t\t\tUniform\t\t0.0,1.0\n')
		f3_p.write('q2_TESS\t\t\tUniform\t\t0.0,1.0\n')
		f3_p.write('mdilution_TESS\t\tFIXED\t\t1.0\n')
		f3_p.write('mflux_TESS\t\tNormal\t\t0.0,1.0\n')
		f3_p.write('sigma_w_TESS\t\tJeffreys\t0.1,10000\n')
		f3_p.write('GP_sigma_TESS\t\tJeffreys\t0.1,10000\n')
		f3_p.write('GP_rho_TESS\t\tUniform\t\t0.001,1e4\n')
		f3_p.write('GP_timescale_TESS\tUniform\t\t0.001,1e4')
		f3_p.close()
		f4_p.write('#Physical parameters of the transiting system:\n')
		f4_p.write('P_p1\t\t\tNormal\t\t' + str(p[i]) + ',' + str(pperr[i]) + '\n')
		f4_p.write('r1_p1\t\t\tUniform\t\t0.0,1.0\n')
		f4_p.write('r2_p1\t\t\tUniform\t\t0.0,1.0\n')
		f4_p.write('a_p1\t\t\tUniform\t\t1,100\n')
		f4_p.write('t0_p1\t\t\tNormal\t\t' + str(t0) + ',0.1\n')
		f4_p.write('ecc_p1\t\t\tFIXED\t\t' + str(ecc[i]) + '\n')
		f4_p.write('omega_p1\t\tFIXED\t\t' + str(ome[i]) + '\n')
		f4_p.write('#Photometric priors for TESS photometry:\n')
		f4_p.write('q1_TESS\t\t\tUniform\t\t0.0,1.0\n')
		f4_p.write('q2_TESS\t\t\tUniform\t\t0.0,1.0\n')
		f4_p.write('mdilution_TESS\t\tFIXED\t\t1.0\n')
		f4_p.write('mflux_TESS\t\tNormal\t\t0.0,1.0\n')
		f4_p.write('sigma_w_TESS\t\tJeffreys\t0.1,10000\n')
		f4_p.write('GP_B_TESS\t\tJeffreys\t1e-5,1000\n')
		f4_p.write('GP_C_TESS\t\tJeffreys\t1e-5,1000\n')
		f4_p.write('GP_L_TESS\t\tJeffreys\t1e-5,1000\n')
		f4_p.write('GP_Prot_TESS\t\tJeffreys\t0.1,30\n')
		f4_p.close()
		utl.move_file(path1 + '/', name[i] + '_sector' + str(i11) + '_priors_exm.dat', path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/priors/', name[i] + '_sector' + str(i11) + '_priors_exm.dat')
		utl.move_file(path1 + '/', name[i] + '_sector' + str(i11) + '_priors_qp.dat', path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) +'/priors/', name[i] + '_sector' + str(i11) + '_priors_qp.dat')
		utl.move_file(path1 + '/', name[i] + '_sector' + str(i11) + '_priors_exm1.dat', path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/priors/', name[i] + '_sector' + str(i11) + '_priors_exm1.dat')
		utl.move_file(path1 + '/', name[i] + '_sector' + str(i11) + '_priors_qp1.dat', path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/priors/', name[i] + '_sector' + str(i11) + '_priors_qp1.dat')
		utl.new_dir(path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/results/qp')
		utl.new_dir(path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/results/exm')
		utl.new_dir(path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) +'/results/qp1')
		utl.new_dir(path1 + '/Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/results/exm1')
		#--------------------------------------
		#-------Created Prior Files------------
		#--------------------------------------
		#--------------------------------------
		#---------Running Juliet---------------
		#--------------------------------------
		#os.system('cd ~')
		#os.system('export LD_LIBRARY_PATH=/home/jayshil/MultiNest/lib/:$LD_LIBRARY_PATH~/.bashrc')
		#os.system('export PATH=$PATH:$HOME/.local/bin/~/.bashrc')
		#os.system('cd ' + path1)
		os.system('python juliet.py -lcfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/data/' + name[i] + '_sector' + str(i11) + '_lc.dat' + ' -lceparamfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/data/' + name[i] + '_sector' + str(i11) + '_lceparam.dat -ldlaw TESS-quadratic -lctimedef TESS-TDB -priorfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/priors/' + name[i] + '_sector' + str(i11) + '_priors_exm.dat -ofolder Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/results/exm -nlive 500')#Ran simulation with Ex-M Kernel
		os.system('python juliet.py -lcfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/data/' + name[i] + '_sector' + str(i11) + '_lc.dat' + ' -lceparamfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/data/' + name[i] + '_sector' + str(i11) + '_lceparam.dat -ldlaw TESS-quadratic -lctimedef TESS-TDB -priorfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/priors/' + name[i] + '_sector' + str(i11) + '_priors_qp.dat -ofolder Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/results/qp -nlive 500')#Ran simulation with QP Kernel
		os.system('python juliet.py -lcfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/data/' + name[i] + '_sector' + str(i11) + '_lc.dat' + ' -lceparamfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/data/' + name[i] + '_sector' + str(i11) + '_lceparam.dat -ldlaw TESS-quadratic -lctimedef TESS-TDB -priorfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/priors/' + name[i] + '_sector' + str(i11) + '_priors_exm1.dat -ofolder Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/results/exm1 -nlive 500')#Ran simulation with Ex-M Kernel and making 'a' uniformally distributed
		os.system('python juliet.py -lcfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/data/' + name[i] + '_sector' + str(i11) + '_lc.dat' + ' -lceparamfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/data/' + name[i] + '_sector' + str(i11) + '_lceparam.dat -ldlaw TESS-quadratic -lctimedef TESS-TDB -priorfile Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/priors/' + name[i] + '_sector' + str(i11) + '_priors_qp1.dat -ofolder Simulations/' + name[i] + '/' + name[i] + '_sector' + str(i11) + '/results/qp1 -nlive 500')#Ran simulation with QP Kernel and making 'a' uniformally distributed
	f_data.write(name[i] + '\t' + str(teff[i]) + '\t' + str(lg[i]) + '\t' + str(mh[i]) + '\t' + str(vturb[i]) + '\t' + str(p[i]) + '\t' + str(pperr[i]) + '\t' + str(pnerr[i]) + '\t' + str(tc[i]) + '\t' + str(aste[i]) + '\t' + str(asteperr[i]) + '\t' + str(astenerr[i]) + '\t' + str(ecc[i]) + '\t' + str(ome[i]) + '\t' + str(rprst[i]) + '\t' + str(rprstperr[i]) + '\t' + str(rprstnerr[i]) + '\t' + str(tce[i]) + '\t' + ra[i] + '\t' + dec[i] + '\n')
	print('******************************************************************************************')
	print('                                                                                          ')
	print('Completed analysis for ' + str(i+1) + ' system(s) / out of ' + str(len(name)) + ' systems\n')
	print('                                                                                          ')
	print('******************************************************************************************')

f_data.close()

print('------------------------------------------------------------------------------------------------------------------')
print('--------------------------------- First Iteration Completed Successfully------------------------------------------')
print('------------------------------------------------------------------------------------------------------------------')

#--------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------Gathering the available data---------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

name2 = np.loadtxt('data2.dat', dtype = str, usecols = 0, unpack = True)
teff2, lg2, mh2, vturb2, p2, pperr2, pnerr2, tc2, aste2, asteperr2, astenerr2, ecc2, ome2, rprst2, rprstperr2, rprstnerr2, tce2 = np.loadtxt('data2.dat', usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), unpack = True)
ra2, dec2 = np.loadtxt('data2.dat', dtype = str, usecols = (18,19), unpack = True)

print('--------------------------------------------------------------------------------------------------------------------------')
print('------------------------ Starting Second Iteration: To Retrieve Values of Various Parameters------------------------------')
print('--------------------------------------------------------------------------------------------------------------------------')

for i in range(len(name2)):
	sec_list = os.listdir(path1 + '/Simulations/' + name2[i])
	sec_list.sort(key = utl.natural_keys)
	ast_tem = []
	per_tem = []
	rp_tem = []
	tc_tem = []
	u1_tem = []
	u2_tem = []
	for j in range(len(sec_list)):
		post_exm1 = pc.load(open(path1 + '/Simulations/' + name2[i] + '/' + sec_list[j] + '/results/exm1/posteriors.pkl', 'rb'), encoding = 'latin1')
		post_qp1 = pc.load(open(path1 + '/Simulations/' + name2[i] + '/' + sec_list[j] + '/results/qp1/posteriors.pkl', 'rb'), encoding = 'latin1')
		lnze1 = post_exm1['lnZ']
		lnzqp1 = post_qp1['lnZ']
		mo1 = np.maximum(lnze1,lnzqp1)
		if lnze1 == mo1:
			model = 'exm1'
		elif lnzqp1 == mo1:
			model = 'qp1'
		post_best = pc.load(open(path1 + '/Simulations/' + name2[i] + '/' + sec_list[j] + '/results/' + model + '/posteriors.pkl', 'rb'), encoding = 'latin1')
		pbe = post_best['posterior_samples']
		ast_tem.append(pbe['a_p1'])
		per_tem.append(pbe['P_p1'])
		tc_tem.append(pbe['t0_p1'])
		bex, p11 = utils.convert_bp(pbe['r1_p1'], pbe['r2_p1'], post_best['pl'], post_best['pu'])
		rp_tem.append(p11)
		u11, u22 = utils.reverse_ld_coeffs('quadratic', pbe['q1_TESS'], pbe['q2_TESS'])
		u1_tem.append(u11)
		u2_tem.append(u22)
	length = np.array([])
	for k in range(len(ast_tem)):
		length = np.hstack((length, len(ast_tem[k])))
	mini1 = np.min(length)
	mini = int(mini1)
	#------------------------------
	#--------- For a/R* -----------
	#------------------------------
	ast_vec = np.random.choice(ast_tem[0], size = mini, replace = False)
	for j1 in range(len(ast_tem)-1):
		xx = np.random.choice(ast_tem[j1+1], size = mini, replace = False)
		ast_vec = np.vstack((ast_vec, xx))
	ast_vec_tran = np.transpose(ast_vec)
	ast_fin = np.array([])
	for j2 in range(mini):
		yy = np.average(ast_vec_tran[j2])
		ast_fin = np.hstack((ast_fin, yy))
	f7.write(name2[i] + '\t' + str(np.median(ast_fin)) + '\t' + str(np.std(ast_fin)) + '\t' + str(np.std(ast_fin)) + '\n')
	#------------------------------
	#----------- For P ------------
	#------------------------------
	per_vec = np.random.choice(per_tem[0], size = mini, replace = False)
	for j1 in range(len(per_tem)-1):
		xx = np.random.choice(per_tem[j1+1], size = mini, replace = False)
		per_vec = np.vstack((per_vec, xx))
	per_vec_tran = np.transpose(per_vec)
	per_fin = np.array([])
	for j2 in range(mini):
		yy = np.average(per_vec_tran[j2])
		per_fin = np.hstack((per_fin, yy))
	f8.write(name2[i] + '\t' + str(np.median(per_fin)) + '\t' + str(np.std(per_fin)) + '\t' + str(np.std(per_fin)) + '\n')
	#------------------------------
	#--------- For Rp--------------
	#------------------------------
	rp_vec = np.random.choice(rp_tem[0], size = mini, replace = False)
	for j1 in range(len(rp_tem)-1):
		xx = np.random.choice(rp_tem[j1+1], size = mini, replace = False)
		rp_vec = np.vstack((rp_vec, xx))
	rp_vec_tran = np.transpose(rp_vec)
	rp_fin = np.array([])
	for j2 in range(mini):
		yy = np.average(rp_vec_tran[j2])
		rp_fin = np.hstack((rp_fin, yy))
	f9.write(name2[i] + '\t' + str(np.median(rp_fin)) + '\t' + str(np.std(rp_fin)) + '\t' + str(np.std(rp_fin)) + '\n')
	#------------------------------
	#--------- For Tc--------------
	#------------------------------
	tc_bjd1 = np.loadtxt(path1 + '/Simulations/' + name2[i] + '/' + sec_list[0] + '/data/' + sec_list[0] + '_lc.dat', usecols = 0, unpack = True)
	phase_to = (tc_bjd1 - tc2[i])/p2[i]
	phase1_to = int(phase_to[0])+1
	tc1 = np.random.normal(tc2[i], tce2[i], 10000)
	p11 = np.random.normal(p2[i], pperr2[i], 10000)
	t11 = tc2[i] + (phase1_to * p11)
	t0_m = np.median(t11)
	t0_s = np.std(t11)
	f1_to.write(name2[i] + '\t' + str(t0_m) + '\t' + str(t0_s) + '\n')
	tc_fin = np.random.choice(tc_tem[0], size = mini, replace = False)
	f10.write(name2[i] + '\t' + str(np.median(tc_fin)) + '\t' + str(np.std(tc_fin)) + '\t' + str(np.std(tc_fin)) + '\n')
	#------------------------------
	#------ For u1 and u2----------
	#------------------------------
	u1_vec = np.random.choice(u1_tem[0], size = mini, replace = False)
	for j1 in range(len(u1_tem)-1):
		xx = np.random.choice(u1_tem[j1+1], size = mini, replace = False)
		u1_vec = np.vstack((u1_vec, xx))
	u1_vec_tran = np.transpose(u1_vec)
	u1_fin = np.array([])
	for j2 in range(mini):
		yy = np.average(u1_vec_tran[j2])
		u1_fin = np.hstack((u1_fin, yy))
	u2_vec = np.random.choice(u2_tem[0], size = mini, replace = False)
	for j1 in range(len(u2_tem)-1):
		xx = np.random.choice(u2_tem[j1+1], size = mini, replace = False)
		u2_vec = np.vstack((u2_vec, xx))
	u2_vec_tran = np.transpose(u2_vec)
	u2_fin = np.array([])
	for j2 in range(mini):
		yy = np.average(u2_vec_tran[j2])
		u2_fin = np.hstack((u2_fin, yy))
	f19.write(name2[i] + '\t' + str(np.median(u1_fin)) + '\t' + str(np.std(u1_fin)) + '\t' + str(np.std(u1_fin)) + '\t' + str(np.median(u2_fin)) + '\t' + str(np.std(u2_fin)) + '\t' + str(np.std(u2_fin)) + '\n')
	print('*****************************************************************************************')
	print('                                                                                         ')
	print('Retrieved values for ' + str(i+1) + ' system(s) / out of ' + str(len(name2)) + ' systems\n')
	print('                                                                                         ')
	print('*****************************************************************************************')

f7.close()
f8.close()
f9.close()
f10.close()
f19.close()
f1_to.close()

print('------------------------------------------------------------------------------------------------------------------')
print('-------------------------------- Second Iteration Completed Successfully------------------------------------------')
print('------------------------------------------------------------------------------------------------------------------')

print('--------------------------------------------------------------------------------------------------------------------------')
print('----------------------------- Starting Third Iteration: To Calculate Theoretical LDCs-------------------------------------')
print('--------------------------------------------------------------------------------------------------------------------------')

for i in range(len(name2)):
	#--------------------------------------
	#------Making input files for----------
	#-------Limb Darkening Code------------
	#--------------------------------------
	fout = open(name2[i] + '.dat','w')
	fout.write('#Name\tTeff\tLog(g)\t[M/H]\tVturb\tRF\t\t\tFT\tmin_w\tmax_w\n')
	fout.write(name2[i] + '\t' + str(teff2[i]) + '\t' + str(lg2[i]) + '\t' + str(mh2[i]) + '\t' + str(vturb2[i]) + '\ttess_res_fun.txt\tA100,P100\t-1\t-1')
	fout.close()
	utl.move_file(path1 + '/', name2[i] + '.dat', path1 + '/limb-darkening-master/input_files/', name2[i] + '.dat')
	#--------------------------------------
	#----Starting Limb-Darkening Code------
	#--------------------------------------
	os.system('python2.7 ' + path1 + '/limb-darkening-master/get_lds.py -ifile ' + path1 +'/limb-darkening-master/input_files/' + name2[i] + '.dat -ofile ' + name2[i] + '_LDC.dat')
	#--------------------------------------
	#---------------ATLAS------------------
	#--------------------------------------
	f5 = open(path1 + '/limb-darkening-master/results/' + name2[i] + '_LDC.dat', 'r')
	con = f5.readlines()
	line = con[24]
	ind = np.array([])
	for i1 in range(len(line)):
		if line[i1] == ',':
			ind = np.hstack((ind,i1))
	u1 = line[int(ind[3]-11):int(ind[3])]
	u2 = line[int(ind[4]-11):int(ind[4])]
	u3 = line[int(ind[5]-11):int(ind[5])]
	u4 = line[int(ind[5]+2):-1]
	f33.write(name2[i] + '\t\t' + u1 + '\t\t' + u2 + '\t\t' + u3 + '\t\t' + u4 + '\n')
	#--------------------------------------
	#-------------PHOENIX------------------
	#--------------------------------------
	line1 = con[36]
	ind1 = np.array([])
	for i2 in range(len(line1)):
		if line[i2] == ',':
			ind1 = np.hstack((ind1,i2))
	u1_p = line1[int(ind1[3]-11):int(ind1[3])]
	u2_p = line1[int(ind1[4]-11):int(ind1[4])]
	u3_p = line1[int(ind1[5]-11):int(ind1[5])]
	u4_p = line1[int(ind1[5]+2):-1]
	f11.write(name2[i] + '\t\t' + u1_p + '\t\t' + u2_p + '\t\t' + u3_p + '\t\t' + u4_p + '\n')
	#--------------------------------------
	#-----Calculating LDCs from PHOENIX----
	#--------from Claret(2017)-------------
	#--------------------------------------
	lg1_p, T1_p, q1_p, q2_p, q3_p, q4_p = np.loadtxt(path1 + '/Phoenix/pho_ldc_claret.dat', usecols=(0,1,4,5,6,7), unpack=True)
	pts_p = np.vstack((lg1_p, T1_p))
	pt_p = np.transpose(pts_p)
	c1_p = inp.griddata(pt_p, q1_p, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	c2_p = inp.griddata(pt_p, q2_p, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	c3_p = inp.griddata(pt_p, q3_p, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	c4_p = inp.griddata(pt_p, q4_p, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	f1.write(name2[i] + '\t' + str(c1_p) + '\t' + str(c2_p) + '\t' + str(c2_p) + '\t' + str(c2_p) + '\n')
	u1_p = ((12*c1_p)/35) + c2_p + ((164*c3_p)/105) + (2*c4_p)
	u2_p = ((10*c1_p)/31) - ((34*c3_p)/63) - c4_p
	f2.write(name2[i] + '\t' + str(u1_p) + '\t' + str(u2_p) + '\n')
	#--------------------------------------
	#-----Calculating LDCs from ATLAS------
	#--------from Claret(2017)-------------
	#--------------------------------------
	lg1_a, T1_a, met_a, q1_a, q2_a, q3_a, q4_a = np.loadtxt(path1 + '/Atlas/atlas_ldc_claret.dat', usecols=(0,1,2,4,5,6,7), unpack=True)
	pts_1 = np.vstack((lg1_a, T1_a))
	pts_a = np.vstack((pts_1,met_a))
	pt_a = np.transpose(pts_a)
	c1_a = inp.griddata(pt_a, q1_a, (lg2[i],teff2[i],mh2[i]), fill_value = 0, method = 'linear')
	c2_a = inp.griddata(pt_a, q2_a, (lg2[i],teff2[i],mh2[i]), fill_value = 0, method = 'linear')
	c3_a = inp.griddata(pt_a, q3_a, (lg2[i],teff2[i],mh2[i]), fill_value = 0, method = 'linear')
	c4_a = inp.griddata(pt_a, q4_a, (lg2[i],teff2[i],mh2[i]), fill_value = 0, method = 'linear')
	f3.write(name2[i] + '\t' + str(c1_a) + '\t' + str(c2_a) + '\t' + str(c2_a) + '\t' + str(c2_a) + '\n')
	u1_a = ((12*c1_a)/35) + c2_a + ((164*c3_a)/105) + (2*c4_a)
	u2_a = ((10*c1_a)/31) - ((34*c3_a)/63) - c4_a
	f4.write(name2[i] + '\t' + str(u1_a) + '\t' + str(u2_a) + '\n')
	print('****************************************************************************************')
	print('                                                                                        ')
	print('Calculated LDCs for ' + str(i+1) + ' system(s) / out of ' + str(len(name2)) + ' systems\n')
	print('                                                                                        ')
	print('****************************************************************************************')


f1.close()
f2.close()
f3.close()
f33.close()
f4.close()
f11.close()

#--------------------------------------------------------------------------------------------------
#----------------Calculating Limiting LDCs from ATLAS (Code)---------------------------------------
#--------------------------------------------------------------------------------------------------

name1 = np.loadtxt(path1 + '/Atlas/code_us_nl_ata.dat',dtype=str,usecols=0,unpack=True)
c1_code_a, c2_code_a, c3_code_a, c4_code_a = np.loadtxt(path1 + '/Atlas/code_us_nl_ata.dat', usecols=(1,2,3,4), unpack=True)

u1_code_a = ((12*c1_code_a)/35) + c2_code_a + ((164*c3_code_a)/105) + (2*c4_code_a)
u2_code_a = ((10*c1_code_a)/31) - ((34*c3_code_a)/63) - c4_code_a

for i in range(len(name1)):
	f44.write(name1[i] + '\t\t' + str(u1_code_a[i]) + '\t\t' + str(u2_code_a[i]) + '\n')

f44.close()

#--------------------------------------------------------------------------------------------------
#----------------Calculating Limiting LDCs from PHOENIX (Code)-------------------------------------
#--------------------------------------------------------------------------------------------------

name11 = np.loadtxt(path1 + '/Phoenix/code_us_nl_pho.dat',dtype=str,usecols=0,unpack=True)
c1_code_p, c2_code_p, c3_code_p, c4_code_p = np.loadtxt(path1 + '/Phoenix/code_us_nl_pho.dat', usecols=(1,2,3,4), unpack=True)

u1_code_p = ((12*c1_code_p)/35) + c2_code_p + ((164*c3_code_p)/105) + (2*c4_code_p)
u2_code_p = ((10*c1_code_p)/31) - ((34*c3_code_p)/63) - c4_code_p

for i in range(len(name11)):
	f22.write(name11[i] + '\t\t' + str(u1_code_p[i]) + '\t\t' + str(u2_code_p[i]) + '\n')

f22.close()

print('--------------------------------------------------------------------------------------------------------------------------')
print('-----------------------------------------Now its time to plot some amazing results----------------------------------------')
print('--------------------------------------------------------------------------------------------------------------------------')

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


#-------------------------------------------------------------------------
print('---------------------Plotting Stellar Parameters-------------------------')
#-------------------------------------------------------------------------

fig_gr = plt.figure(figsize=(27,15.88))
x_gr = np.arange(0, len(name2), 1)
y_gr = np.arange(3.5,5.2,0.1)
plt.errorbar(name2, lg2, fmt='o', elinewidth=1, alpha=0.5, color='darkblue', zorder=5, markersize=16.25)
plt.xticks(x_gr, rotation = 90, fontsize = 20)
plt.yticks(y_gr, fontsize=25)
plt.ylabel('Surface gravity (log(g)) of stellar host', fontsize = 27)
plt.grid()
plt.savefig(path1 + '/Results/stellar_prop/lg.png')
plt.close(fig_gr)

fig_teff = plt.figure(figsize=(27,15.88))
x_teff = np.arange(0, len(name2), 1)
y_teff = np.arange(3000,7500,500)
plt.errorbar(name2, teff2, fmt='o', elinewidth=1, alpha=0.5, color='darkblue', zorder=5, markersize=16.25)
plt.xticks(x_teff, rotation = 90, fontsize = 20)
plt.yticks(y_teff, fontsize = 25)
plt.ylabel('Effective temperature of stellar host', fontsize = 27)
plt.grid()
plt.savefig(path1 + '/Results/stellar_prop/teff.png')
plt.close(fig_teff)

#-------------------------------------------------------------------------
print('----------------Plots to check how good was the fit----------------------')
#-------------------------------------------------------------------------

amax_t = np.max(aste2)
pmax_t = np.max(p2)
rmax_t = np.max(rprst2)
amin_t = np.min(aste2)
pmin_t = np.min(p2)
rmin_t = np.min(rprst2)

#------------------------------------------------
print('-----------------For a/R*-----------------------')
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
	at1 = np.random.normal(aste2[i],asteperr2[i],10000)
	ac1 = np.random.normal(a_j[i], a_jp[i], 10000)
	diff1 = at1 - ac1
	am1 = np.median(diff1)
	ae1 = np.std(diff1)
	diff_a = np.hstack((diff_a,am1))
	diff_ae = np.hstack((diff_ae,ae1))

fig_a = plt.figure(figsize = (8,10))
gs_a = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_a = plt.subplot(gs_a[0])

ax_a.errorbar(aste2, a_j, xerr = [astenerr2, asteperr2], yerr = [a_jn, a_jp], fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xla_j, xua_j])
plt.ylim([xla_j, xua_j])
ax_a.plot(x1a, y1a, 'k--')
ax_a.grid()
plt.ylabel(r'$a/R_*$ (Observed)')
#plt.title('Comparison between literature values and calculated values of a/R*')

ax1_a = plt.subplot(gs_a[1], sharex = ax_a)

ax1_a.errorbar(aste2, diff_a, xerr = [astenerr2, asteperr2], yerr = diff_ae, fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

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
		fa.write(name[i] + '\n')

fa.close()
plt.close(fig_a)

#------------------------------------------------
print('-----------------For Period---------------------')
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
	pt1 = np.random.normal(p2[i],pperr2[i],10000)
	pc1 = np.random.normal(p_j[i], p_jp[i], 10000)
	diff1 = (pt1 - pc1)*86400
	pm1 = np.median(diff1)
	pe1 = np.std(diff1)
	diff_p = np.hstack((diff_p,pm1))
	diff_pe = np.hstack((diff_pe,pe1))

fig_p = plt.figure(figsize = (8,10))
gs_p = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_p = plt.subplot(gs_p[0])

ax_p.errorbar(p2, p_j, xerr = [pnerr2, pperr2], yerr = [p_jn, p_jp], fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xlp_j, xup_j])
plt.ylim([xlp_j, xup_j])
ax_p.plot(x1p, y1p, 'k--')
ax_p.grid()

plt.ylabel('Period (Observed - in days)')
#plt.title('Comparison between literature values and calculated values of period')

ax1_p = plt.subplot(gs_p[1], sharex = ax_p)

ax1_p.errorbar(p2, diff_p, xerr = [pnerr2, pperr2], yerr = diff_pe, fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

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
		fp.write(name[i] + '\n')

fp.close()
plt.close(fig_p)

#------------------------------------------------
print('-----------------For Rp/R*----------------------')
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
	rt1 = np.random.normal(rprst2[i],rprstperr2[i],10000)
	rc1 = np.random.normal(r_j[i], r_jp[i], 10000)
	diff1 = rt1 - rc1
	rm1 = np.median(diff1)
	re1 = np.std(diff1)
	diff_r = np.hstack((diff_r,rm1))
	diff_re = np.hstack((diff_re,re1))

fig_r = plt.figure(figsize = (8,10))
gs_r = gd.GridSpec(2, 1, height_ratios = [4,1])

ax_r = plt.subplot(gs_r[0])

ax_r.errorbar(rprst2, r_j, xerr = [rprstnerr2, rprstperr2], yerr = [r_jn, r_jp], fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xlr_j, xur_j])
plt.ylim([xlr_j, xur_j])
ax_r.plot(x1r, y1r, 'k--')
ax_r.grid()

plt.ylabel(r'$R_p/R_*$ (Observed)')
#plt.title('Comparison between literature values and calculated values of Rp/R*')

ax1_r = plt.subplot(gs_r[1], sharex = ax_r)

ax1_r.errorbar(rprst2, diff_r, xerr = [rprstnerr2, rprstperr2], yerr = diff_re, fmt = '.', elinewidth=1, alpha=0.5, color='black', zorder=5)

plt.xlim([xlr_j, xur_j])
plt.ylim([-0.20,0.10])
ax1_r.plot(x1r, y11r, 'k--')
ax1_r.grid()
plt.ylabel('Residuals')
plt.xlabel('Values of Rp/R* taken from literature')

plt.subplots_adjust(hspace = 0.2)
plt.savefig(path1 + '/Results/comp_a_r_p/r.pdf')

diff_r1 = np.abs(diff_r)
errr1 = 3*(rprstperr2 + r_jp)

fr = open(path1 + '/Results/comp_a_r_p/off_r.dat', 'w')
fr.write('These are the systems which have residulas larger than 3-sigma\n')

for i in range(len(diff_r)):
	if diff_r1[i] > errr1[i]:
		fr.write(name2[i]+ '\n')

fr.close()
plt.close(fig_r)

#------------------------------------------------
print('----------------For Tc--------------------------')
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
plt.ylabel('Residuals (s)')
plt.xlabel('Values of Tc taken from literature')

plt.subplots_adjust(hspace = 0.2)
plt.savefig(path1 + '/Results/comp_a_r_p/tc.pdf')
plt.close(fig_t)

diff_t1 = np.abs(diff_t)
errt1 = 3*(tc_te + tc_jp)*1440

ft = open(path1 + '/Results/comp_a_r_p/off_t.dat', 'w')
ft.write('These are the systems which have residulas larger than 3-sigma\n')

for i in range(len(diff_t)):
	if diff_t1[i] > errt1[i]:
		ft.write(name[i] + '\n')

ft.close()

#-------------------------------------------------------------------------
print('-------------------------Plots to compare LDCs---------------------------')
#-------------------------------------------------------------------------

u1_j, u1_jp, u1_jn, u2_j, u2_jp, u2_jn = np.loadtxt(path1 + '/Results/cal_us_and_evidance/cal_u1_u2.dat', usecols = (1,2,3,4,5,6), unpack = True)

u1max_j = np.max(u1_j)
u1min_j = np.min(u1_j)
u2max_j = np.max(u2_j)
u2min_j = np.min(u2_j)

#----------------------------------------
#--How good are the Claret(2017) LDCs?---
#----------------------------------------

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
#--------u1-------------
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

plt.xlim([xlo, xup])
plt.ylim([xlo, xup])

plt.legend(loc='best')
plt.ylabel(r'$u_1$ (Theoretical)')
plt.title('Values from Claret(2017)')

ax1_u1_c_p = plt.subplot(gs_u1_c_p[1], sharex = ax_u1_c_p)

ax1_u1_c_p.errorbar(u1_j, diff_u1_c_p, xerr = [u1_jn, u1_jp], yerr = diff_u1_c_pe, fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5)
ax1_u1_c_p.errorbar(u1_j, diff_u1_c_a, xerr = [u1_jn, u1_jp], yerr = diff_u1_c_ae, fmt = '.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5)

ax1_u1_c_p.plot(x1u1_c_p, y11u1_c_p, 'k--')
ax1_u1_c_p.grid()
plt.ylabel('Residuals')
plt.xlabel(r'$u_1$ (Observed)')

plt.subplots_adjust(hspace = 0.2)
plt.savefig(path1 + '/Results/cal_us_and_evidance/u1_cla.pdf')
plt.close(fig_u1_c_p)

#-----------------------
#--------u2-------------
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

plt.xlim([xlo, xup])
plt.ylim([xlo, xup])

plt.legend(loc='best')
plt.ylabel(r'$u_2$ (Theoretical)')
plt.title('Values from Claret(2017)')

ax1_u2_c_p = plt.subplot(gs_u2_c_p[1], sharex = ax_u2_c_p)

ax1_u2_c_p.errorbar(u2_j, diff_u2_c_p, xerr = [u2_jn, u2_jp], yerr = diff_u2_c_pe, fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5)
ax1_u2_c_p.errorbar(u2_j, diff_u2_c_a, xerr = [u2_jn, u2_jp], yerr = diff_u2_c_ae, fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5)

ax1_u2_c_p.plot(x1u2_c_p, y11u2_c_p, 'k--')
ax1_u2_c_p.grid()
plt.ylabel('Residuals')
plt.xlabel(r'$u_2$ (Observed)')

plt.subplots_adjust(hspace = 0.2)
plt.savefig(path1 + '/Results/cal_us_and_evidance/u2_cla.pdf')
plt.close(fig_u2_c_p)

#------------------------------------------------
#--------How good are the Code LDCs?-------------
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
#--------u1-------------
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

plt.xlim([xlo, xup])
plt.ylim([xlo, xup])

plt.legend(loc='best')
plt.ylabel(r'$u_1$ (Theoretical)')
plt.title('Values from Espinoza \& Jordan(2015)')

ax1_u1_co_p = plt.subplot(gs_u1_co_p[1], sharex = ax_u1_co_p)

ax1_u1_co_p.errorbar(u1_j, diff_u1_co_p, xerr = [u1_jn, u1_jp], yerr = diff_u1_co_pe, fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5)
ax1_u1_co_p.errorbar(u1_j, diff_u1_co_a, xerr = [u1_jn, u1_jp], yerr = diff_u1_co_ae, fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5)

ax1_u1_co_p.plot(x1u1_co_p, y11u1_co_p, 'k--')
ax1_u1_co_p.grid()
plt.ylabel('Residuals')
plt.xlabel(r'$u_1$ (Observed)')

plt.subplots_adjust(hspace = 0.2)
plt.savefig(path1 + '/Results/cal_us_and_evidance/u1_code.pdf')
plt.close(fig_u1_co_p)

#-----------------------
#--------u2-------------
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

plt.xlim([xlo, xup])
plt.ylim([xlo, xup])

plt.legend(loc='best')
plt.ylabel(r'$u_2$ (Theoretical)')
plt.title('Values from Espinoza \& Jordan(2015)')

ax1_u2_co_p = plt.subplot(gs_u2_co_p[1], sharex = ax_u2_co_p)

ax1_u2_co_p.errorbar(u2_j, diff_u2_co_p, xerr = [u2_jn, u2_jp], yerr = diff_u2_co_pe, fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5)
ax1_u2_co_p.errorbar(u2_j, diff_u2_co_a, xerr = [u2_jn, u2_jp], yerr = diff_u2_co_ae, fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5)

ax1_u2_co_p.plot(x1u2_co_p, y11u2_co_p, 'k--')
ax1_u2_co_p.grid()
plt.ylabel('Residuals')
plt.xlabel(r'$u_2$ (Observed)')

plt.subplots_adjust(hspace = 0.2)
plt.savefig(path1 + '/Results/cal_us_and_evidance/u2_code.pdf')
plt.close(fig_u2_co_p)


#--------------------------------------------------------------------------------------------------
print('----------------------------Making Table of mean offset in LDCs-----------------------------------')
#--------------------------------------------------------------------------------------------------

f101 = open(path1 + '/Results/mean_off.dat', 'w')

#-----------diff_u1_c_p
m1 = np.array([])
diff = np.random.normal(0,0,10000)
mean1 = np.array([])
for i in range(len(diff_u1_c_p)):
	diff1 = np.random.normal(diff_u1_c_p[i], diff_u1_c_pe[i], 10000)
	diff = np.vstack((diff, diff1))
for i in range(10000):
	for j in range(len(diff_u1_c_p)):
		m1 = np.hstack((m1,diff[j+1][i]))
		m2 = np.mean(m1)
		mean1 = np.hstack((mean1,m2))
med1 = np.median(mean1)
std1 = np.std(mean1)

#-----------diff_u2_c_p
m1 = np.array([])
diff = np.random.normal(0,0,10000)
mean2 = np.array([])
for i in range(len(diff_u2_c_p)):
	diff1 = np.random.normal(diff_u2_c_p[i], diff_u2_c_pe[i], 10000)
	diff = np.vstack((diff, diff1))
for i in range(10000):
	for j in range(len(diff_u2_c_p)):
		m1 = np.hstack((m1,diff[j+1][i]))
		m2 = np.mean(m1)
		mean2 = np.hstack((mean2,m2))
med2 = np.median(mean2)
std2 = np.std(mean2)

#-----------diff_u1_c_a
m1 = np.array([])
diff = np.random.normal(0,0,10000)
mean3 = np.array([])
for i in range(len(diff_u1_c_a)):
	diff1 = np.random.normal(diff_u1_c_a[i], diff_u1_c_ae[i], 10000)
	diff = np.vstack((diff, diff1))
for i in range(10000):
	for j in range(len(diff_u1_c_a)):
		m1 = np.hstack((m1,diff[j+1][i]))
		m2 = np.mean(m1)
		mean3 = np.hstack((mean3,m2))
med3 = np.median(mean3)
std3 = np.std(mean3)

#-----------diff_u2_c_a
m1 = np.array([])
diff = np.random.normal(0,0,10000)
mean4 = np.array([])
for i in range(len(diff_u2_c_a)):
	diff1 = np.random.normal(diff_u2_c_a[i], diff_u2_c_ae[i], 10000)
	diff = np.vstack((diff, diff1))
for i in range(10000):
	for j in range(len(diff_u2_c_a)):
		m1 = np.hstack((m1,diff[j+1][i]))
		m2 = np.mean(m1)
		mean4 = np.hstack((mean4,m2))
med4 = np.median(mean4)
std4 = np.std(mean4)

###############################

#-----------diff_u1_co_p
m1 = np.array([])
diff = np.random.normal(0,0,10000)
mean5 = np.array([])
for i in range(len(diff_u1_co_p)):
	diff1 = np.random.normal(diff_u1_co_p[i], diff_u1_co_pe[i], 10000)
	diff = np.vstack((diff, diff1))
for i in range(10000):
	for j in range(len(diff_u1_co_p)):
		m1 = np.hstack((m1,diff[j+1][i]))
		m2 = np.mean(m1)
		mean5 = np.hstack((mean5,m2))
med5 = np.median(mean5)
std5 = np.std(mean5)

#-----------diff_u2_co_p
m1 = np.array([])
diff = np.random.normal(0,0,10000)
mean6 = np.array([])
for i in range(len(diff_u2_co_p)):
	diff1 = np.random.normal(diff_u2_co_p[i], diff_u2_co_pe[i], 10000)
	diff = np.vstack((diff, diff1))
for i in range(10000):
	for j in range(len(diff_u2_c_p)):
		m1 = np.hstack((m1,diff[j+1][i]))
		m2 = np.mean(m1)
		mean6 = np.hstack((mean6,m2))
med6 = np.median(mean6)
std6 = np.std(mean6)

#-----------diff_u1_co_a
m1 = np.array([])
diff = np.random.normal(0,0,10000)
mean7 = np.array([])
for i in range(len(diff_u1_co_a)):
	diff1 = np.random.normal(diff_u1_co_a[i], diff_u1_co_ae[i], 10000)
	diff = np.vstack((diff, diff1))
for i in range(10000):
	for j in range(len(diff_u1_co_a)):
		m1 = np.hstack((m1,diff[j+1][i]))
		m2 = np.mean(m1)
		mean7 = np.hstack((mean7,m2))
med7 = np.median(mean7)
std7 = np.std(mean7)

#-----------diff_u2_co_a
m1 = np.array([])
diff = np.random.normal(0,0,10000)
mean8 = np.array([])
for i in range(len(diff_u2_co_a)):
	diff1 = np.random.normal(diff_u2_co_a[i], diff_u2_co_ae[i], 10000)
	diff = np.vstack((diff, diff1))
for i in range(10000):
	for j in range(len(diff_u2_co_a)):
		m1 = np.hstack((m1,diff[j+1][i]))
		m2 = np.mean(m1)
		mean8 = np.hstack((mean8,m2))
med8 = np.median(mean8)
std8 = np.std(mean8)

f101.write('\t\t\t\tu1\t\t\t\tu2\n')
f101.write('Claret(2017), PHOENIX\t\t' + str(med1) + ' +/- ' + str(std1) + '\t\t' + str(med2) + ' +/- ' + str(std2) + '\n')
f101.write('Claret(2017), ATLAS\t\t' + str(med3) + ' +/- ' + str(std3) + '\t\t' + str(med4) + ' +/- ' + str(std4) + '\n')
f101.write('EJ(2015), PHOENIX\t\t' + str(med5) + ' +/- ' + str(std5) + '\t\t' + str(med6) + ' +/- ' + str(std6) + '\n')
f101.write('EJ(2015), ATLAS\t\t\t' + str(med7) + ' +/- ' + str(std7) + '\t\t' + str(med8) + ' +/- ' + str(std8) + '\n')
f101.close()

#--------------------------------------------------------------------------------------------------
print('-------------------------------Plots of LDCs with temperature variation---------------------------')
#--------------------------------------------------------------------------------------------------

tmin = np.min(teff2)
tmax = np.max(teff2)

x = np.linspace(tmin, tmax, 100)
y = np.zeros(len(x))

def line(x,m,c):
	function = m*x + c
	return function

#------------------------------------------------
#---------u1 - Claret(2017) ---------------------
#------------------------------------------------

fig1 = plt.figure(figsize=(8,6))

plt.fill_between(teff2, med1 + std1, med1 - std1, color = 'red', alpha = 0.5)
plt.errorbar(teff2, diff_u1_c_p, yerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
plt.fill_between(teff2, med3 + std3, med3 - std3, color = 'blue', alpha = 0.5)
plt.errorbar(teff2, diff_u1_c_a, yerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

popt, pcov = cft(line, teff2, diff_u1_c_p)
popt1, pcov1 = cft(line, teff2, diff_u1_c_a)

plt.plot(teff2, line(teff2, *popt), color='orangered', linestyle=':')
plt.plot(teff2, line(teff2, *popt1), color='cornflowerblue', linestyle=':')

plt.plot(x, y, 'k--')
plt.grid()
plt.legend(loc='best')
plt.ylabel(r'Residuals in $u_1$ Claret(2017) LDCs')
plt.xlabel('Effective Temperature')
plt.savefig(path1 + '/Results/variation_with_temp/u1_cla_te.pdf')
plt.close(fig1)

#------------------------------------------------
#---------u2 - Claret(2017) ---------------------
#------------------------------------------------

fig2 = plt.figure(figsize=(8,6))

plt.fill_between(teff2, med2 + std2, med2 - std2, color = 'red', alpha = 0.5)
plt.errorbar(teff2, diff_u2_c_p, yerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
plt.fill_between(teff2, med4 + std4, med4 - std4, color = 'blue', alpha = 0.5)
plt.errorbar(teff2, diff_u2_c_a, yerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

popt2, pcov2 = cft(line, teff2, diff_u2_c_p)
popt3, pcov3 = cft(line, teff2, diff_u2_c_a)

plt.plot(teff2, line(teff2, *popt2), color='orangered', linestyle=':')
plt.plot(teff2, line(teff2, *popt3), color='cornflowerblue', linestyle=':')

plt.plot(x, y, 'k--')
plt.grid()
plt.legend(loc='best')
plt.ylabel(r'Residuals in $u_2$ Claret(2017) LDCs')
plt.xlabel('Effective Temperature')
plt.savefig(path1 + '/Results/variation_with_temp/u2_cla_te.pdf')
plt.close(fig2)

#------------------------------------------------
#-----------u1 - Code ---------------------------
#------------------------------------------------

fig5 = plt.figure()

plt.fill_between(teff2, med5 + std5, med5 - std5, color = 'red', alpha = 0.5)
plt.errorbar(teff2, diff_u1_co_p, yerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
plt.fill_between(teff2, med7 + std7, med7 - std7, color = 'blue', alpha = 0.5)
plt.errorbar(teff2, diff_u1_co_a, yerr = [u1_jn, u1_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

popt4, pcov4 = cft(line, teff2, diff_u1_co_p)
popt5, pcov5 = cft(line, teff2, diff_u1_co_a)

plt.plot(teff2, line(teff2, *popt4), color='orangered', linestyle=':')
plt.plot(teff2, line(teff2, *popt5), color='cornflowerblue', linestyle=':')

plt.plot(x, y, 'k--')
plt.grid()
plt.legend(loc='best')
plt.ylabel(r'Residuals in $u_1$ - EJ(2015) LDCs')
plt.xlabel('Effective Temperature')
plt.savefig(path1 + '/Results/variation_with_temp/u1_code_te.pdf')
plt.close(fig5)

#------------------------------------------------
#------------u2 - Code - PHOENIX-----------------
#------------------------------------------------

fig6 = plt.figure()

plt.fill_between(teff2, med6 + std6, med6 - std6, color = 'red', alpha = 0.5)
plt.errorbar(teff2, diff_u2_co_p, yerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='orangered', zorder=5, label = 'PHOENIX LDCs')
plt.fill_between(teff2, med7 + std7, med7 - std7, color = 'blue', alpha = 0.5)
plt.errorbar(teff2, diff_u2_co_a, yerr = [u2_jn, u2_jp], fmt='.', elinewidth=1, alpha=0.5, color='cornflowerblue',zorder=5, label = 'ATLAS LDCs')

popt6, pcov6 = cft(line, teff2, diff_u2_co_p)
popt7, pcov7 = cft(line, teff2, diff_u2_co_a)

plt.plot(teff2, line(teff2, *popt6), color='orangered', linestyle=':')
plt.plot(teff2, line(teff2, *popt7), color='cornflowerblue', linestyle=':')

plt.plot(x, y, 'k--')
plt.grid()
plt.legend(loc='best')
plt.ylabel(r'Residuals in $u_2$ - EJ(2015) LDCs')
plt.xlabel('Effective Temperature')
plt.savefig(path1 + '/Results/variation_with_temp/u2_code_te.pdf')
plt.close(fig6)

print("------------------------------------------------------------------------")
print("---------------------Your task is complete!!----------------------------")
print("------------------------------------------------------------------------")
