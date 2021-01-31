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
f1 = open(path1 + '/Phoenix/claret_us_nl_pho.dat','w')#-------------------Non-linear-s-------------
f1.write('#Name\t\tc1\t\t\tc2\t\t\tc3\t\t\tc4\n')

f1r = open(path1 + '/Phoenix/claret_us_nl_pho_r.dat','w')#--------------------Non-linear-r------------
f1r.write('#Name\t\tc1\t\t\tc2\t\t\tc3\t\t\tc4\n')

f2 = open(path1 + '/Phoenix/claret_limiting_LDC_pho.dat','w')#-------------Limiting----------------
f2.write('#Name\t\tu1\t\tu2\n')

f2r = open(path1 + '/Phoenix/claret_limiting_LDC_pho_r.dat','w')#------------Limiting-r---------------
f2r.write('#Name\t\tu1\t\tu2\n')

#--------------------------------------------------------------------------------------------------
#--------------------------------Claret (2017) ATLAS LDCs------------------------------------------
#--------------------------------------------------------------------------------------------------
f3 = open(path1 + '/Atlas/claret_us_nl_ata.dat','w')#---------------------Non-linear---------------
f3.write('#Name\t\tc1\t\t\tc2\t\t\tc3\t\t\tc4\n')

f4 = open(path1 + '/Atlas/claret_limiting_LDC_ata.dat','w')##-------------Limiting------------------
f4.write('#Name\t\tu1\t\tu2\n')

#--------------------------------------------------------------------------------------------------
#----------------------------------Code ATLAS LDCs-------------------------------------------------
#--------------------------------------------------------------------------------------------------
f33 = open(path1 + '/Atlas/code_us_nl_ata.dat','w')#-------------------Non-linear------------------
f33.write('#Name\t\tc1\t\t\tc2\t\t\tc3\t\t\tc4\n')

f44 = open(path1 + '/Atlas/code_limiting_LDC_ata.dat','w')##-----------------Limiting---------------
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
		print('******************************************************************************************')
		print('                                                                                          ')
		print('Completed analysis for ' + str(i+1) + ' system(s) / out of ' + str(len(name)) + ' systems')
		print('                                                                                          ')
		print('******************************************************************************************')
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
#ra2, dec2 = np.loadtxt('data2.dat', dtype = str, usecols = (18,19), unpack = True)

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
	#----Calculating LDCs from PHOENIX-s---
	#--------from Claret(2017)-------------
	#--------------------------------------
	lg1_p, T1_p, q1_p, q2_p, q3_p, q4_p = np.loadtxt(path1 + '/Phoenix/pho_ldc_claret.dat', usecols=(0,1,4,5,6,7), unpack=True)
	pts_p = np.vstack((lg1_p, T1_p))
	pt_p = np.transpose(pts_p)
	c1_p = inp.griddata(pt_p, q1_p, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	c2_p = inp.griddata(pt_p, q2_p, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	c3_p = inp.griddata(pt_p, q3_p, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	c4_p = inp.griddata(pt_p, q4_p, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	f1.write(name2[i] + '\t' + str(c1_p) + '\t' + str(c2_p) + '\t' + str(c3_p) + '\t' + str(c4_p) + '\n')
	u1_p = ((12*c1_p)/35) + c2_p + ((164*c3_p)/105) + (2*c4_p)
	u2_p = ((10*c1_p)/21) - ((34*c3_p)/63) - c4_p
	f2.write(name2[i] + '\t' + str(u1_p) + '\t' + str(u2_p) + '\n')
	#--------------------------------------
	#----Calculating LDCs from PHOENIX-r---
	#--------from Claret(2017)-------------
	#--------------------------------------
	lg1_pr, T1_pr, q1_pr, q2_pr, q3_pr, q4_pr = np.loadtxt(path1 + '/Phoenix/pho_ldc_claret_r.dat', usecols=(0,1,4,5,6,7), unpack=True)
	pts_pr = np.vstack((lg1_pr, T1_pr))
	pt_pr = np.transpose(pts_pr)
	c1_pr = inp.griddata(pt_pr, q1_pr, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	c2_pr = inp.griddata(pt_pr, q2_pr, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	c3_pr = inp.griddata(pt_pr, q3_pr, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	c4_pr = inp.griddata(pt_pr, q4_pr, (lg2[i],teff2[i]), fill_value = 0, method = 'cubic')
	f1r.write(name2[i] + '\t' + str(c1_pr) + '\t' + str(c2_pr) + '\t' + str(c3_pr) + '\t' + str(c4_pr) + '\n')
	u1_pr = ((12*c1_pr)/35) + c2_pr + ((164*c3_pr)/105) + (2*c4_pr)
	u2_pr = ((10*c1_pr)/21) - ((34*c3_pr)/63) - c4_pr
	f2r.write(name2[i] + '\t' + str(u1_pr) + '\t' + str(u2_pr) + '\n')
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
	f3.write(name2[i] + '\t' + str(c1_a) + '\t' + str(c2_a) + '\t' + str(c3_a) + '\t' + str(c4_a) + '\n')
	u1_a = ((12*c1_a)/35) + c2_a + ((164*c3_a)/105) + (2*c4_a)
	u2_a = ((10*c1_a)/21) - ((34*c3_a)/63) - c4_a
	f4.write(name2[i] + '\t' + str(u1_a) + '\t' + str(u2_a) + '\n')
	print('****************************************************************************************')
	print('                                                                                        ')
	print('Calculated LDCs for ' + str(i+1) + ' system(s) / out of ' + str(len(name2)) + ' systems\n')
	print('                                                                                        ')
	print('****************************************************************************************')


f1.close()
f2.close()
f1r.close()
f2r.close()
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
u2_code_a = ((10*c1_code_a)/21) - ((34*c3_code_a)/63) - c4_code_a

for i in range(len(name1)):
	f44.write(name1[i] + '\t\t' + str(u1_code_a[i]) + '\t\t' + str(u2_code_a[i]) + '\n')

f44.close()

#--------------------------------------------------------------------------------------------------
#----------------Calculating Limiting LDCs from PHOENIX (Code)-------------------------------------
#--------------------------------------------------------------------------------------------------

name11 = np.loadtxt(path1 + '/Phoenix/code_us_nl_pho.dat',dtype=str,usecols=0,unpack=True)
c1_code_p, c2_code_p, c3_code_p, c4_code_p = np.loadtxt(path1 + '/Phoenix/code_us_nl_pho.dat', usecols=(1,2,3,4), unpack=True)

u1_code_p = ((12*c1_code_p)/35) + c2_code_p + ((164*c3_code_p)/105) + (2*c4_code_p)
u2_code_p = ((10*c1_code_p)/21) - ((34*c3_code_p)/63) - c4_code_p

for i in range(len(name11)):
	f22.write(name11[i] + '\t\t' + str(u1_code_p[i]) + '\t\t' + str(u2_code_p[i]) + '\n')

f22.close()

print("----------------------------------------------------------------------------------")
print("---------------------Your computing task is complete!!----------------------------")
print("----------------------------------------------------------------------------------")

print('----------------------------------------------------------------------------------')
print('------------------Now you can to plot some amazing results------------------------')
print('--------------------Please run plot.py to start plotting--------------------------')
print('----------------------------------------------------------------------------------')
