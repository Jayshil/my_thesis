import numpy as np
import os
import re
from scipy import stats

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------


def new_dir(path):
	os.system('mkdir ' + '-p' + ' ' + path)
	return

def copy_file(in_path, fi_name, out_path):
	if not os.path.exists(out_path):
		new_dir(out_path)
	os.system('cp ' + in_path + fi_name + ' ' + out_path)
	return

def move_file(in_path, fi_name, out_path, new_name):
	if not os.path.exists(out_path):
		new_dir(out_path)
	os.system('mv ' + in_path + fi_name + ' ' + out_path + new_name)
	#os.system('rm ' + in_path + fi_name)
	return

#------------------------------------------------------------------------------------------
#-------------------------------Natural Sorting--------------------------------------------
#------------------------------------------------------------------------------------------
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def convert_bp(r1,r2,pl,pu):
    Ar = (pu - pl)/(2. + pl + pu)
    nsamples = len(r1)
    p = np.zeros(nsamples)
    b = np.zeros(nsamples)
    for i in range(nsamples):
        if r1[i] > Ar:
            b[i],p[i] = (1+pl)*(1. + (r1[i]-1.)/(1.-Ar)),\
                        (1-r2[i])*pl + r2[i]*pu
        else:
            b[i],p[i] = (1. + pl) + np.sqrt(r1[i]/Ar)*r2[i]*(pu-pl),\
                        pu + (pl-pu)*np.sqrt(r1[i]/Ar)*(1.-r2[i])
    return b,p

def convert_ld_coeffs(ld_law, coeff1, coeff2):
    if ld_law == 'quadratic':
        q1 = (coeff1 + coeff2)**2
        q2 = coeff1/(2.*(coeff1+coeff2))
    elif ld_law=='squareroot':
        q1 = (coeff1 + coeff2)**2
        q2 = coeff2/(2.*(coeff1+coeff2))
    elif ld_law=='logarithmic':
        q1 = (1-coeff2)**2
        q2 = (1.-coeff1)/(1.-coeff2)
    return q1,q2

def reverse_ld_coeffs(ld_law, q1, q2):
    if ld_law == 'quadratic':
        coeff1 = 2.*np.sqrt(q1)*q2
        coeff2 = np.sqrt(q1)*(1.-2.*q2)
    elif ld_law=='squareroot':
        coeff1 = np.sqrt(q1)*(1.-2.*q2)
        coeff2 = 2.*np.sqrt(q1)*q2
    elif ld_law=='logarithmic':
        coeff1 = 1.-np.sqrt(q1)*q2
        coeff2 = 1.-np.sqrt(q1)
    elif ld_law=='linear':
        return q1,0.
    return coeff1,coeff2

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def freedman_diaconis(data, returnas="width"):
	"""
	Use Freedman Diaconis rule to compute optimal histogram bin width. 
	``returnas`` can be one of "width" or "bins", indicating whether
	the bin width or number of bins should be returned respectively.

	Parameters
	----------
	data: np.ndarray
		One-dimensional array.

	returnas: {"width", "bins"}
		If "width", return the estimated width for each histogram bin. 
		If "bins", return the number of bins suggested by rule.
	"""
	data = np.asarray(data, dtype=np.float_)
	IQR = stats.iqr(data, rng=(25, 75), scale="raw", nan_policy="omit")
	N = data.size
	bw = (2 * IQR) / np.power(N, 1/3)
	if returnas=="width":
		result = bw
	else:
		datmin, datmax = data.min(), data.max()
		datrng = datmax - datmin
		result = int((datrng / bw) + 1)
	return(result)
