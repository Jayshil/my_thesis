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

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def lowest_bic(a):
	"""
	Find the lowest BIC (Bayesian Inference Crterian) among
	three BICs.
	I know I can just look at them manually but I am a little-
	bit lazy to do so. Hence here I am with a function defined for it.
	
	Parameters
	----------
	a: dict
		Values of BICs in dict.
	n: int
		Degree of highest polynomial of model.
	
	returns: dict
		Containing the best fit model, with its BIC.
		If two models have BIC difference less than 2,
		then the model with a fewer number of parameters
		would be selected.
	"""
	b = np.array([])
	for i in a.values():
		b = np.hstack((b,i))
	def get_key(val, my_dict):
		for key, value in my_dict.items():
			if val == value:
				return key
		return "Key does not exist"
	xx = np.min(b)
	yy = get_key(xx, a)
	ret = {}
	con = np.abs(a[yy] - a['constant'])
	lin = np.abs(a[yy] - a['linear'])
	qua = np.abs(a[yy] - a['quadratic'])
	if yy == 'quadratic':
		if con < 2:
			ret['constant'] = a['constant']
		elif lin < 2 and lin < con:
			ret['linear'] = a['linear']
		else:
			ret['quadratic'] = a['quadratic']
	if yy == 'linear':
		if con < 2:
			ret['constant'] = a['constant']
		else:
			ret['linear'] = a['linear']
	if yy == 'constant':
		ret['constant'] = a['constant']
	return ret

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def binned_data(datax, datay, nos=10, datax_err=None, datay_err=None):
	"""
	This function creates binned array from the given array.
	
	Parameters
	----------
	datax: np.ndarray
		One dimensional array. x-coordinate
	datax_err: np.ndarray
		One dimensional array. Error in x-coordinate
		If not provided then assumes to be a zero matrix
		of length of datax.
	datay: np.ndarray
		One dimensional array. y-coordinate
	datay_err: np.ndarray
		One dimensional array. Error in y-coordinate
		If not provided then assumes to be a zero matrix
		of length of datay.
	nos: float
		Number of binned data points you want to set.
		Default is 10.

	returns: {np.ndarray}
		numpy array of binned data
	"""
	if datax_err is None:
		datax_err = np.zeros(len(datax))
	if datay_err is None:
		datay_err = np.zeros(len(datay))
	aa = []
	for i in range(len(datax)):
		xxx = (datax[i], datax_err[i], datay[i], datay_err[i])
		aa.append(xxx)
	bb = sorted(aa)
	aaa = bbb = ccc = ddd = np.array([])
	for i in range(len(bb)):
		aaa = np.hstack((aaa, bb[i][0]))
		bbb = np.hstack((bbb, bb[i][1]))
		ccc = np.hstack((ccc, bb[i][2]))
		ddd = np.hstack((ddd, bb[i][3]))
	rep = int((len(datax))/(nos-1))
	rem = len(datax) - ((nos-1)*rep)
	bin_datax = bin_datax_err = bin_datay = bin_datay_err = np.array([])
	k = 0
	for i in range(nos-1):
		du_t = np.zeros(1000)
		for j in range(rep):
			abc1 = np.random.normal(aaa[k], bbb[k], 1000)
			du_t = np.vstack((du_t, abc1))
			k = k+1
		bint = np.mean(du_t[1:], axis=0)
		bin_datax = np.hstack((bin_datax, np.median(bint)))
		bin_datax_err = np.hstack((bin_datax_err, np.std(bint)))
	rem_t = np.zeros(1000)
	for i in range(rem):
		abc1 = np.random.normal(aaa[k], bbb[k], 1000)
		rem_t = np.vstack((rem_t, abc1))
	remt = np.mean(rem_t[1:], axis=0)
	bin_datax = np.hstack((bin_datax, np.median(remt)))
	bin_datax_err = np.hstack((bin_datax_err, np.std(remt)))
	k1 = 0
	for i in range(nos-1):
		du_d = np.zeros(1000)
		for j in range(rep):
			abc1 = np.random.normal(ccc[k1], ddd[k1], 1000)
			du_d = np.vstack((du_d, abc1))
			k1 = k1+1
		bind = np.mean(du_d[1:], axis=0)
		bin_datay = np.hstack((bin_datay, np.median(bind)))
		bin_datay_err = np.hstack((bin_datay_err, np.std(bind)))
	rem_d = np.zeros(1000)
	for i in range(rem):
		abc1 = np.random.normal(ccc[k1], ddd[k1], 1000)
		rem_d = np.vstack((rem_d, abc1))
	remd = np.mean(rem_d[1:], axis=0)
	bin_datay = np.hstack((bin_datay, np.median(remd)))
	bin_datay_err = np.hstack((bin_datay_err, np.std(remd)))
	return bin_datax, bin_datax_err, bin_datay, bin_datay_err
