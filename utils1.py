import numpy as np
import os
import re

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

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#------------------This Function is taken from the 'Exotoolbox' by ---------------------------
#--------------------------------Dr Nestor Espinoza-------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def get_quantiles(dist,alpha = 0.68, method = 'median'):
    """
    get_quantiles function

    DESCRIPTION

        This function returns, in the default case, the parameter median and the error% 
        credibility around it. This assumes you give a non-ordered 
        distribution of parameters.

    OUTPUTS

        Median of the parameter,upper credibility bound, lower credibility bound

    """
    ordered_dist = dist[np.argsort(dist)]
    param = 0.0
    # Define the number of samples from posterior
    nsamples = len(dist)
    nsamples_at_each_side = int(nsamples*(alpha/2.)+1)
    if(method == 'median'):
       med_idx = 0
       if(nsamples%2 == 0.0): # Number of points is even
          med_idx_up = int(nsamples/2.)+1
          med_idx_down = med_idx_up-1
          param = (ordered_dist[med_idx_up]+ordered_dist[med_idx_down])/2.
          return param,ordered_dist[med_idx_up+nsamples_at_each_side],\
                 ordered_dist[med_idx_down-nsamples_at_each_side]
       else:
          med_idx = int(nsamples/2.)
          param = ordered_dist[med_idx]
          return param,ordered_dist[med_idx+nsamples_at_each_side],\
                 ordered_dist[med_idx-nsamples_at_each_side]

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

def RA_to_deg(coords):
    """
    Given a RA string in hours (e.g., '11:12:12.11'), returns the corresponding 
    coordinate in degrees.
    """
    hh,mm,ss = coords.split(':')

    hours = np.double(hh) + np.double(mm)/60. + np.double(ss)/3600.
    return hours * 360./24.

def DEC_to_deg(coords):
    """
    Given a DEC string in degrees (e.g., '-30:12:12.11'), returns the corresponding 
    coordinate in degrees.
    """
    dd,mm,ss = coords.split(':')
    if dd[0] == '-':
        return np.double(dd) - np.double(mm)/60. - np.double(ss)/3600.
    else:
        return np.double(dd) + np.double(mm)/60. + np.double(ss)/3600.

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
