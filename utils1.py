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
