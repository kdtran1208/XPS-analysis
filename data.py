# XPS Data Analysis
# Duyen Tran 07-15-2016
#

#========================================================================================#
#                               START OF CODE
#========================================================================================#

import os,sys
import numpy as np
from matplotlib import *
from math import factorial
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as integrate

def readxps(filename):
	'''
	Function to read raw data lines by lines and assign them to variable "lines"
	In: name of the raw data file --> this is what Avantage gives us after overlapping all spectrums of all samples
	Out: one variable storing raw data --> can be accessed lines by lines
	'''
	fo = open(filename)
	lines = fo.readlines()
	#print "Name of the file",fo.name
	return lines

def spectrum_finder(lines):
	''' 
	Function to extract headers of all elemental scans and assign to variable "plots"
	'''
	for i in lines:
		if "Overlay Axis (Pos)" in i:
			hi= i
			plots=hi.strip().split('\t')[1:]
	return plots

def scannum_finder(plots):
	'''
	Function to find the number of spectra plots of each test cases, including survey
	'''
	index = [i for i, x in enumerate(plots) if x == 'Survey']
	scan_num = index[1] - index[0]
	return scan_num

def element_finder(plots):
	index = [i for i, x in enumerate(plots) if x == 'Survey']
	elements = plots[index[0]+1 : index[1]]
	return elements
# ------------------------------------- #
# --- Importing treatment condition --- #
# ------------------------------------- #

def treatment_condition(ID):
	'''Function to sort out treatment condition:
	surface
	type of deposition
	treatment 
	temperature
	time
	surface cleanning
	'''
	''' 
	in: coupon ID excel file
	out: dictionary sorted by surface/deposition/treatment/temp/time
	'''
	# %%%%% extract data %%%%%% #
	# %%%%%%%%%%%%%%%%%%%%%%%%% #
	# ex: dict['Surface']
	# ID must be an excel file followed the following format "Name of file.xlsx"
	
	df = pd.read_excel(ID)
	dict = pd.DataFrame.to_dict(df)
	return dict

# ------------------------------------------------- #
# --- Arrange test condition coupons by coupons --- #
# ------------------------------------------------- #
def sort_condition(dict):
	'''
	Function to rearrange raw treatment condition into coupons by coupons 
	'''
	'''
	In: dictionary created by function treatment_condition
	Out: list of test condition
	'''
	# ex: condition[1] = ['deposition','surface','temp','time','treatment']
	condition =[]
	for i in range(len(dict['surface'])):
		condition.append({'surface':dict['surface'][i],'deposition':dict['deposition'][i],'treatment':dict['treatment'][i],'temp':dict['temp'][i],'time':dict['time'][i],
		'surface cleaning':dict['surface cleaning'][i]})
		#print condition
	return condition


def check(condition,plots,scan_num):
	'''
	Function to double check whether everything matches up
	The number of test condition (number of coupons) must be equal to number of plots/scan_num
	'''
	'''
	In: condition, plots, number of scan_num
	Out: permission to go ahead or warning to re-check raw data file/ coupon ID file
	'''
	if len(condition) != len(plots)/scan_num:
		print "Check the ID file or raw data file again for correct number of samples"
	else:
		print "Everything is correct. Please proceed forward"
		
# ------------------------------------------------- #
# --- Arrange test condition coupons by coupons --- #
# ------------------------------------------------- #

def sort_keynames(condition,plots,scan_num):
	'''
	Function to include treating condition and elemental scan for each test case
	'''
	'''
	In: condition, plots, number of scan
	Out: compiling all necessary information for each test case
		condition: deposition / surface / temp / time / treatment
		scan : survey / Si2p / P2p / C1s / N1s / Ti2p / O1s
	'''
	keynames =[]
	for i in range(len(condition)):
		#for j in range(scan_num):
		keynames.append({'condition':condition[i],'scan':plots[scan_num*i:i*scan_num+scan_num]})
			#keynames['scan']=plots[scan_num*i:i*scan_num+scan_num]
			#print keynames
	return keynames
	
# change element scan from 1 list to 7 lists. Each list is one scan

def compile (keynames,scan_num):
	''' 
	Function to convert elemental scan to convert elemental scan from string to list
	so we can assign x (binding energy) and y (counts per seconds) later for plotting purposes
	'''
	'''
	In: keynames, number of scan
	Out: compiling all necessary information for each test case
		condition: deposition / surface / temp / time / treatment
		scan : survey / Si2p / P2p / C1s / N1s / Ti2p / O1s -------> in form of "list" instead of string
	'''
	new_key = []
	for i in range(len(keynames)):
		new_key.append({'condition': keynames[i]['condition'].copy(),'scan':[]})
		# newkey['scan']=[]
		for j in range(scan_num):
			#newkey = key[i]['condition'].copy()
			#newkey[i]['scan'] = []
			if j==0:
				#new_key.append({'scan':key[i]['scan'][j].split()})
				new_key[i]['scan'].append(keynames[i]['scan'][j].split())
			else:
				#new_key.append({'scan':key[i]['scan'][j].replace(" ","").split()})
				new_key[i]['scan'].append(keynames[i]['scan'][j].replace(" ","").split())
				#data.append(new_key)
	return new_key

def init_data(new_key,scan_num):
	'''
	Function to initialize x and y value for each elemental scan
	'''
	'''
	In: new_key --> elemental scans are stored as list
	Out: assigned initial value to x and y for each elemental scan
	'''
	for i in range(len(new_key)):
		for j in range(scan_num):
				new_key[i]['scan'][j].append({'binding energy (ev)':[],'counts per sec':[]})
	return new_key
	
def stackdata(lines,new_key,scan_num):	
	energy = []
	split = []
	'''
	This loop is to strip all unnecessary information out of variable "lines"
	Will retain info that we want and store in variable "split"
	'''
	for i in range(len(lines)):
		if '#INF0' in lines[i]:
			split.append(lines[i].strip().split('\t'))
	''' 
	This loop is to extract binding energy and assign to the variable "energy" 
	All scans of all samples will have the same x value which is energy 
	'''		
	for i in range(len(split)):
		energy.append (split[i][0])
		
	'''
	This loop is to assign "energy" variable to the x value of all scan
	for all test cases
	'''
	for m in range(len(new_key)):
		for n in range(scan_num):
			new_key[m]['scan'][n][1]['binding energy (eV)'] =  energy
			
	'''
	This loop is to extract signal values ("y") (counts per second) of each scan
	for each test case and assign them in the "new_key" variable
	'''
	for i in range(len(split)):
		#for j in range(len(split[i])):
		for m in range(len(new_key)):
			for n in range(scan_num):
			#		if j == 0 :
			#			continue
			#		elif split[i][j] == '':
			#			continue
			#		else:
						new_key[m]['scan'][n][1]['counts per sec'].append([split[i][scan_num*m+n+2]])
	return new_key

def cutdata(new_key, scan_num):
	'''
	Function to trim the placeholder '-1.#INF0' off the y value
	In: Data with x and y extracted from raw data file
	Out: Trim out the placeholder '-1.#INF0' off some y values
	'''
	for i in range(len(new_key)):
		for j in range(scan_num):
			#for k in range(len(new_key[i]['scan'][j][1]['counts per sec'])):
			index = [m for m, x in enumerate(new_key[i]['scan'][j][1]['counts per sec']) if x!= ['-1.#INF0']]
			new_key[i]['scan'][j][1]['binding energy (eV)'] = new_key[i]['scan'][j][1]['binding energy (eV)'][index[0]:index[-1]]
			new_key[i]['scan'][j][1]['counts per sec'] = new_key[i]['scan'][j][1]['counts per sec'][index[0]:index[-1]]
			# for k in range(len(new_key[i]['scan'][j][1]['binding energy (eV)'])):
				# new_key[i]['scan'][j][1]['binding energy (eV)'][k] = float(new_key[i]['scan'][j][1]['binding energy (eV)'][k])
				# new_key[i]['scan'][j][1]['counts per sec'][k] = float(''.join(new_key[i]['scan'][j][1]['counts per sec'][k]))
	return new_key
	
def convertdata(new_key,scan_num):
	'''
	Function to convert x and y values to 'float' for ease of calculation later
	'''
	for i in range(len(new_key)):
		for j in range(scan_num):
			for k in range(len(new_key[i]['scan'][j][1]['counts per sec'])):
					new_key[i]['scan'][j][1]['binding energy (eV)'][k] = float(new_key[i]['scan'][j][1]['binding energy (eV)'][k])
					new_key[i]['scan'][j][1]['counts per sec'][k] = float(''.join(new_key[i]['scan'][j][1]['counts per sec'][k]))
	return new_key

def energy_finder(new_key,scan_num):
	'''
	Function to extract binding energy of all test case for ease of access/plotting later
	This is just a support function --> don't really need to use
	'''
	energy = []
	for i in range(len(new_key)):
		for j in range(scan_num):
			energy.append(new_key[i]['scan'][j][1]['binding energy (eV)'])
	return energy

def rawsignal_finder(new_key,scan_num):
	'''
	Function to extract 'counts per sec' value of all test case for ease of access/plotting later
	This is just a support function --> don't really need to use
	'''
	rawsignal = []
	for i in range(len(new_key)):
		for j in range(scan_num):
			rawsignal.append(new_key[i]['scan'][j][1]['counts per sec'])
	return rawsignal

def background_finder(new_key, energy , scan_num):
	'''
	Function to find background noise and substract it out
	In: Data with x and y without the placeholder
	Out: 
		Background signal
		Real signal (substracted background)
	'''
	background = []
	coeff = []
	for i in range(len(new_key)):
		for j in range(scan_num):
			x = [ new_key[i]['scan'][j][1]['binding energy (eV)'][0] , new_key[i]['scan'][j][1]['binding energy (eV)'][-1] ]
			x = np.asarray(x)
			y = [ new_key[i]['scan'][j][1]['counts per sec'][0] , new_key[i]['scan'][j][1]['counts per sec'][-1] ]
			y = np.asarray(y)
			coeff.append(np.polyfit(x,y,1))
			#energy.append(new_key[i]['scan'][j][1]['binding energy (eV)'])
	for i in range(len(energy)):
		#for j in range(len(energy[i])):
		yspan = (coeff[i][0]*np.asarray(energy[i]) + coeff[i][1]).tolist()
		background.append(yspan)
	return background
	
def background_substractor(new_key,background,scan_num):
	'''
	Function to substract background signal from raw data and produce real signal
	'''
	realsignal = []
	for i in range(len(new_key)):
		for j in range(scan_num):
			raw = new_key[i]['scan'][j][1]['counts per sec']
			substract = (np.abs(np.asarray(raw) - np.asarray(background[scan_num*i + j]))).tolist()
			realsignal.append(substract)
	return realsignal

def perc_atomic(new_key, scan_num, realsignal , energy ):
	'''
	Function to calculate the atomic percentage of each elements for each test case
	We will calculate area under the curve (trapezoidal rule)
		--> convert that to % atomic
	In:
		real signal --> should equal len(new_key)*len(scan_num))
		energy 		--> should equal len(new_key)*len(scan_num))

	Out:
		atomic percentage of each atom in one test case
		ex: in test case 1: we are interested in 6 element: Si, Ti, O, P, N, C
		--> this function will give the atomic percentage of each element
		based on the ratio of their area under the curve
	'''
	percentage 			= []
	split_area 			= []
	split_percentage 	= []
	area = []
	
	for i in range(len(new_key)):
		for j in range(scan_num)[1:]:
			area.append(integrate.trapz(realsignal[scan_num*i+j]))
		
	for i in range(len(new_key)):
		split_area.append(area[(scan_num-1)*i : (scan_num-1)*(i+1)])
		for j in range(len(split_area[i])):
			percentage.append(split_area[i][j]*100 / sum(split_area[i]))
		split_percentage.append(percentage[(scan_num-1)*i : (scan_num-1)*(i+1)])
	return split_percentage

def exportcsv(new_key, background, realsignal, elements , split_percentage , scan_num):
	csv_1 = {} # this is to store all raw signals
	csv_2 = {} # this is to store atomic percentage of each element
	total = []
	length = []
	totalvalue = []
	split = []
	colname = ['00-surface' , '01-deposition' , '02-treatment' , '03-temp' , '04-time' , '05-surface cleaning' ,
			   '06-scan' , '07-binding energy (eV)' , '08-raw (cps)' ,'09-background (cps)' , '10-real (cps)'] 
	
	for i in range(len(new_key)):
		for j in range(scan_num):	
			length = len(new_key[i]['scan'][j][1]['counts per sec'])
			total.append(length)
		split.append(total[(scan_num*i) : (scan_num*i + scan_num)])
		totalvalue.append(sum(split[i]))
		
	surface = []
	deposition = []
	treatment = []
	temp = []
	time = []
	surfacecleaning = []
	scan = []
	bindingenergy =[]
	raw = []
	control = []
	real = []
	for i in range(len(new_key)):
		#for j in range(scan_num):
			surface.extend([str(new_key[i]['condition']['surface'])] * totalvalue [i])
			deposition.extend([str(new_key[i]['condition']['deposition'])] * totalvalue[i])
			treatment.extend([str(new_key[i]['condition']['treatment'])]*totalvalue[i])
			temp.extend([str(new_key[i]['condition']['temp'])] * totalvalue[i])
			time.extend([str(new_key[i]['condition']['time'])] * totalvalue[i])
			surfacecleaning.extend([str(new_key[i]['condition']['surface cleaning'])] * totalvalue[i])
			for j in range(scan_num):
				scan.extend([new_key[i]['scan'][j][0]]*total[scan_num*i + j])
				bindingenergy.extend(new_key[i]['scan'][j][1]['binding energy (eV)'])
				raw.extend(new_key[i]['scan'][j][1]['counts per sec'])
				control.extend(background[scan_num*i+j])
				real.extend(realsignal[scan_num*i+j])
			
	csv_1.update({colname[0]:surface    , colname[1]:deposition    , colname[2]:treatment , 
				  colname[3]:temp       , colname[4]:time          , colname[5]:surfacecleaning ,
				  colname[6]:scan    	, colname[7]:bindingenergy , colname[8]:raw , 
				  colname[9]:control    , colname[10]:real})
	
	percentage = [[] for i in range(len(elements))]
	for i in range(len(percentage)):
		for j in range(len(split_percentage)):
			percentage[i].append(split_percentage[j][i])
		csv_2.update({elements[i]:percentage[i]})
	
	df1 = pd.DataFrame.from_dict(csv_1)
	df2 = pd.DataFrame.from_dict(csv_2)
	df1.to_csv('XPS_rawdata.csv')
	df2.to_csv('XPS_atomicpercentage.csv')
	
#========================================================================================#
#                               END OF CODE
#========================================================================================#

				
# def savitzky_golay(new_key, window_size, order, deriv=0, rate=1):
    # try:
        # window_size = np.abs(np.int(window_size))
        # order = np.abs(np.int(order))
    # except ValueError, msg:
        # raise ValueError("window_size and order have to be of type int")
    # if window_size % 2 != 1 or window_size < 1:
        # raise TypeError("window_size size must be a positive odd number")
    # if window_size < order + 2:
        # raise TypeError("window_size is too small for the polynomials order")
    # order_range = range(order+1)
    # half_window = (window_size -1) // 2
    # # precompute coefficients
    # b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    # m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # # pad the signal at the extremes with
    # # values taken from the signal itself
	# for i in range(len(new_key)):
		# for j in range(scan_num):
			# y = new_key[i]['scan'][j][1]['counts per sec']
		
	# firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
	# lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
	# y = np.concatenate((firstvals, new_key[i]['scan'][j][1]['counts per sec'], lastvals))
    # return np.convolve( m[::-1], new_key[i]['scan'][j][1]['counts per sec'], mode='valid')

