#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
#from __future__ import division
"""Mass Balance Calculator"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '1.0.0'
__date__ = '06/11/2015'
#
## IMPORTS
import sys
import os
# import glob
# import fnmatch
# import gc
import csv
# import pickle
# import time
from datetime import datetime, timedelta
# import subprocess
# import urllib2
# import os.path
# import zipfile
# import StringIO
#
## SCIENTIFIC
# import math
from math import *
import numpy as np
# import numpy.ma as ma
# import scipy.linalg
# import scipy.stats as stats
from scipy.stats.stats import nanmean
# from scipy.stats import norm
# from pandas import Series, DataFrame
# import pandas as pd
#
## GIS SUPPORT
# import gdal
# import gdalconst
# from geostatsmodels import utilities, kriging, variograms, model, geoplot
# from qgis.core import QgsProject
# from PyQt4.QtCore import QFileInfo
# To run in QGIS enter following in QGIS console:
# execfile(u'/path/to/script.py'.encode('utf-8'))
#
## PLOT TOOLS
# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.dates as dates
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.backends.backend_pdf import PdfPages
# import pylab
# pylab.ion()
# pylab.ioff()
#
#
from standard import *
from spatialfiles import *
from kriging import *
#
## FUNCTIONS
#
def envset(*args):
	'''Set up working folders and files via parameters found in a text file
	Settings file should look like this (omit or # those not needed).
	Paths can be absolute (from root /Users/username) or relative (./path).
	**********************************************************************
	Input Data Folder =./InData/
	Output Folder =./OutData
	Probing = SG_Probe_2013.csv
	Density = SG_Density_2013.csv
	Stakes = SG_Stakes_2013.csv
	DEM=./InData/2010_DEM2.tif
	#Data Files =SG_Probe_2013.csv
	Plot DEM? =y
	Plot results? =y
	#Base =
	Detrend = n'''
	#Called by:
	env = os.environ
	curLoc = env['PWD']
	curList = os.listdir('./')
	if len(args)>0:
		if os.path.exists(args[0]):
			envFile = args[0]
	else:
		print '\nList of files in ', curLoc,'\n'
		for i in curList:
			if not os.path.isdir(i):
				print i
		envFile = './settings.txt'
		while not os.path.exists(envFile):
			envFile = raw_input('\nEnter name of file containing environment settings(txt): ')
	settings = {}
	InFile = open(envFile,'rb')
	# Check contents and set up dictionary
	for row in InFile:
		line = row.strip().split('=')
		print line
		settings[line[0].strip()]=line[1].strip()
	# Folder containing data files
	if 'Input Data Folder' not in settings:
		settings['Input Data Folder'] = './Indata/.'
	inDataDir = settings['Input Data Folder']
	if not os.path.exists(inDataDir):
		sys.exit('Input Data Folder setting incorrect')
	# Root folder for output data
	outDir = settings['Output Folder']
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	# Name namstr used for output subfolder, set first to day_time
	namstr = datetime.now().strftime('%j_%H%M%S')
	destDir = os.path.join(outDir,namstr)
	if 'Data Files' in settings:
		fileList = settings['Data Files'].split(',')
		settings['Data Files'] = fileList
		for inDataName in fileList:
			if not os.path.exists(os.path.join(inDataDir,inDataName)):
				sys.exit('Data File setting incorrect')
	if 'Probing' in settings:
			if not os.path.exists(os.path.join(inDataDir,settings['Probing'])):
				sys.exit('Probing setting incorrect')
	if 'Density' in settings:
		print settings['Density']
		densityList = settings['Density'].split(',')
		settings['Density'] = densityList
		print settings['Density']
		print densityList
		for densityName in densityList:
			if not os.path.exists(os.path.join(inDataDir,densityName)):
				sys.exit('Density setting incorrect')
	if 'Stakes' in settings:
			if not os.path.exists(os.path.join(inDataDir,settings['Stakes'])):
				sys.exit('Stakes setting incorrect')
	# Create output folder for session
	if not os.path.exists(destDir):
		os.makedirs(destDir)
	settings['Write Folder']=destDir
	if 'DEM' in settings:
		demFile = settings['DEM']
		if not os.path.exists(demFile):
			print settings['DEM']
			sys.exit('DEM file setting incorrect')
	if 'MSK' in settings:
		demFile = settings['MSK']
		if not os.path.exists(demFile):
			print settings['MSK']
			sys.exit('MSK (mask) file setting incorrect')
	if 'Pickle' in settings:
		pickleFile = settings['Pickle']
		if not os.path.exists(pickleFile):
			sys.exit('No such pickle file')
	settings['This File'] = envFile
	print '\nSettings are: '
	for i in settings.keys():
		print i,': ',settings[i]
	print '\n'*2
	return settings, envFile
#
def cont():
	'''Simple function to break script flow'''
	#Called by:
	print '\n'*2,'*'*10
	answers = ['c','y','s','e']
	a = ''
	while a not in answers:
		a = raw_input('\nDo you wish to continue (y or c), skip (s) or end (e)?: ')
	if a == 'e':
		sys.exit(0)
	elif a == 's':
		ans = 1
	elif a == 'c' or a == 'y':
		ans = 0
	print '\n'*2
	return ans
#
def readDensity(settings):
	'''1. READ DENSITY PROFILE DATA FROM CSV AND CREATE DENSITY PROFILE FROM IT
	Requires "settings".'''
	#
	print readDensity.__doc__,'\n','_'*20,'\n'
	time_one = datetime.now()
	print time_one.strftime('\n1. started at day:%j %H:%M:%S')
	# Process files into dictionary
	densityVectors = {}
	for file in settings['Density']:
		dfile = os.path.join(settings['Input Data Folder'],file)
		densityVectors[file] = import2vector(dfile)
	# Create least squares fitted functions
	densityFiles = densityVectors.keys()
	#outDir = os.path.join(settings['Write Folder'],'Density')
	outDir = settings['Write Folder']
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	DPS =  open(os.path.join(outDir,'DensityPoly.txt'), 'a')
	print 'Density files, with headers:'
	for i in densityFiles:
		xvec = 'Cumulative Snow'
		yvec = 'Cumulative Mass'
		print 'File: ',i, '\n'
		for j in densityVectors[i].keys():
			print j
		qus1 = '\nUse '+xvec+' & '+yvec+'? (y)'
		ans1 = raw_input(qus1)
		if ans1 != 'y':
			while xvec not in densityVectors[i].keys():
				xvec = raw_input('Give name of X vector: ')
			while yvec not in densityVectors[i].keys():
				yvec = raw_input('Give name of Y vector: ')
		densityVectors[i]['DensFunc1'] = fitFunc(densityVectors[i][xvec],densityVectors[i][yvec],1)
		densityVectors[i]['DensFunc2'] = fitFunc(densityVectors[i][xvec],densityVectors[i][yvec],2)
		densityVectors[i]['DensFunc3'] = fitFunc(densityVectors[i][xvec],densityVectors[i][yvec],3)
		print i
		print densityVectors[i]['DensFunc1']
		print densityVectors[i]['DensFunc2']
		print densityVectors[i]['DensFunc3']
		DPS.write(str(i)+'\n')
		DPS.write(str(densityVectors[i]['DensFunc1'])+'\n')
		DPS.write(str(densityVectors[i]['DensFunc2'])+'\n')
		DPS.write(str(densityVectors[i]['DensFunc3'])+'\n')
	# Save data dictionary to pickle database
	DPS.close()
	picksave(densityVectors,'density',outDir)
	print '\nOutput to ',settings['Write Folder'],'\n'
	return densityVectors, densityFiles
#
#
def convertProbing(settings, densityVectors, densityFiles):
	'''2. READ PROBING DATA AND CONVERT SNOW DEPTH TO M W.E. USING DENSITY PROFILE - ACCUMULATION
	Requires "settings" and densityVectors, densityFiles from "readDensity" '''
	#
	print convertProbing.__doc__,'\n','_'*20,'\n'
	time_two = datetime.now()
	print time_two.strftime('\n2. started at day:%j %H:%M:%S')
	# Read csv file containing probing data
	File = csv.reader(open(os.path.join(settings['Input Data Folder'],settings['Probing']),'rb'), delimiter=',', quotechar='|')
	# Read first line as header into a list
	colNames = File.next()
	colConv = 'Depth'
	# Ask for column to be converted. Give headers present but offer 'Depth' as default
	print 'Column headers in input file.\n',colNames
	qus2 = 'Use '+colConv+'? (y)'
	ans2 = raw_input(qus2)
	if ans2 != 'y':
		colConv = ''
		while colConv not in colNames:
			colConv = raw_input('Which column is to be converted to m w.e.? ')
	colConvLoc = colNames.index(colConv)
	qus3 = 'Conversion factor for depth values in '+colConv+'\n e.g. 1 for '+colConv+' in metres, 100 for centimetres?: '
	ans3 = raw_input(qus3)
	try:
		convFact = float(ans3)
	except:
		convFact = 1
		print '\n Set conversion factor to 1 as given value not valid.\n'
	# Show the choice of density functions available
	print 'Density functions:'
	for i in densityFiles:
		print 'File: ',i
		keys = densityVectors[i].keys()
		fun2=[]
		for fun1 in keys:
			if 'DensFunc' in fun1:
				print fun1
				fun2.append(fun1)
	# Chose the function to use by file and function (for multiple density pits and fits)
	fileChoice = ''
	while fileChoice not in densityFiles:
		fileChoice = raw_input('\nEnter density file name (or "Enter" for first file): ')
		if fileChoice == '': fileChoice = densityFiles[0]
	funcChoice = ''
	while funcChoice not in fun2:
		print 'Functions: ',fun2
		funcChoice = raw_input('\nEnter function or press "Enter" for DensFunc1: ')
		if funcChoice == '': funcChoice = 'DensFunc1'
	func = densityVectors[fileChoice][funcChoice]
	# Write converted data back to csv file
	outName,_,_,_ = namer(settings['Probing'])
	outDataNamestr = os.path.basename(outName + '_e.csv')
	csvOut = open(os.path.join(settings['Write Folder'],outDataNamestr),'wb')
	#colNames.append(colConv +'X' +	 funcChoice)
	colNames.append('mwe')
	Header = ','.join(colNames) + '\n'
	csvOut.write(Header)
	csvOut = csv.writer(open(os.path.join(settings['Write Folder'],outDataNamestr),'ab'), delimiter = ',', lineterminator='\n')
	for inrow in File:
		outrow = inrow
		outrow.append(str("%.2f" % func(float(outrow[colConvLoc])/convFact)))
		csvOut.writerow(outrow)
	# Weird thing of opening again just to close it but seems to be necessary
	csvOut = open(os.path.join(settings['Write Folder'],outDataNamestr),'ab')
	csvOut.close()
	# Pass file names on to next stage
	print 'Output to ',settings['Write Folder'],'\n'
	sourceFolder = settings['Write Folder']
	sourceFile = outDataNamestr
	return sourceFolder, sourceFile, func
#
#
def krigAcc(settings, sourceFolder, sourceFile):
	'''3. KRIG ACCUMULATION FILE ACROSS DEM
	Requires "settings" and sourceFolder, sourceFile'''
	#
	print krigAcc.__doc__,'\n','_'*20,'\n'
	time_three = datetime.now()
	print time_three.strftime('\n3. started at day:%j %H:%M:%S')
	#settings['Input Data Folder'] = sourceFolder
	#settings['Data Files'] = sourceFile
	# Call kriging function
	DEM = settings['DEM']
	dataFile = os.path.join(sourceFolder,sourceFile)
	krigedFile = kriging(dataFile,DEM)
	_,_,Dir,_ = namer(krigedFile)
	print 'Output to ',Dir,'\n'
	return krigedFile
#
#
def sliceDEM(settings,slice=20):
	'''4. SLICE DEM INTO ELEVATION BANDS
	Requires "settings" for DEM file and optional elevation slice. Default is 20m.'''
	#
	print sliceDEM.__doc__,'\n','_'*20,'\n'
	time_four = datetime.now()
	print time_four.strftime('\n4. started at day:%j %H:%M:%S')
	# Set the elevation interval for the DEM bands (m).
	#slice = 20
	demSlcNm,demSlcExt,demSlcPth,demSlcNmEx = namer(settings['DEM'])
	demSlcDir = os.path.join(demSlcPth,demSlcNm)
	if not os.path.exists(demSlcDir):
		os.makedirs(demSlcDir)
	demslicer(settings['DEM'],demSlcDir,slice)
	print 'Output to ',demSlcDir,'\n'
	return demSlcDir
#
#
def acc2bands(settings,krigedFile,demSlcDir):
	'''5. MULTIPLY THE MASS BALANCE LAYER BY THE ELEVATION BANDS
	Requires "settings" and krigedFile,demSlcDir.
	demSlcDir comes from sliceDEM and is a folder containing the elevation bands.'''
	#
	print acc2bands.__doc__,'\n','_'*20,'\n'
	time_five = datetime.now()
	print time_five.strftime('\n5. started at day:%j %H:%M:%S')
	flist = os.listdir(demSlcDir)
	name,ext,path,namefull = namer(krigedFile)
	outname = name+'_bands'
	outfolder = os.path.join(settings['Write Folder'],outname)
	if not os.path.exists(outfolder):
		os.makedirs(outfolder)
	rasXras(flist,demSlcDir,krigedFile,outfolder)
	print 'Output to ',outfolder,'\n'
	return 0
#
#
def readStakes(settings, densFunc):
	'''6. READ STAKE FILE AND CONVERT TO m w.e. USING DENSITY FUNCTION.
	DO NOT USE FOR DATA FROM THE OLDER PROTOCOLS (PRE 2013).
	Requires "settings" and densFunc from convertProbing'''
	#
	print readStakes.__doc__,'\n','_'*20,'\n'
	time_six = datetime.now()
	print time_six.strftime('\n6. started at day:%j %H:%M:%S')
	# Read csv file containing stake data
	File = csv.reader(open(os.path.join(settings['Input Data Folder'],settings['Stakes']),'rb'), delimiter=',', quotechar='|')
	# Read first line as header into a list
	colNames = File.next()
	# Ask for column to be converted.
	print 'Column headers in input file.\n',colNames
	# Snow depth, winter
	Dw = ''
	while Dw not in colNames:
		Dw = raw_input('Which column contains winter snow depth (Dw) readings? ')
	DwLoc = colNames.index(Dw)
	# Stake height, winter
	Hw = ''
	while Hw not in colNames:
		Hw = raw_input('Which column contains winter stake height (Hw) readings? ')
	HwLoc = colNames.index(Hw)
	# Stake height, summer
	Hs = ''
	while Hs not in colNames:
		Hs = raw_input('Which column contains summer stake height (Hs) readings? ')
	HsLoc = colNames.index(Hs)
	# Open data file and append mass balance headers
	outName,_,_,_ = namer(settings['Stakes'])
	outDataNamestr = os.path.basename(outName + '_e.csv')
	csvOut = open(os.path.join(settings['Write Folder'],outDataNamestr),'wb')
	colNames.append('Bw')
	colNames.append('Bs')
	colNames.append('Bn')
	Header = ','.join(colNames) + '\n'
	csvOut.write(Header)
	# Reopen file for appending calculated values
	csvOut = csv.writer(open(os.path.join(settings['Write Folder'],outDataNamestr),'ab'), delimiter = ',', lineterminator='\n')
	for inrow in File:
		# copy row
		outrow = inrow
		# Calculate winter balance from snow depth and density function
		Bw = densFunc(float(outrow[DwLoc]))
		# Get distance to winter's ice surface
		Ds = (float(outrow[HwLoc]) + float(outrow[DwLoc])) - float(outrow[HsLoc])
		# Abl = Acc - Ds * density. Density depends on sign of Ds
		if Ds >= 0 :
			rho = 0.55
		else:
			rho = 0.9
		Bs = Bw - (Ds * rho)
		Bn = Bw - Bs
		outrow.append(str("%.2f" % Bw))
		outrow.append(str("%.2f" % Bs))
		outrow.append(str("%.2f" % Bn))
		csvOut.writerow(outrow)
	# Weird thing of opening again just to close it but seems to be necessary
	csvOut = open(os.path.join(settings['Write Folder'],outDataNamestr),'ab')
	csvOut.close()
	# Pass file names on to next stage
	print 'Output to ',settings['Write Folder'],'\n'
	sourceFolder = settings['Write Folder']
	sourceFile = outDataNamestr
	return sourceFolder, sourceFile
#
#
def krigAbl(settings, sourceFolder, sourceFile):
	'''7. KRIG ABLATION FILE ACROSS DEM
	Requires "settings" and sourceFolder, sourceFile'''
	#
	print krigAcc.__doc__,'\n','_'*20,'\n'
	time_seven = datetime.now()
	print time_seven.strftime('\n7. started at day:%j %H:%M:%S')
	#settings['Input Data Folder'] = sourceFolder
	#settings['Data Files'] = sourceFile
	# Call kriging function
	DEM = settings['DEM']
	dataFile = os.path.join(sourceFolder,sourceFile)
	krigedFile = kriging(dataFile,DEM)
	_,_,Dir,_ = namer(krigedFile)
	print 'Output to ',Dir,'\n'
	return krigedFile
#
#
def abl2bands(settings,krigedFile,demSlcDir):
	'''8. MULTIPLY THE MASS BALANCE LAYER BY THE ELEVATION BANDS
	Requires "settings" and krigedFile,demSlcDir.
	demSlcDir comes from sliceDEM and is a folder containing the elevation bands.'''
	#
	print acc2bands.__doc__,'\n','_'*20,'\n'
	time_eight = datetime.now()
	print time_eight.strftime('\n8. started at day:%j %H:%M:%S')
	flist = os.listdir(demSlcDir)
	name,ext,path,namefull = namer(krigedFile)
	outname = name+'_bands'
	outfolder = os.path.join(settings['Write Folder'],outname)
	if not os.path.exists(outfolder):
		os.makedirs(outfolder)
	rasXras(flist,demSlcDir,krigedFile,outfolder)
	print 'Output to ',outfolder,'\n'
	return 0

#
def main():
	# READ SETTINGS FILE
	settings, envFile = envset()
	#
	# Ask whether to continue
	print readDensity.__doc__,' is next.\n'
	contAns = cont()
	#
	#
	if contAns == 0:
		# READ DENSITY DATA
		densityVectors, densityFiles = readDensity(settings)
	# Ask whether to continue
	print convertProbing.__doc__,'\n is next.\n'
	contAns = cont()
	#
	#
	if contAns == 0:
		# CONVERT PROBING DEPTHS TO M W.E. (ACCUMULATION)
		sourceFolderAcc, sourceFileAcc, densFunc = convertProbing(settings, densityVectors, densityFiles)
	# Ask whether to continue
	print krigAcc.__doc__,'\n is next.\n'
	contAns = cont()
	#
	#
	if contAns == 0:
		# KRIG ACCUMULATION
		accumulationK = krigAcc(settings,sourceFolderAcc, sourceFileAcc)
	# Ask whether to continue
	print sliceDEM.__doc__,'\n is next.\n'
	contAns = cont()
	#
	#
	if contAns == 0:
		# SLICE THE DEM INTO ELEVATION BANDS
		demSlcDir = sliceDEM(settings,20)
	# Ask whether to continue
	print acc2bands.__doc__,'\n is next.\n'
	contAns = cont()
	#
	#
	if contAns == 0:
		# SPLIT THE ACCUMULATION INTO ELEVATION BANDS
		acc2bands(settings,accumulationK,demSlcDir)
	# Ask whether to continue
	print readStakes.__doc__,'\n is next.\n'
	contAns = cont()
	#
	#
	if contAns == 0:
		# READ ABLATION STAKES AND CONVERT TO m w.e.
		sourceFolderAbl, sourceFileAbl = readStakes(settings, densFunc)
	# Ask whether to continue
	print krigAbl.__doc__,'\n is next.\n'
	contAns = cont()
	#
	#
	if contAns == 0:
		# KRIG OR INTERPOLATE ABLATION
		ablationK = krigAbl(settings, sourceFolderAbl, sourceFileAbl)
	# Ask whether to continue
	print abl2bands.__doc__,'\n is next.\n'
	contAns = cont()
	#
	#
	if contAns == 0:
		# SPLIT ABLATION INTO ELEVATION BANDS
		abl2bands(settings,ablationK,demSlcDir)
	# CALCULATE NET BALANCE
	# PRINT RESULTS
	return 0
#
if __name__ == "__main__":
	main() # Calls first function, named "main" which contains main body of programme.

