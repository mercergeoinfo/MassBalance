#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
#from __future__ import division
"""Template file for Python 2.7"""
__author__ = 'Andrew Mercer'
__contact__ = 'mercergeoinfo@gmail.com'
__maintainer__ = 'Andrew Mercer'
__version__ = '1.0.2'
__date__ = '25/08/2015'
#
## IMPORTS
import sys
import os
# import glob
# import fnmatch
# import gc
import csv
import pickle
# import time
# from datetime import datetime, timedelta
# import subprocess
# import urllib2
# import os.path
# import zipfile
# import StringIO
#
## IMAGE MANIPULATION
# import cv2
# import PIL
# from PIL import Image
# from PIL.GifImagePlugin import getheader, getdata
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
## FUNCTIONS
#
def getClipboardData():
	'''Get contents of system clipboard
	To use: data = getClipboardData()'''
	## Called by:
	p = subprocess.Popen(['pbpaste'], stdout=subprocess.PIPE)
	retcode = p.wait()
	data = p.stdout.read()
	return data
#
def setClipboardData(data):
	'''Set contents of system clipboard
	To use: setClipboardData(data)'''
	## Called by:
	p = subprocess.Popen(['pbcopy'], stdin=subprocess.PIPE)
	p.stdin.write(data)
	p.stdin.close()
	retcode = p.wait()
	return 0
#
def namer(pathstring):
	'''Get name and extension of a file
	To use: name,ext,path,namefull = namer(pathstring) '''
	## Called by:
	namefull = os.path.basename(pathstring)
	path = os.path.dirname(pathstring)
	ext = namefull.split('.')[1]
	name = namefull.split('.')[0]
	return name,ext,path,namefull
#
def makeDir(path,name):
	'''Make a directory at given location
	To use: outDir = makeDir(path, name)'''
	## Called by:
	targetDir = os.path.join(path,name)
	if not os.path.exists(targetDir):
		os.makedirs(targetDir)
	return targetDir
#
def makeDirTS(path,name):
	'''Make a directory at given location with current time stamp added to name
	To use: outDir, containDir = makeDirTS(path, name)'''
	## Called by:
	tempDir = os.path.join(path,name)
	namstr = datetime.now().strftime('%j_%H%M%S')
	targetDir = os.path.join(tempDir,namstr)
	if not os.path.exists(targetDir):
		os.makedirs(targetDir)
	return targetDir, tempDir
#
def filelist(folder,ending):
	"""Return list of specific file type in folder other than current folder.
	To use: filelist = filelist(folder,ending)"""
	# Called by:
	matchstring = '*.' + ending
	filelist = fnmatch.filter(os.listdir(folder),matchstring)
	print 'From ', folder, ' the following matched ', matchstring, '\n', filelist
	return filelist
#
def picksave(data,fileName,location):
	'''Pickle data into a pickle database for easy retrieval into Python'''
	# Save data dictionary to pickle database
	fileName = fileName + '.p'
	pickleFile = os.path.join(location,fileName)
	with open(pickleFile, 'wb') as fp:
		pickle.dump(data, fp)
	return pickleFile
#
def pickload(pickleFile):
	'''Load data from a pickle database'''
	with open(pickleFile, 'rb') as fp:
					data = pickle.load(fp)
	return data
#
def import2vector(fileName, dateString = '%d/%m/%y %H:%M:%S'):
	'''Imports the data as vectors in a dictionary. dateString is optional and can be set to match datetime format
	To use: db = import2vector(filename) or db = import2vector(filename, dateString = '%d/%m/%y %H:%M:%S')'''
	# Called by:
	# Open file
	InFile = open(fileName,'rb')
	line = InFile.next()
	Name = os.path.basename(fileName).split('.')[0]
	# Get headers
	Headers = line.strip().split(',')
	# Create dictionary for data
	data = {}
	data['Description'] = {}
	data['Description']['Source'] = Name
	i=0
	# Set up list of Nan values
	nanlist = ['NAN','NaN','nan','NULL','null','-9999',-9999,'']
	# Read through data file
	for i in Headers:
		data[i] = []
	for row in InFile:
		if row != "\n":
			# Split read line of data into list
			dataIn = row.strip().split(',')
			for i in range(len(dataIn)):
				# Check for NaN and empty values
				if dataIn[i] in nanlist:
					 dataIn[i] = np.nan
				else:
					# Try date formatted data conversion
					try:
						dataIn[i] = datetime.strptime(dataIn[i],dateString)
					except:
						# Try converting to float
						try:
							dataIn[i] = float(dataIn[i])
						except:
							# Leave as string
							dataIn[i] = dataIn[i]
				# Add to vector
				data[Headers[i]].append(dataIn[i])
	for i in Headers:
		# Convert to numpy arrays
		data[i] = np.array(data[i])
		try:
			# Create posts containing basic statistics for each numerical column (vector)
			data['Description'][(str(i)+'_min')] = np.nanmin(data[i])
			data['Description'][(str(i)+'_max')] = np.nanmax(data[i])
			data['Description'][(str(i)+'_mean')] = np.nanmean(data[i])
			data['Description'][(str(i)+'_stdDev')] = np.nanstd(data[i])
			data['Description'][(str(i)+'_median')] = np.median(data[i])
		except:
			print "\nStatistics not computable for %s\n" % str(i)
	return data
#
def nancov(m):
	rows = float(np.shape(m)[0])
	print rows
	rowssqr =np.ones([rows,rows])
	mdivn = rowssqr/rows
	#mavg = mdivn.dot(m)
	mavg = matdot(mdivn,m)
	print 'Averages: ',mavg[0]
	mdiff = m-mavg
	mdiffT = mdiff.T
	#mdevscore = mdiffT.dot(mdiff)
	mdevscore = matdot(mdiffT,mdiff)
	print 'Deviation Score:\n',mdevscore
	mcov = mdevscore/(rows-1)
	print 'Covariance matrix:\n',mcov
	return mcov, mdevscore, mavg[0]
#
#
def fitFunc(xvec,yvec,fitOrder=1):
	"""Create least squares fitted function through data"""
	z = np.polyfit(np.array(xvec), np.array(yvec), fitOrder)
	p = np.poly1d(z)
	return p
