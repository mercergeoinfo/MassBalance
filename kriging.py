#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
#from __future__ import division
"""Kriging"""
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
# import csv
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
import numpy.ma as ma
import scipy.linalg
import scipy.stats as stats
from scipy.stats.stats import nanmean
# from scipy.stats import norm
# from pandas import Series, DataFrame
# import pandas as pd
#
## GIS SUPPORT
import gdal
import gdalconst
# from geostatsmodels import utilities, kriging, variograms, model, geoplot
# from qgis.core import QgsProject
# from PyQt4.QtCore import QFileInfo
# To run in QGIS enter following in QGIS console:
# execfile(u'/path/to/script.py'.encode('utf-8'))
#
## PLOT TOOLS
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as dates
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
# import pylab
# pylab.ion()
# pylab.ioff()
#
from standard import *
from spatialfiles import *
#
def kriging(*args,**kwargs):
	'''Main kriging starter function
	If only two arguments passed:  outfile = kriging(Datafile,DEM)
	Otherwise use named arguments: Datafile = *.csv, DEM = *.tif,
	If no DEM is given provide these named arguments for definition of grid:
	llx_ = flt|int, ll_y = flt|int, xres = flt|int, yres = flt|int, gcols = int, grows = int
	The following are optional:
	chc = int, modelSel = int, lagrc = flt|int, lagnrrc = int, sillrc = flt, nuggrc = flt, rangrc = flt, alpha = int
	chc: data column choice
	modelSel = variogram model number, lagrc: lag distance, lagnrrc = no. of lags, sillrc: sill, nuggrc: nugget,
	rangrc: range, alpha: exponent for Gaussian model
	'''
	if len(kwargs)>1:
		try:
			print kwargs.keys()
			if ('Datafile' in kwargs):
				Datafile = kwargs['Datafile']
			if ('DEM' in kwargs):
				DEM = kwargs['DEM']
			# If no DEM provided
			if ('ll_x' in kwargs):
				ll_x = kwargs['ll_x']
			if ('ll_y' in kwargs):
				ll_y = kwargs['ll_y']
			if ('xres' in kwargs):
				xres = kwargs['xres']
			if ('yres' in kwargs):
				yres = kwargs['yres']
			if ('gcols' in kwargs):
				gcols = kwargs['gcols']
			if ('grows' in kwargs):
				grows = kwargs['grows']
			# Optional arguments
			if ('chc' in kwargs):
				chc = kwargs['chc']
			if ('modelSel' in kwargs):
				modelSel = kwargs['modelSel']
			if ('lagrc' in kwargs):
				lagrc = kwargs['lagrc']
			if ('lagnrrc' in kwargs):
				lagnrrc = kwargs['lagnrrc']
			if ('sillrc' in kwargs):
				sillrc = kwargs['sillrc']
			if ('nuggrc' in kwargs):
				nuggrc = kwargs['nuggrc']
			if ('rangrc' in kwargs):
				rangrc = kwargs['rangrc']
			if ('alpha' in kwargs):
				alpha = kwargs['alpha']
			else:
				alpha = 2.0
			print "Run kriging using named arguments"
		except:
			print "Named arguments failed, this may cause a problem"
	else:
		try:
			Datafile = args[0]
			DEM = args[1]
			alpha = 2.0
			print "Run kriging using implied arguments %s\nand %s" %(Datafile, DEM)
		except:
			print "No named arguments and no implied arguments"
			sys.exit()
	#
	# Get DEM
	try:
		demdata,Xg,Yg,Xg1,Yg1,rx,ry,demmeta = DemImport(DEM)
		print "Start kriging with data from %s\nand DEM file %s" %(Datafile,DEM)
	except:
		try:
			xstart = int(round(ll_x))
			xstop = xstart + int(gcols)*int(round(xres))
			ystart = int(round(ll_y))
			ystop = ystart + int(grows)*int(round(yres))
			demRx = range(xstart,xstop,int(round(xres)))
			demRy = range(ystart,ystop,int(round(yres)))
			Xg1,Yg1 = np.meshgrid(demRx,demRy)
			# Convert grid to vectors
			Yg=Yg1.reshape((-1,1))
			Xg=Xg1.reshape((-1,1))
			#
			rx = len(demRx)
			ry = len(demRy)
			print "Start kriging with data from %s\nand given grid definition." %(Datafile)
		except:
			return 1
	# Get data to krige
	try:
		Xin,Yin,Zin = InDataArray(Datafile,chc)
	except:
		Xin,Yin,Zin = InDataArray(Datafile)
	#
	# Set output directory
	namstr,ext_,outDir,full_ = namer(Datafile)
	outDir = os.path.join(outDir,namstr)
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	#
	# Remove NaN and -inf values from data
	checklist =[]
	for i in range(len(Zin)):
		if np.isnan(Zin[i]) | np.isinf(Zin[i]):
			checklist.append(i)
	Znan = np.delete(Zin, checklist)
	X = np.delete(Xin, checklist)
	Y = np.delete(Yin, checklist)
	#
	#
	## _kriging_logtrans__
	# Should data be log transformed?
	Zlog = logTranHist(Znan,outDir,namstr)
	print "*"*10
	print "Log transform of data\n"
	print '\n**********CHECK IN**********\n', outDir, '\n*****FOR DATA HISTOGRAMS*****\n'
	yn = [['y','Y','yes','Yes','YES'],['n','N','no','No','NO']]
	trana = ''
	while (trana not in yn[0])and(trana not in yn[1]):
		print trana
		trana = raw_input('Log transform Z value?: ')
		if trana in yn[0]:
			Z = Zlog
		elif trana in yn[1]:
			Z = Znan
	#
	## __kriging_report1__
	# Report data summary to to user interface
	stat1 = 'Data file: '+full_+'\n'
	stat2 = 'Data ranges\n'
	stat3 = 'X: '+str(np.min(X))+' '+str(np.max(X))+'\n'
	stat4 = 'Y: '+str(np.min(Y))+' '+str(np.max(Y))+'\n'
	stat5 = 'Z: '+str(np.min(Z))+' '+str(np.max(Z))+'\n'
	stat5a = 'Zin: '+str(np.min(Zin))+' '+str(np.max(Zin))+'\n'
	stat5b = 'Zlog: '+str(np.min(Zlog))+' '+str(np.max(Zlog))+'\n'
	print stat1,stat2,stat3,stat4
	print stat5,stat5a,stat5b
	print '*'*10
	#
	## __kriging_meshgrids__
	# Create grids of raw data
	X1,X2 = np.meshgrid(X,X)
	Y1,Y2 = np.meshgrid(Y,Y)
	Z1,Z2 = np.meshgrid(Z,Z)
	#
	## __kriging_Dcalc__
	# Calculate distances in XY plane (pythagoras)
	D = np.sqrt((X1 - X2)**2 + (Y1 -Y2)**2)
	#Dang = np.arctan2((X1 - X2),(Y1 -Y2))*180/np.pi
	#Dang[Dang<0]+=360
	## __kriging_Gcalc__
	# Calculate experimental variogram
	G = 0.5*(Z1 - Z2)**2
	#
	modelGood = ""
	while modelGood not in yn[0]:
		# Calculate lags and maximum lag
		try:
			lag, max_lags = lagCheck(D, lagrc, lagnrrc)
		except:
			lag, max_lags = lagCheck(D)
		#
		# Create experimental variogram
		print "D:\n", np.shape(D)
		print "G:\n", np.shape(G)
		DE,GE,GEest,nuggest,sillest,rangeest,GErsqrd = expvar(D,G,lag,max_lags)
		#
		# If values provided as optional arguments reset here but give option to revert manually below.
		try:
			print 'Sill overridden from %f to %f'%(sillest,sillrc)
			sillest = sillrc
		except:
			pass
		try:
			print 'Nugget overridden from %f to %f'%(nuggest,nuggrc)
			nuggest = nuggrc
		except:
			pass
		try:
			print 'Range overridden from %f to %f'%(rangeest,rangrc)
			rangeest = rangrc
		except:
			pass
		#
		## __kriging_varEstAppr__
		# Plot variogram estimates for appraisal
		nugget = float(nuggest)
		sill = float(sillest)
		range_ = float(rangeest)
		# Create variogram models
		outDirEst = outDir+'/initVar'
		if not os.path.exists(outDirEst):
			os.makedirs(outDirEst)
		G_mod, modelTypes = modelPlotter(nugget,sill,D,range_,lag,max_lags,DE,GE,GEest,GErsqrd,namstr,outDirEst,alpha)
		print 'VARIOGRAM ESTIMATES AND MODEL PLOTS SAVED TO:\n',outDir,'\nPLEASE REVIEW BEFORE PROCEEDING.'
		#
		## __kriging_varModSetup__
		# Ask user to set sill, nugget and lag
		contAns = ''
		while contAns not in yn[0]:
			sillan = []
			while not isinstance(sillan, float):
				print 'Sill estimate: ',sillest
				print 'Sill set at: ',sill
				sillin = raw_input('\nEnter value for sill: ')
				try:
					sillan = float(sillin)
					sill = sillan
					print 'Sill reset at: ',sill
				except:
					if sillan == []:
						sillan = sill
						print 'Sill set at: ',sill
					else:
						sillan = []
			nuggan = []
			while not isinstance(nuggan, float):
				print 'Nugget estimate: ',nuggest
				print 'Nugget set at: ',nugget
				nuggin = raw_input('\nEnter value for nugget: ')
				try:
					nuggan = float(nuggin)
					nugget = nuggan
					print 'Nugget reset at: ',nugget
				except:
					if nuggan == []:
						nuggan = nugget
						print 'Nugget set at: ',nugget
					else:
						nuggan = []
			rangan = []
			while not isinstance(rangan,float):
				print 'Range estimate: ',rangeest
				print 'Range set at: ',range_
				rangin = raw_input('\nEnter value for range: ')
				try:
					rangan = float(rangin)
					range_ = rangan
					print 'Range reset at: ',range_
				except:
					if rangan == []:
						rangan = range_
						print 'Range set at : ',range_
					else:
						rangan = []
			for i in range(len(modelTypes)):
				print i, ': ',modelTypes[i]
			ansa = []
			while ansa not in range(len(modelTypes)):
				ansa = int(raw_input('Which model to use? '))
			print ansa, ' ',modelTypes[ansa]
			modelnr = int(ansa)
			G_mod, modSel = modelSelect(modelnr,nugget,sill,D,range_,alpha)
			G_mod_export,G_mod_export_rsqrd,G_mod_poly1d,G_mod_exp = VarModPrep(D,lag,max_lags,DE,GE,G_mod)
			varestplt(DE,GE,GEest,GErsqrd,nugget,sill,range_,namstr,outDir,G_mod_export,G_mod_export_rsqrd,modSel)
			print 'Variogram model plots saved to:\n',outDir,'\nPLEASE REVIEW BEFORE PROCEEDING.\n\n'
			contAns = raw_input('Enter "yes" to accept this model or any other key to repeat choices: ')
		modelGood = raw_input("Are you sure you want to continue with this model? ")
	# Krige data
	Zkg,G_mod_export,G_mod_export_rsqrd,model = krig(nugget,sill,range_,D,Z,X,Y,Xg,Yg,rx,ry,lag,max_lags,DE,GE,modelnr,alpha)
	print np.nanmin(Zkg), np.nanmean(Zkg), np.nanmax(Zkg)
	# Reverse log transform
	if trana  in yn[0]:
		#for x in np.nditer(Zkg, op_flags=['readwrite']):
			#x[...] = 10 ** x
		Zkg = np.exp(Zkg)
		print np.min(Zkg), np.mean(Zkg), np.max(Zkg)
		Z = np.exp(Z)
	#
	try:
		Z_ = Zkg.reshape(ry,rx)
		#SK = s2_k.reshape(ry,rx)
		DEMmaskedZ_ = np.ma.masked_where(demdata == demmeta[6], Z_)
		plotson = 'y'
	except:
		plotson = 'n'
# Write output to files
#	  outDir = os.path.join(outDir,namstr)
#	  if not os.path.exists(outDir):
#		  os.makedir(outDir)
	name = namstr + '_kriged_data'
	print "Data kriged: ",outDir,name
	if plotson == 'y':
		try:
			outfilename = datawrite(Z_,demdata,demmeta,name,outDir)
		except:
			print "Kriged data could not be written to file"
			sys.exit()
	else:
		print plotson
		outfilename = 'Not kriged'
	# Write summary to file
	txtName = namstr + '_vario.csv'
	outtxtName = os.path.join(outDir,txtName)
	varioFile = open(outtxtName, "ab")
	lagout = 'Lag: '+str(lag)+' nr. of lags: '+str(max_lags)+'\n'
	DEout = 'DE: '+str(DE)+'\n'
	GEout = 'GE: '+str(GE)+'\n'
	GEestout = 'GEest: '+str(GEest)+'\n'
	Nugout = 'Nugget, estimate: '+str(nugget)+' '+str(nuggest)+'\n'
	Sillout = 'Sill, estimate: '+str(sill)+' '+str(sillest)+'\n'
	Rangout = 'Range, estimate: '+str(range_)+' '+str(rangeest)+'\n'
	Alphout = 'Alpha set to:'+str(alpha)
	GErout = 'GE R^2: '+str(GErsqrd)+'\n'
	varioFile.write(stat1)
	varioFile.write(stat2)
	varioFile.write(stat3)
	varioFile.write(stat4)
	varioFile.write(stat5)
	varioFile.write(stat5a)
	varioFile.write(stat5b)
	varioFile.write(lagout)
	varioFile.write(DEout)
	varioFile.write(GEout)
	varioFile.write(GEestout)
	varioFile.write(Nugout)
	varioFile.write(Sillout)
	varioFile.write(Rangout)
	varioFile.write(Alphout)
	varioFile.write(GErout)
	varioFile.close()
	## PLOT RESULTS
	DEMplotReq = 'n'
	if DEMplotReq == 'y':
		DEMmasked = np.ma.masked_where(demdata == demmeta[6], demdata)
		name = 'DEM'
		map3d(Xg1,Yg1,DEMmasked,plt.cm.RdBu,DEMmasked.min(),DEMmasked.max(),name,outDir)
	#
	## Plot figure of data to krige
	name = namstr + '_in'
	maplot(X,Y,Z,'r','o',name,outDir)
	## Plot variogram estimators
	name = namstr + '_data'
	varestplt(DE,GE,GEest,GErsqrd,nugget,sill,range_,name,outDir,G_mod_export,G_mod_export_rsqrd,model)
	## Plot semi-variograms
	name = namstr + '_data'
	semvar(Z,D,G,name,outDir)
	if plotson == 'y':
	## Plot kriged data
		name = namstr + '_data'
		krigplot(Xg1,Yg1,X,Y,Z,DEMmaskedZ_,name,outDir)
	sys.stdout.flush()
	print 'Data column kriged: ', ZcolName
	print '\a'
	return outfilename
#
def lagCheck(D, **kwargs):
	'''Calculate number and size of lags and check with user'''
	if ('lagrc' in kwargs):
		lagrc = kwargs['lagrc']
		print 'lag size: ',lagrc
	if ('lagnrrc' in kwargs):
		lagnrrc = kwargs['lagnrrc']
		print 'lag no.: ',lagnrrc
	#
	lag, max_lags = lagcalc(D)
	#
	try:
		print 'Lag overridden from %f to %f'%(lag,lagrc)
		lag = lagrc
	except:
		pass
	try:
		print 'Number of lags overridden from %f to %f'%(max_lags,lagnrrc)
		max_lags = lagnrrc
	except:
		pass
	## __kriging_lagcheck__
	lagan = []
	while not lagan:
		print 'Lag size set at: ',lag
		lagin = raw_input('\nEnter value for lag: ')
		try:
			lagan = float(lagin)
			lag = lagan
			print 'Lag set at: ',lag
		except:
			if lagan == []:
				lagan = lag
				print 'Lag set at: ',lag
			else:
				lagan = []
	maxlagan = []
	while not maxlagan:
		print 'Lag number set at: ',max_lags
		maxlagin = raw_input('\nEnter number of lags: ')
		try:
			maxlagan = float(maxlagin)
			max_lags = maxlagan
			print 'Lag number set at: ',max_lags
		except:
			if maxlagan == []:
				maxlagan = max_lags
				print 'No of lags set to: ',max_lags
			else:
				maxlagan = []
	return lag, max_lags
#
#
def lagcalc(D1):
	'''Calculate lag (kriging)'''
	print '\nCalculating lag'
	# Replace zeros in diagonal with nan
	D2 = np.copy(D1)
	np.fill_diagonal(D2,np.nan)
	# nanmin needs to know which axis to check along, otherwise gives min of whole
	lag_c =np.mean(np.nanmin(D2,axis=1))
	# Get maximum distance and divide by 2
	hmd_c = np.nanmax(D1)/2
	# Set number of lags equal to 1/2 max distance / lag length
	max_lags_c = np.floor(hmd_c/lag_c)
	print lag_c, ' ',max_lags_c,'\n'
	return lag_c, max_lags_c
#
#
def expvar(D,G,lag,max_lags):
	'''Calculate variogram for kriging'''
	## Called by: kriging
	print '\nCalculating variogram'
	# Set seperation distances and calculate variogram
	LAGS = np.ceil(D/lag)
	DE=[]
	PN=[]
	GE=[]
	for i in range(1, int(max_lags)):
		SEL = (LAGS == i); #Selection matrix
		DE.insert(i, np.mean(np.mean(D[SEL]))) #Mean lag
		PN.insert(i, sum(sum(SEL == 1))/2) #Number of pairs
		GE.insert(i, np.mean(np.mean(G[SEL]))) #Variogram estimator
		#print 'DE: %f, PN: %f, GE: %f = '% (np.mean(np.mean(D[SEL])), sum(sum(SEL == 1))/2, np.mean(np.mean(G[SEL])))
	# polyfit to estimate
	GEpoly = np.polyfit(DE,GE,3,full=True)
	# Create function from polyfitted line
	GEpoly1d = np.poly1d(GEpoly[0])
	GEest = GEpoly1d(DE)
	nuggest = GEpoly1d(0)
	GEpoly1d1 = GEpoly1d.deriv(1)
	GEpoly1d2 = GEpoly1d.deriv(2)
	# Get roots of first derivative (inflection points)
	GEinfl = GEpoly1d1.r
	sillest = GEpoly1d(GEinfl.real[0])
	rangeest = GEinfl.real[0]
	rangeest = abs(rangeest)
	# Calculate R^2 for polyfit
	GErsqrd = rSquared(DE,GE,GEest)
	print 'Nugget: ',nuggest,' sill: ',sillest,' range: ',rangeest,'\n'
	return DE,GE,GEest,nuggest,sillest,rangeest,GErsqrd
#
#
def rSquared(x,y,f):
	'''Calculate coefficient of determination, R^2'''
	x = np.array(x)
	y = np.array(y)
	f = np.array(f)
	ybar = np.mean(y)
	ydif = y - ybar
	ydifsq = ydif * ydif
	SStot = np.sum(ydifsq)
	fydif = f - y
	fydifsq = fydif * fydif
	SSres = np.sum(fydifsq)
	Rsqrd = 1 - (SSres/SStot)
#	 print 'R^2 function:\n'
#	  print 'X: ',x
#	  print '\nY: ',y
#	  print '\nf: ',f
#	  print 'Mean y: ',ybar,'\n'
#	  print 'y-difs: ',ydif,'\n',ydifsq,'\n'
#	  print 'SStot: ',SStot,'\n'
#	  print 'f-difs: ',fydif,'\n',fydifsq,'\n'
#	  print 'SSres: ',SSres,'\n'
	print 'R^2: ',Rsqrd,'\n'
	return Rsqrd
#
#
def krig(nuggest,sillest,rangeest,D,Z,X,Y,Xg,Yg,rx,ry,lag,max_lags,DE,GE,modelnr,alpha=2.0):
	'''Kriging function. Kriging routine based on Trauth "Matlab Recipes for Earth Sciences", 2nd edition'''
	## Called by: kriging
	#####
	nugget = float(nuggest)
	sill = float(sillest)
	range_ = float(rangeest)
	#
	# Create variogram model
	G_mod, modSel = modelSelect(modelnr,nugget,sill,D,range_,alpha,lag)
	print 'Use ',modSel
	#
	# Prepare variogram for printing
	G_mod_export,G_mod_export_rsqrd,G_mod_poly1d,G_mod_exp = VarModPrep(D,lag,max_lags,DE,GE,G_mod)
	#####
	# Create empty matrices to accpet data (krigged estimates and variance)
	Zg = np.empty(Xg.shape)
	#s2_k = np.empty(Xg.shape)
	#
	# Manipulate model for application to data
	n = len(X)
	naddr = np.ones((1,n))
	naddc = np.ones((n+1,1))
	G_mod=np.vstack((G_mod,naddr))
	G_mod=np.hstack((G_mod,naddc))
	G_mod[n,n] = 0
	try:
		G_inv = np.linalg.inv(G_mod)
	except:
		return 'bad model',G_mod_export,G_mod_export_rsqrd, modSel
	## Krig values at grid points
	kto = int(len(Xg))
	for k in range(0, kto):
		DOR = ((X - Xg[k])**2 + (Y - Yg[k])**2)**0.5
		G_R,_ = modelSelect(modelnr,nugget,sill,DOR,range_,alpha,lag)
		G_R = np.append(G_R,1)
		# Lagrange multiplier
		E = np.dot(G_inv,G_R)
		Zg[k] = sum(E[0:n]*Z)
		#s2_k[k] = sum(E[0:n]*G_R[0:n])+E[n]
	Zg.reshape(ry,rx)
	#for i in Zg: print i
	#Z_ = Zg.reshape(ry,rx)
	#SK = s2_k.reshape(ry,rx)
	return Zg,G_mod_export,G_mod_export_rsqrd, modSel
#
def modelPlotter(nugget,sill,D,range_,lag,max_lags,DE,GE,GEest,GErsqrd,namstr,outDir,alpha=2.0):
	'''Call variogram model creators and plot results to file.'''
	## Called by: kriging
	modelTypes = ['TRS','Linear','Exponential with Nugget','Spherical','Gaussian or stable']
	for modSel in modelTypes:
		if modSel == 'Exponential with Nugget':
			print modSel
			G_mod = varMod_expNug(nugget,sill,D,range_)
		elif modSel == 'Spherical':
			print modSel
			G_mod = varMod_spheric(nugget,sill,D,range_)
		elif modSel == 'Gaussian or stable':
			print modSel
			G_mod = varMod_Gauss(nugget,sill,D,range_,alpha)
		elif modSel == 'Linear':
			print modSel
			G_mod = varMod_Linear(nugget,sill,D,range_)
		elif modSel == 'TRS':
			print modSel
			G_mod = varMod_TRS(nugget,sill,D,range_)
		G_mod_export,G_mod_export_rsqrd,G_mod_poly1d,G_mod_exp = VarModPrep(D,lag,max_lags,DE,GE,G_mod)
		varestplt(DE,GE,GEest,GErsqrd,nugget,sill,range_,namstr,outDir,G_mod_export,G_mod_export_rsqrd,modSel)
	return G_mod, modelTypes
#
#
def modelSelect(modelnr,nugget,sill,D,range_,lag,alpha=2.0):
	modelTypes = ['TRS','Linear','Exponential with Nugget','Spherical','Gaussian or stable']
	try:
		modSel = modelTypes[modelnr]
	except:
		modSel = 'Exponential with Nugget'
	if modSel == 'Exponential with Nugget':
		G_mod = varMod_expNug(nugget,sill,D,range_)
	elif modSel == 'Spherical':
		G_mod = varMod_spheric(nugget,sill,D,range_)
	elif modSel == 'Gaussian or stable':
		G_mod = varMod_Gauss(nugget,sill,D,range_,alpha)
	elif modSel == 'Linear':
		G_mod = varMod_Linear(nugget,sill,D,range_)
	elif modSel == 'TRS':
		G_mod = varMod_TRS(nugget,sill,D,range_)
	return G_mod, modSel
#
#
def VarModPrep(D,lag,max_lags,DE,GE,G_mod):
	'''Create variogram model for printing'''
	## Called by: modelPlotter, kriging
	# Export variogram for printing
	LAGS = np.ceil(D/lag)
	G_mod_copy = G_mod[:]
	G_mod_exp = []
	for i in range(1, int(max_lags)):
		SEL = (LAGS == i)
		G_mod_exp.insert(i, np.mean(np.mean(G_mod_copy[SEL])))
	G_mod_poly = np.polyfit(DE,G_mod_exp,3,full=True)
	# Create function from polyfitted line
	G_mod_poly1d = np.poly1d(G_mod_poly[0])
	G_mod_export = G_mod_poly1d(DE)
	# Calculate R^2 for polyfit
	G_mod_export_rsqrd = rSquared(DE,GE,G_mod_export)
	return G_mod_export,G_mod_export_rsqrd,G_mod_poly1d,G_mod_exp
#
#
def varMod_expNug(nugget,sill,D,range_):
	'''Variogram model from Trauth book.
	Variogram: exponential with nugget variance'''
	## Called by: modelPlotter, modelSelect
	#print '\nExponential with nugget\n'
	sill=sill-nugget
	G_mod = (nugget + sill*(1 - np.exp(-3*D/range_)))*(D>0)
	return G_mod
#
#
def varMod_spheric(nugget,sill,D,range_):
	'''Variogram: spherical'''
	## Called by: modelPlotter, modelSelect
	#print '\nSpherical\n'
	sill=sill-nugget
	G_mod_1 = (nugget + sill*(((3*D)/(2*range_)) - 0.5*((D/range_)**3)))*(D>0)*(D<range_)
	G_mod_2 = (nugget + sill) *(D>0)*(D>=range_)
	G_mod = G_mod_1 + G_mod_2
	return G_mod
#
#
def varMod_Gauss(nugget,sill,D,range_,alpha = 2.0):
	'''Variogram: Gaussian or stable. If alpha not passed it is set to 2 and the model is
	Gaussian.'''
	## Called by: modelPlotter, modelSelect
	#print '\nGaussian or stable model\n'
	#print 'Alpha = ',alpha
	sill=sill-nugget
	#G_mod = (nugget + sill*(1 - np.exp(-1*(D**alpha)/(((range_**2)/3)**alpha))))*(D>0)
	#r1 = (range_**2)/3
	#G_mod = (nugget + sill*(1 - np.exp(-1*D**alpha/r1**alpha)))*(D>0)
	G_mod = (nugget + sill*(1 - np.exp(-3*D**alpha/range_**alpha)))*(D>0)
	return G_mod
#
#
def varMod_Linear(nugget,sill,D,range_):
	'''Bounded linear model as used by TRS'''
	## Called by: modelPlotter, modelSelect
	#print 'Bounded linear model as used by TRS'
	#range_ = 500
	#sill = 1
	#sill=sill-nugget
	G_mod_1 = sill*(D/range_)*(D>0)*(D<range_)
	G_mod_2 = sill *(G_mod_1 == 0)
	G_mod = G_mod_1 + G_mod_2
	return G_mod
#
def varMod_TRS(nugget,sill,D,range_):
	'''Bounded linear model as used by TRS.
	According R. Petterssons kriging.m file in the massbalance toolbox for Matlab in lines 189 to 212
	linear model with no nugget, slope of 1 and power of 1. These are default values left unaltered.'''
	## Called by: modelPlotter, modelSelect
	#print 'Bounded linear model as used by TRS'
	range_ = range_
	slope = 1
	G_mod_1 = slope*(D/range_)*(D>0)*(D<range_)
	G_mod_2 = slope*(G_mod_1 == 0)
	G_mod = G_mod_1 + G_mod_2
	return G_mod
#
def polyfit2d(x, y, z, order=2):
	'''Fit a polynomial to a surface (for detrending kriging data)
	 Taken from http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent'''
	print '\nCreate polynomial trend surface'
	ncols = (order + 1)**2
	G = np.zeros((x.size, ncols))
	ij = itertools.product(range(order+1), range(order+1))
	for k, (i,j) in enumerate(ij):
		G[:,k] = x**i * y**j
	#m, mresids,mrank,ms = linalg.lstsq(G, z)
	m,_,_,_ = np.linalg.lstsq(G, z)
	return m
#
def polyval2d(x, y, m):
	'''Apply polynomial to 3D data (for detrending kriging data)
	Taken from http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent'''
	print '\nApply polynomial trend surface'
	order = int(np.sqrt(len(m))) - 1
	ij = itertools.product(range(order+1), range(order+1))
	z = np.zeros_like(x)
	for a, (i,j) in zip(m, ij):
		z += a * x**i * y**j
	return z
#
def detrend(X,Y,Z,Xg,Yg,Xg1):
		## Call function to create trend polynomial
		detrend_m = polyfit2d(np.array(X),np.array(Y),np.array(Z))
		## Call function to detrend data
		Z_trend = polyval2d(np.array(X), np.array(Y), detrend_m)
		Zdtr = Z -Z_trend
		Zdtr = Zdtr.tolist()
		#
		# Set up grid to "retrend" kriged detrended data
		Xll = Xg.tolist()
		Xl = []
		for i in Xll:
			for j in i:
				Xl.append(j)
		Xa=np.array(Xl)
		#
		Yll = Yg.tolist()
		Yl = []
		for i in Yll:
			for j in i:
				Yl.append(j)
		Ya=np.array(Yl)
		#
		Z_retrend_in = polyval2d(Xa, Ya, detrend_m)
		Z_retrend = Z_retrend_in.reshape(np.shape(Xg1))
		return Zdtr, Z_retrend
#
def fitFunc(xvec,yvec,fitOrder=1):
	"""Create least squares fitted function through data"""
	z = np.polyfit(np.array(xvec), np.array(yvec), fitOrder)
	p = np.poly1d(z)
	return p
#
def maplot(x,y,z,c,m,name,outDir):
	plt.figure()
	plt.autoscale(enable=True, axis='both', tight=False)
	plt.plot(x,y,m,color=c,hold=True)
	zrange = [np.min(z),np.max(z)]
	zmess = 'Z Range = '+ str(zrange[0]) + ' to ' + str(zrange[1])
	plt.xlabel(zmess)
	plt.title(name)
	indatplt = name + '_dataplot.png'
	plt.savefig(os.path.join(outDir,indatplt))
	plt.close()
	return 0
#
def map3d(x,y,z,c,valmin,valmax,name,outDir):
	'''Plot 3D map (kriging)'''
	plt.figure(figsize=(16,8),facecolor='w')
	h = plt.pcolor(x,y,z,cmap=c, vmax = valmax, vmin = valmin)
	plt.autoscale(enable=True, axis='both', tight=False)
	plt.savefig(os.path.join(outDir, (name + '_data3d.png')))
	plt.close()
	return 0
#
def logTranHist(Zin,outfolder,name):
	''' Plot histograms of data and log transformed data'''
	## Called by: kriging
	# Should data be log transformed?
	Zlog = np.log(Zin)
	if np.min(Zin) < 0:
		return Zlog
	checklist =[]
	for i in range(len(Zlog)):
		if np.isnan(Zlog[i]) | np.isinf(Zlog[i]):
			checklist.append(i)
	Zlog = np.delete(Zlog, checklist)
	plt.figure(1)
	plt.hist(Zin, bins=20, histtype='stepfilled', normed=True, color='b', label='Untransformed')
	plt.title("Untransformed Data")
	plt.xlabel("Value")
	plt.ylabel("Number")
	plt.legend()
	plt.savefig(os.path.join(outfolder,(name +'Hist.png')))
	#plt.show()
	plt.figure(2)
	plt.hist(Zlog, bins=20, histtype='stepfilled', normed=True, color='r', label='Log transformed')
	plt.title("Log Transform of data")
	plt.xlabel("Value")
	plt.ylabel("Number")
	plt.legend()
	#plt.show()
	plt.savefig(os.path.join(outfolder,(name +'logHist.png')))
	plt.close()
	return Zlog
#
def semvar(Z,D,G,name,outDir):
	'''Plot variogram (kriging)'''
	indx=[]
	for i in range(len(Z)):indx.append(i+1)
	C,R = np.meshgrid(indx,indx)
	I = (R > C)
	Dpl = D * I
	Gpl = G * I
	plt.figure()
	plt.plot(Dpl,Gpl,'.')
	plt.xlabel('Lag Distance')
	titletext = name + ' Semivariogram'
	plt.title(titletext)
	plt.ylabel('Variogram')
	plt.savefig(os.path.join(outDir,(name + '_semivar.png')))
	plt.close()
	return 0
#
def varestplt(DE,GE,GEest,GErsqrd,nuggest,sillest,rangeest,name,outDir,G_mod_export,G_mod_export_rsqrd,model):
	'''Plot variogram estimate (kriging)'''
	## Called by: modelPlotter
	plt.figure()
	plt.plot(DE,GE,'.',hold=True)
	b = [0, max(DE)]
	c = [sillest,sillest]
	plt.plot(b,c,'--r',hold=True)
	y1 = 1.1 * max(GE)
	plt.ylim(0,y1)
	plt.plot(DE,GEest,':b',hold=True)
	plt.plot(DE,G_mod_export,'--k',hold=True)
	plt.xlabel('Averaged distance between observations')
	plt.ylabel('Averaged semivariance')
	titletext = name + ' Variogram Estimator and ' + model + ' model'
	plt.title(titletext)
	nuggtext = 'Nugget estimate =\n %s' % (str(nuggest))
	silltext = 'Sill estimate =\n %s' % (str(sillest))
	rangetext = 'Range estimate = \n %s' % (str(rangeest))
	GErsqrdtext = 'Estimate R^2 =\n %s' % (str(GErsqrd))
	G_modrsqrdtext = 'Model R^2 =\n %s' % (str(G_mod_export_rsqrd))
	plt.text(max(DE)*0.8,max(GE)*0.1,nuggtext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.text(max(DE)*0.8,max(GE)*0.9,silltext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.text(max(DE)*0.8,max(GE)*0.3,rangetext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.text(max(DE)*0.8,max(GE)*0.5,GErsqrdtext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.text(max(DE)*0.8,max(GE)*0.7,G_modrsqrdtext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.savefig(os.path.join(outDir,(name + '_' + model+ '.png')))
	plt.close()
	return 0
#
def krigvarplot(Xg1,Yg1,SK,name,outDir):
	'''Plot variance (kriging)'''
	plt.figure(figsize=(10,5),facecolor='w')
	h = plt.pcolor(Xg1,Xg2,SK)
	plt.title(' Kriging Variance')
	plt.xlabel('x-Coordinates')
	plt.ylabel('y-Coordinates')
	plt.colorbar()
	plt.axis("equal")
	plt.plot(X,Y,'ok',hold=True)
	krigvarplt = name + '_krigvarplot.png'
	plt.savefig(os.path.join(outDir,krigvarplt))
	plt.close()
	return
#
def InDataArray(fileName):
	'''Simplified csv data reader that returns 3 vectors'''
	## Called by: kriging
	print '\n'*2,'*'*10
	print 'Read ',fileName
	InFile = open(fileName,'rb')
	Headers = InFile.next().strip().split(',')
	Xchc = ['Easting','SWEREF99_E','RT90_E','E','X','East','EASTING','EAST']
	Ychc = ['Northing','SWEREF99_N','RT90_N','N','Y','North','NORTHING','NORTH']
	Zchc = ['Elevation','ELEV','Z','mwe','Bw','Bs','Bn','Dw']
	print 'File Headers:\n', Headers
	Xcol = ''
	Ycol = ''
	Zcol = ''
	for a in Xchc:
		if a in Headers:
			print a
			a_ans = raw_input('Use as x? (y/any other key): ')
			if a_ans == 'y':
				Xcol = a
				break
	for b in Ychc:
		if b in Headers:
			print b
			b_ans = raw_input('Use as y? (y/any other key): ')
			if b_ans == 'y':
				Ycol = b
				break
	for c in Zchc:
		if c in Headers:
			print c
			c_ans = raw_input('Use as z? (y/any other key): ')
			if c_ans == 'y':
				Zcol = c
				break
	while Xcol not in Headers:
		Xcol = raw_input('Enter column for "x": ')
	while Ycol not in Headers:
		Ycol = raw_input('Enter column for "y": ')
	while Zcol not in Headers:
		Zcol = raw_input('Enter column for "z": ')
	#
	# Dirty little hack to get name of z for kriging return
	global ZcolName
	ZcolName = Zcol
	#
	inFile = open(fileName,'rb')
	X = []
	Y = []
	Z = []
	j = 1
	for line in csv.DictReader(inFile, delimiter=','):
		catch = 0
		#print j,': ',' x: ',line[Xcol].strip(),' y: ',line[Ycol].strip(),' z: ',line[Zcol].strip()
		try:
			float(line[Zcol].strip())
		except ValueError:
			catch = 1
			print 'Line ',j,': ',line, ' not processed'
			j = j+1
			continue
		#if not (line[Zcol].strip()) or isnan(float(line[c3].strip())) == True:
			#catch = 1
			#print 'Line ',j,': ',line, ' not processed'
		if catch == 0:
			X.append(float(line[Xcol].strip()))
			Y.append(float(line[Ycol].strip()))
			Z.append(float(line[Zcol].strip()))
		j = j+1
	return X,Y,Z
	#
#
def DemImport(demfile):
	'''Import a DEM file and return grid plus meta data. Written for kriging.
	To use: demdata,Xg,Yg,Xg1,Yg1,rx,ry,demmeta = DemImport(demfile)'''
	## Called by: kriging
	time_one = datetime.now()
	print '\nImporting ',demfile
	#print time_one.strftime('at day:%j %H:%M:%S')
	# register all of the GDAL drivers
	gdal.AllRegister()
	# Open file
	dem = gdal.Open(demfile)
	if dem is None:
		print 'Could not open ',demfile,'\n'
		sys.exit(1)
	# Get coordinate system parameters
	projec = dem.GetProjection()
	transf = dem.GetGeoTransform()
	ul_x = transf[0]
	ul_y = transf[3]
	xres = transf[1]
	yres = transf[5]
	# get image size
	demrows = dem.RasterYSize
	demcols = dem.RasterXSize
	# Calculate corners
	ll_x = ul_x
	ll_y = ul_y + (demrows * yres)
	#
	driveshrt = dem.GetDriver().ShortName
	driver = gdal.GetDriverByName(driveshrt)
	#metadata = driver.GetMetadata()
	# Read the dem band to a matrix called band_1
	demband = dem.GetRasterBand(1)
	# Access data in rastern band as array
	demdata = demband.ReadAsArray(0,0,demcols,demrows)
	# gdal interpolation creates "upside down files, hence this
	Yres=yres
	if Yres/np.fabs(Yres)==-1:
		Yres=-1*Yres
		demdata = np.flipud(demdata)
	#
	# get nodata value
	demnandat = demband.GetNoDataValue()
	# get minimum and maximum value
	demmindat = demband.GetMinimum()
	demmaxdat = demband.GetMaximum()
	if demmindat is None or demmaxdat is None:
			(demmindat,demmaxdat) = demband.ComputeRasterMinMax(1)
	#
	## Create grid to krig to.
	xstart = int(round(ll_x))
	xstop = xstart + int(demcols)*int(round(xres))
	ystart = int(round(ll_y))
	ystop = ystart + int(demrows)*int(round(Yres))
	demRx = range(xstart,xstop,int(round(xres)))
	demRy = range(ystart,ystop,int(round(Yres)))
	Xg1,Yg1 = np.meshgrid(demRx,demRy)
	# Convert grid to vectors
	Yg=Yg1.reshape((-1,1))
	Xg=Xg1.reshape((-1,1))
	#
	rx = len(demRx)
	ry = len(demRy)
	# Collect metadata to list for return
	demmeta = []
	demmeta.append(['projection','geotransform','driver','rows','columns','nanvalue','min','max'])
	demmeta.append(projec)
	demmeta.append(transf)
	demmeta.append(driver)
	demmeta.append(demrows)
	demmeta.append(demcols)
	demmeta.append(demnandat)
	demmeta.append(demmindat)
	demmeta.append(demmaxdat)
	#for i in demmeta: print i
	#print 'demmeta:\n',demmeta[0],'\n'
	return demdata,Xg,Yg,Xg1,Yg1,rx,ry,demmeta
#
def varestplt(DE,GE,GEest,GErsqrd,nuggest,sillest,rangeest,name,outDir,G_mod_export,G_mod_export_rsqrd,model):
	'''Plot variogram estimate (kriging)'''
	## Called by: modelPlotter
	plt.figure()
	plt.plot(DE,GE,'.',hold=True)
	b = [0, max(DE)]
	c = [sillest,sillest]
	plt.plot(b,c,'--r',hold=True)
	y1 = 1.1 * max(GE)
	plt.ylim(0,y1)
	plt.plot(DE,GEest,':b',hold=True)
	plt.plot(DE,G_mod_export,'--k',hold=True)
	plt.xlabel('Averaged distance between observations')
	plt.ylabel('Averaged semivariance')
	titletext = name + ' Variogram Estimator and ' + model + ' model'
	plt.title(titletext)
	nuggtext = 'Nugget estimate =\n %s' % (str(nuggest))
	silltext = 'Sill estimate =\n %s' % (str(sillest))
	rangetext = 'Range estimate = \n %s' % (str(rangeest))
	GErsqrdtext = 'Estimate R^2 =\n %s' % (str(GErsqrd))
	G_modrsqrdtext = 'Model R^2 =\n %s' % (str(G_mod_export_rsqrd))
	plt.text(max(DE)*0.8,max(GE)*0.1,nuggtext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.text(max(DE)*0.8,max(GE)*0.9,silltext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.text(max(DE)*0.8,max(GE)*0.3,rangetext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.text(max(DE)*0.8,max(GE)*0.5,GErsqrdtext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.text(max(DE)*0.8,max(GE)*0.7,G_modrsqrdtext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.savefig(os.path.join(outDir,(name + '_' + model+ '.png')))
	plt.close()
	return 0
#
def krigplot(Xg1,Yg1,X,Y,Z,DEMmdZ_,name,outDir):
	'''Plot kriged data as map'''
	# Get range of values for colour scale
	#print name,type(DEMmdZ_)
	#currentmin = DEMmdZ_.min()
	#currentmax = DEMmdZ_.max()
	currentmin = np.nanmin(DEMmdZ_[DEMmdZ_ != np.inf])
	currentmax = np.nanmax(DEMmdZ_[DEMmdZ_ != np.inf])
	print "Kriged values range from {0:.3f} to {1:.3f}\n".format(currentmin, currentmax)
	plt.figure(figsize=(16,8),facecolor='w')
	if currentmin > -1.5 and currentmax < 1.5:
		print "-1.5 to 1.5"
		h = plt.pcolor(Xg1,Yg1,DEMmdZ_,cmap=plt.cm.RdBu, vmax = 1.5, vmin = -1.5)
		plt.scatter(X, Y, c=Z, cmap=plt.cm.RdBu, vmax = 1.5, vmin = -1.5, hold=True)
	elif currentmin > -3.0 and currentmax < 3.0:
		print "-3.0 to 3.0"
		h = plt.pcolor(Xg1,Yg1,DEMmdZ_,cmap=plt.cm.RdBu, vmax = 3.0, vmin = -3.0)
		plt.scatter(X, Y, c=Z, cmap=plt.cm.RdBu, vmax = 3.0, vmin = -3.0, hold=True)
	elif currentmin > 0 and currentmax < 3.0:
		print "0 to 3.0"
		h = plt.pcolor(Xg1,Yg1,DEMmdZ_,cmap=plt.cm.RdBu, vmax = 3.0, vmin = 0)
		plt.scatter(X, Y, c=Z, cmap=plt.cm.RdBu, vmax = 3.0, vmin = 0, hold=True)
	elif currentmin > -5.0 and currentmax < 5.0:
		print "-5.0 to 5.0"
		h = plt.pcolor(Xg1,Yg1,DEMmdZ_,cmap=plt.cm.RdBu, vmax = 5.0, vmin = -5.0)
		plt.scatter(X, Y, c=Z, cmap=plt.cm.RdBu, vmax = 5.0, vmin = -5.0, hold=True)
	elif currentmax < 6:
		print "<6"
		h = plt.pcolor(Xg1,Yg1,DEMmdZ_,cmap=plt.cm.Reds, vmax = 5, vmin = 0)
		plt.scatter(X, Y, c=Z, cmap=plt.cm.RdBu, vmax = 5, vmin = 0, hold=True)
	else:
		print "Otherwise"
		h = plt.pcolor(Xg1,Yg1,DEMmdZ_,cmap=plt.cm.Blues, vmax = int(currentmax), vmin = int(currentmin))
		plt.scatter(X, Y, c=Z, cmap=plt.cm.RdBu, vmax = int(currentmax), vmin = int(currentmin), hold=True)
	titletext = name + ' Kriging Estimate'
	plt.title(titletext)
	plt.xlabel('x-Coordinates')
	plt.ylabel('y-Coordinates')
	plt.colorbar()
	plt.axis("equal")
	#plt.plot(X,Y,'.k',hold=True)
	#plt.scatter(X, Y, c=Z, hold=True)
	krigtext = 'Min = %f\nMax = %f' % (currentmin,currentmax)
	plt.text(min(X)+((max(X)-min(X))*0.9),min(Y)+((max(Y)-min(Y))*0.9),krigtext,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	krigplt = name + '_krigplot.png'
	plt.savefig(os.path.join(outDir,krigplt))
	plt.close()
	return
#
def krigvarplot(Xg1,Yg1,SK,name,outDir):
	'''Plot variance (kriging)'''
	plt.figure(figsize=(10,5),facecolor='w')
	h = plt.pcolor(Xg1,Xg2,SK)
	plt.title(' Kriging Variance')
	plt.xlabel('x-Coordinates')
	plt.ylabel('y-Coordinates')
	plt.colorbar()
	plt.axis("equal")
	plt.plot(X,Y,'ok',hold=True)
	krigvarplt = name + '_krigvarplot.png'
	plt.savefig(os.path.join(outDir,krigvarplt))
	plt.close()
	return
#
