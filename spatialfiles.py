#!/Users/andrew/anaconda/bin/python
# -*- coding: utf-8 -*-
#from __future__ import division
"""Spatial File Functions"""
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
# import scipy.stats as stats
from scipy.stats.stats import nanmean
# from scipy.stats import norm
# from pandas import Series, DataFrame
# import pandas as pd
#
## GIS SUPPORT
import gdal
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
from standard import *
#
## FUNCTIONS
#
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
#
def rasterImport(file):
	'''Import a raster file and return grid plus meta data.
	To use: data, meta, metadata = rasterImport(file)'''
	time_one = datetime.now()
	print '\nImporting ',file
	#print time_one.strftime('at day:%j %H:%M:%S')
	# register all of the GDAL drivers
	gdal.AllRegister()
	# Open file
	raster = gdal.Open(file)
	if raster is None:
		print 'Could not open ',file,'\n'
		sys.exit(1)
	# Get coordinate system parameters
	projec = raster.GetProjection()
	srs=osr.SpatialReference(wkt=projec)
	transf = raster.GetGeoTransform()
	ul_x = transf[0]
	ul_y = transf[3]
	xres = transf[1]
	yres = transf[5]
	# get image size
	rows = raster.RasterYSize
	cols = raster.RasterXSize
	dims = {'xres':xres,'yres':yres,'rows':rows,'cols':cols}
	# Calculate corners
	ll_x = ul_x
	ll_y = ul_y + (rows * yres)
	ur_x = ul_x + (cols * xres)
	ur_y = ul_y
	lr_x = ur_x
	lr_y = ll_y
	corners = {'ll':(ll_x,ll_y),'ul':(ul_x,ul_y),'ur':(ur_x,ur_y),'lr':(lr_x,lr_y)}
	#
	driveshrt = raster.GetDriver().ShortName
	driver = gdal.GetDriverByName(driveshrt)
	metadata = driver.GetMetadata()
	metakeys = metadata.keys()
	# Read the file band to a matrix called band_1
	band = raster.GetRasterBand(1)
	# Access data in rastern band as array
	data = band.ReadAsArray(0,0,cols,rows)
	#
	# gdal interpolation creates "upside down files, hence this
	Yres=yres
	if Yres/np.fabs(Yres)==-1:
		Yres=-1*Yres
		data = np.flipud(data)
	#
	# get nodata value
	nandat = band.GetNoDataValue()
	# get minimum and maximum value
	mindat = band.GetMinimum()
	maxdat = band.GetMaximum()
	if mindat is None or maxdat is None:
			(mindat,maxdat) = band.ComputeRasterMinMax(1)
	dats = [nandat,mindat,maxdat]
	sumvals = data[np.where(np.logical_not(data == nandat))]
	sumval = sumvals.sum()
	numvals = len(sumvals)
	avval = sumval/numvals
	volvals = sumval*(math.fabs((transf[1]*transf[5])))
	meta = {'transform':transf,'projection':projec,'corners':corners,'dimension':dims,'dats':dats,'sumval':sumval,'numvals':numvals,'avval':avval,'volvals':volvals}
	if srs.IsProjected:
		print 'projcs 1: ',srs.GetAttrValue('projcs')
		print 'projcs 2: ',srs.GetAuthorityCode('projcs')
		print 'projcs 3: ',srs.GetAuthorityName('projcs')
	print 'geogcs 1:: ',srs.GetAttrValue('geogcs')
	print 'Sum of pixel values = ',sumval,' in range ',mindat,' to ',maxdat
	print 'From ',numvals,' of ',rows,' X ',cols,' = ',rows*cols,' pixels, or ',((1.0*numvals)/(rows*cols))*100,'%'
	print 'Average then is ',avval
	print 'Pixels size, x:',transf[1],' y:',transf[5]
	print '"Volume" is then the sum * 1 pixel area: ',volvals
	return data, meta, metadata
#
def rasterExport(outdata,meta,filnm):
	'''Export array as single layer GeoTiff. Uses metadata from rasterImport and numpy nan as mask.
	To use: rasterExport(outdata,meta,filnm)'''
	inrows =  meta['dimension']['rows']
	incols =  meta['dimension']['cols']
	datdriver = gdal.GetDriverByName( "GTiff" )
	datout = datdriver.Create(filnm,incols,inrows,1,gdal.GDT_Float32)
	datout.SetGeoTransform(meta['transform'])
	datout.SetProjection(meta['projection'])
	outdata_m = np.flipud(outdata)
	outdata_m = np.where(np.isnan(outdata_m),meta['dats'][0],outdata_m)
	datout.GetRasterBand(1).WriteArray(outdata_m)
	datout.GetRasterBand(1).SetNoDataValue(meta['dats'][0])
	datout.GetRasterBand(1).ComputeStatistics(True)
	datout = None
	return 0
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
	# Vectors of values along grid axes
	demRx = range(xstart,xstop,int(round(xres)))
	demRy = range(ystart,ystop,int(round(Yres)))
	# Make grid 2D arrays
	Xg1,Yg1 = np.meshgrid(demRx,demRy)
	# Convert grids to 1D vectors. NOTE THE STUPID NAMES Xg1 IS NOT 1D
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
def demslicer(file,outfolder,slice=10):
	'''Slice a raster (usually a DEM) into bands according to pixel value
	To use: demslicer(file,outfolder,slice=10) '''
	print slice
	print '\nImporting ',file
	time_one = datetime.now()
	print time_one.strftime('at day:%j %H:%M:%S')
	demdata,Xg,Yg,Xg1,Yg1,rx,ry,demmeta = DemImport(file)
	demmeta.append(['projection','geotransform','driver','rows','columns','nanvalue','min','max'])
	# Get range for slices by rounding range of data to nearest slice size
	# round to nearest whole
	floordat = np.floor(demmeta[7])
	ceildat = np.ceil(demmeta[8])
	# Divide to multiple
	floordat_10 = floordat / slice
	ceildat_10 = ceildat / slice
	# Round to multiple
	floordat_10_f = np.floor(floordat_10)
	ceildat_10_c = np.ceil(ceildat_10)
	# Multiply back to original scale
	floord = floordat_10_f * slice
	ceild = ceildat_10_c * slice
	# Create slice range
	slicerange = range(int(floord),int(ceild),int(slice))
	#
	# Slice the data
	for i in slicerange:
		data0 = np.zeros(demdata.shape, demdata.dtype)
		# Select array members with value within range
		slice1 = np.where(np.logical_and(demdata>=i,demdata<(i+slice)))
		# Set array data0 members indexed by slice1 to 1
		data0[slice1] = 1
		# Create output file
		name = str(i) + 'to' + str(i+slice)
		filnm = datawrite(data0,demdata,demmeta,name,outfolder)
		print filnm
	print '\n'
	return 0
#
def rasXras(flist,flistfolder,Afile,outfolder):
	'''Multiply one raster file with a list of other raster files and calculate "volume".
	To use: filnm = rasXras(flist,flistfolder,Afile,outfolder) '''
	time_one = datetime.now()
	print time_one.strftime('Start rasXras at:%j %H:%M:%S')
#	 print flist
#	 print flistfolder
#	 print Afile
#	 print outfolder
	# Create/Open text file for storing results
	resnam = os.path.join(outfolder,(time_one.strftime('%j_%H%M%S')+'_results.csv'))
	f = open(resnam,'a')
	f.write('file,sum,count,pixel_avg,pixelsize,area,volume \n')
	# Get first raster file
	Adata,AXg,AYg,AXg1,AYg1,Arx,Ary,Ademmeta = DemImport(Afile)
	print Afile,' opened as first file.'
	# Go through list of other raster files
	for B in flist:
		print 'Test: ',B
		Bfile = os.path.join(flistfolder,B)
		Bdata,BXg,BYg,BXg1,BYg1,Brx,Bry,Bdemmeta = DemImport(Bfile)
		print Bfile,' opened.'
		print 'Rows (A,B): ',Arx,Brx,' Columns (A,B): ',Ary,Bry
		print 'xres (A,B): ',Ademmeta[2][1],Bdemmeta[2][1],' yres (A,B): ',Ademmeta[2][5],Bdemmeta[2][5]
		print 'No data (A): ',Ademmeta[6],' No data (B): ',Bdemmeta[6]
		# Check matching resolution
		if Arx != Brx or Ary != Bry:
			print B,' resolution mismatch with ', Afile
			continue
		elif Ademmeta[4] != Bdemmeta[4] or Ademmeta[5] != Bdemmeta[5]:
			print 'Size mismatch between ',B,' and ',Afile
			continue
		# Multiply first file with current file
		AdataMasked = np.ma.masked_where(Adata == Ademmeta[6], Adata)
		BdataMasked = np.ma.masked_where(Bdata == Bdemmeta[6], Bdata)
		outdata = AdataMasked * BdataMasked
		name = namer(B)[0]+'X'+namer(Afile)[0]
		#print "A: ",AdataMasked
		#print "B: ",BdataMasked
		#print "outdata: ", outdata
		#print "Adata: ", Adata
		print "Ademmeta: ",Ademmeta
		print "name: ",name
		print "outfolder: ",outfolder
		filnm = datawrite(outdata,Adata,Ademmeta,name,outfolder)
		sumval = np.sum(outdata)
		print 'Sum of values: ',sumval,'\n'
		avval = np.mean(outdata)
		print 'Mean value ',avval,'\n'
		countval = sumval/avval
		print 'Number of values: ',countval
		pixsize = abs(Ademmeta[2][1] * Ademmeta[2][5])
		area = pixsize * countval
		volume = avval * area
		# write values to text file
		blnk = ','
		endln = '\n'
		sumvalstr = str(sumval)
		countvalstr = str(countval)
		avvalstr = str(avval)
		pixsizestr = str(pixsize)
		areastr = str(area)
		volumestr = str(volume)
		results = B.split('.')[0] + blnk + sumvalstr + blnk + countvalstr + blnk + \
		avvalstr + blnk + pixsizestr + blnk + areastr + blnk + volumestr + endln
		f.write(results)
	f.close()
	print '\n'
	return filnm
#
#
def rasterAverage(flist,flistfolder,outfolder):
	'''Average all rasters together
	To use: filnmsum,filnmmean = rasterAverage(flist,flistfolder,outfolder) '''
	time_one = datetime.now()
	print time_one.strftime('Start rasterAverage at:%j %H:%M:%S')
	# Get first raster file
	Afile = os.path.join(flistfolder,flist[0])
	Adata,AXg,AYg,AXg1,AYg1,Arx,Ary,Ademmeta = DemImport(Afile)
	# Get second for creating mask to deal with NaN
	Bfile = os.path.join(flistfolder,flist[1])
	Bdata,BXg,BYg,BXg1,BYg1,Brx,Bry,Bdemmeta = DemImport(Bfile)
	Adata[Adata==Ademmeta[6]]=np.nan
	Bdata[Bdata==Bdemmeta[6]]=np.nan
	sumdata = np.ma.masked_array(np.nan_to_num(Adata), mask=np.isnan(Adata) & np.isnan(Bdata))
	counter = 1
	print '\n'*3
	print '*'*20
	# Go through list of other raster files
	for B in flist[1:]:
		print 'file ',counter+1,' of ',len(flist)
		Bfile = os.path.join(flistfolder,B)
		Bdata,BXg,BYg,BXg1,BYg1,Brx,Bry,Bdemmeta = DemImport(Bfile)
		print B,' data type: ',type(Bdata)
		print 'Rows (A,B): ',Arx,Brx,' Columns (A,B): ',Ary,Bry
		print 'xres (A,B): ',Ademmeta[2][1],Bdemmeta[2][1],' yres (A,B): ',Ademmeta[2][5],Bdemmeta[2][5]
		print 'No data (A): ',Ademmeta[6],' No data (B): ',Bdemmeta[6]
		# Check matching resolution
		if Arx != Brx or Ary != Bry:
			print B,' resolution mismatch with ', Afile
			continue
		elif Ademmeta[4] != Bdemmeta[4] or Ademmeta[5] != Bdemmeta[5]:
			print 'Size mismatch between ',B,' and ',Afile
			continue
		# Add current file to sum
		Bdata[Bdata==Bdemmeta[6]]=np.nan
		BdataMasked  = np.ma.masked_array(np.nan_to_num(Bdata), mask=sumdata.mask)
		sumdata = (sumdata + BdataMasked).filled(np.nan)
		sumdata = np.ma.masked_array(np.nan_to_num(sumdata), mask=np.isnan(sumdata))
		counter = counter + 1
	meandata = sumdata/counter
	#
	sumname = 'file_sum'
	filnmsum = datawrite(sumdata,Adata,Ademmeta,sumname,outfolder)
	meanname = 'file_mean'
	filnmmean = datawrite(meandata,Adata,Ademmeta,meanname,outfolder)
	return filnmsum,filnmmean
#
def rasterDiff(Afile, Bfile, outfolder):
	'''Subtract one raster from another
	To use: filnm = rasSubras(Afile, Bfile, outfolder) '''
	A,Aext,Apath,Anamefull = namer(Afile)
	B,Bext,Bpath,Bnamefull = namer(Bfile)
	# Get first raster file
	Adata,AXg,AYg,AXg1,AYg1,Arx,Ary,Ademmeta = DemImport(Afile)
	Adata[Adata==Ademmeta[6]]=np.nan
	Avec = np.reshape(Adata,(Arx*Ary))
	AvecNan = Avec[np.logical_not(np.isnan(Avec))]
	# Get second
	Bdata,BXg,BYg,BXg1,BYg1,Brx,Bry,Bdemmeta = DemImport(Bfile)
	print '\n'*2
	print '*'*20
	print A,' data type: ', type(Adata)
	print B,' data type: ',type(Bdata)
	print 'Rows (A,B): ',Arx,Brx,' Columns (A,B): ',Ary,Bry
	print 'xres (A,B): ',Ademmeta[2][1],Bdemmeta[2][1],' yres (A,B): ',Ademmeta[2][5],Bdemmeta[2][5]
	print 'No data (A): ',Ademmeta[6],' No data (B): ',Bdemmeta[6]
	# Check matching resolution
	if Arx != Brx or Ary != Bry:
		print B,' resolution mismatch with ', Afile
	elif Ademmeta[4] != Bdemmeta[4] or Ademmeta[5] != Bdemmeta[5]:
		print 'Size mismatch between ',B,' and ',Afile
		return 1
	# Add current file to sum
	Bdata[Bdata==Bdemmeta[6]]=np.nan
	Bvec = np.reshape(Bdata,(Brx*Bry))
	BvecNan = Bvec[np.logical_not(np.isnan(Bvec))]
	#
	AdataMasked = np.ma.masked_array(np.nan_to_num(Adata), mask=np.isnan(Adata) & np.isnan(Bdata))
	BdataMasked  = np.ma.masked_array(np.nan_to_num(Bdata), mask=AdataMasked.mask)
	#
	ABcovmat = np.cov(AvecNan,BvecNan)
	ABcov = ABcovmat[0,1]
	Amean = np.nanmean(Adata)
	Astd = ABcovmat[0,0]**0.5
	print "Mean of A: %4.4f, %4.4f" % (Amean, Astd)
	Bmean = np.nanmean(Bdata)
	Bstd = ABcovmat[1,1]**0.5
	print "Mean of B: %4.4f, %4.4f" % (Bmean, Bstd)
	print "Covariance: %4.4f" % (ABcov)
	#
	outdata = (BdataMasked - AdataMasked).filled(np.nan)
	outmean = np.nanmean(outdata)
	outstd = np.nanstd(outdata)
	print "Mean of outdata: %4.4f, %4.4f" % (outmean, outstd)
	propagated = (ABcovmat[0,0] + ABcovmat[1,1] - 2*ABcovmat[0,1])**0.5
	print "StdDev of propagated: %4.4f" % (propagated)
	outdata = np.ma.masked_array(np.nan_to_num(outdata), mask=np.isnan(outdata))
	outname = B+'Sub'+A
	filnm = datawrite(outdata,Adata,Ademmeta,outname,outfolder)
	return filnm
#
def rasSubras(flist,flistfolder,Afile,outfolder):
	'''Subtract one file from others in list
	To use: filnm = rasSubras(flist,flistfolder,Afile,outfolder) '''
	#time_one = datetime.now()
	#print time_one.strftime('Start raster subtraction at:%j %H:%M:%S')
	# Get first raster file
	Adata,AXg,AYg,AXg1,AYg1,Arx,Ary,Ademmeta = DemImport(Afile)
	Amean = nanmean(Adata)
	# Get second for creating mask to deal with NaN
#	 Bfile = os.path.join(flistfolder,flist[1])
#	 Bdata,BXg,BYg,BXg1,BYg1,Brx,Bry,Bdemmeta = DemImport(Bfile)
#	 Adata[Adata==Ademmeta[6]]=np.nan
#	 Bdata[Bdata==Bdemmeta[6]]=np.nan
#	 AdataMasked = np.ma.masked_array(np.nan_to_num(Adata), mask=np.isnan(Adata) & np.isnan(Bdata))
	counter = 1
	print '\n'*3
	print '*'*20
	# Go through list of other raster files
	for B in flist:
		print 'file ',counter+1,' of ',len(flist)
		Bfile = os.path.join(flistfolder,B)
		Bdata,BXg,BYg,BXg1,BYg1,Brx,Bry,Bdemmeta = DemImport(Bfile)
		print B,' data type: ',type(Bdata)
		print 'Rows (A,B): ',Arx,Brx,' Columns (A,B): ',Ary,Bry
		print 'xres (A,B): ',Ademmeta[2][1],Bdemmeta[2][1],' yres (A,B): ',Ademmeta[2][5],Bdemmeta[2][5]
		print 'No data (A): ',Ademmeta[6],' No data (B): ',Bdemmeta[6]
		# Check matching resolution
		if Arx != Brx or Ary != Bry:
			print B,' resolution mismatch with ', Afile
			continue
		elif Ademmeta[4] != Bdemmeta[4] or Ademmeta[5] != Bdemmeta[5]:
			print 'Size mismatch between ',B,' and ',Afile
			continue
		# Add current file to sum
		Bdata[Bdata==Bdemmeta[6]]=np.nan
		AdataMasked = np.ma.masked_array(np.nan_to_num(Adata), mask=np.isnan(Adata) & np.isnan(Bdata))
		BdataMasked  = np.ma.masked_array(np.nan_to_num(Bdata), mask=AdataMasked.mask)
		outdata = (BdataMasked - AdataMasked).filled(np.nan)
		outdata = np.ma.masked_array(np.nan_to_num(outdata), mask=np.isnan(outdata))
		counter = counter + 1
		outname = namer(B)[0]+'Sub'+namer(Afile)[0]
		filnm = datawrite(outdata,Adata,Ademmeta,outname,outfolder)
	return filnm
#
def geopixsum(filename):
	'''Sum all the non NaN values in a raster file
	To Use:[sumval, area, average, countval] = geopixsum(filename) '''
	# register all of the GDAL drivers
	gdal.AllRegister()
	sumval = 'No File'
	# open the image
	try:
		inDs = gdal.Open(filename)
	except:
		print 'Could not open ',file,'\n'
	# get image size
	rows = inDs.RasterYSize
	cols = inDs.RasterXSize
	transf = inDs.GetGeoTransform()
	ul_x = transf[0]
	ul_y = transf[3]
	xres = transf[1]
	yres = transf[5]
	#print 'rows = ',rows,' cols = ',cols
	# read band 1 into data
	band1 = inDs.GetRasterBand(1)
	data = band1.ReadAsArray(0,0,cols,rows)
	print np.shape(data)
	# get nodata value
	nandat = band1.GetNoDataValue()
	print "NaN value: ", nandat
	sumvals = data[np.logical_not((np.isnan(data)) + (np.isinf(data)) + (data==nandat))]
	sumval = sumvals.sum()
	countval = len(sumvals)
	average = sumval/countval
	area = countval * abs(xres * yres)
	print "Sum = %2.3f, Area = %2.1f, Average = %2.3f, Number = %d" % (sumval, area, average, countval)
	inDs = None
	return [sumval, area, average, countval]
#
def xydiff(atype, bearing, distance):
	"""Calculate the x and y coordinates (relative) from HD and Bearing"""
	# Convert bearing to radians
	if atype == 'd':
		bearing = math.radians(bearing)
	elif atype == 'g':
		bearing = bearing * (math.pi / 200)
	elif atype == 'r':
		bearing = bearing
	#
	# Calculate difference in Easting
	ediff = math.sin(bearing)*distance
	# Calculate difference in Northing
	ndiff = math.cos(bearing)*distance
	return  ediff, ndiff
#
def hordist(atype,vangle,slpdist):
	"""Calculate the horisontal distance from the sloping distance and the vertical angle"""
	# Convert angle to radians
	if atype == 'd':
		vangle = math.radians(vangle)
	elif atype == 'g':
		vangle = vangle * (math.pi / 200)
	elif atype == 'r':
		vangle = vangle
	#
	# Calculate horisontal distance
	hdist = math.cos(vangle)*slpdist
	# Calculate height difference
	vdist = math.sin(vangle)*slpdist
	#
	return hdist, vdist
#
def pythag(Xdif,Ydif):
	'''Calculate hypotenuse of right angled triangle'''
	# Called by: pythag2p
	hypo = (Xdif**2 + Ydif**2)**0.5
	return hypo
#
def pythag2p(Xa,Ya,Xb,Yb):
	'''Calculate distance in plane between two points'''
	# Called by: inbindning, inskarning, freestation
	Xdif = Xb - Xa
	Ydif = Yb - Ya
	hypo = pythag(Xdif, Ydif)
	return hypo
#
def extender(ax,ay,bx,by,length):
	'''Extend a line between two by a set distance'''
	lenAB = pythag2p(ax,ay,bx,by)
	cx = bx + (bx - ax) / lenAB * length
	cy = by + (by - ay) / lenAB * length
	return cx,cy
#
def mustard(file,kernel):
	'''Run a kernel (window,matrix) over a single layer raster. Uses rasterImport and rasterExport.
	Kernel must be numpy array of odd number dimensions.
	To use: outdata = mustard(file,kernel) '''
	# Example high pass filter
	# kernel = np.array([[-1.0,-1.0,-1.0],[-1.0,9.0,-1.0],[-1.0,-1.0,-1.0]])
	# Import data file
	data, meta, metadata = rasterImport(file)
	# Set no data values to numpy nan
	data[data==meta['dats'][0]]=np.nan
	# Get size of indata array
	inrows =  meta['dimension']['rows']
	incols =  meta['dimension']['cols']
	# Get size of kernel
	krows = np.shape(kernel)[0]
	kcols = np.shape(kernel)[1]
	# Check kernel is smaller than data grid.
	if krows >= inrows or kcols >= incols: sys.exit('Bad kernel. Too large.')
	# Check kernel has a central pixel
	if krows % 2 == 0: sys.exit('Bad kernel. Even number rows.')
	if kcols % 2 == 0: sys.exit('Bad kernel. Even number columns.')
	# Get central pixel location in kernel
	kmidrow = int(krows)/2
	kmidcol = int(kcols)/2
	# Create relative extent of kernel
	rowminext = -1*((krows-1)/2)
	rowmaxext = (krows-1)/2
	rowrange = range(rowminext,rowmaxext+1,1)
	colminext = -1*((kcols-1)/2)
	colmaxext = (kcols-1)/2
	colrange = range(colminext,colmaxext+1,1)
	# Set initial starting location of kernel on grid
	dati = kmidrow
	datj = kmidcol
	# Get number of rows to run kernel over
	gridrows = range(inrows + rowminext*2)
	# Get number of columns to run kernel over
	gridcols = range(incols + colminext*2)
	# Create output array filled with nan
	outdata = np.empty((inrows,incols))
	outdata[:] = np.nan
	# Start loop
	for row in gridrows:
		datj = kmidcol
		rowvec = np.ones((1,krows))*dati + rowrange
		for col in gridcols:
			if np.isnan(data[dati,datj]):
				datj = datj + 1
				outdata[dati,datj] = np.nan
			else:
				colvec = np.ones((1,kcols))*datj + colrange
				extract = np.empty((krows,kcols))
				for i in range(krows):
					for j in range(kcols):
						extract[i,j] = data[int(rowvec[0,i]),int(colvec[0,j])]
				pixval = np.nansum(extract * kernel)
				outdata[dati,datj] = pixval
				datj = datj + 1
		dati = dati + 1
	return outdata,meta
#
def kernelChoice(ans = ''):
	'''Returns a kernel to use with mustard as well as the name of the kernel.
	To use: kernel, ans = kernelChoice(ans = '') '''
	kernelDict ={
	'Average':np.array([[0.111,0.111,0.111],[0.111,0.111,0.111],[0.111,0.111,0.111]]),
	'HighPass1':np.array([[-0.7,-1.0,-0.7],[-1.0,6.8,-1.0],[-0.7,-1.0,-0.7]]),
	'HighPass2':np.array([[-1.0,-1.0,-1.0],[-1.0,9.0,-1.0],[-1.0,-1.0,-1.0]]),
	'PrewittX':np.array([[-1.0,0.0,1.0],[-1.0,0.0,1.0],[-1.0,0.0,1.0]]),
	'PrewittY':np.array([[1.0,1.0,1.0],[0.0,0.0,0.0],[-1.0,-1.0,-1.0]]),
	'SobelX':np.array([[-1.0,0.0,1.0],[-2.0,0.0,2.0],[-1.0,0.0,1.0]]),
	'SobelY':np.array([[1.0,2.0,1.0],[0.0,0.0,0.0],[-1.0,-2.0,-1.0]])
	}
	KerList = kernelDict.keys()
	KerList.sort()
	while ans not in KerList:
		print KerList
		ans = raw_input('Enter name of kernel: ')
	kernel = kernelDict[ans]
	print ans,':\n', kernel
	return kernel, ans
#
def curve(file,kernel):
	'''Copy of mustard, altered to calculate curvature according to Equation 4 from Park, S.; McSweeney, K. & Lowery, B. Identification of the spatial distribution of soils using a process-based terrain characterization Geoderma, 2001, 103, 249-272
	as suggested by http://casoilresource.lawr.ucdavis.edu/drupal/node/937
	Calls: rasterImport
	To Use: outdata, meta = curve(file,kernel) '''
	# Import data file
	data, meta, metadata = rasterImport(file)
	# Set no data values to numpy nan
	data[data==meta['dats'][0]]=np.nan
	# Get size of indata array
	inrows =  meta['dimension']['rows']
	incols =  meta['dimension']['cols']
	cellsize =  meta['dimension']['xres']
	# Get size of kernel
	krows = np.shape(kernel)[0]
	kcols = np.shape(kernel)[1]
	# Check kernel is smaller than data grid.
	if krows >= inrows or kcols >= incols: sys.exit('Bad kernel. Too large.')
	# Check kernel has a central pixel
	if krows % 2 == 0: sys.exit('Bad kernel. Even number rows.')
	if kcols % 2 == 0: sys.exit('Bad kernel. Even number columns.')
	# Get central pixel location in kernel
	kmidrow = int(krows)/2
	kmidcol = int(kcols)/2
	# Create relative extent of kernel
	rowminext = -1*((krows-1)/2)
	rowmaxext = (krows-1)/2
	rowrange = range(rowminext,rowmaxext+1,1)
	colminext = -1*((kcols-1)/2)
	colmaxext = (kcols-1)/2
	colrange = range(colminext,colmaxext+1,1)
	# Set initial starting location of kernel on grid
	dati = kmidrow
	datj = kmidcol
	# Get number of rows to run kernel over
	gridrows = range(inrows + rowminext*2)
	# Get number of columns to run kernel over
	gridcols = range(incols + colminext*2)
	# Create output array filled with nan
	outdata = np.empty((inrows,incols))
	outdata[:] = np.nan
	# Start loop
	for row in gridrows:
		datj = kmidcol
		rowvec = np.ones((1,krows))*dati + rowrange
		for col in gridcols:
			if np.isnan(data[dati,datj]):
				datj = datj + 1
				outdata[dati,datj] = np.nan
			else:
				colvec = np.ones((1,kcols))*datj + colrange
				extract = np.empty((krows,kcols))
				for i in range(krows):
					for j in range(kcols):
						extract[i,j] = data[int(rowvec[0,i]),int(colvec[0,j])]
					D = np.zeros(np.shape(extract))
					for i in range(krows):
						for j in range(kcols):
			# To here same as mustard
							D[i,j] = (((kmidrow - i)**2+(kmidcol - j)**2)**0.5) * cellsize
					D[kmidrow,kmidcol] = np.nan
					E = np.ones(np.shape(extract))
					E = E*extract[kmidrow,kmidcol]
					EmA = E - extract
					Asum1 = EmA/D
					Asum2 = np.nansum(Asum1)
					pixval = Asum2 / (np.size(extract)-1)
				outdata[dati,datj] = pixval
				datj = datj + 1
		dati = dati + 1
	return outdata, meta
#
def curveKernel():
	'''Creates an rXc kernel of ones. Called by curveDem for use with curve.
	Calls:
	To Use: kernel = curveKernel() '''
	print 'Define size of window to calculate curve from.\nMust have odd number of rows and columns'
	Ashape_r = ''
	Ashape_c = ''
	while not Ashape_r:
		ans = raw_input('Enter number of rows: ')
		try:
			Ashape_r = int(ans)
		except:
			Ashape_r = ''
	while not Ashape_c:
		ans = raw_input('Enter number of columns: ')
		try:
			Ashape_c = int(ans)
		except:
			Ashape_c = ''
	kernel = np.ones((Ashape_r,Ashape_c))
	return kernel
#
def curveDem(*infile):
	'''Create a curvature surface from a DEM.
	Calls: curveKernel, curve, namer, rasterExport
	To use: curveDem() '''
	if infile:
		file = infile[0]
	else:
		file = ''
	while os.path.isfile(file) is False:
		file = raw_input('Enter DEM file path and name: ')
	name,ext,path,namefull = namer(file)
	filnm = os.path.join(path,(name + '_curve.tif'))
	kernel = curveKernel()
	outdata, meta = curve(file,kernel)
	rasterExport(outdata,meta,filnm)
	return 0
#
def pt2fmt(pt):
	fmttypes = {
		GDT_Byte: 'B',
		GDT_Int16: 'h',
		GDT_UInt16: 'H',
		GDT_Int32: 'i',
		GDT_UInt32: 'I',
		GDT_Float32: 'f',
		GDT_Float64: 'f'
		}
	return fmttypes.get(pt, 'x')
#
#
def getRasterFile(file):
	'''Get georeferenced raster and open for use
	To use: raster, transf, bandcount = getRasterFile(file)'''
	# Called by:
	# register all of the GDAL drivers
	gdal.AllRegister()
	# Open file
	raster = gdal.Open(file, GA_ReadOnly)
	if raster is None:
		print 'Could not open ',file,'\n'
		sys.exit(1)
	# Get coordinate system parameters
	projec = raster.GetProjection()
	srs=osr.SpatialReference(wkt=projec)
	transf = raster.GetGeoTransform()
	bandcount = raster.RasterCount
	#
	return raster, transf, bandcount
#
def GetRasterVals(x,y, raster, transf, bandcount):
	'''Create vector of pixel values at location in raster
	To use: vals = GetRasterVals(x,y, raster, transf, bandcount)'''
	# Called by:
	# get image size
	#print x, y
	success, transfInv = gdal.InvGeoTransform(transf)
	if not success:
		print "Failed InvGeoTransform()"
		sys.exit(1)
	rows = raster.RasterYSize
	cols = raster.RasterXSize
#	 xdifpix = math.floor((x - transf[0])/transf[1])
#	 ydifpix = math.floor((y - transf[3])/transf[5])
#	 xpix = xdifpix -1
#	 ypix = ydifpix -1
	xpix, ypix = gdal.ApplyGeoTransform(transfInv, x, y)
	# Read the file band to a matrix called band_1
	vals = []
	for i in range(1,bandcount+1):
		band = raster.GetRasterBand(i)
		bandtype = gdal.GetDataTypeName(band.DataType)
		if band is None:
			continue
		# Access data in raster band as array
		#data = band.ReadAsArray(0,0,cols,rows)
		#vals[i] = (data[ypix,xpix])
		structval = band.ReadRaster(int(xpix), int(ypix), 1,1, buf_type = band.DataType )
		fmt = pt2fmt(band.DataType)
		intval = struct.unpack(fmt , structval)
		vals[i] = intval[0]
	return vals
#
def datawrite(outdata,demdata,meta,name,outDir):
	'''Write an array of grid data to a georeferenced raster file'''
	#meta = ['projection','geotransform','driver','rows','columns','nanvalue']
	filnm = os.path.join(outDir,(name + '.tif'))
	datdriver = gdal.GetDriverByName( "GTiff" )
	datout = datdriver.Create(filnm,meta[5],meta[4],1,gdal.GDT_Float32)
	datout.SetGeoTransform(meta[2])
	datout.SetProjection(meta[1])
	nanmask = demdata != meta[6]
	outdata_m = np.flipud(outdata * nanmask)
	outdata_m = np.where(outdata_m==0,-9999,outdata_m)
	datout.GetRasterBand(1).WriteArray(outdata_m)
	datout.GetRasterBand(1).SetNoDataValue(-9999)
	datout.GetRasterBand(1).ComputeStatistics(True)
	datout = None
	print "datawrite returns: ",filnm
	return filnm
#
def vrtMaker(fileName):
	'''Read csv coordinate file and create vrt header
	To use: Zfield, namevrt = vrtMaker(datafile)'''
	print '\n'*2,'*'*10
	print 'Read ',fileName
	InFile = open(fileName,'rb')
	Headers = InFile.next().strip().split(',')
	InFile.close()
	warn = 0
	for i in range(len(Headers)):
		if len(Headers[i]) < 2:
			print '"gdal_grid" cannot parse single letter column headers.\n'
			print 'If you intend to pass this column as x,y or z then change the column name.'
			warn = 1
		if ' ' in Headers[i] or '.' in Headers[i]:
			print '"gdal_grid" cannot parse column headers with spaces or dots.\n'
			print 'If you intend to pass this column as x,y or z then change the column name.\n'
			warn = 1
	print 'File Headers:\n'
	for i in Headers:
		print i
	print '\n'
	if warn == 1:
		print 'THERE ARE BAD COLUMN HEADERS IN YOUR INPUT FILE\n'
		ans = cont()
		if ans == 1:
			print '\nLeaving makevrt function\n'
			return 0
	Xcol = ''
	Ycol = ''
	Zcol = ''
	useAns = ['y','n']
	useFirst = ''
	while useFirst not in useAns:
		useFirst = raw_input('Use first three fields as Easting, Northing and Z? (y/n): ')
	if useFirst == 'n':
		while Xcol not in Headers:
			Xcol = raw_input('Enter column for "Easting": ')
		while Ycol not in Headers:
			Ycol = raw_input('Enter column for "Northing": ')
		while Zcol not in Headers and Zcol !='none':
			print Zcol
			Zcol = raw_input('Enter column for Z or "none": ')
	elif useFirst == 'y':
		Xcol = 'field_1'
		Ycol = 'field_2'
		Zcol = 'field_3'
	epsgList = ['3006','3021','7030']
	nameList = ['SWEREF99TM','RT90 2,5gV','WGS84']
	for i in range(len(nameList)):
		print nameList[i],' = ',epsgList[i]
	epsg = ''
	while epsg not in epsgList:
		epsg = raw_input('Enter epsg code for source file coordinate system: ')
	namevrt = makevrt(fileName,epsg,Xcol,Ycol,Zcol)
	return Zcol, namevrt, epsg
#
def makevrt(fileName,epsg,Xcol,Ycol,Zcol='none'):
	layername,ext,path,namefull = namer(fileName)
	line1 = '<OGRVRTDataSource>\n'
	line2 = '<OGRVRTLayer name="'+layername+'">\n'
	line3 = '<SrcDataSource>'+layername+'.csv</SrcDataSource>\n'
	line4 = '<GeometryType>wkbPoint</GeometryType>\n'
	line4a = '<LayerSRS>EPSG:'+epsg+'</LayerSRS>\n'
	if Zcol != 'none': line5 = '<GeometryField encoding="PointFromColumns" x="'+Xcol+'" y="'+Ycol+'" z="'+Zcol+'"/>\n'
	elif Zcol == 'none': line5 = '<GeometryField encoding="PointFromColumns" x="'+Xcol+'" y="'+Ycol+'"/>\n'
	else: return 1
	line6 = '</OGRVRTLayer>\n'
	line7 = '</OGRVRTDataSource>\n'
	#
	vrtName = layername + '.vrt'
	outName = os.path.join(path,vrtName)
	with open(outName,'ab') as OutFile:
		OutFile.write(line1)
		OutFile.write(line2)
		OutFile.write(line3)
		OutFile.write(line4)
		OutFile.write(line4a)
		OutFile.write(line5)
		OutFile.write(line6)
		OutFile.write(line7)
	return outName
#
def idw(dataFile,DEM):
	'''Implement gdalgrid for IDW'''
	# Set output directory
	namstr,ext_,outDir,full_ = namer(dataFile)
	outDir = os.path.join(outDir,namstr)
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	# Make vrt for extrapolation
	Zfield, namevrt, epsg = vrtMaker(dataFile)
	# Get metadata from DEM for setting extrapolation extents
	_, meta, metadata = rasterImport(DEM)
	txel = meta['corners']['ll'][0]
	txeu = meta['corners']['ur'][0]
	tyel = meta['corners']['ll'][1]
	tyeu = meta['corners']['ur'][1]
	outsider = meta['dimensions']['rows']
	outsidec = meta['dimensions']['cols']
	# Parameters used in IDW
	power = 2.0
	smoothing = 2.0
	# Create gdal command for IDW
	layer = namestr
	idwfile = os.path.join(outDir, namestr + '_idw.tif')
	gridout = "gdal_grid -a invdist:power={0:0.2f}:smoothing={1:0.2f} -txe {2:0.0f} {3:0.0f} -tye {4:0.0f} {5:0.0f} -outsize {5:0.0f} {6:0.0f} -a_srs EPSG:{7} -of GTiff -ot Float64 -l {8} {9} {10}".format(power, smoothing, txel, txeu, tyel, tyeu, outsider, outsidec, epsg, layer, namevrt, idwfile)
	os.system(gridout)
	return idwfile
