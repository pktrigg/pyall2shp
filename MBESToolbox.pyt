import arcpy
from collections import namedtuple
import sys
import time
import os
import fnmatch
import math
from datetime import datetime
from datetime import timedelta
from glob import glob
import pyall
import shapefile
import geodetic

###############################################################################
class Toolbox(object):
	def __init__(self):
		"""Define the toolbox (the name of the toolbox is the name of the
		.pyt file)."""
		self.label = "MBES Toolbox"
		self.alias = ""
		# List of tool classes associated with this toolbox
		self.tools = [all2shp]

###############################################################################
class all2shp(object):
	def __init__(self):
		"""Define the tool (tool name is the name of the class)."""
		self.label = "Kongsberg ALL file coverage extraction V1.17"
		self.description = "Kongsberg .ALL file coverage and trackplot tool"
		self.canRunInBackground = False

###############################################################################
	def getParameterInfo(self):
		# Input Features parameter
		param0 = arcpy.Parameter(
			displayName="Input Folder (containing .all files)",
			name="inputfolder",
			datatype="DEFolder",
			parameterType="Required",
			direction="Input")

		param1 = arcpy.Parameter(
			displayName="Process .all files in sub-folders",
			name="recursive",
			datatype="GPBoolean",
			parameterType="Required",
			direction="Input")
		
		param2 = arcpy.Parameter(
			displayName="Step size between data records.  use 0 for every record which is highly detailed (seconds: default 30)",
			name="step",
			datatype="GPLong",
			parameterType="Required",
			direction="Input")
		
		param3 = arcpy.Parameter(
			displayName="Create Coverage Polygons",
			name="createcoverage",
			datatype="GPBoolean",
			parameterType="Required",
			direction="Input")

		param4 = arcpy.Parameter(
			displayName="Create Track Lines",
			name="createtracklines",
			datatype="GPBoolean",
			parameterType="Required",
			direction="Input")

		param5 = arcpy.Parameter(
			displayName="Create Track Points",
			name="createtrackpoints",
			datatype="GPBoolean",
			parameterType="Required",
			direction="Input")

		# param5 = arcpy.Parameter(
		# 	displayName="Create Track Points",
		# 	name="createtrackpoints",
		# 	datatype="GPBoolean",
		# 	parameterType="Required",
		# 	direction="Input")

		param6 = arcpy.Parameter(
			displayName="Output File Name:",
			name="output",
			datatype="GPString",
			parameterType="Required",
			direction="Input")

		param7 = arcpy.Parameter(
			displayName="Discovered files to process:",
			name="filestoprocess",
			datatype="GPString",
			parameterType="Required",
			direction="Output",
			enabled=False)
			
		# param0.value = 'c:\\development\\python'
		param1.value = False
		param2.value = 10
		param3.value = True
		param4.value = True
		param5.value = True
		# param5.value = True
		param6.value = "coverage"
		
		parameters = [param0, param1, param2, param3, param4, param5, param6, param7]
		return parameters

###############################################################################
	def isLicensed(self):
		"""Set whether tool is licensed to execute."""
		return True

###############################################################################
	def updateParameters(self, parameters):
		"""Modify the values and properties of parameters before internal
		validation is performed.  This method is called whenever a parameter
		has been changed."""

		if not os.path.exists(str(parameters[0].valueAsText)):
			return
		files = findfiles(parameters[0].valueAsText + "\\*.all", parameters[1].value)
		txt = "Count:" + str(len(files)) + ", "
		for f in files:
			txt = txt + f + ", "
			if len(txt) > 128:
				break
		parameters[7].value = txt
		return

###############################################################################
	def updateMessages(self, parameters):
		"""Modify the messages created by internal validation for each tool
		parameter.  This method is called after internal validation."""
		return

###############################################################################
	def execute(self, parameters, messages):
		"""The source code of the tool."""
		inputFolder = parameters[0].valueAsText + "\\*.all"
		recursive = parameters[1].value
		stepsize = int(parameters[2].valueAsText)
		createcoverage = parameters[3].value
		createtrackline = parameters[4].value
		createtrackpoint = parameters[5].value
		# createraster = parameters[5].value
		outfile = parameters[6].value

		arcpy.AddMessage(inputFolder)
		arcpy.AddMessage("Coverage:"+ str(createcoverage))
		arcpy.AddMessage("Lines:"+ str(createtrackline))
		arcpy.AddMessage("Points:"+ str(createtrackpoint))
		# arcpy.AddMessage("raster:"+ str(createraster))

		if createcoverage == False and createtrackline == False and createtrackpoint == False:
			arcpy.AddMessage("Nothing to do, quitting.")

		# we have some inputs, so we can call the processing script.
		# we need to ensure the file is a shp extension
		args = {'inputFolder': inputFolder, 'step': stepsize, 'coverage': createcoverage, 'recursive': recursive, 'trackline':createtrackline, 'trackpoint':createtrackpoint, 'outputFile': outfile, 'csv':False}
		# args = {'inputFile': inputFolder, 'step': stepsize, 'coverage': createcoverage, 'recursive': False, 'trackline':createtrackline, 'trackpoint':createtrackpoint, 'raster':createraster, 'outputFile': outfile, 'csv':False}
		a = namedtuple('GenericDict', args.keys())(**args)
		
		process(a)

		return

###############################################################################
def findfiles(inputFolder, recursive):
	matches = []
	if recursive:
		for root, dirnames, filenames in os.walk(os.path.dirname(inputFolder)):
			for f in fnmatch.filter(filenames, '*.all'):
				matches.append(os.path.join(root, f))
	else:
		if os.path.isfile(inputFolder):
			matches.append (os.path.abspath(inputFolder))
			# arcpy.AddMessage ("2 Adding:" + os.path.abspath(inputFile))
		else:
			for filename in glob(inputFolder):
				if os.path.isfile(filename):
					# arcpy.AddMessage ("Adding:" + filename)
					matches.append(filename)

	arcpy.AddMessage ("Files to process:" + str(matches))
	return matches

###############################################################################
def process(args):
	arcpy.AddMessage("Processing files...")
	matches = findfiles(args.inputFolder, args.recursive)

	# there are no files to process, so quit
	if len(matches) == 0:
		arcpy.AddMessage ("Nothing to process, quitting...")
		exit(0)

	fname, ext = os.path.splitext(os.path.expanduser(args.outputFile))

	outputfolder = os.path.dirname(args.inputFolder)
	# outputfolder = os.path.dirname(os.path.abspath(matches[0]))
	arcpy.AddMessage ("output folder:" + outputfolder)
	trackLineFileName = os.path.join(outputfolder, fname + "_trackLine.shp")
	trackPointFileName = os.path.join(outputfolder, fname + "_trackPoint.shp")
	trackCoverageFileName = os.path.join(outputfolder, fname + "_trackCoverage.shp")
	rasterFileName = os.path.join(outputfolder, fname + "_raster")

	# open the output files once only.
	# create the destination shape files 
	if args.trackpoint:
		TPshp = createSHP(trackPointFileName, shapefile.POINT)
		if len(TPshp.fields) <= 1: #there is a default deletion flag field always set, to we need to accoutn for this.
			TPshp.field("LineName", "C")
			TPshp.field("SurveyDate", "C")
			TPshp.field("SurveyTime", "C")
			TPshp.field("UNIXTime", "N")
			TPshp.field("SpeedKnots", "C")
			TPshp.field("Heading", "C")
			TPshp.field("PortCover", "N")
			TPshp.field("StbdCover", "N")
			TPshp.field("PortWidth", "N")
			TPshp.field("StbdWidth", "N")
			TPshp.field("DepthMode", "C")
			TPshp.field("Absorption", "C")
			TPshp.field("PulseLength", "N")
			TPshp.field("TVG", "N")
			TPshp.field("DualSwath", "C")
			TPshp.field("SpikeFilt", "C")
			TPshp.field("Stabilise", "C")
			TPshp.field("MinZGate", "N")
			TPshp.field("MaxZGate", "N")
			TPshp.field("BeamSpace", "C")
			TPshp.field("MedianDepth", "C")

	if args.trackline:
		TLshp = createSHP(trackLineFileName, shapefile.POLYLINE)
		if len(TLshp.fields) <= 1:
			TLshp.field("LineName", "C")
			TLshp.field("SurveyDate", "C")

	if args.coverage:
		TCshp = createSHP(trackCoverageFileName, shapefile.POLYGON)
		if len(TCshp.fields) <= 1:
			TCshp.field("LineName", "C")
			TCshp.field("SurveyDate", "C")

	fileCounter=0		
	for filename in matches:
		# this is not a .all file so skip
		if not filename.lower().endswith('.all'):
			fileCounter +=1
			continue
		
		reader = pyall.ALLReader(filename)
		start_time = time.time() # time  the process

		# create the track point with a point recoard and metadata per ping
		if args.trackpoint:
			createTrackPoint(reader, TPshp, float(args.step))

		# create the track polyline
		if args.trackline:
			createTrackLine(reader, TLshp, float(args.step))

		# create the coverage polygon
		if args.coverage:
			createCoverage(reader, TCshp, float(args.step))

		# create the csv polygon of coverage
		if args.csv:
			createCSV(reader, csv, float(args.step), filename)

		fileCounter +=1
		reader.close()

	if args.coverage:
		arcpy.AddMessage ("Saving coverage polygon shapefile: %s" % trackCoverageFileName)		
		TCshp.save(trackCoverageFileName)
		# now write out a prj file so the data has a spatial Reference
		prj = open(trackCoverageFileName.replace('.shp','.prj'), 'w')
		prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
		prj.close() # you can omit in most cases as the destructor will call it
		addshapetoArcMap(trackCoverageFileName, False)

	if args.trackline:
		arcpy.AddMessage ("Saving track line shapefile: %s" % trackLineFileName)		
		TLshp.save(trackLineFileName)
		# now write out a prj file so the data has a spatial Reference
		prj = open(trackLineFileName.replace('.shp','.prj'), 'w')
		prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
		prj.close() # you can omit in most cases as the destructor will call it
		addshapetoArcMap(trackLineFileName, False)

	if args.trackpoint:
		arcpy.AddMessage ("Saving track point shapefile: %s" % trackPointFileName)		
		TPshp.save(trackPointFileName)
		# now write out a prj file so the data has a spatial Reference
		prj = open(trackPointFileName.replace('.shp','.prj'), 'w')
		prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
		prj.close() # you can omit in most cases as the destructor will call it
		addshapetoArcMap(trackPointFileName, True)
	return

###############################################################################
def addshapetoArcMap(fileName, refresh=True):
	#if the shape file exists, and we are in arcmap, we can add the new layers to the map document.
	if 'ArcMap' in sys.executable:
		arcpy.AddMessage ("Adding layers to arcmap" )		
		if os.path.isfile(fileName):
			try:
				mxd = arcpy.mapping.MapDocument("CURRENT")
				df = arcpy.mapping.ListDataFrames(mxd,"*")[0]
				newlayer1 = arcpy.mapping.Layer(fileName)
				arcpy.mapping.AddLayer(df, newlayer1, "TOP")  
				if refresh:
					arcpy.RefreshActiveView()
					arcpy.RefreshTOC()
			except Exception as ex:
				print(ex.args[0])
	return

###############################################################################
def createTrackPoint(reader, shp, step):
	lastTimeStamp = 0
	recTime = 0
	recTimeString = ""
	latitude = 0
	longitude = 0
	heading = 0
	points = []
	selectedPositioningSystem = None
	maximumPortCoverageDegrees = 0
	maximumPortWidth = 0
	maximumStbdCoverageDegrees = 0
	maximumStbdWidth = 0
	depthmode = ""
	absorptioncoefficient = 0
	pulselength = 0
	tvg = 0
	dualswath = "N/A"
	spikefilter = "N/A"
	stabilisation = "N/A"
	mindepthgate = 0
	maxdepthgate = 0
	beamspacing = "N/A"
	mediandepth = 0

	# remember the previous records so we can compute the speed
	prevX = 0
	prevY = 0
	prevT = 0
	arcpy.AddMessage("Creating Track Points for:" + reader.fileName)

	reader.rewind() #rewind to the start of the file
	while reader.moreData():
		TypeOfDatagram, datagram = reader.readDatagram()

		if TypeOfDatagram == 'D' or TypeOfDatagram == 'X':
			datagram.read()
			if len(datagram.Depth) == 0:
				continue
			mediandepth = (median(datagram.Depth))

		if (TypeOfDatagram == 'P'):
			datagram.read()
			if (selectedPositioningSystem == None):
				selectedPositioningSystem = datagram.Descriptor
			if (selectedPositioningSystem == datagram.Descriptor):
				latitude = datagram.Latitude
				longitude = datagram.Longitude
				heading = datagram.Heading
				recTime = to_timestamp(reader.currentRecordDateTime())

		if TypeOfDatagram == 'R':
			datagram.read()
			maximumPortCoverageDegrees = datagram.maximumPortCoverageDegrees
			maximumPortWidth = datagram.maximumPortWidth
			maximumStbdCoverageDegrees = datagram.maximumStbdCoverageDegrees
			maximumStbdWidth = datagram.maximumStbdWidth
			depthmode = datagram.DepthMode + "+" + datagram.TXPulseForm
			absorptioncoefficient = datagram.absorptionCoefficient
			pulselength = datagram.transmitPulseLength
			tvg = datagram.tvg
			dualswath = datagram.dualSwathMode
			spikefilter = datagram.filterSetting
			stabilisation = datagram.yawAndPitchStabilisationMode
			mindepthgate = datagram.minimumDepth
			maxdepthgate = datagram.maximumDepth
			beamspacing = datagram.beamSpacingString

		if recTime - lastTimeStamp > step:
			if latitude == 0 or longitude == 0:
				continue
			lastTimeStamp = recTime
			shp.point(longitude,latitude)
			# now add to the shape file.
			recDate = from_timestamp(recTime).strftime("%Y-%m-%d")
			recTimeString = from_timestamp(recTime).strftime("%H:%M:%S")
			# write out the shape file FIELDS data
			# compute the speed as distance/time
			distance = math.sqrt( ((longitude - prevX) **2) + ((latitude - prevY) **2))
			dtime = max(lastTimeStamp - prevT, 0.001)
			speed = int((distance/dtime) * 60.0 * 3600) # need to convert from degrees to knots
			shp.record(os.path.basename(reader.fileName), 
				recDate, 
				recTimeString, 
				int(lastTimeStamp), 
				round(speed,2), 
				round(heading,2), 
				int(maximumPortCoverageDegrees), 
				int(maximumStbdCoverageDegrees), 
				int(maximumPortWidth), 
				int(maximumStbdWidth), 
				depthmode, 
				round(absorptioncoefficient,2), 
				int(pulselength), 
				int(tvg), 
				dualswath, 
				spikefilter, 
				stabilisation, 
				int(mindepthgate), 
				int(maxdepthgate), 
				beamspacing,
				round(mediandepth,2))

			# remember the last update
			prevX = longitude
			prevY = latitude
			prevT = lastTimeStamp

###############################################################################
def median(lst):
	n = len(lst)
	if n < 1:
			return None
	if n % 2 == 1:
			return sorted(lst)[n//2]
	else:
			return sum(sorted(lst)[n//2-1:n//2+1])/2.0

###############################################################################
def normalizeAngle(angle):
	newAngle = angle
	while (newAngle <= -180):
		newAngle += 360
	while (newAngle > 180):
		newAngle -= 360
	return newAngle

###############################################################################
def createTrackLine(reader, trackLine, step):
	lastTimeStamp = 0
	line_parts = []
	line = []
	navigation = reader.loadNavigation()
	arcpy.AddMessage("Creating Track Line for:" + reader.fileName)

	# create the trackline shape file
	for update in navigation:
		if update[0] - lastTimeStamp >= step:
			line.append([float(update[2]),float(update[1])])
			lastTimeStamp = update[0]
	# now add the very last update
	line.append([float(navigation[-1][2]),float(navigation[-1][1])])
		
	line_parts.append(line)
	trackLine.line(parts=line_parts)
	# now add to the shape file.
	recDate = from_timestamp(navigation[0][0]).strftime("%Y-%m-%d")
	# write out the shape file FIELDS data
	trackLine.record(os.path.basename(reader.fileName), recDate) 

###############################################################################
def createCoverage(reader, coveragePoly, step):
	lastTimeStamp = 0
	lastheading = None # remember the heading so we can chop up into smaller polygons.  This is important if the vessel is running in circles rather than lines.
	left = []
	right = []
	selectedPositioningSystem = None
	latitude = 0;
	longitude = 0
	heading = []
	leftside = [] #sliding window
	rightside = [] #sliding window
	window = step #sliding window size in number of pings, to smooth the data
	pendingrecord = False

 	arcpy.AddMessage("Creating Coverage for:" + reader.fileName)

	reader.rewind()
	while reader.moreData():
		TypeOfDatagram, datagram = reader.readDatagram()
		if (TypeOfDatagram == 'P'):
			datagram.read()
			if (selectedPositioningSystem == None):
				selectedPositioningSystem = datagram.Descriptor
			if (selectedPositioningSystem == datagram.Descriptor):
				latitude = datagram.Latitude
				longitude = datagram.Longitude
				if lastheading == None:
					lastheading = datagram.Heading
				pendingrecord = True #performance upgrade

		if pendingrecord == True:
			if TypeOfDatagram == 'D' or TypeOfDatagram == 'X':
				datagram.read()
				if len(datagram.AcrossTrackDistance) == 0:
					continue
				if len(datagram.Depth) == 0:
					continue
				nadirbeamno = math.floor(len(datagram.Depth)/2)
				if (math.fabs(datagram.AcrossTrackDistance[0]) > 0) and (math.fabs(datagram.AcrossTrackDistance[-1]) > 0):
					leftside.append(min(datagram.AcrossTrackDistance))
					rightside.append(max(datagram.AcrossTrackDistance))
					heading.append(datagram.Heading)
					pendingrecord = False #performance upgrade - only decode teh bathy if needed.

		# add to the shape file at the user required interval
		if to_timestamp(reader.currentRecordDateTime()) - lastTimeStamp >= step:
			if len(leftside) == 0 or len(rightside) == 0:
				continue
			if longitude == 0 or latitude == 0:
				continue
			leftextent = median(leftside)
			# leftextent = statistics.median(leftside)
			rightextent = median(rightside)
			# rightextent = statistics.median(rightside)
			hdg = median(heading)
			# hdg = statistics.median(heading)
			leftLatitude, leftLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, hdg - 90 , math.fabs(leftextent))
			rightLatitude, rightLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, hdg + 90, math.fabs(rightextent))
			left.append([leftLongitude,leftLatitude])
			right.append([rightLongitude,rightLatitude])
			lastTimeStamp = to_timestamp(reader.currentRecordDateTime())

			# if the heading within the line has changed by more than 30 degrees, write out a polygon and continue.
			# gracefully deal with a rotation about 360.0/0.000
			if (180 - abs (abs(hdg - lastheading) -180)) > 30: 
				writepolygon(coveragePoly, left, right, lastTimeStamp, reader.fileName)
				lastheading = hdg
				
		if len(leftside) > window:
			leftside.pop(0)
			rightside.pop(0)
			heading.pop(0)
	# now write out the last part of the polygon to ensure full coverage
	coveragePoly = writepolygon(coveragePoly, left, right, lastTimeStamp, reader.fileName)
	return coveragePoly

###############################################################################
def writepolygon(coveragePoly, left, right, lastTimeStamp, fileName):
	parts = []
	poly = []
	for p in left:
		poly.append(p)

	for p in reversed(right):
		poly.append(p)
		
	parts.append(poly)
	coveragePoly.poly(parts=parts)
	recDate = from_timestamp(lastTimeStamp).strftime("%Y-%m-%d")
	# write out the shape file FIELDS data
	coveragePoly.record(os.path.basename(fileName), recDate) 
	# we have added the record, so now pop everything except the last record
	while len(left) > 1:
		left.pop(0)
		right.pop(0)

	return coveragePoly

###############################################################################
def createSHP(fileName, geometrytype=shapefile.POLYLINE):
	'''open for append or create the shape files. This can be a polyline <false> or polygon '''
	if os.path.isfile(fileName):
		try:
			# Create a shapefile reader
			r = shapefile.Reader(fileName)
			# Create a shapefile writer
			# using the same shape type
			# as our reader
			writer = shapefile.Writer(r.shapeType)
			# Copy over the existing dbf fields
			writer.fields = list(r.fields)
			# Copy over the existing polygons
			writer._shapes.extend(r.shapes())
			# Copy over the existing dbf records
			writer.records.extend(r.records())
		except shapefile.error:
			print ("Problem opening existing shape file, aborting!")
			exit()
	else:
		writer = shapefile.Writer(geometrytype)
		writer.autoBalance = 1
	return writer

###############################################################################
def from_timestamp(unixtime):
	return datetime(1970, 1 ,1) + timedelta(seconds=unixtime)

###############################################################################
def to_timestamp(recordDate):
	return (recordDate - datetime(1970, 1, 1)).total_seconds()

###############################################################################
def createOutputFileName(path):
	'''Create a valid output filename. if the name of the file already exists the file name is auto-incremented.'''
	path	  = os.path.expanduser(path)

	if not os.path.exists(os.path.dirname(path)):
		os.makedirs(os.path.dirname(path))

	if not os.path.exists(path):
		return path

	root, ext = os.path.splitext(os.path.expanduser(path))
	dir	   = os.path.dirname(root)
	fname	 = os.path.basename(root)
	candidate = fname+ext
	index	 = 1
	ls		= set(os.listdir(dir))
	while candidate in ls:
			candidate = "{}_{}{}".format(fname,index,ext)
			index	+= 1

	return os.path.join(dir, candidate)