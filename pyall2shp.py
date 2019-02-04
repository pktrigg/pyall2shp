import sys
import time
import os
import tempfile
import fnmatch
import math
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from datetime import datetime
from datetime import timedelta
from glob import glob
import multiprocessing as mp
import pyproj
import logging

# local imports
# import pyall

# local from the shared area...
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'shared'))
import fileutils
import geodetic
import SSDM
import shapefile

# we need to do this as airflow and regular python cmdline interpreter differ.
localpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(localpath)

# import shapefile
import pyall


# sys.path.append("../shared")
# import Common
# from ..Shared import SSDM
##############################################################################
# def main():
def main(*opargs, **kwargs):
	logging.info('Running PYALL2SHP')

	parser = ArgumentParser(description='Read Kongsberg ALL file and create an ESRI shape file of the trackplot.',
			epilog='Example: \n To process a single file use -i c:/temp/myfile.all \n to mass process every file in a folder use -i c:/temp/*.all\n To convert all .all files recursively in a folder, use -r -i c:/temp \n To convert all .all files recursively from the current folder, use -r -i ./*.all \n', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', 		dest='inputFile', 	action='store', 		default='',		help='input ALL filename to image. It can also be a wildcard, e.g. *.all')
	parser.add_argument('-o', 		dest='outputFile', 	action='store', 		default='', 	help='output filename to create. e.g. coverage.shp [Default: coverage.shp]')
	parser.add_argument('-s', 		dest='step', 		action='store', 		default='30', 	help='decimate the data to reduce the output size. [Default: 30]')
	parser.add_argument('-ta', 		dest='trackall', 	action='store_true', 	default=False, 	help='create trackline, point and coverage polygon shapefiles.')
	parser.add_argument('-tc', 		dest='trackcoverage',action='store_true', 	default=False,  help='create coverage polygon shapefile.')
	parser.add_argument('-tl', 		dest='trackline',	action='store_true', 	default=False,  help='create track polyline shapefile.')
	parser.add_argument('-tp', 		dest='trackpoint',	action='store_true', 	default=False,  help='create track point shapefile, with runtime information per ping.')
	parser.add_argument('-odir', 	dest='odir', 		action='store', 		default="", 	help='Specify a relative output folder e.g. -odir GIS')
	parser.add_argument('-opath', 	dest='opath', 		action='store', 		default="", 	help='Specify an output path e.g. -opath c:/temp')
	parser.add_argument('-odix', 	dest='odix', 		action='store', 		default="", 	help='Specify an output filename appendage e.g. -odix _coverage')
	parser.add_argument('-r', 		dest='recursive',	action='store_true', 	default=False,  help='search recursively.')
	parser.add_argument('-epsg', 	dest='epsg', 		action='store', 		default="0", 	help='Specify an output EPSG code for transforming from WGS84 to East,North,e.g. -epsg 4326')
	parser.add_argument('-cores', 	dest='cores', 		action='store', 		default="0", 	help='Specify the number of cores to use for processing. By default all cores shall be used, e.g. -cores 4')

	# logging.info("OPARGS",opargs)

	# need to handle args from command line AND from Airflow.  they come from different streams
	if len(opargs) > 0:
		args = parser.parse_args(opargs)
	else:
		args = parser.parse_args()

	# logging.info("XXX",args)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)

	# if len(sys.argv)==1:
	# 	parser.print_help()
	# 	sys.exit(1)

	# args = parser.parse_args()
	process(args)

###############################################################################
def process(args):

	matches = fileutils.findFiles(args.recursive, args.inputFile, "*.all")

	if len(args.outputFile) == 0:
		fname, ext = os.path.splitext(os.path.expanduser(matches[0]))
		args.outputFile = fname
		# args.outputFile = "Track"
	if len(args.opath) == 0:
		args.opath = os.path.dirname(os.path.abspath(args.outputFile))

	# trackLineFileName = os.path.join(os.path.dirname(os.path.abspath(args.outputFile)), fname + "_MBESLine.shp")
	trackLineFileName = os.path.join(args.opath, args.odir, os.path.basename(args.outputFile) + "Survey_TrackLines.shp")
	trackLineFileName  = addFileNameAppendage(trackLineFileName, args.odix)
	trackLineFileName = fileutils.createOutputFileName(trackLineFileName)

	# trackPointFileName = os.path.join(os.path.dirname(os.path.abspath(args.outputFile)), fname + "_MBESPoint.shp")
	trackPointFileName = os.path.join(args.opath, args.odir, os.path.basename(args.outputFile) + "Survey_TrackPoint.shp")
	trackPointFileName  = addFileNameAppendage(trackPointFileName, args.odix)
	trackPointFileName = fileutils.createOutputFileName(trackPointFileName)

	# trackCoverageFileName = os.path.join(os.path.dirname(os.path.abspath(args.outputFile)), fname + "_trackCoverage.shp")
	trackCoverageFileName = os.path.join(args.opath, args.odir, os.path.basename(args.outputFile) + "Survey_TrackCoverage.shp")
	trackCoverageFileName  = addFileNameAppendage(trackCoverageFileName, args.odix)
	trackCoverageFileName = fileutils.createOutputFileName(trackCoverageFileName)

	#load the python proj projection object library if the user has requested it
	if int(args.epsg) == 4326:
		args.epsg = "0"
	if len(args.epsg) > 0:
		projection = geodetic.loadProj(args.epsg)
	else:
		projection = None

	# open the output files once only.
	# create the destination shape files
	TPshp = None
	TLshp = None
	TCshp = None

	if args.trackall:
		args.trackpoint = True
		args.trackline = True
		args.trackcoverage = True


	if args.trackpoint:
		TPshp = SSDM.createPointShapeFile(trackPointFileName)
	if args.trackline:
		TLshp = SSDM.createSurveyTracklineSHP(trackLineFileName)
	if args.trackcoverage:
		TCshp = SSDM.createCoverageSHP(trackCoverageFileName)

	MULTICORE = True

	if MULTICORE:
		# # for multiprocessing we call the process tool based on the number of cores.
		if int(args.cores) > 0:
			cores = int(args.cores)
		else:
			cores = mp.cpu_count()
		pool = mp.Pool(processes=int(cores))

		if args.trackline:
			results = [pool.apply(processTrackLine, args=(filename, float(args.step), trackLineFileName, projection)) for filename in matches]
			#merge all the results
			for r in results:
				for shape in r._shapes:
					TLshp._shapes.append(shape)
				# now add the fields
				for fields in r.records:
					TLshp.records.append(fields)

		if args.trackcoverage:
			results = [pool.apply(processTrackCoverage, args=(filename, float(args.step), trackCoverageFileName, projection)) for filename in matches]
			#merge all the results
			for r in results:
				for shape in r._shapes:
					TCshp._shapes.append(shape)
				# now add the fields
				for fields in r.records:
					TCshp.records.append(fields)

		if args.trackpoint:
			results = [pool.apply(processTrackPoint, args=(filename, float(args.step), trackPointFileName, projection)) for filename in matches]
			#merge all the results
			for r in results:
				for shape in r._shapes:
					TPshp._shapes.append(shape)
				# now add the fields
				for fields in r.records:
					TPshp.records.append(fields)
	else:
		for filename in matches:
			reader = pyall.ALLReader(filename)
			if args.trackpoint:
				# TPshp = processTrackPoint(filename, float(args.step),trackPointFileName, projection)
				createTrackPoint(reader, TPshp, float(args.step), projection)
			if args.trackline:
				# TLshp = processTrackLine(filename, float(args.step), trackLineFileName, projection)
				totalDistanceRun = createTrackLine(reader, TLshp, float(args.step), projection)
				# print ("%s Trackplot Length: %.3f" % (filename, totalDistanceRun))
			if args.trackcoverage:
				# TCshp = processTrackCoverage(filename, float(args.step), trackCoverageFileName, projection)
				createTrackCoverage(reader, TCshp, float(args.step), projection)

	# 	update_progress("Processed: %s (%d/%d)" % (filename, fileCounter, len(matches)), (fileCounter/len(matches)))
	# 	fileCounter +=1

	#now we can write out the results to a shape file...
	# update_progress("Process Complete: ", (fileCounter/len(matches)))
	if args.trackpoint:
		print ("Saving track point: %s" % trackPointFileName)
		TPshp.save(trackPointFileName)
		# now write out a prj file so the data has a spatial Reference
		filename = trackPointFileName.replace('.shp','.prj')
		geodetic.writePRJ(filename, args.epsg)
	if args.trackline:
		print ("Saving track line: %s" % trackLineFileName)
		TLshp.save(trackLineFileName)
		# now write out a prj file so the data has a spatial Reference
		filename = trackLineFileName.replace('.shp','.prj')
		geodetic.writePRJ(filename, args.epsg)
	if args.trackcoverage:
		print ("Saving coverage polygon: %s" % trackCoverageFileName)
		TCshp.save(trackCoverageFileName)
		# now write out a prj file so the data has a spatial Reference
		filename = trackCoverageFileName.replace('.shp','.prj')
		geodetic.writePRJ(filename, args.epsg)

def processTrackPoint(filename, step, trackPointFileName, projection):
	# dirname, basename = os.path.split(filename)
	# trackPointFileName = tempfile.NamedTemporaryFile(prefix=basename, dir=dirname)
	shp = SSDM.createPointShapeFile(trackPointFileName)
	reader = pyall.ALLReader(filename)
	# create the track point with a point recoard and metadata per ping
	createTrackPoint(reader, shp, step, projection)
	print ("Trackpoint created for: %s" % (filename))
	reader.close()
	return shp

def processTrackLine(filename, step, trackLineFileName, projection):
	shp = SSDM.createSurveyTracklineSHP(trackLineFileName)
	reader = pyall.ALLReader(filename)
	# create the track polyline
	totalDistanceRun = createTrackLine(reader, shp, step, projection)
	print ("Trackplot created for: %s, Length: %.3f" % (filename, totalDistanceRun))
	reader.close()
	return shp

def processTrackCoverage(filename, step, trackCoverageFileName, projection):
	shp = SSDM.createCoverageSHP(trackCoverageFileName)
	reader = pyall.ALLReader(filename)
	# create the coverage polygon
	createTrackCoverage(reader, shp, step, projection)
	print ("Trackcoverage created for: %s" % (filename))
	reader.close()
	return shp

###############################################################################
def createTrackPoint(reader, shp, step, projection):
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
	# arcpy.AddMessage("Creating Track Points for:" + reader.fileName)

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
			if projection is not None:
				x,y = projection(float(longitude),float(latitude))
				shp.point(x,y)
			else:
				shp.point(longitude,latitude)
			# now add to the shape file.
			recDate = from_timestamp(recTime).strftime("%Y%m%d")
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
def createTrackLine(reader, trackLine, step, projection):
	lastTimeStamp = 0
	line_parts = []
	line = []
	navigation = reader.loadNavigation()
	totalDistanceRun = 0
	if len(navigation) == 0: #trap out empty .all files.
		return totalDistanceRun
	prevX =  navigation[0][2]
	prevY = navigation[0][1]

	for update in navigation:
		longitude = update[2]
		latitude = update[1]
		distance = geodetic.est_dist(latitude, longitude, prevY, prevX)
		# distance = math.sqrt( ((longitude - prevX) **2) + ((latitude - prevY) **2))
		totalDistanceRun += distance
		prevX = longitude
		prevY = latitude

	# create the trackline shape file
	for update in navigation:
		if update[0] - lastTimeStamp >= step:
			if projection is not None:
				x,y = projection(float(update[2]),float(update[1]))
				line.append([x,y])
			else:
				line.append([float(update[2]),float(update[1])])
			lastTimeStamp = update[0]
	# now add the very last update
	if projection is not None:
		x,y = projection(float(navigation[-1][2]),float(navigation[-1][1]))
		line.append([x,y])
	else:
		line.append([float(navigation[-1][2]),float(navigation[-1][1])])
	# line.append([float(navigation[-1][2]),float(navigation[-1][1])])
	# print("Points added to track: %d" % (len(line)))
	line_parts.append(line)
	trackLine.line(parts=line_parts)
	# now add to the shape file.
	recDate = from_timestamp(navigation[0][0]).strftime("%Y%m%d")
	# write out the shape file FIELDS data
	# preparedDate = datetime.now()
	userName = os.getenv('username')
	filename = os.path.basename(reader.fileName)
	trackLine.record(LAST_UPDATE = recDate,
		LAST_UPDATE_BY=userName,
		FEATURE_ID=0,
		SURVEY_ID=0,
		SURVEY_ID_REF=filename,
		REMARKS=reader.fileName,
		LINE_ID=0,
		LINE_NAME=filename[:40],
		LINE_DIRECTION=1.123,
		SYMBOLOGY_CODE=0,
		LAST_SEIS_PT_ID=0,
		DATA_SOURCE=filename,
		CONTRACTOR_NAME= userName,
		LINE_LENGTH = distance,
		FIRST_SEIS_PT_ID= 0,
		HIRES_SEISMIC_EQL_URL= "",
		OTHER_DATA_URL = "",
		LAYER= "",
		SHAPE_Length= distance
		)
	return totalDistanceRun

###############################################################################
def createTrackCoverage(reader, coveragePoly, step, projection):
	lastTimeStamp = 0
	lastheading = None # remember the heading so we can chop up into smaller polygons.  This is important if the vessel is running in circles rather than lines.
	left = []
	right = []
	selectedPositioningSystem = None
	latitude = 0
	longitude = 0
	heading = []
	leftside = [] #sliding window
	rightside = [] #sliding window
	window = step #sliding window size in number of pings, to smooth the data
	pendingrecord = False

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
				pendingrecord = True #performance upgrade. only make the coverage at the same rate as the navigation. no oint any faster

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

			if projection is not None:
				x,y = projection(float(leftLongitude),float(leftLatitude))
				left.append([x,y])
				x,y = projection(float(rightLongitude),float(rightLatitude))
				right.append([x,y])
			else:
				left.append([leftLongitude,leftLatitude])
				right.append([rightLongitude,rightLatitude])

			lastTimeStamp = to_timestamp(reader.currentRecordDateTime())

			# if the heading within the line has changed by more than 30 degrees, write out a polygon and continue.
			# if abs (hdg - lastheading) > 30:
			if (180 - abs (abs(hdg - lastheading) -180)) > 30:
				coveragePoly = writepolygon(coveragePoly, left, right, lastTimeStamp, reader.fileName)
				lastheading = hdg

		if len(leftside) > window:
			leftside.pop(0)
			rightside.pop(0)
			heading.pop(0)
	# now write out the last part of the polygon to ensure full coverage
	leftLatitude, leftLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, hdg - 90 , math.fabs(leftextent))
	rightLatitude, rightLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, hdg + 90, math.fabs(rightextent))
	if projection is not None:
		x,y = projection(float(leftLongitude),float(leftLatitude))
		left.append([x,y])
		x,y = projection(float(rightLongitude),float(rightLatitude))
		right.append([x,y])
	else:
		left.append([leftLongitude,leftLatitude])
		right.append([rightLongitude,rightLatitude])
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
	recDate = from_timestamp(lastTimeStamp).strftime("%Y%m%d")
	# write out the shape file FIELDS data

	coveragePoly.record(os.path.basename(fileName), recDate)
	# we have added the record, so now pop everything except the last record
	while len(left) > 1:
		left.pop(0)
		right.pop(0)

	return coveragePoly


###############################################################################
def from_timestamp(unixtime):
	return datetime(1970, 1 ,1) + timedelta(seconds=unixtime)

###############################################################################
def to_timestamp(recordDate):
	return (recordDate - datetime(1970, 1, 1)).total_seconds()

###############################################################################
def update_progress(job_title, progress):
	length = 20 # modify this to change the length
	block = int(round(length*progress))
	msg = "\r{0}: [{1}] {2}%".format(job_title, "#"*block + "-"*(length-block), round(progress*100, 2))
	if progress >= 1: msg += " DONE\r\n"
	sys.stdout.write(msg)
	sys.stdout.flush()

###############################################################################
def addFileNameAppendage(path, appendage):
	'''Create a valid output filename. if the name of the file already exists the file name is auto-incremented.'''
	path = os.path.expanduser(path)

	if not os.path.exists(os.path.dirname(path)):
		os.makedirs(os.path.dirname(path))

	# if not os.path.exists(path):
	# 	return path

	root, ext = os.path.splitext(os.path.expanduser(path))
	dir	   = os.path.dirname(root)
	fname	 = os.path.basename(root)
	candidate = "{}{}{}".format(fname, appendage, ext)

	return os.path.join(dir, candidate)

# ###############################################################################
# def createCSV(reader, csv, step, filename):
# 	lastTimeStamp = 0
# 	swathwidth = [] #a list of swath widths so we can smooth
# 	depths = [] # a list of nadir depths so we can smooth a little
# 	ping = 0
# 	left = 0
# 	right = 0
# 	window = 100 #sliding window size in number of pings, to smooth the data
# 	maximumPortCoverageDegrees = 0
# 	maximumPortWidth = 0
# 	maximumStbdCoverageDegrees = 0
# 	maximumStbdWidth = 0
# 	mode = ''
# 	acrosstrackresolution = 0
# 	alongtrackresolution = 0
# 	nbeams = 0
# 	speed = 0
# 	pingtimes = []
# 	prevX = 0
# 	prevY = 0
# 	prevT = 0

# 	reader.rewind()
# 	while reader.moreData():
# 		TypeOfDatagram, datagram = reader.readDatagram()
# 		if TypeOfDatagram == 'R':
# 			datagram.read()
# 			maximumPortCoverageDegrees = datagram.maximumPortCoverageDegrees
# 			maximumPortWidth = datagram.maximumPortWidth
# 			maximumStbdCoverageDegrees = datagram.maximumStbdCoverageDegrees
# 			maximumStbdWidth = datagram.maximumStbdWidth
# 			depthmode = datagram.DepthMode

# 		if TypeOfDatagram == 'P':
# 			datagram.read()
# 			speed = datagram.SpeedOverGround
# 			distance = math.sqrt( ((datagram.Longitude - prevX) **2) + ((datagram.Latitude - prevY) **2))
# 			dtime = max(to_timestamp(reader.currentRecordDateTime()) - prevT, 0.001)
# 			speed = (distance/dtime) * 60.0 * 3600 # need to convert from degrees to knots
# 			prevX = datagram.Longitude
# 			prevY = datagram.Latitude
# 			prevT = to_timestamp(reader.currentRecordDateTime())

# 		if TypeOfDatagram == 'D' or TypeOfDatagram == 'X':
# 			datagram.read()
# 			if len(datagram.AcrossTrackDistance) == 0:
# 				continue
# 			if len(datagram.Depth) == 0:
# 				continue
# 			nadirbeamno = math.floor(len(datagram.Depth)/2)
# 			if (math.fabs(datagram.AcrossTrackDistance[0]) > 0) and (math.fabs(datagram.AcrossTrackDistance[-1]) > 0):
# 				left = min(datagram.AcrossTrackDistance)
# 				right = max(datagram.AcrossTrackDistance)
# 				swathwidth.append(math.fabs(left) + math.fabs(right))
# 				depths.append(median(datagram.Depth))
# 				# depths.append(statistics.median(datagram.Depth))
# 				# depths.append(datagram.Depth[nadirbeamno])
# 				pingtimes.append(to_timestamp(reader.currentRecordDateTime()))
# 				ping = datagram.Counter
# 				nbeams = datagram.NBeams

# 		# add to the shape file at the user required interval
# 		if to_timestamp(reader.currentRecordDateTime()) - lastTimeStamp >= step:
# 			if len(depths) == 0 or len(swathwidth) == 0:
# 				continue
# 			# we use the maximum depth as the spikes are not real and are typically always shallower
# 			# we use the median swath width to get a representative swath
# 			swath = median(swathwidth)
# 			# swath = statistics.median(swathwidth)
# 			acrosstrackresolution = swath / nbeams
# 			duration = ( pingtimes[-1] - pingtimes[0] ) / len(pingtimes) #compute the time between pings over the moving window so we see some stability
# 			alongtrackresolution = (speed * (1852/3600)) * duration #convert speed to metres/second
# 			csv.write("%s,%s,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%s,%.3f,%.3f,%.3f\n" % (filename, str(reader.currentRecordDateTime()), ping, max(depths), swath, maximumPortWidth, maximumStbdWidth, maximumPortCoverageDegrees, maximumStbdCoverageDegrees, depthmode, speed, acrosstrackresolution, alongtrackresolution))
# 			lastTimeStamp = to_timestamp(reader.currentRecordDateTime())

# 		if len(depths) > window:
# 			depths.pop(0)
# 			swathwidth.pop(0)
# 			pingtimes.pop(0)
# 	return

if __name__ == "__main__":
	main()
