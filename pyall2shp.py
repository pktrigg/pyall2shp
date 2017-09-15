import sys
import time
import os
import fnmatch
import math
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from datetime import datetime
from datetime import timedelta
from glob import glob

# local imports
import pyall
import shapefile
import geodetic
# import statistics

##############################################################################
def main():
	parser = ArgumentParser(description='Read Kongsberg ALL file and create an ESRI shape file of the trackplot.',
			epilog='Example: \n To process a single file use -i c:/temp/myfile.all \n to mass process every file in a folder use -i c:/temp/*.all\n To convert all .all files recursively in a folder, use -r -i c:/temp \n To convert all .all files recursively from the current folder, use -r -i ./*.all \n', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', dest='inputFile', action='store', help='-i <ALLfilename> : input ALL filename to image. It can also be a wildcard, e.g. *.all')
	parser.add_argument('-o', dest='outputFile', action='store', default='coverage.shp', help='-o <SHPfilename.shp> : output filename to create. e.g. coverage.shp [Default: coverage.shp]')
	parser.add_argument('-s', dest='step', action='store', default='30', help='-s <step size in seconds> : decimate the data to reduce the output size. [Default: 30]')
	parser.add_argument('-c', action='store_true', default=False, dest='coverage', help='-c : create coverage polygon shapefile.')
	parser.add_argument('-tl', action='store_true', default=False, dest='trackline', help='-tl : create track polyline shapefile.')
	parser.add_argument('-tp', action='store_true', default=False, dest='trackpoint', help='-tp : create track point shapefile, with runtime information per ping.')
	parser.add_argument('-csv', action='store_true', default=False, dest='csv', help='-cv : create CSV coverage file, with runtime information per ping.')
	parser.add_argument('-r', action='store_true', default=False, dest='recursive', help='-r : search recursively.')
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
		
	args = parser.parse_args()
	process(args)

###############################################################################
def process(args):
	matches = []
	if args.recursive:
		for root, dirnames, filenames in os.walk(os.path.dirname(args.inputFile)):
			for f in fnmatch.filter(filenames, '*.all'):
				matches.append(os.path.join(root, f))
				print (matches[-1])
	else:
		if os.path.exists(args.inputFile):
			matches.append (os.path.abspath(args.inputFile))
		else:
			for filename in glob(args.inputFile):
				matches.append(filename)
	if len(matches) == 0:
		print ("Nothing found to convert, quitting")
		exit()
	print (matches)

	fname, ext = os.path.splitext(os.path.expanduser(args.outputFile))
	trackLineFileName = os.path.join(os.path.dirname(os.path.abspath(args.outputFile)), fname + "_trackLine.shp")
	trackPointFileName = os.path.join(os.path.dirname(os.path.abspath(args.outputFile)), fname + "_trackPoint.shp")
	trackCoverageFileName = os.path.join(os.path.dirname(os.path.abspath(args.outputFile)), fname + "_trackCoverage.shp")

	# open the output files once only.
	# create the destination shape files 
	fileCounter=0
	if args.trackpoint:
		TPshp = createSHP(trackPointFileName, shapefile.POINT)
		if len(TPshp.fields) <= 1: #there is a default deletion flag field always set, to we need to accoutn for this.
			TPshp.field("LineName", "C")
			TPshp.field("SurveyDate", "D")
			TPshp.field("SurveyTime", "C")
			TPshp.field("UNIXTime", "N")
			TPshp.field("SpeedKnots", "N")
			TPshp.field("PortCover", "N")
			TPshp.field("StbdCover", "N")
			TPshp.field("PortWidth", "N")
			TPshp.field("StbdWidth", "N")
			TPshp.field("DepthMode", "C")
			TPshp.field("Absorption", "N")
			TPshp.field("PulseLength", "N")
			TPshp.field("TVG", "N")
			TPshp.field("DualSwath", "C")
			TPshp.field("SpikeFilt", "C")
			TPshp.field("Stabilise", "C")
			TPshp.field("MinZGate", "N")
			TPshp.field("MaxZGate", "N")
			TPshp.field("BeamSpace", "C")

	if args.trackline:
		TLshp = createSHP(trackLineFileName, shapefile.POLYLINE)
		if len(TLshp.fields) <= 1:
			TLshp.field("LineName", "C")
			TLshp.field("SurveyDate", "D")
	if args.coverage:
		TCshp = createSHP(trackCoverageFileName, shapefile.POLYGON)
		if len(TCshp.fields) <= 1:
			TCshp.field("LineName", "C")
			TCshp.field("SurveyDate", "D")
	if args.csv:
		csv = open(args.outputFile, 'a')
		csv.write("Filename,Date,Ping,Depth(m),Swathwidth(m),maximumPortWidth(m),maximumStbdWidth(m),maximumPortCoverageDegrees,maximumStbdCoverageDegrees,DepthMode,Speed(Kts),AcrossTrackResolution(m),AlongTrackResolution(m)\n")

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

		update_progress("Processed: %s (%d/%d)" % (filename, fileCounter, len(matches)), (fileCounter/len(matches)))
		fileCounter +=1
		reader.close()

	update_progress("Process Complete: ", (fileCounter/len(matches)))
	if args.trackpoint:
		print ("Saving track point shapefile: %s" % trackPointFileName)		
		TPshp.save(trackPointFileName)
		# now write out a prj file so the data has a spatial Reference
		prj = open(trackPointFileName.replace('.shp','.prj'), 'w')
		prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
		prj.close() # you can omit in most cases as the destructor will call it
	if args.trackline:
		print ("Saving track line shapefile: %s" % trackLineFileName)		
		TLshp.save(trackLineFileName)
		# now write out a prj file so the data has a spatial Reference
		prj = open(trackLineFileName.replace('.shp','.prj'), 'w')
		prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
		prj.close() # you can omit in most cases as the destructor will call it
	if args.coverage:
		print ("Saving coverage polygon shapefile: %s" % trackCoverageFileName)		
		TCshp.save(trackCoverageFileName)
		# now write out a prj file so the data has a spatial Reference
		prj = open(trackCoverageFileName.replace('.shp','.prj'), 'w')
		prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
		prj.close() # you can omit in most cases as the destructor will call it

###############################################################################
def createTrackPoint(reader, shp, step):
	lastTimeStamp = 0
	recTime = 0
	recTimeString = ""
	latitude = 0
	longitude = 0
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
	# remember the previous records so we can compute the speed
	prevX = 0
	prevY = 0
	prevT = 0
	# arcpy.AddMessage("Creating Track Points for:" + reader.fileName)

	reader.rewind() #rewind to the start of the file
	while reader.moreData():
		TypeOfDatagram, datagram = reader.readDatagram()
		if (TypeOfDatagram == 'P'):
			datagram.read()
			if (selectedPositioningSystem == None):
				selectedPositioningSystem = datagram.Descriptor
			if (selectedPositioningSystem == datagram.Descriptor):
				latitude = datagram.Latitude
				longitude = datagram.Longitude
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
				speed, 
				int(maximumPortCoverageDegrees), 
				int(maximumStbdCoverageDegrees), 
				int(maximumPortWidth), 
				int(maximumStbdWidth), 
				depthmode, 
				absorptioncoefficient, 
				pulselength, 
				tvg, 
				dualswath, 
				spikefilter, 
				stabilisation, 
				mindepthgate, 
				maxdepthgate, 
				beamspacing)
			
			# remember the last update
			prevX = longitude
			prevY = latitude
			prevT = lastTimeStamp


# def createTrackPoint(reader, shp, step):
# 	lastTimeStamp = 0
# 	points = []
# 	navigation = reader.loadNavigation()
# 	# remember the previous records so we can compute the speed
# 	prevX = navigation[0][2]
# 	prevY = navigation[0][1]
# 	prevT = navigation[0][0] - 0.001 

# 	# create the trackpoint shape file
# 	for update in navigation:
# 		if update[0] - lastTimeStamp >= step:
# 			shp.point(float(update[2]),float(update[1]))
# 			# now add to the shape file.
# 			recDate = from_timestamp(navigation[0][0]).strftime("%Y%m%d")
# 			# write out the shape file FIELDS data
# 			# compute the speed as distance/time
# 			distance = math.sqrt( ((update[2]-prevX) **2) + ((update[1]-prevY) **2))
# 			dtime = max(update[0] - prevT, 0.001)
# 			speed = int((distance/dtime) * 60.0 * 3600 * 100) # need to convert from degrees to centi knots
# 			shp.record(os.path.basename(reader.fileName), int(navigation[0][0]), recDate, speed) 
# 			lastTimeStamp = update[0]
			
# 			# remember the last update
# 			prevX = update[2]
# 			prevY = update[1]
# 			prevT = update[0]

def median(lst):
	n = len(lst)
	if n < 1:
			return None
	if n % 2 == 1:
			return sorted(lst)[n//2]
	else:
			return sum(sorted(lst)[n//2-1:n//2+1])/2.0

def createTrackLine(reader, trackLine, step):
	lastTimeStamp = 0
	line_parts = []
	line = []
	navigation = reader.loadNavigation()

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
	recDate = from_timestamp(navigation[0][0]).strftime("%Y%m%d")
	# write out the shape file FIELDS data
	speed = 0
	trackLine.record(os.path.basename(reader.fileName), int(navigation[0][0]), recDate, speed) 

###############################################################################
def createCSV(reader, csv, step, filename):
	lastTimeStamp = 0
	swathwidth = [] #a list of swath widths so we can smooth
	depths = [] # a list of nadir depths so we can smooth a little
	ping = 0
	left = 0
	right = 0
	window = 100 #sliding window size in number of pings, to smooth the data
	maximumPortCoverageDegrees = 0
	maximumPortWidth = 0
	maximumStbdCoverageDegrees = 0
	maximumStbdWidth = 0
	mode = ''
	acrosstrackresolution = 0
	alongtrackresolution = 0
	nbeams = 0
	speed = 0
	pingtimes = []
	prevX = 0
	prevY = 0
	prevT = 0

	reader.rewind()
	while reader.moreData():
		TypeOfDatagram, datagram = reader.readDatagram()
		if TypeOfDatagram == 'R':
			datagram.read()
			maximumPortCoverageDegrees = datagram.maximumPortCoverageDegrees
			maximumPortWidth = datagram.maximumPortWidth
			maximumStbdCoverageDegrees = datagram.maximumStbdCoverageDegrees
			maximumStbdWidth = datagram.maximumStbdWidth
			depthmode = datagram.DepthMode

		if TypeOfDatagram == 'P':
			datagram.read()
			speed = datagram.SpeedOverGround
			distance = math.sqrt( ((datagram.Longitude - prevX) **2) + ((datagram.Latitude - prevY) **2))
			dtime = max(to_timestamp(reader.currentRecordDateTime()) - prevT, 0.001)
			speed = (distance/dtime) * 60.0 * 3600 # need to convert from degrees to knots
			prevX = datagram.Longitude
			prevY = datagram.Latitude
			prevT = to_timestamp(reader.currentRecordDateTime())
		
		if TypeOfDatagram == 'D' or TypeOfDatagram == 'X':
			datagram.read()
			if len(datagram.AcrossTrackDistance) == 0:
				continue
			if len(datagram.Depth) == 0:
				continue
			nadirbeamno = math.floor(len(datagram.Depth)/2)
			if (math.fabs(datagram.AcrossTrackDistance[0]) > 0) and (math.fabs(datagram.AcrossTrackDistance[-1]) > 0):
				left = min(datagram.AcrossTrackDistance) 
				right = max(datagram.AcrossTrackDistance)
				swathwidth.append(math.fabs(left) + math.fabs(right))
				depths.append(median(datagram.Depth))
				# depths.append(statistics.median(datagram.Depth))
				# depths.append(datagram.Depth[nadirbeamno])
				pingtimes.append(to_timestamp(reader.currentRecordDateTime()))
				ping = datagram.Counter
				nbeams = datagram.NBeams

		# add to the shape file at the user required interval
		if to_timestamp(reader.currentRecordDateTime()) - lastTimeStamp >= step:
			if len(depths) == 0 or len(swathwidth) == 0:
				continue
			# we use the maximum depth as the spikes are not real and are typically always shallower
			# we use the median swath width to get a representative swath
			swath = median(swathwidth)
			# swath = statistics.median(swathwidth)
			acrosstrackresolution = swath / nbeams
			duration = ( pingtimes[-1] - pingtimes[0] ) / len(pingtimes) #compute the time between pings over the moving window so we see some stability
			alongtrackresolution = (speed * (1852/3600)) * duration #convert speed to metres/second
			csv.write("%s,%s,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%s,%.3f,%.3f,%.3f\n" % (filename, str(reader.currentRecordDateTime()), ping, max(depths), swath, maximumPortWidth, maximumStbdWidth, maximumPortCoverageDegrees, maximumStbdCoverageDegrees, depthmode, speed, acrosstrackresolution, alongtrackresolution))
			lastTimeStamp = to_timestamp(reader.currentRecordDateTime())

		if len(depths) > window:
			depths.pop(0)
			swathwidth.pop(0)
			pingtimes.pop(0)
	return
	
def createCoverage(reader, coveragePoly, step):
	lastTimeStamp = 0
	parts = []
	left = []
	right = []
	selectedPositioningSystem = None
	latitude = 0;
	longitude = 0
	# leftLatitude = 0;
	# leftLongitude = 0
	# rightLatitude = 0;
	# rightLongitude = 0

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

			# datagram.read()
			# if len(datagram.AcrossTrackDistance) == 0:
			#	 continue
			# if (math.fabs(datagram.AcrossTrackDistance[0]) > 0) and (math.fabs(datagram.AcrossTrackDistance[-1]) > 0):
			#	 if (longitude > 0):
					# leftLatitude, leftLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, datagram.Heading - 90 , math.fabs(datagram.AcrossTrackDistance[0]))
					# rightLatitude, rightLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, datagram.Heading + 90, math.fabs(datagram.AcrossTrackDistance[-1]))

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
			right.insert(0,[rightLongitude,rightLatitude])		
			lastTimeStamp = to_timestamp(reader.currentRecordDateTime())

		if len(leftside) > window:
			leftside.pop(0)
			rightside.pop(0)
			heading.pop(0)
	
	poly = []
	for p in left:
		poly.append(p)
	for p in right:
		poly.append(p)
	parts.append(poly)
	coveragePoly.poly(parts=parts)
	recDate = from_timestamp(lastTimeStamp).strftime("%Y%m%d")
	# write out the shape file FIELDS data
	speed = 0
	coveragePoly.record(os.path.basename(reader.fileName), int(lastTimeStamp), recDate, speed) 
		
	return coveragePoly

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
		# writer.field("LineName", "C")
		# # w.field("WaterDepth", "N")
		# writer.field("UNIXTime", "N")
		# writer.field("SurveyDate", "D")
		# writer.field("SpeedCentiKnots", "N")
	return writer

def from_timestamp(unixtime):
	return datetime(1970, 1 ,1) + timedelta(seconds=unixtime)

def to_timestamp(recordDate):
	return (recordDate - datetime(1970, 1, 1)).total_seconds()

def update_progress(job_title, progress):
	length = 20 # modify this to change the length
	block = int(round(length*progress))
	msg = "\r{0}: [{1}] {2}%".format(job_title, "#"*block + "-"*(length-block), round(progress*100, 2))
	if progress >= 1: msg += " DONE\r\n"
	sys.stdout.write(msg)
	sys.stdout.flush()

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

if __name__ == "__main__":
	main()

