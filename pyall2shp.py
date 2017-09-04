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
import statistics

def main():

    matches = []

    parser = ArgumentParser(description='Read Kongsberg ALL file and create an ESRI shape file of the trackplot.',
            epilog='Example: \n To convert a single file use -i c:/temp/myfile.all \n to convert all files in a folder use -i c:/temp/*.all\n To convert all .all files recursively in a folder, use -r -i c:/temp \n To convert all .all files recursively from the current folder, use -r -i ./ \n', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', dest='inputFile', action='store', help='-i <ALLfilename> : input ALL filename to image. It can also be a wildcard, e.g. *.all')
    parser.add_argument('-o', dest='outputFile', action='store', default='', help='-o <SHPfilename.shp> : output filename to create. e.g. trackplot.shp [Default: track.shp]')
    parser.add_argument('-s', dest='step', action='store', default='30', help='-s <step size in seconds> : decimate the position datagrams to reduce the shapefile size.  Some systems record at 100Hz.  [Default: 30]')
    parser.add_argument('-c', action='store_true', default=False, dest='coverage', help='-c : create coverage polygon shapefile.  [Default: False]')
    parser.add_argument('-tl', action='store_true', default=False, dest='trackline', help='-tl : create track polyline shapefile.  [Default: False]')
    parser.add_argument('-tp', action='store_true', default=False, dest='trackpoint', help='-tp : create track point shapefile, with runtime information per ping  [Default: False]')
    parser.add_argument('-csv', action='store_true', default=False, dest='csv', help='-cv : create CSV coverage file, with runtime information per ping  [Default: False]')
    parser.add_argument('-r', action='store_true', default=False, dest='recursive', help='-r : search recursively.  [Default: False]')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    # we need to ensure the file is a shp extension

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
        print (matches)

    if len(args.outputFile) == 0:
        trackLineFileName = os.path.join(os.path.dirname(os.path.abspath(matches[0])), "trackline.shp")
        trackLineFileName = createOutputFileName(trackLineFileName)
        trackPointFileName = os.path.join(os.path.dirname(os.path.abspath(matches[0])), "trackpoint.shp")
        trackPointFileName = createOutputFileName(trackPointFileName)
        trackCoverageFileName = os.path.join(os.path.dirname(os.path.abspath(matches[0])), "trackcoverage.shp")
        trackCoverageFileName = createOutputFileName(trackCoverageFileName)
    else:
        trackLineFileName = args.outputFile 
    # if not trackLineFileName.lower().endswith('.shp'):
    #     trackLineFileName += '.shp'
    #     trackPointFileName += '.shp'
    #     trackCoverageFileName += '.shp'

    fileCounter=0
    matches = []
        
    if args.recursive:
        for root, dirnames, filenames in os.walk(os.path.dirname(args.inputFile)):
            for f in fnmatch.filter(filenames, '*.all'):
                matches.append(os.path.join(root, f))
                print (matches[-1])
    else:
        for filename in glob(args.inputFile):
            matches.append(filename)
        print (matches)
    if len(matches) == 0:
        print ("Nothing found to convert, quitting")
        exit()

    # open the output files once only.
    # create the destination shape files 
    if args.trackpoint:
        shp = createSHP(trackPointFileName, shapefile.POINT)
    if args.trackline:
        shp = createSHP(trackLineFileName, shapefile.POLYLINE)
    if args.coverage:
        shp = createSHP(trackCoverageFileName, shapefile.POLYGON)
    if args.coverage:
        shp = createSHP(coverageFileName, True)
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
            createTrackPoint(reader, shp, float(args.step))

        # create the track polyline
        if args.trackline:
            createTrackLine(reader, shp, float(args.step))

        # create the coverage polygon
        if args.coverage:
            createCoverage(reader, shp, float(args.step))

        # create the csv polygon of coverage
        if args.csv:
            createCSV(reader, csv, float(args.step), filename)

        update_progress("Processed: %s (%d/%d)" % (filename, fileCounter, len(matches)), (fileCounter/len(matches)))
        fileCounter +=1
        reader.close()

    update_progress("Process Complete: ", (fileCounter/len(matches)))
    if args.trackpoint:
        print ("Saving track point shapefile: %s" % trackPointFileName)        
        shp.save(trackPointFileName)
        # now write out a prj file so the data has a spatial Reference
        prj = open(trackPointFileName.replace('.shp','.prj'), 'w')
        prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
        prj.close() # you can omit in most cases as the destructor will call it
    if args.trackline:
        print ("Saving track line shapefile: %s" % trackLineFileName)        
        shp.save(trackLineFileName)
        # now write out a prj file so the data has a spatial Reference
        prj = open(trackLineFileName.replace('.shp','.prj'), 'w')
        prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
        prj.close() # you can omit in most cases as the destructor will call it
    if args.coverage:
        print ("Saving coverage polygon shapefile: %s" % coverageFileName)        
        shp.save(coverageFileName)
        # now write out a prj file so the data has a spatial Reference
        prj = open(coverageFileName.replace('.shp','.prj'), 'w')
        prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
        prj.close() # you can omit in most cases as the destructor will call it

def createTrackPoint(reader, shp, step):
    lastTimeStamp = 0
    points = []
    navigation = reader.loadNavigation()
    # remember the previous records so we can compute the speed
    prevX = navigation[0][2]
    prevY = navigation[0][1]
    prevT = navigation[0][0] - 0.001 

    # create the trackpoint shape file
    for update in navigation:
        if update[0] - lastTimeStamp >= step:
            shp.point(float(update[2]),float(update[1]))
            # now add to the shape file.
            recDate = from_timestamp(navigation[0][0]).strftime("%Y%m%d")
            # write out the shape file FIELDS data
            # compute the speed as distance/time
            distance = math.sqrt( ((update[2]-prevX) **2) + ((update[1]-prevY) **2))
            dtime = max(update[0] - prevT, 0.001)
            speed = (distance/dtime) * 60.0 * 3600 # need to convert from degrees to knots
            shp.record(os.path.basename(reader.fileName), int(navigation[0][0]), recDate, speed) 
            
            # remember the last update
            prevX = update[2]
            prevY = update[1]
            prevT = update[0]

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
                depths.append(statistics.median(datagram.Depth))
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
            swath = statistics.median(swathwidth)
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
    leftLatitude = 0;
    leftLongitude = 0
    rightLatitude = 0;
    rightLongitude = 0

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
        if TypeOfDatagram == 'D':
            datagram.read()
            if len(datagram.AcrossTrackDistance) == 0:
                continue
            if (math.fabs(datagram.AcrossTrackDistance[0]) > 0) and (math.fabs(datagram.AcrossTrackDistance[-1]) > 0):
                if (longitude > 0):
                    leftLatitude, leftLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, datagram.Heading - 90 , math.fabs(datagram.AcrossTrackDistance[0]))

                    rightLatitude, rightLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, datagram.Heading + 90, math.fabs(datagram.AcrossTrackDistance[-1]))

        if TypeOfDatagram == 'X':
            datagram.read()
            if len(datagram.AcrossTrackDistance) == 0:
                continue
            if (math.fabs(datagram.AcrossTrackDistance[0]) > 0) and (math.fabs(datagram.AcrossTrackDistance[-1]) > 0) :
                if (longitude > 0):
                    leftLatitude, leftLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, datagram.Heading - 90, math.fabs(datagram.AcrossTrackDistance[0]))

                    rightLatitude, rightLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, datagram.Heading + 90,  math.fabs(datagram.AcrossTrackDistance[-1]))

        # add to the shape file at the user required interval
        if to_timestamp(reader.currentRecordDateTime()) - lastTimeStamp >= step:
            if (leftLongitude > 0) and (longitude > 0):
                left.append([leftLongitude,leftLatitude])        
                right.insert(0,[rightLongitude,rightLatitude])        
                lastTimeStamp = to_timestamp(reader.currentRecordDateTime())
    
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
            # Copy over the existing dbf records
            writer.records.extend(r.records())
            # Copy over the existing polygons
            writer._shapes.extend(r.shapes())
        except shapefile.error:
            print ("Problem opening existing shape file, aborting!")
            exit()
    else:
        writer = shapefile.Writer(geometrytype)
        writer.autoBalance = 1
        writer.field("LineName", "C")
        # w.field("WaterDepth", "N")
        writer.field("UNIXTime", "N")
        writer.field("SurveyDate", "D")
        writer.field("SpeedKnots", "N")
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
     path      = os.path.expanduser(path)

     if not os.path.exists(os.path.dirname(path)):
         os.makedirs(os.path.dirname(path))

     if not os.path.exists(path):
        return path

     root, ext = os.path.splitext(os.path.expanduser(path))
     dir       = os.path.dirname(root)
     fname     = os.path.basename(root)
     candidate = fname+ext
     index     = 1
     ls        = set(os.listdir(dir))
     while candidate in ls:
             candidate = "{}_{}{}".format(fname,index,ext)
             index    += 1

     return os.path.join(dir, candidate)


if __name__ == "__main__":
    main()

