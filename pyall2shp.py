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

def main():
    parser = ArgumentParser(description='Read Kongsberg ALL file and create an ESRI shape file of the trackplot.',
            epilog='Example: \n To convert a single file use -i c:/temp/myfile.all \n to convert all files in a folder use -i c:/temp/*.all\n To convert all .all files recursively in a folder, use -r -i c:/temp \n To convert all .all files recursively from the current folder, use -r -i ./ \n', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', dest='inputFile', action='store', help='-i <ALLfilename> : input ALL filename to image. It can also be a wildcard, e.g. *.all')
    parser.add_argument('-o', dest='outputFile', action='store', default='', help='-o <SHPfilename.shp> : output filename to create. e.g. trackplot.shp [Default: track.shp]')
    parser.add_argument('-s', dest='step', action='store', default='30', help='-s <step size in seconds> : decimate the position datagrams to reduce the shapefile size.  Some systems record at 100Hz.  [Default: 30]')
    parser.add_argument('-c', action='store_true', default=False, dest='coverage', help='-c : create coverage polygon shapefile.  [Default: False]')
    parser.add_argument('-t', action='store_true', default=True, dest='track', help='-t : create track polyline shapefile.  [Default: True]')
    parser.add_argument('-r', action='store_true', default=False, dest='recursive', help='-r : search recursively.  [Default: False]')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    # we need to ensure the file is a shp extension
    if len(args.outputFile) == 0:
        trackFileName = os.path.dirname(args.inputFile) + '/track.shp'
    else:
        trackFileName = args.outputFile 
    if not trackFileName.lower().endswith('.shp'):
        trackFileName += '.shp'

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
    if args.track:
        # create the destination shape files 
        trackLine = createSHP(trackFileName, False)

    if args.coverage:
        # we need to make a coverage shp file, so auto create a name based on the trackline 
        # create a name for the coverage shape file
        coverageFileName = trackFileName.lower().replace('.shp', '_coverage.shp')
        coveragePoly = createSHP(coverageFileName, True)

    for filename in matches:
        # print ("processing file: %s" % filename)
        
        reader = pyall.ALLReader(filename)
        start_time = time.time() # time  the process

        if args.track:
            createTrack(reader, trackLine, float(args.step))

        # create the polygon
        if args.coverage:
            createCoverage(reader, coveragePoly, float(args.step))

        update_progress("Processed: %s (%d/%d)" % (filename, fileCounter, len(matches)), (fileCounter/len(matches)))
        # lastTimeStamp = update[0]
        fileCounter +=1
        reader.close()

    update_progress("Process Complete: ", (fileCounter/len(matches)))
    if args.track:
        print ("Saving track line shapefile: %s" % trackFileName)        
        trackLine.save(trackFileName)
    if args.coverage:
        print ("Saving coverage polygon shapefile: %s" % coverageFileName)        
        coveragePoly.save(coverageFileName)

    # now write out a prj file so the data has a spatial Reference
    prj = open(trackFileName.replace('.shp','.prj'), 'w')
    prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
    prj.close() # you can omit in most cases as the destructor will call it

def createTrack(reader, trackLine, step):
    lastTimeStamp = 0
    # trackRecordCount = 0
    line_parts = []
    line = []
    navigation = reader.loadNavigation()

    # create the trackline shape file
    for update in navigation:
        if update[0] - lastTimeStamp >= step:
            line.append([float(update[2]),float(update[1])])
            # trackRecordCount += 1
            lastTimeStamp = update[0]
    # now add the very last update
    line.append([float(navigation[-1][2]),float(navigation[-1][1])])
        
    line_parts.append(line)
    trackLine.line(parts=line_parts)
    # now add to the shape file.
    recDate = from_timestamp(navigation[0][0]).strftime("%Y%m%d")
    trackLine.record(os.path.basename(reader.fileName), int(navigation[0][0]), recDate) 

    
def createCoverage(reader, coveragePoly, step):
    lastTimeStamp = 0
    parts = []
    left = []
    right = []
    selectedPositioningSystem = None
    Latitude = 0;
    Longitude = 0
    leftLatitude = 0;
    leftLongitude = 0
    rightLatitude = 0;
    rightLongitude = 0

    reader.rewind()
    while reader.moreData():
        TypeOfDatagram, datagram = reader.readDatagram()
        if (TypeOfDatagram == 'P'):
            datagram.read()
            # recDate = self.currentRecordDateTime()
            if (selectedPositioningSystem == None):
                selectedPositioningSystem = datagram.Descriptor
            if (selectedPositioningSystem == datagram.Descriptor):
                latitude = datagram.Latitude
                longitude = datagram.Longitude
        if TypeOfDatagram == 'D':
            datagram.read()
            leftLatitude, leftLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, datagram.Heading - 90 , math.fabs(datagram.AcrossTrackDistance[0]))

            rightLatitude, rightLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, datagram.Heading + 90, math.fabs(datagram.AcrossTrackDistance[-1]))

        if TypeOfDatagram == 'X':
            datagram.read()
            if (math.fabs(datagram.AcrossTrackDistance[0]) > 0) and (math.fabs(datagram.AcrossTrackDistance[-1]) > 0) :
                leftLatitude, leftLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, datagram.Heading - 90, math.fabs(datagram.AcrossTrackDistance[0]))

                rightLatitude, rightLongitude, leftAz = geodetic.calculateGeographicalPositionFromRangeBearing(latitude, longitude, datagram.Heading + 90,  math.fabs(datagram.AcrossTrackDistance[-1]))

        # add to the shape file at the user required interval
        if to_timestamp(reader.currentRecordDateTime()) - lastTimeStamp >= step:
            if leftLongitude > 0:
                left.append([leftLongitude,leftLatitude])        
                right.insert(0,[rightLongitude,rightLatitude])        
                lastTimeStamp = to_timestamp(reader.currentRecordDateTime())
    
    poly = []
    for p in left:
        poly.append(p)
    for p in right:
        poly.append(p)
    # left.append(right)
    # parts.append(left)
    # reverse the list on the right side so we get a nice polygon.
    # right.reverse()
    parts.append(poly)
    coveragePoly.poly(parts=parts)
    recDate = from_timestamp(lastTimeStamp).strftime("%Y%m%d")
    coveragePoly.record(os.path.basename(reader.fileName), int(lastTimeStamp), recDate) 
        
    return coveragePoly

def createSHP(fileName, createPolyGon):
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
        if createPolyGon:
            writer = shapefile.Writer(shapefile.POLYGON)
        else:
            writer = shapefile.Writer(shapefile.POLYLINE)
        writer.autoBalance = 1
        writer.field("LineName", "C")
        # w.field("WaterDepth", "N")
        writer.field("UNIXTime", "N")
        writer.field("SurveyDate", "D")
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

if __name__ == "__main__":
    main()

