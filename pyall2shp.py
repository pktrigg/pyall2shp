import sys
sys.path.append("C:/development/Python/pyall")

import argparse
from datetime import datetime
import geodetic
from glob import glob, iglob
import pyall
import time
import os
import os.path
import warnings
from xml.etree.ElementTree import Element, SubElement, Comment, ElementTree
from xml.etree import ElementTree
from xml.dom import minidom
import shapefile
import fnmatch

def main():
    parser = argparse.ArgumentParser(description='Read Kongsberg ALL file and create an ESRI shape file of the trackplot.')
    parser.add_argument('-i', dest='inputFile', action='store', help='-i <ALLfilename> : input ALL filename to image. It can also be a wildcard, e.g. *.all')
    parser.add_argument('-r', action='store_true', default=False, dest='recursive', help='-r : search recursively.  [Default: False]')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    # we need to remember the previous record so we only create uniq values, not duplicates
    fileOut = "track.shp"
    fileCounter=0
    matches = []
    if os.path.isfile(fileOut):
        # Create a shapefile reader
        r = shapefile.Reader(fileOut)
        # Create a shapefile writer
        # using the same shape type
        # as our reader
        w = shapefile.Writer(r.shapeType)
        # Copy over the existing dbf fields
        w.fields = list(r.fields)
        # Copy over the existing dbf records
        w.records.extend(r.records())
        # Copy over the existing polygons
        w._shapes.extend(r.shapes())
    else:
        # w = shapefile.Writer(shapefile.POINTZ)
        w = shapefile.Writer(shapefile.POLYLINE)
        w.autoBalance = 1
        w.field("LineName", "C")
        
    if args.recursive:
        for root, dirnames, filenames in os.walk(args.inputFile):
            for f in fnmatch.filter(filenames, '*.all'):
                matches.append(os.path.join(root, f))
                print (matches[-1])
    else:
        for filename in glob(args.inputFile):
            matches.append(filename)
        print (matches)
    for filename in matches:
        print ("processing file: %s" % filename)
        lastTimeStamp = 0
        trackRecordCount = 0
        line_parts = []
        line = []
        
        r = pyall.ALLReader(filename)
        start_time = time.time() # time  the process
        navigation = r.loadNavigation()
        for update in navigation:
            if update[0] - lastTimeStamp > 30:
                line.append([float(update[2]),float(update[1])])
                # trackRecordCount += 1
                lastTimeStamp = update[0]
        # now add the very last update
        line.append([float(navigation[-1][2]),float(navigation[-1][1])])
            
        line_parts.append(line)
        w.line(parts=line_parts)
        # now add to the shape file.
        w.record(os.path.basename(filename))

        # update_progress("Processed file: %s TrackRecords: %d" % (filename, trackRecordCount), (fileCounter/len(matches)))
        lastTimeStamp = update[0]
        fileCounter +=1
        r.close()

    print ("Saving shapefile: %s" % fileOut)        
    w.save(fileOut)

def update_progress(job_title, progress):
    length = 20 # modify this to change the length
    block = int(round(length*progress))
    msg = "\r{0}: [{1}] {2}%".format(job_title, "#"*block + "-"*(length-block), round(progress*100, 2))
    if progress >= 1: msg += " DONE\r\n"
    sys.stdout.write(msg)
    sys.stdout.flush()

if __name__ == "__main__":
    main()

