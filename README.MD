pyall
=====
* created:          November 2016
* version           1.30
* by:               p.kennedy@fugro.com
* description:      python module to read an Kongsberg .ALL file and create a shape file 
* notes:            See main at end of script for example how to use this
* based on ALL file version October 2013 
* developed for Python version 3.4 

Done
====


2Do
===
import shapefile
# Polygon shapefile we are updating.
# We must include a file extension in
# this case because the file name
# has multiple dots and pyshp would get
# confused otherwise.
file_name = "ep202009.026_5day_pgn.shp"
# Create a shapefile reader
r = shapefile.Reader(file_name)
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
# Add a new polygon
w.poly(parts=[[[-104,24],[-104,25],[-103,25],[-103,24],[-104,24]]])
# Add a new dbf record for our polygon making sure we include
# all of the fields in the original file (r.fields)
w.record("STANLEY","TD","091022/1500","27","21","48","ep")
# Overwrite the old shapefile or change the name and make a copy 
w.save(file_name)

Example Reader
-------
```
    r = ALLReader("C:/development/all/sampledata/EM2040/GeoFocusEM2040400kHzdual-Rx0.5degx1degPitchStabilised.all")
    start_time = time.time() # time the process

    while r.moreData():
        # read a datagram.  If we support it, return the datagram type and aclass for that datagram
        # The user then needs to call the read() method for the class to undertake a fileread and binary decode.  This keeps the read super quick.
        TypeOfDatagram, datagram = r.readDatagram()

        if TypeOfDatagram == 'P':
            datagram.read()
            print ("Lat: %.5f Lon: %.5f" % (datagram.Latitude, datagram.Longitude))

        if TypeOfDatagram == 'X':
            datagram.read()

    print("Read Duration: %.3f seconds" % (time.time() - start_time)) # print the processing time. It is handy to keep an eye on processing performance.

    r.rewind()
    print("Complete reading ALL file :-)")
    r.close()    
    # see the test code in main() at the end of pyall for more details.  Have Fun
```
```
```