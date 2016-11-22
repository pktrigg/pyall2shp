pyall
=====
* created:          November 2016
* version           1.30
* by:               p.kennedy@fugro.com
* description:      python module to read an Kongsberg .ALL file and create a shape file containing the trackplot in WGS84
* notes:            See main at end of script for example how to use this
* based on ALL file version October 2013 
* developed for Python version 3.4 

Done
====
# basic implementation.  It uses the pyall built in function to extract all navigation into a list.  this makes it a small script
# uses pyshp to create a shape file
# if the shp file exists, it appends
# if the shp file does not exist, it creates a new one
# adds the .all filename to an attribute fields


2Do
===
#write in SSDM schema format

Example
-------

to see the 
'''
python pyall2shp.py
'''
python pyall2shp.py -i <filename.all>  
python pyall2shp.py -i c:\temp -r  

output
------
'''
<?xml version="1.0"?>
<HIPSVesselConfig Version="2.0">
  <VesselShape>
    <PlanCoordinates>
      <Entry X="-2.000000" Y="-0.000000"/>
      <Entry X="2.000000" Y="-0.000000"/>

```
```