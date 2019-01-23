pyall2shp
=====
* created:          November 2016
* version           2.14
* by:               p.kennedy@fugro.com
* description:      python module to read an Kongsberg .ALL file and create a shape file containing the trackplot in WGS84
* notes:            See main at end of script for example how to use this
* based on ALL file version October 2013
* developed for Python version 3.4+ and 2.7 (for ESRI compatibilty)

2DO
===

Done
====
* moved shared code to ../shared
* added option to run in non-multicore.  this is good for debugging
* update shp file to be SSDM compliant trackplot **note shape file field length is contstrained, so not all fieldsnames comply as they are truncated.
* add  support for user to specify EPSG code
* update to use pyproj so we can create East/North files
* add support for multiprocessing using many cores
* https://sebastianraschka.com/Articles/2014_multiprocessing.html#multi-threading-vs-multi-processing
* update shp file to be SSDM compliant trackplot **note shape file field length is contstrained, so not all fieldsnames comply as they are truncated.
* fixed bug in date field on coverage (thanks Rhiannon)
* fixed bug if we open to append to a shape file instead of creating a new one.
* now adds the hapes as layers to arcmap if the tool is run withing arcmap
* improved coverage polygon if the vessel is running circular lines
* tested against python 2.7
* added wrapper so the tool works as an ESRI Arc Tool in ArcGIS 10.4
* improved auto filenaming convention
* coverage shp file is now smoother
* decode now uses R_Runtime datagrams
* Added WGS84 prj file so ESRI has a spatial reference.  .ALL files are always WGS, so this is an easy one.
* Added UNIXTime and SurveyDate field for Geosciences Australia
* Added backward compatibility for ArcGIS 2.7. It now works on Python 2.7 and 3.4
* Basic implementation.  It uses the pyall built in function to extract all navigation into a list.  this makes it a small script
* added option to export coverage details to CSV so we can analyse performance.
* Uses pyshp to create a shape file
* If the shp file exists, it appends
* If the shp file does not exist, it creates a new one
* Adds the .all filename to an attribute fields
* Now supports wildcard for recursive searches
* Improved user feedbad with progress text

Example
-------

To see the command line help...
```
python pyall2shp.py
```
To create a simple shape file from 1 .all file...
```
python pyall2shp.py -i <filename.all>
```
To create a single shape file from every .all file in the c:\temp folder, including subfolders...
```
python pyall2shp.py -i c:\temp\*.all -r
```
To create a single shape file from every .all file in the current folder, including subfolders...
```
python pyall2shp.py -i .\ -r
```

output
------
```
A nice polyline shape file with a trackplot at first, last and 30 second intervals.  The attributes include the .all filename, UNIXTime and SurveyDate.  More can be added as needed.
```