pyall
=====
* created:          November 2016
* version           1.1
* by:               p.kennedy@fugro.com
* description:      python module to read an Kongsberg .ALL file and create a shape file containing the trackplot in WGS84
* notes:            See main at end of script for example how to use this
* based on ALL file version October 2013 
* developed for Python version 3.4 

Done
====
* Added UNIXTime and SurveyDate field for Geosciences Australia
* Added backward compatibility for ArcGIS 2.7. It now works on Python 2.7 and 3.4
* Basic implementation.  It uses the pyall built in function to extract all navigation into a list.  this makes it a small script
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