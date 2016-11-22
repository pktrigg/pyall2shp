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
* basic implementation.  It uses the pyall built in function to extract all navigation into a list.  this makes it a small script
* uses pyshp to create a shape file
* if the shp file exists, it appends
* if the shp file does not exist, it creates a new one
* adds the .all filename to an attribute fields


2Do
===
* write in SSDM schema format

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
python pyall2shp.py -i c:\temp -r  
```

output
------
```

```