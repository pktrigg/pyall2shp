# import arcpy
from collections import namedtuple

# arcpy.AddMessage("eat me...")
args = {'inputFile': "inputFolder", 'step': "stepsize", 'coverage': "createcoverage", 'recursive': "False", 'trackline':"createtrackline", 'trackpoint':"createtrackpoint", 'outputFile': 'coverage.shp', 'csv':"False"}
a = namedtuple('GenericDict', args.keys())(**args)
print (a.recursive)		

