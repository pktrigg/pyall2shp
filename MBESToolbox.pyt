import arcpy
from collections import namedtuple
# from glob import glob
# import pyall
# import shapefile
# import geodetic
import pyall2shp2
import apple

class Toolbox(object):
	def __init__(self):
		"""Define the toolbox (the name of the toolbox is the name of the
		.pyt file)."""
		self.label = "MBES Toolbox"
		self.alias = ""

		# List of tool classes associated with this toolbox
		self.tools = [all2shp]


class all2shp(object):
	def __init__(self):
		"""Define the tool (tool name is the name of the class)."""
		self.label = "all file coverage"
		self.description = "Kongsberg all file coverage and trackplot tool"
		self.canRunInBackground = False

	def getParameterInfo(self):
		#Define parameter definitions

		# Input Features parameter
		param0 = arcpy.Parameter(
			displayName="Input Folder (containing .all files)",
			name="inputfolder",
			datatype="DEFolder",
			parameterType="Required",
			direction="Input")
		
		param1 = arcpy.Parameter(
			displayName="Step size (seconds: default 30)",
			name="step",
			datatype="GPLong",
			parameterType="Optional",
			direction="Input")
		
		param2 = arcpy.Parameter(
			displayName="Create Coverage Polygons",
			name="createcoverage",
			datatype="GPBoolean",
			parameterType="Required",
			direction="Input")

		param3 = arcpy.Parameter(
			displayName="Create Track Lines",
			name="createtracklines",
			datatype="GPBoolean",
			parameterType="Required",
			direction="Input")

		param4 = arcpy.Parameter(
			displayName="Create Track Points",
			name="createtrackpoints",
			datatype="GPBoolean",
			parameterType="Required",
			direction="Input")
			
		param0.value = 'c:\\development\\python'
		param1.value = 10
		param2.value = True
		param3.value = True
		param4.value = False
		
		parameters = [param0, param1, param2, param3, param4]		
		return parameters

	def isLicensed(self):
		"""Set whether tool is licensed to execute."""
		return True

	def updateParameters(self, parameters):
		"""Modify the values and properties of parameters before internal
		validation is performed.  This method is called whenever a parameter
		has been changed."""
		return

	def updateMessages(self, parameters):
		"""Modify the messages created by internal validation for each tool
		parameter.  This method is called after internal validation."""
		return

	def execute(self, parameters, messages):
		"""The source code of the tool."""
		inputFolder = parameters[0].valueAsText + "\\*.all"
		stepsize = int(parameters[1].valueAsText)
		createcoverage = parameters[2].value
		createtrackline = parameters[3].value
		createtrackpoint = parameters[4].value

		# arcpy.AddMessage(inputFolder)
		# arcpy.AddMessage("Coverage:"+ str(createcoverage))
		# arcpy.AddMessage("Lines:"+ str(createtracklines))
		# arcpy.AddMessage("Points:"+ str(createtrackpoints))

		if createcoverage==False and createtracklines==False and createtrackpoints==False:
			arcpy.AddMessage("Nothing to do, quitting.")

		# we have some inputs, so we can call the processing script.
			# we need to ensure the file is a shp extension
		args = {'inputFile': inputFolder, 'step': stepsize, 'coverage': createcoverage, 'recursive': False, 'trackline':createtrackline, 'trackpoint':createtrackpoint, 'outputFile': 'coverage.shp', 'csv':False}
		a = namedtuple('GenericDict', args.keys())(**args)
		
		arcpy.AddMessage("calling3..." + str(a))
		# apple.eat()

		pyall2shp2.process(a)
		arcpy.AddMessage("done...")

		return
