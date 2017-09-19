from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import os
from collections import defaultdict
import geodetic
import math
import numpy as np
import pyall
import sys
import time
import os.path
import warnings

# ignore numpy NaN warnings when applying a mask to the images.
warnings.filterwarnings('ignore')

##############################################################################
def main():
	parser = ArgumentParser(description='Read Kongsberg ALL file and create an ESRI shape file of the trackplot.',
			epilog='Example: \n To process a single file use -i c:/temp/myfile.all \n to mass process every file in a folder use -i c:/temp/*.all\n To convert all .all files recursively in a folder, use -r -i c:/temp \n To convert all .all files recursively from the current folder, use -r -i ./*.all \n', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', dest='inputFile', action='store', help='-i <ALLfilename> : input ALL filename to image. It can also be a wildcard, e.g. *.all')
	parser.add_argument('-o', dest='outputFile', action='store', default='coverageraster', help='-o <rasterfilename> : output raster filename to create. e.g. coverageraster [Default: coverageraster]')
	parser.add_argument('-g', dest='gridsize', action='store', default='1', help='-g <grid cell size> : output raster cell size in metres to create. e.g. [Default: 1]')
	
	args = parser.parse_args()
	start_time = time.time() # time the process
	convert(args.inputFile, args.outputFile)
	# grid = defaultdict(dict)

	# for x in range (1000):
	# 	for y in range (5000, 6000):
	# 		grid[x][y] = x * y


	# # for key, value in d.iteritems(): 2.7
	# for x, row in grid.items():
	# 	for y, value in row.items():
	# 		print (x, y, value)
	# print (grid[0,-45])
	# print (grid[0,-45])

	# print("Duration: %.3f seconds, Count" % (time.time() - start_time)) # print the processing time. It is handy to keep an eye on processing performance.

def convert(fileName, outputFile):	
	recCount = 0

	npx = np.array([])
	npy = np.array([])
	npz = np.array([])

	r = pyall.ALLReader(fileName)
	eprint("loading navigation...")
	navigation = r.loadNavigation()
	# eprint("done.")

	arr = np.array(navigation)
	times = arr[:,0]
	latitudes = arr[:,1]
	longitudes = arr[:,2]

	start_time = time.time() # time the process

	ptsx = []
	ptsy = []
	ptsz = []
	while r.moreData():
		TypeOfDatagram, datagram = r.readDatagram()
		if (TypeOfDatagram == 'X') or (TypeOfDatagram == 'D'):
			datagram.read()
			recDate = r.currentRecordDateTime()

			if datagram.NBeams > 1:
				# interpolate so we know where the ping is located
				lat = np.interp(pyall.to_timestamp(recDate), times, latitudes, left=None, right=None)
				lon = np.interp(pyall.to_timestamp(recDate), times, longitudes, left=None, right=None)
				latRad = math.radians(lat)
				lonRad = math.radians(lon)

				# needed for an optimised algorithm
				localradius = calculateradiusFromLatitude(lat)
				
				# for each beam in the ping, compute the real world position
				for i in range(len(datagram.Depth)):
					#native python version are faster than numpy
					# given the Dx,Dy soundings, compute a range, bearing so we can correccttly map out the soundings
					brg = 90 - ( (180 / math.pi) * math.atan2(datagram.AlongTrackDistance[i], datagram.AcrossTrackDistance[i]) )
					rng = math.sqrt( (datagram.AcrossTrackDistance[i]**2) + (datagram.AlongTrackDistance[i]**2) )

					x,y = destinationPoint(lat, lon, rng, brg + datagram.Heading, localradius)
					z = datagram.Depth[i] + datagram.TransducerDepth
					ptsx.append(x)
					ptsy.append(y)
					ptsz.append(z)
					# print ("%.10f, %.10f, %.3f" % (x, y, z))
			# if len(ptsx) > 100000:
			# 	break

			recCount = recCount + 1
			if recCount % 100 == 0:
				print (recCount)

	npx = np.append(npx, ptsx)
	npy = np.append(npy, ptsy)
	npz = np.append(npz, ptsz)
	print("Duration %.3fs" % (time.time() - start_time )) # time the process

	import matplotlib.pyplot as plt
	import matplotlib.mlab as ml

	binsize_m = 10
	cellsize = binsize_m / (60 * 1852) 
	xi = np.linspace(np.min(npx), np.max(npx, cellsize))
	yi = np.linspace(np.min(npy), np.max(npy, cellsize))
	grid, bins, binloc = griddata2(npx, npy, npz, binsize=cellsize)

	# minimum values for colorbar. filter our nans which are in the grid
	zmin    = grid[np.where(np.isnan(grid) == False)].min()
	zmax    = grid[np.where(np.isnan(grid) == False)].max()

	# colorbar stuff
	palette = plt.matplotlib.colors.LinearSegmentedColormap('jet3',plt.cm.datad['jet'],2048)
	palette.set_under(alpha=0.0)

	# plot the results.  first plot is x, y vs z, where z is a filled level plot.
	extent = (npx.min(), npx.max(), npy.min(), npy.max()) # extent of the plot
	# plt.subplot(1, 2, 1)
	plt.imshow(grid, extent=extent, cmap = plt.get_cmap('rainbow'), origin='lower', vmin=zmin, vmax=zmax, aspect='auto', interpolation='bilinear')
	# plt.xlabel('X values')
	# plt.ylabel('Y values')
	# plt.title('Z = F(X, Y)')
	plt.colorbar()
	plt.show()

	# zi = ml.griddata(npx, npy, npz, xi, yi, interp='linear')
	# xi, yi, zi = np.meshgrid(npx, npy, npz, sparse=True) 
	# plt.contour(xi, yi, zi, 15, linewidths = 0.5, colors = 'k')
	# plt.pcolormesh(xi, yi, zi, cmap = plt.get_cmap('rainbow'))
	# plt.pcolormesh(grid[0], grid[1], bins, cmap = plt.get_cmap('rainbow'))
	# plt.colorbar() 
	# plt.scatter(x, y, marker = 'o', c = 'b', s = 5, zorder = 10)
	# plt.xlim(np.min(npx), np.max(npx))
	# plt.ylim(np.min(npy), np.max(npy))
	# plt.show()

	r.close()
	print("Duration %.3fs" % (time.time() - start_time )) # time the process

	return navigation

def destinationPoint(lat1, lon1, distance, bearing, radius):
	'''
	given a latitude, longitude in degrees, a bearing in metres, a distance in metres
	compute a new geopgraphical position in an efficient manner.
	'''
	radius = 6371000

	# // sinlatitude = sinlat1⋅cosangulardist + coslat1⋅sinangulardist⋅cosbearing
	# // tanangulardistλ = sinbearing⋅sinangulardist⋅coslat1 / cosangulardist−sinlat1⋅sinlatitude
	# // see http://williams.best.vwh.net/avform.htm#LL

	angulardist = distance / radius # angular distance in radians
	bearing = math.radians(bearing)

	lat1 = math.radians(lat1) 
	lon1 = math.radians(lon1) 

	sinlat1 = math.sin(lat1)
	coslat1 = math.cos(lat1)
	sinangulardist = math.sin(angulardist)
	cosangulardist = math.cos(angulardist)
	sinbearing = math.sin(bearing)
	cosbearing = math.cos(bearing)

	sinlatitude = sinlat1*cosangulardist + coslat1*sinangulardist*cosbearing
	latitude = math.asin(sinlatitude)
	y = sinbearing * sinangulardist * coslat1
	x = cosangulardist - sinlat1 * sinlatitude
	longitude = lon1 + math.atan2(y, x)

	return ((math.degrees(longitude)+540) % 360-180, (math.degrees(latitude)+540) % 360-180) # normalise to −180..+180°
	# return new LatLon(latitude.toDegrees(), (longitude.toDegrees()+540)%360-180); //

##################################################################################
def calculateradiusFromLatitude(lat):
	'''
	given a latitude compute a localised earth radius in metres using wgs84 ellipsoid 
	https://rechneronline.de/earth-radius/
	'''
	r = 6378.137 # semi major axis for wgs84
	rp = 6356.752 # semi minor axis for wgs 84
	B = math.radians(lat)
	cosB = math.cos(B)
	sinB = math.sin(B) 

	R = (((r**2) * cosB)**2 + ((rp**2) * sinB)**2) / ((r * cosB)**2 + (rp * sinB)**2)
	R = math.sqrt(R)
	return R * 1000

##################################################################################

def griddata2(x, y, z, binsize=0.01, retbin=True, retloc=True):
    """
    Place unevenly spaced 2D data on a grid by 2D binning (nearest
    neighbor interpolation).
    
    Parameters
    ----------
    x : ndarray (1D)
        The idependent data x-axis of the grid.
    y : ndarray (1D)
        The idependent data y-axis of the grid.
    z : ndarray (1D)
        The dependent data in the form z = f(x,y).
    binsize : scalar, optional
        The full width and height of each bin on the grid.  If each
        bin is a cube, then this is the x and y dimension.  This is
        the step in both directions, x and y. Defaults to 0.01.
    retbin : boolean, optional
        Function returns `bins` variable (see below for description)
        if set to True.  Defaults to True.
    retloc : boolean, optional
        Function returns `wherebins` variable (see below for description)
        if set to True.  Defaults to True.
   
    Returns
    -------
    grid : ndarray (2D)
        The evenly gridded data.  The value of each cell is the median
        value of the contents of the bin.
    bins : ndarray (2D)
        A grid the same shape as `grid`, except the value of each cell
        is the number of points in that bin.  Returns only if
        `retbin` is set to True.
    wherebin : list (2D)
        A 2D list the same shape as `grid` and `bins` where each cell
        contains the indicies of `z` which contain the values stored
        in the particular bin.

    Revisions
    ---------
    2010-07-11  ccampo  Initial version
    """
    # get extrema values.
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()

    # make coordinate arrays.
    xi      = np.arange(xmin, xmax+binsize, binsize)
    yi      = np.arange(ymin, ymax+binsize, binsize)
    xi, yi = np.meshgrid(xi,yi)

    # make the grid.
    grid           = np.zeros(xi.shape, dtype=x.dtype)
    nrow, ncol = grid.shape
    if retbin: bins = np.copy(grid)

    # create list in same shape as grid to store indices
    if retloc:
        wherebin = np.copy(grid)
        wherebin = wherebin.tolist()

    # fill in the grid.
    for row in range(nrow):
        for col in range(ncol):
            xc = xi[row, col]    # x coordinate.
            yc = yi[row, col]    # y coordinate.

            # find the position that xc and yc correspond to.
            posx = np.abs(x - xc)
            posy = np.abs(y - yc)
            ibin = np.logical_and(posx < binsize/2., posy < binsize/2.)
            ind  = np.where(ibin == True)[0]

            # fill the bin.
            bin = z[ibin]
            if retloc: wherebin[row][col] = ind
            if retbin: bins[row, col] = bin.size
            if bin.size != 0:
                binval         = np.median(bin)
                grid[row, col] = binval
            else:
                grid[row, col] = np.nan   # fill empty bins with nans.

    # return the grid
    if retbin:
        if retloc:
            return grid, bins, wherebin
        else:
            return grid, bins
    else:
        if retloc:
            return grid, wherebin
        else:
            return grid
##############################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

##############################################################################
if __name__ == "__main__":
	main()

