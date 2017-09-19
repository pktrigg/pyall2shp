import os
import time
from collections import defaultdict

##############################################################################
def main():
	start_time = time.time() # time the process
	grid = defaultdict(dict)

	for x in range (1000):
		for y in range (5000, 6000):
			grid[x][y] = x * y

	print("Duration: %.3f seconds, Count" % (time.time() - start_time)) # print the processing time. It is handy to keep an eye on processing performance.

	# for key, value in d.iteritems(): 2.7
	for x, row in grid.items():
		for y, value in row.items():
			print (x, y, value)
	print (grid[0,-45])
	print (grid[0,-45])

	# >>> a[23]
	# {11: 1}
	# >>> a[23][11]
	# 1

##############################################################################
if __name__ == "__main__":
	main()

