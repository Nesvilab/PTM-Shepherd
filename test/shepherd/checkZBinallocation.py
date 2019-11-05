import sys.argv

peps = []

with open(sys.agv[1]) as fin:
	next(fin)
	for line in fin:
		line = line.split('\t')
		if peps
