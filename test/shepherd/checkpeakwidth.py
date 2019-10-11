import sys

with open(sys.argv[1]) as fin:
    next(fin)
    for line in fin:
        line = line.split('\t')
        print(float(line[1]) - float(line[2]))
