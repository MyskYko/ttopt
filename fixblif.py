import sys
import os

filename = sys.argv[1]
gilename = sys.argv[2]
f = open(filename)
g = open(gilename, "w")
line = f.readline()
while not line.startswith('.model'):
    line = f.readline()
g.write(line)
line = f.readline()
while not line.startswith('.model'):
    g.write(line)
    line = f.readline()

g.write('.end\n')

g.write(line)
line = f.readline()
while line:
    g.write(line)
    line = f.readline()
