import os

f = open('dct.old.blif')
g = open('dct.blif', 'w')

line = f.readline()
while line:
    line = f.readline()
    if '-' not in line:
        g.write(line)
        continue
    words = line.split()
    cnt = words[0].count('-')
    print(cnt)
    for i in range(2**cnt):
        for c in words[0]:
            if c != '-':
                g.write(c)
                continue
            g.write(str(i % 2))
            i = i >> 1
        g.write(' ')
        g.write(words[1])
        g.write('\n')

g.flush()
        
os.system('yosys -p "read_blif dct.blif; write_blif dct.blif"')
