import sys
print(sys.argv)

dirname = sys.argv[1]
layerid = sys.argv[2]

layer = open(dirname + f'/layer{layerid}.v')
line = layer.readline()
innum = int(line[line.find('[')+1:line.find(':')]) + 1
line = line[line.find('output'):]
outnum = int(line[line.find('[')+1:line.find(':')]) + 1

blif = open(f'layer{layerid}.blif', 'w')
blif.write(f'.model layer{layerid}\n')
blif.write('.inputs')
for i in range(innum):
    blif.write(f' M0[{i}]')
blif.write('\n')
blif.write('.outputs')
for i in range(outnum):
    blif.write(f' M1[{i}]')
blif.write('\n')

ttid = 0
curoutnum = 0
while layer:
    line = layer.readline()
    words = line.split()
    if words and words[0] == 'endmodule':
        break
    if words and words[0] == 'wire':
        inputs = line[line.find('{')+1:line.find('}')].split()
        inputs = ''.join(inputs).split(',')
        pats = []
        tt = open(dirname + f'/layer{layerid}_N{ttid}.v')
        words = ''
        while not words or words[0] != 'case':
            words = tt.readline().split()
        while True:
            line = tt.readline()
            if not line.split():
                break
            inpat = line[line.find('b')+1:line.find(':')]
            line = line[line.find('='):]
            outpat = line[line.find('b')+1:line.find(';')]
            pats.append((inpat, outpat))
        ttoutnum = len(pats[0][1])
        for j in range(ttoutnum):
            blif.write('.names ')
            blif.write(' '.join(inputs))
            blif.write(f' M1[{curoutnum}]\n')
            for pat in pats:
                if pat[1][ttoutnum-j-1] == '1':
                    blif.write(pat[0])
                    blif.write(' 1\n')
            curoutnum += 1
        ttid += 1
        
blif.write('.end')
