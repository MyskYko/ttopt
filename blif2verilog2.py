f = open('dct_sop_fx.blif')

nnode = 0
line = f.readline()
while line:
    words = line.split()
    for word in words:
        if word.startswith('new_n'):
            word = word[5:-1]
            if nnode < int(word):
                nnode = int(word)
    line = f.readline()
print(nnode)
f = open('dct_sop_fx.blif')
g = open('dct_sop_fx.v', 'w')

g.write('module dct(input [2:0] x, y, u, v, output reg [31:0] coef);\n')
g.write('reg ')
delim = ''
for i in range(nnode + 1):
    g.write(f'{delim}new_n{i}_')
    delim = ', '
g.write(';\n')

g.write('always @(*) begin\n')

out = ''
line = f.readline()
while line:
    line = f.readline().strip()
    while len(line) and line[-1] == '\\':
        line = line[0:-1] + f.readline().strip()
    if line.startswith('.names'):
        if out != '':
            g.write('default : ' + out + f' = 1\'b{1-int(words[1])};\n')
            g.write('endcase\n')
#            g.write('end\n')
        words = line.split()
        out = words[-1]
        nin = len(words) - 2
#        g.write('always @(*) begin\n')
        g.write('casez ({' + ', '.join(words[1:-1]) + '})\n')
    if line.startswith('0') or line.startswith('1') or line.startswith('-'):
        words = line.split()
        words[0] = words[0].replace('-', '?')
        g.write(f'{nin}\'b' + words[0] + ' : ' + out + f' = 1\'b{words[1]};\n')

g.write('default : ' + out + f' = 1\'b{1-int(words[1])};\n')
g.write('endcase\n')
g.write('end\n')
g.write('endmodule\n')
