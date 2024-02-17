import sys
import os

projdirname = sys.argv[1]
layerid = int(sys.argv[2])
outperneuron = int(sys.argv[3])
blifdirname = projdirname + '/blif'
blifpartdirname = projdirname + '/blifpart'

layer = open(blifdirname + f'/layer{layerid}.blif')
line = layer.readline()
while not line.startswith('.inputs'):
    line = layer.readline()
inputs = line.strip().split()[1:]
print(inputs)

while not line.startswith('.outputs'):
    line = layer.readline()
outputs = line.strip().split()[1:]
print(outputs)

topmodel = open(blifpartdirname + f'/layer{layerid}.blif', "w")
topmodel.write(f'.model layer{layerid}\n')
topmodel.write(f".inputs {' '.join(inputs)}\n")
topmodel.write(f".outputs {' '.join(outputs)}\n")

while not line.startswith('.names'):
    line = layer.readline()

for i in range(int(len(outputs) / outperneuron)):
    modelname = f'layer{layerid}_{i}'
    neuron = open(blifpartdirname + f'/{modelname}.blif', "w")
    neuron.write(f'.model {modelname}\n')
    neuron.write(f".inputs {' '.join(inputs)}\n")
    neuronoutputs = outputs[i * outperneuron: (i + 1) * outperneuron]
    neuron.write(f".outputs {' '.join(neuronoutputs)}\n")
    for j in range(outperneuron):
        neuron.write(line)
        line = layer.readline()
        while not line.startswith('.names') and not line.startswith('.end'):
            neuron.write(line)
            line = layer.readline()
    neuron.write('.end\n')

    topmodel.write(f'.subckt {modelname}')
    for inp in inputs:
        topmodel.write(f' {inp}={inp}')
    for oup in neuronoutputs:
        topmodel.write(f' {oup}={oup}')
    topmodel.write('\n')

topmodel.write('.end\n')
