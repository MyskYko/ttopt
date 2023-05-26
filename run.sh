f="dct.blif"
./build/ttopt ${f} ${f}.ttopt.blif dummy 32 0
abc -c "read ${f}.ttopt.blif; strash; &get; &lnetmap -I 12 -O 32; write_blif ${f}.ttopt.blif"

./postsynth.sh $f.ttopt.blif
