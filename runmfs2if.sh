f="dct.blif"
abc -c "read_blif ${f}; mfs2; &get; &if -K 6 -a; &put; write_blif ${f}.abc.blif"

./postsynth.sh $f.abc.blif
