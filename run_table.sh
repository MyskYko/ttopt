dirname=$1
nin=$2
nout=$3
rarity=$4
ninl=$5
noutl=$6
ninf=$7
if [ -z "$dirname" ] || [ -z "$nin" ] || [ -z "$nout" ] || [ -z "$rarity" ]
then
    echo "specify project-dir num-neuron-inputs num-neuron-outputs rarity"
    exit 1
fi
if [ -z "$ninl" ]
then
    ninl=$nin
fi
if [ -z "$noutl" ]
then
   noutl=$nout
fi
if [ -z "$ninf" ]
then
    ninf=$nin
fi
nfile=`ls ${dirname}/blif | wc -l`

i=0
f="layer${i}.blif"
./build/ttopt_table ${dirname}/blif/${f} ${dirname}/blifopt/${f} ${dirname}/train${i}.sim $nout $rarity
#abc -c "read ${dirname}/blifopt/${f}; strash; &get; &lnetmap -I ${ninf} -O ${nout}; write_blif ${dirname}/blifmap/${f}; write_verilog -fm ${dirname}/veropt/${f}.v"

for i in `seq 1 $(($nfile - 2))`
do
    f="layer${i}.blif"
    ./build/ttopt_table ${dirname}/blif/${f} ${dirname}/blifopt/${f} ${dirname}/train${i}.sim $nout $rarity
#    abc -c "read ${dirname}/blifopt/${f}; strash; &get; &lnetmap -I ${nin} -O ${nout}; write_blif ${dirname}/blifmap/${f}; write_verilog -fm ${dirname}/veropt/${f}.v"
done

i=$(($nfile - 1))
f="layer${i}.blif"
./build/ttopt_table ${dirname}/blif/${f} ${dirname}/blifopt/${f} ${dirname}/train${i}.sim $noutl $rarity
#abc -c "read ${dirname}/blifopt/${f}; strash; &get; &lnetmap -I ${ninl} -O ${noutl}; write_blif ${dirname}/blifmap/${f}; write_verilog -fm ${dirname}/veropt/${f}.v"

#./postsynth.sh $dirname
