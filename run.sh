dirname=$1
nin=$2
nout=$3
rarity=$4
if [ -z "$dirname" ] || [ -z "$nin" ] || [ -z "$nout" ] || [ -z "$rarity" ]
then
    echo "specify project-dir num-neuron-inputs num-neuron-outputs rarity"
    exit 1
fi
i=0
for f in `ls ${dirname}/blif`
do
    ./build/ttopt ${dirname}/blif/${f} ${dirname}/blifopt/${f} ${dirname}/train${i}.sim $nout $rarity
    abc -c "read ${dirname}/blifopt/${f}; strash; &get; &lnetmap -I ${nin} -O ${nout}; write_blif ${dirname}/blifmap/${f}; write_verilog -fm ${dirname}/veropt/${f}.v"
    i=$(($i + 1))
done
