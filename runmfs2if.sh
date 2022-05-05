dirname=$1
if [ -z "$dirname" ]
then
    echo "specify project-dir"
    exit 1
fi
i=0
nfile=`ls ${dirname}/blif | wc -l`
for i in `seq 0 $(($nfile - 1))`
do
    f="layer${i}.blif"
    abc -c "read ${dirname}/blif/${f}; mfs2; &get; &if -K 6 -a; &put; write_blif ${dirname}/blifmap/${f}; write_verilog -fm ${dirname}/veropt/${f}.v"
done

./postsynth.sh $dirname
