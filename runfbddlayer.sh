dirname=$1
i=$2
nout=$3

python3 part.py $dirname $i $nout
npart=`ls ${dirname}/blifpart/layer${i}_* | wc -l`
cp ${dirname}/blifpart/layer${i}.blif ${dirname}/blifopt/.
cp ${dirname}/blifpart/layer${i}_*.blif ${dirname}/blifopt/.
for j in `seq 0 $(($npart - 1))`
do
    timeout 1m ./fbdd ${dirname}/blifopt/layer${i}_${j}.blif > ${dirname}/blifopt/layer${i}_${j}.log
    if [ -e ${dirname}/blifopt/layer${i}_${j}_opt.blif ]
    then
        python3 fixblif.py ${dirname}/blifopt/layer${i}_${j}_opt.blif ${dirname}/blifopt/layer${i}_${j}_fixed.blif
        abc -c "read ${dirname}/blifopt/layer${i}_${j}_fixed.blif; strash; write ${dirname}/blifopt/layer${i}_${j}_flat.blif"
        cat ${dirname}/blifopt/layer${i}_${j}_flat.blif >> ${dirname}/blifopt/layer${i}.blif
    else
        cat ${dirname}/blifopt/layer${i}_${j}.blif >> ${dirname}/blifopt/layer${i}.blif
    fi
done
