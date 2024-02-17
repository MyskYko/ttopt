dirname=$1
if [ -z "$dirname" ]
then
    echo "specify project directory"
    exit 1
fi
mkdir "${dirname}/blif"
mkdir "${dirname}/blifopt"
mkdir "${dirname}/blifmap"
mkdir "${dirname}/blifpart"
mkdir "${dirname}/veropt"
python3 conv.py $dirname
abc -c "&lnetread ${dirname}/train_input.txt ${dirname}/train0.sim; &lnetread ${dirname}/test_input.txt  ${dirname}/test.sim"
i=0
for f in `ls ${dirname}/blif`
do
    abc -c "read ${dirname}/blif/${f}; strash; &get; &lnetsim ${dirname}/train${i}.sim  ${dirname}/train$(($i + 1)).sim"
    i=$(($i + 1))
done
