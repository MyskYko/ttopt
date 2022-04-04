dirname=$1
nout=$2
if [ -z "$dirname" ] || [ -z "$nout" ]
then
    echo "specify project-dir num-neuron-outputs"
    exit 1
fi
l=`ls ${dirname}/blifmap/*`
abc -c "putontop ${l}; sw; ps; strash; &get; &lnetsim ${dirname}/train0.sim ${dirname}/train.simo; &lneteval -O ${nout} ${dirname}/train.simo ${dirname}/train_output.txt; &lnetsim ${dirname}/test.sim ${dirname}/test.simo; &lneteval -O ${nout} ${dirname}/test.simo ${dirname}/test_output.txt"
