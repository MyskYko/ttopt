dirname=$1
if [ -z "$dirname" ]
then
    echo "specify project-dir"
    exit 1
fi
l=`ls ${dirname}/blifmap/*`
f=${dirname}/all.blif
fnew=${dirname}/all_new.blif
best=`abc -c "putontop ${l}; sw; write_blif $f; ps" | grep -Po "nd =[ ]*\K[0-9]*"`
echo $best

new=`abc -c "read $f; mfs2; write_blif $fnew; ps" | grep -Po "nd =[ ]*\K[0-9]*"`
echo $new

while (( $new < $best ))
do
    best=$new
    mv $fnew $f
    new=`abc -c "read $f; if -K 6 -a; mfs2; write_blif $fnew; ps" | grep -Po "nd =[ ]*\K[0-9]*"`
    echo $new
done

rm $fnew


