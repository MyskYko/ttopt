f=$1
if [ -z "$f" ]
then
    echo "specify blif name"
    exit 1
fi
fnew=${f}.new.blif
best=`abc -c "read_blif $f; ps" | grep -Po "nd =[ ]*\K[0-9]*"`
#echo $best

new=`abc -c "read $f; mfs2; write_blif $fnew; ps" | grep -Po "nd =[ ]*\K[0-9]*"`
#echo $new

while (( $new < $best ))
do
    best=$new
    mv $fnew $f
    new=`abc -c "read $f; if -K 6 -a; mfs2; write_blif $fnew; ps" | grep -Po "nd =[ ]*\K[0-9]*"`
#    echo $new
done

rm $fnew

echo $best
