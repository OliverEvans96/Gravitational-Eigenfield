#!/bin/bash
imgdir=$1
outfile=$2
nloops=$3
rm -rf tmp
mkdir tmp -p
k=0
echo "Copying files"
for i in `seq 1 $nloops`
do
	# Skip first image
	for f in `ls $imgdir | awk '$1 !~ /001/{print}'`
	do
		name=`printf 'img%04d.png' $k`
		k=$((k+1))
		rm -rf tmp/$name
		cp $imgdir/$f tmp/$name
	done
done

echo "Converting to movie"
ffmpeg -i tmp/img%04d.png -c:v libx264 -crf 24 -pix_fmt yuv420p -r 30 -framerate 30 -y $outfile
rm -rf tmp

echo "Done!"
