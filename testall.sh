#! /bin/bash

$path=results
if [ -d "$path" ]; then
	datetime=`date '+%Y%m%d-%H%M%S'`
	newpath="$path.$datetime"
	mv $path $newpath
fi

./test2d.sh $path
./test3d.sh $path

