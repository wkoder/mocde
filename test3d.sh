#! /bin/bash

path=$1
if [ -z "$path" ]; then
	path="mocde"

	if [ -d "$path" ]; then
		datetime=`date '+%Y%m%d-%H%M%S'`
		newpath="$path.$datetime"
		mv $path $newpath
	fi
fi

./test.sh dtlz1 12 $path
./test.sh dtlz2 12 $path
#./test.sh r_dtlz2 10 $path
./test.sh dtlz3 12 $path
#./test.sh dtlz5im 10 $path
./test.sh dtlz7 12 $path
./test.sh uf8 30 $path
./test.sh uf9 30 $path
./test.sh uf10 30 $path

