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

./test.sh deb2 2 $path
./test.sh deb3 2 $path
./test.sh fonseca2 5 $path
./test.sh kursawe 3 $path
./test.sh zdt1 30 $path
./test.sh zdt2 30 $path
./test.sh zdt3 30 $path
#./test.sh wfg1 $path
#./test.sh wfg2 $path
#./test.sh wfg6 $path
./test.sh uf1 30 $path
./test.sh uf2 30 $path
./test.sh uf3 30 $path
./test.sh uf4 30 $path
./test.sh uf5 30 $path
./test.sh uf6 30 $path
./test.sh uf7 30 $path

