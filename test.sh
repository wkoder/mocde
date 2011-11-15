#! /bin/bash

path=$3
mkdir -p $path
varfile=$path/$1"_var.out"
objfile=$path/$1"_front.out"
cmd="./mocderunner $1 $2"
echo "----------------------------------"
echo "Running '$cmd'"
time $cmd $varfile $objfile --silent < test.in
echo "----------------------------------"

