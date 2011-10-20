#! /bin/bash

path=results
mkdir -p $path
varfile=$path/$1"_var.out"
objfile=$path/$1"_front.out"
./mocderunner $1 $2 $varfile $objfile < test.in

