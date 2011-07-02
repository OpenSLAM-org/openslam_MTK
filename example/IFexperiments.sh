#!/bin/bash

# run experiments for Information Fusion paper
# this script needs to be run from main directory of slom, with build-directory build

# data sources:
# ../iSAM/data/sphere400.txt
# ../iSAM/data/sphere2500.txt
# ../g2o-prerelease-udo/datasets/3D/stanford-garage/parking-garage.g2o

# 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0


for DS in "box .graph Q" "parking-garage .g2o Q" "sphere400 .txt E" "sphere2500 .txt E"
do
	set -- $DS
	echo $1
	# angle noise in milli-radian (with leading zeroes)
	for noise in 200 100 050 020 010 005 002 001 000
	do
		for orient in M E S Q
		do
			for lambda in + -
			do
				build/example/relation -v 0 -d 3  -a$3 -s -100 -N 0.$noise --seed 1337 -m logs/$1$2 -l${lambda}1e-3 -T$orient -S$1-$noise-$orient$lambda.dat
			done
		done
	done
done
