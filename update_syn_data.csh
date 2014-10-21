#!/bin/csh

csh cleandata.csh
rm -rf sacdata/*
foreach event (`ls ../gsdf_measure/CSmeasure/*.mat | cut -c 27-38`)
	echo $event
	cp -r ../MINOS/SAC/$event sacdata/
end
#cp -r ../MINOS/SAC/* sacdata/
cd sacdata
ls -d 2* > eventlist
cd ..
