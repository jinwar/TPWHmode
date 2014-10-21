#!/bin/csh

foreach event (`ls ../gsdf_measure/CSmeasure/*.mat | cut -c 27-38`)
	echo $event
	cp pa5_01_allevent_eventmat/$event* eventmat/
end
