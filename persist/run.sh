#!/bin/bash

echo "radius based connectivity"

for radius in {1..9}
do
	
	echo "sphere radius: " $radius " mm"
	./persist CC003.ply "cc03_"$radius".ply" $radius

done

echo "done!"
