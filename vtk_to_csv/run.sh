#!/bin/bash

echo "calculating curvetures"

for file in *; do
	
	echo "file " $file
	./curvatures $file curv_$file

done

echo "done!"
