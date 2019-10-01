#!/bin/bash

a=$1

for folder in $a/*; do

	
	echo $folder
	ndatasets=`wc -l $folder/*datasets*`
	nmv=`ls -l $folder/snv*/*pvm* | wc -l`; 
	nerror=`tail -n2 $folder/*err | grep "ERROR" | wc -l`; 
	nmemory=`tail -n2 $folder/*err | grep "CANCELLED" | wc -l`; 
	echo $ndatasets $nerror $nmv $nmemory; 
	grep "CANCELLED" $folder/*err
	grep "ERROR" $folder/*out

done > summary26julio_2.txt




for folder in $a/*; do

        echo $folder
        ls -l $folder/snv*/*pvm*
    
done > summary26julio_Analysed_2.txt 

grep mapping summary26julio_2.txt | cut -f2,3 -d ":" > a.txt
grep mapping summary26julio_2.txt | cut -f3 -d "/" | cut -f3 -d"_" | awk -F"." '{print $1}' > b.txt
paste b.txt a.txt > summary26julio_ERROR2.txt
rm b.txt
rm a.txt
		
