#!/bin/bash


#usage  ./obtain_shift.sh cluster_path  

cluster_name=$1
win=8
shifttime=3
awk -v win=${win} -v shifttime=${shifttime} '{ 


if(NR==1) {reference=$3; 
    print "  ./bin/corr1d  1 " reference , $3 , win , shifttime}
else{
print "  ./bin/corr1d  1 " reference , $3 , win , shifttime

}

}' $cluster_name | sh
