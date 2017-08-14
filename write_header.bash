#!/bin/bash
# delete without Ppick
cd data_sac
saclst a f *.SAC | awk '{if($2==-12345){print "rm -rf ",$1}}' |sh
# mv Ppick header to t1
for i in `ls *.SAC`
do
    t1=`saclst a  f $i |awk '{print $2}'  ` 
    echo "r $i" > "a.m"
    echo "ch t1 $t1" >>"a.m"
    echo "wh" >>"a.m"
    echo "q"  >>"a.m"
sac a.m
done




