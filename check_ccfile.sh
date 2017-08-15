#!/bin/bash
# this subroutine is used to check the ccfile
# correct ccfile would be used to cluster
# Date:2017-08-14
# HJ

NR=` more  ./output/ccfile | awk '{if(NR>NF){total=NR} }END{if(total<NR){print total}else{print 0}}' `
cat ./output/index.list | awk -v numth=${NR} '{if(NR==numth){print "line=" ,NR ,"file: ",$2 ," t1-shift <b or t1+shift<e"}}'


