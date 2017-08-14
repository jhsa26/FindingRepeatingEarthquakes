#!/bin/bash
rm -rf corr1d cc2clusterHierarchical

gcc -o corr1d   corr1d.c -lm
gcc -o cc2clusterHierarchical  cc2clusterHierarchical.c

cp corr1d cc2clusterHierarchical ../bin
