#/bin/bash
#Date:2017-07-07
#HuJing @USTC
#email:jhsa26@mail.ustc.edu.cn

# you just need to modify the getcc1D.pl script to define the windows and shift 

# and the following two parameters ############################
nevent=380
ccthreshold=0.8

########################################
#step1
test -d output || mkdir output
#rm -rf output/*
#1.create the file list and index list
ls -1 data_sac/*.SAC | gawk '{print NR,"./"$0}'  >./output/index.list_old
ls -1 data_sac/*.SAC | gawk '{print "./"$0}'  >./output/file.list_old
head -n${nevent} ./output/file.list_old >./output/file.list
head -n${nevent} ./output/index.list_old >./output/index.list
#step2
######edit the getcorr1d.pl to define the windows and shift
#2.calculate the cc matrix
############################################################
./getcc1D.pl ./output/index.list  >./output/ccfile

#exit

#step3 hierarchical clustering 
#always we want to use the maximum method to make event whithin
./bin/cc2clusterHierarchical  ./output/ccfile $nevent $ccthreshold m 
#step4 output the final figure

./bin/cc2clusterHierarchical  ./output/ccfile $nevent $ccthreshold m |paste - ./output/file.list |\
    ./getcluster.awk | ./pltcluster.awk -v ps="./output/output.ps" |sh
