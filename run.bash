#/bin/bash
#Date:2017-07-07
#HuJing @USTC
#email:jhsa26@mail.ustc.edu.cn

# you just need to modify the getcc1D.pl script to define the windows and shift 

# and the following two parameters ############################
nevent=5000  # only the first nevent events in file.list_old are used to do cluster
ccthreshold=0.85
minev=5      # only cluster with more than 5 repeating earthquakes is reserved           

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


#step3 hierarchical clustering 
#always we want to use the maximum method to make event whithin
./bin/cc2clusterHierarchical  ./output/ccfile $nevent $ccthreshold m 
#step4 output the final figure

./bin/cc2clusterHierarchical  ./output/ccfile $nevent $ccthreshold m |paste - ./output/file.list |\
    ./getcluster.awk  | uniq -w11 -D   > "all_cluster"

rm -rf Cluster/*
awk -v minev="$minev" 'BEGIN{count=0;}
{
old_cluster=$1"_"$2"_cluster"
if(NR==1){
    count=1;
    #print old_cluster,count,NR
    print $0 >>old_cluster; 
temp_cluster=old_cluster}

else{
#print old_cluster,temp_cluster
if(old_cluster==temp_cluster){
count = count +1
print $0 >> old_cluster
#print old_cluster,count,NR
}
else{
#print old_cluster,count,NR
    if(count < minev){print "rm -rf ",temp_cluster }
    print $0 >> old_cluster
    count =1
}
temp_cluster=old_cluster
}

}
END{
if(count<minev) {print "rm -rf ",temp_cluster}
}
' "all_cluster" | sh
mv cluster_*_cluster Cluster


for i in `ls Cluster/cluster_*_cluster`
do
    ./plot_repeating.sh $i
done

exit



rm -rf output/*.ps
rm -rf output/*.jpg
for i in `ls Cluster/cluster_*_cluster`
do
    temp=`echo $i |awk 'BEGIN{FS="/"}{print $2}' `
    echo $temp
    cat $i | ./pltcluster.awk -v ps="./output/${temp}.ps" |sh
done
