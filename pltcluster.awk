#!/usr/bin/gawk -f
# Usage: read out put from getcluster.awk whose format is 
#        cluster  clusternm filenm(single component)
#        plot the clusters.
# the events within one cluster will have  the same color,and the 
# color between clusters will change from red blue black in turn.
# just for fun and esay display.
BEGIN{
     #set time parameters
     time="t1";
     btime=-1; #before time,also the min of Rx
     atime=10;  #after time,also the max of Rx
     cut1=btime;  #cut begin before time
     cut2=atime;   #cut end after time
     #set the output file
     #ps="swarm.ps"
     #set the GMT layout and paper size
     print "gmtset PAPER_MEDIA a4";
     ncluster=0;
     yshift=0;
     colorarray[1]="black";
     colorarray[0]="blue";
     colorarray[2]="red";
}
{
   n=NR;    
   f[$2,++cluster[$2]]=$3;
   n++;
 }
END{
    print n
    print "psbasemap -JX5i/10i -R"btime"/"atime"/0/"n" -Ba1/a10WSen -P -K >"ps
#find the optimal size.We will shift Y later,so make this as accurate as possible.
    size=sprintf("%.6f",10/n);
    for(i in cluster){
       ncluster++;
# set color here
    ncolor=ncluster%3; 
    color=colorarray[ncolor];
    printf "pssac2 "
    for(j=1;j<=cluster[i];j++){ #note:n is number_of_file+1
        printf "%s ",f[i,j];
       }
    printf " -J -R -En"time" -M"size"i -C"cut1"/"cut2" -W1p/"color" -Y"yshift"i  -P -O -K >>"ps"\n"
# set -Y here
    yshift=size*cluster[i];
    sumyshift+=yshift;
      }
#end of each cluster
#now shift all back for print the event name number.  
    lastyshift=size*cluster[i];
    sumyshift-=lastyshift;
#now create the file for pstext
    print "cat << EOF |pstext -J -R -D0.2/0 -N -Y-"sumyshift"i -P -O >>"ps
       for(ii in cluster){
          for(jj=1;jj<=cluster[ii];jj++){
          ncount++;
          print atime,ncount,10,0,0,"LM",substr(f[ii,jj],8,12); 
          }
     }
    print "EOF"
#end of the gmt scripts
    print "gs "ps" &"
}
