#!/usr/bin/gawk -f
#BEGIN{
#   if(ARGC!=2){
#     print ARGC
#     print "./getcluster.awk infofile";
#     print "INFO Format:"
#     print "event enentnumber cluster clusternumer eventname"
#     exit;
#}
#} 
{ event[$4,++cluster[$4]]=$5;}
END{
   n=1;
   for(i in cluster){
     # print i,cluster[i];
      for(j=1;j<=cluster[i];j++){
          printf("%s %3d %s\n","cluster",n,event[i,j]);
         }
   n++;
     }   
      
}
