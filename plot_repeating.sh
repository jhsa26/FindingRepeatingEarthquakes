#!/bin/bash

filename=$1
J="X4i/10i"
xitval=2
yitval=1
btime=-1
atime=10
cut1=-5
cut2=15
outmap=${filename}".ps"
nev=`cat $filename | wc -l | awk '{print $1+1}'`
size=`echo ${nev} |awk '{print 10/$1}'`
size=0.5  #`echo ${nev} |awk '{print 10/$1}'`
Ysize=`echo ${size} ${nev} | awk '{print $1*$2}'`
R=${btime}/${atime}/0/${nev}
J="X4i/${Ysize}i"
psbasemap -J${J} -R${R} -Ba1:"T(s)":/a1wSen -P -K >${outmap}
count=0
awk '{print "saclst kevnm kzdate kztime mag f " $3}' ${filename} | sh  >"evinfo"

./obtain_shift_time.sh ${filename} | paste - ${filename} | paste - evinfo | sort -k 10 -r  |  awk -v cut1=${cut1} -v cut2=${cut2} -v  size=${size} -v outmap=${outmap} -v atime=${atime} 'BEGIN{yshift=0}
{
   dt=0-$4
   print "pssac2  "  $7 " -J  -R -C"cut1"/"cut2  " -Ent1 "    "-S"dt " -M"size"i" " -Y"yshift"i"   " -K  -O -P >>" outmap
   yshift=size
   print "echo ",atime,1," 10 0 0 LM ", $9 "|pstext -J -R -D0.2/0 -N  -P -O -K>>"outmap
   print "echo ",atime+1.2,1," 10 0 0 LM ", $10 "|pstext -J -R -D0.2/0 -N  -P -O -K>>"outmap
   print "echo ",atime+3.5,1," 10 0 0 LM ", $11 "|pstext -J -R -D0.2/0 -N  -P -O -K>>"outmap
   print "echo ",atime+6.8,1," 10 0 0 CM ", $12 "|pstext -J -R -D0.2/0 -N  -P -O -K>>"outmap
}' | sh 

echo "$atime 1.8 11 0 0 LM Event    Date               Time        Magnitude "  | pstext -J -R -D0.2/0 -N  -P -O  >>${outmap}
ps2raster -A -Tj  -E300 ${outmap}
rm -rf ${outmap}

