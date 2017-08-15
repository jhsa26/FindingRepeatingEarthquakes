#!/usr/bin/perl 
#use this script to get the left-triangle matrix of coeff of the files in list;
#first use ls xxx.sac|gawk '{print NR,$1}' >list to get the list.
use strict;
#check the arguments
if(@ARGV!=1)
{ 
  print ("Usage:getcc.pl list \n");
  exit;}

open(LIST,"$ARGV[0]");
my @content;
@content=<LIST>;
chomp @content;
close(LIST);

my @fileid;
my @filenm;
my $i;
for($i=0;$i<@content;$i++){
   ($fileid[$i],$filenm[$i])=split(/\ /,$content[$i]);
#  print "$fileid[$i]\t$filenm[$i]\n";
}
my $j;
for($i=0;$i<@filenm;$i++){
#print "$filenm[$i]";
    for($j=0;$j<=$i;$j++){ 
       my $out;
#print "$filenm[$i]\t$filenm[$j]\n";
       #my $out=`./bin/corr1d 1 $filenm[$i] $filenm[$j] 5 2`;
       # for gofar data                 window length(s) shifttime (s)
       my $out=`./bin/corr1d 1 $filenm[$i] $filenm[$j] 8 3`;
       my @outcc=split(/\t/,$out,4);
       print "$outcc[1]\t";
 }
print "\n";
}
