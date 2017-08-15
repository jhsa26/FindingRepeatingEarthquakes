/* Programe to calculate the nomalized cross correlation coefficent of two sac file.
   (Part of them ).
   First cut the syn array from file1.windows is from hdr.tn to hdr.tn+length.
   With the same winlength,calculate the coefficent with the waveform in file2.
   Cut wave of file2.time from hdr2.tn-dt to hdr2.tn+length+dt.shift from -dt to dt
   to calculate the coefficent,then output the max value and maxindex,then convert
   to the shift time;
   Output the max coefficent and the timeshift.
   NOTICE: 
   Usage:returnCC timeflag file1 file2 winlenth dt
          Initial:mark t1 as the first arrival time in each record.
          Input: file1---event1
                 file2---event2 for compare
                 length---windows length in second
                 dt---   time shift in file2 around file2.tn
          Output:1.winlength---simple print out
                 2.corr-coefficent-----the max value within the shift
                 3.the time shift from the file.o.

                            wangwt@2008
                            revised 2008.10.07*/
/*############################################################################*/
#include<stdio.h>
#include<stdlib.h>
#include"sac.h"
#include <math.h>
float corr1d(float *w1,int idx1,float *w2,int idx2,int len,int norm);
int main(int argc,char*argv[])
{  FILE *file1,*file2; 
   SACHEAD hd1,hd2;
   int flag;
   int npts1,npts2;
   float flagtime1,flagtime2;
   int i,maxindex,minindex;
   int ibegin1,ibegin2,ndt,ncc,nlength;
   float length,dt,delta;
   float maxval,minval,deltatime;
   float *array1,*array2,*cc;
  //usage
   if(argc!=6){
        printf("Usage:corr1d flag file1 file2 windowlengh shiftlengh\n\n");
        printf("Details:\n");
	printf("    flag is the timeID for the beginning of the window\n");
	printf("    -5->header.b -4->header.e -3->header.o -2->header.a\n");
	printf("    n->header.tn -1->I dont know\n");
	printf("    windowlength is the length for crosscorrelation which\n");
	printf("    will shift from -dt to dt around time defined by flag\n");
        printf("Meaning of shift time:\n");
        printf("    1.it is nothing to file1.tn,file1 is ONLY a pattern.\n");
        printf("    2.it means if using the pattern in file1,the real f2.tn should be\n");
        printf("      F2.tn_accurate=F2.tn_old+shifttime\n");
        printf("    3.if we take the time marker in file1 is right,then the S-P time in file1\n");
        printf("      TimeDist1=F1.t_s-F1.t_p\n");
        printf("      while using file1's P and S as pattern,the time distance in file2 should be\n");
        printf("      TimeDist2=(F2.ts+s_shift)-(F2.tp+p_shift)\n");
        printf("               =F2.ts-F2.tp+s_shift-p_shift\n");
        exit(-1);}
    /*number to stand for the time marker in file:
      -5 stands for h.b,
      -4 stands for h.e
      -3 stands for h.o
      -2 stands for h.a
      -1 stands for h. ??     a,o,b,e,internal1
       n stands for h.tn
       You can transform the hdr into int and the value is 
       nhdr+10+number.  Just see below.*/
//now initial the argv
   flag=atoi(argv[1]);
   file1=fopen(argv[2],"rb");
   if(file1==NULL){
       printf("open file %s error!\n",argv[2]);
       exit(-1);
	 }
   file2=fopen(argv[3],"rb");
   if(file2==NULL){
       printf("open file %s error!\n",argv[3]);
       exit(-1);
	 }
   length=atof(argv[4]);
   dt=atof(argv[5]);
//test
//printf("file1:%s\tfile:%s\tlength:%f\tdt:%f\n\n",argv[1],argv[2],length,dt);
//now read in the input file
   if(fread(&hd1,sizeof(SACHEAD),1,file1)!=1){
	printf("read header of file:%s error!\n",argv[2]);
	exit(-1);
      }
   if(fread(&hd2,sizeof(SACHEAD),1,file2)!=1){
	printf("read header of file:%s error!\n",argv[3]);
	exit(-1);
      }
   if(hd1.delta!=hd2.delta){
	printf("The two file should has same delta\n");
	exit(-1);
      }
//now check the time id
   flagtime1=*((float*)(&hd1)+10+flag);
   flagtime2=*((float*)(&hd2)+10+flag);
//test
// printf("timeflag in file1 is %.6f\ttimeflag in file2 is%.6f\n",flagtime1,flagtime2);
   if(fabs(flagtime2-hd2.b)<dt||fabs(hd2.e-length-flagtime2)<dt){
      //printf("the dt value seems too large for file:%s %10.5f\n",argv[3],flagtime1);
      printf("coeff1d\t%f\tshifttime\t%f\n",0.0,-999.0);
      exit(-1);
    }
    npts1=hd1.npts;
    npts2=hd2.npts;
    delta=hd1.delta;
//    printf("npts1=%d\tnpts2=%d\tdelta=%f\n",npts1,npts2,delta);
//now read in the data
    array1=(float*)malloc(sizeof(float)*npts1);
    array2=(float*)malloc(sizeof(float)*npts2);
    if(fread(array1,sizeof(float)*npts1,1,file1)!=1){
	printf("read in data in %s error!\n",argv[2]);
	exit(-1);
     }
    if(fread(array2,sizeof(float)*npts2,1,file2)!=1){
	printf("read in data in %s error!\n",argv[3]);
	exit(-1);
     }
    ibegin1=(flagtime1-hd1.b)/delta+0.5;
    ibegin2=(flagtime2-dt-hd2.b)/delta+0.5;
    nlength=length/delta+1+0.5;
    ndt=dt/delta;
    ncc=ndt*2+1;
//printf("ibegin1=%d\tibegin2=%d\tnlength=%d\tndt=%d\n",ibegin1,ibegin2,nlength,ndt);
    cc=(float*)malloc(sizeof(float)*ncc);
    maxval=-1e20;
    minval=1e20;
    //now calculate the coef array
    for(i=0;i<ncc;i++){
      cc[i]=corr1d(array1,ibegin1,array2,ibegin2+i,nlength,1);
  //    printf("cc[%d]==%f\n",i,cc[i]);
      if(cc[i]>maxval){maxval=cc[i];maxindex=i-ndt;}
      if(cc[i]<minval){minval=cc[i];minindex=i-ndt;}
    }
    if(fabs(maxval)>fabs(minval)){
      printf("coeff1d\t%f\tshifttime\t%f\n",maxval,maxindex*delta);
      }
    else {
      printf("coeff1d\t%f\tshifttime\t%f\n",minval,minindex*delta);
    }
    fclose(file1);
    fclose(file2);
    free(array1);
    free(array2);
    free(cc);
}
//the sub func of corr
float corr1d(float *w1,int idx1,float *w2,int idx2,int len,int norm)
{
    int i;
    double aw1=0,aw2=0,sw1=0,sw2=0,sw12=0;
    for(i=0;i<len;i++){
        aw1+=w1[idx1+i];
        aw2+=w2[idx2+i];
        }
        aw1/=len;
        aw2/=len;
    for(i=0;i<len;i++){
        sw1+=(w1[idx1+i]-aw1)*(w1[idx1+i]-aw1);
        sw2+=(w2[idx2+i]-aw2)*(w2[idx2+i]-aw2);
        sw12+=(w1[idx1+i]-aw1)*(w2[idx2+i]-aw2);
      }
    //printf("s1=%f\ts2=%f\ts12=%f\n",sw1,sw2,sw12);
    if(norm) sw12/=sqrt(sw1*sw2+1e-20);
    return (float)(sw12);
}

