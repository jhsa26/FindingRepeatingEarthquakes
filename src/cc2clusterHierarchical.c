/****************************************************************************/
/*Hierarchical clustering using the normalized cross correlation value matrix.
  Subrountines revised from Cluster3.0     
                                            wangwt@2009@yunnanEarthquakeAdmn */

/*Below is the original license in Cluser Routine*/
/* The C clustering library.
 * Copyright (C) 2002 Michiel Jan Laurens de Hoon.
 *
 * This library was written at the Laboratory of DNA Information Analysis,
 * Human Genome Center, Institute of Medical Science, University of Tokyo,
 * 4-6-1 Shirokanedai, Minato-ku, Tokyo 108-8639, Japan.
 * Contact: mdehoon 'AT' gsc.riken.jp
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */
/****************************************************************************/
#include<stdio.h>
#include<stdlib.h>
/*****************PreDefination**********************************************/
#ifndef min
#define min(x, y)   ((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define max(x, y)   ((x) > (y) ? (x) : (y))
#endif
typedef struct {int left; int right; double distance;} Node;
/****************Subroutines for Hierarchical Clustering *********************/
/****************subroutines prototype****************************************/
double find_closest_pair(int n, double** distmatrix, int* ip, int* jp);
Node* pmlcluster (int nelements, double** distmatrix);
Node* palcluster (int nelements, double** distmatrix);
Node* pslcluster (int nelements, double** distmatrix);
void cuttree (int nelements, Node* tree, int nclusters, int clusterid[]);

/***************Find closest Pair ********************************************/
/******double find_closest_pair(int n, double** distmatrix, int* ip, int* jp)*/
/*
This function searches the distance matrix to find the pair with the shortest
distance between them. The indices of the pair are returned in ip and jp; the
distance itself is returned by the function.
INPUT:
             n:int
                The number of elements in the distance matrix.
    distmatrix:double**
                A ragged array containing the distance matrix. The number of columns in each
                row is one less than the row index.
OUTPUT:
            ip: int* 
                A pointer to the integer that is to receive the first index of the pair with
                the shortest distance.

            jp: int*
                A pointer to the integer that is to receive the second index of the pair with
                the shortest distance.
*/
double find_closest_pair(int n, double** distmatrix, int* ip, int* jp)
    {   int i, j;
        double temp;
        double distance = distmatrix[1][0];
        *ip = 1;
        *jp = 0;
        for (i = 1; i < n; i++)
            { for (j = 0; j < i; j++) 
                { temp = distmatrix[i][j];
                if (temp<distance)
                    { distance = temp;
                      *ip = i;
                      *jp = j;
                    }
                }
            }
    return distance;
}
/***************Maximum Linkage Pairwize Method ******************************/
/*We need use the maximum linkage method if we want either pair of elements or
  cells are within some threshold distance   wangwt@2009@yunnanAirPort */
/************** Node* pmlcluster (int nelements, double** distmatrix)********/
/* 
Purpose
=======
The pmlcluster routine performs clustering using pairwise maximum- (complete-)
linking on the given distance matrix. 
Arguments
=========
INPUT:
      nelements:  int
                The number of elements to be clustered.
     distmatrix: double**
                The distance matrix, with nelements rows, each row being filled up to the
                diagonal. The elements on the diagonal are not used, as they are assumed to be
                zero. The distance matrix will be modified by this routine.
Return value
============
        A pointer to a newly allocated array of Node structs, describing the
        hierarchical clustering solution consisting of nelements-1 nodes. Depending on
        whether genes (rows) or microarrays (columns) were clustered, nelements is
        equal to nrows or ncolumns. See src/cluster.h for a description of the Node
        structure.
        If a memory error occurs, pmlcluster returns NULL.
========================================================================
*/
Node* pmlcluster (int nelements, double** distmatrix)
{ 
    int j;
    int n;
    int* clusterid;
    Node* result;

    clusterid=malloc(nelements*sizeof(int));
    if(!clusterid) return NULL;
    result = malloc((nelements-1)*sizeof(Node));
    if (!result)
       { free(clusterid);
         return NULL;
        }
    /* Setup a list specifying to which cluster an elements belongs */
    for (j = 0; j < nelements; j++) clusterid[j] = j;

    for (n = nelements; n > 1; n--)
        { int is = 1;
          int js = 0;
          result[nelements-n].distance = find_closest_pair(n, distmatrix, &is, &js);
 
     /* Fix the distances */
    for (j = 0; j < js; j++)
        distmatrix[js][j] = max(distmatrix[is][j],distmatrix[js][j]);
    for (j = js+1; j < is; j++)
        distmatrix[j][js] = max(distmatrix[is][j],distmatrix[j][js]);
    for (j = is+1; j < n; j++)
        distmatrix[j][js] = max(distmatrix[j][is],distmatrix[j][js]);
    for (j = 0; j < is; j++) distmatrix[is][j] = distmatrix[n-1][j];
    for (j = is+1; j < n-1; j++) distmatrix[j][is] = distmatrix[n-1][j];
 
     /* Update clusterids */
    result[nelements-n].left = clusterid[is];
    result[nelements-n].right = clusterid[js];
    clusterid[js] = n-nelements-1;
    clusterid[is] = clusterid[n-1];
         }
    free(clusterid);
    return result;
}
/***************Average Linkage Pairwize Method ******************************/
/***************Node* palcluster (int nelements, double** distmatrix)*******/
/*
Purpose
=======
The palcluster routine performs clustering using pairwise average
linking on the given distance matrix.
Arguments
=========
INPUT:
        nelements:  int
            The number of elements to be clustered.

        distmatrix: double**
            The distance matrix, with nelements rows, each row being filled up to the
            diagonal. The elements on the diagonal are not used, as they are assumed to be
            zero. The distance matrix will be modified by this routine.
Return value
============
            A pointer to a newly allocated array of Node structs, describing the
            hierarchical clustering solution consisting of nelements-1 nodes. Depending on
            whether genes (rows) or microarrays (columns) were clustered, nelements is
            equal to nrows or ncolumns. See src/cluster.h for a description of the Node
            structure.
If a memory error occurs, palcluster returns NULL.
========================================================================
*/
Node* palcluster (int nelements, double** distmatrix)
    {
        int j;
        int n;
        int* clusterid;
        int* number;
        Node* result;

        clusterid = malloc(nelements*sizeof(int));
        if(!clusterid) return NULL;
        number = malloc(nelements*sizeof(int));
        if(!number)
           { 
            free(clusterid);
            return NULL;
            }
        result = malloc((nelements-1)*sizeof(Node));
        if (!result)
           {
            free(clusterid);
            free(number);
             return NULL;
            }

  /* Setup a list specifying to which cluster a gene belongs, and keep track
   * of the number of elements in each cluster (needed to calculate the
   * average). */
        for (j = 0; j < nelements; j++)
            { number[j] = 1;
              clusterid[j] = j;
                }

        for (n = nelements; n > 1; n--)
            { int sum;
              int is = 1;
              int js = 0;
              result[nelements-n].distance = find_closest_pair(n, distmatrix, &is, &js);

    /* Save result */
              result[nelements-n].left = clusterid[is];
              result[nelements-n].right = clusterid[js];

    /* Fix the distances */
              sum = number[is] + number[js];
        for (j = 0; j < js; j++)
            { distmatrix[js][j] = distmatrix[is][j]*number[is]
                        + distmatrix[js][j]*number[js];
              distmatrix[js][j] /= sum;
                }
        for (j = js+1; j < is; j++)
            { distmatrix[j][js] = distmatrix[is][j]*number[is]
                        + distmatrix[j][js]*number[js];
              distmatrix[j][js] /= sum;
                }
        for (j = is+1; j < n; j++)
            { distmatrix[j][js] = distmatrix[j][is]*number[is]
                        + distmatrix[j][js]*number[js];
              distmatrix[j][js] /= sum;
                }

        for (j = 0; j < is; j++) distmatrix[is][j] = distmatrix[n-1][j];
        for (j = is+1; j < n-1; j++) distmatrix[j][is] = distmatrix[n-1][j];

    /* Update number of elements in the clusters */
        number[js] = sum;
        number[is] = number[n-1];

    /* Update clusterids */
        clusterid[js] = n-nelements-1;
        clusterid[is] = clusterid[n-1];
        }
        free(clusterid);
        free(number);

        return result;
}
/* ****************Node Compare Sub required by pslcluster***************/
/* ******************************************************************** */
int nodecompare(const void* a, const void* b)
/* Helper function for qsort. */
{
    const Node* node1 = (const Node*)a;
    const Node* node2 = (const Node*)b;
    const double term1 = node1->distance;
    const double term2 = node2->distance;
      if (term1 < term2) return -1;
      if (term1 > term2) return +1;
    return 0;
}
/* ---------------------------------------------------------------------- */
/***************Single Linkage Pairwize Method ******************************/
//Node* pslcluster (int nrows, int ncolumns, double** data, int** mask,
/*In the original cluster library,the the single linkage method using the ori
  data to caltulate the distance in some case.But we have got the cross 
  correlation value matrix when we do the clustering,so I change the prototype
  of that into
    Node* pslcluster (int nelements,double** distmatrix)
                                   wangwt@beijing@20090107 */
/*
Purpose
=======
The pslcluster routine performs single-linkage hierarchical clustering, using
either the distance matrix directly, if available, or by calculating the
distances from the data array.(HERE WE USE THE PRECALCULATED DISTANCEMATRIX)
This implementation is based on the SLINK algorithm, described in:
Sibson, R. (1973). SLINK: An optimally efficient algorithm for the single-link
cluster method. The Computer Journal, 16(1): 30-34.
The output of this algorithm is identical to conventional single-linkage
hierarchical clustering, but is much more memory-efficient and faster. Hence,
it can be applied to large data sets, for which the conventional single-
linkage algorithm fails due to lack of memory.
Arguments
=========
INPUT:
        nelements: int
                  the number of events  we want to cluster
        distance : double **
                  the nelements by nelements distance matrix

Return value
============

        A pointer to a newly allocated array of Node structs, describing the
        hierarchical clustering solution consisting of nelements-1 nodes. Depending on
        whether genes (rows) or microarrays (columns) were clustered, nelements is
        equal to nrows or ncolumns. See src/cluster.h for a description of the Node
        structure.If a memory error occurs, pslcluster returns NULL.
========================================================================
*/
Node* pslcluster (int nelements,double** distmatrix)
{
     double DISTMAX=5.00;
     int i, j, k;
     int nnodes = nelements - 1;
     int* vector;
     double* temp;
     int* index;
     Node* result;
     temp = malloc(nnodes*sizeof(double));
     if(!temp) return NULL;
     index = malloc(nelements*sizeof(int));
     if(!index)
        { free(temp);
          return NULL;
        }
     vector = malloc(nnodes*sizeof(int));
     if(!vector)
        { free(index);
          free(temp);
          return NULL;
        }
     result = malloc(nelements*sizeof(Node));
     if(!result)
        { free(vector);
          free(index);
          free(temp);
          return NULL;
        }

     for (i = 0; i < nnodes; i++) vector[i] = i;

     //for (i = 0; i < nrows; i++)
     for (i = 0; i < nelements; i++)
        { result[i].distance =DISTMAX;
            for (j = 0; j < i; j++) temp[j] = distmatrix[i][j];
            for (j = 0; j < i; j++)
                { k = vector[j];
                  if (result[j].distance >= temp[j])
                    { if (result[j].distance < temp[k]) temp[k] = result[j].distance;
                          result[j].distance = temp[j];
                          vector[j] = i;
                        }
                  else if (temp[j] < temp[k]) temp[k] = temp[j];
                    }
     for (j = 0; j < i; j++)
          {
            if (result[j].distance >= result[vector[j]].distance) vector[j] = i;
            }
        }
     free(temp);

     for (i = 0; i < nnodes; i++) result[i].left = i;
     qsort(result, nnodes, sizeof(Node), nodecompare);

     for (i = 0; i < nelements; i++) index[i] = i;
     for (i = 0; i < nnodes; i++)
        { j = result[i].left;
          k = vector[j];
          result[i].left = index[j];
          result[i].right = index[k];
          index[k] = -i-1;
        }
     free(vector);
     free(index);

     result = realloc(result, nnodes*sizeof(Node));

     return result;
}

/*****************Subroutine to cut tree************************************* */
/***void cuttree (int nelements, Node* tree, int nclusters, int clusterid[])***/
/*
Purpose
=======
The cuttree routine takes the output of a hierarchical clustering routine, and
divides the elements in the tree structure into clusters based on the
hierarchical clustering result. The number of clusters is specified by the user.
Arguments
=========
INPUT:
         nelements: int
                The number of elements that were clustered.
              tree: Node[nelements-1]
                The clustering solution. Each node in the array describes one linking event,
                with tree[i].left and tree[i].right representig the elements that were joined.
                The original elements are numbered 0..nelements-1, nodes are numbered
                 -1..-(nelements-1).
        nclusters: int
                The number of clusters to be formed.
OUTPUT:
         clusterid : int[nelements]
                The number of the cluster to which each element was assigned. Space for this
                array should be allocated before calling the cuttree routine. If a memory
                error occured, all elements in clusterid are set to -1.
========================================================================
*/
void cuttree (int nelements, Node* tree, int nclusters, int clusterid[])
{
     int i, j, k;
     int icluster = 0;
     const int n = nelements-nclusters; /* number of nodes to join */
     int* nodeid;
     for (i = nelements-2; i >= n; i--)
        { k = tree[i].left;
            if (k>=0)
                { clusterid[k] = icluster;
                  icluster++;
                }
            k = tree[i].right;
            if (k>=0)
                { clusterid[k] = icluster;
                  icluster++;
                }
        }
     nodeid = malloc(n*sizeof(int));
     if(!nodeid)
        { for (i = 0; i < nelements; i++) clusterid[i] = -1;
          return;
        }
     for (i = 0; i < n; i++) nodeid[i] = -1;
     for (i = n-1; i >= 0; i--)
        { if(nodeid[i]<0)
            { j = icluster;
              nodeid[i] = j;
              icluster++;
            }
          else j = nodeid[i];
          k = tree[i].left;
          if (k<0) nodeid[-k-1] = j; else clusterid[k] = j;
          k = tree[i].right;
          if (k<0) nodeid[-k-1] = j; else clusterid[k] = j;
        }
      free(nodeid);
      return;
}
/******************************************************************************/
/****usage*/
void usage()
{
      printf("\n Usage:cc2clusterHierarchical ccfile nevent ccthreshold method[m|a|s]\n\n");
      printf(" This program using Hierarchical Clustering Method to find");
      printf(" similar earthquakes in a whole data set\n\n");
      printf(" ccfile is the lower trianglar matrix with the cc value\n");
      printf(" nevent is the number of event to be clustered\n");
      printf(" ccthreshold is the threshold to define the similar eartquake,such as 0.8\n");
      printf(" method is the method used to clustering.m or a or s only,see below.\n\n");
      printf(" The distance matrix is defined as dist(a,b)=1-corr(a,b)\n");
      printf(" Three different cross-subcluster-distance calculation:\n");
      printf(" Maximum linkage pairwise method: method=m,the defaults\n");
      printf("      If we want each pair events within one cluster have cc value >threshold\n");
      printf("      Use the maximum method to perfom the clustering.\n");
      printf(" Single  linkage pairwise method: method=s.\n");
      printf(" Average linkage pairwise method: method=a.\n\n");
      printf(" Output format:\n");
      printf(" event eventid cluster clusterid\n\n");
      printf(" Later,use\n");
      printf(" ./exe ccmatrixfile nevents threshold method|paste - file.list|getcluster.awk\n");
      printf(" To get the input for plotcluster.awk to create the figure by GMT\n");
      
}
/******************************************************************************/

/**************************NOW the Main ***************************************/


//#include"mine.h"
int main(int argc,char **argv)
{ 
    FILE *ccfile;
    double **imatrix;
    int nevent,ncluster;
    double ccvalue,mindist;
    int i,j;
    int count=0;
    int *clusterid;
    Node* tree;
    char method;
    //method='m';
/************************check arguments*************************/
    if(argc!=5){   usage(); exit(-1);}
/************************initial arguments***********************/
    ccfile=fopen(argv[1],"r");
    if(ccfile==NULL){
       printf("open file %s error!\n",argv[1]);
       exit(-1);}
    nevent=atoi(argv[2]);
    ccvalue=atof(argv[3]);
   // printf("the threshold is %f\n",ccvalue);
    sscanf(argv[4],"%c",&method);
   // printf("the method flag is %c\n",method);
   //printf("Totally %d events will be classified!\n",nevent);
/*In the original Cluster3.0,the matrix is not one square matrix.But it does not 
  matter if we use a square matrix because the distance along the trace is awlways
  zero and never be used.*/
    imatrix=(double**)malloc(nevent*sizeof(double*));
    for(i=0;i<nevent;i++){
          imatrix[i]=(double*)malloc(sizeof(double)*(i+1));
            for(j=0;j<=i;j++){
                 fscanf(ccfile,"%lf",&imatrix[i][j]);
                 imatrix[i][j]=1.0-imatrix[i][j];
                }
        }
//test output
//   for(i=0;i<nevent;i++){
//      for(j=0;j<=i;j++){
//         printf("%lf\t",imatrix[i][j]);
//             }
//      printf("\n");
//  }
/************now do the cluster using different linkage method ******************/
    switch(method)
    {
       case 'm': //maximum linkage 
                tree=pmlcluster(nevent,imatrix);
                break;
       case 's': //single linkage method
                tree=pslcluster(nevent,imatrix);
                break;
       case 'a': //averaged  linkage method
                tree=palcluster(nevent,imatrix);
                break; 
       defaults :
            printf("The method should be and can only be 'm','s','a'\n");
            exit(-1);
    }
//tree=pmlcluster(nevent,distmatrix);
/*************Now output the results for test ********************************/
//printf("Node     Item 1   Item 2    Distance\n");
//for(i=0; i<nevent-1; i++) printf("%3d:%9d%9d      %g\n",-i-1, tree[i].left, tree[i].right, tree[i].distance);
//printf("\n");
/***********using the mindist to get parameter for cuttree*************/
     mindist=1-ccvalue;
//printf("ccvalue is %lf\n",ccvalue);
//printf("mindist is %lf\n",mindist);
    for(i=0;i<nevent-1;i++){
        if(tree[i].distance<mindist)  count++;
            }
//printf("Count is %d\n",count);
     ncluster=nevent-count;
//printf("ncluster is %d\n",ncluster);

//printf("=============== Cutting a hierarchical clustering tree ==========\n");
    clusterid =(int*)malloc(nevent*sizeof(int));
    cuttree (nevent, tree, ncluster, clusterid);  
/*************Now out put the final results************8***/
    for(i=0; i<nevent; i++) printf("eventid \t%4d\tclusterid\t%d\n",i+1, clusterid[i]);

/******deallocate memory*/
    free(tree);
    for(i=0;i<nevent;i++)  free(imatrix[i]);    
    free(imatrix);
    fclose(ccfile);
    return(1);
} 
