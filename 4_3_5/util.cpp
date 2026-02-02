#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> 
#include <math.h>
#include "util.h"
#include "bnsch.h"

#define NR_END 1

/*===============
 * max |A - B|
 * ===============*/

float max_norm(float ***a, float ***b,
        int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j,k;
    float x = 0.0;
    
    for(i = xl; i <= xr; i++)
        for(j = yl; j <= yr; j++)
            for(k = zl; k <=zr; k++){
        
        if (fabs(a[i][j][k] - b[i][j][k]) > x)
            x = fabs(a[i][j][k] - b[i][j][k]);
            }
    return x;
}



float ***cube(int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i,j,nrow=xr-xl+1, ncol=yr-yl+1, ndep=zr-zl+1;
    float ***t;
    
    t = (float ***) malloc(((nrow+1)*sizeof(float**)));
    t += 1;
    t -= xl;
    
    t[xl]=(float **) malloc(((nrow*ncol+1)*sizeof(float *)));
    t[xl] += 1;
    t[xl] -= yl;
    
    t[xl][yl] = (float *) malloc(((nrow*ncol*ndep+1)*sizeof(float)));
    t[xl][yl] += 1;
    t[xl][yl] -= zl;
    
    for(j=yl+1; j<=yr; j++)
        t[xl][j]=t[xl][j-1]+ndep;
    
    for(i=xl+1; i<=xr; i++)
    {
        t[i] = t[i-1]+ncol;
        t[i][yl] = t[i-1][yl]+ncol*ndep;
        for(j=yl+1; j<=yr; j++)
            t[i][j] = t[i][j-1]+ndep;
    }
    
    return t;
}


void free_cube(float ***t, int xl, int xr, int yl, int yr, int zl, int zr)
{
    free((char *) (t[xl][yl]+zl-1));
    free((char *) (t[xl]+yl-1));
    free((char *) (t+xl-1));
}



void zero_cube(float ***a, int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++)
        for (j=yl; j<=yr; j++)
            for (k=zl; k<=zr; k++){
        
        a[i][j][k] = 0.0;
        
            }
    
    return;
}




void cube_add(float ***a, float ***b, float ***c,
        int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++)
        for (j=yl; j<=yr; j++)
            for (k=zl; k<=zr; k++){
        
        a[i][j][k]=b[i][j][k]+c[i][j][k];
        
            }
    
    return;
}

void cube_add2(float ***a, float ***b, float ***c,
        float ***a2, float ***b2, float ***c2,
        int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++)
        for (j=yl; j<=yr; j++)
            for (k=zl; k<=zr; k++){
        
        a[i][j][k]=b[i][j][k]+c[i][j][k];
        a2[i][j][k]=b2[i][j][k]+c2[i][j][k];
        
            }
    
    return;
}








void cube_sub(float ***a, float ***b, float ***c,
        int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++)
        for (j=yl; j<=yr; j++)
            for (k=zl; k<=zr; k++){
        
        a[i][j][k]=b[i][j][k]-c[i][j][k];
        
            }
    
    return;
}


void cube_sub2(float ***a, float ***b, float ***c,
        float ***a2, float ***b2, float ***c2,
        int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++)
        for (j=yl; j<=yr; j++)
            for (k=zl; k<=zr; k++){
        
        a[i][j][k]=b[i][j][k]-c[i][j][k];
        a2[i][j][k]=b2[i][j][k]-c2[i][j][k];
        
            }
    
    return;
}



void cube_copy(float ***a, float ***b,
        int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++)
        for (j=yl; j<=yr; j++)
            for (k=zl; k<=zr; k++){
        
        a[i][j][k]=b[i][j][k];
        
            }
    
    return;
}


void cube_copy2(float  ***a, float  ***b,
        float  ***a2, float  ***b2,
        int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    
    for (i=xl; i<=xr; i++)
        for (j=yl; j<=yr; j++)
            for (k=zl; k<=zr; k++){
        
        a[i][j][k]=b[i][j][k];
        a2[i][j][k]=b2[i][j][k];
        
            }
    
    return;
}


float cube_max(float ***a,
        int xl, int xr, int yl, int yr, int zl, int zr)
{
    int i, j, k;
    float x = 0.0;
    
    for (i=xl; i<=xr; i++)
        for (j=yl; j<=yr; j++)
            for (k=zl; k<=zr; k++){
        
        if (fabs(a[i][j][k]) > x)
            x = fabs(a[i][j][k]);
            }
    
    return x;
}


void print_cube(FILE *fptr,
        float ***a, int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
    int i, j, k;
    
    for(k = ndl; k <= ndh; k++) {
        for(j = ncl; j <= nch; j++) {
            for(i = nrl; i <= nrh; i++) {
                fprintf(fptr, " %f", a[i][j][k]);
                
                fprintf(fptr, "\n");
            }}}
    
    return;
    
    
}


void print_data_sf( float ***fx, float ***fy,float ***fz,int kk, int index)
{
    extern int nx, ny, nz;
    
    char bufferfx[100], bufferfy[100],bufferfz[100];
    
    FILE  *sfx, *sfy, *sfz;
    
    sprintf(bufferfx,"data%d/sfx%d.m",index, kk);
    sprintf(bufferfy,"data%d/sfy%d.m",index, kk);
    sprintf(bufferfz,"data%d/sfz%d.m",index, kk);
    
    
    sfx = fopen(bufferfx,"w");
    sfy = fopen(bufferfy,"w");
    sfz = fopen(bufferfz,"w");
    
    print_cube(sfx, fx, 0, nx, 1, ny, 1, nz);
    print_cube(sfy, fy, 1, nx, 0, ny, 1, nz);
    print_cube(sfz, fz, 1, nx, 1, ny, 0, nz);
    
    fclose(sfx);
    fclose(sfy);
    fclose(sfz);
    
    return;
}


void print_data(float ***phi,int kk, int index)
{
    int i,j,k;
    char buffer[20];
    FILE *fc;
    sprintf(buffer,"data%d/phi%d.m",index, kk);
    fc = fopen(buffer,"w");
    ijkloop
            fprintf(fc, " %10.7f\n", phi[i][j][k]);
    fclose(fc);
    
    return;
}



