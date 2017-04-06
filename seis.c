/* Computes peak ground velocity from AWP-ODC output.
   Improved version using MPI-IO.
   Daniel Roten, March 2015 <droten@sdsc.edu> */

#include "awp_data.h"
#include "fd3dparam.h"
#include "xapiir.h"
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char *argv[] ){
   FD3D_param par;
   read_settings(&par, "IN3D.out");
   char *xfile, *yfile, *zfile;
   int nx=1, ny=1, nt=1, nz=1;
   float **x, **y, **z;
   float *bufx, *bufy, *bufz;
   float *PH, *PZ, vh, vz;
   float *DX, *DY, *DZ;
   MPI_File xfid, yfid, zfid;
   int nchunks = 88, csize;
   int rank, nprocs;
   int k, l, n;
   MPI_Offset off;
   int s0, xi, yi, zi;

   /*parameters for xapiir*/
   int iord=6, npas=1;
   float trbndw=0., a=0.; /* chebyshev parameters */
   char *aproto="BU";
   char *ftype="LP";
   float hp=0., lp=15.00; /*lp=1.0*/

   float dt2;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   /* determine dimensions from these parameters */
   nx=(int) floorf( (par.nedx-par.nbgx+1)/par.nskpx); 
   ny=(int) floorf( (par.nedy-par.nbgy+1)/par.nskpy); 
   nz=(int) floorf( (par.nedz-par.nbgz+1)/par.nskpz); 
   nt=(int) floorf( (par.tmax/par.dt/par.ntiskp));

   
   /* output parameters for debug */
   fprintf(stdout, "  >> nx=%d, ny=%d, nz=%d, nt=%d," 
                   " ivelocity=%d, ntiskp=%d, writestep=%d\n",
		    nx, ny, nz, nt, par.itype, par.ntiskp, par.writestep);

   xfile=(char*) calloc(100, sizeof(char) );
   yfile=(char*) calloc(100, sizeof(char) );
   zfile=(char*) calloc(100, sizeof(char) );
   sscanf(par.sxrgo, "%s", xfile);
   sscanf(par.syrgo, "%s", yfile);
   sscanf(par.szrgo, "%s", zfile);
   /* this removes the remaining ' at the end of the string, if any */
   xfile=strsep(&xfile, "\'");
   yfile=strsep(&yfile, "\'");
   zfile=strsep(&zfile, "\'");

   /* number of points to process per read operation */
   if (((nx*ny*nz) % (nprocs * nchunks)) != 0) {
      fprintf(stdout, "number of points not divisible by number of cpus * buffer size\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
   }
 
   csize = nx*ny*nz / nprocs / nchunks;

   /* buffers for unsorted velocities */
   bufx = (float*) calloc(csize*nt, sizeof(float));
   bufy = (float*) calloc(csize*nt, sizeof(float));
   bufz = (float*) calloc(csize*nt, sizeof(float));

   /* 2D arrays for velocities */
   x = (float**) calloc (csize, sizeof(float*));
   y = (float**) calloc (csize, sizeof(float*));
   z = (float**) calloc (csize, sizeof(float*));

   for (l=0; l<csize; l++){
      x[l] = (float*) calloc(nt, sizeof(float));
      y[l] = (float*) calloc(nt, sizeof(float));
      z[l] = (float*) calloc(nt, sizeof(float));
   }
   if (z[csize-1]==NULL) perror("error using calloc()");

   MPI_File_open(MPI_COMM_WORLD, "x_velocity.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &xfid);
   MPI_File_open(MPI_COMM_WORLD, "y_velocity.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &yfid);
   MPI_File_open(MPI_COMM_WORLD, "z_velocity.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, 
       MPI_INFO_NULL, &zfid);

   for (k=0; k<nchunks; k++){
      s0=rank*nchunks*csize + k*csize;
      zi=s0 / (nx*ny);
      yi=(s0 % (nx*ny)) / nx;
      xi=s0 % nx;
      //fprintf(stdout, "(%d): s0=%d, xi=%d, yi=%d, zi=%d\n", rank, s0, xi, yi, zi);
      read_awp_timeseries_multi(xfile, nx, ny, nz, xi, yi, zi, nt, 
         par.writestep, par.ntiskp, csize, bufx);
      read_awp_timeseries_multi(yfile, nx, ny, nz, xi, yi, zi, nt, 
         par.writestep, par.ntiskp, csize, bufy);
      read_awp_timeseries_multi(zfile, nx, ny, nz, xi, yi, zi, nt, 
         par.writestep, par.ntiskp, csize, bufz);

      MPI_Barrier(MPI_COMM_WORLD);

      for (n=0; n<nt; n++){
         for (l=0; l<csize; l++){
            x[l][n]=bufx[n*csize+l];
            y[l][n]=bufy[n*csize+l];
            z[l][n]=bufz[n*csize+l];
         }
         //if (k==0) fprintf(stdout, "%e\n", x[3851][n]);
         //if (k==1) fprintf(stdout, "%e\n", x[3871][n]);
         //if (k==87) fprintf(stdout, "%e\n", x[87][n]);
      }

      for (l=0; l<csize; l++) {
         /*apply a filter here if desired */
         dt2=par.dt*par.ntiskp;
         xapiir_(x[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
         xapiir_(y[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);
         xapiir_(z[l], &nt, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt2, &npas);

         /* write seismograms */
         MPI_Barrier(MPI_COMM_WORLD);
         off = (MPI_Offset) nt * s0 * sizeof(float);
         MPI_File_write_at_all(xfid, off, PH, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
         MPI_File_write_at_all(yfid, off, PZ, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
         MPI_File_write_at_all(zfid, off, DX, csize, MPI_FLOAT, MPI_STATUS_IGNORE);
      }
   }
   MPI_File_close(&hfid); 
   MPI_File_close(&zfid); 
   MPI_File_close(&dxfid); 
   MPI_File_close(&dyfid); 
   MPI_File_close(&dzfid); 

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(0);
}

