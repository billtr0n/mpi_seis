#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

/* reads several time series from contiguous sites (non-multiplexed) */
void read_awp_timeseries(char *fname, int nx, int ny, int nz, int xi, int yi, int zi, 
    int nt, int nsites, float *buf){
    MPI_File fh;
    MPI_Datatype filetype;
    MPI_Offset disp;

    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    MPI_Type_vector(nt, nsites, nx*ny*nz, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    disp = sizeof(MPI_FLOAT) * (nx*ny*zi + nx*yi + xi);
    MPI_File_set_view(fh, disp, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);

    MPI_File_read_all(fh, buf, nt*nsites, MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
}

/* reads several time series from contiguous sites (multiplexed) */
void read_awp_timeseries_multi(char *fbase, int nx, int ny, int nz, int xi, int yi, int zi, 
    int nt, int wstep, int nskip, int nsites, float *buf){
   int nfiles, k;
   MPI_Offset bpos;
   char fname[200];

   nfiles = nt / wstep;
   for (k=0; k<nfiles; k++){
      sprintf(fname, "%s%07d", fbase, (k+1)*wstep*nskip);
      bpos=k*wstep*nsites;
      read_awp_timeseries(fname, nx, ny, nz, xi, yi, zi, wstep, nsites, buf+bpos);
   }
};

/*
int main (int argc, char **argv){
    int rank, nprocs;
    float *buf, **x;
    int k, l;
    int nsites=20, nt=200;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    buf=(float*) calloc(nt*nsites, sizeof(float));
    read_awp_timeseries_multi("output_sfc/SX", 88, 88, 44, 67, 43, 0, 100, 
        100, 10, nsites, buf);

    x=(float**) calloc(nsites, sizeof(float*));
    for (l=0; l<nsites; l++) x[l] = (float*) calloc(nt, sizeof(float));

    for (k=0; k<nt; k++){
       for (l=0; l<nsites; l++){
          x[l][k]=buf[k*nsites+l];
       }
    }

    free(buf);
    for (k=0; k<nt; k++) fprintf(stdout, "%e %e\n", x[0][k], x[19][k]);
    return(0);

}*/
