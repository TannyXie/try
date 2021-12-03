/*
 * Simple parallel program to test for percolation of a cluster
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "percolate.h"
#include "arralloc.h"
#include "subroutine.h"

int L,M,N;
int MPROC, NPROC;
int main(int argc, char *argv[])
{
  clock_t start_total, end_total, start_calc, end_calc;
  start_total = clock();

  /*
   *  Variables that define the simulation, the initial value is default value
   */

  int seed = 8759, max_size = 2;
  double rho = 0.4064;
  L = 480;

  
  MPI_Init(&argc, &argv);


  int opt;
  int exit_code = -1;
  char* ptr;

  while ((opt = getopt(argc, argv, ":hd:p:l:s:m:r:")) != -1){
    switch (opt){
    case 's':
      seed = atoi(optarg);
      printf("Recognized seed : %d\n", seed);
      break;
    case 'm':
      max_size = atoi(optarg);
      printf("Recognized max_size : %d\n", max_size);
      break;
    case 'l':	
      L = atoi(optarg);
      printf("Recognized L : %d\n", L);
      break;
    case 'r':	
      rho = strtod(optarg, &ptr);
      printf("Recognized rho : %f\n", rho);
      break;
    case '?':
      printf("Unknown option: %c\n", optopt);
      exit_code = 1;
      break;
    }
  }
  for(; optind < argc; optind++){
    printf("Extra argument: %s\n", argv[optind]);
    exit_code = 1;
  }

  if(exit_code != -1) {
    MPI_Finalize();
    return 0;
  }

  /*
   *  Local variables
   */

  int i, j, k, w, step, maxstep, oldval, newval;
  int nchangelocal, nchange, printfreq;
  long nmaplocal, nmap;
  int itop, ibot, perc;

  /*
   *  MPI variables and Cartesian variables
   */

  MPI_Comm comm = MPI_COMM_WORLD, newcomm;
  MPI_Status status, statuses[8];

  int size, rank, left, right, up, down, temp;
  int coord[2], offset[2], neighbours[4], mesh[2], periodic[2], boundaryflag[2];
  int tag = 1;


  if (MPI_Get_info(&size, &rank, coord, offset, neighbours, mesh, periodic, boundaryflag, comm, &newcomm) == 1) {
    printf("Cannot create a cartesian\n");
    MPI_Finalize();
    return 0;
  }
  comm = newcomm;

  
  /*
   *  Define the main arrays for the simulation
   */

  int ** old = arralloc(sizeof(int), 2, M+2, N+2);
  int ** new = arralloc(sizeof(int), 2, M+2, N+2);


  int ** map = arralloc(sizeof(int), 2, L, L);
  int ** smallmap = arralloc(sizeof(int), 2, M, N);

  maxstep = 5*L*L;
  printfreq = 100;

  if (rank == 0) {
    initializeMap(map, seed, rho, size, maxstep);
  }


  scatterMap(smallmap, map, offset, comm);

  initializeOldMap(old, smallmap);


  step = 1;
  nchange = 1;

  start_calc = clock();

  while (step <= maxstep && nchange != 0) {
    transmit(rank, old, boundaryflag[0], boundaryflag[1], periodic[0], 
          periodic[1], neighbours[UP], neighbours[DOWN], neighbours[LEFT], 
          neighbours[RIGHT], tag, comm);

    /**
     * Update old map and store to new map.
     * Return the number of elements begin updated.
     **/
    nchangelocal = update(old, new);
    nmaplocal = 0;
    nmap = 0;
    for(i=0;i<M;++i){
      for(j=0;j<N;++j){
        nmaplocal += new[i+1][j+1];
      }
    }
    nmaplocal /= (M*N);

    end_calc = clock();
    
    /**
     * Compute global number of changes on rank 0
     **/
    MPI_Allreduce(&nchangelocal, &nchange, 1, MPI_INT, MPI_SUM, comm);
    MPI_Reduce(&nmaplocal, &nmap, 1, MPI_INT, MPI_SUM, 0, comm);
    /**
     * Report progress every now and then
     **/
    if (step % printfreq == 0 && rank == 0){
          printf("percolate: changes on step %d is %d, average is %f\n",
                step, nchange, ((double)nmap) / (MPROC*NPROC));
    }

    /**
     * Copy back in preparation for next step, omitting halos
     **/
    for (i=1; i<=M; i++){
      for (j=1; j<=N; j++){
          old[i][j] = new[i][j];
        }
    }

    step++;
  }
  if(rank == 0)
    printf("This is sync3 over from rank %d\n\n", rank);
  /*
   *  We set a maximum number of steps to ensure the algorithm always
   *  terminates. However, if we hit this limit before the algorithm
   *  has finished then there must have been a problem (e.g. the value
   *  of maxstep is too small)
   */

  
  if (rank == 0 && nchange != 0){
            printf("percolate: WARNING max steps = %d reached \
                 but nchange != 0\n",
        maxstep);
  }

  /*
   *  Copy the centre of old, excluding the halos, into smallmap
   */
  reduceOldMaps(old, map, offset, comm);
  
  /*
   *  Now gather smallmap back to map
   */
  
  if(rank == 0)
    printf("This is sync4 over\n\n");
  

  if (rank == 0){
    checkPercolate(map, max_size);
  }
  MPI_Barrier(comm);

  end_total = clock();

  if(rank == 0) {
    printf("The whole program part takes time of %f\nthe parallelled calculation part takes time %f\n", \
            (double)(end_total - start_total), (double)(end_calc - start_calc));
  }

  MPI_Finalize();
  return 0;
}
