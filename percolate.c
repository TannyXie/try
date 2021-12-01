/*
 * Simple parallel program to test for percolation of a cluster
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <string.h>
#include <math.h>

#include "percolate.h"
#include "arralloc.h"

int L,M,N;
int MPROC, NPROC;
enum direction {LEFT, RIGHT, UP, DOWN};


void rank_v(int* up, int* down, int* left, int* right, int size)  {
  if (!rank_valid(*left, size))
    {
      *left = MPI_PROC_NULL;
    }
  if (!rank_valid(*right, size))
    {
      *right = MPI_PROC_NULL;
    }
  if (!rank_valid(*down, size))
    {
      *down = MPI_PROC_NULL;
    }
  if (!rank_valid(*up, size))
    {
      *up = MPI_PROC_NULL;
    }
}

int MPI_Get_info(int* s, int* r, int coord[], int neighbours[], int mesh[], int periodic[], int boundaryflag[], MPI_Comm comm) {
  MPI_Comm_size(comm, s);
  MPI_Comm_rank(comm, r);
  // TODO:
  mesh[0] = (int)(sqrt(*s));
  mesh[1] = (int)(*s / sqrt(*s));
  MPROC = mesh[0];
  NPROC = mesh[1];
  M = L / MPROC;
  N = L / NPROC;
  int rank = *r;
  int size = *s;
  coord[0] = rank / NPROC;
  coord[1] = rank % NPROC;
  if(rank == 0)
    printf("NPROC:%d, MPROC:%d\n", NPROC, MPROC);
  neighbours[LEFT] = (rank / NPROC)*NPROC + (rank - 1 + NPROC) % NPROC;
  neighbours[RIGHT] = (rank / NPROC)*NPROC + (rank + 1 + NPROC) % NPROC;
  neighbours[UP] = (rank - NPROC + size) % size;
  neighbours[DOWN] = (rank + NPROC + size) % size;
  printf("size: %d\n", size);
  if(rank == 1) {
    printf("rank has down process rank %d\n", neighbours[DOWN]);
  }

  int temp;
  periodic[0] = getNonperiodicLeft(rank);
  temp = rank;
  rank = NPROC - 1 - coord[1];
  periodic[1] = getNonperiodicLeft(rank);
  rank = temp;
  boundaryflag[0] = (rank / NPROC);
  boundaryflag[1] = (size-1-rank) / NPROC;
  printf("rank : %d, left_nonperiodic:%d, right_nonperiodic: %d\n", rank, periodic[0], periodic[1]);

  rank_v(&neighbours[0], &neighbours[1], &neighbours[2], &neighbours[3], size);

  if (NPROC * MPROC != size) {
    if (rank == 0) {
      printf("percolate: ERROR, NPROC = %d but running on %d\n",
      NPROC, size);
    }

    return 1;
  }
  return 0;
}

void initializeMap(int** map, int seed, int size, int maxstep) {
  double rho, r;
  int i, j, nhole;
  printf("percolate: running on %d process(es)\n", size);

  /*
    *  Set most important value: the rock density rho (between 0 and 1)
    */

  rho = 0.4064;

  /*
    *  Set the randum number seed and initialise the generator
    */

  printf("percolate: L = %d, L = %d, rho = %f, seed = %d, maxstep = %d\n",
    L, L, rho, seed, maxstep);

  rinit(seed);

  /*
    *  Initialise map with density rho. Zero indicates rock, a positive
    *  value indicates a hole. For the algorithm to work, all the holes
    *  must be initialised with a unique integer
    */

  nhole = 0;

  for (i=0; i < L; i++)
  {
    for (j=0; j < L; j++)
    {
      r=uni();
  
      if(r < rho)
      {
        map[i][j] = 0;
      }
      else
      {
        nhole++;
        map[i][j] = nhole;
      }
    }
  }
  printf("The initial map is:\n");
  for( i = 0; i < L; ++i) {
    for(j = 0; j < L; ++j) {
      printf("%3d ", map[i][j]);
    }
    printf("\n");
  }

  printf("percolate: rho = %f, actual density = %f\n",
    rho, 1.0 - ((double) nhole)/((double) L*L) );
}

/**
 * Map is initialized only in rank 0 and rank 0 needs to scatter is
 * to other ranks.
 **/
void scatterMap(int** smallmap, int** map, int rank, MPI_Comm comm) {
  MPI_Bcast(&map[0][0], L * L, MPI_INT, 0, comm);
  for(int i=0; i < M; ++i) {
    for(int j=0; j < N; ++j){
      int coordx = rank / NPROC, coordy = rank % NPROC;
      smallmap[i][j] = map[coordx * M + i][coordy * N + j];
    }
  }
}

/**
 * Print the content of the old to debug.
 **/
void printOldMap(int** old, int rank, int step, MPI_Comm comm) {
  if(rank == 0){
    printf("In rank 0, after receiving data in step: %d, the small map looks like\n", step);
    for(int i=0; i < M+2; ++i) {
      for(int j=0; j < N+2; ++j) {
        printf("%2d ", old[i][j]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);

  if(rank == 1){
    printf("In rank 1, after receiving data in step: %d, the small map looks like\n", step);
    for(int i=0; i < M+2; ++i) {
      for(int j=0; j < N+2; ++j) {
        printf("%2d ", old[i][j]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);

  if(rank == 2){
    printf("In rank 2, after receiving data in step: %d, the small map looks like\n", step);
    for(int i=0; i < M+2; ++i) {
      for(int j=0; j < N+2; ++j) {
        printf("%2d ", old[i][j]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);

  if(rank == 3){
    printf("In rank 3, after receiving data in step: %d, the small map looks like\n", step);
    for(int i=0; i < M+2; ++i) {
      for(int j=0; j < N+2; ++j) {
        printf("%2d ", old[i][j]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);
}

/**
 * update the whole map after receiving the data from its neighbours
 * updated result is in new map
 */
int update(int** old, int** new) {
  int newval;
  int oldval;
  int nchangelocal = 0;
  for (int i=1; i<=M; i++){
    for (int j=1; j<=N; j++)
      {
        oldval = old[i][j];
        newval = oldval;

        /*
        * Set new[i][j] to be the maximum value of old[i][j]
        * and its four nearest neighbours
        */

        if (oldval != 0)
        {
          if (old[i][j-1] > newval) newval = old[i][j-1];
          if (old[i][j+1] > newval) newval = old[i][j+1];
          if (old[i-1][j] > newval) newval = old[i-1][j];
          if (old[i+1][j] > newval) newval = old[i+1][j];

          if (newval != oldval) {
            ++nchangelocal;
          }
        }

        new[i][j] = newval;
      }
  }
  return nchangelocal;
}

int getNonperiodicLeft(int rank) {
  int left_nonperiodic;
  int ratio = 6;
  if((rank % NPROC + 1) * N <= L/ratio){
    left_nonperiodic = N;
  }
  else if((rank % NPROC) * N >= L/ratio){
    left_nonperiodic = 0;
  } else {
    left_nonperiodic = L/ratio - (rank % NPROC) * N;
  }
  return left_nonperiodic;
}

/**
 * Initialise the old array: copy the LxL array smallmap to the centre of
 * old, and set the halo values to zero.
 **/
void initializeOldMap(int** old, int** smallmap) {
  int i,j;
  for (i=1; i <= M; i++)
  {
    for (j=1; j <= N; j++)
    {
      old[i][j] = smallmap[i-1][j-1];
    }
  }

  for (i=0; i <= M+1; i++)  // zero the bottom and top halos
  {
    old[i][0]   = 0;
    old[i][N+1] = 0;
  }

  for (j=0; j <= N+1; j++)  // zero the left and right halos
  {
    old[0][j]   = 0;
    old[M+1][j] = 0;
  }
}

/**
 * Swap data with its neighbours
 **/
void transmit(int** old, int left_nonperiodic, int right_nonperiodic, int upmost, int downmost, int up, int down, int left, int right, int tag, MPI_Comm comm) {
  MPI_Request requests[8],request_s,request_r; // for send and receive for up down, left right
  MPI_Status statuses[8];
  MPI_Datatype column_type;
  MPI_Type_vector(M, 1, N + 2, MPI_INT, &column_type);
  MPI_Type_commit(&column_type);

  int i,j;

  MPI_Issend(&old[M][1], N, MPI_INT, down, tag, comm, &requests[0]);
  MPI_Irecv(&old[0][1], N, MPI_INT, up, tag, comm, &requests[1]);
  MPI_Issend(&old[1][1], N, MPI_INT, up, tag, comm, &requests[2]);
  MPI_Irecv(&old[M+1][1], N, MPI_INT, down, tag, comm, &requests[3]);

  MPI_Issend(&old[1][1], 1, column_type, left, tag, comm, &requests[4]);
  MPI_Irecv(&old[1][N+1], 1, column_type, right, tag, comm, &requests[5]);
  MPI_Issend(&old[1][N], 1, column_type, right, tag, comm, &requests[6]);
  MPI_Irecv(&old[1][0], 1, column_type, left, tag, comm, &requests[7]);
  
  MPI_Waitall(8, requests, statuses);


  if(upmost == 0){
    for(i = 0; i < left_nonperiodic; ++i) {
      old[0][i+1] = 0;
    }
    for(i = 0; i < right_nonperiodic; ++i) {
      old[0][M-i] = 0;
    }
  }
  else if(downmost == 0) {
    for(i = 0; i < left_nonperiodic; ++i) {
      old[M+1][i+1] = 0;
    }
    for(i = 0; i < right_nonperiodic; ++i) {
      old[M+1][M-i] = 0;
    }
  }
}

void reduceOldMaps(int** old, int** map, int rank, MPI_Comm comm) {
  int** maptp = arralloc(sizeof(int), 2, L, L);
  int i, j;
  for (i=0; i<M; i++) {
    for (j=0; j<N; j++) {
      maptp[i+rank/NPROC*M][j+rank%NPROC*N] = old[i+1][j+1];
    }
  }
  //this method can replace the following code block
  MPI_Reduce(&maptp[0][0], &map[0][0], L*L, MPI_INT, MPI_SUM, 0, comm);
}

/*
  *  Test to see if percolation occurred by looking for positive numbers
  *  that appear on both the top and bottom edges
  */
void checkPercolate(int** map, int max_size){
  int perc = 0, itop, ibot, i, j;

  for (itop=0; itop < L; itop++){
    if (map[itop][L-1] <= 0) continue;
    for (ibot=0; ibot < L; ibot++){
      if (map[ibot][0] == map[itop][L-1]){
        perc = 1;
      }
    }
  }
  printf("After the gather:\n");
  for( i = 0; i < L; ++i) {
    for(j = 0; j < L; ++j) {
      printf("%3d ", map[i][j]);
    }
    printf("\n");
  }

  if (perc != 0)
  {
    printf("percolate: cluster DOES percolate\n");
  }
  else
  {
    printf("percolate: cluster DOES NOT percolate\n");
  }

  /*
    *  Write the map to the file "map.pgm", displaying the two
    *  largest clusters. If the last argument here was 3, it would
    *  display the three largest clusters etc. The picture looks
    *  cleanest with only a single cluster, but multiple clusters
    *  are useful for debugging.
    */

  //mapwrite("map.pgm", map, max_size);
  mapwritedynamic("map.pgm", map, L, max_size);
}

int main(int argc, char *argv[])
{

  /*
   *  Variables that define the simulation
   */

  int seed = 2, max_size = 2;
  L = 8;

  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */
  MPI_Init(&argc, &argv);



//seed be set to a initial value

  int opt;
  int exit_code = -1;

  while ((opt = getopt(argc, argv, ":hd:p:l:s:m:r:")) != -1){
    switch (opt){
    case 's':
      seed = atoi(optarg);
      break;
    case 'm':
      max_size = atoi(optarg);
      break;
    case 'l':	
      L = atoi(optarg);
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

  /*
   *  Local variables
   */

  int i, j, k, w, step, maxstep, oldval, newval;
  int nchangelocal, nchange, printfreq;
  int itop, ibot, perc;

  /*
   *  MPI variables
   */

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Status status, statuses[8];

  int size, rank, left, right, up, down, temp;
  int coord[2], neighbours[4], mesh[2], periodic[2], boundaryflag[2];
  int tag = 1;


  if (MPI_Get_info(&size, &rank, coord, neighbours, mesh, periodic, boundaryflag, comm) == 1) {
    MPI_Finalize();
    return 0;
  }
  if(rank == 0) {
    for(int i=0; i < 4; ++i) {
      printf("neighbour %d: %d, ", i, neighbours[i]);
    }
    printf("\n");
  }
  
  /*
   *  Define the main arrays for the simulation
   */

  int ** old = arralloc(sizeof(int), 2, M+2, N+2);
  int ** new = arralloc(sizeof(int), 2, M+2, N+2);


  int ** map = arralloc(sizeof(int), 2, L, L);
  int ** smallmap = arralloc(sizeof(int), 2, M, N);

  maxstep = 5*L*L;
  printfreq = 10; 

  if (rank == 0) {
    initializeMap(map, atoi(argv[1]), size, maxstep);
  }

  MPI_Barrier(comm);

  // TODO:
  scatterMap(smallmap, map, rank, comm);

  initializeOldMap(old, smallmap);

  printOldMap(old, rank, -1, comm);

  step = 1;
  nchange = 1;

  while (step <= maxstep && nchange != 0) {
    // TODO:
    transmit(old, boundaryflag[0], boundaryflag[1], periodic[0], 
          periodic[1], neighbours[UP], neighbours[DOWN], neighbours[LEFT], 
          neighbours[RIGHT], tag, comm);
    printOldMap(old, rank, step, comm);

    /**
     * Update old map and store to new map.
     * Return the number of elements begin updated.
     **/
    nchangelocal = update(old, new);
    printf("rank %d, nchangelocal : %d\n", rank, nchangelocal);
    
    /**
     * Compute global number of changes on rank 0
     **/
    MPI_Allreduce(&nchangelocal, &nchange, 1, MPI_INT, MPI_SUM, comm);
    /**
     * Report progress every now and then
     **/
    if (step % printfreq == 0 && rank == 0){
          printf("percolate: changes on step %d is %d, average is %f\n",
                step, nchange, ((float)nchange)/(L*L));
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
  MPI_Barrier(MPI_COMM_WORLD);
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
  reduceOldMaps(old, map, rank, comm);
  
  /*
   *  Now gather smallmap back to map
   */
  // TODO:
  MPI_Barrier(comm);
  if(rank == 0)
    printf("This is sync4 over\n\n");

  if (rank == 0){
    checkPercolate(map, max_size);
  }

  MPI_Finalize();
  return 0;
}
