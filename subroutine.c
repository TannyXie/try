#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <string.h>
#include <math.h>

#include "percolate.h"
#include "arralloc.h"


extern int L,M,N;
extern int MPROC, NPROC;
enum direction {LEFT, RIGHT, UP, DOWN};

int getNonperiodicLeft(int coord[2]);

int MPI_Get_info(int* s, int* r, int coord[], int neighbours[], int mesh[], int periodic[], int boundaryflag[], MPI_Comm comm, MPI_Comm* newcomm) {

  MPI_Comm_size(comm, s);
  MPI_Comm_rank(comm, r);
  mesh[0] = mesh[1] = 0;
  int periodical[] = {1, 0};
  MPI_Dims_create(*s, 2, mesh);
  MPI_Cart_create(comm, 2, mesh, periodical, 1, newcomm);
  comm = *newcomm;
  MPI_Comm_size(comm, s);
  MPI_Comm_rank(comm, r);
  MPI_Cart_coords(comm, *r, 2, coord);
  MPI_Cart_shift(comm, 1, 1, &neighbours[LEFT], &neighbours[RIGHT]);
  MPI_Cart_shift(comm, 0, 1, &neighbours[UP], &neighbours[DOWN]);
  int rank = *r;
  int size = *s;
  MPROC = mesh[0];
  NPROC = mesh[1];
  M = L / MPROC;
  N = L / NPROC;

  if(rank == 1) {
    printf("rank has down process rank %d\n", neighbours[DOWN]);
  }


  int temp;
  periodic[0] = getNonperiodicLeft(coord);
  temp = coord[1];
  coord[1] = NPROC - 1 - coord[1];
  periodic[1] = getNonperiodicLeft(coord);
  coord[1] = temp;

  boundaryflag[0] = (coord[0] == 0);
  boundaryflag[1] = (coord[0] == MPROC - 1);

  printf("rank : %d, left_nonperiodic:%d, right_nonperiodic: %d\n", rank, periodic[0], periodic[1]);

  if (NPROC * MPROC != size) {
    if (rank == 0) {
      printf("percolate: ERROR, NPROC = %d but running on %d\n",
      NPROC, size);
    }

    return 1;
  }
  return 0;
}

void initializeMap(int** map, int seed, double rho, int size, int maxstep) {
  double r;
  int i, j, nhole;
  printf("percolate: running on %d process(es)\n", size);


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
  /*
  for( i = 0; i < L; ++i) {
    for(j = 0; j < L; ++j) {
      printf("%3d ", map[i][j]);
    }
    printf("\n");
  }
  */

  printf("percolate: rho = %f, actual density = %f\n",
    rho, 1.0 - ((double) nhole)/((double) L*L) );
}

/**
 * Map is initialized only in rank 0 and rank 0 needs to scatter is
 * to other ranks.
 **/
void scatterMap(int** smallmap, int** map, int coord[2], MPI_Comm comm) {
  MPI_Bcast(&map[0][0], L * L, MPI_INT, 0, comm);
  for(int i=0; i < M; ++i) {
    for(int j=0; j < N; ++j){
      int coordx = coord[0], coordy = coord[1];
      smallmap[i][j] = map[coordx * M + i][coordy * N + j];
    }
  }
}

/**
 * Print the content of the old to debug.
 **/
void printOldMap(int** old, int rank, int step, MPI_Comm comm) {
  printf("In rank %d, after receiving data in step: %d, the small map looks like\n", rank, step);
  for(int i=0; i < M+2; ++i) {
    printf("line %d ", i);
    for(int j=0; j < N+2; ++j) {
      printf("%2d ", old[i][j]);
    }
    printf("\n");
  }
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

int getNonperiodicLeft(int coord[2]) {
  int left_nonperiodic = 0;
  int ratio = 6;
  if((coord[1] + 1) * N <= L/ratio){
    left_nonperiodic = N;
  }
  else if(coord[1] * N >= L/ratio){
    left_nonperiodic = 0;
  } else {
    left_nonperiodic = L/ratio - coord[1] * N;
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
void transmit(int rank, int** old, int upmost, int downmost, int left_nonperiodic, int right_nonperiodic, int up, int down, int left, int right, int tag, MPI_Comm comm) {
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

  if(upmost == 1){
    for(i = 0; i < left_nonperiodic; ++i) {
      old[0][i+1] = 0;
    }
    for(i = 0; i < right_nonperiodic; ++i) {
      old[0][M-i] = 0;
    }
  }
  else if(downmost == 1) {
    for(i = 0; i < left_nonperiodic; ++i) {
      old[M+1][i+1] = 0;
    }
    for(i = 0; i < right_nonperiodic; ++i) {
      old[M+1][M-i] = 0;
    }
  }
  //printOldMap(old, rank, -1, comm);
}


/**
 * Reduce the old maps to the 
 **/
void reduceOldMaps(int** old, int** map, int coord[], MPI_Comm comm) {
  int** maptp = arralloc(sizeof(int), 2, L, L);
  int i, j;
  for (i=0; i<L; i++) {
    for (j=0; j<L; j++) {
      maptp[i][j] = 0;
    }
  }
  
  for (i=0; i<M; i++) {
    for (j=0; j<N; j++) {
      maptp[i+coord[0]*M][j+coord[1]*N] = old[i+1][j+1];
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
  /*
  for( i = 0; i < L; ++i) {
    for(j = 0; j < L; ++j) {
      printf("%3d ", map[i][j]);
    }
    printf("\n");
  }
  */
  
  
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
