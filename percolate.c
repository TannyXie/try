/*
 * Simple parallel program to test for percolation of a cluster
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "percolate.h"
int rank_v(int rank, int size)  {
  return (rank >=0 ) && (rank < size);
}
/*
 * Simple parallel program to test for percolation of a cluster.
 */

int update(int old[][N+2], int new[][N+2]) {
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

          if (newval != oldval){
            ++nchangelocal;
          }
        }

        new[i][j] = newval;
      }
  }
  return nchangelocal;
}

int main(int argc, char *argv[])
{
  /*
   *  Define the main arrays for the simulation
   */

  int old[M+2][N+2] = {0};
  int new[M+2][N+2] = {0};

  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */

  int map[L][L] = {0};

  /*
   *  Array to store local part of map
   */

  int smallmap[M][N] = {0};  

  /*
   *  Variables that define the simulation
   */

  int seed;
  double rho;

  /*
   *  Local variables
   */

  int i, j, k, w, nhole, step, maxstep, oldval, newval;
  int nchangelocal, nchange, printfreq;
  int itop, ibot, perc;
  double r;

  /*
   *  MPI variables
   */

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Status status;

  int size, rank, left, right, up, down;
  int tag = 1;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
// TODO:
  left = (rank / NPROC)*NPROC + (rank - 1 + NPROC) % NPROC;
  right = (rank / NPROC)*NPROC + (rank + 1 + NPROC) % NPROC;
  down = rank + NPROC;
  up = rank - NPROC;
  printf("size: %d\n", size);
  if(rank == 1) {
    printf("rank has down process rank %d\n", down);
  }

  /*
   * Non-periodic boundary conditions
   *
   * Note that the special rank of MPI_PROC_NULL is a "black hole" for
   * communications. Using this value for processes off the edges of the
   * image means there is no additional logic needed to ensure processes
   * at the edges do not attempt to send to or receive from invalid
   * ranks (i.e. rank = -1 and rank = NPROC).
   *
   * Proper solution would compute neighbours with a Cartesian topology
   * and MPI_Cart_shift, where MPI_PROC_NULL is assigned automatically.
   */

  if (!rank_valid(left, size))
    {
      left = MPI_PROC_NULL;
    }
  if (!rank_valid(right, size))
    {
      right = MPI_PROC_NULL;
    }
  if (!rank_v(down, size))
    {
      down = MPI_PROC_NULL;
    }
  if (!rank_v(up, size))
    {
      up = MPI_PROC_NULL;
    }

  if (NPROC * MPROC != size)
    {
      if (rank == 0)
      {
        printf("percolate: ERROR, NPROC = %d but running on %d\n",
        NPROC, size);
      }

      MPI_Finalize();
      return 0;
  }

  if (argc != 2)
    {
      if (rank == 0)
      {
        printf("Usage: percolate <seed>\n");
      }

      MPI_Finalize();
      return 0;
    }

  /*
   *  Update for a fixed number of steps, periodically report progress
   */

  maxstep = 5*L*L;
  printfreq = 10; 
  printf("Before Init\n");

  if (rank == 0) {
    printf("percolate: running on %d process(es)\n", size);

    /*
      *  Set most important value: the rock density rho (between 0 and 1)
      */

    rho = 0.4064;

    /*
      *  Set the randum number seed and initialise the generator
      */

    seed = atoi(argv[1]);

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

  MPI_Barrier(comm);

  /*
   *  Now scatter map to smallmap
   */

  // TODO:
  //MPI_Scatter(map, M*N, MPI_INT, smallmap, M*N, MPI_INT, 0, comm);
  /**
   * For MPROC * NPROC processes, the rank 0 process needs to scatter the map to 
   * them. So there would be two loops, range from 0 to MPROC, and from 0 to NPROC.
   */

  if (rank != 0) {
    for (i = 0; i < M; ++i) {
      MPI_Recv(&smallmap[i][0], N, MPI_INT, 0, 0, comm, &status);
      //int a[N];
      //MPI_Recv(a, N, MPI_INT, 0, 0, comm, &status);
    }
    //printf("Rank %d recv over\n", rank);
  }
  else {
    for(i = 0; i < MPROC; ++i) {
      for(j = 0; j < NPROC; ++j) {
        if (i == 0 && j == 0) continue;
        for(k = 0; k < M; ++k) {
          //int a[N] = {0};
          MPI_Ssend(&map[k+i*M][j*N], N, MPI_INT, i*MPROC+j, 0, comm);
          //MPI_Ssend(a, N, MPI_INT, i*MPROC+j, 0, comm);
        }
        printf("Rank %d send over\n", i*MPROC+j);
      }
    }
    for(i = 0; i < M; ++i) {
      for(j = 0; j < N; ++j) {
        smallmap[i][j] = map[i][j];
      }
    }
  }

  if(rank == 1){
    printf("In rank 1, the small map looks like\n");
    for(int i=0; i < M; ++i) {
      for(int j=0; j < N; ++j) {
        printf("%2d ", smallmap[i][j]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);

  if(rank == 2){
    printf("In rank 2, the small map looks like\n");
    for(int i=0; i < M; ++i) {
      for(int j=0; j < N; ++j) {
        printf("%2d ", smallmap[i][j]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);

  if(rank == 3){
    printf("In rank 3, the small map looks like\n");
    for(int i=0; i < M; ++i) {
      for(int j=0; j < N; ++j) {
        printf("%2d ", smallmap[i][j]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);

  if(rank == 0){
    printf("In rank 0, the small map looks like\n");
    for(int i=0; i < M; ++i) {
      for(int j=0; j < N; ++j) {
        printf("%2d ", smallmap[i][j]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);
  printf("This is sync1 over from rank %d\n", rank);


  /*
   * Initialise the old array: copy the LxL array smallmap to the centre of
   * old, and set the halo values to zero.
   */

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

  step = 1;
  nchange = 1;

  while (step <= maxstep && nchange != 0) {
    /*
      *  Swap halos up and down
      */

    /*
      * Communications is done using the sendrecv routine; a proper
      * solution would use non-blocking communications (e.g. some
      * combination of issend/recv or ssend/irecv)
      */
    // TODO:
    MPI_Request request_s, request_r;
    
    MPI_Sendrecv(&old[M][1], N, MPI_INT, down, tag,
      &old[0][1], N, MPI_INT, up, tag,
      comm, &status);
    MPI_Sendrecv(&old[1][1], N, MPI_INT, up, tag,
      &old[M+1][1], N, MPI_INT, down, tag,
      comm, &status);
    MPI_Barrier(comm);
    if(rank == 0 && step == 1)
      printf("This is sync21 over\n");
    MPI_Issend(&old[1][1], N, MPI_INT, down, tag, comm, &request_s);
    MPI_Irecv(&old[M+1][1], N, MPI_INT, up, tag, comm, &request_r);
    MPI_Wait(&request_s, &status);
    MPI_Wait(&request_r, &status);
    MPI_Barrier(comm);
    if(rank == 0 && step == 1)
      printf("This is sync22 over\n");
    int temp_send_1[M], temp_send_N[M], temp_recv_0[M], temp_recv_Np1[M];
    for(i = 0; i < M; ++i) {
      temp_send_1[i] = old[i+1][1];
      temp_send_N[i] = old[i+1][N];
    }
    MPI_Issend(temp_send_1, M, MPI_INT, left, tag, comm, &request_s);
    MPI_Irecv(temp_recv_Np1, N, MPI_INT, right, tag, comm, &request_r);
    MPI_Wait(&request_s, &status);
    MPI_Wait(&request_r, &status);
    MPI_Barrier(comm);
    if(rank == 0 && step == 1)
      printf("This is sync23 over\n");
    MPI_Issend(temp_send_N, M, MPI_INT, right, tag, comm, &request_s);
    MPI_Irecv(temp_recv_0, N, MPI_INT, left, tag, comm, &request_r);
    MPI_Wait(&request_s, &status);
    MPI_Wait(&request_r, &status);
    MPI_Barrier(comm);
    if(rank == 0 && step == 1)
      printf("This is sync24 over\n");
    
    for(i = 0; i < M; ++i){
      old[i+1][0] = temp_recv_0[i];
      old[i+1][N+1] = temp_recv_Np1[i];
    }

    nchangelocal = 0;
  if(rank == 1){
    printf("In rank 1, after receiving data in step: %d, the small map looks like\n", step);
    for(int i=0; i < M+2; ++i) {
      for(int j=0; j < N+2; ++j) {
        printf("%2d ", old[i][j]);
      }
      printf("\n");
    }
  }

    nchangelocal = update(old, new);
      printf("rank %d, nchangelocal : %d\n", rank, nchangelocal);
    
 /*
  if(rank == 1){
    printf("In rank 1, after updating, in step: %d, the small map looks like\n", step);
    for(int i=0; i < M+2; ++i) {
      for(int j=0; j < N+2; ++j) {
        printf("%2d ", old[i][j]);
      }
      printf("\n");
    }
  }
  */
    /*
      *  Compute global number of changes on rank 0
      */

    MPI_Allreduce(&nchangelocal, &nchange, 1, MPI_INT, MPI_SUM, comm);
    /*
      *  Report progress every now and then
      */

    if (step % printfreq == 0)
    {
      if (rank == 0)
        {
                printf("percolate: changes on step %d is %d\n",
                      step, nchange);
        }
    }

    /*
      *  Copy back in preparation for next step, omitting halos
      */

    for (i=1; i<=M; i++)
    {
      for (j=1; j<=N; j++)
        {
          old[i][j] = new[i][j];
        }
    }
    /*
printf("---------step:%d---------\n", step);
  if(rank == 0){
    printf("In rank 0, step: %d, the small map looks like\n", step);
    for(int i=0; i < M; ++i) {
      for(int j=0; j < N; ++j) {
        printf("%2d ", old[i-1][j-1]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);
  
  if(rank == 1){
    printf("In rank 1, step: %d, the small map looks like\n", step);
    for(int i=0; i < M; ++i) {
      for(int j=0; j < N; ++j) {
        printf("%2d ", old[i-1][j-1]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);

  if(rank == 2){
    printf("In rank 2, step: %d, the small map looks like\n", step);
    for(int i=0; i < M; ++i) {
      for(int j=0; j < N; ++j) {
        printf("%2d ", old[i-1][j-1]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);

  if(rank == 3){
    printf("In rank 3, step: %d, the small map looks like\n", step);
    for(int i=0; i < M; ++i) {
      for(int j=0; j < N; ++j) {
        printf("%2d ", old[i-1][j-1]);
      }
      printf("\n");
    }
  }
  MPI_Barrier(comm);
*/

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

  if (rank == 0)
    {
      if (nchange != 0)
      {
              printf("percolate: WARNING max steps = %d reached but nchange != 0\n",
          maxstep);
      }
    }

  /*
   *  Copy the centre of old, excluding the halos, into smallmap
   */
  
  for (i=1; i<=M; i++)
    {
      for (j=1; j<=N; j++)
      {
        smallmap[i-1][j-1] = old[i][j];
      }
    }
  /*
   *  Now gather smallmap back to map
   */
  // TODO:
  //MPI_Gather(smallmap, M*N, MPI_INT, map, M*N, MPI_INT, 0, comm);
  for(i=0; i<L; ++i ) {
    for(j=0;j<N;++j){
      map[i][j] = 0;
    }
  }
  if(rank == 0) {
    for(i = 0; i < MPROC; ++i) {
      for(j = 0; j < NPROC; ++j) {
        // skip rank 0 itself
        if (i == 0 && j == 0) continue;
        for(k = 0; k < M; ++k) {
          MPI_Recv(&map[k+i*M][j*N], N, MPI_INT, i*NPROC+j, tag, comm, &status);
        }
        printf("After receiving from process %d\n", i*NPROC+j);
        for(k = 0; k < L; ++k) {
          for(int p = 0; p < L; ++p)
            printf("%3d ", map[k][p]);
          printf("\n");
        }
      }
    }
    for(i = 0; i < M; ++i) {
      for(j = 0; j < N; ++j) {
        map[i][j] = smallmap[i][j];
      }
    }
  }
  else {
    for (int i = 0; i < M; ++i) {
      MPI_Ssend(&smallmap[i][0], N, MPI_INT, 0, tag, comm);
    }
  }
  MPI_Barrier(comm);
  if(rank == 0)
    printf("This is sync4 over\n\n");

  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */

  if (rank == 0){
    perc = 0;

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
    int map_t[L][L] = {0};

    mapwrite("map.pgm", map, 2);
  }

  MPI_Finalize();

  return 0;
}
