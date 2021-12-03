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

int MPI_Get_info(int* s, int* r, int coord[], int neighbours[], int mesh[], int periodic[], int boundaryflag[], MPI_Comm comm, MPI_Comm* newcomm);
void initializeMap(int** map, int seed, double rho, int size, int maxstep);
void scatterMap(int** smallmap, int** map, int coord[2], MPI_Comm comm);
void printOldMap(int** old, int rank, int step, MPI_Comm comm);
int update(int** old, int** new);
int getNonperiodicLeft(int coord[2]);
void initializeOldMap(int** old, int** smallmap);
void transmit(int rank, int** old, int upmost, int downmost, int left_nonperiodic, int right_nonperiodic, int up, int down, int left, int right, int tag, MPI_Comm comm);
void reduceOldMaps(int** old, int** map, int coord[], MPI_Comm comm);
void checkPercolate(int** map, int max_size);