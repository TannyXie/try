/*
 *  Main header file for percolation code.
 */

/*
 *  System size L L
 */

#define L 8


/*
 *  Use 1D decomposition over NPROC processes across first dimension
 *  For an LxL simulation, the local arrays are of size MxN
 */

// TODO:
#define NPROC 2
#define MPROC 2

#define M L/MPROC
#define N L/NPROC

/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */

void mapwrite(char *percfile, int map[L][L], int ncluster);
void mapwritedynamic(char *percfile, int **map, int l, int ncluster);

/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);

/**
 * utils.c
 */
int rank_valid(int rank, int size);