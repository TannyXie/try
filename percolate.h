/*
 *  Main header file for percolation code.
 */

/*
 *  System size L
 */

#define LM 60
#define LN 40


/*
 *  Use 1D decomposition over NPROC processes across first dimension
 *  For an LxL simulation, the local arrays are of size MxN
 */

// TODO:
#define NPROC 4
#define MPROC 4

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
