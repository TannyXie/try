/*
 *  Main header file for percolation code.
 */

/*
 *  System size LM LN
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

#define M LM/MPROC
#define N LN/NPROC

/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */

void mapwrite(char *percfile, int map[LM][LN], int ncluster);
void mapwritedynamic(char *percfile, int **map, int l, int ncluster);

/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);
