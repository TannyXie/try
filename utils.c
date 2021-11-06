#include <stdlib.h>
#include "percolate.h"

int rank_valid(int rank, int size)  {
  return (rank >=0 ) && (rank < size);
}