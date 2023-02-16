#include "../include/envi.h"

int main(int argc, char* argv[])
{
  uOctree octreeIn = createOctreeF({10.0,10.0,10.0}, 20.0);

  size_t SIZE = (argc>1)? atol(argv[1]) : 1000;
  for(int rep=0; rep<4; rep++){
    printf("vuelta %d\n", rep);
    // uOctree octreeIn = createOctreeF({10.0,10.0,10.0}, 20.0);
    Lpoint* pointer = (Lpoint*)malloc(SIZE*sizeof(Lpoint));
    for(int i = 0; i < SIZE; i++){
      pointer[i].id = i;
      pointer[i].x = i%40;
      pointer[i].y = i%50;
      pointer[i].z = i%100;
      // insertPointMinRadius(&pointer[i], octreeIn, minRadius);
      insertPointF(&pointer[i], octreeIn, 1.0, 0);
      // printf("punto %d\n", i);
    }

    free(pointer);
    deleteOctree(octreeIn);
  }

  printf("FIN\n");

  return 0;
}
