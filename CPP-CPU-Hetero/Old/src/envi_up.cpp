// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <float.h>
// #include <math.h>
// #include <omp.h>
#include "../include/envi_up.h"

static const int REALLOC_INCREMENT = 256;

// Octree_t::Octree_t(Vector2D c, float r) : center(c), radius(r) {
//   for(int i = 0; i < 4; i++)
//     octants[i] = NULL;
// }

double round2d(double z){
  return round(z*100.0)/100.0;
}

// Make a box with center the point and the specified radius
void makeBox(Lpoint *point, float radius, Vector3D *min, Vector3D *max)
{
    // printf("Radio: %.2f\n", radius);
    // printf("Centro: [ %.2lf, %.2lf]\n",point->x,point->y );
    min->x = point->x - radius;
    min->y = point->y - radius;
    min->z = point->z - radius;
    max->x = point->x + radius;
    max->y = point->y + radius;
    max->z = point->z + radius;
}

int insideBox2D(Lpoint *point, Vector3D min, Vector3D max)
{
    if(point->x > min.x && point->y > min.y)
    {
        if(point->x < max.x && point->y < max.y)
        {
            return 1;
        }
    }

    return 0;
}

int boxOverlap2D(Vector3D boxMin, Vector3D boxMax, Octree oct)
{
    if(oct->center.x + oct->radius < boxMin.x ||
       oct->center.y + oct->radius < boxMin.y)
        return 0;

    if(oct->center.x - oct->radius > boxMax.x ||
       oct->center.y - oct->radius > boxMax.y)
        return 0;

    return 1;
}


unsigned int findMin(Lpoint** neighbors, unsigned int cellPoints) {
  unsigned int idmin=0;
  double zzmin = neighbors[0]->z;
  for(int i=1; i<cellPoints; i++)
    if(neighbors[i]->z < zzmin){
      zzmin = neighbors[i]->z;
      idmin=i;
    }
  return idmin;
}

/** Calculate the radius in each axis and save the max radius of the bounding box */
Vector3D getRadius(Vector3D min, Vector3D max, float *maxRadius)
{
    Vector3D radii;

    radii.x = (max.x - min.x) / 2.0;
    radii.y = (max.y - min.y) / 2.0;
    radii.z = (max.z - min.z) / 2.0;

    if(radii.x >= radii.y && radii.x >= radii.z)
    {
        *maxRadius = radii.x;
    }
    else if(radii.y >= radii.x && radii.y >= radii.z)
    {
        *maxRadius = radii.y;
    }
    else
    {
        *maxRadius = radii.z;
    }

    return radii;
}

/** Calculate the center of the bounding box */
Vector3D getCenter(Vector3D min, Vector3D radius)
{
    Vector3D center;

    center.x = min.x + radius.x;
    center.y = min.y + radius.y;
    center.z = min.z + radius.z;

    return center;
}

int isLeaf(Octree oct)
{
    return oct->octants[0] == NULL;
}

// int isEmpty(Octree oct)
// {
//     return oct->points.size() == 0;
// }


/** Wrapper to handle malloc() nicely */
void* mallocWrap(size_t size)
{
    void *ptr = malloc(size);
    if(!ptr)
    {
        fprintf(stderr, "Not enough memory\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void* reallocWrap(void *ptr, size_t size)
{
    void *tmp = realloc(ptr, size);
    if(tmp == NULL)
    {
        fprintf(stderr, "Error in realloc() of size %zu\n", size);
        exit(EXIT_FAILURE);
    }
    else
    {
        ptr = tmp;
    }
    return ptr;
}


void createOctantsF(Octree oct)
{
    int i = 0;
    Vector2D newCenter;
    float newRadius = oct->radius * 0.5;

    for(i = 0; i < 4; i++)
    {
        newCenter = oct->center;
        newCenter.x += oct->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += oct->radius * (i&1 ? 0.5f : -0.5f);
        // newCenter.z += oct->radius * (i&1 ? 0.5f : -0.5f);
        // printf("(%lf,%lf) ", newCenter.x, newCenter.y);

        // oct->octants[i] = new Octree_t(oct, newCenter, newRadius);
        // oct->octants[i] = new Octree_t(newCenter, newRadius);
        (oct->octants[i]).reset( new nLeaf(newCenter, oct->radius * 0.5) );
    }
    // printf("\n");
}

// Find the child corresponding a given point
int octantIdxF(Lpoint *point, Octree octree)
{
    int child = 0;

    if(point->x >= octree->center.x) child |= 2;
    if(point->y >= octree->center.y) child |= 1;
    // if(point->z >= octree->center.z) child |= 1;

    return child;
}


void insertPointF(Lpoint *point, uOctree& octree, float minRadius)
{
    int idx = 0;

    if(isLeaf(octree.get()))
    {
        // printf("octante hoja nivel %d\n",nivel);
        if(octree->radius * 0.5 > minRadius)    // still divisible -> divide
        {
          octree.reset( new nInter(octree->center, octree->radius) );
          createOctantsF(octree.get());
          // fillOctants(octree);
          idx = octantIdxF(point, octree.get());
          insertPointF(point, octree->octants[idx], minRadius);

        } else {
          // octree->points.push_back(point);
          // vector_t points = octree->getVector();
          octree->getVector().push_back(point);
        }
    }
    else                                // No leaf -> search the correct one
    {
      // printf("octante intermedio nivel %d\n", nivel);
      idx = octantIdxF(point, octree.get());
      insertPointF(point, octree->octants[idx], minRadius);
    }
}




Lpoint** neighbors2D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *ptsInside_size, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
      vector_t points = octree->getVector();
      uint64_t size = points.size();

      for(i = 0; i < size; i++)
      {
          if(insideBox2D(points[i], boxMin, boxMax))
          {
              if (*numInside >= *ptsInside_size) {
                  (*ptsInside_size) += REALLOC_INCREMENT;
                  ptsInside = (Lpoint**)reallocWrap(ptsInside, *ptsInside_size * sizeof(Lpoint*));
              }
              ptsInside[(*numInside)++] = points[i];
          }
      }

    }
    else
    {
        for(i = 0; i < 4; i++)
        {
            if(!boxOverlap2D(boxMin, boxMax, (octree->octants[i]).get() ))
            {
                continue;
            }
            else
            {
                ptsInside = neighbors2D(point, boxMin, boxMax, (octree->octants[i]).get(), ptsInside, ptsInside_size, numInside);
            }
        }
    }

    return ptsInside;
}

Lpoint** searchNeighbors2D(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int ptsInside_size = 0;


    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);
    ptsInside = neighbors2D(point, boxMin, boxMax, octree, ptsInside, &ptsInside_size, numNeighs);

    return ptsInside;
}


void deleteOctree(Octree octree)
{
    int i;
    if(isLeaf(octree))
    {
      // octree->points.clear();
      // vector_t points = octree->getVector();
      octree->getVector().clear();
    } else {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(i = 0; i < 4; i++) {
            // Check
            deleteOctree((octree->octants[i]).get());
            // delete(octree->octants[i]);
            (octree->octants[i]).reset(NULL); // free memory for unique_ptr
        }

    }

    return;
}


unsigned int stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min){

  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  Lpoint cellCenter = {0,0.0,0.0,0.0};

  // Lpoint** neighbors = NULL;

  unsigned int countMin = 0;

  int cellPoints = 0;

  int ii,jj;

  /*  #pragma omp parallel for shared(countMin) firstprivate(minIDs,octreeIn,Wsize,Overlap,Displace,initX,initY) \
                            private(cellCenter,neighbors,cellPoints,idmin) schedule(dynamic,block_size)
    */
  #pragma omp parallel
  {
      // printf("Thread: %d\n", omp_get_thread_num());
      // printf("Thread num: %d\n", omp_get_num_threads());
      for( int jj = omp_get_thread_num() ; jj < Ccol ; jj+=omp_get_num_threads() ){
      // for( jj = 1 ; jj < Ccol ; jj++ ){

          // cellCenter.y = initY + jj*Displace;
          // printf("\nCeldas no descartadas thread %d:   %d\n",omp_get_thread_num(), countMin);

          /*          #pragma omp parallel shared(countMin,minIDs) firstprivate(Crow,cellCenter,octreeIn,Wsize,Overlap,Displace,initX) \
                                    private(ii,neighbors,cellPoints,idmin)
          {
              for( ii = omp_get_thread_num() ; ii < Crow ; ii+=omp_get_num_threads() ){
            */
              #pragma omp parallel for private(cellPoints) schedule(dynamic,1)
              for( int ii=0 ; ii < Crow ; ii++ ){

                  // cellCenter.x = initX + ii*Displace;
                  Lpoint cellCenter = {0, initX + ii*Displace, initY + jj*Displace, 0.0};

                  // printf("Centro de %d: %.2f %.2f\n",omp_get_thread_num(), cellCenter.x, cellCenter.y);
                  // printf("Busco los vecinos\n");
                  Lpoint** neighbors = searchNeighbors2D(&cellCenter, octreeIn, Wsize/2, &cellPoints);
                  // printf("Numero de elementos de la celda: %d\n", cellPoints );
                  if(cellPoints >= minNumPoints ){
                      // printf("Numero de elementos de la celda: %d\n", cellPoints );
                      int idmin = findMin(neighbors, cellPoints);
                      #pragma omp critical
                      {
                          minIDs[countMin] = neighbors[idmin]->id;
                          countMin++;
                      }
                      // printf("El mínimo %d de la celda es: %.2f\n",countMin ,neighbors[idmin]->z);

                  }

                  free(neighbors);
                  neighbors = NULL;
              }
          }
      // }
  }

  return countMin;
}

unsigned int stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min){

  // Tamaño del bloque del scheduler de OMP
  // unsigned short block_size = 1;

  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  Lpoint cellCenter = {0,0.0,0.0,0.0};

  Lpoint** neighbors = NULL;

  unsigned int countMin = 0, idmin;

  int cellPoints = 0;

  int ii,jj;

      for( jj = 0 ; jj < Ccol ; jj++ ){

          cellCenter.y = initY + jj*Displace;
          // printf("\nCeldas no descartadas thread %d:   %d\n",omp_get_thread_num(), countMin);

              for( ii = 0 ; ii < Crow ; ii++ ){

                  cellCenter.x = initX + ii*Displace;
                  // printf("Centro de %d: %.2f %.2f\n",omp_get_thread_num(), cellCenter.x, cellCenter.y);
                  // printf("Busco los vecinos\n");
                  neighbors = searchNeighbors2D(&cellCenter, octreeIn, Wsize/2, &cellPoints);
                  // printf("Numero de elementos de la celda: %d\n", cellPoints );
                  if(cellPoints >= minNumPoints ){
                      // printf("Numero de elementos de la celda: %d\n", cellPoints );
                      idmin = findMin(neighbors, cellPoints);
                      // #pragma omp critical
                      // {
                          minIDs[countMin] = neighbors[idmin]->id;
                          countMin++;
                      // }
                      // printf("El mínimo %d de la celda es: %.2f\n",countMin ,neighbors[idmin]->z);

                  }

                  free(neighbors);
                  neighbors = NULL;
              }
          }
      // }
  // }

  return countMin;
}


unsigned int stage2(unsigned int countMin, int* minIDs){

  unsigned int index = 0;

  int ii,jj,id;

  for( ii=0 ; ii<countMin ; ii=jj ){
      id = minIDs[ii];

      for( jj=ii+1 ; id==minIDs[jj] && jj<countMin ; jj++ );

      /* esta es la forma de descartar las posiciones no utilizadas, inicializadas a -1 */
      if(jj-ii > 1){
          minIDs[index]=id;
          index++;
      }
  }

  return index;
}

unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min){

    unsigned int addMin = 0, idmin = 999;

    int cellPoints = 0;

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    Lpoint** neighbors = NULL;

    Lpoint** minimos = NULL;

    int ii,jj;

    #pragma omp parallel private(ii,jj,cellCenter,neighbors,minimos,cellPoints,idmin)
    {

        for( jj = omp_get_thread_num() ; jj < Ccol ; jj+=omp_get_num_threads() ){
        // for( jj = 0 ; jj < Ccol ; jj++ ){
           cellCenter.y = min.y + Bsize/2 + jj*Bsize;

           for( ii = 0 ; ii < Crow ; ii++ ){
               cellCenter.x = min.x + Bsize/2 + ii*Bsize;

               minimos = searchNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);
               // printf("Numero de elementos de la celda: %d\n", cellPoints );
               //Tengo que hacerlo porque el algoritmo no lo hace
               if(cellPoints == 0){
                   // printf("    Tengo una celda vacia en la malla\n" );
                   neighbors = searchNeighbors2D(&cellCenter, octreeIn, Bsize/2, &cellPoints);
                   if(cellPoints>0){

                       idmin = findMin(neighbors, cellPoints);
                       #pragma omp critical
                       {
                          minGridIDs[addMin] = neighbors[idmin]->id;
                          addMin++;
                       }
                       // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
                   }
                   free(neighbors);
                   neighbors = NULL;

               }
               free(minimos);
               minimos=NULL;
           }
        }
    }

    return addMin;

}

unsigned int stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min){

    unsigned int addMin = 0, idmin;

    int cellPoints = 0;

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    Lpoint** neighbors = NULL;

    Lpoint** minimos = NULL;

    int ii,jj;

        for( jj = 0 ; jj < Ccol ; jj++ ){
           cellCenter.y = min.y + Bsize/2 + jj*Bsize;

           for( ii = 0 ; ii < Crow ; ii++ ){
               cellCenter.x = min.x + Bsize/2 + ii*Bsize;

               minimos = searchNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);

               //Tengo que hacerlo porque el algoritmo no lo hace
               if(cellPoints == 0){

                   neighbors = searchNeighbors2D(&cellCenter, octreeIn, Bsize/2, &cellPoints);
                   if(cellPoints>0){

                      idmin = findMin(neighbors, cellPoints);
                      minGridIDs[addMin] = neighbors[idmin]->id;
                      addMin++;
                   }
                   free(neighbors);
                   neighbors = NULL;

               }
               free(minimos);
               minimos=NULL;
           }
        }
    // }

    return addMin;

}
