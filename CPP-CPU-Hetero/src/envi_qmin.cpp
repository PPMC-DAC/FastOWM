// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <float.h>
// #include <math.h>
// #include <omp.h>
#include "../include/envi_qmin.h"

static const int REALLOC_INCREMENT = 256;

Octree_t::Octree_t(Vector2D c, float r) : center(c), radius(r) {
  for(int i = 0; i < 4; i++)
    octants[i] = NULL;
}

// Octree_t::Octree_t(Octree oct, Vector2D c, float r) : oparent(oct), center(c), radius(r) {
//   for(int i = 0; i < 4; i++)
//     octants[i] = NULL;
// }

double round2d(double z){
  return round(z*100.0)/100.0;
}

/** Calculate the radius in each axis and save the max radius of the bounding box */
Vector2D getRadius(Vector2D min, Vector2D max, float *maxRadius)
{
    Vector2D radii;

    radii.x = (max.x - min.x) / 2.0;
    radii.y = (max.y - min.y) / 2.0;
    // radii.z = (max.z - min.z) / 2.0;

    if(radii.x >= radii.y)
    {
        *maxRadius = radii.x;
    }
    else
    {
        *maxRadius = radii.y;
    }

    return radii;
}

/** Calculate the center of the bounding box */
Vector2D getCenter(Vector2D min, Vector2D radius)
{
    Vector2D center;

    center.x = min.x + radius.x;
    center.y = min.y + radius.y;
    // center.z = min.z + radius.z;

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
// void* mallocWrap(size_t size)
// {
//     void *ptr = malloc(size);
//     if(!ptr)
//     {
//         fprintf(stderr, "Not enough memory\n");
//         exit(EXIT_FAILURE);
//     }
//     return ptr;
// }

// void* reallocWrap(void *ptr, size_t size)
// {
//     void *tmp = realloc(ptr, size);
//     if(tmp == NULL)
//     {
//         fprintf(stderr, "Error in realloc() of size %zu\n", size);
//         exit(EXIT_FAILURE);
//     }
//     else
//     {
//         ptr = tmp;
//     }
//     return ptr;
// }


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
        oct->octants[i] = new Octree_t(newCenter, newRadius);
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


void insertPointF(Lpoint *point, Octree octree, float minRadius)
{
    int idx = 0;

    if(isLeaf(octree))
    {
        // printf("octante hoja nivel %d\n",nivel);
        if(octree->radius * 0.5 > minRadius)    // still divisible -> divide
        {
          createOctantsF(octree);
          // fillOctants(octree);
          idx = octantIdxF(point, octree);
          insertPointF(point, octree->octants[idx], minRadius);

        } else {
          octree->points.push_back(point);
        }
    }
    else                                // No leaf -> search the correct one
    {
      // printf("octante intermedio nivel %d\n", nivel);
      idx = octantIdxF(point, octree);
      insertPointF(point, octree->octants[idx], minRadius);
    }
}

// Make a box with center the point and the specified radius
void makeBox(Vector2D *point, float radius, Vector2D *min, Vector2D *max)
{
    // printf("Radio: %.2f\n", radius);
    // printf("Centro: [ %.2lf, %.2lf]\n",point->x,point->y );
    min->x = point->x - radius;
    min->y = point->y - radius;
    // min->z = point->z - radius;
    max->x = point->x + radius;
    max->y = point->y + radius;
    // max->z = point->z + radius;
}

int insideBox2D(Lpoint *point, Vector2D min, Vector2D max)
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

int boxInside2D(Vector2D boxMin, Vector2D boxMax, Octree oct)
{
    if(oct->center.x + oct->radius > boxMax.x ||
       oct->center.y + oct->radius > boxMax.y)
        return 0;

    if(oct->center.x - oct->radius < boxMin.x ||
       oct->center.y - oct->radius < boxMin.y)
        return 0;

    return 1;
}

int boxOverlap2D(Vector2D boxMin, Vector2D boxMax, Octree oct)
{
    if(oct->center.x + oct->radius < boxMin.x ||
       oct->center.y + oct->radius < boxMin.y)
        return 0;

    if(oct->center.x - oct->radius > boxMax.x ||
       oct->center.y - oct->radius > boxMax.y)
        return 0;

    return 1;
}

int boxTotalOverlap2D(Vector2D boxMin, Vector2D boxMax, Octree oct)
{
    if(oct->center.x + oct->radius < boxMax.x ||
       oct->center.y + oct->radius < boxMax.y)
        return 0;

    if(oct->center.x - oct->radius > boxMin.x ||
       oct->center.y - oct->radius > boxMin.y)
        return 0;

    return 1;
}

void deleteOctree(Octree octree)
{
    int i;
    if(isLeaf(octree))
    {
      octree->points.clear();
    } else {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(i = 0; i < 4; i++) {
            // Check
            deleteOctree(octree->octants[i]);
            delete(octree->octants[i]);
            // (octree->octants[i]).reset(NULL); // free memory for unique_ptr
        }

    }

    return;
}

Lpoint findValidMin(Octree octree, Vector2D* boxMin, Vector2D* boxMax, int* numInside)
{
    // Lpoint tmp, min = nomin;
    Lpoint tmp, min = {0,0.0,0.0,99999.0};
    int i;

    if(isLeaf(octree))
    {

      if(boxInside2D(*boxMin, *boxMax, octree)){
        for(Lpoint* p : octree->points) {
          if (p->z < min.z) {
              min = *p;
          }
          (*numInside)++;
        }
      } else {
        for(Lpoint* p : octree->points) {
          if(insideBox2D(p, *boxMin, *boxMax))
          {
            if (p->z < min.z) {
                min = *p;
            }
            (*numInside)++;
          }
        }
      }

    } else {

        for(i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(*boxMin, *boxMax, octree->octants[i]))
                continue;
            else {
                tmp = findValidMin(octree->octants[i], boxMin, boxMax, numInside);
                if (tmp.z < min.z) {
                    min = tmp;
                }
            }

        }

    }

    return min;
}

Lpoint searchNeighborsMin(Vector2D* point, Octree octree, float radius, int* numNeighs)
{
    Vector2D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    return findValidMin(octree, &boxMin, &boxMax, numNeighs);
}


void countNeighbors(Octree octree, Vector2D* boxMin, Vector2D* boxMax, int* numInside)
{
    int i;

    if(isLeaf(octree))
    {
        // uint64_t size = octree->points.size();
        // if(size > 0)
        // {
        //     for (i = 0; i < size; i++)
        //     {
        //       if(insideBox2D(octree->points[i], *boxMin, *boxMax))
        //       {
        //         (*numInside)++;
        //       }
        //     }
        //
        // }

        for(Lpoint* p : octree->points) {
          if(insideBox2D(p, *boxMin, *boxMax))
          {
            (*numInside)++;
          }
        }

    } else {

        for(i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(*boxMin, *boxMax, octree->octants[i]))
              continue;
            else {
              countNeighbors(octree->octants[i], boxMin, boxMax, numInside);
            }

        }

    }

    return;
}

void countNeighbors2D(Vector2D* point, Octree octree, float radius, int* numNeighs)
{
    Vector2D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    countNeighbors(octree, &boxMin, &boxMax, numNeighs);

    return;
}

Lpoint findValidMinPath(Octree octree, Vector2D* boxMin, Vector2D* boxMax,
                    int* numInside, std::vector<Octree>& path)
{
  Lpoint tmp, min = {0,0.0,0.0,99999.0};
  int i;

  if(isLeaf(octree))
  {

    for(Lpoint* p : octree->points) {
      if(insideBox2D(p, *boxMin, *boxMax))
      {
        if (p->z < min.z) {
            min = *p;
        }
        (*numInside)++;
      }
    }

  } else {
    std::vector<Octree> local_path, best_path;

    for(i = 0; i < 4; i++) {
      // Check
      if(!boxOverlap2D(*boxMin, *boxMax, octree->octants[i]))
        continue;
      else {
        local_path.push_back(octree->octants[i]);
        tmp = findValidMinPath(octree->octants[i], boxMin, boxMax, numInside, local_path);
        if (tmp.z < min.z) {
          min = tmp;
          best_path = local_path;
        }
        local_path.clear();
      }
    }

    path.insert(path.end(), best_path.begin(), best_path.end());

  }

    return min;
}

Lpoint searchNeighborsMinPath(Vector2D* point, Octree octree, float radius, int* numNeighs, std::vector<Octree>& path)
{
    Vector2D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    return findValidMinPath(octree, &boxMin, &boxMax, numNeighs, path);
}


// Lpoint findValidMinParent(Octree octree, Vector2D* boxMin, Vector2D* boxMax,
//                     int* numInside, Octree& minparent)
// {
//   Lpoint tmp, min = {0,0.0,0.0,99999.0};
//   int i;
//
//   if(isLeaf(octree))
//   {
//     if(boxInside2D(*boxMin, *boxMax, octree)){
//       for(Lpoint* p : octree->points) {
//         if (p->z < min.z) {
//             min = *p;
//         }
//         (*numInside)++;
//       }
//     } else {
//       for(Lpoint* p : octree->points) {
//         if(insideBox2D(p, *boxMin, *boxMax))
//         {
//           if (p->z < min.z) {
//               min = *p;
//           }
//           (*numInside)++;
//         }
//       }
//     }
//
//   } else {
//     for(i = 0; i < 4; i++) {
//       /*Tengo que seguir haciendo el chequeo*/
//       if(!boxOverlap2D(*boxMin, *boxMax, octree->octants[i]))
//         continue;
//       else {
//         tmp = findValidMinParent(octree->octants[i], boxMin, boxMax, numInside, minparent);
//         if (tmp.z < min.z) {
//           min = tmp;
//           if(isLeaf(octree->octants[i]))
//             minparent = octree->oparent;
//         }
//       }
//     }
//   }
//
//     return min;
// }
//
//
// Lpoint searchNeighborsMinParent(Vector2D* point, Octree& oparent, float radius, int* numNeighs)
// {
//     Vector2D boxMin, boxMax;
//     Octree actual_parent = oparent;
//
//     *numNeighs = 0;
//     makeBox(point, radius, &boxMin, &boxMax);
//
//     while(actual_parent->oparent != NULL && !boxTotalOverlap2D(boxMin, boxMax, actual_parent)){
//     // while(actual_parent->oparent != NULL){
//       actual_parent = actual_parent->oparent;
//     }
//
//     // if(actual_parent->oparent != NULL) {
//     //   // printf("my size: %zu\n", (*oparent)->points.size());
//     //   actual_parent = actual_parent->oparent;
//     //   // printf("parent size: %zu\n", (*oparent)->points.size());
//     //   // if(actual_parent->oparent != NULL)
//     //   //   actual_parent = actual_parent->oparent;
//     // }
//
//     return findValidMinParent(actual_parent, &boxMin, &boxMax, numNeighs, oparent);
// }

// unsigned int stage1cppParent(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
//   unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector2D min)
// {
//
//     // double initX = min.x + Wsize/2 - Wsize*Overlap;
//     // double initY = min.y + Wsize/2 - Wsize*Overlap;
//     double initX = min.x - Wsize/2 + Wsize*Overlap;
//     double initY = min.y - Wsize/2 + Wsize*Overlap;
//
//     double Displace = round2d(Wsize*(1-Overlap));
//
//     Vector2D cellCenter;
//
//     // Lpoint** neighbors = NULL;
//     Lpoint newmin;
//
//     unsigned int countMin = 0;
//
//     int cellPoints;
//
//     // std::vector<int> vcount(Ccol,0);
//     Octree mypointer;
//     // std::vector<Octree> path;
//
//     #pragma omp parallel for private(cellCenter,cellPoints,newmin,mypointer) schedule(dynamic)
//     for(int jj = 0 ; jj < Ccol ; jj++ ){
//         /* lo renuevo para cada fila porque cada vez que cambio de fila voy
//         a estar lejos de donde estaba antes en el mapa*/
//         mypointer = octreeIn;
//         // path.clear();
//
//         cellCenter.y = initY + jj*Displace;
//         // printf("\nCeldas no descartadas thread %d:   %d\n",omp_get_thread_num(), countMin);
//
//         for(int ii = 0 ; ii < Crow ; ii++ ){
//
//             cellCenter.x = initX + ii*Displace;
//
//             newmin = searchNeighborsMinParent(&cellCenter, mypointer, Wsize/2, &cellPoints);
//             // newmin = searchNeighborsMinPath(&cellCenter, octreeIn, Wsize/2, &cellPoints, path);
//             // printf("Numero de elementos de la celda: %d\n", cellPoints );
//             if(cellPoints >= minNumPoints ){
//               #pragma omp critical
//               {
//                   minIDs[countMin] = newmin.id;
//                   countMin++;
//               }
//             }
//         }
//         // vcount[jj] = countMin;
//         // countMin=0;
//     }
//
//     // countMin = std::accumulate(vcount.begin(), vcount.end(), 0);
//
//     return countMin;
// }



unsigned int stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector2D min){

  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  Vector2D cellCenter;

  Lpoint newmin;

  unsigned int countMin = 0;

  int cellPoints = 0;

  #pragma omp parallel for private(cellCenter,cellPoints,newmin) schedule(dynamic)
  for( int jj = 0 ; jj < Ccol ; jj++ ){

    cellCenter.y =  initY + jj*Displace;

      for( int ii = 0 ; ii < Crow ; ii++ ){

          cellCenter.x = initX + ii*Displace;

          newmin = searchNeighborsMin(&cellCenter, octreeIn, Wsize/2, &cellPoints);

          if(cellPoints >= minNumPoints ){
              #pragma omp critical
              {
                  minIDs[countMin] = newmin.id;
                  countMin++;
              }
          }
      }
  }


  return countMin;
}

unsigned int stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector2D min){

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    unsigned int countMin = 0;

    int cellPoints;

        for( int jj = 0 ; jj < Ccol ; jj++ ){

                for( int ii=0 ; ii < Crow ; ii++ ){

                    Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

                    Lpoint newmin = searchNeighborsMin(&cellCenter, octreeIn, Wsize/2, &cellPoints);

                    if(cellPoints >= minNumPoints ){
                        #pragma omp critical
                        {
                            minIDs[countMin] = newmin.id;
                            countMin++;
                        }
                    }

                }
            }

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
          int* minGridIDs, Octree octreeIn, Octree grid, Vector2D min){

    unsigned int addMin = 0;

    Vector2D cellCenter;

    int cellPoints;

    #pragma omp parallel for private(cellCenter,cellPoints) schedule(dynamic)
    for( int jj = 0 ; jj < Ccol ; jj++ ){
       cellCenter.y = min.y + Bsize/2 + jj*Bsize;

       for( int ii = 0 ; ii < Crow ; ii++ ){
           cellCenter.x = min.x + Bsize/2 + ii*Bsize;

           countNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);

           //Tengo que hacerlo porque el algoritmo no lo hace
           if(cellPoints == 0){
               // printf("    Tengo una celda vacia en la malla\n" );
               Lpoint newmin = searchNeighborsMin(&cellCenter, octreeIn, Bsize/2, &cellPoints);
               if(cellPoints>0){
                   #pragma omp critical
                   {
                      minGridIDs[addMin] = newmin.id;
                      addMin++;
                   }
                   // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
               }
           }
       }
    }


    return addMin;

}

unsigned int stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector2D min){

    unsigned int addMin = 0;

    Vector2D cellCenter;

    int cellPoints;

    for( int jj = 0 ; jj < Ccol ; jj++ ){
       cellCenter.y = min.y + Bsize/2 + jj*Bsize;

       for( int ii = 0 ; ii < Crow ; ii++ ){
           cellCenter.x = min.x + Bsize/2 + ii*Bsize;

           countNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);

           //Tengo que hacerlo porque el algoritmo no lo hace
           if(cellPoints == 0){
               // printf("    Tengo una celda vacia en la malla\n" );
               Lpoint newmin = searchNeighborsMin(&cellCenter, octreeIn, Bsize/2, &cellPoints);
               if(cellPoints>0){
                   #pragma omp critical
                   {
                      minGridIDs[addMin] = newmin.id;
                      addMin++;
                   }
                   // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
               }
           }
       }
    }

    return addMin;

}
