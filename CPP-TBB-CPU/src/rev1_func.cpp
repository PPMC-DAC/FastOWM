// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <float.h>
// #include <math.h>
// #include <omp.h>
#include "../include/rev1_func.h"

static const int REALLOC_INCREMENT = 256;

Qtree_t::Qtree_t(Vector2D c, float r) : center(c), radius(r) {
  for(int i = 0; i < 4; i++)
    quadrants[i] = NULL;
}

unsigned int stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Qtree qtreeIn, Vector3D min, const unsigned chunk){

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    unsigned int countMin = 0;

#ifdef COLLAPSE
#ifdef STATIC
    #pragma omp parallel for collapse(2) schedule(static,chunk)
#elif defined(GUIDED)
    #pragma omp parallel for collapse(2) schedule(guided,chunk)
#else
    #pragma omp parallel for collapse(2) schedule(dynamic,chunk)
#endif
#else //COLLAPSE
#ifdef STATIC
    #pragma omp parallel for schedule(static,chunk)
#elif defined(GUIDED)
    #pragma omp parallel for schedule(guided,chunk)
#else
    #pragma omp parallel for schedule(dynamic,chunk)
#endif
#endif //COLLAPSE

    for( int jj = 0 ; jj < Ccol ; jj++ ){
        for( int ii=0 ; ii < Crow ; ii++ ){
            Lpoint cellCenter = {0, initX + ii*Displace, initY + jj*Displace, 0.0};
            int cellPoints = 0;
            Lpoint** neighbors = searchNeighbors2D(&cellCenter, qtreeIn, Wsize/2, &cellPoints);
            if(cellPoints >= minNumPoints ){
                int idmin = findMin(neighbors, cellPoints);
                #pragma omp critical
                {
                    minIDs[countMin] = neighbors[idmin]->id;
                    countMin++;
                }
            }

            free(neighbors);
            neighbors = NULL;
        }
    }

    return countMin;
}
unsigned int stage1t(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Qtree qtreeIn, Vector3D min, const unsigned chunk){

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    unsigned int countMin = 0;

    #pragma omp parallel
    {
        #pragma omp single
        {
            for( int jj = 0 ; jj < Ccol ; jj++ ){
                for( int ii=0 ; ii < Crow ; ii++ ){
                    #pragma omp task firstprivate(ii,jj)
                    {
                        Lpoint cellCenter = {0, initX + ii*Displace, initY + jj*Displace, 0.0};
                        int cellPoints = 0;
                        Lpoint** neighbors = searchNeighbors2D(&cellCenter, qtreeIn, Wsize/2, &cellPoints);
                        if(cellPoints >= minNumPoints ){
                            int idmin = findMin(neighbors, cellPoints);
                            #pragma omp critical
                            {
                                minIDs[countMin] = neighbors[idmin]->id;
                                countMin++;
                            }
                        }

                        free(neighbors);
                        neighbors = NULL;
                    }
                } // for ii
            } // for jj
        } // omp single
    } // omp parallel

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
          int* minGridIDs, Qtree qtreeIn, Qtree grid, Vector3D min){

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
                   neighbors = searchNeighbors2D(&cellCenter, qtreeIn, Bsize/2, &cellPoints);
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

double round2d(double z){
  return round(z*100.0)/100.0;
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

int isLeaf(Qtree quad)
{
    return quad->quadrants[0] == NULL;
}

int isEmpty(Qtree quad)
{
    return quad->points.size() == 0;
}


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


void createQuadrantsF(Qtree quad)
{
    int i = 0;
    Vector2D newCenter;
    float newRadius = quad->radius * 0.5;

    for(i = 0; i < 4; i++)
    {
        newCenter = quad->center;
        newCenter.x += quad->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += quad->radius * (i&1 ? 0.5f : -0.5f);
        // newCenter.z += quad->radius * (i&1 ? 0.5f : -0.5f);
        // printf("(%lf,%lf) ", newCenter.x, newCenter.y);

        // quad->quadrants[i] = new Qtree_t(quad, newCenter, newRadius);
        quad->quadrants[i] = new Qtree_t(newCenter, newRadius);
    }
    // printf("\n");
}

// Find the child corresponding a given point
int quadrantIdxF(Lpoint *point, Qtree qtree)
{
    int child = 0;

    if(point->x >= qtree->center.x) child |= 2;
    if(point->y >= qtree->center.y) child |= 1;
    // if(point->z >= qtree->center.z) child |= 1;

    return child;
}


void insertPointF(Lpoint *point, Qtree qtree, float minRadius)
{
    int idx = 0;

    if(isLeaf(qtree))
    {
        // printf("quadrante hoja nivel %d\n",nivel);
        if(qtree->radius * 0.5 > minRadius)    // still divisible -> divide
        {
          createQuadrantsF(qtree);
          // fillQuadrants(qtree);
          idx = quadrantIdxF(point, qtree);
          insertPointF(point, qtree->quadrants[idx], minRadius);

        } else {
          qtree->points.push_back(point);
        }
    }
    else                                // No leaf -> search the correct one
    {
      // printf("quadrante intermedio nivel %d\n", nivel);
      idx = quadrantIdxF(point, qtree);
      insertPointF(point, qtree->quadrants[idx], minRadius);
    }
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

int boxOverlap2D(Vector3D boxMin, Vector3D boxMax, Qtree quad)
{
    if(quad->center.x + quad->radius < boxMin.x ||
       quad->center.y + quad->radius < boxMin.y)
        return 0;

    if(quad->center.x - quad->radius > boxMax.x ||
       quad->center.y - quad->radius > boxMax.y)
        return 0;

    return 1;
}

Lpoint** neighbors2D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Qtree qtree, Lpoint **ptsInside, int *ptsInside_size, int *numInside)
{
    int i = 0;

    if(isLeaf(qtree))
    {
        if(!isEmpty(qtree))
        {
          size_t mysize = qtree->points.size();
            for(i = 0; i < mysize; i++)
            {
                if(insideBox2D(qtree->points[i], boxMin, boxMax))
                {
                    if (*numInside >= *ptsInside_size) {
                        (*ptsInside_size) += REALLOC_INCREMENT;
                        ptsInside = (Lpoint**)reallocWrap(ptsInside, *ptsInside_size * sizeof(Lpoint*));
                    }
                    ptsInside[(*numInside)++] = qtree->points[i];
                }
            }
        }
    }
    else
    {
        for(i = 0; i < 4; i++)
        {
            if(!boxOverlap2D(boxMin, boxMax, qtree->quadrants[i]))
            {
                continue;
            }
            else
            {
                ptsInside = neighbors2D(point, boxMin, boxMax, qtree->quadrants[i], ptsInside, ptsInside_size, numInside);
            }
        }
    }

    return ptsInside;
}

Lpoint** searchNeighbors2D(Lpoint *point, Qtree qtree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int ptsInside_size = 0;


    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);
    ptsInside = neighbors2D(point, boxMin, boxMax, qtree, ptsInside, &ptsInside_size, numNeighs);

    return ptsInside;
}


void deleteQtree(Qtree qtree)
{
    int i;
    if(isLeaf(qtree))
    {
      qtree->points.clear();
    } else {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(i = 0; i < 4; i++) {
            // Check
            deleteQtree(qtree->quadrants[i]);
            delete(qtree->quadrants[i]);
            // (qtree->quadrants[i]).reset(NULL); // free memory for unique_ptr
        }

    }

    return;
}
