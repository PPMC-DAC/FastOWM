#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include "../include/OWM_functions.h"

#define MIN_RADIUS 0.10 //Minimum size of the leaf bounding box

static const int REALLOC_INCREMENT = 256;

//This function returns the point with the minimum z-value of a Sliding Window (SW)
//The minium has to be selected among at least minNumPoints found inside the SW
//It returns the number of valid minimums found and the array minIDs that contains the IDs of the valid minimums

unsigned int stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min){

  double Displace = round2d(Wsize*(1-Overlap)); //Displacement between SW (sliding window)

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  unsigned int countMin = 0;

#ifdef PARALLEL
    #pragma omp parallel for schedule(dynamic,1)
#endif
      for(int jj = 0 ; jj < Ccol ; jj++ ){
        for( int ii=0 ; ii < Crow ; ii++ ){

            Lpoint cellCenter = {0, initX + ii*Displace, initY + jj*Displace, 0.0};

            // printf("Center of SW processed by thread %d: %.2f %.2f\n",omp_get_thread_num(), cellCenter.x, cellCenter.y);
            // printf("Search points inside the SW: neighbors\n");
            int cellPoints = 0; //Number of points inside the cell/SW
            Lpoint** neighbors = searchNeighbors2D(&cellCenter, octreeIn, Wsize/2, &cellPoints);
            if(cellPoints >= minNumPoints ){ //If there are enough points inside the cell/SW the minimum is considered
                //printf("Number of points found in the SW: %d\n", cellPoints );
                int idmin = findMin(neighbors, cellPoints);
#ifdef PARALLEL
                #pragma omp critical
#endif
                {
                    minIDs[countMin] = neighbors[idmin]->id;
                    countMin++;
                }
            //printf("Th:%d; Minimum at SW %04d,%04d has count:%d  and z-value:%.2f\n",omp_get_thread_num(),jj,ii,countMin ,neighbors[idmin]->z);  
            }

            free(neighbors);
            neighbors = NULL;
        }
     }

  return countMin; //Number of valid minimums found
}
//Sequential version of stage1
unsigned int stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min){

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
          
              for( ii = 0 ; ii < Crow ; ii++ ){

                  cellCenter.x = initX + ii*Displace;
                  neighbors = searchNeighbors2D(&cellCenter, octreeIn, Wsize/2, &cellPoints);
                  if(cellPoints >= minNumPoints ){
                      idmin = findMin(neighbors, cellPoints);
                      minIDs[countMin] = neighbors[idmin]->id;
                      countMin++;
                    //printf("Minimum at SW %04d,%04d has count:%d  and z-value:%.2f\n",jj,ii,countMin ,neighbors[idmin]->z);  
                  }
                  
                  free(neighbors);
                  neighbors = NULL;
              }
          }

  return countMin;
}

//This function returns in array minIDs the LLPs and index (the number of LLPs)
//An LLP (Local Lowest Point) is a valid minimum that is found more than once in stage1
unsigned int stage2(unsigned int countMin, int* minIDs){

  unsigned int index = 0;

  int ii,jj,id;

  for( ii=0 ; ii<countMin ; ii=jj ){
      id = minIDs[ii];
    // skipping repeated values
      for( jj=ii+1 ; id==minIDs[jj] && jj<countMin ; jj++ );

    //Check if there were repeated values
      if(jj-ii > 1){
          minIDs[index]=id;
          index++;
      }
  }

  return index;
}

//This function fills empty Bsize x Bsize cells that do not contain any LLP
unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min){

    unsigned int addMin = 0, idmin = 999;

    int cellPoints = 0;

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    Lpoint** neighbors = NULL;
    Lpoint** minimos = NULL;

    int ii,jj;
#ifdef PARALLEL
    #pragma omp parallel private(ii,jj,cellCenter,neighbors,minimos,cellPoints,idmin)
#endif
    {
#ifdef PARALLEL
        for( jj = omp_get_thread_num() ; jj < Ccol ; jj+=omp_get_num_threads() ){
#else        
        for( jj = 0 ; jj < Ccol ; jj++ ){
#endif
           cellCenter.y = min.y + Bsize/2 + jj*Bsize;

           for( ii = 0 ; ii < Crow ; ii++ ){
               cellCenter.x = min.x + Bsize/2 + ii*Bsize;

               minimos = searchNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);
               
               if(cellPoints == 0){ //If there are no LLPs in the cell
                   neighbors = searchNeighbors2D(&cellCenter, octreeIn, Bsize/2, &cellPoints);
                   if(cellPoints>0){
                       idmin = findMin(neighbors, cellPoints); //Find the mimimum in the cell (not an LLP, just the minimum z-value)
#ifdef PARALLEL
                       #pragma omp critical
#endif
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

    return addMin; //Number of minimums added

}
//Sequential version of stage 3
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

int isLeaf(Octree oct)
{
    return oct->octants[0] == NULL;
}

int isEmpty(Octree oct)
{
    return oct->numPts == 0;
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

/** Create an octree node with the given center and radius */
Octree createOctree(Vector3D center, float radius)
{
    int i = 0;
    Octree oct = NULL;

    oct = mallocWrap(sizeof(struct Octree));
    oct->center = center;
    oct->radius = radius;
    oct->points = NULL;
    oct->numPts = 0;
    // oct->plane[0] = 0.0;
    // oct->plane[1] = 0.0;
    // oct->plane[2] = 0.0;
    // oct->plane[3] = 0.0;
    for(i = 0; i < 8; i++)
        oct->octants[i] = NULL;

    return oct;
}

void createOctants(Octree oct)
{
    int i = 0;
    Vector3D newCenter;

    for(i = 0; i < 8; i++)
    {
        newCenter = oct->center;
        newCenter.x += oct->radius * (i&4 ? 0.5f : -0.5f);
        newCenter.y += oct->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.z += oct->radius * (i&1 ? 0.5f : -0.5f);
        oct->octants[i] = createOctree(newCenter, oct->radius * 0.5);
    }
}

// Find the child corresponding a given point
int octantIdx(Lpoint *point, Octree octree)
{
    int child = 0;

    if(point->x >= octree->center.x) child |= 4;
    if(point->y >= octree->center.y) child |= 2;
    if(point->z >= octree->center.z) child |= 1;

    return child;
}

/** Move the points inside the octree to the corresponding octants */
void fillOctants(Octree octree)
{
    int i = 0, idx = 0;

    for(i = 0; i < octree->numPts; i++)
    {
        idx = octantIdx(octree->points[i], octree);
        insertPoint(octree->points[i], octree->octants[idx]);
    }
    octree->numPts = 0;
    octree->points = NULL;
}

// Insert a point in the octree creating the appropiate childs
void insertPoint(Lpoint *point, Octree octree)
{
    int idx = 0;

    if(isLeaf(octree))
    {
        if(isEmpty(octree))             // Empty leaf -> insert point
        {
            octree->points = mallocWrap(sizeof(Lpoint*));
            octree->points[0] = point;
            octree->numPts = 1;
        }
        else                            // Not empty but still divisible -> divide
        {
            //if(octree->numPts > MAX_POINTS)
            // MIN_RADIUS defined as constat in environment.c
            if(octree->radius / 2.0 > MIN_RADIUS)
            {
                createOctants(octree);
                fillOctants(octree);
                idx = octantIdx(point, octree);
                insertPoint(point, octree->octants[idx]);

            }
            else                         // Not empty and isn't divisible -> insert point
            {
                octree->points = reallocWrap(octree->points, ++octree->numPts * sizeof(Lpoint*));
                octree->points[octree->numPts-1] = point;
            }
        }
    }
    else                                // No leaf -> search the correct one
    {
        idx = octantIdx(point, octree);
        insertPoint(point, octree->octants[idx]);
    }
}


// Insert a point in the octree creating the appropiate childs. Keep dividing until reaching minRadius
void insertPointMinRadius(Lpoint *point, Octree octree, float minRadius)
{
    int idx = 0;

    if(isLeaf(octree))
    {
        if(isEmpty(octree))             // Empty leaf -> insert point
        {
            octree->points = mallocWrap(sizeof(Lpoint*));
            octree->points[0] = point;
            octree->numPts = 1;
        }
        else                            // Not empty but still divisible -> divide
        {
            if(octree->radius / 2.0 > minRadius)
            {
                createOctants(octree);
                fillOctants(octree);
                idx = octantIdx(point, octree);
                insertPointMinRadius(point, octree->octants[idx],minRadius);

            }
            else                         // Not empty and isn't divisible -> insert point
            {
                octree->points = reallocWrap(octree->points, ++octree->numPts * sizeof(Lpoint*));
                octree->points[octree->numPts-1] = point;
            }
        }
    }
    else                                // No leaf -> search the correct one
    {
        idx = octantIdx(point, octree);
        insertPointMinRadius(point, octree->octants[idx],minRadius);
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

//This function returns the points inside a bounding box (boxMin, boxMax) in the octree
//The points are returned in array ptsInside and the number of points in the array is in numInside
//ptsInside_size is the size allocated for ptsInside (it can grow if needed)
Lpoint** neighbors2D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *ptsInside_size, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideBox2D(octree->points[i], boxMin, boxMax))
                {
                    if (*numInside >= *ptsInside_size) {
                        (*ptsInside_size) += REALLOC_INCREMENT;
                        ptsInside = reallocWrap(ptsInside, *ptsInside_size * sizeof(Lpoint*));
                    }
                    ptsInside[(*numInside)++] = octree->points[i];
                }
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            if(!boxOverlap2D(boxMin, boxMax, octree->octants[i]))
            {
                continue;
            }
            else
            {
                ptsInside = neighbors2D(point, boxMin, boxMax, octree->octants[i], ptsInside, ptsInside_size, numInside);
            }
        }
    }

    return ptsInside;
}

//Wrapper function for neighbors2D that receives the center of the SW (point) and the radius of the box
//With that, makeBox computes the corresponding BBox of the SW and calls neghbors2D
Lpoint** searchNeighbors2D(Lpoint *point, Octree octree, float radius, int *numInside)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int ptsInside_size = 0;

    *numInside = 0;
    makeBox(point, radius, &boxMin, &boxMax);
    ptsInside = neighbors2D(point, boxMin, boxMax, octree, ptsInside, &ptsInside_size, numInside);

    return ptsInside;
}
