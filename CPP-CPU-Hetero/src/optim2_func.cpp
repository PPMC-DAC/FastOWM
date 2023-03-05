// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <float.h>
// #include <math.h>
// #include <omp.h>
#include "../include/optim2_func.hpp"

static const int REALLOC_INCREMENT = 256;

Qtree_t::Qtree_t(Vector2D c, float r) : center(c), radius(r) {
  for(int i = 0; i < 4; i++)
    quadrants[i] = NULL;
}


double round2d(double z){
  return round(z*100.0)/100.0;
}

/** Calculate the radius in each axis and save the max radius of the bounding box */
Vector2D getRadius(Vector2D &min, Vector2D &max, float *maxRadius)
{
    Vector2D radii;

    radii.x = (max.x - min.x) / 2.0;
    radii.y = (max.y - min.y) / 2.0;
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
Vector2D getCenter(Vector2D &min, Vector2D &radius)
{
    Vector2D center;

    center.x = min.x + radius.x;
    center.y = min.y + radius.y;
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

        quad->quadrants[i] = new Qtree_t(newCenter, newRadius);
    }
}

// Find the child corresponding a given point
int quadrantIdxF(Lpoint *point, Qtree qtree)
{
    int child = 0;

    if(point->x >= qtree->center.x) child |= 2;
    if(point->y >= qtree->center.y) child |= 1;

    return child;
}


void insertPointF(Lpoint *point, Qtree qtree, float minRadius)
{
    int idx = 0;

    if(isLeaf(qtree))
    {
        if(qtree->radius * 0.5 > minRadius)    // still divisible -> divide
        {
          createQuadrantsF(qtree);
          idx = quadrantIdxF(point, qtree);
          insertPointF(point, qtree->quadrants[idx], minRadius);

        } else {
          qtree->points.push_back(point);
        }
    }
    else                                // No leaf -> search the correct one
    {
      idx = quadrantIdxF(point, qtree);
      insertPointF(point, qtree->quadrants[idx], minRadius);
    }
}

// Make a box with center the point and the specified radius
void makeBox(Vector2D &point, float radius, Vector2D &min, Vector2D &max)
{
    min.x = point.x - radius;
    min.y = point.y - radius;
    max.x = point.x + radius;
    max.y = point.y + radius;
}

int insideBox2D(Lpoint* point, Vector2D &min, Vector2D &max)
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

int boxInside2D(Vector2D &boxMin, Vector2D &boxMax, Qtree quad)
{
    if(quad->center.x + quad->radius > boxMax.x ||
       quad->center.y + quad->radius > boxMax.y)
        return 0;

    if(quad->center.x - quad->radius < boxMin.x ||
       quad->center.y - quad->radius < boxMin.y)
        return 0;

    return 1;
}

int boxOverlap2D(Vector2D &boxMin, Vector2D &boxMax, Qtree quad)
{
    if(quad->center.x + quad->radius < boxMin.x ||
       quad->center.y + quad->radius < boxMin.y)
        return 0;

    if(quad->center.x - quad->radius > boxMax.x ||
       quad->center.y - quad->radius > boxMax.y)
        return 0;

    return 1;
}

int boxTotalOverlap2D(Vector2D &boxMin, Vector2D &boxMax, Qtree &quad)
{
    if(quad->center.x + quad->radius < boxMax.x ||
       quad->center.y + quad->radius < boxMax.y)
        return 0;

    if(quad->center.x - quad->radius > boxMin.x ||
       quad->center.y - quad->radius > boxMin.y)
        return 0;

    return 1;
}

void deleteQtree(Qtree qtree)
{
    int i;
    if(isLeaf(qtree))
    {
      qtree->points.clear();
    } else {
        for(i = 0; i < 4; i++) {
            // Check
            deleteQtree(qtree->quadrants[i]);
            delete(qtree->quadrants[i]);
        }

    }

    return;
}

void findValidMin(Qtree qtree, Vector2D &boxMin, Vector2D &boxMax, int &numInside, Lpoint * &minptr)
{
    if(isLeaf(qtree))
    {
      if(boxInside2D(boxMin, boxMax, qtree)){
        for(Lpoint* p : qtree->points) {
          if (p->z < minptr->z) {
              minptr = p;
          }
          numInside++; //passed by reference
        }
      } else {
        for(Lpoint* p : qtree->points) {
          if(insideBox2D(p, boxMin, boxMax))
          {
            if (p->z < minptr->z) {
                minptr = p;
            }
            numInside++;
          }
        }
      }

    } else {
        for(int i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(boxMin, boxMax, qtree->quadrants[i]))
                continue;
            else {
                findValidMin(qtree->quadrants[i], boxMin, boxMax, numInside, minptr);
            }
        }
    }
}

Lpoint searchNeighborsMin(Vector2D &SW_center, Qtree qtree, float radius, int &numInside)
{
    Vector2D boxMin, boxMax;
    Lpoint temp{0, 0.0, 0.0, std::numeric_limits<double>::max()};
    Lpoint *minptr = &temp; 
    numInside = 0;
    makeBox(SW_center, radius, boxMin, boxMax); //updates boxMin,boxMax to be de BBox of SW with center and radius

    findValidMin(qtree, boxMin, boxMax, numInside, minptr);
    return *minptr; 
}


void countNeighbors(Qtree qtree, Vector2D &boxMin, Vector2D &boxMax, int &numInside)
{
    int i;

    if(isLeaf(qtree))
    {
        for(Lpoint* p : qtree->points) {
          if(insideBox2D(p, boxMin, boxMax))
          {
            numInside++;
          }
        }
    } else {
        for(i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(boxMin, boxMax, qtree->quadrants[i]))
              continue;
            else {
              countNeighbors(qtree->quadrants[i], boxMin, boxMax, numInside);
            }
        }
    }
    return;
}

void countNeighbors2D(Vector2D &point, Qtree qtree, float radius, int &numInside)
{
    Vector2D boxMin, boxMax;

    numInside = 0;
    makeBox(point, radius, boxMin, boxMax);

    countNeighbors(qtree, boxMin, boxMax, numInside);

    return;
}


Lpoint** neighbors2D(Vector2D &point, Vector2D &boxMin, Vector2D &boxMax, Qtree qtree, Lpoint **ptsInside, int &ptsInside_size, int &numInside)
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
                    if (numInside >= ptsInside_size) {
                        ptsInside_size += REALLOC_INCREMENT;
                        ptsInside = (Lpoint**)reallocWrap(ptsInside, ptsInside_size * sizeof(Lpoint*));
                    }
                    ptsInside[numInside++] = qtree->points[i];
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

Lpoint** searchNeighbors2D(Vector2D &point, Qtree qtree, float radius, int &numInside)
{
    Vector2D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int ptsInside_size = 0;


    numInside = 0;
    makeBox(point, radius, boxMin, boxMax);
    ptsInside = neighbors2D(point, boxMin, boxMax, qtree, ptsInside, ptsInside_size, numInside);

    return ptsInside;
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

void stage1(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Qtree qtreeIn, Vector2D min){

  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  #pragma omp parallel for schedule(dynamic,1)
  for( int step = 0 ; step < Crow*Ccol ; step++ ){
          int ii=step/Ccol, jj=step%Ccol;
          Vector2D cellCenter={initX + ii*Displace, initY + jj*Displace};
          int cellPoints = 0;
//New method          
          Lpoint newmin = searchNeighborsMin(cellCenter, qtreeIn, Wsize/2, cellPoints);
          //printf("Step: %d.%d; Min id: %.2f; cellPoints: %d\n",ii,jj,newmin.id, cellPoints);
//Old method
#ifdef DEBUG
          Vector2D cellCenter_org = {initX + ii*Displace, initY + jj*Displace};
          int cellPoints_org=0;
          Lpoint** neighbors = searchNeighbors2D(cellCenter_org, qtreeIn, Wsize/2, cellPoints_org);
#endif
          if(cellPoints >= minNumPoints ){
#ifdef DEBUG
              int idmin = findMin(neighbors, cellPoints_org);
              if(neighbors[idmin]->id != newmin.id && neighbors[idmin]->z < newmin.z){
                printf("Step: %d.%d; Center of SW: %.2f %.2f\n",ii,jj,cellCenter_org.x, cellCenter_org.y);
                printf("ERROR (old,new) ids: (%d, %d), z: (%.3f %.3f); count:(%d, %d)\n", 
                  neighbors[idmin]->id, newmin.id,neighbors[idmin]->z, newmin.z, cellPoints_org, cellPoints);
              }
#endif
              minIDs[step] = newmin.id;
          }
#ifdef DEBUG
          free(neighbors);
          neighbors = NULL;
#endif
  }
}

void stage1tbb(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Qtree qtreeIn, Vector2D min){

  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  tbb::parallel_for(tbb::blocked_range<int>{0,Crow*Ccol},
                      [&](tbb::blocked_range<int> r ) {
        for( int step = r.begin() ; step < r.end() ; step++ ){
                int ii=step/Ccol, jj=step%Ccol;
                Vector2D cellCenter={initX + ii*Displace, initY + jj*Displace};
                int cellPoints = 0;
        //New method          
                Lpoint newmin = searchNeighborsMin(cellCenter, qtreeIn, Wsize/2, cellPoints);
                //printf("Step: %d.%d; Min id: %.2f; cellPoints: %d\n",ii,jj,newmin.id, cellPoints);
        //Old method
        #ifdef DEBUG
                Vector2D cellCenter_org = {initX + ii*Displace, initY + jj*Displace};
                int cellPoints_org=0;
                Lpoint** neighbors = searchNeighbors2D(cellCenter_org, qtreeIn, Wsize/2, cellPoints_org);
        #endif
                if(cellPoints >= minNumPoints ){
        #ifdef DEBUG
                    int idmin = findMin(neighbors, cellPoints_org);
                    if(neighbors[idmin]->id != newmin.id && neighbors[idmin]->z < newmin.z){
                        printf("Step: %d.%d; Center of SW: %.2f %.2f\n",ii,jj,cellCenter_org.x, cellCenter_org.y);
                        printf("ERROR (old,new) ids: (%d, %d), z: (%.3f %.3f); count:(%d, %d)\n", 
                        neighbors[idmin]->id, newmin.id,neighbors[idmin]->z, newmin.z, cellPoints_org, cellPoints);
                    }
        #endif
                    minIDs[step] = newmin.id;
                }
                else minIDs[step] = -1;
        #ifdef DEBUG
                free(neighbors);
                neighbors = NULL;
        #endif
        }
  });
}


//Receives a sorted list of minIDs with -1s at the beginning due to the SWs that do not have an LLP
unsigned int stage2(unsigned int Ncells, int* minIDs){

  int i=0;
  while(minIDs[i]<0 && i<Ncells) i++; //Skip -1s
  unsigned int counter = 0;
  int jj=0;
  for( int ii=i ; ii<Ncells ; ii=jj ){
      int id = minIDs[ii];
      for( jj=ii+1 ; id==minIDs[jj] && jj<Ncells ; jj++ );
      if(jj-ii > 1){
          minIDs[counter]=id;
          counter++;
      }
  }
  return counter;
}

unsigned int stage3(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Qtree qtreeIn, Qtree grid, Vector2D min){

    unsigned int addMin = 0;

    #pragma omp parallel for schedule(dynamic)
    for( int jj = 0 ; jj < Ccol ; jj++ ){
       Vector2D cellCenter;
       cellCenter.y = min.y + Bsize/2 + jj*Bsize;

       for( int ii = 0 ; ii < Crow ; ii++ ){
           cellCenter.x = min.x + Bsize/2 + ii*Bsize;
           int cellPoints = 0;
           countNeighbors2D(cellCenter, grid, Bsize/2, cellPoints);
           if(cellPoints == 0){
               Lpoint newmin = searchNeighborsMin(cellCenter, qtreeIn, Bsize/2, cellPoints);
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
