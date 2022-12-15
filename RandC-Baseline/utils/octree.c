#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "octree.h"
#include "file.h"
#include "util.h"
#include "plane.h"
#include "math.h"
//#include "laspec.h"
#include "wrapper.h"

//#define MIN_RADIUS 1.5  //!< Stop dividing the octree when childs have this radius
#define MIN_RADIUS 0.10 //For the discretization
#define MAX_POINTS 100  //!< Keep dividing the octree while octants have more points than this
#define SORT 1

static const int REALLOC_INCREMENT = 1024;

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


/** Create a octree with the given center and radius */
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

/** Simply insert points in the octree */
void insertPoints(Lpoint *points, int numPts, Octree oct)
{
    int i = 0;

    for(i = 0; i < numPts; i++)
        insertPoint(&points[i], oct);
    printf("%d points inserted\n", numPts);
}

int getPointGroup(Lpoint point, Octree octree)
/*Función que devuelve el grupo al que pertenece un punto arbitrario.
Se calcula asociando al punto point el grupo del primer vecino encontrado en un radio determinado.*/
{
  int i, numNeighs, gid;
  float currDistance = 0, distance;
  float radius = 1;
  Lpoint **neighs = NULL;
  neighs = searchNeighbors2D(&point, octree, radius, &numNeighs);

  if (numNeighs == 0) return(-1);

  distance = sqrt((point.x-neighs[0]->x) * (point.x-neighs[0]->x) + (point.y-neighs[0]->y) * (point.y-neighs[0]->y));
  gid = neighs[0]->gId;

  for (i = 1; i < numNeighs; i++)
  {
    currDistance = sqrt((point.x-neighs[i]->x) * (point.x-neighs[i]->x) + (point.y-neighs[i]->y) * (point.y-neighs[i]->y));
    if (currDistance < distance)
    {
      distance = currDistance;
      gid = neighs[i]->gId;
    }
  }
  return(gid);
}

///** Create an octree storing the point's pointers */
//Octree buildOctree(Lpoint *points, int numPts, Vector3D *octCenter, float *octRadius, int label)
//{
//    FILE *fp = NULL;
//    Octree oct = NULL;
//
//    printf("Building octree...\n");
//    *octCenter = mbb(points, numPts, octRadius); // @todo: This consider all labels
//    oct = createOctree(*octCenter, *octRadius, NULL);
//    if(label == -1)
//        insertPoints(points, numPts, oct, fp);
//    else
//        insertPointsAvoidLabel(points, numPts, oct, label, fp);
//#if DEBUG
//    fp = fopen(OCTREE_FILE, "w");
//    //printOctreeCenters(oct, fp, 0);
//    writeLeafQuads(oct, fp, 1);
//#endif
//
//    return oct;
//}

// Find the child corresponding a given point
int octantIdx(Lpoint *point, Octree octree)
{
    int child = 0;

    if(point->x >= octree->center.x) child |= 4;
    if(point->y >= octree->center.y) child |= 2;
    if(point->z >= octree->center.z) child |= 1;

    return child;
}

int isLeaf(Octree oct)
{
    return oct->octants[0] == NULL;
}

int isEmpty(Octree oct)
{
    return oct->numPts == 0;
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

void mergeOctants(Octree dst, Octree src)
{
    int newSize = 0;

    newSize = dst->numPts + src->numPts;
    dst->points = reallocWrap(dst->points, newSize * sizeof(Lpoint*));
    memcpy(&dst->points[dst->numPts], src->points, src->numPts * sizeof(Lpoint*));
    dst->numPts = newSize;
    removeOctree(src);
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
    int i = 0, idx = 0;

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

// Insert a point in the octree creating the appropiate childs. Keep dividing until reaching radius
void insertPointMinRadius(Lpoint *point, Octree octree, float minRadius)
{
    int i = 0, idx = 0;

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

// Insert a point in the octree creating the appropiate childs. Keep dividing until reaching max points in leaf
void insertPointMaxPoints(Lpoint *point, Octree octree, int maxPoints)
{
    int i = 0, idx = 0;

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
            if(octree->numPts > MAX_POINTS)
            {
                createOctants(octree);
                fillOctants(octree);
                idx = octantIdx(point, octree);
                insertPointMaxPoints(point, octree->octants[idx],maxPoints);

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
        insertPointMaxPoints(point, octree->octants[idx],maxPoints);
    }
}


void octreeDiscretization(Octree octree, Octree octreeDisc, int *totalPoints, int *newPoints, float minRadius)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            // First point
            //writePointRaw(file, octree->points[0]);
            // Random point
            //int index = rand() % octree->numPts;
            //writePointRaw(file, octree->points[index]);
            // Lowest point
//            int j = 0, minIdx = 0;
//            double minZ = octree->points[0]->z;
//            for(j = 1; j < octree->numPts; j++)
//                if(octree->points[j]->z < minZ)
//                {
//                    minZ = octree->points[j]->z;
//                    minIdx = j;
//                }
//            //Highest point
//            int j = 0, maxIdx = 0;
//            double maxZ = octree->points[0]->z;
//            for(j = 1; j < octree->numPts; j++)
//                if(octree->points[j]->z > maxZ)
//                {
//                    maxZ = octree->points[j]->z;
//                    maxIdx = j;
//                }
            //Lowest and Highest point
            int j = 0, maxIdx = 0, minIdx = 0;
            double maxZ = octree->points[0]->z;
            double minZ = octree->points[0]->z;
            *totalPoints += octree->numPts;
            for(j = 1; j < octree->numPts; j++){
                if(octree->points[j]->z > maxZ)
                {
                    maxZ = octree->points[j]->z;
                    maxIdx = j;
                }else if(octree->points[j]->z < minZ)
                {
                    minZ = octree->points[j]->z;
                    minIdx = j;
                }
            }
            if(minIdx!=maxIdx){
                insertPointMinRadius(octree->points[maxIdx], octreeDisc,minRadius);
                insertPointMinRadius(octree->points[minIdx], octreeDisc,minRadius);
                *newPoints += 2;
            }else{
                insertPointMinRadius(octree->points[maxIdx], octreeDisc,minRadius);
                *newPoints += 1;
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            octreeDiscretization(octree->octants[i], octreeDisc, totalPoints, newPoints,minRadius);
        }
    }
}


void insertPointNoSplit(Lpoint *point, Octree octree)
{
    if(isEmpty(octree))
    {
        octree->points = mallocWrap(sizeof(Lpoint*));
        octree->points[0] = point;
        octree->numPts = 1;
    }
    else
    {
        octree->points = reallocWrap(octree->points, ++octree->numPts * sizeof(Lpoint*));
        octree->points[octree->numPts-1] = point;
    }
}

void populateChilds(Octree oct)
{
    int i = 0, newOct = 0;

    for(i = 0; i < oct->numPts; i++)
    {
        newOct = octantIdx(oct->points[i], oct);
        insertPointNoSplit(oct->points[i], oct->octants[newOct]);
    }

    oct->numPts = 0;
    free(oct->points);
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

int insideBox3D(Lpoint *point, Vector3D min, Vector3D max)
{
    if(point->x > min.x && point->y > min.y  && point->z > min.z)
    {
        if(point->x < max.x && point->y < max.y && point->z < max.z)
        {
            return 1;
        }
    }

    return 0;
}

int octantInsideBox(Octree oct, Vector3D min, Vector3D max)
{
    if(oct->center.x - oct->radius > min.x && oct->center.y - oct->radius > min.y  && oct->center.z - oct->radius > min.z)
    {
        if(oct->center.x + oct->radius < max.x && oct->center.y + oct->radius < max.y && oct->center.z + oct->radius < max.z)
        {
            return 1;
        }
    }

    return 0;
}

Vector3D subVector(Vector3D v1, float number)
{
    Vector3D result;

    result.x = v1.x - number;
    result.y = v1.y - number;
    result.z = v1.z - number;

    return result;
}

Vector3D addVector(Vector3D v1, float number)
{
    Vector3D result;

    result.x = v1.x + number;
    result.y = v1.y + number;
    result.z = v1.z + number;

    return result;
}

int boxOverlap3DOld(Vector3D boxMin, Vector3D boxMax, Octree oct)
{
    Vector3D octMin = subVector(oct->center, oct->radius);
    Vector3D octMax = addVector(oct->center, oct->radius);

    if(octMax.x < boxMin.x || octMax.y < boxMin.y || octMax.z < boxMin.z)  return 0;
    if(octMin.x > boxMax.x || octMin.y > boxMax.y || octMin.z > boxMax.z)  return 0;

    return 1;
}

int boxOverlap3D(Vector3D boxMin, Vector3D boxMax, Octree oct)  // TODO CPU: Precomput octant min and max boxes (need more memory)
{
    if(oct->center.x + oct->radius < boxMin.x ||
       oct->center.y + oct->radius < boxMin.y ||
       oct->center.z + oct->radius < boxMin.z)
        return 0;

    if(oct->center.x - oct->radius > boxMax.x ||
       oct->center.y - oct->radius > boxMax.y ||
       oct->center.z - oct->radius > boxMax.z)
        return 0;

    return 1;
}

int boxOverlap2DOld(Vector3D boxMin, Vector3D boxMax, Octree octree)
{
    Vector3D octMin = subVector(octree->center, octree->radius);
    Vector3D octMax = addVector(octree->center, octree->radius);

    if(octMax.x < boxMin.x || octMax.y < boxMin.y)  return 0;
    if(octMin.x > boxMax.x || octMin.y > boxMax.y)  return 0;

    return 1;
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

Lpoint** neighbors2D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *ptsInside_size, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideBox2D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id)
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

// Distinc group id
Lpoint** neighbors2DDG(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideBox2D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id && point->gId != octree->points[i]->gId)
                {
                    ptsInside = reallocWrap(ptsInside, ++*numInside * sizeof(Lpoint*));
                    ptsInside[*numInside-1] = octree->points[i];

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
                ptsInside = neighbors2DDG(point, boxMin, boxMax, octree->octants[i], ptsInside, numInside);
            }
        }
    }

    return ptsInside;
}

Lpoint** neighbors2DReset(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideBox2D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id)
                {
                    ptsInside = reallocWrap(ptsInside, ++*numInside * sizeof(Lpoint*));
                    octree->points[i]->id  = *numInside - 1;
                    ptsInside[*numInside-1] = octree->points[i];
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
                ptsInside = neighbors2DReset(point, boxMin, boxMax, octree->octants[i], ptsInside, numInside);
            }
        }
    }

    return ptsInside;
}


Lpoint** neighbors3D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *ptsInside_size, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            if(octantInsideBox(octree, boxMin, boxMax))
            {
                //ptsInside = reallocWrap(ptsInside, (*numInside + octree->numPts) * sizeof(Lpoint*));
                for(i = 0; i < octree->numPts; i++)
                {
                    if(point->id != octree->points[i]->id)
                    {
                       if (*numInside >= *ptsInside_size) {
                            (*ptsInside_size) += REALLOC_INCREMENT;
                            ptsInside = reallocWrap(ptsInside, *ptsInside_size * sizeof(Lpoint*));
                        }
                       ptsInside[(*numInside)++] = octree->points[i];
                    }
                }
            }
            else
            {
                for(i = 0; i < octree->numPts; i++)
                {
                    if(insideBox3D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id)
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
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            if(!boxOverlap3D(boxMin, boxMax, octree->octants[i]))
            {
                continue;
            }
            else
            {
                ptsInside = neighbors3D(point, boxMin, boxMax, octree->octants[i], ptsInside, ptsInside_size, numInside);
            }
        }
    }

    return ptsInside;
}

Lpoint** neighbors2DNoGroupLowerZ(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int gId, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideBox2D(octree->points[i], boxMin, boxMax) && octree->points[i]->z < point->z && point->id != octree->points[i]->id && octree->points[i]->gId != gId)
                {
                    ptsInside = reallocWrap(ptsInside, ++*numInside * sizeof(Lpoint*));
                    ptsInside[*numInside-1] = octree->points[i];
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
                ptsInside = neighbors2DNoGroupLowerZ(point, boxMin, boxMax, octree->octants[i], ptsInside, gId, numInside);
            }
        }
    }

    return ptsInside;
}

Lpoint** neighbors3DNoGroup(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int gId, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideBox3D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id &&  octree->points[i]->gId != gId)
                {
                    ptsInside = reallocWrap(ptsInside, ++*numInside * sizeof(Lpoint*));
                    ptsInside[*numInside-1] = octree->points[i];
                }
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            if(!boxOverlap3D(boxMin, boxMax, octree->octants[i]))
            {
                continue;
            }
            else
            {
                ptsInside = neighbors3DNoGroup(point, boxMin, boxMax, octree->octants[i], ptsInside, gId, numInside);
            }
        }
    }

    return ptsInside;
}


Lpoint** neighbors3DOld(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *ptsInside_size,  int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideBox3D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id)
                {
                    //ptsInside = reallocWrap(ptsInside, ++*numInside * sizeof(Lpoint*));
                    //ptsInside[*numInside-1] = octree->points[i];
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
            if(!boxOverlap3D(boxMin, boxMax, octree->octants[i]))
            {
                continue;
            }
            else
            {
                ptsInside = neighbors3DOld(point, boxMin, boxMax, octree->octants[i], ptsInside, ptsInside_size, numInside);
            }
        }
    }

    return ptsInside;
}

Lpoint** circleNeighbors(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *ptsInside_size, int *numInside, float circleRadius)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideCircle(octree->points[i], point, circleRadius) && point->id != octree->points[i]->id)
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
            if(!boxOverlap3D(boxMin, boxMax, octree->octants[i]))
            {
                continue;
            }
            else
            {
                ptsInside = circleNeighbors(point, boxMin, boxMax, octree->octants[i], ptsInside, ptsInside_size, numInside, circleRadius);
            }
        }
    }

    return ptsInside;
}

int numNeighbors3D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideBox3D(octree->points[i], boxMin, boxMax) && octree->points[i]->intensity != -1.0)
                {
                    ++*numInside;
                }
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            if(!boxOverlap3D(boxMin, boxMax, octree->octants[i]))
            {
                continue;
            }
            else
            {
                *numInside = numNeighbors3D(point, boxMin, boxMax, octree->octants[i], numInside);
            }
        }
    }

    return *numInside;
}

int numNeighbors2DDG(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, int *numInside) // WARNING: Avoid same point id?
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(octree->points[i]->gId != point->gId && insideBox2D(octree->points[i], boxMin, boxMax))
                {
                    ++*numInside;
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
                *numInside = numNeighbors2DDG(point, boxMin, boxMax, octree->octants[i], numInside);
            }
        }
    }

    return *numInside;
}

int numNeighbors2D(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, int *numInside)
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
                    ++*numInside;
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
                *numInside = numNeighbors2D(point, boxMin, boxMax, octree->octants[i], numInside);
            }
        }
    }

    return *numInside;
}

int numNeighbors2DDiffGroup(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideBox2D(octree->points[i], boxMin, boxMax)
                        && octree->points[i]->gId != point->gId)
                {
                    ++*numInside;
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
                *numInside = numNeighbors2DDiffGroup(point, boxMin, boxMax, octree->octants[i], numInside);
            }
        }
    }

    return *numInside;
}

// Return also the octrees containing inside points
Lpoint** neighborsNBH(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint** ptsInside, Octree **octrees, int *numInside)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {

            if(octantInsideBox(octree, boxMin, boxMax))
            {
                ptsInside = reallocWrap(ptsInside, (*numInside + octree->numPts) * sizeof(Lpoint*));
                *octrees = reallocWrap(*octrees, (*numInside + octree->numPts) * sizeof(Octree));
                memcpy(&ptsInside[*numInside], octree->points, octree->numPts * sizeof(Lpoint*));
                for(i = 0; i < octree->numPts; i++)
                      (*octrees)[i+*numInside] = octree;
                *numInside += octree->numPts;
            }

            else
            {
                for(i = 0; i < octree->numPts; i++)
                {   //printf("[%d] %.2f\n", i, octree->points[i]->x);
                    if(insideBox3D(octree->points[i], boxMin, boxMax) && point->id != octree->points[i]->id)
                    {
                        *numInside += 1;
                        ptsInside = reallocWrap(ptsInside, *numInside * sizeof(Lpoint*));
                        ptsInside[*numInside-1] = octree->points[i];
                        *octrees = reallocWrap(*octrees, *numInside * sizeof(Octree));
                        (*octrees)[*numInside-1] = octree;
                    }
                }
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            if(!boxOverlap3D(boxMin, boxMax, octree->octants[i]))
            {
                continue;
            }
            else
            {
                ptsInside = neighborsNBH(point, boxMin, boxMax, octree->octants[i], ptsInside, octrees, numInside);
            }
        }
    }

    return ptsInside;
}

void neighborsOctant(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree oct, Octree **octrees, int *numInside)
{
    int i = 0;

    if(isLeaf(oct))
    {
        if(!isEmpty(oct))
        {
            if(boxOverlap3D(boxMin, boxMax, oct))
            {
                *octrees = reallocWrap(*octrees, ++*numInside * sizeof(Octree));
                (*octrees)[*numInside-1] = oct;
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            if(!boxOverlap3D(boxMin, boxMax, oct->octants[i]))
            {
                continue;
            }
            else
            {
                neighborsOctant(point, boxMin, boxMax, oct->octants[i], octrees, numInside);
            }
        }
    }
}

void neighborsOctantFelipe(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree oct, Octree **octrees, int *numInside)
{
    int i = 0;

    if(isLeaf(oct))
    {
        if(!isEmpty(oct))
        {
            if(boxOverlap2D(boxMin, boxMax, oct))
            {
                *octrees = reallocWrap(*octrees, ++*numInside * sizeof(Octree));
                (*octrees)[*numInside-1] = oct;
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            if(!boxOverlap2D(boxMin, boxMax, oct->octants[i]))
            {
                continue;
            }
            else
            {
                neighborsOctantFelipe(point, boxMin, boxMax, oct->octants[i], octrees, numInside);
            }
        }
    }
}

Octree* searchOctNeighbors(Octree oct, Octree octree, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint center;
    Octree *neighs = NULL;

    *numNeighs = 0;
    center.x = oct->center.x;
    center.y = oct->center.y;
    center.z = oct->center.z;
    makeBox(&center, oct->radius + 0.01, &boxMin, &boxMax); // Enlarge octant a bit
    neighborsOctant(&center, boxMin, boxMax, octree, &neighs, numNeighs);

    return neighs;
}

Octree* searchOctNeighborsFelipe(Lpoint *center, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Octree *neighs = NULL;

    *numNeighs = 0;
    // center.x = oct->center.x;
    // center.y = oct->center.y;
    // center.z = oct->center.z;
    makeBox(center, radius, &boxMin, &boxMax); // Enlarge octant a bit
    neighborsOctantFelipe(center, boxMin, boxMax, octree, &neighs, numNeighs);

    return neighs;
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

void makeBoxWrap(Lpoint *point, float radiusX, float radiusY, Vector3D *min, Vector3D *max)
{
    min->x = point->x - radiusX;
    min->y = point->y - radiusY;
    max->x = point->x + radiusX;
    max->y = point->y + radiusY;
}

// Make the box that represents an octree
void ownBox(Octree octree, Vector3D *min, Vector3D *max)
{
    min->x = octree->center.x - octree->radius;
    min->y = octree->center.y - octree->radius;
    min->z = octree->center.z - octree->radius;
    max->x = octree->center.x + octree->radius;
    max->y = octree->center.y + octree->radius;
    max->z = octree->center.z + octree->radius;
}

int numOctreeNeighbors3D(Lpoint *point, Octree octree, float radius)
{
    Vector3D boxMin, boxMax;
    int numInside = 0;

    makeBox(point, radius, &boxMin, &boxMax);
    numInside = numNeighbors3D(point, boxMin, boxMax, octree, &numInside);

    return numInside;
}

int numOctreeNeighbors2D(Lpoint *point, Octree octree, float radius)
{
    Vector3D boxMin, boxMax;
    int numInside = 0;

    makeBox(point, radius, &boxMin, &boxMax);
    numInside = numNeighbors2D(point, boxMin, boxMax, octree, &numInside);

    return numInside;
}

int numOctreeNeighbors2DDiffGroup(Lpoint *point, Octree octree, float radius)
{
    Vector3D boxMin, boxMax;
    int numInside = 0;

    makeBox(point, radius, &boxMin, &boxMax);
    numInside = numNeighbors2DDiffGroup(point, boxMin, boxMax, octree, &numInside);

    return numInside;
}

Lpoint** searchNeighbors2DGroupDistinct(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);
    ptsInside = neighbors2DDG(point, boxMin, boxMax, octree, ptsInside, numNeighs);

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

Lpoint** searchNeighbors2DWrap(Lpoint *point, Octree octree, float radiusX, float radiusY, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int ptsInside_size = 0;


    *numNeighs = 0;
    makeBoxWrap(point, radiusX, radiusY, &boxMin, &boxMax);
    ptsInside = neighbors2D(point, boxMin, boxMax, octree, ptsInside, &ptsInside_size, numNeighs);

    return ptsInside;
}

Lpoint** searchNeighbors2DResetId(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);
    ptsInside = neighbors2DReset(point, boxMin, boxMax, octree, ptsInside, numNeighs);

    return ptsInside;
}

Lpoint** searchNeighbors3D(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int ptsInside_size = 0;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);
    ptsInside = neighbors3D(point, boxMin, boxMax, octree, ptsInside,&ptsInside_size, numNeighs);
/*
    if(SORT)
    {
        int *newIdxs = NULL;
        newIdxs = nearestIndices(point, ptsInside, *numNeighs);
        sortPtsByFeature(ptsInside, *numNeighs, newIdxs);
    }
*/
    return ptsInside;
}

///** Inside a sphere & include query point */
//Lpoint** searchSphereNeighborhood(Lpoint *point, Octree octree, float radius, int *numNeighs)
//{
//    int i = 0, count = 0;
//    Vector3D boxMin, boxMax;
//    Lpoint **ptsInside = NULL, **ptsInside2 = NULL;
//
//    *numNeighs = 0;
//    makeBox(point, radius, &boxMin, &boxMax);
//    ptsInside = neighbors3D(point, boxMin, boxMax, octree, ptsInside, numNeighs);
//    ptsInside2 = malloc((*numNeighs + 1) * sizeof(Lpoint*));
//    ptsInside2[0] = point;
//    count = 1;
//    for(i = 0; i < *numNeighs; i++)
//        if(pointDist3D(point, ptsInside[i]) < radius)
//            ptsInside2[count++] = ptsInside[i];
//    free(ptsInside);
//    *numNeighs = count;
//    return ptsInside2;
//}

///** Inside a sphere */
//Lpoint** searchSphereNeighbors(Lpoint *point, Octree octree, float radius, int *numNeighs)
//{
//    int i = 0, count = 0;
//    Vector3D boxMin, boxMax;
//    Lpoint **ptsInsideVoxel = NULL, **ptsInsideSphere = NULL;
//
//    *numNeighs = 0;
//    makeBox(point, radius, &boxMin, &boxMax);
//    ptsInsideVoxel = neighbors3D(point, boxMin, boxMax, octree, ptsInsideVoxel, numNeighs);
//
//    ptsInsideSphere = malloc(*numNeighs * sizeof(Lpoint*));
//    for(i = 0; i < *numNeighs; i++)
//        if(pointDist3D(point, ptsInsideVoxel[i]) < radius)
//            ptsInsideSphere[count++] = ptsInsideVoxel[i];
//
//    ptsInsideSphere = reallocWrap(ptsInsideSphere, count * sizeof(Lpoint*));
//    *numNeighs = count;
//    free(ptsInsideVoxel);
//
//    return ptsInsideSphere;
//}

/** Ignore points inside the specified group */
Lpoint** searchNeighbors2DNoGroupLowerZ(Lpoint *point, Octree octree, float radius, int gId, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);
    ptsInside = neighbors2DNoGroupLowerZ(point, boxMin, boxMax, octree, ptsInside, gId, numNeighs);
    return ptsInside;
}

/** Ignore points inside the specified group */
Lpoint** searchNeighbors3DNoGroup(Lpoint *point, Octree octree, float radius, int gId, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);
    ptsInside = neighbors3DNoGroup(point, boxMin, boxMax, octree, ptsInside, gId, numNeighs);
    return ptsInside;
}

Lpoint** searchCircleNeighbors(Lpoint *point, Octree octree, float radius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int numInside = 0;
    int ptsInside_size = 0;

    makeBox(point, radius, &boxMin, &boxMax);
    ptsInside = circleNeighbors(point, boxMin, boxMax, octree, ptsInside,&ptsInside_size, &numInside, radius);
    *numNeighs = numInside;

    return ptsInside;
}


/** Exchange the point with the last one and decrement the octree size */
void hideOctPoint(Lpoint *p, Octree oct)
{
    Lpoint aux, *last = NULL;

    last = oct->points[oct->numPts-1];
    aux = *p;
    *p = *last;
    *last = aux;
    oct->numPts--;
}

void hideOctPointByIndex(int index, Octree oct)
{
    Lpoint *p = NULL, aux, *last = NULL;

    p = oct->points[index];
    last = oct->points[oct->numPts-1];
    aux = *p;
    *p = *last;
    *last = aux;
    oct->numPts--;
}

/** Undo the exchange done by the hideOctPoint() */
void recoverOctPoint(Lpoint *p, Lpoint *orig, Octree oct)
{
    Lpoint aux;

    aux = *p;
    *p = *orig;
    *orig = aux;
    oct->numPts++;
}

Octree searchOctree(Lpoint *point, Octree octree, Octree *leaf)
{
    int i = 0;
    Vector3D boxMin, boxMax;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(octree->points[i]->id == point->id)
                {
                    *leaf = octree;
                }
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            ownBox(octree->octants[i], &boxMin, &boxMax);

            if(!insideBox3D(point, boxMin, boxMax))
            {
                continue;
            }
            else
            {
                *leaf = searchOctree(point, octree->octants[i], leaf);
            }
        }
    }

    return *leaf;
}

void printOctreeCenters(Octree octree, FILE *file, int alsoEmpty)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            fprintf(file, "%.2f\t%.2f\t%.2f\t%.2f\n", octree->center.x, octree->center.y, octree->center.z, octree->radius);
            return;
        }
    }
    else
    {
        if(alsoEmpty) fprintf(file, "%.2f\t%.2f\t%.2f\t%.2f\n", octree->center.x, octree->center.y, octree->center.z, octree->radius);

        for(i = 0; i < 8; i++)
        {
            printOctreeCenters(octree->octants[i], file, alsoEmpty);
        }
    }
}

void printOctreePoints(Octree octree, FILE *file)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(file == NULL) printf("%d\t%.2f\t%.2f\t%.2f\n", octree->points[i]->id, octree->points[i]->x, octree->points[i]->y, octree->points[i]->z);
                else fprintf(file, "%d\t%.2f\t%.2f\t%.2f\n", octree->points[i]->id, octree->points[i]->x, octree->points[i]->y, octree->points[i]->z);
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            printOctreePoints(octree->octants[i], file);
        }
    }
}

// void writeOctreeNormals(Octree octree, FILE *file)
// {
//     int i = 0;
//
//     if(isLeaf(octree))
//     {
//         if(!isEmpty(octree))
//         {
//             //if(octree->numPoints >= 10)
//             //{
//                 //fprintf(file, "%.2f\t%.2f\t%.2f\n", octree->center.x, octree->center.y, octree->points[octree->numPoints/2]->z);
//                 fprintf(file, "%.2f\t%.2f\t%.2f\n", octree->points[octree->numPts/2]->x, octree->points[octree->numPts/2]->y, octree->points[octree->numPts/2]->z);
//                 fprintf(file, "%.2f\t%.2f\t%.2f\n", octree->plane[0], octree->plane[1], octree->plane[2]);
//             //}
//         }
//     }
//     else
//     {
//         for(i = 0; i < 8; i++)
//         {
//             writeOctreeNormals(octree->octants[i], file);
//         }
//     }
// }

int numOctreePoints(Octree octree, int *num)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                ++*num;
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            *num = numOctreePoints(octree->octants[i], num);
        }
    }

    return *num;
}

void searchMinMaxOct(Octree octree, double *min, double *max)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(octree->points[i]->intensity > *max)
                {
                    *max = octree->points[i]->intensity;
                }
                if(octree->points[i]->intensity < *min)
                {
                    *min = octree->points[i]->intensity;
                }
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            searchMinMaxOct(octree->octants[i], min ,max);
        }
    }
}

void searchMinFelipe(Octree octree, double* min, unsigned int* id)
{
    //Devuelvo el índice
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(octree->points[i]->z < *min)
                {
                    *min = octree->points[i]->z;
                    *id = octree->points[i]->id;
                }
            }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            searchMinFelipe(octree->octants[i], min, id);
        }
    }
}

// void normalizeIntOct(Octree octree, double min, double max)
// {
//     int i = 0;
//     double range = 0;
//
//     range = max - min; //printf("Range %.2lf Min %.2lf Max %.2lf\n", range, min, max);
//
//     if(isLeaf(octree))
//     {
//         if(!isEmpty(octree))
//         {
//            for(i = 0; i < octree->numPts; i++)
//            {
//                octree->points[i]->intensity = ((octree->points[i]->intensity - min) / range);
//            }
//         }
//     }
//     else
//     {
//         for(i = 0; i < 8; i++)
//         {
//             normalizeIntOct(octree->octants[i], min ,max);
//         }
//     }
// }

void cpyOctNoDiv(Octree dst, Octree src)
{
    int i = 0;

    if(isLeaf(src))
    {
        if(!isEmpty(src))
        {
           for(i = 0; i < src->numPts; i++)
           {
                insertPointNoSplit(src->points[i], dst);
           }
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            cpyOctNoDiv(dst, src);
        }
    }
}

// void normalizeOctIntensity(Octree octree)
// {
//     double min = DBL_MAX, max = -DBL_MAX;
//
//     searchMinMaxOct(octree, &min, &max); printf("Min %f max %f\n", min, max);
//     normalizeIntOct(octree, min, max);
// }

void sortOctsByFeature(Octree *octrees, int num, int* sortIdxs)
{
    int i = 0;
    Octree *aux = NULL;

    aux = mallocWrap(num * sizeof(Lpoint**));

    for(i = 0; i < num; i++)
    {
        aux[i] = octrees[sortIdxs[i]];
    }

    memcpy(octrees, aux, num * sizeof(Lpoint*));
    free(aux);
}

// Get centroid point of octant
Lpoint octCentroid (Octree oct)
{
    int i = 0;
    double x = 0, y = 0, z = 0, inten = 0;
    Lpoint centroid;

    for(i = 0; i < oct->numPts; i++)
    {
        x += oct->points[i]->x;
        y += oct->points[i]->y;
        z += oct->points[i]->z;
        inten += oct->points[i]->intensity;
    }

    x /= oct->numPts;
    y /= oct->numPts;
    z /= oct->numPts;
    inten /= oct->numPts;

    centroid.x = x;
    centroid.y = y;
    centroid.z = z;
    centroid.intensity = inten;

    return centroid;
}

void octPadding() // On i5 x64 tends to pad until divisible by 8 (currently struct is 136 bytes w/o padding)
    {
        int membersSizeDou = sizeof(Octree) * 8 + sizeof(double) * 7 + sizeof(Lpoint**) + sizeof(float) + sizeof(int);
        int membersSizeVec = sizeof(Octree) * 8 + sizeof(Vector3D) + sizeof(Lpoint**) + sizeof(float) + sizeof(int) + sizeof(double) * 4;
        int size = sizeof(struct Octree);
        printf("SIZE %d bytes DBL PAD: %d V3C PAD: %d\n", size, size - membersSizeDou, size - membersSizeVec);
    }

void octreeSpacing(Octree oct, int *numPoints, double *area, FILE *file)
{
    int i = 0;

    if(isLeaf(oct))
    {
        if(!isEmpty(oct))
        {
           *numPoints += oct->numPts;
           *area += (oct->radius * 2) * (oct->radius * 2);
           fprintf(file, "%.2f\t%.2f\t%.2f\t%.2f\n", oct->center.x, oct->center.y, oct->center.z, oct->radius);
        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            octreeSpacing(oct->octants[i], numPoints, area, file);
        }
    }
}

void removeOctree(Octree oct)
{
    free(oct->points);
    oct->points = NULL;
    oct->numPts = 0;
}

//// Create and populate an octree without splitting the points in octants
//Octree buildOctreeUnsplit(Vector3D center, float radius, Lpoint *points, int numPts, char label)
//{
//    int i = 0;
//    Octree oct = NULL;
//
//    printf("Building octree unsplit ...\n");
//    oct = createOctree(center, radius, NULL);
//    for(i = 0; i < numPts; i++)
//        if(label == -1 || points[i].class == label)
//            insertPointNoSplit(&points[i], oct);
//
//    return oct;
//}

// void gatherLeafs(Octree oct, Octree **fitLeafs, int *numFit, Octree **unFitLeafs, int *numUnFit)
// {
//     int i = 0;
//
//     if(isLeaf(oct))
//     {
//         if(!isEmpty(oct))
//         {
//             if(oct->plane[3] == -1.0)
//             {
//                 *unFitLeafs = reallocWrap(*unFitLeafs, ++*numUnFit * sizeof(Octree));
//                 (*unFitLeafs)[*numUnFit-1] = oct;
//             }
//             else if(oct->numPts > 2)
//             {
//                 *fitLeafs = reallocWrap(*fitLeafs, ++*numFit * sizeof(Octree));
//                 (*fitLeafs)[*numFit-1] = oct;
//             }
//         }
//     }
//     else
//     {
//         for(i = 0; i < 8; i++)
//         {
//             gatherLeafs(oct->octants[i], fitLeafs, numFit, unFitLeafs, numUnFit);
//         }
//     }
// }

//!< Check if a point is inside an octree
int insideOct(Lpoint *p, Octree oct)
{
    Vector3D boxMin, boxMax;
    Lpoint octCenter;

    octCenter.x = oct->center.x;
    octCenter.y = oct->center.y;
    octCenter.z = oct->center.z;
    makeBox(&octCenter, oct->radius, &boxMin, &boxMax);

    return insideBox3D(p, boxMin, boxMax);
}

/*
 LANDING
 */

Lpoint** donutNeighbors(Lpoint *point, Vector3D boxMin, Vector3D boxMax, Octree octree, Lpoint **ptsInside, int *ptsInside_size, int *numInside, float extRadius, float intRadius)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
            for(i = 0; i < octree->numPts; i++)
            {
                if(insideCircle(octree->points[i], point, extRadius) & (point->id != octree->points[i]->id) )
                {
			if(insideCircle(octree->points[i], point, intRadius) & (point->id != octree->points[i]->id) )
			{
				//do nothing
			}else{
                            if (*numInside >= *ptsInside_size) {
                                (*ptsInside_size) += REALLOC_INCREMENT;
                                ptsInside = reallocWrap(ptsInside, *ptsInside_size * sizeof(Lpoint*));
                            }
                            ptsInside[(*numInside)++] = octree->points[i];
			}
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
                ptsInside = donutNeighbors(point, boxMin, boxMax, octree->octants[i], ptsInside, ptsInside_size, numInside, extRadius, intRadius);
            }
        }
    }

    return ptsInside;
}


Lpoint** searchDonutNeighbors(Lpoint *point, Octree octree, float extRadius, float intRadius, int *numNeighs)
{
    Vector3D boxMin, boxMax;
    Lpoint **ptsInside = NULL;
    int numInside = 0;
    int ptsInside_size = 0;

    makeBox(point, extRadius, &boxMin, &boxMax);
    ptsInside = donutNeighbors(point, boxMin, boxMax, octree, ptsInside,&ptsInside_size, &numInside, extRadius, intRadius);
    *numNeighs = numInside;

    return ptsInside;
}

void getLeafMinRadius(Octree octree, float *radius)
{
    int i = 0;

    if(isLeaf(octree))
    {
        if(!isEmpty(octree))
        {
           if(octree->radius < *radius)
                {
                    *radius = octree->radius ;
                }

        }
    }
    else
    {
        for(i = 0; i < 8; i++)
        {
            getLeafMinRadius(octree->octants[i], radius);
        }
    }
}

void printLeafRadius(Octree octree){

    float radius=octree->radius;
    getLeafMinRadius(octree,&radius);

    printf("Octree min leaf size: %.3f\n",radius);
}
