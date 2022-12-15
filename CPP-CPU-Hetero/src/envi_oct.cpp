#include "../include/envi.h"
#include "tbb/parallel_for.h"
#include <numeric>

double round2d(double z){
  return round(z*100)/100;
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

// int boxOverlap2D(Vector3D boxMin, Vector3D boxMax, Octree_t* oct)
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

int boxTotalOverlap2D(Vector3D boxMin, Vector3D boxMax, Octree oct)
{
    if(oct->center.x + oct->radius < boxMax.x ||
       oct->center.y + oct->radius < boxMax.y)
        return 0;

    if(oct->center.x - oct->radius > boxMin.x ||
       oct->center.y - oct->radius > boxMin.y)
        return 0;

    return 1;
}

int isLeaf(Octree oct)
{
    return oct->octants[0] == NULL;
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

// Find the child corresponding a given point
int octantIdx(Lpoint *point, Octree octree)
{
    int child = 0;

    if(point->x >= octree->center.x) child |= 4;
    if(point->y >= octree->center.y) child |= 2;
    if(point->z >= octree->center.z) child |= 1;

    return child;
}


// Lpoint findValidMin(Octree octree, const std::vector<std::function<bool(Vector3D, float)>>& conditions)
Lpoint findValidMin(Octree octree, Vector3D* boxMin, Vector3D* boxMax, int* numInside)
{
    // Lpoint tmp, min = nomin;
    Lpoint tmp, min = {0,0.0,0.0,99999.0};
    int i;

    if(isLeaf(octree))
    {
        // printf("octree with %d points\n", octree->numPts);
        // printf("nodo hoja con radio %g\n", octree->radius/2.0);
        uint64_t size = octree->points.size();
        if(size > 0)
        {
            for (i = 0; i < size; i++) {
              if(insideBox2D(octree->points[i], *boxMin, *boxMax))
              {
                tmp = *(octree->points[i]);
                if (tmp.z < min.z) {
                    min = tmp;
                }
                (*numInside)++;
              }
            }

        }

    } else {

        for(i = 0; i < 8; i++) {
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

Lpoint searchNeighborsMin(Lpoint* point, Octree octree, float radius, int* numNeighs)
{
    Vector3D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    return findValidMin(octree, &boxMin, &boxMax, numNeighs);
}

void countNeighbors(Octree octree, Vector3D* boxMin, Vector3D* boxMax, int* numInside)
{
    int i;

    if(isLeaf(octree))
    {
        uint64_t size = octree->points.size();
        if(size > 0)
        {
            for (i = 0; i < size; i++)
            {
              if(insideBox2D(octree->points[i], *boxMin, *boxMax))
              {
                (*numInside)++;
              }
            }

        }

    } else {

        for(i = 0; i < 8; i++) {
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

void countNeighbors2D(Lpoint* point, Octree octree, float radius, int* numNeighs)
{
    Vector3D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    countNeighbors(octree, &boxMin, &boxMax, numNeighs);

    return;
}

Lpoint findValidMinPath(Octree octree, Vector3D* boxMin, Vector3D* boxMax,
                    int* numInside, std::vector<Octree>& path)
{
  Lpoint tmp, min = {0,0.0,0.0,99999.0};
  int i;

  if(isLeaf(octree))
  {
    // printf("octree with %d points\n", octree->numPts);
    // printf("nodo hoja con radio %g\n", octree->radius/2.0);
    uint64_t size = octree->points.size();
    if(size > 0)
    {
      for (i = 0; i < size; i++) {
        if(insideBox2D(octree->points[i], *boxMin, *boxMax))
        {
          tmp = *(octree->points[i]);
          if (tmp.z < min.z) {
              min = tmp;
          }
          (*numInside)++;
        }
      }
    }

  } else {
    std::vector<Octree> local_path, best_path;

    for(i = 0; i < 8; i++) {
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

Lpoint searchNeighborsMinPath(Lpoint* point, Octree octree, float radius, int* numNeighs, std::vector<Octree>& path)
{
    Vector3D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    return findValidMinPath(octree, &boxMin, &boxMax, numNeighs, path);
}


Lpoint findValidMinParent(Octree octree, Vector3D* boxMin, Vector3D* boxMax,
                    int* numInside, Octree* minparent, int nivel)
{
  Lpoint tmp, min = {0,0.0,0.0,99999.0};
  int i;

  if(isLeaf(octree))
  {
    // printf("octree with %d points\n", octree->numPts);
    // printf("nodo hoja con radio %g\n", octree->radius/2.0);
    uint64_t size = octree->points.size();
    if(size > 0)
    {
      for (i = 0; i < size; i++) {
        if(insideBox2D(octree->points[i], *boxMin, *boxMax))
        {
          tmp = *(octree->points[i]);
          if (tmp.z < min.z) {
              min = tmp;
              // printf("min level %d\n", nivel);
          }
          (*numInside)++;
        }
      }
    }

  } else {
    for(i = 0; i < 8; i++) {
      // Check
      if(!boxOverlap2D(*boxMin, *boxMax, octree->octants[i]))
        continue;
      else {
        tmp = findValidMinParent(octree->octants[i], boxMin, boxMax, numInside, minparent, nivel+1);
        if (tmp.z < min.z) {
          min = tmp;
          if(isLeaf(octree->octants[i]))
            *minparent = octree->octants[i];
        }
      }
    }
  }

    return min;
}


Lpoint searchNeighborsMinParent(Lpoint* point, Octree* oparent, float radius, int* numNeighs)
{
    Vector3D boxMin, boxMax;
    Octree actual_parent = *oparent;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    while(actual_parent->oparent != NULL && !boxTotalOverlap2D(boxMin, boxMax, actual_parent)){
    // while(actual_parent->oparent != NULL){
      actual_parent = actual_parent->oparent;
    }

    // if(actual_parent->oparent != NULL) {
    //   // printf("my size: %zu\n", (*oparent)->points.size());
    //   actual_parent = actual_parent->oparent;
    //   // printf("parent size: %zu\n", (*oparent)->points.size());
    //   // if(actual_parent->oparent != NULL)
    //   //   actual_parent = actual_parent->oparent;
    // }

    return findValidMinParent(actual_parent, &boxMin, &boxMax, numNeighs, oparent, 0);
}


/** Create a octree with the given center and radius */
// Octree createOctreeI(Vector3D center, float radius)
// {
//     int i = 0;
//     Octree oct = NULL;
//
//     oct = (Octree_t*)mallocWrap(sizeof(Octree_t));
//     // oct = std::unique_ptr<Octree_t>{ (Octree_t*)mallocWrap(sizeof(Octree_t)) };
//     oct->center = center;
//     oct->radius = radius;
//     // oct->points = NULL;
//     // oct->numPts = 0;
//     for(i = 0; i < 8; i++)
//       oct->octants[i] = NULL;
//
//     return oct;
// }

// uOctree createOctreeF(Vector3D center, float radius)
// {
//     int i = 0;
//     // uOctree oct = uOctree{(Octree_t*)malloc(sizeof(Octree_t))};
//     // uOctree oct( new Octree_t );
//     Octree oct( new Octree_t );
//
//     // oct = (Octree_t*)mallocWrap(sizeof(Octree_t));
//     // oct = std::unique_ptr<Octree_t>{ (Octree_t*)mallocWrap(sizeof(Octree_t)) };
//     // oct = std::unique_ptr<Octree_t>{ new Octree_t };
//     oct->center = center;
//     oct->radius = radius;
//     oct->points = NULL;
//     oct->numPts = 0;
//     for(i = 0; i < 8; i++)
//       oct->octants[i] = NULL;
//
//     return oct;
// }


void createOctantsF(Octree oct)
{
    int i = 0;
    Vector3D newCenter;
    float newRadius = oct->radius * 0.5;

    for(i = 0; i < 8; i++)
    {
        newCenter = oct->center;
        newCenter.x += oct->radius * (i&4 ? 0.5f : -0.5f);
        newCenter.y += oct->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.z += oct->radius * (i&1 ? 0.5f : -0.5f);
        // printf("(%lf,%lf,%lf) ", newCenter.x, newCenter.y, newCenter.z);
        // oct->octants[i] = std::move(createOctreeF(newCenter, oct->radius * 0.5));
        // oct->octants[i] = createOctreeF(newCenter, oct->radius * 0.5);
        // oct->octants[i] = uOctree{ new nHoja2(newCenter, oct->radius * 0.5) };
        // (oct->octants[i]).reset( new nHoja2(newCenter, oct->radius * 0.5) );
        oct->octants[i] = new Octree_t(oct, newCenter, newRadius);
    }
    // printf("\n");
}


// Insert a point in the octree creating the appropiate childs
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
          idx = octantIdx(point, octree);
          insertPointF(point, octree->octants[idx], minRadius);

        } else {
          octree->points.push_back(point);
        }
    }
    else                                // No leaf -> search the correct one
    {
      // printf("octante intermedio nivel %d\n", nivel);
      idx = octantIdx(point, octree);
      insertPointF(point, octree->octants[idx], minRadius);
      // insertPoint(point, octree->octants[idx]);
    }
}

void deleteOctree(Octree octree)
{
    int i;
    if(isLeaf(octree))
    {
      octree->points.clear();
    } else {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(i = 0; i < 8; i++) {
            // Check
            deleteOctree(octree->octants[i]);
            delete(octree->octants[i]);
            // (octree->octants[i]).reset(NULL); // free memory for unique_ptr
        }

    }

    return;
}

unsigned int stage1cpp(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min)
{


    double initX = min.x + Wsize/2 - Wsize*Overlap;
    double initY = min.y + Wsize/2 - Wsize*Overlap;

    double Displace = round2d(Wsize*(1-Overlap));

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    // Lpoint** neighbors = NULL;
    Lpoint newmin;

    unsigned int countMin;

    int cellPoints;

    // std::vector<int> vcount(Ccol,0);

    #pragma omp parallel for private(cellCenter,cellPoints,newmin) schedule(dynamic) reduction(+:countMin)
    for(int jj = 0 ; jj < Ccol ; jj++ ){

        // countMin=0; // counts mins per row

        cellCenter.y = initY + jj*Displace;
        // printf("\nCeldas no descartadas thread %d:   %d\n",omp_get_thread_num(), countMin);

        for(int ii=0 ; ii < Crow ; ii++ ){

            // std::vector<std::function<bool(Vector3D,float)>> conditions = {
            //     [&diag,&cellCenter] (Vector3D c, float r) {
            //         double d = sqrt(pow(cellCenter.x - c.x,2) + pow(cellCenter.y - c.y,2));
            //         return d < (diag+r*sqrt(2)); // suma de las diagonales
            //     },
            // };

            cellCenter.x = initX + ii*Displace;
            // printf("Centro de %d: %.2f %.2f\n",omp_get_thread_num(), cellCenter.x, cellCenter.y);
            // printf("Busco los vecinos\n");
            newmin = searchNeighborsMin(&cellCenter, octreeIn, Wsize/2, &cellPoints);
            // printf("Numero de elementos de la celda: %d\n", cellPoints );
            if(cellPoints >= minNumPoints ){

                minIDs[jj*Crow+ii] = newmin.id; //los voy guardando por posición de celda; esto hace que tenga que modificar stage2
                countMin++;
            }
        }
        // vcount[jj] = countMin;
        // countMin=0;
    }

    // countMin = std::accumulate(vcount.begin(), vcount.end(), 0);

    return countMin;
}

unsigned int stage1cppParent(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min)
{


    double initX = min.x + Wsize/2 - Wsize*Overlap;
    double initY = min.y + Wsize/2 - Wsize*Overlap;

    double Displace = round2d(Wsize*(1-Overlap));

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    // Lpoint** neighbors = NULL;
    Lpoint newmin;

    unsigned int countMin = 0;

    int cellPoints;

    // std::vector<int> vcount(Ccol,0);
    Octree mypointer;

    #pragma omp parallel for private(cellCenter,cellPoints,newmin,mypointer) schedule(dynamic) reduction(+:countMin)
    for(int jj = 0 ; jj < Ccol ; jj++ ){
        /* lo renuevo para cada fila porque cada vez que cambio de fila voy
        a estar lejos de donde estaba antes en el mapa*/
        mypointer = octreeIn;

        cellCenter.y = initY + jj*Displace;
        // printf("\nCeldas no descartadas thread %d:   %d\n",omp_get_thread_num(), countMin);

        for(int ii=0 ; ii < Crow ; ii++ ){

            // std::vector<std::function<bool(Vector3D,float)>> conditions = {
            //     [&diag,&cellCenter] (Vector3D c, float r) {
            //         double d = sqrt(pow(cellCenter.x - c.x,2) + pow(cellCenter.y - c.y,2));
            //         return d < (diag+r*sqrt(2)); // suma de las diagonales
            //     },
            // };

            cellCenter.x = initX + ii*Displace;
            // printf("Centro de %d: %.2f %.2f\n",omp_get_thread_num(), cellCenter.x, cellCenter.y);
            // printf("Busco los vecinos\n");
            // newmin = searchNeighborsMin(&cellCenter, octreeIn, Wsize/2, &cellPoints);
            newmin = searchNeighborsMinParent(&cellCenter, &mypointer, Wsize/2, &cellPoints);
            // printf("Numero de elementos de la celda: %d\n", cellPoints );
            if(cellPoints >= minNumPoints ){

                minIDs[jj*Crow+ii] = newmin.id; //los voy guardando por posición de celda; esto hace que tenga que modificar stage2
                countMin++;
            }
        }
        // vcount[jj] = countMin;
        // countMin=0;
    }

    // countMin = std::accumulate(vcount.begin(), vcount.end(), 0);

    return countMin;
}

int stage1tbb(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minIDs, Octree octreeIn, Vector3D min)
{


    double initX = min.x + Wsize/2 - Wsize*Overlap;
    double initY = min.y + Wsize/2 - Wsize*Overlap;

    double Displace = round2d(Wsize*(1-Overlap));

    // std::vector<int> vcount(Ccol,0);

    tbb::parallel_for(0, static_cast<int>(Ccol), 1,
        [&](int jj) {

        /* lo renuevo para cada fila porque cada vez que cambio de fila voy
        a estar lejos de donde estaba antes en el mapa*/
        Lpoint newmin;
        Lpoint cellCenter;
        int cellPoints;
        Octree mypointer = octreeIn;

        cellCenter.y = initY + jj*Displace;

        for(int ii=0 ; ii < Crow ; ii++ ){

          cellCenter.x = initX + ii*Displace;
          newmin = searchNeighborsMinParent(&cellCenter, &mypointer, Wsize/2, &cellPoints);
          if(cellPoints >= minNumPoints ){
              minIDs[jj*Crow+ii] = newmin.id; //los voy guardando por posición de celda; esto hace que tenga que modificar stage2
              // countMin++;
          }
        }
    });

    // finalcount = std::accumulate(vcount.begin(), vcount.end(), 0);

    return 0;
}



unsigned int stage2(unsigned int countMin, int* minIDs){

    unsigned int index = 0;

    int ii,jj,id;

    for( ii=0 ; ii<countMin ; ii=jj ){
        id = minIDs[ii];

        for( jj=ii+1 ; id==minIDs[jj] && jj<countMin ; jj++ );

        /* esta es la forma de descartar las posiciones no utilizadas, inicializadas a -1 */
        if((id >= 0) && (jj-ii > 1)){
            // if(jj-ii > 1){
            minIDs[index]=id;
            index++;
            // }
        }
    }

    return index;
}


unsigned int stage3cpp(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          int* minGridIDs, Octree octreeIn, Octree grid, Vector3D min){

    unsigned int addMin = 0;

    int cellPoints;

    Lpoint cellCenter = {0,0.0,0.0,0.0};

    // Lpoint** neighbors = NULL;

    // Lpoint** minimos = NULL;

    Lpoint newmin;

    int ii,jj;

    Octree mypointer;

    // #pragma omp parallel for private(cellCenter,cellPoints,newmin) schedule(static)
    #pragma omp parallel for private(cellCenter,cellPoints,newmin,mypointer) schedule(dynamic) reduction(+:addMin)
    for(int jj = 0 ; jj < Ccol ; jj++ ){

      mypointer = octreeIn;

      cellCenter.y = min.y + Bsize/2 + jj*Bsize;

      for(int ii = 0 ; ii < Crow ; ii++ ){
          cellCenter.x = min.x + Bsize/2 + ii*Bsize;

          // Candidata a VOID porque no necesito el mínimo, solo el número de puntos encontrados.
          // newmin = searchNeighborsMin(&cellCenter, grid, Bsize/2, &cellPoints);
          countNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);
          // printf("Numero de elementos de la celda: %d\n", cellPoints );
          // Tengo que hacerlo porque el algoritmo no lo hace
          if(cellPoints == 0){
              // printf("    Tengo una celda vacia en la malla\n" );
              newmin = searchNeighborsMin(&cellCenter, octreeIn, Bsize/2, &cellPoints);
              // newmin = searchNeighborsMinParent(&cellCenter, &mypointer, Bsize/2, &cellPoints);
              if(cellPoints>0){

                  // idmin = findMin(neighbors, cellPoints);
                  // #pragma omp critical
                  // {
                  //     minGridIDs[addMin] = newmin.id;
                  //     addMin++;
                  // }
                  minGridIDs[jj*Crow+ii] = newmin.id; //los voy guardando por posición de celda; esto hace que tenga que modificar stage2
                  addMin++;
                  // printf("%d: %.2f , %.2f , %.2f\n", addMin, pointer[neighbors[idmin]->id].x, pointer[neighbors[idmin]->id].y,pointer[neighbors[idmin]->id].z);
              }
          }
      }
    }

    return addMin;

}


// void prueba1( uOctree2& octree )
// {
//   printf("myfloat %g\n", octree->radius);
//   octree.reset( new nHoja2({0.0,0.0,0.0}, 33.0) );
//   // printf("myfloat %g\n", (octree)->radius);
//   // octree.reset();
//   return;
// }
