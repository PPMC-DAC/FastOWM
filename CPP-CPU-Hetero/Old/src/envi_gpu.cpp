#include "../include/envi_gpu.h"

#include "tbb/parallel_reduce.h"
#include "tbb/concurrent_vector.h"
#include "tbb/blocked_range2d.h"

#include <atomic>

#include <fstream>
#include <sstream>


extern unsigned short  minNumPoints, Crow, Ccol, Crowg, Ccolg;

extern unsigned short Wsize;

extern Qtree qtreeIn;

extern double Displace;


double round2d(double z){
  return round(z*100.0)/100.0;
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

void makeBox(Vector2D *point, double radiusX, double radiusY, Vector2D *min, Vector2D *max)
{

    min->x = point->x - radiusX;
    min->y = point->y - radiusY;

    max->x = point->x + radiusX;
    max->y = point->y + radiusY;

}

int boxInside2D(Vector2D boxMin, Vector2D boxMax, Qtree qt)
{
    if(qt->center.x + qt->radius > boxMax.x ||
       qt->center.y + qt->radius > boxMax.y)
        return 0;

    if(qt->center.x - qt->radius < boxMin.x ||
       qt->center.y - qt->radius < boxMin.y)
        return 0;

    return 1;
}

int insideBox2D(Vector3D *point, Vector2D min, Vector2D max)
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


// int boxOverlap2D(Vector2D boxMin, Vector2D boxMax, Qtree_t* qt)
int boxOverlap2D(Vector2D boxMin, Vector2D boxMax, Qtree qt)
{
    if(qt->center.x + qt->radius < boxMin.x ||
       qt->center.y + qt->radius < boxMin.y)
        return 0;

    if(qt->center.x - qt->radius > boxMax.x ||
       qt->center.y - qt->radius > boxMax.y)
        return 0;

    return 1;
}

int boxOverlap2D(Vector2D boxMin, Vector2D boxMax, QtreeG qt)
{
    if(qt->center.x + qt->radius < boxMin.x ||
       qt->center.y + qt->radius < boxMin.y)
        return 0;

    if(qt->center.x - qt->radius > boxMax.x ||
       qt->center.y - qt->radius > boxMax.y)
        return 0;

    return 1;
}

int boxTotalOverlap2D(Vector2D boxMin, Vector2D boxMax, Qtree qt)
{
    if(qt->center.x + qt->radius < boxMax.x ||
       qt->center.y + qt->radius < boxMax.y)
        return 0;

    if(qt->center.x - qt->radius > boxMin.x ||
       qt->center.y - qt->radius > boxMin.y)
        return 0;

    return 1;
}

int isLeaf(Qtree qt)
{
    return qt->quadrants[0] == NULL;
}

int isLeaf(QtreeG qt)
{
    return qt->quadrants[0] == NULL;
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

// Find the child corresponding a given point
int quadrantIdx(Lpoint *point, Qtree qtree)
{
    int child = 0;

    if(point->x >= qtree->center.x) child |= 2;
    if(point->y >= qtree->center.y) child |= 1;
    // if(point->z >= qtree->center.z) child |= 1;

    return child;
}

int quadrantIdx(Lpoint *point, QtreeG qtree)
{
    int child = 0;

    if(point->x >= qtree->center.x) child |= 2;
    if(point->y >= qtree->center.y) child |= 1;
    // if(point->z >= qtree->center.z) child |= 1;

    return child;
}


// Lpoint findValidMin(Qtree qtree, const std::vector<std::function<bool(Vector2D, float)>>& conditions)
Lpoint findValidMin(Qtree qtree, Vector2D* boxMin, Vector2D* boxMax, int* numInside)
{
    // Lpoint tmp, min = nomin;
    Lpoint tmp, min = {0,0.0,0.0,99999.0};

    if(isLeaf(qtree))
    {

      if(boxInside2D(*boxMin, *boxMax, qtree)){
        for(Lpoint* p : qtree->points) {
          if (p->z < min.z) {
              min = *p;
          }
          (*numInside)++;
        }
      } else {
        for(Lpoint* p : qtree->points) {
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

        for(int i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(*boxMin, *boxMax, qtree->quadrants[i]))
                continue;
            else {
                tmp = findValidMin(qtree->quadrants[i], boxMin, boxMax, numInside);
                if (tmp.z < min.z) {
                    min = tmp;
                }
            }

        }

    }

    return min;
}

Lpoint searchNeighborsMin(Vector2D* point, Qtree qtree, float radius, int* numNeighs)
{
    Vector2D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    return findValidMin(qtree, &boxMin, &boxMax, numNeighs);
}


Lpoint searchOverlap(Vector2D* point, Qtree qtree, double radiusX, double radiusY, int* numNeighs)
{
    Vector2D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radiusX, radiusY, &boxMin, &boxMax);

    return findValidMin(qtree, &boxMin, &boxMax, numNeighs);
}

void countNeighbors(Qtree qtree, Vector2D* boxMin, Vector2D* boxMax, int* numInside)
{
    int i;

    if(isLeaf(qtree))
    {
      if(boxInside2D(*boxMin, *boxMax, qtree)){

          (*numInside) += qtree->points.size();

      } else {
        for(Lpoint* p : qtree->points) {
          if(insideBox2D(p, *boxMin, *boxMax))
          {
            (*numInside)++;
          }
        }
      }

    } else {

        for(i = 0; i < 4; i++) {
            // Check
            if(!boxOverlap2D(*boxMin, *boxMax, qtree->quadrants[i]))
              continue;
            else {
              countNeighbors(qtree->quadrants[i], boxMin, boxMax, numInside);
            }

        }

    }

    return;
}

void countNeighbors2D(Vector2D* point, Qtree qtree, float radius, int* numNeighs)
{
    Vector2D boxMin, boxMax;

    *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    countNeighbors(qtree, &boxMin, &boxMax, numNeighs);

    return;
}

Qtree createQtreeF(Qtree parent, Vector2D center, float radius)
{
    Qtree qt = new Qtree_t;

    qt->center = center;
    qt->radius = radius;
    qt->parent = parent;

    qt->quadrants.reserve(4);
    
    for( int i = 0; i < 4; i++)
      qt->quadrants[i] = NULL;

    return qt;
}


void createQuadrantsF(Qtree qt)
{
    Vector2D newCenter;
    float newRadius = qt->radius * 0.5;

    for( int i = 0; i < 4; i++)
    {
        newCenter = qt->center;
        newCenter.x += qt->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += qt->radius * (i&1 ? 0.5f : -0.5f);

        qt->quadrants[i] = createQtreeF(qt, newCenter, newRadius);

    }
    // printf("\n");
}


// Insert a point in the qtree creating the appropiate childs
void insertPointF(Lpoint *point, Qtree qtree, float minRadius)
{
    int idx = 0;

    if(isLeaf(qtree))
    {
        // printf("octante hoja nivel %d\n",nivel);
        if(qtree->radius * 0.5 > minRadius)    // still divisible -> divide
        {
          createQuadrantsF(qtree);
          // fillOctants(qtree);
          idx = quadrantIdx(point, qtree);
          insertPointF(point, qtree->quadrants[idx], minRadius);

        } else {
          qtree->points.push_back(point);
        }
    }
    else                                // No leaf -> search the correct one
    {
      // printf("octante intermedio nivel %d\n", nivel);
      idx = quadrantIdx(point, qtree);
      insertPointF(point, qtree->quadrants[idx], minRadius);
    }
}

void fillQuadrants(Qtree qtree, float minRadius)
{
    int idx;

    for(Lpoint* p : qtree->points)
    {
      idx = quadrantIdx(p, qtree);
      insertPointF(p, qtree->quadrants[idx], minRadius);
    }
    qtree->points.clear();
}

void deleteQtree(Qtree qtree)
{
    if(isLeaf(qtree))
    {
      qtree->points.clear();
      qtree->quadrants.clear();
      // delete(qtree);
    } else {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(int i = 0; i < 4; i++) {
            // Check
            deleteQtree(qtree->quadrants[i]);
            delete(qtree->quadrants[i]);
            // qtree->quadrants.clear();
            // (qtree->quadrants[i]).reset(NULL); // free memory for unique_ptr
        }
        qtree->quadrants.clear();
        // delete(qtree);
    }

    return;
}

void deleteQtreeGPU(QtreeG qtree)
{
    if(isLeaf(qtree))
    {
      qtree->points.clear();
      free(qtree->quadrants, q);
    } 
    else 
    {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(int i = 0; i < 4; ++i) {
            // Check
            deleteQtreeGPU(qtree->quadrants[i]);
            free(qtree->quadrants[i], q);
            // delete(qtree->quadrants[0]);
            // qtree->quadrants.clear();
            // (qtree->quadrants[i]).reset(NULL); // free memory for unique_ptr
        }
        free(qtree->quadrants, q);
    }

    return;
}

void deleteQtreeGPU2(QtreeG qtree) // Este tengo que utilizarlo cuando creo los nodos contiguos
{
    if(isLeaf(qtree))
    {
      qtree->points.clear();
      free(qtree->quadrants, q);
    } 
    else 
    {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(int i = 0; i < 4; ++i) {
            // Check
            deleteQtreeGPU2(qtree->quadrants[i]);
            // free(qtree->quadrants[i], q);
            // delete(qtree->quadrants[0]);
            // qtree->quadrants.clear();
            // (qtree->quadrants[i]).reset(NULL); // free memory for unique_ptr
        }
        free(qtree->quadrants[0], q);
        free(qtree->quadrants, q);
    }

    return;
}

void deleteQtreeGPU3(QtreeG qtree) // Este tengo que utilizarlo cuando creo los nodos contiguos
{
    if(isLeaf(qtree))
    {
      qtree->points.clear();
      free(qtree->quadrants, q);
    } 
    else 
    {
        // aqui no borro porque los nodos intermedios no van a tener vector
        for(int i = 0; i < 4; ++i) {
            // Check
            // printf("octante intermedio %d\n", i);
            // fflush(stdout);
            deleteQtreeGPU3(qtree->quadrants[i]);
            // free(qtree->quadrants[i], q);
            // delete(qtree->quadrants[0]);
            // qtree->quadrants.clear();
            // (qtree->quadrants[i]).reset(NULL); // free memory for unique_ptr
        }
        // free(qtree->quadrants[0], q);
        free(qtree->quadrants, q);
    }

    return;
}


std::vector<Qtree>::iterator findIterator2(Qtree aux, Qtree current, int& idx) 
{
  std::vector<Qtree>::iterator iterador = std::find(aux->quadrants.begin(), aux->quadrants.end(), current);

  int distance = std::distance(aux->quadrants.begin(), iterador);

  if(distance<3){

    std::advance(iterador, 1);
    idx = distance+1;

  } else { //termino

    iterador = aux->quadrants.begin();
    std::advance(iterador, 3);
    idx = 4;

  }

  // current = *itt;

  return iterador;
}

std::vector<Qtree>::iterator findIterator(Qtree aux, Qtree current, int& idx) 
{
  std::vector<Qtree>::iterator iterador = aux->quadrants.begin();

  if(aux->quadrants[0] == current){
    std::advance(iterador, 1);
    idx=1;
  }
  else if(aux->quadrants[1] == current){
    std::advance(iterador, 2);
    idx=2;
  }
  else if(aux->quadrants[2] == current){
    std::advance(iterador, 3);
    idx=3;
  }
  else {
    //lo pongo apuntando al último
    std::advance(iterador, 3);
    idx = 4; // termino, pero dejo listo current
  }

  return iterador;
}

Lpoint cpuSearchNeighborsMin(Vector2D* point, Qtree qtree, float radius, int* numNeighs)
{
    int idx = 0;
    int numInside = 0;
    Vector2D boxMin, boxMax;

    // *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    Lpoint min = {0,0.0,0.0,99999.0};

    Qtree aux = qtree;
    auto itt = qtree->quadrants.begin();
    Qtree current = *itt;
    // Qtree current = qtree->quadrants[0];

    // int nivel = 1;

    // auto iterador = std::find(std::begin(qtree->quadrants), std::end(qtree->quadrants), current);

    // if(iterador != std::end(qtree->quadrants))
    //   printf("\t------ENCONTRADO!-------\n");
    // else
    //   printf("\t------NO ENCONTRADO!-------\n");


    while(current != qtree) {

      // printf("ITER ");
      // fflush(stdout);

      if(idx > 3) {

        // printf(" termino en nodo %d, nivel %d; ", idx, nivel);
        // fflush(stdout);

        // nivel--;

        // printf(" PADRE... ");
        // fflush(stdout);

        current = current->parent;

        if(current != qtree){

          aux = current->parent;

          itt = findIterator(aux, current, idx);

          current = *itt;

          // printf("anterior idx: %d, NIVEL: %d\n", idx, nivel);
        }
        
      } else {

        if(isLeaf(current)) {

          // printf("\thoja %d\n", idx);
          // fflush(stdout);
          if(boxOverlap2D(boxMin, boxMax, current)) {

            if(boxInside2D(boxMin, boxMax, current)) {
              for(Lpoint* p : current->points) {
                if (p->z < min.z) {
                    min = *p;
                }
                numInside++;
              }
            } else {
              for(Lpoint* p : current->points) {
                if(insideBox2D(p, boxMin, boxMax))
                {
                  if (p->z < min.z) {
                      min = *p;
                  }
                  numInside++;
                }
              }
            }

          }

          idx++;
          if(idx < 4) {
            std::advance(itt, 1);
            current = *itt;
          }

        } else {

          // printf("solape? ");
          // fflush(stdout);
          
          if(!boxOverlap2D(boxMin, boxMax, current)) {
            // printf("nodo %d ", idx);
            idx++;
            // printf("no -> siguiente; idx: %d\n", idx);
            if(idx < 4) {
              std::advance(itt, 1);
              current = *itt;
            }
          }
          else {
            // printf("nodo %d ", idx);
            idx = 0;
            itt = current->quadrants.begin();
            current = *itt;
            // nivel++;
            // printf("sí -> sigo; NIVEL: %d\n", nivel);
          }
        }
      }

    }

    (*numNeighs) = numInside;

    return min;
}

void stage1s(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, std::vector<int>& minIDs, Qtree qtreeIn, Vector2D min)
{

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    Vector2D cellCenter;

    Lpoint newmin;

    int cellPoints;

    for(int jj = 0 ; jj < Ccol ; jj++ ){

        cellCenter.y = initY + jj*Displace;

        for(int ii = 0 ; ii < Crow ; ii++ ){

          // printf("CELDA %d,%d <-------------------------------------------------------\n",jj,ii);

            cellCenter.x = initX + ii*Displace;

            // newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);
            newmin = cpuSearchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);
            // newmin = gpuSearchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

            if(cellPoints >= minNumPoints ){
              // printf("minimo: %d", newmin.id);
              minIDs.push_back(newmin.id);
            }
        }

    }

    return;
}

std::vector<int> stage1reduce(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
                        unsigned short minNumPoints, Qtree qtreeIn, Vector2D min)
{
  double Displace = round2d(Wsize*(1-Overlap));

  double initX = min.x - Wsize/2 + Displace;
  double initY = min.y - Wsize/2 + Displace;

  std::vector<int> init;

  return tbb::parallel_reduce(
      tbb::blocked_range2d<int,int>{0,Ccol,0,Crow},
      init,
      [&](const tbb::blocked_range2d<int,int>& r, std::vector<int> v) {

        Lpoint newmin;
        Vector2D cellCenter;
        int cellPoints;
        int je = r.rows().end();
        int ie = r.cols().end();

        for(int jj = r.rows().begin(); jj < je; ++jj) {

            cellCenter.y = initY + jj*Displace;

            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

                cellCenter.x = initX + ii*Displace;

                newmin = searchNeighborsMin(&cellCenter, qtreeIn, Wsize/2, &cellPoints);

                if( cellPoints >= minNumPoints ){
                    v.push_back(newmin.id);
                }
            }
        }
        return v;
      },
      [](std::vector<int> a, std::vector<int> b) {
         a.insert(a.end(), b.begin(), b.end());
         return a;
      }
  );

}

unsigned int stage2(unsigned int countMin, std::vector<int>& minIDs){

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



void stage3s(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
          std::vector<int>& minGridIDs, Qtree qtreeIn, Qtree grid, Vector2D min){

    int cellPoints;

    Vector2D cellCenter;

    Lpoint newmin;

    double initX = min.x + Bsize/2;
    double initY = min.y + Bsize/2;

    for(int jj = 0 ; jj < Ccol ; jj++ ){

      cellCenter.y = initY + jj*Bsize;

      for(int ii = 0 ; ii < Crow ; ii++ ){

        cellCenter.x = initX + ii*Bsize;

        countNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);

        // Tengo que hacerlo porque el algoritmo no lo hace
        if(cellPoints == 0){

          newmin = searchNeighborsMin(&cellCenter, qtreeIn, Bsize/2, &cellPoints);

          if(cellPoints>0){
            minGridIDs.push_back(newmin.id);
          }
        }
      }
    }

    return;
}

std::vector<int> stage3reduce(unsigned short Bsize, unsigned short Crow, unsigned short Ccol,
                            Qtree qtreeIn, Qtree grid, Vector2D min)
{
  double initX =  min.x + Bsize/2;
  double initY = min.y + Bsize/2;

  std::vector<int> init;

  return tbb::parallel_reduce(
      tbb::blocked_range2d<int,int>{0,Ccol,0,Crow},
      init,
      [&](const tbb::blocked_range2d<int,int>& r, std::vector<int> v) {

        Lpoint newmin;
        Vector2D cellCenter;
        int cellPoints;
        int je = r.rows().end();
        int ie = r.cols().end();

        for(int jj = r.rows().begin(); jj < je; ++jj) {

            cellCenter.y = initY + jj*Bsize;

            for(int ii = r.cols().begin() ; ii < ie ; ii++ ){

                cellCenter.x = initX + ii*Bsize;

                countNeighbors2D(&cellCenter, grid, Bsize/2, &cellPoints);

                if( cellPoints == 0 ){

                    newmin = searchNeighborsMin(&cellCenter, qtreeIn, Bsize/2, &cellPoints);

                    if( cellPoints > 0 ){
                        v.push_back(newmin.id);
                    }
                }
            }
        }
        return v;
      },
      [](std::vector<int> a, std::vector<int> b) {
         a.insert(a.end(), b.begin(), b.end());
         return a;
      }
  );

}

int save_file(std::string file_name, std::vector<int>& ids, Lpoint** pointer)
{
  // unsigned int size = pointList.size();
  std::ofstream out;
  out.precision(std::numeric_limits<double>::digits10);
  out.open(file_name);
  // Point aux;
  if(!out.is_open())
    return -1;
  for( int i : ids ){
    out << (*pointer)[i].x << " " << (*pointer)[i].y << " " << (*pointer)[i].z << std::endl;
    // fprintf(out, "%.2f %.2f %.2f\n", (*pointer)[i].x, (*pointer)[i].y, (*pointer)[i].z);
  }
  out.close();
  if(out.is_open())
    return -1;
  return 0;
}

int check_results(std::string filename, std::vector<int>& ids, Lpoint** pointer, float Displace)
{
  std::string line;
  // char delim;
  Vector3D aux;
  unsigned int npoints = 0;
  std::vector<Vector3D> vgold, v;

  std::ifstream input(filename);
  if (input.fail()) {
      std::cout << "The GOLD point list doesn't exits" << std::endl;
      return 0;
  }

  while(getline(input,line)) {
    std::stringstream s(line);
    // s >> aux.x >> delim >> aux.y >> delim >> aux.z;
    s >> aux.x >> aux.y >> aux.z;
    vgold.push_back(aux);
    npoints++;
  }

  (input).close();
  if((input).is_open())
    return 0;

  std::sort(vgold.begin(), vgold.end(), [](const Vector3D a, const Vector3D b){
    return a.z < b.z;
  });

  // for(auto p1 : v){
  //   printf("%lf %lf %lf\n", p1.x, p1.y, p1.z);
  // }

  for( int i : ids ){
    aux.x = (*pointer)[i].x;
    aux.y = (*pointer)[i].y;
    aux.z = (*pointer)[i].z;
    v.push_back(aux);
  }

  std::sort(v.begin(), v.end(), [](const Vector3D a, const Vector3D b){
    return a.z < b.z;
  });

  int count = 0;
  // float Displace = 1.0;
  Vector2D boxMin, boxMax;
  for(auto p : v){
    Vector2D aux = {p.x, p.y};
    makeBox(&aux, Displace*0.5, &boxMin, &boxMax);
    for(auto pg : vgold){
      if(insideBox2D(&pg,boxMin,boxMax)){
        // if(fabs(p.z - pg.z) < 0.01){
          count++;
          break;
        }
    }
  }

  double rate = count/((double)(npoints))*100.0;
  printf("%d points correct; %.2f%%\n", count, rate);

  return 0;

  // return std::equal(v.begin(), v.end(), v2.begin(), [](const Vector3D a, const Vector3D b){
  //   return fabs(a.z - b.z) < 0.01;
  // });
}

uint64_t read_points(std::string filename, Lpoint** point_cloud) 
{
  std::ifstream input(filename);
  if (input.fail())
    return -1;

  std::string line;
  uint64_t id = 0;
  
  while(getline(input,line)) {
    std::stringstream s(line); 
    (*point_cloud)[id].id = id;
    s >> (*point_cloud)[id].x >> (*point_cloud)[id].y >> (*point_cloud)[id].z;
    id++;
  }

  input.close();
  if(input.is_open())
    return -1;

  return id;
}


int boxInside2D(Vector2D boxMin, Vector2D boxMax, QtreeG qt)
{
    if(qt->center.x + qt->radius > boxMax.x ||
       qt->center.y + qt->radius > boxMax.y)
        return 0;

    if(qt->center.x - qt->radius < boxMin.x ||
       qt->center.y - qt->radius < boxMin.y)
        return 0;

    return 1;
}

/** Wrapper to handle malloc() nicely */
void* mallocWrap(size_t size)
{
    void *ptr = malloc_host(size, q);
    if(!ptr)
    {
        fprintf(stderr, "Not enough memory\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

QtreeG createQtreeGPU(QtreeG parent, Vector2D center, float radius)
{
    QtreeG qt = static_cast<QtreeG>(mallocWrap(sizeof(QtreeG_t)));

    new (qt) QtreeG_t(p_alloc);

    qt->quadrants = static_cast<QtreeG*>(mallocWrap(4 * sizeof(QtreeG)));

    qt->parent = parent;
    qt->center = center;
    qt->radius = radius;

    // usm_allocator<int, usm::alloc::shared> p_alloc{q};
    // qt->points = std::vector<Lpoint*, usm_allocator<Lpoint*, usm::alloc::shared>> vector_prueba(p_alloc);
    // qt->points.reserve(100);
    // qt->numPts = 0;
    // qt->plane[0] = 0.0;
    // qt->plane[1] = 0.0;
    // qt->plane[2] = 0.0;
    // qt->plane[3] = 0.0;
    for(int i = 0; i < 4; i++){
        qt->quadrants[i] = NULL;
    }

    return qt;
}

QtreeG createQtreeGPU2(QtreeG parent, Vector2D center, float radius)
{
    QtreeG_t newNode(parent, center, radius, p_alloc);

    // newNode.quadrants = static_cast<QtreeG*>(mallocWrap(4 * sizeof(QtreeG)));

    // for(int i = 0; i < 4; i++){
    //     newNode.quadrants[i] = NULL;
    // }

    all_nodes.push_back(newNode);

    return &(all_nodes.back());
}

void createQuadrantsGPU(QtreeG qt)
{
    Vector2D newCenter;
    float newRadius = qt->radius * 0.5;

    for(int i = 0; i < 4; i++)
    {
        newCenter = qt->center;
        newCenter.x += qt->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += qt->radius * (i&1 ? 0.5f : -0.5f);
        // newCenter.z += qt->radius * (i&1 ? 0.5f : -0.5f);
        // printf("(%lf,%lf) ", newCenter.x, newCenter.y);

        // qt->quadrants[i] = new Qtree_t(qt, newCenter, newRadius);
        qt->quadrants[i] = createQtreeGPU(qt, newCenter, newRadius);
    }
    // printf("\n");
}

void createQuadrantsGPU2(QtreeG qt)
{
    Vector2D newCenter;
    float newRadius = qt->radius * 0.5;

    QtreeG quadrants = static_cast<QtreeG>(mallocWrap(4 * sizeof(QtreeG_t)));

    for(int i = 0; i < 4; i++)
    {
        newCenter = qt->center;
        newCenter.x += qt->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += qt->radius * (i&1 ? 0.5f : -0.5f);
        // newCenter.z += qt->radius * (i&1 ? 0.5f : -0.5f);
        // printf("(%lf,%lf) ", newCenter.x, newCenter.y);

        new (quadrants + i) QtreeG_t(qt, newCenter, newRadius, p_alloc);

        // qt->quadrants[i] = new Qtree_t(qt, newCenter, newRadius);
        // qt->quadrants[i] = createQtreeGPU(qt, newCenter, newRadius);

        qt->quadrants[i] = quadrants + i;
    }
    // printf("\n");
}

void createQuadrantsGPU3(QtreeG qt)
{
    Vector2D newCenter;
    float newRadius = qt->radius * 0.5;

    for(int i = 0; i < 4; i++)
    {
        // printf("Parent pos: %p; ", qt);
        newCenter = qt->center;
        newCenter.x += qt->radius * (i&2 ? 0.5f : -0.5f);
        newCenter.y += qt->radius * (i&1 ? 0.5f : -0.5f);
        // newCenter.z += qt->radius * (i&1 ? 0.5f : -0.5f);
        // printf("(%lf,%lf) ", newCenter.x, newCenter.y);

        QtreeG_t newNode(qt, newCenter, newRadius, p_alloc);

        // if(qt->quadrants[i] != NULL)
        //   printf("ERROR, %d esta lleno! ", i);

        // printf(" %d",i);
        // fflush(stdout);
        all_nodes.push_back(newNode);
        // printf(".");
        // fflush(stdout);
        qt->quadrants[i] = &(all_nodes.back());
        // printf(".");
        // fflush(stdout);
    }
    // printf("\n");
}

// void createQuadrantsGPU4(QtreeG qt)
// {
//     Vector2D newCenter;
//     float newRadius = qt->radius * 0.5;

//     for(int i = 0; i < 4; i++)
//     {
//         // printf("Parent pos: %p; ", qt);
//         newCenter = qt->center;
//         newCenter.x += qt->radius * (i&2 ? 0.5f : -0.5f);
//         newCenter.y += qt->radius * (i&1 ? 0.5f : -0.5f);

//         QtreeG_t newNode(qt, newCenter, newRadius, p_alloc);

//         all_nodes.push_back(newNode);

//         all_addresses.push_back(&(all_nodes.back()));

//         qt->quadrants[i] = &(all_nodes.back());

//     }
//     // printf("\n");
// }

// template<typename nT>
// void insertPointGPU(Lpoint *point, nT qtree, float minRadius)
// {
//     int idx = 0;

//     if(isLeaf(qtree))
//     {
//         // printf("  octante hoja\n");
//         // fflush(stdout);
//         if(qtree->radius * 0.5 > minRadius)    // still divisible -> divide
//         {
//           // printf("    divido; ");
//           // fflush(stdout);
//           // createQuadrantsGPU(qtree);
//           // createQuadrantsGPU2(qtree);
//           createQuadrantsGPU3(qtree);
//           // createQuadrantsGPU4(qtree);

//           // printf("; Ahora hay %zu nodos\n", all_nodes.size());
//           // fflush(stdout);

//           // fillOctants(qtree);
//           idx = quadrantIdx(point, qtree);
//           insertPointGPU(point, qtree->quadrants[idx], minRadius);

//         } else {
//           // printf("    inserto\n");
//           // fflush(stdout);
//           qtree->points.push_back(point);
//         }
//     }
//     else                                // No leaf -> search the correct one
//     {
//       idx = quadrantIdx(point, qtree);
//       // printf("octante intermedio %d\n", idx);
//       // fflush(stdout);
//       insertPointGPU(point, qtree->quadrants[idx], minRadius);
//     }
// }

QtreeG findPosition(QtreeG parent, QtreeG current, int& idx) 
{

  if(parent->quadrants[0] == current){
    // std::advance(iterador, 1);
    idx=1;
    return parent->quadrants[1];
  }
  else if(parent->quadrants[1] == current){
    // std::advance(iterador, 2);
    idx=2;
    return parent->quadrants[2];
  }
  else if(parent->quadrants[2] == current){
    // std::advance(iterador, 3);
    idx=3;
    return parent->quadrants[3];
  } 
  else {
    //lo pongo apuntando al último
    // std::advance(iterador, 3);
    idx = 4; // termino, pero dejo listo current
    return parent->quadrants[3];
  }

}


Lpoint gpuSearchNeighborsMin(Vector2D* point, QtreeG qtree, float radius, int* numNeighs)
{
    int idx = 0;
    int numInside = 0;
    Vector2D boxMin, boxMax;

    // *numNeighs = 0;
    makeBox(point, radius, &boxMin, &boxMax);

    Lpoint min = {0,0.0,0.0,99999.0};

    QtreeG parent = qtree;
    QtreeG current = qtree->quadrants[0];
    // QtreeG current = qtree->quadrants[0];

    // int nivel = 1;

    // auto iterador = std::find(std::begin(qtree->quadrants), std::end(qtree->quadrants), current);

    // if(iterador != std::end(qtree->quadrants))
    //   printf("\t------ENCONTRADO!-------\n");
    // else
    //   printf("\t------NO ENCONTRADO!-------\n");


    while(current != qtree) {

      // printf("ITER ");
      // fflush(stdout);

      if(idx > 3) {

        // printf(" termino en nodo %d, nivel %d; ", idx, nivel);
        // fflush(stdout);

        // nivel--;

        // printf(" PADRE... ");
        // fflush(stdout);

        current = current->parent;

        if(current != qtree){

          parent = current->parent;

          current = findPosition(parent, current, idx);

          // current = *itt;

          // printf("anterior idx: %d, NIVEL: %d\n", idx, nivel);
        }
        
      } else {

        if(isLeaf(current)) {

          // printf("\thoja %d\n", idx);
          // fflush(stdout);
          if(boxOverlap2D(boxMin, boxMax, current)) {

            if(boxInside2D(boxMin, boxMax, current)) {
              for(Lpoint* p : current->points) {
                if (p->z < min.z) {
                    min = *p;
                }
                numInside++;
              }
            } else {
              for(Lpoint* p : current->points) {
                if(insideBox2D(p, boxMin, boxMax))
                {
                  if (p->z < min.z) {
                      min = *p;
                  }
                  numInside++;
                }
              }
            }

          }

          idx++;
          if(idx < 4) {
            current = parent->quadrants[idx];
            // std::advance(itt, 1);
            // current = *itt;
          }

        } else {

          // printf("solape? ");
          // fflush(stdout);
          
          if(!boxOverlap2D(boxMin, boxMax, current)) {
            // printf("nodo %d ", idx);
            idx++;
            // printf("no -> siguiente; idx: %d\n", idx);
            if(idx < 4) {
              current = parent->quadrants[idx];
              // std::advance(itt, 1);
              // current = *itt;
            }
          }
          else {
            // printf("nodo %d ", idx);
            idx = 0;
            parent = current;
            current = current->quadrants[0];
            // itt = current->quadrants.begin();
            // current = *itt;
            // nivel++;
            // printf("sí -> sigo; NIVEL: %d\n", nivel);
          }
        }
      }

    }

    (*numNeighs) = numInside;

    return min;
}


void stage1gpu(unsigned short Wsize, double Overlap, unsigned short Crow, unsigned short Ccol,
  unsigned short minNumPoints, int* minGPU, QtreeG qtree, Vector2D min)
{
    // usm_allocator<int, usm::alloc::shared> new_alloc{q};

    // tbb::concurrent_vector<int, usm_allocator<int, usm::alloc::host>> v(new_alloc);
    // std::vector<int, usm_allocator<int, usm::alloc::shared>> v(new_alloc);

    // v.reserve(10);

    std::cout << " Using device: " << q.get_device().get_info<info::device::name>() << std::endl;

    double Displace = round2d(Wsize*(1-Overlap));

    double initX = min.x - Wsize/2 + Displace;
    double initY = min.y - Wsize/2 + Displace;

    // std::vector<int> values(Ccol);
    // std::iota(std::begin(values), std::end(values), 0);

    // q.parallel_for(range<1>(Ccol), [=](id<1> j) { 

    //   int jj = static_cast<int>(j);

    //   for(int ii = 0. ; ii < Crow ; ii++ ){  

    //     Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

    //     int cellPoints;

    //     Lpoint newmin = gpuSearchNeighborsMin(&cellCenter, qtree, Wsize/2, &cellPoints);


    //     if( cellPoints >= minNumPoints ){
    //       minGPU[jj*Crow+ii] = newmin.id;
    //     }

    //   }

    // }).wait();

    // std::vector<int> prueba;

    // prueba.reserve(100);

    // buffer<int, 1> bufout_vect(prueba.data(), range<1>(100));

    q.submit([&](handler &h) {
      //# setup sycl stream class to print standard output from device code
      auto out = stream(4096, 1024, h);

      // auto V3 = bufout_vect.get_access<access::mode::write>(h);

      //# nd-range kernel
      h.parallel_for(range<1>(Ccol), [=](id<1> j) { 

        int jj = static_cast<int>(j);

        for(int ii = 0. ; ii < Crow ; ii++ ){  

          Vector2D cellCenter = {initX + ii*Displace, initY + jj*Displace};

          int cellPoints;

          Lpoint newmin = gpuSearchNeighborsMin(&cellCenter, qtree, Wsize/2, &cellPoints);

          // out << " (" << Ccol << ", " << Crow << ") " <<  "minimo: (" << ii << ", " << jj << ") " << endl;

          if( cellPoints >= minNumPoints ){

              // out << "minimo: (" << ii << ", " << jj << ") " << endl;
              // out << "minimo: " << newmin.id << endl;
              minGPU[jj*Crow+ii] = newmin.id;

              // v.push_back(newmin.id);
              // printf("Tengo un minimo!\n");
              // std::cout << "Tengo un minimo\n";

          }
        }

      });

    }).wait();


    return;
}


unsigned int stage2GPU(unsigned int countMin, int* minIDs){

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
