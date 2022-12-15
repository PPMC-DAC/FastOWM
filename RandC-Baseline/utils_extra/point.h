#ifndef POINT_H

#define	POINT_H

typedef struct 
{
    unsigned int id;
    unsigned int gId;   
    double x;
    double y;
    double z;
    unsigned int intensity;
    float normal[3];

} Lpoint; // LiDAR point data structure

Lpoint createPoint(unsigned int id, double x, double y, double z, unsigned int intensity, double normal[3]);

void copyCoords(Lpoint *dest, Lpoint *source);

Lpoint midPoint(Lpoint *p1, Lpoint *p2);

Lpoint avgPoint(Lpoint *points, int num);

Lpoint avgPointP(Lpoint **points, int num);

Lpoint avgPointInt(Lpoint **points, int num);

Lpoint* mergePoints(Lpoint *pts1, int num1, Lpoint *pts2, int num2);

Lpoint** distinctPoints(Lpoint **pts1, int numPts1, Lpoint **pts2, int numPts2, int *distNum);

void invertPoints(Lpoint *points, int num);

void invertPointsP(Lpoint **points, int num);

void rotatePoint(Lpoint *p, double sin, double cos, Lpoint *orig);

double pointLineDistance(Lpoint *p, Lpoint* lp1, Lpoint* lp2);

void linePoints(Lpoint *p, double *lineCoefs, Lpoint *p1, Lpoint *p2);

int lineIntersection(Lpoint *p1, Lpoint *p2, Lpoint *p3, Lpoint *p4, Lpoint *intersect);

Lpoint* curvePoints(Lpoint *dataPts, int numData, double *coefs, int *numCurve); 

Lpoint* curvePoints3(Lpoint *dataPts, int numData, double *coefs, int *numCurve);

double* lineCoeficients(Lpoint *p1, Lpoint *p2);

double boundingCubeWidth(Lpoint *points, int numPoints, Lpoint center);

#endif	/* POINT_H */