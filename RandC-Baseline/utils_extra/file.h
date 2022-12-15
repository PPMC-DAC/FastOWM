#ifndef FILE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vector.h"

#define	FILE_H


int readPoint(char *string, Lpoint *point);

int writePointsWithClass(Lpoint *points, int numPoints, char *fileName);

#endif  /* FILE_H */

