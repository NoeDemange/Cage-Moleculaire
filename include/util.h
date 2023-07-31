#ifndef __UTIL_H
#define __UTIL_H

#include "structure.h"

float radianToDegre(float);
float degreToRadian(float);
Point_t initPoint(float);
Point_t addPoint(Point_t, Point_t);
Point_t subPoint(Point_t, Point_t);
Point_t mulPoint(Point_t, float);
Point_t divPoint(Point_t, float);
Point_t merPoint(Point_t, Point_t);
int EqualPoint(Point_t A, Point_t B);
float dist(Point_t, Point_t);
int distInf(Point_t A, Point_t B, float dist);
float dist_manhattan(Point_t, Point_t);
float dist_cercle(Point_t, Point_t);
float dist_obstacle(Point_t A, Point_t B, Molecule_t* sub);
Point_t normalization(Point_t, float);
float angle(Point_t, Point_t, Point_t);
Point_t vector(Point_t, Point_t);
Point_t addThirdPoint(Point_t, Point_t, Point_t, float);
Point_t planNormal(Point_t, Point_t, Point_t);
Point_t rotation(Point_t, float, Point_t);

Point_t AX1E1(Point_t, Point_t, float);
Point_t AX2E1(Point_t, Point_t, Point_t, float);
Point_t AX1E2(Point_t, Point_t, Point_t, float);
Point_t AX3E1(Point_t, Point_t, Point_t, Point_t, float);
Point_t AX2E2(Point_t, Point_t, Point_t, float);
Point_t AX1E3(Point_t, Point_t, Point_t, float);

#endif
