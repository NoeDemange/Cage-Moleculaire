#ifndef __PATHFINDING_H
#define __PATHFINDING_H

#include "structure.h"
#include "voxelization.h"

float dijkstra(Point3D start, Point3D goal, VOXELGRID voxelGrid);
float aStarPathfinding(Point3D start, Point3D goal, VOXELGRID voxelGrid);

#endif