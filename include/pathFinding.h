#ifndef __PATHFINDING_H
#define __PATHFINDING_H

#include "structure.h"
#include "voxelization.h"

float dijkstra(Point3D, Point3D, VOXELGRID, VMap***, NodeHeap);
float aStarPathfinding(Point3D, Point3D, VOXELGRID, VMap***, NodeHeap);
float distWithObstacles(Point_t startPos, Point_t endPos, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap);

#endif