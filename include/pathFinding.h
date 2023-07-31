#ifndef __PATHFINDING_H
#define __PATHFINDING_H

#include "structure.h"
#include "voxelization.h"

/**
 * @file pathfinding.h
 * @brief Pathfinding Header File
 *
 * This file contains function for pathfinding algorithms on the voxel grid.
 */

/**
 * @brief Perform Dijkstra's algorithm on the voxel grid to find the shortest path.
 *
 * @param start The starting point.
 * @param goal The goal point.
 * @param voxelGrid The 3D voxel grid.
 * @param vMap The voxel map for storing distances and index heap information.
 * @param nodeHeap The node heap for storing nodes during the algorithm.
 * @return The distance of the shortest path from the start to the goal, or -1 if no path is found.
 */
float dijkstra(Point3D, Point3D, VOXELGRID, VMap***, NodeHeap);

/**
 * @brief Perform A* pathfinding on the voxel grid to find the shortest path.
 *
 * @param start The starting point.
 * @param goal The goal point.
 * @param voxelGrid The 3D voxel grid.
 * @param vMap The voxel map for storing distances and index heap information.
 * @param nodeHeap The node heap for storing nodes during the algorithm.
 * @return The distance of the shortest path from the start to the goal, or -1 if no path is found.
 */
float aStarPathfinding(Point3D, Point3D, VOXELGRID, VMap***, NodeHeap);

/**
 * @brief Compute the distance between two points considering obstacles.
 *
 * This function calculates the distance between the start and end points, considering obstacles
 * in the voxel grid using A* pathfinding algorithm. It uses voxel distance heuristics to compute
 * the path distance.
 *
 * @param startPos The starting position.
 * @param endPos The ending position.
 * @param voxelGrid The 3D voxel grid.
 * @param vMap The voxel map for storing distances and index heap information.
 * @param nodeHeap The node heap for storing nodes during the algorithm.
 * @return The computed distance between the start and end points, considering obstacles.
 */
float distWithObstacles(Point_t startPos, Point_t endPos, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap);

#endif