#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "AStar.h"
#include "constant.h"
#include "util.h"

// Helper function to calculate the Manhattan distance between two points
float manhattanDistance(Point3D a, Point3D b) {
    return abs((START_GRID_X + a.x * LENGTH_GRID_X) - (START_GRID_X + b.x * LENGTH_GRID_X)) 
        + abs((START_GRID_Y + a.y * LENGTH_GRID_Y) - (START_GRID_Y + b.y * LENGTH_GRID_Y)) 
        + abs((START_GRID_Z + a.z * LENGTH_GRID_Z) - (START_GRID_Z + b.z * LENGTH_GRID_Z));
}

// Helper function to calculate the distance between two points,  src: The Jump Point Search Pathfinding System in 3D
/*The three-dimensional distance between a node n and node n′ on a 26-connected grid is composed of distances ∆x,∆y and ∆z 
in each respective axis. Let dmax=max(∆x,∆y,∆z), dmin=min(∆x,∆y,∆z), and dmid={∆x,∆y,∆z}\{dmax,dmin}.
  The voxel distance between nodes n and n′ is subsequently calculated as:h(n,n′) = (√3−√2)dmin+ (√2−1)dmid+dmax
*/
float dist26C(Point3D a, Point3D b) {
    float deltaX = abs((START_GRID_X + a.x * LENGTH_GRID_X) - (START_GRID_X + b.x * LENGTH_GRID_X));
    float deltaY = abs((START_GRID_Y + a.y * LENGTH_GRID_Y) - (START_GRID_Y + b.y * LENGTH_GRID_Y));
    float deltaZ = abs((START_GRID_Z + a.z * LENGTH_GRID_Z) - (START_GRID_Z + b.z * LENGTH_GRID_Z));
    float dmax = (deltaX > deltaY) ? ((deltaX > deltaZ) ? deltaX : deltaZ) : ((deltaY > deltaZ) ? deltaY : deltaZ);
    float dmin = (deltaX < deltaY) ? ((deltaX < deltaZ) ? deltaX : deltaZ) : ((deltaY < deltaZ) ? deltaY : deltaZ);
    float dmid = deltaX + deltaY + deltaZ - dmax - dmin;
    return ((sqrt(3) - sqrt(2)) * dmin + (sqrt(2) - 1) * dmid + dmax);
}

// Helper function to check if a point is within the grid bounds
int isPointValid(Point3D point) {
    return (point.x >= 0 && point.x < GRID_SIZE_X && point.y >= 0 && point.y < GRID_SIZE_Y && point.z >= 0 && point.z < GRID_SIZE_Z);
}

// Helper function to check if a point is traversable (not an obstacle) in the voxel grid
int isPointTraversable(Point3D point, VOXELGRID voxelGrid) {
    if(voxelGrid[point.x][point.y][point.z] == 0 ) return 1;
    return 0;
}

// Helper function to check if a point is the goal point
int isPointGoal(Point3D point, Point3D goal) {
    return (point.x == goal.x && point.y == goal.y && point.z == goal.z);
}

void freeGScores(float *** gScores){
    // Free dynamically allocated memory
    for (int x = 0; x < GRID_SIZE_X; x++) {
        for (int y = 0; y < GRID_SIZE_Y; y++) {
            free(gScores[x][y]);
        }
        free(gScores[x]);
    }
    free(gScores);
}

void freeCloseSet(int*** closedSet){
    // Free the dynamically allocated memory
    for (int x = 0; x < GRID_SIZE_X; x++) {
        for (int y = 0; y < GRID_SIZE_Y; y++) {
            free(closedSet[x][y]);
        }
        free(closedSet[x]);
    }
    free(closedSet);
}

// A* algorithm for pathfinding in voxel grid
float aStarPathfinding(Point3D start, Point3D goal, VOXELGRID voxelGrid) {
    float*** gScores = (float***)malloc(GRID_SIZE_X * sizeof(float**));
    for (int x = 0; x < GRID_SIZE_X; x++) {
        gScores[x] = (float**)malloc(GRID_SIZE_Y * sizeof(float*));
        for (int y = 0; y < GRID_SIZE_Y; y++) {
            gScores[x][y] = (float*)malloc(GRID_SIZE_Z * sizeof(float));
            for (int z = 0; z < GRID_SIZE_Z; z++) {
                gScores[x][y][z] = __FLT_MAX__;
            }
        }
    }
    Node startNode = createNode(start, 0.0, dist26C(start, goal));
    //printf("x : %d, y : %d, z : %d\n",start.x,start.y, start.z);
    gScores[start.x][start.y][start.z] = 0.0;

    Node* openSet = (Node*)malloc(GRID_SIZE_X * GRID_SIZE_Y * GRID_SIZE_Z * sizeof(Node));
    int openSetSize = 0;
    openSet[openSetSize++] = startNode;

    int*** closedSet = (int***)calloc(GRID_SIZE_X, sizeof(int**));
    for (int x = 0; x < GRID_SIZE_X; x++) {
        closedSet[x] = (int**)calloc(GRID_SIZE_Y, sizeof(int*));
        for (int y = 0; y < GRID_SIZE_Y; y++) {
            closedSet[x][y] = (int*)calloc(GRID_SIZE_Z, sizeof(int));
        }
    }

    while (openSetSize > 0) {
        // Find the node with the lowest f value
        int currentIndex = 0;
        for (int i = 1; i < openSetSize; i++) {
            if (openSet[i].f < openSet[currentIndex].f) {
                currentIndex = i;
            }
        }
        Node current = openSet[currentIndex];
        Point3D currentPoint = current.point;

        if (isPointGoal(currentPoint, goal)) {
            freeGScores(gScores);
            free(openSet);
            freeCloseSet(closedSet);
            return current.g; // Return the distance of the minimal path found
        }

        // Remove the current node from the open set
        openSet[currentIndex] = openSet[openSetSize - 1];
        openSetSize--;

        closedSet[currentPoint.x][currentPoint.y][currentPoint.z] = 1;

        // Generate the neighbors
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (dx == 0 && dy == 0 && dz == 0) {
                        continue;
                    }
                    Point3D neighbor = { currentPoint.x + dx, currentPoint.y + dy, currentPoint.z + dz };
                    if (!isPointValid(neighbor) || closedSet[neighbor.x][neighbor.y][neighbor.z] || !isPointTraversable(neighbor, voxelGrid)) {
                        continue;
                    }
                    float tentativeGScore = current.g + dist26C(currentPoint,neighbor);/*1;*///A changer 
                    if (tentativeGScore < gScores[neighbor.x][neighbor.y][neighbor.z]) {
                        gScores[neighbor.x][neighbor.y][neighbor.z] = tentativeGScore;
                        float hScore = dist26C(neighbor, goal);
                        Node neighborNode = createNode(neighbor, tentativeGScore, hScore);
                        openSet[openSetSize++] = neighborNode;
                    }
                }
            }
        }
    }
    freeGScores(gScores);
    free(openSet);
    freeCloseSet(closedSet);
    return -1; // Return -1 if path not found
}