#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pathFinding.h"

/**************************************/
/* util  ******************************/
/**************************************/

int isPointGoal(Point3D point, Point3D goal) {
    return (point.x == goal.x && point.y == goal.y && point.z == goal.z);
}

int isPointValid(Point3D point) {
    return (point.x >= 0 && point.x < GRID_SIZE && point.y >= 0 && point.y < GRID_SIZE && point.z >= 0 && point.z < GRID_SIZE);
}

int isPointTraversable(Point3D point, VOXELGRID voxelGrid) {
    if(voxelGrid[point.x][point.y][point.z] == 0 ) return 1;
    return 0;
}

void freeVMap(VMap *** vMap){
    // Free dynamically allocated memory
    for (int x = 0; x < GRID_SIZE; x++) {
        for (int y = 0; y < GRID_SIZE; y++) {
            free(vMap[x][y]);
        }
        free(vMap[x]);
    }
    free(vMap);
}

/**************************************/
/* heuristics and cost  ***************/
/**************************************/

float cost(Point3D a, Point3D b) {
    int deltaX = abs(a.x-b.x);
    int deltaY = abs(a.y-b.y);
    int deltaZ = abs(a.z-b.z);
    int sum = deltaX + deltaY + deltaZ;
    if(sum == 3) return sqrt(3)*LENGTH_GRID;
    if(sum == 2) return sqrt(2)*LENGTH_GRID;
    return LENGTH_GRID;
}

float manhattanDistance(Point3D a, Point3D b) {
    return fabs((START_GRID + a.x * LENGTH_GRID) - (START_GRID + b.x * LENGTH_GRID)) 
        + fabs((START_GRID + a.y * LENGTH_GRID) - (START_GRID + b.y * LENGTH_GRID)) 
        + fabs((START_GRID + a.z * LENGTH_GRID) - (START_GRID + b.z * LENGTH_GRID));
}

//src: The Jump Point Search Pathfinding System in 3D
/*The three-dimensional distance between a node n and node n′ on a 26-connected grid is composed of distances ∆x,∆y and ∆z 
in each respective axis. Let dmax=max(∆x,∆y,∆z), dmin=min(∆x,∆y,∆z), and dmid={∆x,∆y,∆z}\{dmax,dmin}.
  The voxel distance between nodes n and n′ is subsequently calculated as:h(n,n′) = (√3−√2)dmin+ (√2−1)dmid+dmax
*/
float voxelDist(Point3D a, Point3D b) {
    float deltaX = fabs((START_GRID + a.x * LENGTH_GRID) - (START_GRID + b.x * LENGTH_GRID));
    float deltaY = fabs((START_GRID + a.y * LENGTH_GRID) - (START_GRID + b.y * LENGTH_GRID));
    float deltaZ = fabs((START_GRID + a.z * LENGTH_GRID) - (START_GRID + b.z * LENGTH_GRID));
    float dmax = (deltaX > deltaY) ? ((deltaX > deltaZ) ? deltaX : deltaZ) : ((deltaY > deltaZ) ? deltaY : deltaZ);
    float dmin = (deltaX < deltaY) ? ((deltaX < deltaZ) ? deltaX : deltaZ) : ((deltaY < deltaZ) ? deltaY : deltaZ);
    float dmid = deltaX + deltaY + deltaZ - dmax - dmin;
    return ((sqrt(3) - sqrt(2)) * dmin + (sqrt(2) - 1) * dmid + dmax);
}


/**************************************/
/* Algorithm **************************/
/**************************************/

float dijkstra(Point3D start, Point3D goal, VOXELGRID voxelGrid){
    NodeHeap nodeHeap = NH_initAlloc(GRID_SIZE*GRID_SIZE*GRID_SIZE);
    Node node;
    VMap*** vMap = (VMap***)malloc(GRID_SIZE * sizeof(VMap**));
    for (int x = 0; x < GRID_SIZE; x++) {
        vMap[x] = (VMap**)malloc(GRID_SIZE * sizeof(VMap*));
        for (int y = 0; y < GRID_SIZE; y++) {
            vMap[x][y] = (VMap*)malloc(GRID_SIZE * sizeof(VMap));
            for (int z = 0; z < GRID_SIZE; z++) {
                    vMap[x][y][z].dist = __FLT_MAX__;
                    vMap[x][y][z].indexHeap = -__INT_MAX__;
            }
        }
    }

    vMap[start.x][start.y][start.z].dist = 0.0;
    node = createNode(start, 0.0, 0.0);
    NH_insert(&nodeHeap,node,vMap);
    #ifdef DEBUGDij 
        int nbD = 0;
    #endif
    while(nodeHeap.size > 0){
        Node minN = NH_extractMin(&nodeHeap,vMap);
        #ifdef DEBUGDij 
            nbD++;
        #endif
        if(isPointGoal(minN.point,goal)){ // Goal achieved
            #ifdef DEBUGDij
                printf("Dijkstra nombre sommet vu : %d\n",nbD);
            #endif
            freeVMap(vMap);
            NH_free(nodeHeap);
            return minN.g;
        }

        // Find the neighbors
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (dx == 0 && dy == 0 && dz == 0) {
                        continue;
                    }
                    Point3D neighbor = { minN.point.x + dx, minN.point.y + dy, minN.point.z + dz };
                    if (!isPointValid(neighbor) || vMap[minN.point.x + dx][minN.point.y + dy][minN.point.z + dz].indexHeap==-1 || !isPointTraversable(neighbor, voxelGrid)) {
                        continue;
                    }
                    float tentativeGScore = minN.g + cost(minN.point,neighbor);
                    if (tentativeGScore < vMap[neighbor.x][neighbor.y][neighbor.z].dist) {
                        vMap[neighbor.x][neighbor.y][neighbor.z].dist = tentativeGScore;
                        if(vMap[neighbor.x][neighbor.y][neighbor.z].indexHeap == -__INT_MAX__){
                            Node neighborNode = createNode(neighbor, tentativeGScore, 0.0);
                            NH_insert(&nodeHeap,neighborNode,vMap);
                        } else{
                            NH_decrease_priority(&nodeHeap, vMap, neighbor, tentativeGScore);
                        }  
                    }
                }
            }
        }
    }

    freeVMap(vMap);
    NH_free(nodeHeap);
    return -1; // Return -1 if path not found
}

float aStarPathfinding(Point3D start, Point3D goal, VOXELGRID voxelGrid){
    NodeHeap nodeHeap = NH_initAlloc(GRID_SIZE*GRID_SIZE*GRID_SIZE);
    Node node;
    VMap*** vMap = (VMap***)malloc(GRID_SIZE * sizeof(VMap**));
    for (int x = 0; x < GRID_SIZE; x++) {
        vMap[x] = (VMap**)malloc(GRID_SIZE * sizeof(VMap*));
        for (int y = 0; y < GRID_SIZE; y++) {
            vMap[x][y] = (VMap*)malloc(GRID_SIZE * sizeof(VMap));
            for (int z = 0; z < GRID_SIZE; z++) {
                    vMap[x][y][z].dist = __FLT_MAX__;
                    vMap[x][y][z].indexHeap = -__INT_MAX__;
            }
        }
    }

    vMap[start.x][start.y][start.z].dist = 0.0;
    node = createNode(start, 0.0, voxelDist(start,goal));
    NH_insert(&nodeHeap,node,vMap);
    #ifdef DEBUGAstar 
        int nb = 0;
    #endif
    while(nodeHeap.size > 0){
        Node minN = NH_extractMin(&nodeHeap,vMap);
        #ifdef DEBUGAstar
            nb++;
        #endif
        if(isPointGoal(minN.point,goal)){ // Goal achieved
            #ifdef DEBUGAstar
                printf("A* nombre sommet vu : %d\n",nb);
            #endif
            freeVMap(vMap);
            NH_free(nodeHeap);
            return minN.g;
        }

        // Find the neighbors
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (dx == 0 && dy == 0 && dz == 0) {
                        continue;
                    }
                    Point3D neighbor = { minN.point.x + dx, minN.point.y + dy, minN.point.z + dz };
                    if (!isPointValid(neighbor) || vMap[minN.point.x + dx][minN.point.y + dy][minN.point.z + dz].indexHeap==-1 || !isPointTraversable(neighbor, voxelGrid)) {
                        continue;
                    }
                    float tentativeGScore = minN.g + cost(minN.point,neighbor);
                    if (tentativeGScore < vMap[neighbor.x][neighbor.y][neighbor.z].dist) {
                        vMap[neighbor.x][neighbor.y][neighbor.z].dist = tentativeGScore;
                        if(vMap[neighbor.x][neighbor.y][neighbor.z].indexHeap == -__INT_MAX__){
                            Node neighborNode = createNode(neighbor, tentativeGScore, voxelDist(neighbor,goal));
                            NH_insert(&nodeHeap,neighborNode,vMap);
                        } else{
                            NH_decrease_priority(&nodeHeap, vMap, neighbor, tentativeGScore);
                        }  
                    }
                }
            }
        }
    }

    freeVMap(vMap);
    NH_free(nodeHeap);
    return -1; // Return -1 if path not found
}