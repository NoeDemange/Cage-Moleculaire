#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pathFinding.h"
#include "util.h"

/**************************************/
/* heuristics  ************************/
/**************************************/

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


//peu performant
float dijkstra(Point3D start, Point3D goal, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap){
    nodeHeap.size = 0;
    Node node;
    for (int x = 0; x < GRID_SIZE; x++) {
        for (int y = 0; y < GRID_SIZE; y++) {
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
        if(minN.point.x == goal.x && minN.point.y == goal.y && minN.point.z == goal.z){ // Goal achieved
            #ifdef DEBUGDij
                printf("Dijkstra nombre sommet vu : %d\n",nbD);
            #endif
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
                    if (!(neighbor.x >= 0 && neighbor.x < GRID_SIZE && neighbor.y >= 0 && neighbor.y < GRID_SIZE && neighbor.z >= 0 && neighbor.z < GRID_SIZE) 
                    || vMap[minN.point.x + dx][minN.point.y + dy][minN.point.z + dz].indexHeap==-1 
                    || voxelGrid[neighbor.x][neighbor.y][neighbor.z] == 1) {
                        continue;
                    }
                    float tentativeGScore = minN.g + sqrt(abs(minN.point.x-neighbor.x)+abs(minN.point.y-neighbor.y)+abs(minN.point.z-neighbor.z))*LENGTH_GRID;
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
    return -1; // Return -1 if path not found
}

float aStarPathfinding(Point3D start, Point3D goal, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap){ //soit on déclare vMap et nodeheap dans la fonction mais trop de malloc sinon on déclare une seule fois mais on ne peut plus paralléliser facilement
    nodeHeap.size = 0;
    Node node;
    for (int x = 0; x < GRID_SIZE; x++) {
        for (int y = 0; y < GRID_SIZE; y++) {
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
        //static int call = 0;
        //call++;
        int nb = 0;
    #endif
    while(nodeHeap.size > 0){
        Node minN = NH_extractMin(&nodeHeap,vMap);
        #ifdef DEBUGAstar
            nb++;
        #endif
        if(minN.point.x == goal.x && minN.point.y == goal.y && minN.point.z == goal.z){ // Goal achieved
            #ifdef DEBUGAstar
                //printf("A* nombre sommet vu : %d\n",nb);
                //printf("call : %d\n",call);
            #endif
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
                    if (!(neighbor.x >= 0 && neighbor.x < GRID_SIZE && neighbor.y >= 0 && neighbor.y < GRID_SIZE && neighbor.z >= 0 && neighbor.z < GRID_SIZE) 
                    || vMap[minN.point.x + dx][minN.point.y + dy][minN.point.z + dz].indexHeap==-1 
                    || voxelGrid[neighbor.x][neighbor.y][neighbor.z] == 1) {
                        continue;
                    }
                    float tentativeGScore = minN.g + sqrt(abs(minN.point.x-neighbor.x)+abs(minN.point.y-neighbor.y)+abs(minN.point.z-neighbor.z))*LENGTH_GRID;
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
    return -1; // Return -1 if path not found
}

/**************************************/
/* Util  ******************************/
/**************************************/

float distWithObstacles(Point_t startPos, Point_t endPos, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap){
    Point3D endPoint = createPoint3D(endPos);
	float EndDist = dist(endPos, createPoint_t(endPoint));
	Point3D startPoint = createPoint3D(startPos);
	float startDist = dist(startPos, createPoint_t(startPoint));
	float computedDist = aStarPathfinding(startPoint, endPoint, voxelGrid, vMap, nodeHeap);
	return computedDist + EndDist + startDist;
}