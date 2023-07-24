#include <stdio.h>
#include "voxelization.h"
#include "constant.h"
#include "output.h"

int isPointInsideSphere(Point_t point, Point_t center, float radius) {
    float dx = point.x - center.x;
    float dy = point.y - center.y;
    float dz = point.z - center.z;
    float distanceSquared = dx * dx + dy * dy + dz * dz;
    float radiusSquared = radius * radius;
    if(distanceSquared <= radiusSquared) return 1;
    return 0;
}

VOXELGRID initVoxelGrid(){
    VOXELGRID voxelGrid = (int***)calloc(GRID_SIZE, sizeof(int**));
    if(voxelGrid == NULL){
        printf("problem calloc\n");
        exit(EXIT_FAILURE);
    }
    for (int x = 0; x < GRID_SIZE; x++) {
        voxelGrid[x] = (int**)calloc(GRID_SIZE, sizeof(int*));
        for (int y = 0; y < GRID_SIZE; y++) {
            voxelGrid[x][y] = (int*)calloc(GRID_SIZE, sizeof(int));
        }
    }
    return voxelGrid;
}

VOXELGRID voxelization(Molecule_t* sub) {
    VOXELGRID voxelGrid = initVoxelGrid();
    for (int i = 0; i < size(sub); i++) {
        float radius = DIST_GAP_SUBSTRATE;
        Point_t center = coords(atom(sub,i));
        for (int x = (int)((abs(START_GRID) + (center.x - radius))/(LENGTH_GRID)); x <= (int)((abs(START_GRID) + (center.x + radius))/(LENGTH_GRID))+1; x++) {
            for (int y = (int)((abs(START_GRID) + (center.y - radius))/(LENGTH_GRID)); y <= (int)((abs(START_GRID) + (center.y + radius))/(LENGTH_GRID))+1; y++) {
                for (int z = (int)((abs(START_GRID) + (center.z - radius))/(LENGTH_GRID)); z <= (int)((abs(START_GRID) + (center.z + radius))/(LENGTH_GRID))+1; z++) {
                    Point_t current = {(START_GRID + LENGTH_GRID * x), (START_GRID + LENGTH_GRID * y), (START_GRID + LENGTH_GRID * z)};
                    if (x >= 0 && x < GRID_SIZE &&
                        y >= 0 && y < GRID_SIZE &&
                        z >= 0 && z < GRID_SIZE &&
                        isPointInsideSphere(current,center,radius)) {
                        voxelGrid[x][y][z] = 1;
                    }
                }
            }
        }
    }
    return voxelGrid;
}


void printVoxelGrid(VOXELGRID voxelGrid){
    Shell_t* voxelVisu = SHL_create();
    for(int x = 0; x<GRID_SIZE; x++){
        for(int y = 0; y<GRID_SIZE; y++){
            for(int z = 0; z<GRID_SIZE; z++){
                if(voxelGrid[x][y][z]==1){
                    Point_t Ptest = {START_GRID + ((LENGTH_GRID) * x), START_GRID + ((LENGTH_GRID) * y), START_GRID + ((LENGTH_GRID) * z)};
                    SHL_addAtom(voxelVisu,Ptest,-1);
                }
            }
        }
    }
    SHL_writeMol2("voxelTest_neg.mol2",voxelVisu);
}

void freeVoxelGrid(VOXELGRID voxelGrid){
    // Free the dynamically allocated memory
    for (int x = 0; x < GRID_SIZE; x++) {
        for (int y = 0; y < GRID_SIZE; y++) {
            free(voxelGrid[x][y]);
        }
        free(voxelGrid[x]);
    }
    free(voxelGrid);
}
