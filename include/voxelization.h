#ifndef __VOXELIZATION_H
#define __VOXELIZATION_H

#include "structure.h"

VOXELGRID initVoxelGrid();
VOXELGRID voxelization(Molecule_t* sub);
void freeVoxelGrid(VOXELGRID voxelGrid);
void printVoxelGrid(VOXELGRID voxelGrid);

#endif