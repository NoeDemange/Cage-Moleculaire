#ifndef __VOXELIZATION_H
#define __VOXELIZATION_H

#include "structure.h"

/**
 * @file voxelization.h
 * @brief Voxelization Header File
 *
 * This file contains functions for voxelization-related operations.
 */

/**
 * @brief Initialize a 3D voxel grid.
 *
 * This function initializes a 3D voxel grid by allocating memory for the grid and setting all
 * voxels to zero.
 *
 * @return The initialized 3D voxel grid (VOXELGRID).
 */
VOXELGRID initVoxelGrid();

/**
 * @brief Perform voxelization of a molecular structure.
 *
 * This function performs voxelization of a given molecule by marking the voxels that contain
 * atoms from the molecular structure. It uses the atom coordinates and radius to determine
 * the voxelization.
 *
 * @param sub The Molecule_t pointer representing the molecular structure.
 * @return The voxelized 3D grid (VOXELGRID) with voxels containing atoms marked as 1.
 */
VOXELGRID voxelization(Molecule_t* sub);

/**
 * @brief Free the dynamically allocated memory of a voxel grid.
 *
 * This function frees the memory dynamically allocated for the voxel grid.
 *
 * @param voxelGrid The voxelized 3D grid (VOXELGRID) to be freed.
 */
void freeVoxelGrid(VOXELGRID voxelGrid);

/**
 * @brief Print the voxel grid for visualization.
 *
 * This function creates a shell structure representing the voxel grid and writes it to a .mol2 file
 * for visualization purposes.
 *
 * @param voxelGrid The voxelized 3D grid (VOXELGRID) to be visualized.
 */
void printVoxelGrid(VOXELGRID voxelGrid);

#endif