#ifndef __INPUT_H
#define __INPUT_H

#include "structure.h"

/**
 * @file input.h
 * @brief Input Handling
 *
 * This file contains the function for handling input data.
 */

/**
 * @brief Read a molecule from a file in XYZ format.
 * 
 * This function reads molecular data from a file in the XYZ format and initializes a molecule with it.
 * The XYZ file should contain information about the atoms and their coordinates in 3D space.
 * 
 * @param filename The name of the file containing the molecular data in the XYZ format.
 * @return A pointer to the initialized molecule.
 */
Molecule_t* readInput_xyz(char*);

/**
 * @brief Read covalent radii of atoms from a file.
 * 
 * This function retrieves the covalent radii of atoms from the "rdc.dat" file and stores them in the
 * corresponding atoms of the molecule.
 * 
 * @param m Pointer to the molecule to which the covalent radii will be added.
 */
void readCovalence(Molecule_t*);

#endif
