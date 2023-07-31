#ifndef __INITIALIZATION_H
#define __INITIALIZATION_H

#include "structure.h"

/**
 * @file initialization.h
 * @brief Initialization of Molecules
 *
 * This file contains the function for initializing molecules.
 * It provides a function to create a molecule from a file by reading the molecule's coordinates and attributes.
 */

/**
 * @brief Initialize a molecule from a file.
 * 
 * This function initializes a molecule by reading its data from a file with the specified filename.
 * The file should be in the XYZ format. It loads the molecular data and creates the atoms and bonds 
 * required to represent the molecule.
 * 
 * @param filename The name of the file containing the molecular data in the XYZ format.
 * @return A pointer to the initialized molecule.
 */
Molecule_t* initMolecule(char*);

#endif
