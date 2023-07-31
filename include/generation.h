#ifndef __GENERATION_H
#define __GENERATION_H

#include "structure.h"

/**
 * @file generation.h
 * @brief Pathless Cage Generation
 *
 * This file contains the function for generating pathless cages. 
 * It adds aromatic rings and hydrogen patterns to the envelope to create pathless cages.
 */

/**
 * @brief Generate pathless cages by adding aromatic rings and hydrogen patterns to the envelope.
 * 
 * This function takes a grouping of the main structures, which includes the substrate and envelope.
 * It generates pathless cages by adding aromatic rings and hydrogen patterns to the envelope.
 * The function modifies the envelope to form pathless cages without creating any path between atoms.
 * 
 * @param m Grouping of the main structures (substrate and envelope).
 */
void generatePathlessCages(Main_t*);

#endif
