#ifndef __EXPANSION_H
#define __EXPANSION_H

#include "structure.h"

/** @file expansion.h
 *  @brief Functions for shell creation and expansion.
 */

Shell_t* createShell(Molecule_t* m, double alpha);
void expansion(Molecule_t* m, Shell_t* s);

#endif
