#ifndef __EXPANSION_H
#define __EXPANSION_H

#include "structure.h"

Shell_t* createShell(Molecule_t*, double);
void expansion(Molecule_t* m, Shell_t* s);

#endif
