#ifndef __INTERFACE_H
#define __INTERFACE_H

#include "structure.h"

Ashape_t* Cashape3d(Shell_t*, double);
Ashape_t* Cashape3d2(Shell_t* s, double alpha, Ashape_t** as3d2);
int* Cinashape3d(Ashape_t*, double*, int);

#endif
