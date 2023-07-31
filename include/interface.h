#ifndef __INTERFACE_H
#define __INTERFACE_H

#include "structure.h"

/**
 * @file interface.h
 * @brief Interface Functions
 *
 * This file contains function for interfacing with R.
 */

/**
 * @brief Call to R for generating 3D Alpha Shapes.
 *
 * This function takes an envelope with a cloud of points and calls R to generate
 * 3D Alpha Shapes based on the provided aphashape parameter.
 * 
 * @param s Pointer to the envelope that already has a cloud of points.
 * @param alpha The aphashape parameter used in generating 3D Alpha Shapes.
 * @return A pointer to the generated 3D Alpha Shape.
 */
Ashape_t* Cashape3d(Shell_t*, double alpha);

#endif
