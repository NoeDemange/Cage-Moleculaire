#ifndef __ASSEMBLY_H
#define __ASSEMBLY_H

#include "structure.h"
#include "main.h"

/** @file assembly.h
 *  @brief Functions related to generating whole cages.
 */

/** @brief Generates whole cages.
 *  
 *  This function generates whole cages by assembling atoms according to the specified options.
 * 
 *  @param mainData A pointer to the Main_t structure containing program options and data.
 *  @param options Options_t structure containing cage assembly options.
 */
void generateWholeCages(Main_t*, Options_t);

#endif
