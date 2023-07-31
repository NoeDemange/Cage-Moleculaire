#include "structure.h"

/**
 * @file structureAsp.c
 * @brief Ashape_t functions for allocating and deallocating memory.
 *
 * This file contains functions for allocating and deallocating memory for an Ashape_t structure.
 */

/**************************************/
/* ASHAPE3D ***************************/
/**************************************/

/**
 * @brief Allocates memory for a new Ashape_t structure.
 *
 * This function allocates memory for a new Ashape_t structure and initializes its fields.
 *
 * @return A pointer to the newly allocated Ashape_t structure.
 */
Ashape_t* ASP_create() {
	Ashape_t* as3d = malloc(sizeof(Ashape_t));

	as3d->nb_triang = 0;
	as3d->nb_edge = 0;
	as3d->nb_vertex = 0;
	as3d->nb_x = 0;
	as3d->nb_alpha = 0;
	
	as3d->triang = NULL;
	as3d->edge = NULL;
	as3d->vertex = NULL;
	as3d->x = NULL;
	as3d->alpha = NULL;

	return as3d;
}

/**
 * @brief Frees memory used by the Ashape_t structure.
 *
 * This function frees the memory used by the Ashape_t structure and its associated arrays.
 *
 * @param as3d Pointer to the Ashape_t structure to be destroyed.
 */
void ASP_delete(Ashape_t* as3d) {

	if (as3d != NULL) {
		free(as3d->triang);
		free(as3d->edge);
		free(as3d->vertex);
		free(as3d->x);
		free(as3d->alpha);
	}

	free(as3d);
}