#include "structure.h"
#include "output.h"

/**
 * @file structureMN.c
 * @brief Structure Main Operations
 * 
 * This file contains the declarations of functions related to main operations
 * in the application. These operations involve initializing, creating,
 * deleting, and managing the Main_t structure and its elements (Shell_t, Molecule_t).
 * It also includes functions for managing the allocation and copying of Shell_t elements.
 */

/**
 * @brief Initializes a new Shell_t (moc).
 *
 * This function initializes a new Shell_t (moc) by setting its size to 0 and
 * its atoms, cycle, and bond pointers to NULL.
 *
 * @param s Pointer to the Shell_t (moc) to be initialized.
 */
void MN_initMoc(Shell_t* s) {

	s->size = 0;
	s->atoms = NULL;
	s->cycle = NULL;
	s->bond = NULL;
}

/**
 * @brief Adds memory allocation for mocs in the Main_t.
 *
 * This function adds memory allocation for additional mocs in the Main_t structure.
 * The size parameter specifies the number of additional mocs to be allocated.
 *
 * @param m Pointer to the Main_t structure.
 * @param size Number of additional mocs to be allocated.
 */
void MN_addAlloc(Main_t* m, unsigned size) {

	int i;

	m->mocs = realloc(m->mocs, (mocSize(m)+size)*sizeof(Shell_t*));

	for (i=0; i<size; i++) {
		moc(m,mocSize(m)+i) = NULL;
	}

	mocSize(m) += size;
}

/**
 * @brief Gets the index of the next free moc in the Main_t.
 *
 * This function searches for the next free moc in the Main_t structure and
 * returns its index. If there are no free mocs, it allocates additional memory
 * to accommodate more mocs.
 *
 * @param m Pointer to the Main_t structure.
 * @return Index of the next free moc in the mocs array.
 */
unsigned MN_getIndiceFree(Main_t* m) {

	int i;

	for (i=0; i<mocSize(m); i++) {
		if (moc(m,i) == NULL)
			return i;
	}
	MN_addAlloc(m, REALLOCSIZE);
	return i;
}

/**
 * @brief Gets the index of the next free moc in the Main_t (adding only one moc).
 *
 * This function searches for the next free moc in the Main_t structure and
 * returns its index. If there are no free mocs, it allocates memory for one more moc.
 *
 * @param m Pointer to the Main_t structure.
 * @return Index of the next free moc in the mocs array.
 */
unsigned MN_getIndiceFree2(Main_t* m) {

	int i;

	for (i=0; i<mocSize(m); i++) {
		if (moc(m,i) == NULL)
			return i;
	}
	MN_addAlloc(m, 1);
	return i;
}

/**
 * @brief Copies a Shell_t (moc) and adds it to the Main_t.
 *
 * This function copies a given Shell_t (moc) and adds it to the Main_t structure.
 * It returns the index of the added moc in the mocs array.
 *
 * @param m Pointer to the Main_t structure.
 * @param s Pointer to the Shell_t (moc) to be copied and added.
 * @return Index of the added moc in the mocs array.
 */
unsigned MN_copyMoc(Main_t* m, Shell_t* s) {

	unsigned indice = MN_getIndiceFree(m);
	moc(m, indice) = SHL_copy(s);

	return indice;
}

/**
 * @brief Creates a new Main_t structure.
 *
 * This function allocates memory for a new Main_t structure and initializes its
 * elements (substrate, envelope, envarom, mocs) to NULL or zero.
 *
 * @return Pointer to the newly allocated Main_t structure.
 */
Main_t* MN_create() {

	Main_t* m = malloc(sizeof(Main_t));
	m->substrat = NULL;
	m->envelope = NULL;
	m->envarom = NULL;
	m->mocs = NULL;
	m->mocSize = 0;

	return m;
}

/**
 * @brief Deletes a Main_t structure and its associated elements.
 *
 * This function deletes a Main_t structure and frees the memory associated with its
 * elements (substrate, envelope, envarom, mocs).
 *
 * @param m Pointer to the Main_t structure to be deleted.
 */
void MN_delete(Main_t* m) {

	int i;

	if (substrat(m) != NULL)
		MOL_delete(substrat(m));

	if (envelope(m) != NULL)
		SHL_delete(envelope(m));

	if (envarom(m) != NULL)
		SHL_delete(envarom(m));

	if (m->mocs != NULL) {
		for (i=0; i<mocSize(m); i++) {
			if (moc(m,i) != NULL)
				SHL_delete(moc(m,i));
		}
		free(m->mocs);
	}

	free(m);
}
