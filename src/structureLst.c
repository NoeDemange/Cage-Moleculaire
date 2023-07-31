#include "structure.h"
#include "util.h"
#include "pathFinding.h"
#include <limits.h>

/**
 * @file structureLst.c
 * @brief List operations functions for manipulating linked lists.
 *
 * This file contains functions for creating, manipulating, and checking linked lists.
 */

/**************************************/
/* LISTE ******************************/
/**************************************/
/**
 * @brief Initializes a new list.
 *
 * This function initializes a new list by setting its elements to NULL and size to 0.
 *
 * @param l Pointer to the list to initialize.
 */
void LST_init(List_t* l) {
	l->elts = NULL;
	l->size = 0;
}

/**
 * @brief Adds allocated memory to the list.
 *
 * This function adds allocated memory to the list to accommodate new elements.
 *
 * @param l Pointer to the list to which memory is added.
 */
void LST_addAlloc(List_t* l) {
	int i;

	l->elts = realloc(l->elts, (size(l)+REALLOCSIZE)*sizeof(int));

	for (i=0; i<REALLOCSIZE; i++)
		l->elts[l->size+i] = -1;

	size(l) += REALLOCSIZE;
}

/**
 * @brief Returns the number of elements in the list.
 *
 * This function returns the number of elements in the list.
 *
 * @param l Pointer to the list.
 * @return The number of elements in the list.
 */
unsigned LST_nbElements(List_t* l) {
	int cpt;

	for (cpt=size(l); cpt>0 && elts(l,cpt-1) == -1; cpt--);

	return cpt;
}

/**
 * @brief Returns the index of the first available free element in the list.
 *
 * This function returns the index of the first available free element in the list.
 * If no free element is available, it adds and allocates memory for more elements.
 *
 * @param l Pointer to the list.
 * @return The index of the first available free element in the list.
 */
unsigned LST_getIndiceFree(List_t* l) {
	int i;

	for (i=0; i<size(l); i++)
		if (elts(l,i) == -1)
			return i;

	LST_addAlloc(l);
	return i;
}

/**
 * @brief Returns the index of the element with the specified identifier in the list.
 *
 * This function returns the index of the element with the specified identifier in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to search for.
 * @return The index of the element with the specified identifier, or -1 if not found.
 */
unsigned LST_getIndice(List_t* l, unsigned id) {
	int i;

	for (i=0; i<size(l) && elts(l,i) != -1; i++)
		if (elts(l,i) == id)
			return i;

	return -1;
}

/**
 * @brief Checks if an element with the specified identifier exists in the list.
 *
 * This function checks if an element with the specified identifier exists in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to check for.
 * @return 1 if the element with the specified identifier exists, 0 otherwise.
 */

unsigned LST_check(List_t* l, unsigned id) {

	if (LST_getIndice(l,id) == -1)
		return 0;
	return 1;
}

/**
 * @brief Adds an element to the list if it does not exist.
 *
 * This function adds an element to the list if it does not already exist.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be added.
 */
void LST_addElement(List_t* l, unsigned id) {
	
	int i;
	if (LST_getIndice(l, id) == -1) {
		i = LST_getIndiceFree(l);
		l->elts[i] = id;
	}
}

/**
 * @brief Removes an element with the specified identifier from the list.
 *
 * This function removes an element with the specified identifier from the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be removed.
 */
void LST_removeElement(List_t* l, unsigned id) {

	int i = LST_getIndice(l,id);

	if (i != -1) {
		while (i < size(l)-1 && elts(l,i) != -1) {
			elts(l,i) = elts(l,i+1);
			i++;
		}
		elts(l,i) = -1;
	}
}

/**
 * @brief Creates a new integer list.
 *
 * This function allocates memory for a new integer list and initializes its fields.
 *
 * @return Pointer to the newly allocated integer list.
 */
List_t* LST_create() {

	List_t* l = malloc(sizeof(List_t));
	LST_init(l);

	return l;
}

/**
 * @brief Creates a copy of an integer list.
 *
 * This function creates a copy of an integer list by copying its elements.
 *
 * @param l Pointer to the integer list to be copied.
 * @return Pointer to the newly created copy of the integer list.
 */
List_t* LST_copy(List_t* l) {

	int i;
	List_t* copy = LST_create();

	size(copy) = LST_nbElements(l);

	copy->elts = malloc(size(copy)*sizeof(int));

	for (i=0; i<size(copy); i++)
		elts(copy,i) = elts(l,i);

	return copy;
}

/**
 * @brief Copies the list with an offset of the numbering of elements 
 * according to the value of the array passed as an argument.
 * 
 * @param l List to copy.
 * @param shifts Array of offsets.
 * @return (List_t*) Copied list with shifts. 
 */
List_t* LST_copyWithShift(List_t* l, int* shifts) {

	int i;
	List_t* copy = LST_create();

	size(copy) = LST_nbElements(l);

	copy->elts = malloc(size(copy) * sizeof(int));

	for (i = 0; i < size(copy); i++)
		elts(copy,i) = elts(l,i) - shifts[elts(l,i)];

	return copy;
}

/**
 * @brief Merges two integer lists.
 *
 * This function merges two integer lists, adding elements from the second list to the first.
 * The function then deletes both input lists and returns the merged list.
 *
 * @param l1 First integer list.
 * @param l2 Second integer list.
 * @return Merged integer list.
 */
List_t* LST_addList(List_t* l1, List_t* l2) {

	int i;

	List_t* out = LST_copy(l1);

	for (i=0; i<size(l2) && elts(l2,i)!=-1; i++)
		LST_addElement(out, elts(l2,i));

	LST_delete(l2);
	LST_delete(l1);
	return out;
}

/**
 * @brief Deletes an integer list and frees the associated memory.
 *
 * This function deletes an integer list and frees the memory associated with its elements.
 *
 * @param l Pointer to the integer list to be deleted.
 */
void LST_delete(List_t* l) {

	free(l->elts);
	free(l);
}

/******************************/
/**
 * @brief Initializes a new Element list.
 *
 * This function initializes a new Element list by setting its first element to NULL.
 *
 * @return Pointer to the newly allocated Element list.
 */
Element* LST_pairs_init(void) {

	Element* list = malloc(sizeof(Element));
	list = NULL;
	return list;
}

// Add at the 
/**
 * @brief Adds an Element to the beginning of the list.
 *
 * This function adds an Element to the beginning of the list.
 *
 * @param list Pointer to the Element list.
 * @param start Start index of the element to add.
 * @param end End index of the element to add.
 */
void LST_pairs_addElement(Element** list, int start, int end) {
	
	Element* elem = malloc(sizeof(Element));
	elem->start = start;
	elem->end = end;
	elem->next = *list;

	*list = elem;
}

/**
 * @brief Adds the pairs of atoms to connect by sorting them 
 * by ascending A* distance.
 * 
 * @param s Cage in progress.
 * @param list List of pairs of atoms to connect.
 * @param start Index of the starting atom of the path in the cage.
 * @param end Index of the ending atom of the path in the cage.
 * @param voxelGrid Grid of voxelization.
 */
void LST_pairs_addElementInOrder(Shell_t* s, Element** list, int start, int end, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap) {
	Point_t endPos = coords(atom(s, end));
	Point_t startPos = coords(atom(s, start));
	float computedDist = distWithObstacles(startPos,endPos,voxelGrid,vMap,nodeHeap);
	Element* currentElem = *list;
	Element* previousElem = NULL;
	while(currentElem){
		if(currentElem->distance < computedDist){
			previousElem = currentElem;
			currentElem = currentElem->next;
		}
		else{break;}
	}

	Element* elem = malloc(sizeof(Element));
	elem->start = start;
	elem->end = end;
	elem->distance = computedDist;
	elem->next = currentElem;
	if (previousElem) {
		previousElem->next = elem;
	}
	else {
		*list = elem;
	}
}

/**
 * @brief Removes the first Element from the list.
 *
 * This function removes the first Element from the list and frees the associated memory.
 *
 * @param list Pointer to the Element list.
 */
void LST_pairs_removeFirst(Element* list) {
	
	if (list) {
		Element* delete = list;
		list = list->next;
		free(delete);
	}
	else {
		fprintf(stderr, "Can't remove first, list is empty.\n");
	}
}

/**
 * @brief Deletes the Element list and frees the associated memory.
 *
 * This function deletes the Element list and frees the memory associated with its elements.
 *
 * @param list Pointer to the Element list to be deleted.
 */
void LST_pairs_delete(Element* list) {

	if (list) {
		while (list) {
			Element* delete = list;
			list = list->next;
			free(delete);
		}
	}
	else {
		fprintf(stderr, "List is already deleted.\n");
	}
}
/******************************/

/**
 * @brief Initializes a new List_m (Shell_t) list.
 *
 * This function initializes a new List_m (used with Shell_t elements) by setting its first element to NULL.
 *
 * @return Pointer to the newly allocated List_m.
 */
List_m* LSTm_init() {
	
	List_m* list = malloc(sizeof(List_m));
	list->first = NULL;
	
	return list;
}

// Ajout au début
/**
 * @brief Adds a Shell_t element at the beginning of the list.
 *
 * This function adds a Shell_t element to the list.
 *
 * @param list Pointer to the list.
 * @param moc Pointer to the Shell_t element to add.
 */
void LSTm_addElement(List_m* list, Shell_t* moc) {
	
	Elem* elem = malloc(sizeof(Elem));
	
	elem->moc = moc;
	elem->next = list->first;
	
	list->first = elem;
}

/**
 * @brief Removes the first Shell_t element from the list.
 *
 * This function removes the first Shell_t element from the list and frees the associated memory.
 *
 * @param list Pointer to the list.
 */
void LSTm_removeFirst(List_m* list) {
	
	Elem* suppr = list->first;
	list->first = list->first->next;
	if(suppr->moc) SHL_delete(suppr->moc);
	free(suppr);
}

/**
 * @brief Deletes the List_m list and frees the associated memory.
 *
 * This function deletes the List_m list and frees the memory associated with its elements.
 *
 * @param list Pointer to the List_m list to be deleted.
 */
void LSTm_delete(List_m* list) {

	while (list->first)
	{
		LSTm_removeFirst(list);
	}
	free(list);
}

/******************************/

/**
 * @brief Initializes a new List_s (Point_t) list.
 *
 * This function initializes a new List_s (used with Point_t elements) by setting its first element to NULL.
 *
 * @return Pointer to the newly allocated List_s.
 */
List_s* LSTs_init() {
	
	List_s* list = malloc(sizeof(List_s));
	list->first = NULL;
	
	return list;
}

// Ajout au début
/**
 * @brief Adds a Point_t element at the beginning of the List_s (Point_t) list.
 *
 * This function adds a Point_t element to the List_s (Point_t) list.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param sommet Point_t element to add.
 */
void LSTs_addElement(List_s* list, Point_t sommet) {
	
	Elem_s* elem = malloc(sizeof(Elem_s));
	
	elem->position.x = sommet.x;
	elem->position.y = sommet.y;
	elem->position.z = sommet.z;
	elem->next = list->first;
	
	list->first = elem;
}

/**
 * @brief Adds a Point_t element to the List_s (Point_t) list in ascending order of A* distance.
 *
 * This function adds a Point_t element to the List_s (Point_t) list in ascending order of A* distance.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param startPos Position of the atom.
 * @param endPos Position of the objective atom.
 * @param voxelGrid Grid of voxelization.
 */
void LSTs_addElementInOrder(List_s* list, Point_t startPos, Point_t endPos, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap) {
	float computedDist = distWithObstacles(startPos,endPos,voxelGrid,vMap,nodeHeap);
	Elem_s* currentElem = list->first;
	Elem_s* previousElem = NULL;
	while(currentElem){
		if(currentElem->distance < computedDist){
			previousElem = currentElem;
			currentElem = currentElem->next;
		}
		else{break;}
	}

	Elem_s* elem = malloc(sizeof(Elem_s));
	elem->position = startPos;
	elem->distance = computedDist;
	elem->next = currentElem;
	if (previousElem) {
		previousElem->next = elem;
	}
	else {
		list->first = elem;
	}
}

/**
 * @brief Removes the first Point_t element from the List_s (Point_t) list.
 *
 * This function removes the first Point_t element from the List_s (Point_t) list and frees the associated memory.
 *
 * @param list Pointer to the List_s (Point_t) list.
 */
void LSTs_removeFirst(List_s* list) {
	
	Elem_s* suppr = list->first;
	list->first = list->first->next;
	free(suppr);
}

/**
 * @brief Deletes the List_s (Point_t) list and frees the associated memory.
 *
 * This function deletes the List_s (Point_t) list and frees the memory associated with its elements.
 *
 * @param list Pointer to the List_s (Point_t) list to be deleted.
 */
void LSTs_delete(List_s* list) {

	while (list->first)
	{
		LSTs_removeFirst(list);
	}
	free(list);
}

/**
 * @brief Removes a Point_t element from the List_s (Point_t) list.
 *
 * This function removes a Point_t element from the List_s (Point_t) list.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param p Point_t element to be removed.
 */
void LSTs_removeElement(List_s* list, Point_t p) {
	
	Elem_s* cursor = list->first;
	Elem_s* suppr = NULL;
	if (cursor)
	{
		if (cursor->position.x == p.x && cursor->position.y == p.y && cursor->position.z == p.z)
		{
			LSTs_removeFirst(list);
		}
		else
		{
			while (cursor->next && !suppr)
			{
				if (cursor->next->position.x == p.x && cursor->next->position.y == p.y && cursor->next->position.z == p.z)
				{
					suppr = cursor->next;
					cursor->next = cursor->next->next;
				}
				else cursor = cursor->next;
			}
		}
		
	}
	if (suppr)
	{
		free(suppr);
	}
}

/*// Retourne le point de la liste le plus proche du point en argument
Point_t minDist(List_s* list, Point_t p) {
	Point_t min = PT_init(0);
	float minDist = __FLT_MAX__; //si déclarer comme un int permet un choix entre le plus optimal et un sous-opitmal (marche mieux que float pour l'instant)
	float computedDist;
	Elem_s* l = list->first;
	while (l) {
		computedDist = dist(l->position, p);
		if (computedDist < minDist) {
			min = l->position;
			minDist = computedDist;
		}
		l = l->next;
	}
	return min;
}*/

/******************************/

/**
 * @brief Initializes a new List_d (integer) list.
 *
 * This function initializes a new List_d (used with integers) by setting its first element to NULL.
 *
 * @return Pointer to the newly allocated List_d.
 */
List_d* LSTd_init() {
	
	List_d* list = malloc(sizeof(List_d));
	list->first = NULL;
	
	return list;
}

/**
 * @brief Adds an integer element at the beginning of the List_d (integer) list.
 *
 * This function adds an integer element to the List_d (integer) list.
 *
 * @param list Pointer to the List_d (integer) list.
 * @param sommet Integer element to add.
 */
void LSTd_addElement(List_d* list, int sommet) {
	
	Elem_d* elem = malloc(sizeof(Elem_d));
	
	elem->idAtom = sommet;
	elem->next = list->first;
	
	list->first = elem;
}

/**
 * @brief Removes the first integer element from the List_d (integer) list.
 *
 * This function removes the first integer element from the List_d (integer) list and frees the associated memory.
 *
 * @param list Pointer to the List_d (integer) list.
 */
void LSTd_removeFirst(List_d* list) {
	
	Elem_d* suppr = list->first;
	list->first = list->first->next;
	free(suppr);
}

/**
 * @brief Removes a specific integer element from the List_d (integer) list.
 *
 * This function removes a specific integer element from the List_d (integer) list.
 *
 * @param list Pointer to the List_d (integer) list.
 * @param sommet Integer element to be removed.
 */
void LSTd_removeSommet(List_d* list, int sommet) {
	Elem_d* cursor = list->first;
	Elem_d* suppr = NULL;
	if (cursor)
	{
		if (cursor->idAtom == sommet)
		{
			LSTd_removeFirst(list);
		}
		else
		{
			while (cursor->next && !suppr)
			{
				if (cursor->next->idAtom == sommet)
				{
					suppr = cursor->next;
					cursor->next = cursor->next->next;
				}
				else cursor = cursor->next;
			}
		}
		
	}
	if (suppr)
	{
		free(suppr);
	}
	
}

/**
 * @brief Deletes the List_d (integer) list and frees the associated memory.
 *
 * This function deletes the List_d (integer) list and frees the memory associated with its elements.
 *
 * @param list Pointer to the List_d (integer) list to be deleted.
 */
void LSTd_delete(List_d* list) {

	while (list->first)
	{
		LSTd_removeFirst(list);
	}
	free(list);
}
