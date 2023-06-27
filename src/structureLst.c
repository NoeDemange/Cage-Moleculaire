#include "structure.h"
#include "util.h"
#include <limits.h>

/**************************************/
/* LISTE ******************************/
/**************************************/
void LST_init(List_t* l) {
	l->elts = NULL;
	l->size = 0;
}

void LST_addAlloc(List_t* l) {
	int i;

	l->elts = realloc(l->elts, (size(l)+REALLOCSIZE)*sizeof(int));

	for (i=0; i<REALLOCSIZE; i++)
		l->elts[l->size+i] = -1;

	size(l) += REALLOCSIZE;
}

unsigned LST_nbElements(List_t* l) {
	int cpt;

	for (cpt=size(l); cpt>0 && elts(l,cpt-1) == -1; cpt--);

	return cpt;
}

unsigned LST_getIndiceFree(List_t* l) {
	int i;

	for (i=0; i<size(l); i++)
		if (elts(l,i) == -1)
			return i;

	LST_addAlloc(l);
	return i;
}

unsigned LST_getIndice(List_t* l, unsigned id) {
	int i;

	for (i=0; i<size(l) && elts(l,i) != -1; i++)
		if (elts(l,i) == id)
			return i;

	return -1;
}

unsigned LST_check(List_t* l, unsigned id) {

	if (LST_getIndice(l,id) == -1)
		return 0;
	return 1;
}

void LST_addElement(List_t* l, unsigned id) {
	
	int i;
	if (LST_getIndice(l, id) == -1) {
		i = LST_getIndiceFree(l);
		l->elts[i] = id;
	}
}

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


List_t* LST_create() {

	List_t* l = malloc(sizeof(List_t));
	LST_init(l);

	return l;
}

List_t* LST_copy(List_t* l) {

	int i;
	List_t* copy = LST_create();

	size(copy) = LST_nbElements(l);

	copy->elts = malloc(size(copy)*sizeof(int));

	for (i=0; i<size(copy); i++)
		elts(copy,i) = elts(l,i);

	return copy;
}

List_t* LST_copyWithShift(List_t* l, int* mod_pos_nei) { //copie de la liste avec décalage de la numérotation de l'éléments en fonction de la valeur du taleau mod_pos_nei

	int i;
	List_t* copy = LST_create();

	size(copy) = LST_nbElements(l);

	copy->elts = malloc(size(copy)*sizeof(int));

	for (i=0; i<size(copy); i++)
		elts(copy,i) = elts(l,i)-mod_pos_nei[elts(l,i)];

	return copy;
}

List_t* LST_addList(List_t* l1, List_t* l2) {

	int i;

	List_t* out = LST_copy(l1);

	for (i=0; i<size(l2) && elts(l2,i)!=-1; i++)
		LST_addElement(out, elts(l2,i));

	LST_delete(l2);
	LST_delete(l1);
	return out;
}

void LST_delete(List_t* l) {

	free(l->elts);
	free(l);
}

/******************************/

List_p* LST2_init() {
	
	List_p* list = malloc(sizeof(List_p));
	list->first = NULL;
	
	return list;
}

// Ajout au début
void LST2_addElement(List_p* list, int depart, int arrivee) {
	
	Element* elem = malloc(sizeof(Element));
	
	elem->start = depart;
	elem->end = arrivee;
	elem->next = list->first;
	
	list->first = elem;
}

void LST2_removeFirst(List_p* list) {
	
	Element* suppr = list->first;
	list->first = list->first->next;
	free(suppr);
}

void LST2_delete(List_p* list) {

	while (list->first)
	{
		LST2_removeFirst(list);
	}
	free(list);
}

/******************************/

List_m* LSTm_init() {
	
	List_m* list = malloc(sizeof(List_m));
	list->first = NULL;
	
	return list;
}

// Ajout au début
void LSTm_addElement(List_m* list, Shell_t* moc) {
	
	Elem* elem = malloc(sizeof(Elem));
	
	elem->moc = moc;
	elem->next = list->first;
	
	list->first = elem;
}

void LSTm_removeFirst(List_m* list) {
	
	Elem* suppr = list->first;
	list->first = list->first->next;
	if(suppr->moc) SHL_delete(suppr->moc);
	free(suppr);
}

void LSTm_delete(List_m* list) {

	while (list->first)
	{
		LSTm_removeFirst(list);
	}
	free(list);
}

/******************************/

List_s* LSTs_init() {
	
	List_s* list = malloc(sizeof(List_s));
	list->first = NULL;
	
	return list;
}

// Ajout au début
void LSTs_addElement(List_s* list, Point_t sommet) {
	
	Elem_s* elem = malloc(sizeof(Elem_s));
	
	elem->position.x = sommet.x;
	elem->position.y = sommet.y;
	elem->position.z = sommet.z;
	elem->next = list->first;
	
	list->first = elem;
}

void LSTs_removeFirst(List_s* list) {
	
	Elem_s* suppr = list->first;
	list->first = list->first->next;
	free(suppr);
}

void LSTs_delete(List_s* list) {

	while (list->first)
	{
		LSTs_removeFirst(list);
	}
	free(list);
}

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

// Retourne le point de la liste le plus proche du point en argument
Point_t minDist(List_s* list, Point_t p) {
	Point_t min = PT_init();
	int distanceMin = INT_MAX;
	Elem_s* l = list->first;
	
	while (l)
	{
		if (dist(l->position, p) < distanceMin)
		{
			min = l->position;
			distanceMin = dist(l->position, p);
		}
		l = l->next;
	}
	return min;
}

/******************************/

List_d* LSTd_init() {
	
	List_d* list = malloc(sizeof(List_d));
	list->first = NULL;
	
	return list;
}

// Ajout au début
void LSTd_addElement(List_d* list, int sommet) {
	
	Elem_d* elem = malloc(sizeof(Elem_d));
	
	elem->idAtom = sommet;
	elem->next = list->first;
	
	list->first = elem;
}

void LSTd_removeFirst(List_d* list) {
	
	Elem_d* suppr = list->first;
	list->first = list->first->next;
	free(suppr);
}

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

void LSTd_delete(List_d* list) {

	while (list->first)
	{
		LSTd_removeFirst(list);
	}
	free(list);
}
