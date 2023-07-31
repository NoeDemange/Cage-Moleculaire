#include "structure.h"
#include "util.h"
#include "output.h"

/**
 * @file structureShl.c
 * @brief Shell Data Structure Implementation
 *
 * This file contains the implementation of the Shell data structure and related operations.
 * The Shell data structure represents a cage or a molecular shell with atoms and their connections.
 */

/**************************************/
/* SHELL ******************************/
/**************************************/

//#define SHL_addVertex(s, id) GPH_addVertex(bond(s), id)

/**
 * @brief Initializes an AtomShl_t with default values.
 *
 * This function initializes an AtomShl_t with default values for its attributes.
 *
 * @param a Pointer to the AtomShl_t to be initialized.
 */
void SHL_initAtom(AtomShl_t* a) {

	flag(a) = NOT_DEF_F;

	atomX(a) = 0;
	atomY(a) = 0;
	atomZ(a) = 0;

	parentAtom(a) = -1;

	neighborhood(a) = LST_create();
}

/**
 * @brief Gets the number of neighbors in the neighborhood of an AtomShl_t.
 *
 * This function returns the number of neighbors in the neighborhood of the given AtomShl_t.
 *
 * @param a Pointer to the AtomShl_t.
 * @return The number of neighbors in the AtomShl_t's neighborhood.
 */
int SHL_nbNeighborhood(AtomShl_t* a) {

	return LST_nbElements(neighborhood(a));
}

/**
 * @brief Gets the index of a free neighbor in the neighborhood of an AtomShl_t.
 *
 * This function returns the index of a free neighbor in the neighborhood of the given AtomShl_t.
 *
 * @param a Pointer to the AtomShl_t.
 * @return The index of a free neighbor in the AtomShl_t's neighborhood.
 */
int SHL_getIndiceFreeNeighbor(AtomShl_t* a) {

	return LST_getIndiceFree(a->neighborhood);
}

/**
 * @brief Get the index of an atom in the neighborhood list.
 *
 * This function retrieves the index of an atom with the given ID in the neighborhood list
 * of the specified AtomShl_t structure.
 *
 * @param a Pointer to the AtomShl_t structure.
 * @param id The ID of the atom to search for in the neighborhood list.
 * @return The index of the atom in the neighborhood list, or -1 if not found.
 */
int SHL_getIndice(AtomShl_t* a, unsigned id) {

	return LST_getIndice(a->neighborhood, id);
}

/**
 * @brief Add a neighbor to the atom's neighborhood.
 *
 * This function adds a neighbor with the given ID to the neighborhood list of the specified AtomShl_t structure.
 *
 * @param a Pointer to the AtomShl_t structure.
 * @param id The ID of the neighbor atom to add to the neighborhood list.
 */
void SHL_addNeighbor(AtomShl_t* a, unsigned id) {

	LST_addElement(a->neighborhood, id);
}

/**
 * @brief Remove a neighbor from the atom's neighborhood.
 *
 * This function removes a neighbor with the given ID from the neighborhood list of the specified AtomShl_t structure.
 *
 * @param a Pointer to the AtomShl_t structure.
 * @param id The ID of the neighbor atom to remove from the neighborhood list.
 */
void SHL_removeNeighbor(AtomShl_t* a, unsigned id) {
	
	LST_removeElement(a->neighborhood, id);
}

/**
 * @brief Allocate and add new atoms to the Shell structure.
 *
 * This function allocates memory for new AtomShl_t structures and adds them to the Shell structure.
 * It dynamically reallocates memory for the atoms array in the Shell structure to accommodate the new atoms.
 * If memory reallocation fails, the function prints an error message and exits the program.
 *
 * @param s Pointer to the Shell_t structure.
 */
void SHL_addAllocAtom(Shell_t* s) {

	AtomShl_t* tmp = realloc(s->atoms, (size(s) + REALLOCSIZE) * sizeof(AtomShl_t));
	if (tmp == NULL) {       
    fprintf(stderr, "A problem occurred during the reallocation (structureShl.c:53).\n");
    exit(EXIT_FAILURE);
	}
	else {
    s->atoms = tmp;
	}

	for (int i = 0; i < REALLOCSIZE; i++) {
		SHL_initAtom(atom(s,size(s) + i));
	}

	size(s) += REALLOCSIZE;
}

/**
 * @brief Get the number of atoms in the Shell structure.
 *
 * This function returns the number of atoms in the specified Shell_t structure.
 *
 * @param s Pointer to the Shell_t structure.
 * @return The number of atoms in the Shell structure.
 */
int SHL_nbAtom(Shell_t* s) {
	int i, cpt = 0;

	for (i=0; i<size(s); i++)
		if (flag(atom(s,i)) != NOT_DEF_F)
			cpt++;

	return cpt;
}

/**
 * @brief Get the number of edges in the Shell structure.
 *
 * This function returns the number of edges (connections between atoms) in the specified Shell_t structure.
 *
 * @param s Pointer to the Shell_t structure.
 * @return The number of edges in the Shell structure.
 */
int SHL_nbEdges(Shell_t* s) {
	int i, cpt = 0;

	for (i=0; i<size(s); i++)
		cpt += SHL_nbNeighborhood(atom(s, i));

	return cpt/2;
}

/**
 * @brief Get the index of a free atom in the Shell structure.
 *
 * This function searches for a free (not defined) atom in the specified Shell_t structure and returns its index.
 * If no free atom is found, it allocates new atoms using SHL_addAllocAtom and returns the index of the first newly added atom.
 *
 * @param s Pointer to the Shell_t structure.
 * @return The index of the free atom in the Shell structure.
 */
int SHL_getIndiceFreeAtom(Shell_t* s) {
	int i;

	for (i=0; i<size(s); i++)
		if (flag(atom(s,i)) == NOT_DEF_F)
			return i;

	SHL_addAllocAtom(s);
	return i;
}

/**
 * @brief Adds a vertex to the Shell with the specified ID.
 *
 * This function adds a vertex with the given ID to the Shell data structure.
 *
 * @param s Pointer to the Shell data structure.
 * @param id The ID of the vertex to be added.
 * @return The ID of the added vertex.
 */
unsigned SHL_addVertex(Shell_t* s, unsigned id) {

	return GPH_addVertex(bond(s), id);
}

/**
 * @brief Removes a vertex from the Shell with the specified ID.
 *
 * This function removes a vertex with the given ID from the Shell data structure.
 *
 * @param s Pointer to the Shell data structure.
 * @param id The ID of the vertex to be removed.
 */
void SHL_removeVertex(Shell_t* s, unsigned id) {

	GPH_removeVertex(bond(s), id);
}

/**
 * @brief Adds a bond between two vertices with the specified IDs.
 *
 * This function adds a bond between two vertices (atoms) with the given IDs to the Shell data structure.
 *
 * @param s Pointer to the Shell data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void SHL_addBond(Shell_t* s, unsigned id1, unsigned id2) {

	GPH_addEdge(bond(s), id1, id2);
}

/**
 * @brief Removes a bond between two vertices with the specified IDs.
 *
 * This function removes a bond between two vertices (atoms) with the given IDs from the Shell data structure.
 *
 * @param s Pointer to the Shell data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void SHL_removeBond(Shell_t* s, unsigned id1, unsigned id2) {

	GPH_removeEdge(bond(s), id1, id2);
}

/**
 * @brief Adds an edge (bond) between two vertices with the specified IDs.
 *
 * This function adds an edge (bond) between two vertices (atoms) with the given IDs to the Shell data structure.
 * It also adds the reverse edge to create an undirected graph representation.
 *
 * @param s Pointer to the Shell data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void SHL_addEdge(Shell_t* s, unsigned id1, unsigned id2) {

	if (id1 < size(s) && id2 <size(s) && id1 != id2) {

		SHL_addNeighbor(atom(s, id1), id2);
		SHL_addNeighbor(atom(s, id2), id1);
	}
}

/**
 * @brief Removes an edge (bond) between two vertices with the specified IDs.
 *
 * This function removes an edge (bond) between two vertices (atoms) with the given IDs from the Shell data structure.
 * It also removes the reverse edge to maintain the undirected graph representation.
 *
 * @param s Pointer to the Shell data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void SHL_removeEdge(Shell_t* s, unsigned id1, unsigned id2) {

	if (id1 < size(s) && id2 <size(s)) {
		
		SHL_removeNeighbor(atom(s, id1), id2);
		SHL_removeNeighbor(atom(s, id2), id1);
	}
}

/**
 * @brief Adds a new atom to the Shell with the specified coordinates and parent atom ID.
 *
 * This function adds a new atom to the Shell data structure with the specified coordinates and parent atom ID.
 * It also assigns a unique ID (indice) to the new atom.
 *
 * @param s Pointer to the Shell data structure.
 * @param coords The coordinates (Point_t) of the new atom.
 * @param parent The ID of the parent atom for the new atom.
 * @return The unique ID (indice) assigned to the newly added atom.
 */
unsigned SHL_addAtom(Shell_t* s, Point_t coords, unsigned parent) {

	unsigned indice = SHL_getIndiceFreeAtom(s);

	flag(atom(s,indice)) = SHELL_F;
	coords(atom(s,indice)) = coords;
	parentAtom(atom(s,indice)) = parent;

	return indice;
}

/**
 * @brief Removes an atom from the Shell with the specified ID.
 *
 * This function removes an atom with the given ID from the Shell data structure.
 * It also updates the neighborhood and other properties accordingly.
 *
 * @param s Pointer to the Shell data structure.
 * @param id The ID of the atom to be removed.
 */
void SHL_removeAtom(Shell_t* s, unsigned id) {

	int i;

	if (id < size(s)) {
		AtomShl_t* a = atom(s,id);

		if (cycle(s,id))
			LST_removeElement(s->cycle, id);

		for (i=0; i<neighborhoodSize(a); i++)
			if (neighbor(a, i) != -1)
				SHL_removeNeighbor(atom(s, neighbor(a,i)), id);

		LST_delete(neighborhood(a));

		if (checkVertex(s,id))
			SHL_removeVertex(s, id);
		SHL_initAtom(a);
	}
}

/**
 * @brief Seeks the border atoms in the Shell for a given atom ID.
 *
 * This function seeks the border atoms in the Shell for a given atom ID and returns them as a list.
 *
 * @param s Pointer to the Shell data structure.
 * @param in Pointer to the input list.
 * @param id The ID of the atom to start seeking from.
 * @return A list of border atoms (border of the Shell) starting from the specified atom ID.
 */
List_t* SHL_seekBorder(Shell_t* s, List_t* in, unsigned id) {

	int i;
	List_t* out = LST_create();
	AtomShl_t* a = atom(s, id);

	if (flag(a) == SHELL_F || flag(a) == LINKABLE_F) {
		LST_addElement(out, id);
		return out;
	}

	LST_addElement(in, id);
	for (i=0; i<neighborhoodSize(a) && neighbor(a,i)!=-1; i++) {
		if (!LST_check(in, neighbor(a,i)))
			out = LST_addList(out, SHL_seekBorder(s, in, neighbor(a,i)));
	}

	return out;
}

/*void SHL_linkBorder(Shell_t* s, unsigned id, List_t* l) {

	int i, j, indiceMin;
	float distMin, dis;
	List_t* tmp = LST_create();
	List_t* border = SHL_seekBorder(s, tmp, id);

	LST_delete(tmp);

	for (i=0; i<size(l) && elts(l,i)!=-1; i++) {

		indiceMin = elts(border,0);
		distMin = dist(coords(atom(s, elts(l,i))), coords(atom(s, elts(border,0))));
		for (j=1; j<size(border) && elts(border,j)!=-1; j++) {

			dis = dist(coords(atom(s, elts(l,i))), coords(atom(s, elts(border,j))));
			if (dis < distMin) {
				indiceMin = elts(border,j);
				distMin = dis;
			}
		}

		SHL_addEdge(s, elts(l,i), indiceMin);
	}

	LST_delete(border);
}*/

/**
 * @brief Adds a cycle to the Shell with the specified atom ID.
 *
 * This function adds a cycle to the Shell data structure with the specified atom ID.
 *
 * @param s Pointer to the Shell data structure.
 * @param id The ID of the atom to be added to the cycle.
 */
void SHL_addCycle(Shell_t* s, unsigned id) {

	LST_addElement(s->cycle, id);
}

/**
 * @brief Merges two atoms in the Shell.
 *
 * This function merges (combines) two atoms in the Shell data structure.
 * It moves all the edges and properties from the eaten atom to the eater atom,
 * and then removes the eaten atom from the Shell.
 *
 * @param s Pointer to the Shell data structure.
 * @param eater The ID of the atom that will "eat" (absorb) the other atom.
 * @param eaten The ID of the atom that will be "eaten" (absorbed) by the other atom.
 */
void SHL_mergeAtom(Shell_t* s, unsigned eater, unsigned eaten) {

	int i;
	AtomShl_t* a;

	if (eater < size(s) && eaten < size(s) && eater != eaten) {

		a = atom(s, eaten);
		
		for (i=0; i<neighborhoodSize(a); i++) {
			if (neighbor(a, i) != -1 ) {
				SHL_addEdge(s, eater, neighbor(a,i));
			}
		}
		parentAtom(atom(s, eater)) = parentAtom(a);

		if (cycle(s, eaten))
			SHL_addCycle(s, eater);

		if (flag(atom(s, eater)) < flag(atom(s, eaten)))
			flag(atom(s, eater)) = flag(atom(s, eaten));

		SHL_removeAtom(s, eaten);
	}
}

/*void SHL_mergeAtom2(Shell_t* s, unsigned id1, unsigned id2) {

	AtomShl_t* a1, *a2;

	if (id1 != id2 && id1 < size(s) && id2 < size(s)) {

		a1 = atom(s, id1);

		a2 = atom(s, id2);

		while (neighbor(a2,0) != -1){
			SHL_addEdge(s, id1, neighbor(a2,0));
			SHL_removeEdge(s, id2, neighbor(a2,0));
			if (flag(a1) < flag(a2))
				flag(a1) = flag(a2);
		}

		coords(a1) = PT_merge(coords(a1), coords(a2));

		SHL_removeAtom(s, id2);
	}
}*/

/**
 * @brief Creates a new Shell data structure.
 *
 * This function creates and initializes a new Shell data structure.
 *
 * @return A pointer to the newly created Shell data structure.
 */
Shell_t* SHL_create() {

	Shell_t *a = malloc(sizeof(Shell_t));

	a->size = 0;
	a->atoms = NULL;
	a->cycle = LST_create();
	a->bond = GPH_create();

	return a;
}

/**
 * @brief Creates a deep copy of a Shell data structure.
 *
 * This function creates a deep copy of the given Shell data structure.
 *
 * @param s Pointer to the original Shell data structure to be copied.
 * @return A pointer to the newly created deep copy of the Shell data structure.
 */
Shell_t* SHL_copy(Shell_t* s) {

	int i;
	Shell_t* copy = malloc(sizeof(Shell_t));

	size(copy) = size(s);
	copy->atoms = malloc(size(copy)*sizeof(AtomShl_t));
	copy->cycle = LST_copy(s->cycle);
	copy->bond = GPH_copy(s->bond);


	for (i=0; i<size(s); i++) {
		
		flag(atom(copy,i)) = flag(atom(s,i));
		coords(atom(copy,i)) = coords(atom(s,i));
		parentAtom(atom(copy,i)) = parentAtom(atom(s,i));
		neighborhood(atom(copy,i)) = LST_copy(neighborhood(atom(s,i)));
	}

	return copy;
}

/**
 * @brief Creates a trimmed copy of a Shell data structure.
 *
 * This function creates a trimmed copy of the given Shell data structure by removing unused atoms or atoms belonging to the envelope.
 *
 * @param s Pointer to the original Shell data structure containing unused atoms.
 * @return A pointer to the newly created trimmed copy of the Shell data structure.
 */
Shell_t* SHL_copyCageAtoms(Shell_t* s) {

	int notDefsCounter = 0;
	int index = 0;
	Shell_t* copy = malloc(sizeof(Shell_t));
	int* relativeEmptyPositions = malloc(size(s) * sizeof(int));

	// Remove the envelope's atoms (change them to unused).
	for (int j = 0; j < size(s); j++) {
		if (flag(atom(s, j)) == SHELL_F) {
			SHL_removeAtom(s, j);
		}
	}

	size(copy) = SHL_nbAtom(s);
	copy->atoms = malloc(size(copy) * sizeof(AtomShl_t));
	copy->cycle = LST_copy(s->cycle);
	copy->bond = GPH_copy(s->bond);

	// Count the offset in position for each kept atom (to copy the list of neighbors).
	for (int i = 0; i < size(s); i++) {
		if ((flag(atom(s, i)) == NOT_DEF_F)) {
			notDefsCounter++;
		}
		else {
			relativeEmptyPositions[i] = notDefsCounter;
			index++;
		}
	}

	// Copy only used atoms.
	index = 0;
	for(int i = 0; i < size(s); i++) {
		if ((flag(atom(s,i)) != NOT_DEF_F)) {
			flag(atom(copy,index)) = flag(atom(s,i));
			coords(atom(copy,index)) = coords(atom(s,i));
			parentAtom(atom(copy,index)) = parentAtom(atom(s,i));
			neighborhood(atom(copy,index)) = LST_copyWithShift(neighborhood(atom(s,i)), relativeEmptyPositions);
			index++;
		}
	}
	free(relativeEmptyPositions);
	return copy;
}

/*void SHL_testDis(Shell_t* s) {

int i, j;
	for (i = 0; i < size(s); i++) {
		if (flag(atom(s,i)) != NOT_DEF_F && flag(atom(s,i)) != CYCLE_F && flag(atom(s,i)) != SHELL_F) {
			for (j = i + 1; j < size(s); j++) {
				if (flag(atom(s,j)) != NOT_DEF_F && flag(atom(s,j)) != CYCLE_F && flag(atom(s,j)) != SHELL_F) {
					if (dist(coords(atom(s,i)), coords(atom(s,j))) < MINDIS) {
						if (flag(atom(s,i)) == LINKABLE_F && flag(atom(s,j)) == HYDRO_PATTERN_F) {
							SHL_removeAtom(s, i);
						}
						else if (flag(atom(s,i)) == HYDRO_PATTERN_F && flag(atom(s,j)) == LINKABLE_F) {
							SHL_removeAtom(s, j);
						}
						else {
							if (flag(atom(s,i)) == LINKABLE_F && flag(atom(s,j)) == LINKABLE_F) {
								flag(atom(s,i)) = CARBON_F; // Change the flag when both linkable to prevent choosing them as the starting or ending atom.
								SHL_mergeAtom2(s, i, j);

								Point_t hydrogen1 = AX2E2(coords(atom(s,i)), coords(atom(s,neighbor(atom(s,i),0))), coords(atom(s,neighbor(atom(s,i),1))), DIST_ATOM_H);
								Point_t hydrogen2 = AX3E1(coords(atom(s,i)), coords(atom(s,neighbor(atom(s,i),0))), coords(atom(s,neighbor(atom(s,i),1))), hydrogen1, DIST_ATOM_H);
								int idHydrogen = SHL_addAtom(s, hydrogen1, -1);
								flag(atom(s, idHydrogen)) = HYDROGEN_F;
								SHL_addEdge(s, i, idHydrogen);
								idHydrogen = SHL_addAtom(s, hydrogen2, -1);
								flag(atom(s, idHydrogen)) = HYDROGEN_F;
								SHL_addEdge(s, i, idHydrogen);
							}
							else{
								SHL_mergeAtom2(s, i, j);
							}
						}
					}
				}
			}
		}
	}
}*/

/**
 * @brief Deletes an atom in the Shell.
 *
 * This function deletes an atom (node) in the Shell data structure.
 *
 * @param a Pointer to the atom to be deleted.
 */
void SHL_deleteAtom(AtomShl_t* a) {

	LST_delete(neighborhood(a));
}

/**
 * @brief Deletes a Shell data structure.
 *
 * This function deletes the entire Shell data structure, including all atoms and the graph representation.
 *
 * @param s Pointer to the Shell data structure to be deleted.
 */
void SHL_delete(Shell_t* s) {
	
	int i;

	if (s->atoms != NULL) {
		for (i=0; i<size(s); i++)
			SHL_deleteAtom(atom(s,i));
		free(s->atoms);
	}

	if (s->cycle != NULL)
		LST_delete(s->cycle);

	if (s->bond != NULL)
		GPH_delete(s->bond);

	free(s);
}
