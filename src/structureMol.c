#include "structure.h"
#include "util.h"

/**
 * @file structureMol.c
 * @brief Molecule Operations
 *
 * This file contains functions related to molecule operations in the application.
 * These operations involve creating, deleting, and managing
 * molecules, atoms, and their properties.
 */

/**************************************/
/* MOLECULE ***************************/
/**************************************/

/**
 * @brief Adds the identifier of a neighbor to an atom.
 *
 * This function adds the identifier of a neighbor to the neighborhood list of an atom.
 * The identifier is added only if it does not already exist in the list.
 *
 * @param a Pointer to the Atom_t representing the atom.
 * @param id Identifier of the neighbor to be added.
 */
void MOL_addNeighbor(Atom_t* a, unsigned id) {

	if (LST_getIndice(neighborhood(a), id) == -1) {
		int indice = LST_getIndiceFree(neighborhood(a));

		if (indice == -1) {
			printf("MOL_addNeighbor : -1 \n");
			exit(1);
		}

		neighbor(a, indice) = id;
	}
}

/**
 * @brief Removes a neighbor from an atom.
 *
 * This function removes a neighbor from the neighborhood list of an atom.
 * If the neighbor with the given identifier does not exist in the list, nothing happens.
 *
 * @param a Pointer to the Atom_t representing the atom.
 * @param id Identifier of the neighbor to be removed.
 */
void MOL_removeNeighbor(Atom_t* a, unsigned id) {

	int indice = LST_getIndice(neighborhood(a), id);

	if (indice != -1)
		neighbor(a, indice) = -1;
}

/**
 * @brief Calculates the number of ligands (doublets liants) of an atom.
 *
 * The number of ligands of an atom is equal to the number of neighbors it has in the molecule.
 *
 * @param a Pointer to the Atom_t representing the atom.
 */
void MOL_nbLigands(Atom_t* a) {
	
	int cpt = 0;
	for (int i = 0; i < 4; i++) {
		if (neighbor(a,i) != -1) {
			cpt++;
		}
	}
	ligands(a) = cpt;
}

/**
 * @brief Calculates the number of lone pairs (doublets non liants) of an atom.
 *
 * This function calculates the number of lone pairs of an atom based on its properties,
 * neighboring atoms' properties, and whether the atom belongs to a cycle.
 *
 * @param a Pointer to the Atom_t representing the atom.
 * @param alpha Average angle formed by the neighbors of the atom.
 * @param stericNeighbor Number of doublets of the neighbor of the atom (used when it has only one neighbor).
 * @param cycle Boolean indicating if the atom belongs to a cycle.
 */
void MOL_nbLonePairs(Atom_t* a, float alpha, int stericNeighbor, unsigned cycle) {

	if (ligands(a) == 1) {
		if (!strcmp(symbol(a), "H")) {
			lonePairs(a) = 1;
		}
		else if(!strcmp(symbol(a), "Cl") || !strcmp(symbol(a), "Br")
			|| !strcmp(symbol(a), "F") || !strcmp(symbol(a), "I")) {
			lonePairs(a) = 3;
		}
		else {
				lonePairs(a) = stericNeighbor - 1;
		}
	}
	else {
		if (ligands(a) == 4)
			lonePairs(a) = 0;
		
		else if (abs(120-alpha) < 4)
			lonePairs(a) = 3 - ligands(a);

		else if (alpha-109 < 7)
			if (cycle == 1 && stericNeighbor == 3)
				lonePairs(a) = 3 - ligands(a);
			else if (cycle == 1 && stericNeighbor == -1)
				lonePairs(a) = -1;
			else
				lonePairs(a) = 4 - ligands(a);
		else 
			lonePairs(a) = 0;
	}
}

/**
 * @brief Finds all vertices belonging to a cycle in a molecule.
 *
 * This function finds all vertices (atoms) belonging to a cycle in a molecule
 * and updates the molecule's cycle list accordingly.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 */
void MOL_seekCycle(Molecule_t* m) {

	Graph_t* g = MolToGph(m);
	m->cycle = GPH_seekCycle(g);
	GPH_delete(g);
}

/**
 * @brief Calculates the number of edges (bonds) in a molecule.
 *
 * The number of edges (bonds) in a molecule is equal to half of the total number of ligands
 * (doublets liants) in all atoms of the molecule.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 * @return The number of edges (bonds) in the molecule.
 */
int MOL_nbEdges(Molecule_t* m) {
	int cpt = 0;

	for (int i = 0; i < size(m); i++) {
		cpt += ligands(atom(m,i));
	}
	return cpt/2;
}

/**
 * @brief Adds an edge (bond) between two vertices (atoms) in a molecule.
 *
 * This function adds a bond between two vertices (atoms) in a molecule.
 * The edge is added by updating the neighborhood lists of the two atoms to include each other.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 * @param id1 Identifier of the first vertex (atom).
 * @param id2 Identifier of the second vertex (atom).
 */
void MOL_addEdge(Molecule_t* m, unsigned id1, unsigned id2) {
	if (id1 != id2) {
		MOL_addNeighbor(atom(m,id1), id2);
		MOL_addNeighbor(atom(m,id2), id1);
	}
}

/**
 * @brief Removes an edge (bond) between two vertices (atoms) in a molecule.
 *
 * This function removes a bond between two vertices (atoms) in a molecule.
 * The edge is removed by updating the neighborhood lists of the two atoms to exclude each other.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 * @param id1 Identifier of the first vertex (atom).
 * @param id2 Identifier of the second vertex (atom).
 */
void MOL_removeEdge(Molecule_t* m, unsigned id1, unsigned id2) {
	
	MOL_removeNeighbor(atom(m,id1), id2);
	MOL_removeNeighbor(atom(m,id2), id1);
}

/**
 * @brief Initializes a new atom with default values.
 *
 * This function initializes a new atom by setting its symbol, radius, number of ligands,
 * number of lone pairs, coordinates, and neighborhood list to default values.
 *
 * @param a Pointer to the Atom_t representing the atom to be initialized.
 */
void MOL_createAtom(Atom_t* a) {

	symbol(a)[0] = '\0';
	radius(a) = -1;
	ligands(a) = -1;
	lonePairs(a) = -1;

	atomX(a) = 0;
	atomY(a) = 0;
	atomZ(a) = 0;

	neighborhood(a) = LST_create();
}

/**
 * @brief Creates the dependency graph of a molecule.
 *
 * The dependency graph groups atoms in the molecule that can participate in hydrogen bonds.
 * An edge exists between two atoms if they cannot participate in hydrogen bonds simultaneously.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 */
void MOL_createBond(Molecule_t* m) {

	int i, j, k;
	Atom_t* a;
	List_t* lh;

	for (i=0; i<size(m); i++) {
		a = atom(m, i);

		if (!strcmp(symbol(a), "O") || !strcmp(symbol(a), "N")
			||  !strcmp(symbol(a), "F")
			/*&& (ligands(a) != 4 && lonePairs(a) != 0)*/) {
			lh = LST_create();
			if (lonePairs(a) > 0 /*sauf (2,2) && (1,3)*/) {
				LST_addElement(lh, i);
				GPH_addVertex(bond(m), i);
			}

			for (j=0; j<neighborhoodSize(a) && neighbor(a,j) != -1; j++)
				if (!strcmp(symbol(atom(m, neighbor(a,j))), "H")) {
					LST_addElement(lh, neighbor(a,j));
					GPH_addVertex(bond(m), neighbor(a,j));
				}

			for (j=0; j<size(lh) && elts(lh,j)!=-1; j++)
				for (k=j+1; k<size(lh) && elts(lh,k)!=-1; k++)
					GPH_addEdge(bond(m), elts(lh,j), elts(lh,k));

			LST_delete(lh);
		}
	}
}

/**
 * @brief Finds the normal vector of a vertex (atom) in a molecule.
 *
 * This function finds the normal vector of a vertex (atom) in a molecule based on its neighbors.
 * It recursively calculates the normal vector for the vertex if it has only one neighbor.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 * @param ida Identifier of the vertex (atom) for which to find the normal vector.
 * @param dad Identifier of the parent vertex (atom) in the recursive search (used for atoms with one neighbor).
 * @return The normal vector of the vertex (atom) in the molecule.
 */
Point_t MOL_seekNormal(Molecule_t* m, unsigned ida, unsigned dad){

	Atom_t* a = atom(m,ida);

	if (ligands(a) == 1) {
		return MOL_seekNormal(m, neighbor(a,0), ida);
	}

	if (steric(a) == 2) {
		if (neighbor(a,0) != dad) {
			return MOL_seekNormal(m, neighbor(a,0), ida);
		}
		else {
			return MOL_seekNormal(m, neighbor(a,1), ida);
		}
	}

	return planNormal(coords(a),
					coords(atom(m,neighbor(a,0))),
					coords(atom(m,neighbor(a,1))));
}

/**
 * @brief Creates and initializes a new molecule.
 *
 * This function creates and initializes a new molecule with the specified number of vertices (atoms).
 * The atoms are created with default values, and the molecule's cycle and bond properties are set to NULL.
 *
 * @param size Number of vertices (atoms) in the molecule.
 * @return Pointer to the newly created molecule.
 */
Molecule_t* MOL_create(unsigned size) {

	Molecule_t *m = malloc(sizeof(Molecule_t));

	size(m) = size;
	m->atoms = malloc(size*sizeof(Atom_t));

	for (int i = 0; i < size; i++) {
		MOL_createAtom(atom(m,i));
	}

	m->cycle = NULL;
	m->bond = GPH_create();

	return m;
}

/**
 * @brief Deletes an atom from a molecule.
 *
 * Deleting an atom involves freeing the memory used by its neighborhood list.
 *
 * @param a Atom to be deleted.
 */
void MOC_deleteAtom(Atom_t* a) {

	LST_delete(neighborhood(a));
}

/**
 * @brief Deletes a molecule and its associated elements.
 *
 * Deleting a molecule involves freeing the memory used by its atoms, cycle list, and bond graph.
 *
 * @param m Pointer to the Molecule_t representing the molecule to be deleted.
 */
void MOL_delete(Molecule_t* m) {

	for (int i = 0; i < size(m); i++) {
		MOC_deleteAtom(atom(m,i));
	}

	free(m->atoms);
	LST_delete(m->cycle);
	GPH_delete(m->bond);
	free(m);
}