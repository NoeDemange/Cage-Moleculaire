#include "initialization.h"
#include "input.h"
#include "util.h"
#include "output.h"

/**
 * @file initialization.c
 * @brief This file contains functions for computing edges, ligands, angles, and lone pairs in a molecule.
 */

/**
 * Compute the edges of the molecule from the atoms' coordinates and their covalent radius.
 *
 * @param m Molecule.
 */
void computeEdges(Molecule_t* m) {
  float distBetweenTwoAtoms;

  // For each pair of atoms.
  for (int i = 0; i < size(m); i++) {
    for (int j = i + 1; j < size(m); j++) {
      // mult. by 100 to adapt the metric to the theoretical covalent radius.
      distBetweenTwoAtoms = dist(coords(atom(m,i)), coords(atom(m,j))) * 100;

      if (distBetweenTwoAtoms <= 20 + (radius(atom(m,i)) + radius(atom(m,j)))) {
        MOL_addEdge(m, i, j);
      }
    }
  }
}

/**
 * Count the edges.
 *
 * @param m Molecule.
 */
void computeLigands(Molecule_t* m) {

	computeEdges(m);

	for (int i = 0; i < size(m); i++) {
		MOL_nbLigands(atom(m,i));
	}
}

/**
 * Compute the mean angle between the different pairs of neighbors of the atom.
 *
 * @param m Molecule.
 * @param id Identifier of the atom in the molecule.
 */
float angleAvg(Molecule_t* m, unsigned id) {
	float alpha = 0;
	Atom_t* a = atom(m, id);

	if (ligands(a) == 1)
		return 0;

	for (int i = 0; i < 4 && neighbor(a, i) != -1; i++) {
		for (int j = i + 1; j < 4 && neighbor(a, j) != -1; j++) {
			alpha += angle(coords(a), coords(atom(m,neighbor(a,i))), coords(atom(m, neighbor(a,j))));
		}
	}
	return (2*alpha) / (ligands(a)*(ligands(a)-1));
}

/** 
 * Compute the number of lone pairs.
 *
 * @param m Molecule.
 */
void computeLonePairs(Molecule_t* m) {

	computeLigands(m);
	MOL_seekCycle(m);

	for (int i = 0; i < size(m); i++) {
		if (ligands(atom(m,i)) != 1) {
			MOL_nbLonePairs(atom(m,i), angleAvg(m,i), -1, cycle(m,i));
		}
	}

	int stericNumberOfNeighbor;
	for (int i = 0; i < size(m); i++) {
		stericNumberOfNeighbor = -1;
		for (int j = 0; j < ligands(atom(m,i)); j++) {
			if (stericNumberOfNeighbor != 3) {
				stericNumberOfNeighbor = steric(atom(m,neighbor(atom(m,i),j)));
			}
		}
		MOL_nbLonePairs(atom(m,i), angleAvg(m,i), stericNumberOfNeighbor, cycle(m, i));
	}
}

/**
 * Initialize the whole molecule.
 *
 * @param name File name of the molecule.
 * @return (Molecule_t) adress of the molecule. 
 */
Molecule_t* initMolecule(char* name) {
	
	Molecule_t* m = readInput_xyz(name);

	readCovalence(m);
	computeLonePairs(m);

	MOL_createBond(m);

	//printf("Graphe de dépendance du sustrat.\n");
	//GPH_write(bond(m));
  
  return m;
}
