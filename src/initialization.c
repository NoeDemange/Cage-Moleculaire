#include "initialization.h"
#include "input.h"
#include "utile.h"
#include "output.h"

/**
 * Calcule les arêtes de la molécule à partir :
 * - des coordonnées des sommets.
 * - des rayons de covalence des atomes.
 *
 * @param    m     Molécule.
 */
void computeEdges(Molecule_t* m) {
  float AB;

  //pour chaque couple de sommets
  for (int i = 0; i < size(m); i++) {
    for (int j = i + 1; j < size(m); j++) {
      //Calcul de la distance entre les deux points.
      //*100 pour adapter la métrique à celle du rayon de covalence théorique.
      AB = dist(coords(atom(m,i)), coords(atom(m,j))) * 100;

      if (AB <= 20 + (radius(atom(m,i)) + radius(atom(m,j)))) {
        MOL_addEdge(m, i, j);
      }
    }
  }
}

/**
 * Compte les arêtes.
 *
 * @param    m     Molécule.
 */
void computeLigands(Molecule_t* m) {

	computeEdges(m);

	for (int i = 0; i < size(m); i++) {
		MOL_nbLigands(atom(m,i));
	}
}

/**
 * Calcule l'angle moyen entre les différentes paires de voisins de l'atome.
 *
 * @param    m     Molécule.
 * @param		id 		Identifiant de l'atome dans la molécule.
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
 * Calcule le nombre de doublets non liants.
 *
 * @param m Molécule.
 */
void computeLonePairs(Molecule_t* m) {

	computeLigands(m);
	MOL_seekCycle(m);

	for (int i = 0; i < size(m); i++) {
		if (ligands(atom(m,i)) != 1) {
			MOL_nbLonePairs(atom(m,i), angleAvg(m,i), -1, cycle(m,i));
		}
	}

	int neighborStericNumber;
	for (int i = 0; i < size(m); i++) {
		neighborStericNumber = -1;
		for (int j = 0; j < ligands(atom(m,i)); j++) {
			if (neighborStericNumber != 3) {
				neighborStericNumber = steric(atom(m,neighbor(atom(m,i),j)));
			}
		}
		MOL_nbLonePairs(atom(m,i), angleAvg(m,i), neighborStericNumber, cycle(m, i));
	}
}

/**
 * Initialise l'ensemble de la molécule.
 *
 * @param 	name	Nom du fichier dans lequel est sauvegardée la molécule.
 * @return  		  Adresse de la Molecule_t.
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
