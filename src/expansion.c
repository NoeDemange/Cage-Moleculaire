#include "expansion.h"
#include "util.h"
#include "interface.h"

#include "output.h"

/**
 * @file expansion.c
 * @brief This file contains functions related to the expansion of atomic steric groupings and the generation of the molecule's envelope.
 */

/**
 * Expansion of atomic steric grouping 2. 
 * Add a point in the envelope.
 * 
 * @param m Molecule (input).
 * @param s Envelope (output).
 * @param id Identifier of the processed atom in the molecule.
 */
void expansion_AX1E1(Molecule_t* m, Shell_t* s, unsigned id) {

	unsigned indice;
	Point_t newCoords;
	Atom_t* a = atom(m,id), *x1 = atom(m,neighbor(a,0));

  newCoords = AX1E1(coords(a), coords(x1), DIST_HYDRO);

  indice = SHL_addAtom(s, newCoords, id);

  if (checkVertex(m,id))
			SHL_addVertex(s, indice);
}

/**
 * @brief Expansion of atomic steric grouping 3. 
 * Add one to four points in the envelope.
 *
 * @param m Molecule (input).
 * @param s Envelope (output).
 * @param id Identifier of the processed atom in the molecule.
 */
void expansion_steric3(Molecule_t* m, Shell_t* s, unsigned id) {

	unsigned indice;
	Point_t a, x1, x2, x3, normal;
	a = coords(atom(m,id));
	x1 = coords(atom(m,neighbor(atom(m,id),0)));

	//Ajout d'un point si le groupement possède moins de deux doublets.
	//X2 devient soit le point ajouté à l'enveloppe, soit le second voisin dans la molécule.
	if (ligands(atom(m,id)) < 2) {

		if (neighbor(atom(m,neighbor(atom(m,id),0)),0) == id)
			x2 = coords(atom(m,neighbor(atom(m,neighbor(atom(m,id),0)),1)));
		else
			x2 = coords(atom(m,neighbor(atom(m,neighbor(atom(m,id),0)),0)));

		normal = planNormal(a, x1, x2);
		x2 = AX1E2(a, x1, normal, DIST_HYDRO);
		indice = SHL_addAtom(s, x2, id);

		if (checkVertex(m,id))
			SHL_addVertex(s, indice);

	}
	else
		x2 = coords(atom(m,neighbor(atom(m,id),1)));

	//Ajout d'un point si le groupement possède moins de trois doublets.
	//X3 devient soit le point ajouté à l'enveloppe, soit le troisième voisin dans la molécule.
	if (ligands(atom(m,id)) < 3) {
		x3 = AX2E1(a, x1, x2, DIST_HYDRO);
		indice = SHL_addAtom(s, x3, id);

		if (checkVertex(m,id))
			SHL_addVertex(s, indice);
	}
	else
		x3 = coords(atom(m,neighbor(atom(m,id),2)));

	//Ajout des deux extensions perpendiculaires.
	normal = normalization(planNormal(x1, x2, x3), DIST_HYDRO);

	if (cycle(m,id)) {
		SHL_addCycle(s, SHL_addAtom(s, PT_add(a, normal), id));
		SHL_addCycle(s, SHL_addAtom(s, PT_sub(a, normal), id));
	}
	else {
		SHL_addAtom(s, PT_add(a, normal), id);
		SHL_addAtom(s, PT_sub(a, normal), id);
	}
}

/**
 * @brief Expansion of atomic steric grouping 4. 
 * Add zero to three points in the envelope.
 *
 * @param m Molecule (input).
 * @param s Envelope (output).
 * @param id Identifier of the processed atom in the molecule.
 */
void expansion_steric4(Molecule_t* m, Shell_t* s, unsigned id) {

	unsigned indice;
	Point_t a, x1, x2, x3, normal;

	a = coords(atom(m,id));
	x1  = coords(atom(m, neighbor(atom(m,id), 0)));

	//Ajout d'un point si le groupement possède moins de deux doublets.
	//X2 devient soit le point ajouté à l'enveloppe, soit le second voisin dans la molécule.
	if (ligands(atom(m,id)) < 2) {

		if (neighbor(atom(m,neighbor(atom(m,id),0)),0) == id)
			x2 = coords(atom(m,neighbor(atom(m,neighbor(atom(m,id),0)),1)));
		else
			x2 = coords(atom(m,neighbor(atom(m,neighbor(atom(m,id),0)),0)));

		//à voir si on le garde
		x2 = PT_add(a, planNormal(a, x1, x2));
		normal = planNormal(a, x1, x2);
		x2 = AX1E3(a, x1, normal, DIST_HYDRO);
		indice = SHL_addAtom(s, x2, id);

		if (checkVertex(m,id))
			SHL_addVertex(s, indice);
	}
	else
		x2 = coords(atom(m, neighbor(atom(m,id), 1)));

	//Ajout d'un point si le groupement possède moins de trois doublets.
	//X3 devient soit le point ajouté à l'enveloppe, soit le troisième voisin dans la molécule.
	if (ligands(atom(m,id)) < 3) {

		x3 = AX2E2(a, x1, x2, DIST_HYDRO);
		indice = SHL_addAtom(s, x3, id);

		if (checkVertex(m,id))
			SHL_addVertex(s, indice);
	}
	else
		x3 = coords(atom(m, neighbor(atom(m,id), 2)));

	if (ligands(atom(m,id)) < 4) {

		indice = SHL_addAtom(s, AX3E1(a, x1, x2, x3, DIST_HYDRO), id);

		if (checkVertex(m,id))
			SHL_addVertex(s, indice);
	}
}

/********************* A refaire ************************/
/**
 * Expansion of atomic group with an angle of 180°.
 * Add four points in the envelope.
 * 
 * @param m Molecule (input).
 * @param s Envelope (output).
 * @param id Identifier of the processed atom in the molecule.
 */
void expansion_AX2E0(Molecule_t* m, Shell_t* s, unsigned id) {
	Point_t normal, newCoords;
	Atom_t* a = atom(m,id),	*x1 = atom(m,neighbor(a,0)),
	*x2 = atom(m, neighbor(a,1));

	normal = normalization(vector(coords(a), coords(x1)),1);

	newCoords = normalization(
		planNormal(coords(a), coords(x1), coords(x2)), DIST_HYDRO);

	SHL_addAtom(s, PT_add(coords(a), newCoords), id);

	SHL_addAtom(s, PT_add(coords(a), normalization(
		rotation(normal, 90, newCoords), DIST_HYDRO)), id);

	SHL_addAtom(s, PT_add(coords(a), normalization(
		rotation(normal, 180, newCoords), DIST_HYDRO)), id);

	SHL_addAtom(s, PT_add(coords(a), normalization(
		rotation(normal, -90, newCoords), DIST_HYDRO)), id);
}

/**
 * Expansion of the whole molecule.
 * Add zero to three points in the envelope.
 * 
 * @param m Molecule (input) already instantiated.
 * @param s Envelope (output) created.
 */
void expansion(Molecule_t* m, Shell_t* s) {
	int i, j;

	//Création du nuage de points constituant l'enveloppe.
	/*rajouter une liste pour être sur qu'ils soient fait dans le bon ordre*/
	/*vérifier que toutes les donnés ont été calculé*/
	for (i = 0; i < size(m); i++) {
		if (ligands(atom(m,i)) == 1 && lonePairs(atom(m,i)) == 1)
			expansion_AX1E1(m, s, i);
		else if (steric(atom(m,i)) == 3)
			expansion_steric3(m, s, i);
		else if (steric(atom(m,i)) == 4)
			expansion_steric4(m, s, i);
		else if (ligands(atom(m,i)) == 2 && lonePairs(atom(m,i)) == 0)
			expansion_AX2E0(m, s, i);
	}

	//Graph edges construction.
	for (i = 0; i < size(bond(s)) && id(vertex(bond(s),i)) != -1; i++)
		for (j = i + 1; j < size(bond(s)) && id(vertex(bond(s),j)) != -1; j++) {
			if (parentAtom(atom(s,id(vertex(bond(s),i)))) == parentAtom(atom(s,id(vertex(bond(s),j))))
					|| checkBond(m,
						parentAtom(atom(s,id(vertex(bond(s),i)))), 
						parentAtom(atom(s,id(vertex(bond(s),j)))))
					)
				SHL_addBond(s, id(vertex(bond(s),i)), id(vertex(bond(s),j)));
		}

}

/**
 * Construction of the edges of the envelope.
 * Call to R.
 * 
 * @param s				Envelope that already has a cloud of points.
 * @param alpha 	Paramètre de la sphère. (3 est souvent une bonne mesure, 4 sinon)
 */
void alphaShape(Shell_t* s, double alpha) {
	int i;
	Ashape_t* as3d = Cashape3d(s, alpha);

	for (i = 0; i < (as3d->nb_edge/2); i++) {
		SHL_addEdge(s, as3d->edge[i]-1, as3d->edge[i+as3d->nb_edge/2]-1);
	}

	for (i = 0; i < size(s); i++) {
		if (neighborhoodSize(atom(s,i)) == 0) {
			SHL_removeAtom(s, i);
			GPH_removeVertex(bond(s),i);
		}
	}
	ASP_delete(as3d);
}

/**
 *  Creation of the whole envelope from the molecule.
 *
 * @param m Molecule (input) already instantiated.
 * @param s Envelope (output) created.
 */
Shell_t* createShell(Molecule_t* m, double alpha) {

	printf("\n####### Start of the envelope and binding patterns generation #######\n");
	Shell_t* s = SHL_create();
	expansion(m, s);
	//SHL_writeMol2("../results/vec.mol2", s);
	alphaShape(s, alpha);
	//printf("Graphe de dépendance de l'enveloppe.\n");
	//GPH_write(bond(s));	

	return s;
}
