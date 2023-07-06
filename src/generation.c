#include "generation.h"
#include "util.h"
#include "output.h"


/**
 * @brief Creation of the beginning of cages according to the dependencies.
 * 
 * @param m Grouping of the main structures (substrate and envelope).
 */
void generateDependancies(Main_t* m) {

	unsigned i, j, copy;
	Shell_t* mo;
	Vertex_t* v;

	//Création des différents moc en fonction des dépendances.
	//mocSize change au fur et à mesure des itérations
	for (i=0; i< /*mocSize(m)*/1; i++) { //Pour tous les mocs
		mo = moc(m,i);
		if (size(moc(m,i)) != 0) {
			for (j=0; j<size(bond(moc(m,i))); j++) {//Pour tous les sommets
				v = vertex(bond(moc(m,i)),j);

				if (id(v) != -1 && nbNeighbors(v) > 0) {

					while (nbNeighbors(v) != 0) {
						copy = MN_copyMoc(m, mo);
						//printf("copy %d \n", copy);
						//printf("number copy %d\n", copy);
						SHL_removeVertex(moc(m, copy), id(v));
						//printf("neighbor(v,0) = %d\n", neighbor(v,0));
						SHL_removeVertex(mo, neighbor(v,0));
					}
				}
			}
		}
		
	}
	//printf("mocSize %d\n", mocSize(m));
}

/*void checkInsertVertex(Shell_t* m, List_t* l, unsigned idv) {

	int i, index = idv;
	float min, distance;
	AtomShl_t *v = atom(m,idv), *s;

	if (size(l) > 0 && elts(l,0) != -1) {
		min = dist(coords(v), coords(atom(m, elts(l,0))));
		index = 0;
	}

	for (i = 1; forEachElement(l, i); i++) {
		distance = dist(coords(v), coords(atom(m, elts(l,i))));
		if (min > distance) {
			min = distance;
			index = i;
		}
	}

	if (index != idv && min < MINDIS) {
		s = atom(m, index);

		while (neighbor(s,0) != -1){
			LST_addElement(l,neighbor(s,0));
			SHL_removeEdge(m, idv, neighbor(s,0));
		}

		flag(v) = flag(s);
		SHL_removeVertex(m, idv);
	}

}*/

/**
 * @brief Insert an acceptor hydrogen pattern with a triangular geometry.
 * 
 * @param m Envelope with the beginning of the cage.
 * @param idcenter Index of the heteroatom involved in the H-bond.
 * @param normal Normal vector starting point.
 * @param dir Direction (terminal point) of the normal vector.
 */
void insertAcceptor1(Shell_t* m, unsigned idcenter, Point_t normal, Point_t dir) {
	//La position de v est la position du premier.
	int index;
	Point_t x1,x2;
	AtomShl_t* center = atom(m,idcenter);
	List_t* l = LST_create();

	flag(center) = HYDRO_PATTERN_F;
	dir = subPoint(initPoint(0), normalization(dir, DIST_SIMPLE));

	// Remove the edges between center and its neighbors and add them to the list.
	while (neighbor(center,0) != -1) {
		LST_addElement(l,neighbor(center,0));
		SHL_removeEdge(m, idcenter, neighbor(center,0));
	}

	//Position du deuxième
	x1 = addPoint(coords(center), rotation(normal, 120, dir));
	//Position du troisième
	x2 = addPoint(coords(center), rotation(normal, -120, dir));

	for (int i = 0; i < size(m); i++) {
		if ((flag(atom(m, i)) != NOT_DEF_F) && (flag(atom(m,i)) != SHELL_F) 
			&& ((dist(coords(atom(m, i)),x1) < DIST_GAP_CAGE)|| (dist(coords(atom(m, i)),x2) < DIST_GAP_CAGE))) {
			flag(center) = SHELL_F;
			return;
		}
	}

	//Ajout du deuxième
	index = SHL_addAtom(m, x1, -1);
	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idcenter, index);
	if (flag(atom(m,index)) < LINKABLE_F)
		flag(atom(m,index)) = LINKABLE_F;
	
	//Ajout du troisième
	index = SHL_addAtom(m, x2, -1);
	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idcenter, index);
	if (flag(atom(m,index)) < LINKABLE_F)
		flag(atom(m,index)) = LINKABLE_F;

	//Rattacher les nouveaux sommets à ceux de la liste l.
	//SHL_linkBorder(m, idcenter, l);

	LST_delete(l);
}

/**
 * @brief Insert a donor hydrogen pattern with a triangular geometry.  
 * 
 * @param m Envelope with the beginning of the cage.
 * @param idhydro Index of the hydrogen involved in the H-bond.
 * @param normal Normal vector starting point.
 * @param dir Direction (terminal point) of the normal vector.
 */
void insertDonor1(Shell_t* m, unsigned idhydro, Point_t normal, Point_t dir) {
	//La position de v est la position du premier.
	int index, idc;
	Point_t center,x1,x2;
	AtomShl_t* hydro = atom(m, idhydro);
	List_t* neighborsFirstAtom = LST_create();

	flag(hydro) = HYDRO_PATTERN_F;

	// Remove the edges between hydro and its neighbors and add them to the list.
	while (neighbor(hydro,0) != -1) {
		LST_addElement(neighborsFirstAtom, neighbor(hydro,0));
		SHL_removeEdge(m, idhydro, neighbor(hydro,0));
	}

	//Position du deuxième sommet : centre du motif
	//Hydrogène+taille d'une liaison simple moyenne.
	dir = normalization(dir, DIST_ATOM_H);
	center = addPoint(coords(hydro), dir);
	//Position du deuxième :
	dir = subPoint(initPoint(0), normalization(dir, DIST_SIMPLE));
	x1 = addPoint(center, rotation(normal, 120, dir));
	//Position du troisième :
	x2 = addPoint(center, rotation(normal, -120, dir));


	for (int i = 0; i < size(m); i++) {
		if ((flag(atom(m, i)) != NOT_DEF_F) && (flag(atom(m,i)) != SHELL_F) 
			&& ((dist(coords(atom(m, i)),center) < DIST_GAP_CAGE)||(dist(coords(atom(m, i)),x1) < DIST_GAP_CAGE)||(dist(coords(atom(m, i)),x2) < DIST_GAP_CAGE))) {
			flag(hydro) = SHELL_F;
			return;
		}
	}

	//Ajout du centre
	idc = SHL_addAtom(m, center, -1);
	//checkInsertVertex(m, l, idc);
	SHL_addEdge(m, idhydro, idc);
	flag(atom(m,idc)) = HYDRO_PATTERN_F;

	//Ajout du deuxième
	index = SHL_addAtom(m, x1, -1);

	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idc, index);
	if (flag(atom(m,index)) < LINKABLE_F)
		flag(atom(m,index)) = LINKABLE_F;

	//Ajout du troisième
	index = SHL_addAtom(m, x2, -1);
	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idc, index);
	if (flag(atom(m,index)) < LINKABLE_F)
		flag(atom(m,index)) = LINKABLE_F;

	//Rattacher les nouvaux sommets à ceux de la liste l.
	//SHL_linkBorder(m, idc, l);
	
	LST_delete(neighborsFirstAtom);
}

/**
 * @brief Insert a donor hydrogen pattern with a tetrahedral geometry.  
 * 
 * @param m Envelope with the beginning of the cage.
 * @param idhydro Index of the hydrogen involved in the H-bond.
 * @param normal Normal vector starting point.
 * @param dir Direction (terminal point) of the normal vector.
 */
void insertDonor2(Shell_t* m, unsigned idhydro, Point_t normal, Point_t dir) {
	//La position de v est la position du premier. Premier sommet du tétraèdre
	int index, idc;
	Point_t x1, center, x2, x3, x4;
	AtomShl_t* hydro = atom(m, idhydro), *c;
	List_t* l = LST_create();

	flag(hydro) = HYDRO_PATTERN_F;

	//Insérer tous les voisins de hydro dans la liste.
	while (neighbor(hydro,0) != -1) {
		LST_addElement(l,neighbor(hydro,0));
		SHL_removeEdge(m, idhydro, neighbor(hydro,0));
	}

	x1 = coords(hydro);

	//Position du centre du motif
	//Hydrogène+taille d'une liaison simple moyenne.
	dir = normalization(dir, DIST_ATOM_H);
	center = addPoint(x1, dir);
	//Deuxième sommet du tétraèdre
	normal = planNormal(center, x1, normal);
	normal = normalization(normal,1);
	x2 = AX1E3(center, x1, normal, DIST_SIMPLE);
	//Troisième sommet du tétraèdre
	x3 = AX2E2(center, x1, x2, DIST_SIMPLE);
	//Quatrième sommet du tétraèdre
	x4 = AX3E1(center, x1, x2, x3, DIST_SIMPLE);
	for (int l = 0; l < size(m); l++) {
		if((flag(atom(m, l)) != NOT_DEF_F) && (flag(atom(m,l)) != SHELL_F) 
		&& ((dist(coords(atom(m,l)),center) < DIST_GAP_CAGE) || (dist(coords(atom(m,l)),x2) < DIST_GAP_CAGE)
		|| (dist(coords(atom(m,l)),x3) < DIST_GAP_CAGE) || (dist(coords(atom(m,l)),x4) < DIST_GAP_CAGE))){
			flag(hydro) = SHELL_F;
			return;
		}
	}
	//Ajout centre du motif
	idc = SHL_addAtom(m, center, -1);
	//checkInsertVertex(m, l, idc);
	SHL_addEdge(m, idhydro, idc);
	c = atom(m,idc);
	flag(c) = HYDRO_PATTERN_F;

	//Ajout deuxième sommet du tétraèdre
	index = SHL_addAtom(m, x2, -1);
	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idc, index);
	flag(atom(m,index)) = LINKABLE_F;

	//Ajout Troisième sommet du tétraèdre
	index = SHL_addAtom(m, x3, -1);
	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idc, index);
	flag(atom(m,index)) = LINKABLE_F;

	//Ajout quatrième sommet du tétraèdre
	index = SHL_addAtom(m, x4, -1);
	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idc, index);
	flag(atom(m,index)) = LINKABLE_F;

	//Rattacher les nouvaux sommets à ceux de la liste l.
	//SHL_linkBorder(m, idc, l);
	
	LST_delete(l);
}

/**
 * @brief Add donor or acceptor hydrogen patterns to the beginning of the cage.
 * 
 * @param m Grouping of the main structures (substrate and envelope).
 */
void generateHydrogenPattern(Main_t* m) {

	AtomShl_t *atomCage;
	Atom_t *parentAtom; // Corresponding atom in the substrate.

	//TODO choose to keep the loop to make the dependencies graph
	for (int i = 0; i < /*mocSize(m)*/1; i++) {
		for (int j = 0; j < size(bond(moc(m, i))); j++) {
			// Get the atom id from the dependency graph
			int idAtomCage = id(vertex(bond(moc(m, i)), j));

			if (idAtomCage != -1) {
				atomCage = atom(moc(m,i), idAtomCage);
				int idParentAtom = parentAtom(atomCage);
				int tooClose = 0;
				for (int k = 0; k < size(moc(m,i)); k++) {
					if((k != idAtomCage) && (flag(atom(moc(m, i), k)) != NOT_DEF_F) && (flag(atom(moc(m, i), k)) != SHELL_F) && (dist(coords(atom(moc(m, i), k)), coords(atomCage)) < DIST_GAP_CAGE)) {
						tooClose = 1;
						break;
					}
				}
				if(!tooClose) {
					parentAtom = atom(substrat(m), idParentAtom);
					if (!strcmp(symbol(parentAtom), "H")) {
						insertAcceptor1(moc(m,i), idAtomCage, MOL_seekNormal(substrat(m), idParentAtom, -1), 
							vector(coords(parentAtom), coords(atomCage)));
					}
					else {
						int haveTriangularGeometry = (steric(parentAtom) == 3);
						if (haveTriangularGeometry) {
							insertDonor1(moc(m,i), idAtomCage, MOL_seekNormal(substrat(m), idParentAtom, -1), 
								vector(coords(parentAtom), coords(atomCage)));
						}
						else {
							insertDonor2(moc(m,i), idAtomCage, MOL_seekNormal(substrat(m), idParentAtom, -1), 
								vector(coords(parentAtom), coords(atomCage)));
						}
					}
				}
			}
		}
		//SHL_testDis(moc(m,i));
	}
}

/**
 * @brief Add aromatic rings to the envelope.
 * 
 * The function find for each atom involved in a cycle
 * if it can be part of a triangular pattern,
 * then add the pattern to the envelpe.
 * 
 * @param s Envelope of the substrate.
 */
void generateCycle(Shell_t* s) {
	AtomShl_t* atom;
	List_t* neighborsNotInCycle;
	List_t* atomsInCycle = LST_create();

	// Find the atoms of the shell involved in a cycle.
	for (int i = 0; i < size(s); i++) {
		if (flag(atom(s, i)) != NOT_DEF_F) {
			if (cycle(s, i)) {	
				LST_addElement(atomsInCycle, i);
			}
		}
	}

	for (int i = 0; forEachElement(atomsInCycle, i); i++) {
		atom = atom(s, elts(atomsInCycle,i));
		neighborsNotInCycle = LST_create();

		// Find the neighbors of the atom not involved in a cycle.
		for (int j = 0; forEachNeighbor(atom, j); j++) {
			if (!LST_check(atomsInCycle, neighbor(atom, j)) ||
			dist(coords(atom), coords(atom(s, neighbor(atom, j)))) > MAXDIS_CYCLE) {
				LST_addElement(neighborsNotInCycle,neighbor(atom, j));
			}
		}

		int haveEnoughNeighborsInCycle = (SHL_nbNeighborhood(atom) - LST_nbElements(neighborsNotInCycle) >= 2);
		if (haveEnoughNeighborsInCycle) {
			// Remove old link between the atom (i) and its neighbors (j).
			for (int j = 0; forEachElement(neighborsNotInCycle, j); j++) {
				SHL_removeEdge(s, elts(atomsInCycle,i), elts(neighborsNotInCycle,j));
			}
			flag(atom) = CYCLE_F;
			if (SHL_nbNeighborhood(atom) == 2) {
				int idNewNeigbor = -1;
				Point_t newNeighborPoint = addThirdPoint(coords(atom), coords(atom(s, neighbor(atom, 0))),
							coords(atom(s, neighbor(atom, 1))), SIMPLE_CYCLE);
				idNewNeigbor = SHL_addAtom(s, newNeighborPoint, -1);

				for (int j = 0; j < size(s); j++) {
					if (flag(atom(s, j)) != NOT_DEF_F && dist(newNeighborPoint, coords(atom(s, j))) < MINDIS_CYCLE) {
						if (LST_check(neighborsNotInCycle, j)) {
							LST_removeElement(neighborsNotInCycle, j);
						}
						if (cycle(s, j)) {
							//TODO remove this condition unless needed
							LST_addElement(atomsInCycle, idNewNeigbor);
						}
						if (flag(atom(s, j)) != SHELL_F) {
							SHL_mergeAtom(s, idNewNeigbor, j);
						}
					}
				}
				if (flag(atom(s, idNewNeigbor)) <= LINKABLE_F) {
					flag(atom(s, idNewNeigbor)) = LINKABLE_F;
				}
				SHL_addEdge(s, elts(atomsInCycle, i), idNewNeigbor);
			}
			//SHL_linkBorder(s, elts(atomT,i), nei);
		}
		LST_delete(neighborsNotInCycle);
	}
	LST_delete(atomsInCycle);
}

void generatePathlessCages(Main_t* m) {

	envarom(m) = SHL_copy(envelope(m));
	//generateDependancies(m);
	printf("### Aromatic rings generation ###\n");
	generateCycle(envarom(m));
	//SHL_write(envarom(m));

	printf("### Hydrogen patterns generation ###\n");
	MN_copyMoc(m, envarom(m));
	generateDependancies(m);
	generateHydrogenPattern(m);

}
