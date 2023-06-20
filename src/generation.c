#include "generation.h"
#include "utile.h"
#include "output.h"


/**
 * @brief Creation of the beginning of cages according to the dependencies.
 * 
 * @param m Gathering of the main elements (substrate and envelope).
 */
void generationDep(Main_t* m) {

	unsigned i, j, copy;
	Shell_t* mo;
	Vertex_t* v;

	//Création des différents moc en fonction des dépendances.
	//mocSize change au fur et à mesure des itérations
	for (i=0; i</*mocSize(m)*/1	; i++) { //Pour tous les mocs
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
 * @brief Insert a hydrogen donor pattern with a triangular geometry.
 * 
 * @param m Envelope with the beginning of the cage.
 * @param idv Index of the hydrogen atom.
 * @param normal Normal vector starting point.
 * @param dir Direction (terminal point) of the normal vector.
 */
void insertDonor1(Shell_t* m, unsigned idv, Point_t normal, Point_t dir) {
	//La position de v est la position du premier.
	int index;
	Point_t new_coords;
	AtomShl_t* v = atom(m,idv);
	List_t* l = LST_create();

	flag(v) = HYDRO_BOND_F;
	dir = subPoint(initPoint(0), normalization(dir, DIST_SIMPLE));

	// Remove the edges between v and its neighbors and add them to the list.
	while (neighbor(v,0) != -1) {
		LST_addElement(l,neighbor(v,0));
		SHL_removeEdge(m, idv, neighbor(v,0));
	}

	//Position du deuxième : (rotation de normal, 120, -dir) + coords(v)
	new_coords = addPoint(coords(v), rotation(normal, 120, dir));
	index = SHL_addAtom(m, new_coords, -1);

	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idv, index);
	if (flag(atom(m,index)) < LINKABLE_F)
		flag(atom(m,index)) = LINKABLE_F;
	
	v = atom(m,idv);
	
	//Position du troisième : (rotation de normal, -120, -dir) + coords(v)
	new_coords = addPoint(coords(v), rotation(normal, -120, dir));
	index = SHL_addAtom(m, new_coords, -1);

	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idv, index);
	if (flag(atom(m,index)) < LINKABLE_F)
		flag(atom(m,index)) = LINKABLE_F;

	//Rattacher les nouveaux sommets à ceux de la liste l.
	//SHL_linkBorder(m, idv, l);

	LST_delete(l);
}

/**
 * @brief Insert a hydrogen acceptor pattern with a triangular geometry.  
 * 
 * @param m Envelope with the beginning of the cage.
 * @param idv Index of the heteroatom.
 * @param normal Normal vector starting point.
 * @param dir Direction (terminal point) of the normal vector.
 */
void insertAcceptor1(Shell_t* m, unsigned idv, Point_t normal, Point_t dir) {
	//La position de v est la position du premier.
	int index, idc;
	Point_t new_coords;
	AtomShl_t* v = atom(m, idv), *c;
	List_t* neighborsFirstAtom = LST_create();

	flag(v) = HYDRO_BOND_F;

	// Remove the edges between v and its neighbors and add them to the list.
	while (neighbor(v,0) != -1) {
		LST_addElement(neighborsFirstAtom, neighbor(v,0));
		SHL_removeEdge(m, idv, neighbor(v,0));
	}

	//Position du deuxième sommet : centre du motif
	//Hydrogène+taille d'une liaison simple moyenne.
	dir = normalization(dir, (DIST_SIMPLE / 2) + (MINDIS / 2));
	new_coords = addPoint(coords(v), dir);
	for (int i = 0; i < size(m); i++) {
		if (flag(atom(m, i)) == HYDRO_BOND_F && dist(coords(atom(m, i)),new_coords) < MINDIS) {
			flag(v) = SHELL_F;
			return;
		}
	}
	idc = SHL_addAtom(m, new_coords, -1);

	//checkInsertVertex(m, l, idc);
	SHL_addEdge(m, idv, idc);
	c = atom(m,idc);
	flag(c) = HYDRO_BOND_F;

	//Position du deuxième : (rotation de normal, 120, -dir) + coords(v)
	dir = subPoint(initPoint(0), normalization(dir, DIST_SIMPLE));
	new_coords = addPoint(coords(c), rotation(normal, 120, dir));
	index = SHL_addAtom(m, new_coords, -1);

	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idc, index);
	if (flag(atom(m,index)) < LINKABLE_F)
		flag(atom(m,index)) = LINKABLE_F;

	//Position du troisième : (rotation de normal, 120, -dir) + coords(v)
	new_coords = addPoint(coords(c), rotation(normal, -120, dir));
	index = SHL_addAtom(m, new_coords, -1);

	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idc, index);
	if (flag(atom(m,index)) < LINKABLE_F)
		flag(atom(m,index)) = LINKABLE_F;

	//Rattacher les nouvaux sommets à ceux de la liste l.
	//SHL_linkBorder(m, idc, l);
	
	LST_delete(neighborsFirstAtom);
}

/**
 * @brief Insert a hydrogen acceptor pattern with a tetrahedral geometry.  
 * 
 * @param m Envelope with the beginning of the cage.
 * @param idv Index of the heteroatom.
 * @param normal Normal vector starting point.
 * @param dir Point in the direction of the normal vector.
 */
void insertAcceptor2(Shell_t* m, unsigned idv, Point_t normal, Point_t dir) {
	//La position de v est la position du premier.
	int index, idc;
	Point_t x1, x2, x3, x4;
	AtomShl_t* v = atom(m, idv), *c;
	List_t* l = LST_create();

	flag(v) = HYDRO_BOND_F;

	//Insérer tous les voisins de v dans la liste.
	while (neighbor(v,0) != -1) {
		LST_addElement(l,neighbor(v,0));
		SHL_removeEdge(m, idv, neighbor(v,0));
	}

	x1 = coords(v);

	//Position du deuxième sommet : centre du motif
	//Hydrogène+taille d'une liaison simple moyenne.
	dir = normalization(dir, (DIST_SIMPLE/2)+(MINDIS/2));
	x2 = addPoint(coords(v), dir);
	for (int l = 0; l< size(m); l++){
		if(flag(atom(m,l)) == HYDRO_BOND_F && dist(coords(atom(m,l)),x2) < MINDIS){
			flag(v) = SHELL_F;
			return;
		}
	}
	idc = SHL_addAtom(m, x2, -1);

	//checkInsertVertex(m, l, idc);
	SHL_addEdge(m, idv, idc);
	c = atom(m,idc);
	flag(c) = HYDRO_BOND_F;

	//Deuxième sommet du tétraèdre
	//x2 = AX1E3(coords(c), x1, normal, DIST_SIMPLE);
	x2 = normalization(vector(x1, coords(c)), 1);
	x2 = rotation(normal, 109.47, x2);
	x2 = normalization(x2, DIST_SIMPLE);
	x2 = addPoint(coords(c), x2);
	index = SHL_addAtom(m, x2, -1);

	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idc, index);
	flag(atom(m,index)) = LINKABLE_F;

	//Troisième sommet du tétraèdre
	x3 = AX2E2(coords(c), x1, x2, DIST_SIMPLE);
	index = SHL_addAtom(m, x3, -1);

	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idc, index);
	flag(atom(m,index)) = LINKABLE_F;

	//Troisième sommet du tétraèdre
	x4 = AX3E1(coords(c), x1, x2, x3, DIST_SIMPLE);
	index = SHL_addAtom(m, x4, -1);

	//checkInsertVertex(m, l, index);
	SHL_addEdge(m, idc, index);
	flag(atom(m,index)) = LINKABLE_F;

	//Rattacher les nouvaux sommets à ceux de la liste l.
	//SHL_linkBorder(m, idc, l);
	
	LST_delete(l);
}

/**
 * @brief Add hydrogen patterns (donor or acceptor) to the beginning of the cage.
 * 
 * @param m Gathering of the main elements (substrate and envelope).
 */
void generationHydro(Main_t* m) {

	AtomShl_t *atomShell;
	Atom_t *parentAtomSub;

	//TODO choose to keep the loop to make the dependency graph
	for (int i = 0; i < 1/*mocSize(m)*/; i++) {
		for (int j = 0; j < size(bond(moc(m, i))); j++) {
			// Get the atom id from the dependency graph
			int idAtomShell = id(vertex(bond(moc(m, i)), j));

			if (idAtomShell != -1) {
				atomShell = atom(moc(m,i), idAtomShell);
				int tooClose = 1;
				for (int l = 0; l< size(moc(m,i)); l++) {
					if(flag(atom(moc(m, i), l)) == HYDRO_BOND_F && dist(coords(atom(moc(m, i), l)), coords(atomShell)) < MINDIS) {
						tooClose = 0;
						break;
					}
				}
				if(tooClose) {
					parentAtomSub = atom(substrat(m), parentAtom(atomShell));
					if (!strcmp(symbol(parentAtomSub), "H")) {
						insertDonor1(moc(m,i), idAtomShell, MOL_seekNormal(substrat(m), parentAtom(atomShell), -1), 
							vector(coords(parentAtomSub), coords(atomShell)));
					}
					else {
						int haveTriangularGeometry = (steric(parentAtomSub) == 3);
						if (haveTriangularGeometry) {
							insertAcceptor1(moc(m,i), idAtomShell, MOL_seekNormal(substrat(m), parentAtom(atomShell), -1), 
								vector(coords(parentAtomSub), coords(atomShell)));
						}
						else {
							insertAcceptor2(moc(m,i), idAtomShell, MOL_seekNormal(substrat(m), parentAtom(atomShell), -1), 
								vector(coords(parentAtomSub), coords(atomShell)));
						}
					}
				}
			}
		}
		SHL_testDis(moc(m,i));
	}
}


/**
 * @brief Add aromatic rings to the shell.
 * 
 * The function find for each atom involved in a cycle
 * if it can be part of a triangular pattern,
 * then add the pattern to the envelpe.
 * 
 * @param s Envelope of the substrate.
 */
void generationCycle(Shell_t* s) {
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

void generationMoc(Main_t* m) {

	envarom(m) = SHL_copy(envelope(m));
	//generationDep(m);
	printf("### Génération des motifs aromatiques ###\n");
	generationCycle(envarom(m));
	//SHL_write(envarom(m));

	printf("### Génération des motifs hydrogènes ###\n");
	MN_copyMoc(m, envarom(m));
	generationDep(m);
	generationHydro(m);

}
