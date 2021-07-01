#include "assembly.h"
#include "interface.h"
#include "expansion.h"
#include <float.h>

#define NB_MOTIF 3
#define MIN_DIST 3

// Donne le type de l'atome inserer
int typeInsert(int numMotif){
	if (numMotif == 0) // Oxygene
	{
		return 1;
	}
	else if (numMotif == 1) // Azote
	{
		return 3;
	}
	else //Carbone
	{
		return 4;
	}
}

// Ajout l'atome projeté a l'enveloppe
void ajoutProjection(Shell_t* mocTraite, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Point_t positionNvDprt) {
	Shell_t* moc = SHL_copy(mocTraite);
	
	int id = SHL_addAtom(moc, positionNvDprt, -1);
	flag(atom(moc, id)) = typeInsert(numMotif);
	SHL_addEdge(moc, depart, id);
	
	LSTm_addElement(mocAtt, moc);
	LSTd_addElement(nvDepart, id);
	
	// Visualisation
	Point_t v = positionNvDprt;
	v.z += 0.5;
	int id2 = SHL_addAtom(moc, v, -1);
	flag(atom(moc, id2)) = 2;
}

// Determine la position de l'atome a inserer
// et l'ajoute a l'enveloppe si c'est possible
void projection(Molecule_t* mol, Shell_t* mocTraite, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, int arrivee) {
	
	Point_t dpt = coords(atom(mocTraite, depart));
	Point_t arv = coords(atom(mocTraite, arrivee));
	
	//Emplacement a refaire
	Point_t positionNvDprt = addPoint(dpt, normalization(vector(dpt, arv), SIMPLE));
	printf("%d %d\n", parentAtom(atom(mocTraite, arrivee)), size(mol));
	if (parentAtom(atom(mocTraite, arrivee)) != -1) // Derniere arrivee
	{
		Point_t parent = coords( atom(mol, parentAtom(atom(mocTraite, arrivee))) );
		positionNvDprt = addPoint(positionNvDprt, normalization(vector(parent, arv), 0.5)); // Décalage vers exterieur
		positionNvDprt = addPoint(dpt, normalization(vector(dpt, positionNvDprt), SIMPLE)); // Remet à la taille SIMPLE
	}
	
	// Si le point n'est pas dans l'enveloppe
	ajoutProjection(mocTraite, mocAtt, depart, nvDepart, numMotif, positionNvDprt); // Ajout a l'enveloppe
}

// Projection pour un atome avec 1 voisin
void projectionOCN_AX1E3(Shell_t* moc, List_m* mocAtt, int depart, int arrivee, List_d* nvDepart, int numMotif) {
	
	List_s* positions = LSTs_init();
	Point_t dpt = coords(atom(moc, depart));
	Point_t arv = coords(atom(moc, arrivee));
	Point_t v1 = coords(atom(moc, neighbor(atom(moc, depart), 0)));
	
	Point_t normal = normalization(PT_alea(), 1);//normalPerpendiculaire(dpt, v1, PT_alea(), 1);
	
	Point_t positionNvDprt = AX1E3(dpt, v1, normal, SIMPLE);
	//ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt); // Ajout a l'enveloppe
	
	for (int i = 0; i < 12; i++) // Rotation a 360
	{
		// Rotation de 30° de la normal
		normal = rotation(normalization(vector(dpt, v1), 1),  30, normal);
		//normal = rotation30(dpt, normal, addPoint(dpt, vector(dpt, v1)), 1);
		positionNvDprt = AX1E3(dpt, v1, normal, SIMPLE);
		
		// Si le point n'est pas dans l'enveloppe
		LSTs_addElement(positions, positionNvDprt);
	}
	
	for (int i = 0; i < 3; i++) // 3 positions les mieux placés 
	{
		positionNvDprt = distMin(positions, arv); 
		LSTs_removeElement(positions, positionNvDprt);
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt); // Ajout a l'enveloppe
	}
	
	/*// Rotation de 30° du point
	normal = normalization(vector(atom(moc, depart), neighbor(atom(moc, depart), 0)), 1);
	positionNvDprt = rotation( addPoint(atom(moc, depart), normal , 30, positionNvDprt );
	positionNvDprt = rotation30(coords(atom(moc, depart)), positionNvDprt, normal, SIMPLE);
	*/
}

// Projection pour un azote avec 2 voisins
void projectionN_AX2E2(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif) {
	
	Point_t positionNvDprt = AX2E2(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), SIMPLE);
	// Si le point n'est pas dans l'enveloppe
	ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt); // Ajout a l'enveloppe
}

// Projection pour un carbone avec 2 voisins dont un oxygene
void projectionC_AX2E1(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif) {
	
	Point_t positionNvDprt = AX2E1(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), SIMPLE);
	// Si le point n'est pas dans l'enveloppe
	ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt); // Ajout a l'enveloppe
}

// Projection pour un carbone avec 2 voisins
void projectionC_AX2E2(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif) {
	
	Point_t positionNvDprt = AX2E2(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), SIMPLE);
	// Si le point n'est pas dans l'enveloppe
	ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt); // Ajout a l'enveloppe
	
	Point_t positionNvDprt2 = AX3E1(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), positionNvDprt, SIMPLE);
	// Si le point n'est pas dans l'enveloppe
	ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt2); // Ajout a l'enveloppe
}

// Projection pour un carbone avec 3 voisins
void projectionC_AX3E1(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif) {
	
	Point_t positionNvDprt = AX3E1(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), coords(atom(moc, neighbor(atom(moc, depart), 2))), SIMPLE);
	// Si le point n'est pas dans l'enveloppe
	ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt); // Ajout a l'enveloppe
}

// Insertion du motif passé en argument
void insererMotif(Molecule_t* mol, Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, int arrivee){
	
	if ( LST_nbElements(neighborhood(atom(moc, depart))) == 1 ) // Oxygene ou Azote ou Carbone avec 1 voisin
	{
		//Projection
		//Diff rotations
		//projection(mol, moc, mocAtt, depart, nvDepart, numMotif, arrivee);
		projectionOCN_AX1E3(moc, mocAtt, depart, arrivee, nvDepart, numMotif);
	}
	else if (flag(atom(moc, depart)) == 3 && LST_nbElements(neighborhood(atom(moc, depart))) == 2) // Azote avec 2 voisins
	{
		//Projection
		//projection(mol, moc, mocAtt, depart, nvDepart, numMotif, arrivee);
		projectionN_AX2E2(moc, mocAtt, depart, nvDepart, numMotif);
	}
	else if (flag(atom(moc, depart)) == 4) // Carbone
	{
		if (LST_nbElements(neighborhood(atom(moc, depart))) == 2) // 2 voisins
		{
			if (flag(atom(moc, neighbor(atom(moc, depart), 0))) == 1 || flag(atom(moc, neighbor(atom(moc, depart), 1))) == 1) // Si 1 des 2 voisins est un oxygene
			{
				//Projection
				projectionC_AX2E1(moc, mocAtt, depart, nvDepart, numMotif);
			}
			else
			{
				// 2 Projections
				projectionC_AX2E2(moc, mocAtt, depart, nvDepart, numMotif);
			}
		}
		else // 3 voisins
		{
			//Projection
			projectionC_AX3E1(moc, mocAtt, depart, nvDepart, numMotif);
		}
		
	}
	
	// Si motif 4 choix position O en plus
}

// Calcule si le nouveau depart est plus loin de l'arrivée que l'ancien départ
int eloigne(Point_t depart, Point_t nvDepart, Point_t arrivee){
	
	float d1 = dist(depart, arrivee);
	float d2 = dist(nvDepart, arrivee);
	
	if (d1 > d2)
	{
		return 0; // L'ancien est plus eloigné
	}
	else
	{
		return 1; // Le nouveau est plus eloigné
	}
	
}

// Génère le chemin entre 2 groupements de motifs
void genererChemin(Molecule_t* mol, List_m* mocAtt, Shell_t* mocTraite, int depart, int arrivee, Elem_d* sommetInter){
	
	for (int i = 0; i < 1/*NB_MOTIF*/; i++)
	{
		List_m* moc = LSTm_init();
		List_d* nvDepart = LSTd_init();
		
		insererMotif(mol, mocTraite, moc, depart, nvDepart, i, sommetInter->sommet);
		
		while (moc->premier)
		{
			if (eloigne( coords(atom(mocTraite, depart)), coords(atom(moc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, arrivee)) )) // Si le nv depart est plus éloigné 
			{
				if (sommetInter->suivant != NULL) // Si ce n'est pas la derniere arrivee
				{
					genererChemin(mol, mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, sommetInter->suivant);
				}
				else // C'est la derniere arrivee
				{
					if (dist( coords(atom(moc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, arrivee)) ) < MIN_DIST) // Proche de l'arrivée 
					{
						SHL_addEdge(moc->premier->moc, nvDepart->premier->sommet, arrivee);// Ajout lien entre dernier sommet du chemin et arrivee
						LSTm_addElement(mocAtt, SHL_copy(moc->premier->moc)); // Ajout dans la liste a traiter
					}
					else
					{
						printf("Modifier angles\n");
						// Modifier angles
						SHL_addEdge(moc->premier->moc, nvDepart->premier->sommet, arrivee);
						LSTm_addElement(mocAtt, SHL_copy(moc->premier->moc)); // Ajout dans la liste a traiter
					}
				}
			}
			else
			{
				genererChemin(mol, mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, sommetInter);
			}
			LSTm_removeFirst(moc);
			LSTd_removeFirst(nvDepart);
		}
	}
}

// Génère le chemin entre 2 groupements de motifs
void genererChemin2(Molecule_t* mol, List_m* mocAtt, Shell_t* mocTraite, int depart, int arrivee, Elem_d* sommetInter){
	
	for (int i = 0; i < 1/*NB_MOTIF*/; i++)
	{
		List_m* moc = LSTm_init();
		List_d* nvDepart = LSTd_init();
		
		insererMotif(mol, mocTraite, moc, depart, nvDepart, i, sommetInter->sommet);
		
		if (moc->premier)
		{
			if (eloigne( coords(atom(mocTraite, depart)), coords(atom(moc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, arrivee)) )) // Si le nv depart est plus éloigné 
			{
				if (sommetInter->suivant != NULL) // Si ce n'est pas la derniere arrivee
				{
					genererChemin2(mol, mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, sommetInter->suivant);
				}
				else // C'est la derniere arrivee
				{
					if (dist( coords(atom(moc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, arrivee)) ) < MIN_DIST) // Proche de l'arrivée 
					{
						SHL_addEdge(moc->premier->moc, nvDepart->premier->sommet, arrivee);// Ajout lien entre dernier sommet du chemin et arrivee
						LSTm_addElement(mocAtt, SHL_copy(moc->premier->moc)); // Ajout dans la liste a traiter
						return;
					}
					else
					{
						//printf("Modifier angles\n");
						// Modifier angles
						SHL_addEdge(moc->premier->moc, nvDepart->premier->sommet, arrivee);
						LSTm_addElement(mocAtt, SHL_copy(moc->premier->moc)); // Ajout dans la liste a traiter
						return;
					}
				}
			}
			else
			{
				genererChemin2(mol, mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, sommetInter);
			}
			LSTm_removeFirst(moc);
			LSTd_removeFirst(nvDepart);
		}
	}
}

void initDijkstra(Shell_t* s, int depart, int arrivee, float** dist, int** predecesseur, List_d** Q) {
	
	*dist = malloc(sizeof(float) * SHL_nbAtom(s));
	//printf("%d : %d : %d\n", size(s), depart, SHL_nbAtom(s));
	
	for (int i = 0; i < SHL_nbAtom(s); i++)
	{
		//printf("Atome : %d\n", flag(atom(s, i)));
		(*dist)[i] = FLT_MAX;
	}
	(*dist)[depart] = 0;
	
	//printf("Dist\n");
	*predecesseur = malloc(sizeof(int) * SHL_nbAtom(s));
	for (int i = 0; i < SHL_nbAtom(s); i++)
	{
		(*predecesseur)[i] = -1;
	}
	//printf("Pred\n");
	*Q = LSTd_init();
	for (int i = 0; i < SHL_nbAtom(s) ; i++)
	{
		if (flag(atom(s, i)) == 0 || depart == i || arrivee == i)
		{
			LSTd_addElement(*Q, i);
		}
	}
	//printf("Q\n");
}

int trouveMin(float* dist, List_d* Q) {
	float mini = FLT_MAX;
	int sommet = -5;
	Elem_d* cursor = Q->premier;
	while (cursor)
	{
		//printf("%f %d \n", dist[cursor->sommet], cursor->sommet);
		
		if (dist[cursor->sommet] <= mini)
		{
			mini = dist[cursor->sommet];
			sommet = cursor->sommet;
		}
		cursor = cursor->suivant;
	}
	
	return sommet;
}

// Verifie si le point passer en argument est trop interne
// dans ce cas a enlever de la liste
int ptInterne(Shell_t* envelope2, Point_t p) {
	//printf("Nombre : %d\n", SHL_nbAtom(envelope2));
	for (int i = 0; i < SHL_nbAtom(envelope2) ; i++)
	{
		if (p.x == atomX(atom(envelope2, i)) && p.y == atomY(atom(envelope2, i)) && p.z == atomZ(atom(envelope2, i)))
		{
			return 0; // Il est sur l'enveloppe grossiere donc pas interne
		}
	}
	return 1; // Il n'est pas sur l'enveloppe grossiere donc trop interne
} 

void majDistances(Shell_t* envelope2, Shell_t* s, float* dist, int* predecesseur, int s1, int s2, int arrivee) {
	float poids = PT_distance(coords(atom(s, s1)), coords(atom(s, s2)));
	if ( ptInterne(envelope2, coords(atom(s,s2))) && s1!=arrivee && s2!=arrivee ) // Si le point est trop interne le passer
	{
		poids += 1000;
	}
	
	if (dist[s2] > dist[s1] + poids)
	{
		dist[s2] = dist[s1] + poids;
		predecesseur[s2] = s1;
	}
}

// Plus court chemin
int* dijkstra(Shell_t* envelope2, Shell_t* s, int depart, int arrivee) {
	float* dist = NULL;
	int* predecesseur = NULL;
	List_d* Q = NULL;
	//printf("\nDijkstra %d\n", depart);
	initDijkstra(s, depart, arrivee, &dist, &predecesseur, &Q);
	//printf("%p %p %p \n", dist, predecesseur, Q);
	
	while (Q->premier) {
		//printf("TMin\n");
		int s1 = trouveMin(dist, Q);
		//printf("\nS1 : %d\n",s1);
		//printf("RmvSommet\n");
		LSTd_removeSommet(Q, s1);
		//printf("for\n");
		//printf("%p\n", neighborhood(atom(s,s1)));
		/*for (int i = 0; i < neighborhoodSize(atom(s,s1)); i++)
		{
			printf("%d\n", neighbor(atom(s,s1), i));
		}*/
		
		for (int i = 0; i < LST_nbElements(neighborhood(atom(s,s1))); i++)
		{
			//printf("Voisin %d\n", i); 
			int s2 = neighbor(atom(s,s1), i);
			//printf("MAJD\n");
			majDistances(envelope2, s, dist, predecesseur, s1, s2, arrivee);
		}
	}
	
	free(dist);
	LSTd_delete(Q);
	
	return predecesseur;
}

// Determine les sommets intermédiaires du chemin 
List_d* sommetIntermediaire(Main_t* m, Shell_t* s, int depart, int arrivee) {
	
	double alpha2 = 20.0;
	Shell_t* envelope2 = createShell(substrat(m), alpha2); // Enveloppe grossiere
	//printf("ENVELOPPE2\n");
	int* predecesseur = dijkstra(envelope2, s, depart, arrivee);
	//printf("Pred : %p\n", predecesseur);
	
	List_d* sommets = LSTd_init();
	
	int si = arrivee;
	while (si != depart)
	{
		flag(atom(s, si)) = 2; // Visualisation 
		LSTd_addElement(sommets, si);
		si = predecesseur[si];
	}
	
	free(predecesseur);
	SHL_delete(envelope2);
	
	return sommets;
}

// Vérifie si le sommet se situe en bordure de motif
int bordureCheck(Shell_t* s, AtomShl_t* sommet) {
	
	for (int i = 0; i < neighborhoodSize(sommet); i++) // Pour tous les voisins du sommet
	{
		// Si ce sommet est dans un motif et qu'un de ses voisins est de l'enveloppe
		if (flag(atom(s, neighbor(sommet, i))) == 0 && flag(sommet) != 0 )
		{
			return 1; 
		}
	}
	
	return 0;
}

// Parcours en profondeur en fonction des indices sur les sommets des motifs uniquement
int parcours(Shell_t* s, List_t* marquer, int indice1, int indice2) {
	
	/*printf("marquer : ");
	for (int i = 0; i < size(marquer); i++)
	{
		printf("%d ", elts(marquer,i));

	}
	printf("\n");
	*/
	AtomShl_t* a = atom(s, indice1);
	LST_addElement(marquer, indice1);
	
	/*printf("Voisins de %d :\n", indice1);
	for (int i = 0; i < neighborhoodSize(a); i++) {
        printf("%d, ", neighbor(a, i));
    }
    printf("\n");
    */
	if (neighborhoodSize(a) == 0)
	{
		//printf("Pas de voisin\n");
		return 0;
	}
	else
	{
		for (int i = 0; i < neighborhoodSize(a) && neighbor(a, i) != -1; i++) // Pour tous les voisins de a
		{
			if (flag(atom(s, neighbor(a, i))) != 0) // Si le sommet est dans un motif donc de priorité != 0
			{
				if (neighbor(a, i) == indice2) // Si l'identifiant recherché est trouvé
				{
					//printf("%d atteint\n",neighbor(a, i));
					return 1;
				}
				else
				{
					if (!LST_check(marquer, neighbor(a, i))) // Si l'identifiant de ce sommet n'est pas déjà marqué
					{
						//printf("\n %d -> %d \n", indice1, neighbor(a, i));
						int valide = parcours(s, marquer, neighbor(a, i), indice2);
						if (valide)
						{
							return 1;
						}
					}
				}
			}
			
			
		}
	}
	return 0;
}

// Vérifie s'il existe un chemin passant seulement par des sommets qui appartiennent aux motifs donc de priorité != 0
// donc si les 2 sommets sont du même groupement
int existeChemin(Shell_t* s, int indice1, int indice2){
	
	List_t* marquer = NULL;
	marquer = LST_create();
	
	int existe = parcours(s, marquer, indice1, indice2);
	
	LST_delete(marquer);
	
	return existe;
}

// Génère tous les couples de sommets à relier possible entre des groupements
List_p* choixSommets(Shell_t* s){
	
	List_p* sommets = LST2_init();
	
	for (int i = 0; i < SHL_nbAtom(s) - 1; i++) //Pour tous les sommets en bordure
	{
		if ( bordureCheck(s, atom(s, i)) )
		{
			for (int j = i+1; j < SHL_nbAtom(s); j++)
			{
				if ( bordureCheck(s, atom(s, j)) )
				{
					if (!existeChemin(s, i, j)) // Si i et j sont de groupements différents
					{
						LST2_addElement(sommets, i, j);
					}
				}
			}
		}
	}
	
	/*printf("%d\n", SHL_nbAtom(s));
	printf("%d\n", size(s));
	Element* e = sommets->premier;
	while (e)
	{
		printf("\n%d %d\n", e->depart+1, e->arrivee+1);
		e = e->suivant;
	}*/
	
	return sommets;
}

void affichage(Shell_t* s) {
	
	for (int i = 0; i < SHL_nbAtom(s); i++)
	{
		printf("Atome %d : ", i+1);
		for (int j = 0; j < neighborhoodSize(atom(s,i)); j++)
		{
			printf("%d ",neighbor(atom(s,i),j)+1);
		}
		printf("flag %d\n",flag(atom(s,i)));
	}
}

// Cree une liste de moc a traiter et vide le tableau de solutions finales
List_m* initMocAtt(Main_t* m){
	List_m* mocAtt = LSTm_init();
	
	for (int i=0; i</*mocSize(m)*/1; i++) //Pour tous les mocs
	{
		if (moc(m,i) != NULL)
		{
			LSTm_addElement(mocAtt, moc(m, i)); // Les mettre dans la liste a traiter
			moc(m, i) = NULL; // Les supprimer du tableau de solutions finales
			
			/*printf("\n");
			for (int j = 0; j <mocSize(m) ; j++)
			{
				printf("Mocs : %p \n", moc(m, j));
			}
			printf("\nTaille moc : %d \n", mocSize(m));
			
			Elem* e = mocAtt->premier;
			while (e)
			{
				printf("%p ", e->moc);
				e = e->suivant;
			}*/
			
		}
	}
	
	/*Elem* e = mocAtt->premier;
	while (e)
	{
		printf("%p ", e->moc);
		e = e->suivant;
	}*/
	free(m->mocs);
	m->mocs = NULL;
	mocSize(m) = 0;
	
	return mocAtt;
}

// Fonction principale
void assemblage(Main_t* m){
	List_m* mocAtt = initMocAtt(m);
	
	/*while (mocAtt->premier) // Tant qu'il existe un moc a traiter
	{
		List_p* sommets = choixSommets(mocAtt->premier->moc);
		
		if (!sommets->premier) // S'il n'y qu'un groupement de motifs
		{
			m->mocs[MN_getIndiceFree(m)] = mocAtt->premier->moc; // Ajout au tableau des solutions finales
			LSTm_removeFirst(mocAtt); // Suppression dans la liste a traiter
		}
		else // S'il y a au moins 2 groupements de motifs
		{
			Shell_t* mocTraite = mocAtt->premier->moc; // Copie le moc a traiter
			mocAtt->premier = mocAtt->premier->suivant; // Supprime de la liste à traiter
			
			while (sommets->premier) // Pour tous les couples de sommets à relier
			{
				Shell_t* mocTraite2 = SHL_copy(mocTraite); // Cree un nouveau moc dans la liste a traiter
				printf("33333333");
				int depart = sommets->premier->depart;
				int arrivee = sommets->premier->arrivee;
				printf("444444444");
				LST2_removeFirst(sommets);
				printf("555555555");
					
				// Choix sommets intermédiaire et generer directions
				// genererChemin();
				
				LSTm_addElement(mocAtt, mocTraite2); // Ajout dans la liste a traiter
			}
			SHL_delete(mocTraite);
		}
	}
	free(mocAtt);*/
}

void assemblage2(Main_t* m, int alpha){
	List_m* mocAtt = initMocAtt(m); // ! Prend le premier moc seulement
	int sol = 1;
	if (mocAtt->premier && sol) // Tant qu'il existe un moc a traiter
	{
		List_p* sommets = choixSommets(mocAtt->premier->moc);
		printf("111111");
		if (!sommets->premier) // S'il n'y qu'un groupement de motifs
		{
			printf("222222");
			sol = 0;
			m->mocs[MN_getIndiceFree(m)] = mocAtt->premier->moc; // Ajout au tableau des solutions finales
			mocAtt->premier = mocAtt->premier->suivant; // Suppression dans la liste a traiter
		}
		else // S'il y a au moins 2 groupements de motifs
		{
			Shell_t* mocTraite = mocAtt->premier->moc; // Copie le moc a traiter
			mocAtt->premier = mocAtt->premier->suivant; // Supprime de la liste à traiter
			
			if (sommets->premier) // Pour tous les couples de sommets à relier
			{
				Shell_t* mocTraite2 = SHL_copy(mocTraite); // Cree un nouveau moc dans la liste a traiter
				printf("33333333");
				int depart = sommets->premier->depart;
				int arrivee = sommets->premier->arrivee;
				printf("444444444");
				LST2_removeFirst(sommets);
				printf("555555555");
					
				List_d* sommetInter = sommetIntermediaire(m, mocTraite2, depart, arrivee); // Choix sommets intermédiaires
				
				/*Ashape_t* as3d = Cashape3d(envelope(m), alpha);
				double* point = malloc(2 * 3* sizeof(double));
				point[0] = 5.2408;     
				point[2] = -2.1055;
				point[4] = 2.4206;
				point[1] = 7;
				point[3] = -4;
				point[5] = 2;
				Point_t p = PT_init();
				p.x = point[0];
				p.y = point[2];
				p.z = point[4];
				Point_t p2 = PT_init();
				p2.x = point[1];
				p2.y = point[3];
				p2.z = point[5];
				int id = SHL_addAtom(moc(m,0), p, -1);
				flag(atom(moc(m,0),id)) = 2;
				id = SHL_addAtom(moc(m,0), p2, -1);
				flag(atom(moc(m,0),id)) = 2;
				Cinashape3d(as3d, point, 6);
				*/
				
				printf("666666");
				for (int i = 0; i < LST_nbElements(neighborhood(atom(mocTraite2, depart))); i++) // Retire les voisins enveloppe de l'atome de départ (bordure)
				{
					if (flag(atom(mocTraite2, neighbor(atom(mocTraite2, depart), i))) == 0)
					{
						SHL_removeEdge(mocTraite2, depart, i);
					}
				}
				
				for (int i = 1; i < 5 ; i++) // Atome de départ (sommet en bordure) est de n'importe quel type
				{
					if (i != 2)
					{
						flag(atom(mocTraite2, depart)) = i;
						genererChemin2(substrat(m), mocAtt, mocTraite2, depart, arrivee, sommetInter->premier);
					}
				}
				
				printf("777777777\n");
				
				LSTd_delete(sommetInter); // Supprime la liste des sommets intermediaires
			}
			//SHL_delete(mocTraite);
		}
	}
	int id = MN_getIndiceFree(m);
	//printf("%p %d\n", m->mocs, id);
	m->mocs[id] = mocAtt->premier->moc;
	
	free(mocAtt);
}
