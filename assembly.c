#include "assembly.h"
#include "interface.h"
#include <float.h>

#define NB_MOTIF 4
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
}

// Determine la position de l'atome a inserer
// et l'ajoute a l'enveloppe si c'est possible
int projection(Shell_t* mocTraite, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Point_t arrivee) {
	
	Point_t positionNvDprt = PT_init();
	//Emplacement a refaire
	positionNvDprt = vector(arrivee, coords(atom(mocTraite, depart)));
	positionNvDprt = normalization(positionNvDprt, SIMPLE);
	
	// Si le point n'est pas dans l'enveloppe
	ajoutProjection(mocTraite, mocAtt, depart, nvDepart, numMotif, positionNvDprt); // Ajout a l'enveloppe
	
	return 1; //-1; // Si l'ajout n'a pas pu se faire
}

// Insertion du motif passé en argument
void insererMotif(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Point_t arrivee){
	
	if ( flag(atom(moc, depart)) == 1 || ( flag(atom(moc, depart)) == 4 && LST_nbElements(neighborhood(atom(moc, depart))) == 1 ) ) // Oxygene ou Carbone avec 1 voisin
	{
		//Projection
		//Diff rotations
		projection(moc, mocAtt, depart, nvDepart, numMotif, arrivee);
		
	}
	else if (flag(atom(moc, depart)) == 3) // Azote
	{
		if (LST_nbElements(neighborhood(atom(moc, depart))) == 1) // 1 voisin
		{
			// 2 Projections
			projection(moc, mocAtt, depart, nvDepart, numMotif, arrivee);
		}
		else // 2 voisins
		{
			//Projection
			projection(moc, mocAtt, depart, nvDepart, numMotif, arrivee);
		}
	}
	else if (flag(atom(moc, depart)) == 4) // Carbone
	{
		if (LST_nbElements(neighborhood(atom(moc, depart))) == 2) // 2 voisins
		{
			if (flag(atom(moc, neighbor(atom(moc, depart), 0))) == 1 || flag(atom(moc, neighbor(atom(moc, depart), 1))) == 1) // Si 1 des 2 voisins est un oxygene
			{
				//Projection
				projection(moc, mocAtt, depart, nvDepart, numMotif, arrivee);
			}
			else
			{
				// 2 Projections
				projection(moc, mocAtt, depart, nvDepart, numMotif, arrivee);
			}
		}
		else // 3 voisins
		{
			//Projection
			projection(moc, mocAtt, depart, nvDepart, numMotif, arrivee);
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
void genererChemin(List_m* mocAtt, Shell_t* mocTraite, int depart, int arrivee, Elem_s* sommetInter){
	
	for (int i = 0; i < NB_MOTIF; i++)
	{
		List_m* moc = LSTm_init();
		List_d* nvDepart = LSTd_init();
		
		insererMotif(mocTraite, moc, depart, nvDepart, i, sommetInter->sommet);
		
		while (moc->premier)
		{
			if (eloigne( coords(atom(mocTraite, depart)), coords(atom(moc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, arrivee)) )) // Si le nv depart est plus éloigné 
			{
				if (sommetInter->suivant != NULL) // Si ce n'est pas la derniere arrivee
				{
					genererChemin(mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, sommetInter->suivant);
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
						//printf("Modifier angles\n");
						// Modifier angles
						LSTm_addElement(mocAtt, SHL_copy(moc->premier->moc)); // Ajout dans la liste a traiter
					}
				}
			}
			else
			{
				genererChemin(mocAtt, moc->premier->moc, nvDepart->premier->sommet, arrivee, sommetInter);
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
 
void majDistances(Shell_t* s, float* dist, int* predecesseur, int s1, int s2) {
	float poids = PT_distance(coords(atom(s, s1)), coords(atom(s, s2)));
	if (dist[s2] > dist[s1] + poids)
	{
		dist[s2] = dist[s1] + poids;
		predecesseur[s2] = s1;
	}
}

// Plus court chemin
int* dijkstra(Shell_t* s, int depart, int arrivee) {
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
			majDistances(s, dist, predecesseur, s1, s2);
		}
	}
	
	free(dist);
	LSTd_delete(Q);
	
	return predecesseur;
}

// Determine les sommets intermédiaires du chemin 
List_s* sommetIntermediaire(Shell_t* s, int depart, int arrivee) {
	
	int* predecesseur = dijkstra(s, depart, arrivee);
	//printf("Pred : %p\n", predecesseur);
	
	List_s* sommets = LSTs_init();
	
	int si = arrivee;
	while (si != depart)
	{
		// flag(atom(s, si)) = 2; // Visualisation 
		LSTs_addElement(sommets, coords(atom(s, si)));
		si = predecesseur[si];
	}
	
	free(predecesseur);
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
	
	while (mocAtt->premier) // Tant qu'il existe un moc a traiter
	{
		List_p* sommets = choixSommets(mocAtt->premier->moc);
		printf("111111");
		if (!sommets->premier) // S'il n'y qu'un groupement de motifs
		{
			printf("222222");
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
					
				List_s* sommetInter = sommetIntermediaire(mocTraite2, depart, arrivee); // Choix sommets intermédiaires
				
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
				genererChemin(mocAtt, mocTraite2, depart, arrivee, sommetInter->premier);
				printf("777777777");
				
				LSTm_addElement(mocAtt, mocTraite2); // Ajout dans la liste a traiter
				LSTs_delete(sommetInter); // Supprime la liste des sommets intermediaires
			}
			//SHL_delete(mocTraite);
		}
	}
	free(mocAtt);
}
