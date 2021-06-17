#include "assembly.h"
#include <float.h>

// Insertion du motif passé en argument
void insererMotif(){
	
}

// Génère le chemin entre 2 groupements de motifs
void genererChemin(){
	
}

void initDijkstra(Shell_t* s, int depart, float** dist, int** predecesseur, List_d** Q) {
	
	*dist = malloc(sizeof(float) * size(s));
	for (int i = 0; i < size(s); i++)
	{
		*dist[i] = FLT_MAX;
	}
	*dist[depart] = 0;
	
	predecesseur = malloc(sizeof(int) * size(s));
	for (int i = 0; i < size(s); i++)
	{
		*predecesseur[i] = -1;
	}
	
	*Q = LSTd_init();
	for (int i = 0; i < size(s) ; i++)
	{
		if (flag(atom(s, i)) == 0)
		{
			LSTd_addElement(*Q, i);
		}
	}
}

int trouveMin(float* dist, List_d* Q) {
	float mini = FLT_MAX;
	int sommet = -1;
	Elem_d* cursor = Q->premier;
	while (cursor)
	{
		if (dist[cursor->sommet] < mini)
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
int* dijkstra(Shell_t* s, int depart) {
	float* dist;
	int* predecesseur;
	List_d* Q;
	
	initDijkstra(s, depart, &dist, &predecesseur, &Q);
	
	while (Q->premier) {
		int s1 = trouveMin(dist, Q);
		LSTd_removeSommet(Q, s1);
		
		for (int i = 0; i < neighborhoodSize(atom(s,s1)); i++)
		{
			int s2 = neighbor(atom(s,s1), i);
			majDistances(s, dist, predecesseur, s1, s2);
		}
	}
	
	free(dist);
	LSTd_delete(Q);
	
	return predecesseur;
}

// Determine les sommets intermédiaires du chemin 
List_s* sommetIntermediaire(Shell_t* s, int depart, int arrivee) {
	
	int* predecesseur = dijkstra(s, depart);
	
	List_s* sommets = LSTs_init();
	
	int si = arrivee;
	while (si != depart)
	{
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
		if (flag(atom(s, neighbor(sommet, i))) == 0 && flag(sommet) != 0 ) /*&& flag(sommet) != 4*/
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
	
	/*Element* e = sommets->premier;
	while (e)
	{
		printf("\n%d %d\n", e->depart+1, e->arrivee+1);
		e = e->suivant;
	}
	*/
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
	
	for (int i=0; i<mocSize(m); i++) //Pour tous les mocs
	{
		if (moc(m,i) != NULL)
		{
			LSTm_addElement(mocAtt, moc(m, i)); // Les mettre dans la liste a traiter
			//moc(m, i) = NULL; // Les supprimer du tableau de solutions finales
			
			/*for (int j = 0; j <mocSize(m) ; j++)
			{
				printf("\nMocs : %p", moc(m, j));
			}
			printf("\nTaille moc : %d \n", mocSize(m));*/
		}
	}
	//free(m->mocs);
	//m->mocs = NULL;
	//mocSize(m) = 0;
	
	return mocAtt;
}

// Fonction principale
void assemblage(Main_t* m, Shell_t* envelope2){
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
