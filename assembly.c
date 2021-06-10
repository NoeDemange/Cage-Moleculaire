#include "assembly.h"

// Insertion du motif passé en argument
void insererMotif(){
	
}

// Génère le chemin entre 2 groupements de motifs
void genererChemin(){
	
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

// Vérifier s'il existe plusieurs groupements de motifs non relié
int checkGroupement(Shell_t* s){
	
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
						return 1; // Il existe plusieurs groupements
					}
				}
			}
		}
	}
	
	return 0;
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

// Fonction principale
void assemblage(Main_t* m, Shell_t* envelope2){
	
	for (int i=0; i</*mocSize(m)*/1	; i++) //Pour tous les mocs
	{
		checkGroupement(moc(m, i));
		/*while (checkGroupement(moc(m, i))) // Tant qu'il existe au moins 2 groupements de motifs
		{
			printf("1111111");
			List_p* sommets = choixSommets(moc(m, i));
			printf("2222222");
			printf("\n %p \n",sommets);
			while (sommets->premier) // Pour tous les couples de sommets à relier
			{
				printf("33333333");
				int depart = sommets->premier->depart;
				int arrivee = sommets->premier->arrivee;
				printf("444444444");
				LST2_removeFirst(sommets);
				printf("555555555");
				// Choix sommets intermédiaire et generer directions
				// genererChemin();
				
			}
			
		}*/
	}
}

void assemblage2(Main_t* m, Shell_t* envelope2){
	
	for (int i=0; i</*mocSize(m)*/1	; i++) //Pour tous les mocs
	{
		List_p* sommets = choixSommets(moc(m, i));
		
		while (sommets->premier) // Tant qu'il existe au moins 2 groupements de motifs
		{
			while (sommets->premier) // Pour tous les couples de sommets à relier
			{
				printf("33333333");
				int depart = sommets->premier->depart;
				int arrivee = sommets->premier->arrivee;
				printf("444444444");
				LST2_removeFirst(sommets);
				printf("555555555");
				
				// Choix sommets intermédiaire et generer directions
				// genererChemin();
			}
			sommets = choixSommets(moc(m, i));
		}
	}
}
