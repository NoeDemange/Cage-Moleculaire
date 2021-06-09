#include "assembly.h"

// Insertion du motif passé en argument
void insererMotif(){
	
}

// Génère le chemin entre 2 groupements de motifs
void genererChemin(){
	
}

// Parcours en profondeur en fonction des indices
int parcours(Shell_t* s, List_t* marquer, int indice1, int indice2) {
	printf("DEBUT");
	AtomShl_t* a = atom(s, indice1);
	if (neighborhoodSize(a) == 0)
	{
		return 0;
	}
	else
	{
		for (int i = 0; i < neighborhoodSize(a); i++) // Pour tous les voisins de a
		{
			if (neighbor(a, i) == indice2) // S'il l'identifiant recherché est trouvé
			{
				return 1;
			}
			else
			{
				if (LST_getIndice(marquer, neighbor(a, i) != -1 )) // Si l'identifiant de ce sommet est déjà marqué
				{
					return 0;
				}
				else
				{
					printf("\n %d -> %d \n", indice1, indice2);
					LST_addElement(marquer, neighbor(a, i));
					int valide = parcours(s, marquer, neighbor(a, i), indice2);
					if (valide)
					{
						return 1;
					}
				}
			}
		}
	}
	return 0;
}

// Vérifie s'il existe un chemin passant seulement par des sommets de priorité 1 ou 3
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
	
	for (int i = 0; i < SHL_nbAtom(s) - 1; i++) //Pour tous les sommets de priorité 1
	{
		printf("\n %d -> flag : %d\n", i, flag(atom(s, i)));
		if (atom(s, i)->flag == 1)
		{
			for (int j = i+1; j < SHL_nbAtom(s); j++)
			{
				if (atom(s, j)->flag == 1)
				{
					if (!existeChemin(s, i, j))
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
	
	for (int i = 0; i < SHL_nbAtom(s) - 1; i++) //Pour tous les sommets de priorité 1
	{
		if (atom(s, i)->flag == 1)
		{
			for (int j = i+1; j < SHL_nbAtom(s); j++)
			{
				if (atom(s, j)->flag == 1)
				{
					if (!existeChemin(s, i, j)) // Si i et j sont de groupements différents
					{
						LST2_addElement(sommets, i, j);
					}
				}
			}
		}
	}
	
	Element* e = sommets->premier;
	while (e)
	{
		printf("\n%d %d\n", e->depart+1, e->arrivee+1);
		e = e->suivant;
	}
	
}

// Fonction principale
void assemblage(Main_t* m, Shell_t* envelope2){
	
	for (int i=0; i</*mocSize(m)*/1	; i++) //Pour tous les mocs
	{
		while (checkGroupement(moc(m, i))) // Tant qu'il existe au moins 2 groupements de motifs
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
			
		}
	}
}
