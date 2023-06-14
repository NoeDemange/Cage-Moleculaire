#include "assembly.h"
#include "interface.h"
#include "expansion.h"
#include "output.h"
#include <float.h>

#define NB_MOTIF 5
#define MIN_DIST 2
#define DIST_ECART_CAGE 1.5
//#define DIST_ECART_ENVELOPPE 2
#define DIST_ECART_SUBSTRAT 2

void affichage(Shell_t* s) {
	
	for (int i = 0; i < SHL_nbAtom(s); i++)
	{
		printf("Atome %d : ", i/*+1*/);
		for (int j = 0; j < neighborhoodSize(atom(s,i)); j++)
		{
			printf("%d ",neighbor(atom(s,i),j)/*+1*/);
		}
		printf("flag %d\n",flag(atom(s,i)));
	}
}

// Vérifie si le point passé en argument est assez éloigné des autres atomes 
// Retourne 1 s'il n'est pas assez éloigné (<DIST_ECART_CAGE) et 0 sinon
int encombrement(Shell_t* moc, Molecule_t* sub, Point_t p){
	for(int i = 0; i < size(moc); i++){
		if(flag(atom(moc,i))!=0 && flag(atom(moc,i))!=-1){
			Point_t A = coords(atom(moc, i));
			/*if(flag(atom(moc,i))==0){
				if(dist(A, p) < DIST_ECART_ENVELOPPE)
					return 1;
			}
			else{*/
			if(dist(A, p) < DIST_ECART_CAGE)
				return 1;
		}
		//}
	}
	for(int i = 0; i < size(sub); i++){
		Point_t A = coords(atom(sub, i));
		if(dist(A, p) < DIST_ECART_SUBSTRAT) 
			return 1;
	}
			return 0;
}

/**************************************/
/* Ajout des motifs *******************/
/**************************************/

// Donne le type de l'atome inseré
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

// Ajout du motif 4 (cycle aromatique) perpendiculaire au plan avec son voisin
void ajoutMotif4(Shell_t* mocTraite, List_m* mocAtt, int depart, List_d* nvDepart, Point_t posNvDprt, Molecule_t* sub) {
	
	for (int i = 0; i < neighborhoodSize(atom(mocTraite,depart)); i++) // Pour tous les plans possibles avec les voisins
	{
		if (neighbor(atom(mocTraite, depart), i) != -1)
		{
			int idSuiv, idCycle;
			int inASh = 0;
			Point_t posSuiv;
			Point_t positionNvDprt = posNvDprt;
			Shell_t* moc = SHL_copy(mocTraite);
			
			Point_t v1 = coords(atom(moc, neighbor(atom(moc, depart), i)));
			Point_t dpt = coords(atom(moc, depart));
			
			// Ajouter le premier atome du cycle dont on a déjà la position
			int id = SHL_addAtom(moc, positionNvDprt, -1);
			flag(atom(moc, id)) = 4;
			SHL_addEdge(moc, depart, id);

			// Cherche la normal pour positionné le cycle
			Point_t normal = planNormal(positionNvDprt, dpt, v1);
			normal = rotation(normalization(vector(positionNvDprt, dpt), 1),  90, normal); // Perpendiculaire
			
			// Positionner les autres atomes du cycle
			v1 = AX1E2(positionNvDprt, coords(atom(moc, depart)), normal, SIMPLE); // Voisin
			positionNvDprt = AX2E1(positionNvDprt, coords(atom(moc, depart)), v1, SIMPLE); 
			if (encombrement(moc, sub, positionNvDprt) == 1) // Si le point n'est pas éloigné de 1.5 des autres atomes
			{
				inASh++;
			}

			idCycle = SHL_addAtom(moc, positionNvDprt, -1);
			flag(atom(moc, idCycle)) = 4;
			SHL_addEdge(moc, id, idCycle);
			
			int idVoisin = idCycle;
			for (int i = 0; i < 4; i++) 
			{
				v1 = coords(atom(moc, neighbor(atom(moc, idVoisin), 0)));
				positionNvDprt = AX1E2(positionNvDprt, v1, normal, SIMPLE);
				if (encombrement(moc, sub, positionNvDprt) == 1) // Si le point n'est pas éloigné de 1.5 des autres atomes
				{
					inASh++;
				}
				idCycle = SHL_addAtom(moc, positionNvDprt, -1);
				flag(atom(moc, idCycle)) = 4;
				SHL_addEdge(moc, idVoisin, idCycle);
				
				if ( i == 1 ) // Position du prochain depart pour continuer le chemin
				{
					posSuiv = positionNvDprt;
					idSuiv = idCycle;
				}
				
				idVoisin = idCycle;
			}
			
			SHL_addEdge(moc, id, idCycle);
			
			// Positionner atome qui suit le cycle
			v1 = coords(atom(moc, neighbor(atom(moc, idSuiv), 0)));
			Point_t v2 = coords(atom(moc, neighbor(atom(moc, idSuiv), 1)));
			positionNvDprt = AX2E1(posSuiv, v1, v2, SIMPLE); 
			if (encombrement(moc, sub, positionNvDprt) == 1) // Si le point n'est pas éloigné de 1.5 des autres atomes
			{
				inASh++;
			}
			int idSuiv2 = SHL_addAtom(moc, positionNvDprt, -1);
			flag(atom(moc, idSuiv2)) = 4;
			SHL_addEdge(moc, idSuiv, idSuiv2);
			
			if (inASh == 0) // Si tous les atomes du motif sont hors de l'enveloppe
			{
				LSTm_addElement(mocAtt, moc);
				LSTd_addElement(nvDepart, idSuiv2);
			}
			else
			{
				SHL_delete(moc);
			}
			
		}
	}
}

// Ajout seulement du 0 sur un C d'un motif 3
List_m* ajoutOMotif3(Shell_t* mocTraite, int depart, Molecule_t* sub) {
	
	List_m* mAtt = LSTm_init();
	Point_t dpt = coords(atom(mocTraite, depart));
	int voisin1 = neighbor(atom(mocTraite, depart), 0); // Voisin
	Point_t v1 = coords(atom(mocTraite, voisin1));
	
	for (int i = 0; i < neighborhoodSize(atom(mocTraite,voisin1)); i++) // Pour tous les plans possibles avec les voisins
	{
		if (neighbor(atom(mocTraite, voisin1), i) != -1 && neighbor(atom(mocTraite, voisin1), i) != depart)
		{
			Shell_t* moc = SHL_copy(mocTraite);
			Shell_t* moc2 = SHL_copy(mocTraite);
			
			Point_t v2 = coords(atom(moc, neighbor(atom(moc, voisin1), i)));
						
			// Cherche la normal pour positionné le O
			Point_t normal = planNormal(dpt, v1, v2);
			
			// Atome O
			//Premier position
			Point_t positionO = AX1E2(dpt, v1, normal, SIMPLE);
			
			if (encombrement(moc, sub, positionO) == 0) // Si le point est éloigné de 1.5 des autres atomes
			{
				int id3 = SHL_addAtom(moc, positionO, -1);
				flag(atom(moc, id3)) = 1;
				SHL_addEdge(moc, depart, id3);
				
				LSTm_addElement(mAtt, moc);
			}
			else
			{
				SHL_delete(moc);
			}
						
			//Seconde position
			positionO = AX2E1(dpt, v1, positionO, SIMPLE);
			
			if (encombrement(moc2, sub, positionO) == 0) // Si le point est éloigné de 1.5 des autres atomes
			{
				int id4 = SHL_addAtom(moc2, positionO, -1);
				flag(atom(moc2, id4)) = 1;
				SHL_addEdge(moc2, depart, id4);
				
				LSTm_addElement(mAtt, moc2);
			}
			else
			{
				SHL_delete(moc2);
			}
			
		}
	}
	
	return mAtt;
}

// Ajout du motif 3 (C = O) dans le plan avec son voisin
void ajoutMotif3(Shell_t* mocTraite, List_m* mocAtt, int depart, List_d* nvDepart, Point_t positionNvDprt, Molecule_t* sub) {
	
	for (int i = 0; i < neighborhoodSize(atom(mocTraite,depart)); i++) // Pour tous les plans possibles avec les voisins
	{
		if (neighbor(atom(mocTraite, depart), i) != -1)
		{
			Shell_t* moc = SHL_copy(mocTraite);
			Shell_t* moc2 = SHL_copy(mocTraite);
			
			Point_t v1 = coords(atom(moc, neighbor(atom(moc, depart), i)));
			Point_t dpt = coords(atom(moc, depart));
			
			// Atome C
			int id = SHL_addAtom(moc, positionNvDprt, -1);
			flag(atom(moc, id)) = 4;
			SHL_addEdge(moc, depart, id);
			
			int id2 = SHL_addAtom(moc2, positionNvDprt, -1);
			flag(atom(moc2, id2)) = 4;
			SHL_addEdge(moc2, depart, id2);
			
			// Cherche la normal pour positionné le O
			Point_t normal = planNormal(positionNvDprt, dpt, v1);
			
			// Atome O
			//Premier position
			Point_t positionO = AX1E2(positionNvDprt, dpt, normal, SIMPLE);
			
			if (encombrement(moc, sub, positionO) == 0) // Si le point est éloigné de 1.5 des autres atomes
			{
				int id3 = SHL_addAtom(moc, positionO, -1);
				flag(atom(moc, id3)) = 1;
				SHL_addEdge(moc, id, id3);
				
				LSTm_addElement(mocAtt, moc);
				LSTd_addElement(nvDepart, id);
			}
			else
			{
				SHL_delete(moc);
			}
			
			//Seconde position
			positionO = AX2E1(positionNvDprt, dpt, positionO, SIMPLE);
			
			if (encombrement(moc2, sub, positionO) == 0) // Si le point est éloigné de 1.5 des autres atomes
			{
				int id4 = SHL_addAtom(moc2, positionO, -1);
				flag(atom(moc2, id4)) = 1;
				SHL_addEdge(moc2, id2, id4);
				
				LSTm_addElement(mocAtt, moc2);
				LSTd_addElement(nvDepart, id2);
			}
			else
			{
				SHL_delete(moc2);
			}
			
		}
	}
}

// Ajout l'atome projeté a l'enveloppe
void ajoutProjection(Shell_t* mocTraite, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Point_t positionNvDprt, Molecule_t* sub) {
	
	if (numMotif == 3) // C = 0
	{
		ajoutMotif3(mocTraite, mocAtt, depart, nvDepart, positionNvDprt, sub);
	}
	else if (numMotif == 4)
	{
		ajoutMotif4(mocTraite, mocAtt, depart, nvDepart, positionNvDprt, sub);
	}
	else
	{
		Shell_t* moc = SHL_copy(mocTraite);
	
		int id = SHL_addAtom(moc, positionNvDprt, -1);
		flag(atom(moc, id)) = typeInsert(numMotif);
		SHL_addEdge(moc, depart, id);
		
		LSTm_addElement(mocAtt, moc);
		LSTd_addElement(nvDepart, id);
		
	}
	
}

/**************************************/
/* Projection emplacement *************/
/**************************************/

// Projection pour un atome avec 1 voisin
void projectionOCN_AX1E3(Shell_t* moc, List_m* mocAtt, int depart, int arrivee, List_d* nvDepart, int numMotif, Molecule_t* sub) {
	
	List_s* positions = LSTs_init();
	Point_t dpt = coords(atom(moc, depart));
	Point_t arv = coords(atom(moc, arrivee));
	Point_t v1 = coords(atom(moc, neighbor(atom(moc, depart), 0)));
	Point_t x2; // Voisin du voisin
	
	if (neighbor(atom(moc,neighbor(atom(moc,depart),0)),0) == depart)
		x2 = coords(atom(moc,neighbor(atom(moc,neighbor(atom(moc,depart),0)),1)));
	else
		x2 = coords(atom(moc,neighbor(atom(moc,neighbor(atom(moc,depart),0)),0)));

	Point_t normal = planNormal(dpt, v1, x2);
	
	Point_t positionNvDprt = AX1E3(dpt, v1, normal, SIMPLE);
	LSTs_addElement(positions, positionNvDprt);
	
	for (int i = 0; i < 11; i++) // Rotation a 360°
	{
		normal = rotation(normalization(vector(dpt, v1), 1),  30, normal); // Rotation de 30° de la normal
		positionNvDprt = AX1E3(dpt, v1, normal, SIMPLE);
		
		if (encombrement(moc, sub, positionNvDprt) == 0) // Si le point est éloigné de 1.5 des autres atomes
		{	
			LSTs_addElement(positions, positionNvDprt);
		}
	}
	
	for (int i = 0; i < 1 && positions->premier; i++) // 1 positions les mieux placés (distance min avec arrivée)
	{
		positionNvDprt = distMin(positions, arv); 
		LSTs_removeElement(positions, positionNvDprt);
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt, sub); // Ajout a l'enveloppe
	}
	
	LSTs_delete(positions);
}

// Projection pour un azote avec 2 voisins
void projectionN_AX2E2(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Molecule_t* sub) {
	
	Point_t positionNvDprt = AX2E2(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), SIMPLE);
	
	if (encombrement(moc, sub, positionNvDprt) == 0) // Si le point est éloigné de 1.5 des autres atomes
	{
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt, sub); // Ajout a l'enveloppe
	}
}

// Projection pour un carbone avec 2 voisins dont un oxygene
void projectionC_AX2E1(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Molecule_t* sub) {
	
	Point_t positionNvDprt = AX2E1(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), SIMPLE);
	
	if (encombrement(moc, sub, positionNvDprt) == 0) // Si le point est éloigné de 1.5 des autres atomes
	{
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt, sub); // Ajout a l'enveloppe
	}
}

// Projection pour un carbone avec 2 voisins
void projectionC_AX2E2(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, Molecule_t* sub) {
	
	Point_t positionNvDprt = AX2E2(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), SIMPLE);
	
	if (encombrement(moc, sub, positionNvDprt) == 0) // Si le point est éloigné de 1.5 des autres atomes
	{
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt, sub); // Ajout a l'enveloppe
	}
	
	Point_t positionNvDprt2 = AX3E1(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), positionNvDprt, SIMPLE);
	
	if (encombrement(moc, sub, positionNvDprt2) == 0) // Si le point est éloigné de 1.5 des autres atomes
	{
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt2, sub); // Ajout a l'enveloppe
	}
}

// Projection pour un carbone avec 3 voisins
void projectionC_AX3E1(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif,Molecule_t* sub) {
	
	Point_t positionNvDprt = AX3E1(coords(atom(moc, depart)), coords(atom(moc, neighbor(atom(moc, depart), 0))), coords(atom(moc, neighbor(atom(moc, depart), 1))), coords(atom(moc, neighbor(atom(moc, depart), 2))), SIMPLE);
	
	if (encombrement(moc, sub, positionNvDprt) == 0) // Si le point est éloigné de 1.5 des autres atomes
	{
		ajoutProjection(moc, mocAtt, depart, nvDepart, numMotif, positionNvDprt, sub); // Ajout a l'enveloppe
	}
}

/**************************************/
/* Générer chemin *********************/
/**************************************/

// Insertion du motif passé en argument
void insererMotif(Shell_t* moc, List_m* mocAtt, int depart, List_d* nvDepart, int numMotif, int arrivee, Molecule_t* sub){
	
	if ( LST_nbElements(neighborhood(atom(moc, depart))) == 1 ) // Oxygene ou Azote ou Carbone avec 1 voisin
	{
		//Projections
		//Diff rotations
		projectionOCN_AX1E3(moc, mocAtt, depart, arrivee, nvDepart, numMotif, sub);
	}
	else if (flag(atom(moc, depart)) == 3 && LST_nbElements(neighborhood(atom(moc, depart))) == 2) // Azote avec 2 voisins
	{
		//Projection
		projectionN_AX2E2(moc, mocAtt, depart, nvDepart, numMotif, sub);
	}
	else if (flag(atom(moc, depart)) == 4) // Carbone
	{
		if (LST_nbElements(neighborhood(atom(moc, depart))) == 2) // 2 voisins
		{
			if (flag(atom(moc, neighbor(atom(moc, depart), 0))) == 1 || flag(atom(moc, neighbor(atom(moc, depart), 1))) == 1) // Si 1 des 2 voisins est un oxygene
			{
				//Projection
				projectionC_AX2E1(moc, mocAtt, depart, nvDepart, numMotif, sub);
			}
			else
			{
				// 2 Projections
				projectionC_AX2E2(moc, mocAtt, depart, nvDepart, numMotif, sub);
			}
		}
		else // 3 voisins
		{
			//Projection
			projectionC_AX3E1(moc, mocAtt, depart, nvDepart, numMotif, sub);
		}
		
	}
}

// Génère le chemin entre 2 groupements de motifs
// Sans sommets intermediaires
void genererChemin(Main_t* m, List_m* mocAtt, Shell_t* mocTraite, int depart, int arrivee, int nbMotif3, int nbMotif4, char* InputFile, int tailleMax, int tailleMocDep){

/*************************************************/
/* Vérification distance interAtome***********/
/*************************************************/
	// Vérifier que la position d'ajout est éloigné de 1.5 de l'enveloppe
	Point_t B = coords(atom(mocTraite, depart));
	for(int i = 0; i < size(mocTraite)/*SHL_nbAtom(mocTraite)*/; i++){
		if(i!=depart && flag(atom(mocTraite, i))!=-1 && flag(atom(mocTraite, i))!=0 ){
		Point_t A = coords(atom(mocTraite, i));
			/*if(flag(atom(mocTraite, i))==0){
				if(dist(A, B) < DIST_ECART_ENVELOPPE) return;
			}
			else {*/
			if(dist(A, B) < 1/*DIST_ECART_CAGE*/) {
				return;
			}//}
		}
	}
	// Vérifier que la position d'ajout est éloigné de 2 du substrat
	for(int i = 0; i < size(substrat(m)); i++){
		Point_t A = coords(atom(substrat(m), i));
			if(dist(A, B) < 1/*DIST_ECART_SUBSTRAT*/) {
				return;
			}
	}
/*************************************************/

	for (int i = 2; i < NB_MOTIF; i++)
	{
		if(i==3) i++;
		List_m* lMoc = LSTm_init();
		List_d* nvDepart = LSTd_init();
		
		insererMotif(mocTraite, lMoc, depart, nvDepart, i, arrivee, substrat(m));
		
		while (lMoc->premier) // Pour toutes les solutions générées en générant le chemin / Diff rotations
		{
			// Compte le nombre de motif 3 d'affilée (C = O)
			if (i == 3)
			{
				nbMotif3++;
			}
			else
			{
				nbMotif3 = 0;
			}
			
			// Compte le nombre de motif 4 (Cycle) 
			if (i == 4) 
			{
				nbMotif4++;
			}
			
			if(tailleMax>=(SHL_nbAtom(lMoc->premier->moc)-tailleMocDep)){
				if (dist( coords(atom(lMoc->premier->moc, nvDepart->premier->sommet)), coords(atom(mocTraite, arrivee)) ) < MIN_DIST /*&& nbMotif4>0*/) // Proche de l'arrivée 
				{
					if(nbMotif4>0){//que s'il y a un cycle dans le chemin
						SHL_addEdge(lMoc->premier->moc, nvDepart->premier->sommet, arrivee); // Ajout lien entre dernier sommet du chemin et arrivee
						LSTm_addElement(mocAtt, SHL_copy(lMoc->premier->moc));// Ajout dans la liste a traiter
					}
				}
				else if (nbMotif3 < 5 && nbMotif4 < 3) // Maximum 4 motifs 3 d'affilée et 2 motifs 4 en tout
				{
					genererChemin(m, mocAtt, lMoc->premier->moc, nvDepart->premier->sommet, arrivee, nbMotif3, nbMotif4, InputFile, tailleMax, tailleMocDep);
				}
			}
			LSTm_removeFirst(lMoc);
			LSTd_removeFirst(nvDepart);
		}
		LSTm_delete(lMoc);
		LSTd_delete(nvDepart);
	}
}

/*************************************************/
/* Sommets départ et arrivée du chemin ***********/
/*************************************************/

// Vérifie si le sommet se situe en bordure de motif
int bordureCheck(Shell_t* s, AtomShl_t* sommet) {
	
	int vertexNumberOfNeighbors = LST_nbElements(neighborhood(sommet));
	for (int i = 0; i < vertexNumberOfNeighbors; i++) // Pour tous les voisins du sommet
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
	
	AtomShl_t* a = atom(s, indice1);
	LST_addElement(marquer, indice1);
	
	if (neighborhoodSize(a) == 0)
	{
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
					return 1;
				}
				else
				{
					if (!LST_check(marquer, neighbor(a, i))) // Si l'identifiant de ce sommet n'est pas déjà marqué
					{
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
	
	return sommets;
}

// Crée une liste de moc a traiter et vide le tableau de solutions finales
List_m* initMocAtt(Main_t* m){
	List_m* mocAtt = LSTm_init();
	
	for (int i=0; i<mocSize(m); i++) //Pour tous les mocs
	{
		if (moc(m,i) != NULL)
		{
			if (i == 0) // Traite juste le premier mocs
			{
				LSTm_addElement(mocAtt, moc(m, i)); // Les mettre dans la liste a traiter
			}
			if (i != 0)
			{
				SHL_delete(moc(m,i)); // Les supprimer du tableau de solutions finales
			}
		}
	}
	
	free(m->mocs);
	m->mocs = NULL;
	mocSize(m) = 0;
	
	return mocAtt;
}

/**************************************/
/* Fonction principale ****************/
/**************************************/

// Fonction principale
void assemblage(char* InputFile, Main_t* m, double alpha, int tailleMax){

	List_m* mocAtt = initMocAtt(m); // ! Prend le premier moc seulement
	int tailleMocInit = SHL_nbAtom(mocAtt->premier->moc); //permet de récupérer la taille avant l'ajout des chemins, fonctionne car on garde qu'un moc ligne d'avant (à modifier sinon)
	while (mocAtt->premier) // Tant qu'il existe un moc a traiter
	{	
		int tailleMocDep = SHL_nbAtom(mocAtt->premier->moc);
		List_p* sommets = choixSommets(mocAtt->premier->moc);
		
		if (!sommets->premier) // S'il n'y a plus qu'un groupement de motifs
		{
			
			//outputShell(InputFile, mocAtt->premier->moc);
			outputShell2(InputFile, mocAtt->premier->moc, tailleMocInit); // Ecriture de la sortie
			LSTm_removeFirst(mocAtt); // Suppression dans la liste a traiter
			
		}
		else // S'il y a au moins 2 groupements de motifs
		{
			Shell_t* mocTraite = mocAtt->premier->moc; // Copie le moc a traiter
			mocAtt->premier->moc = NULL;
			LSTm_removeFirst(mocAtt);
			
			// Quand on ne veut traiter qu'un couple de sommet, permet de changer de sommet à relier
			/*for (int i = 0; i < 30; i++)
			{
				LST2_removeFirst(sommets);
			}*/
			
			while (sommets->premier) // Pour tous les couples de sommets à relier
			{
				Shell_t* mocTraite2 = SHL_copy(mocTraite); // Crée un nouveau moc dans la liste a traiter
				int depart = sommets->premier->depart;
				int arrivee = sommets->premier->arrivee;
				LST2_removeFirst(sommets);

				for (int i = 0; i < LST_nbElements(neighborhood(atom(mocTraite2, depart))); i++) // Retire les voisins enveloppe de l'atome de départ (bordure)
				{
					if (flag(atom(mocTraite2, neighbor(atom(mocTraite2, depart), i))) == 0)
					{
						SHL_removeEdge(mocTraite2, depart, neighbor(atom(mocTraite2, depart), i));
						i--;
					}
				}
#pragma omp parallel for
				for (int i = 2; i < 4 ; i++) // Attribution de tous les types a l'atome de départ (sommet en bordure)
				{
					flag(atom(mocTraite2, depart)) = typeInsert(i);

					if (i == 3) 
					{
						if (LST_nbElements(neighborhood(atom(mocTraite2, depart))) == 1) // Motif possible que si le depart n'a qu'un voisin
						{
							List_m* mAtt = ajoutOMotif3(mocTraite2, depart,substrat(m)); // Ajout du O du motif 3 ( C = O )

							while (mAtt->premier) // Traiter tous les mocs générés par cet ajout
							{
								genererChemin(m, mocAtt, mAtt->premier->moc, depart, arrivee, 0, 0, InputFile, tailleMax, tailleMocDep);
								LSTm_removeFirst(mAtt);
							}
							LSTm_delete(mAtt);
						}
					}
					else
					{	
						genererChemin(m, mocAtt, mocTraite2, depart, arrivee, 0, 0, InputFile, tailleMax, tailleMocDep);

					}
					
				}
				
				SHL_delete(mocTraite2);
			}
			SHL_delete(mocTraite);
		}
		LST2_delete(sommets);
	}
	
	free(mocAtt);
}
