#ifndef __CONSTANTE_H
#define __CONSTANTE_H

//Main
#define OPTSTR "i:a:s:h"
#define USAGE_FMT  "%s usage : [-i inputfile] [-a alpha (default : %1.f)] [-s sizemax (default : %d)] [-h]\n"
#define DEFLT_ALPHA 3.
#define DEFLT_SIZEMAX 5

#define PATHNAME "alphashape.R"

//Structure
#define REALLOCSIZE 4

//distance
#define HYDRO 1.8 //taille liaison hydrogène
#define SIMPLE 1.5 //taille simple d'une liaison covalente
#define DIST_ERROR 0.5 //taille de l'erreur acceptable pour la liaison covalente
#define MINDIS 0.75
//#define MAX_DIST_ARRIVAL 2 //Distance pour relier notre chemin au sommet d'arrivée
#define DIST_GAP_CAGE 1.34//1.5  //distance entre 2 atomes de la cage
#define DIST_GAP_SUBSTRATE 1.8//2 //distance entre 1 atome de la cage et 1 atome du substrat


//Chemin
#define NB_MOTIF 5 //nombre de motifs pour les chemins

//Flag chemin
#define CARBONE 6
#define AZOTE 5
#define OXYGENE 4

#define NB_RESULTAT 10 //nombre de résultats


#endif