#include "structure.h"
#include "initialization.h"
#include "expansion.h"
#include "generation.h"
#include "output.h"
#include "utile.h"
#include "assembly.h"

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>

#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define PATHNAME "alphashape.R"

void source (const char* name) {
	SEXP e;

	PROTECT(e = lang2(install("source"), mkString(name)));
  R_tryEval(e, R_GlobalEnv, NULL);
  UNPROTECT(1);
}

int main(int argc, char** argv) {
	
	time_t debut = time(NULL);
	
	/************ Initialisation de l'environnement R ************/
	//char * oldPath;
	int r_argc = 2;

	//oldPath = getenv ("R_HOME");
	setenv ("R_HOME", "/usr/lib/R", 1);

	char *r_argv [] = {"R", "--silent"};
	Rf_initEmbeddedR (r_argc, r_argv);

	source (PATHNAME);

	//Vérification qu'une entrée est passée en paramètre
	if (argc < 4) {
		printf("Veuillez rentrer le nom du fichier de la molécule et/ou l'alpha et/ou la taille maximale d'un chemin en nombre d'atomes et taille max du chemin\n");
		exit(1);
	}

	char* name = argv[1];
	double alpha = atof(argv[2]);
	int tailleMax = atoi(argv[3]);

	Main_t* m = MN_create();
	
	substrat(m) = initMolecule(name);
	MOL_write(substrat(m));
	envelope(m) = createShell(substrat(m), alpha);

	SHL_write(envelope(m));
	printf("alpha = %0.1f, Nb sommets env = %d\n", alpha, SHL_nbAtom(envelope(m)));
	generationMoc(m);
	
	/********** Assemblage des motifs **********/
	
	assemblage(name, m, alpha, tailleMax);

	/********** Écriture des résultats dans des fichiers **********/

	output(name, m);
	
	MN_delete(m);

	/************** Fermeture de l'environnement R ***************/
	Rf_endEmbeddedR (0);
	//setenv ("R_HOME", oldPath, 1);
		
	time_t fin = time(NULL);
	long secondes = (long) difftime(fin, debut);
		
	int heure = secondes / 3600;
	secondes -= heure * 3600;
	int minute = secondes / 60;
	secondes -= minute * 60;
	printf("Temps d'execution : %d heure(s) %d minute(s) %ld seconde(s)\n", heure, minute, secondes);
	
	return 0;
}
