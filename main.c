#include "structure.h"
#include "initialization.h"
#include "expansion.h"
#include "generation.h"
#include "output.h"
#include "utile.h"

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>

#include <stdlib.h>

#define PATHNAME "alphashape.R"

void source (const char* name) {
	SEXP e;

	PROTECT(e = lang2(install("source"), mkString(name)));
  R_tryEval(e, R_GlobalEnv, NULL);
  UNPROTECT(1);
}

int main(int argc, char** argv) {
	/************ Initialisation de l'environnement R ************/
	//char * oldPath;
	int r_argc = 2;

	//oldPath = getenv ("R_HOME");
	setenv ("R_HOME", "/usr/lib/R", 1);

	char *r_argv [] = {"R", "--silent"};
	Rf_initEmbeddedR (r_argc, r_argv);

	source (PATHNAME);

	//Vérification qu'une entrée est passée en paramètre
	if (argc < 2) {
		printf("Veillez rentrer le nom du fichier de la molécule et/ou l'alpha\n");
		exit(1);
	}

	char* name = argv[1];
	double alpha = atof(argv[2]);

	Main_t* m = MN_create();

	substrat(m) = initMolecule(name);
	MOL_write(substrat(m));
	envelope(m) = createShell(substrat(m), alpha);
	SHL_write(envelope(m));
	printf("alpha = %0.1f, Nb sommets env = %d\n", alpha, SHL_nbAtom(envelope(m)));
	generationMoc(m);

	
	//SHL_write(envelope);
	//SHL_write(moc(m,0));

	/********** Écriture des résultats dans des fichiers **********/

	//Création du dossier de sortie.
	/*char* name = getBasename (name);
  char* dirName = createDir(name);
  char outputname[512];

  printf("Écriture de la molécule %s.\n", name);
  sprintf(outputname, "%s/%s.mol2", dirName, name);
  MOL_writeMol2(outputname, substrat(m));

  printf("Écriture de l'enveloppe.\n");
  sprintf(outputname, "%s/%s_shell.mol2", dirName, name);
  SHL_writeMol2(outputname, envelope(m));

  printf("Écriture de l'enveloppe après insertion des motifs aromatique.\n");
  sprintf(outputname, "%s/%s_aro.mol2", dirName, name);
  SHL_writeMol2(outputname, envarom(m));

  //Suppression des stuctures.
	free(name);
	free(dirName);*/

	output(name, m);

	

	MN_delete(m);

	/************** Fermeture de l'environnement R ***************/
	Rf_endEmbeddedR (0);
	//setenv ("R_HOME", oldPath, 1);

	return 0;
}
