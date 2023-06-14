#include "structure.h"
#include "initialization.h"
#include "expansion.h"
#include "generation.h"
#include "output.h"
#include "utile.h"
#include "main.h"
#include "assembly.h"

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <libgen.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>

int main(int argc, char** argv) {

	time_t start = time(NULL);

	int opt;
  options_t options = { NULL, DEFLT_ALPHA, DEFLT_SIZEMAX };

  while ((opt = getopt(argc, argv, OPTSTR)) != EOF) {
    switch(opt) {
      case 'i':
				options.input = optarg;
        break;

			case 'a':
        options.alpha = atof(optarg);
        break;

			case 's':
        options.sizeMax = atoi(optarg);
        break;

      case 'h':
      default:
        usage(argv[0]);
        break;
    }
	}

	if (options.input == NULL) {
		fprintf(stderr, "Fichier .xyz du substrat manquant.\n");
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	printf("\n####### Informations #######\n");
	printf("  - Substrat : %s\n  - Alpha : %.1f\n  - Taille maximum d'un chemin (en atomes) : %d\n",
					 options.input, options.alpha, options.sizeMax);

	Main_t* m = MN_create();
	substrat(m) = initMolecule(options.input);

	printf("\n####### Initialisation de l'environnement R  #######\n");
	
	int r_argc = 2;
	char *r_argv [] = {"R", "--silent"};

	setenv("R_HOME", "/usr/lib/R", 1);
	Rf_initEmbeddedR(r_argc, r_argv);

	setWorkingDirectory("src");
	source(PATHNAME);

	printf("\n####### Début de génération de l'enveloppe avec motifs liants #######\n");

	envelope(m) = createShell(substrat(m), options.alpha);
	generationMoc(m);

	Rf_endEmbeddedR(0);

	printf("\n####### Ecriture du substrat et de l'enveloppe dans le dossier de résultat #######\n");

	output(options.input, m);
	
	printf("\n####### Début de génération des chemins #######\n");
	assemblage(options.input, m, options.alpha, options.sizeMax);
	
	MN_delete(m);
		
	time_t end = time(NULL);
	long seconds = (long) difftime(end, start);
		
	int hours = seconds / 3600;
	seconds -= hours * 3600;
	int minuts = seconds / 60;
	seconds -= minuts * 60;
	printf("\nTemps d'execution : %d heure(s) %d minute(s) %ld seconde(s)\n", hours, minuts, seconds);
	
	return EXIT_SUCCESS;
}

void usage(char *argv0) {
	if (!argv0) {
    fprintf(stderr, "Argument manquant pour la fonction 'usage'.\n");
    exit(EXIT_FAILURE);
  }
	fprintf(stderr, USAGE_FMT, basename(argv0), DEFLT_ALPHA, DEFLT_SIZEMAX);
	exit(EXIT_FAILURE);
}

void source(const char* name) {
	SEXP e;
	int errorOccurred;

	PROTECT(e = lang2(install("source"), mkString(name)));
  R_tryEval(e, R_GlobalEnv, &errorOccurred);

	if (errorOccurred) {
		fprintf(stderr, "Une erreur est survenue lors de l'utilisation de R.\n");
		exit(EXIT_FAILURE);
  }
  UNPROTECT(1);
}

void setWorkingDirectory(char* dir) {
	SEXP e;
	int errorOccurred;

	PROTECT(e = lang2(install("setwd"), mkString(dir)));
  R_tryEval(e, R_GlobalEnv, &errorOccurred);

	if (errorOccurred) {
		fprintf(stderr, "Une erreur est survenue lors de l'utilisation de R.\n");
		exit(EXIT_FAILURE);
  }
  UNPROTECT(1);
}
