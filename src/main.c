#include "structure.h"
#include "initialization.h"
#include "expansion.h"
#include "generation.h"
#include "output.h"
#include "util.h"
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

	/********************************* Options *****/
	int opt;
  Options_t options = { NULL, DEFLT_ALPHA, DEFLT_SIZEMAX, DEFLT_MAX_RESULTS };

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

			case 'r':
        options.maxResults = atoi(optarg);
        break;

      case 'h':
      default:
        usage();
        break;
    }
	}

	if (options.input == NULL) {
		fprintf(stderr, "Input filename .xyz of the substrate is missing.\n");
		usage();
		exit(EXIT_FAILURE);
	}

	/*********************************** Infos *****/

	printf("\n####### Informations #######\n");
	printf("  - Substrate : %s\n  - Alpha : %.1f\n  - Maximum size of a path (in patterns) : %d\n  - Maximum number of results : %d\n",
					 options.input, options.alpha, options.sizeMax, options.maxResults);

	Main_t* m = MN_create();
	substrat(m) = initMolecule(options.input);

	/*************************************** R *****/

	printf("\n####### R environment initialization #######\n");
	
	int r_argc = 2;
	char *r_argv [] = {"R", "--silent"};

	setenv("R_HOME", "/usr/lib/R", 1);
	Rf_initEmbeddedR(r_argc, r_argv);

	setWorkingDirectory("src");
	source(PATHNAME);

	/*********** Envelope and binding patterns *****/

	envelope(m) = createShell(substrat(m), options.alpha);
	generatePathlessCages(m);

	Rf_endEmbeddedR(0);

	/***************************** Whole cages *****/

	writeMainOutput(options.input, m);
	
	generateWholeCages(m, options);
	
	MN_delete(m);

	/************************************ Time *****/
		
	time_t end = time(NULL);
	long seconds = (long) difftime(end, start);
		
	int hours = seconds / 3600;
	seconds -= hours * 3600;
	int minutes = seconds / 60;
	seconds -= minutes * 60;
	printf("\nExecution time : %d hour(s) %d minute(s) %ld second(s)\n", hours, minutes, seconds);
	
	return EXIT_SUCCESS;
}

void usage() {
	fprintf(stderr, USAGE_FMT, DEFLT_ALPHA, DEFLT_SIZEMAX, DEFLT_MAX_RESULTS);
	exit(EXIT_FAILURE);
}

void source(const char* name) {
	SEXP e;
	int errorOccurred;

	PROTECT(e = lang2(install("source"), mkString(name)));
  R_tryEval(e, R_GlobalEnv, &errorOccurred);

	if (errorOccurred) {
		fprintf(stderr, "An error occurred while using R.\n");
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
		fprintf(stderr, "An error occurred while using R.\n");
		exit(EXIT_FAILURE);
  }
  UNPROTECT(1);
}
