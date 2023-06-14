#ifndef _MAIN_H
#define _MAIN_H

#define OPTSTR "i:a:s:h"
#define USAGE_FMT  "%s usage : [-i inputfile] [-a alpha (default : %1.f)] [-s sizemax (default : %d)] [-h]\n"
#define DEFLT_ALPHA 3.
#define DEFLT_SIZEMAX 13

#define PATHNAME "alphashape.R"

typedef struct {
  char*         input;
  double				alpha;
	int					sizeMax;
} options_t;

void usage(char *);
void source(const char*);
void setWorkingDirectory(char*);

#endif