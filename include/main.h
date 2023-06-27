#ifndef _MAIN_H
#define _MAIN_H

typedef struct {
  char*         input;
  double				alpha;
	int					sizeMax;
  int      maxResults;
} Options_t;

void usage();
void source(const char*);
void setWorkingDirectory(char*);

#endif