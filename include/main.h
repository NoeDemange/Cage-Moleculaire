#ifndef _MAIN_H
#define _MAIN_H

typedef struct {
  char*         input;
  double				alpha;
	int					sizeMax;
} options_t;

void usage(char *);
void source(const char*);
void setWorkingDirectory(char*);

#endif