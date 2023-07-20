#include "structure.h"

/*         H     H      
            \___/       
   nei    _C/   \__       ...   H    H     
  /   \  /  \___/  \  /          \  /
nei1  start /   \    C            end
      /  \ H     H /  \
     H    H       H    H

              +-----+-----+-----+-----+--  -+-----+
numPattern    |  0  |  0  |  1  |  0  | ... |  0  |
              +-----+-----+-----+-----+--  -+-----+
              +-----+-----+-----+-----+--  -+-----+
patternsBuffer|     |  H  |  H  |  H  |     |  H  |
              | nei1|  H  |  H  |  H  | ... |  H  |
              | nei |start|  C  |  C  |     | end |
             ...             C1
                             C2
*/

Path_t* PTH_init(int size, Element* elem) {
  Path_t* path = malloc(sizeof(Path_t));
  //path->positionsBuffer = malloc((360 / ROTATION_ANGLE_AX1E3) * sizeof(Point_t)); // Unused for the moment. TODO make use of it.
  // 3 additional columns for neighbor of starting atom, starting atom and ending atom.
  path->patternsBuffer = malloc((size + 3) * sizeof(Point_t*));
	path->numPatternArray = malloc((size + 3) * sizeof(int));
	for (int i = 0; i < size; i++) {
		path->patternsBuffer[i] = malloc(MAX_ATOMS_NB_IN_PATTERN * sizeof(Point_t));
	}
  path->idStart = elem->start;
  path->idEnd = elem->end;
  path->index = 0;

  return path;
}

void PTH_delete(Path_t* path, int size) {
  if (path) {
    for (int i = 0; i < size; i++) {
		  free(path->patternsBuffer[i]);
	  }
	  free(path->patternsBuffer);
	  free(path->positionsBuffer);
	  free(path->numPatternArray);
  }
}

void PTH_fillFirstColumn(Path_t* path, Shell_t moc) {
  
}

Point_t** keptPositions_init() {

  int size = NUMBER_POSITION_AX1E3 < 2 ? 2 : NUMBER_POSITION_AX1E3;
  Point_t** table = malloc(size * sizeof(Point_t*));
  for (int i = 0; i < size; i++) {
		table[i] = malloc(size * sizeof(Point_t));
	}
  return table;
}

void keptPositions_delete(Point_t** table) {

  if (table) {
    int size = NUMBER_POSITION_AX1E3 < 2 ? 2 : NUMBER_POSITION_AX1E3;
    for (int i = 0; i < size; i++) {
		  free(table[i]);
	  }
    free(table);
  }
}