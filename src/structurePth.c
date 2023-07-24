#include "structure.h"

/*  
Example - for 2 kept positions (pointers) and 4 different orientations of a cycle. 
             H     H      
              \___/       
   nei    C0__/   \__C1       ...   H    H     
      \  /    \___/    \  /          \  /
      start   /   \     C2            end
      /  \   H     H   /  \             \
     H    H           H    H            nei

              +-----+-----+-----+-----+-----+--  -+
orientations  | -1  | -1  | -1  |[0,3]| -1  | ... |
              +-----+-----+-----+-----+-----+--  -+
              +-----+-----+-----+-----+-----+--  -+
position num  |  0  |  0  |[0,1]|[0,1]|[0,1]|[0,1]|
              +-----+-----+-----+-----+-----+--  -+
              +-----+-----+-----+-----+-----+--  -+
 pattern num  |  0  |  0  |  0  |  1  |  0  |[0,1]|
              +-----+-----+-----+-----+-----+--  -+
              +-----+-----+-----+-----+-----+--  -+
patterns      | nei |start|C0.0 |C1.0 |C2.0 |     |
              +-----+-----|-----+-----+-----+--  -+
              |     |     |C0.1 |C1.1 |     |     |
              +-----+-----+-----+-----+-----+--  -+
*/

/**
 * @brief Create a path object.
 * 
 * @param size Maximum number of patterns in a path.
 * @param pair Pair of starting and ending atoms positions.
 * @return (Path_t*) Initialized path.
 */
Path_t* PTH_init(int size, Pair_t* pair) {
  Path_t* path = malloc(sizeof(Path_t));
  //path->positionsBuffer = malloc((360 / ROTATION_ANGLE_AX1E3) * sizeof(Point_t)); // Unused for the moment. TODO make use of it.
  // 2 additional columns for starting atom and its neighbor.
  int sizeMax = size + 2;
  path->patterns = malloc(sizeMax * sizeof(Point_t*));
	path->orientations = malloc(sizeMax * sizeof(int));
  path->positionNum = malloc(sizeMax * sizeof(int));
  path->patternNum = malloc(sizeMax * sizeof(int));
  path->maxPositions = malloc(sizeMax * sizeof(int));

  int sizePosition = (NUMBER_POSITION_AX1E3 < 2 ? 2 : NUMBER_POSITION_AX1E3);
	for (int i = 0; i < sizeMax; i++) {
		path->patterns[i] = malloc(sizePosition * sizeof(Point_t));
    for (int j = 0; j < sizePosition; j++) {
      path->patterns[i][j] = PT_init();
    }
    path->orientations[i] = -1;
    path->positionNum[i] = 0;
    path->patternNum[i] = 0;
    path->maxPositions[i] = -1;
	}
  path->idStart = pair->start;
  path->idEnd = pair->end;
  path->sizeMax = sizeMax - 1; // Because indexing starts at 0.
  path->size = 0;

  return path;
}

void PTH_delete(Path_t* path) {
  if (path) {
    for (int i = 0; i <= path->sizeMax; i++) {
		  free(path->patterns[i]);
	  }
	  free(path->patterns);
	  //free(path->positionsBuffer);
	  free(path->orientations);
    free(path->positionNum);
    free(path->patternNum);
  }
}

int PTH_countAroRings(Path_t* path) {
  int counter = 0;
  for (int i = 0; i <= path->sizeMax; i++) {
    if (path->patternNum[i] == CYCLE_PATTERN) {
      counter++;
    }
  }
  return counter;
}

void PTH_printTopBottom(int size) {
  printf("+");
  for (int i = 0; i < size; i++) {
    printf("--------+");
  }
   printf("\n");
}

void PTH_printPointers(int* array, int size, char* name) {
  
  printf("%s\n", name);
  PTH_printTopBottom(size);
  printf("|");
  for (int i = 0; i < size; i++) {
    if (array[i] < 0) {
      printf("  %d    |", array[i]);
    } else {
      printf("   %d    |", array[i]);
    }
  }
  printf("\n");
  PTH_printTopBottom(size);
}

void PTH_printPatterns(Path_t* path) {
  
  printf("Patterns\n");
  int size = (NUMBER_POSITION_AX1E3 < 2 ? 2 : NUMBER_POSITION_AX1E3);
  PTH_printTopBottom(path->sizeMax + 1);
  for (int i = 0; i < size; i++) {
    printf("|");
    for (int j = 0; j <= path->sizeMax; j++) {
      Point_t point = path->patterns[j][i];
      if (point.x < 0) {
        printf("%.0f,", point.x);
      } else {
        printf(" %.0f,", point.x);
      }
      if (point.y < 0) {
        printf("%.0f,", point.y);
      } else {
        printf(" %.0f,", point.y);
      }
      if (point.z < 0) {
        printf("%.0f|", point.z);
      } else {
        printf(" %.0f|", point.z);
      }
    }
    printf("\n");
    PTH_printTopBottom(path->sizeMax + 1);
  }
}

void PTH_printPath(Path_t* path) {
  PTH_printPointers(path->orientations, path->sizeMax + 1, "Orientations");
  PTH_printPointers(path->positionNum, path->sizeMax + 1, "Position number");
  PTH_printPointers(path->patternNum, path->sizeMax + 1, "Pattern number");
  PTH_printPointers(path->maxPositions, path->sizeMax + 1, "Max index of positions");
  PTH_printPatterns(path);
  printf("Path size : %d (max : %d)\n", path->size, path->sizeMax);
}