#include "structure.h"
#include "util.h"

/*  
Example - for 2 kept positions.                  atoms in cycle number

       H3 H4 H5   H6 H9 H10                          H5      H6     
         \/   \___/   \/                              \     /
   nei    C0__/   \__C1       ...   H    H            Cb___Cc       
      \  /    \___/    \  /          \  /            /       \
      start   /   \     C2            end       C0__Ca       Cd__C1
      /  \   H7    H8  /  \             \            \       / 
     H1   H2          H11  H12          nei           Cf___Ce
                                                      /     \
                                                    H7      H8

              +-----+-----+-----+-----+-----+--  -+
position num  |  0  |  0  |[0,1]|[0,1]|[0,1]|[0,1]|
              +-----+-----+-----+-----+-----+--  -+
              +-----+-----+-----+-----+-----+--  -+
 pattern num  |  0  |  0  |  0  |  1  |  0  |[0,1]|
              +-----+-----+-----+-----+-----+--  -+
              +-----+-----+-----+-----+-----+--  -+
max positions | -1  |  0  |  1  |  1  |  0  | -1  |
              +-----+-----+-----+-----+-----+--  -+              
              +-----+-----+-----+-----+-----+--  -+
patterns      | nei |start|C0.0 |C1.0 |C2.0 |     |
              |     |     |H1.0 |H3.0 |H9.0 |     |
              |     |     |H2.0 |H4.0 |H10.0|     |
              |     |     |     |Ca.0 |     |     |
              |     |     |     |Cb.0 |     |     |
              |     |     |     |H5.0 |     |     |
              |     |     |     |Cc.0 |     | ... |
              |     |     |     |H6.0 |     |     |
              |     |     |     |Cd.0 |     |     |
              |     |     |     |Ce.0 |     |     |
              |     |     |     |H8.0 |     |     |
              |     |     |     |Cf.0 |     |     |
              |     |     |     |H7.0 |     |     |
              +-----+-----|-----+-----+-----+--  -+
              |     |     |C0.1 |C1.1 |     |     |
              |     |     |H1.1 |h3.1 |     |     |
              |     |     |H2.1 |h4.1 |     |     |
                               ...
              |     |     |     |     |     |     |
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
  //path->positionsBuffer = malloc((360 / ROTATION_ANGLE_AX1E3) * sizeof(Point_t)); // Unused for the moment. TODO use it to store positions computed in projectionAX1E3.
  // 2 additional columns for starting atom and its neighbor.
  int sizeMax = size + 2;
  path->patterns = malloc(sizeMax * sizeof(Point_t**));
  path->positionNum = malloc(sizeMax * sizeof(int));
  path->patternNum = malloc(sizeMax * sizeof(int));
  path->maxPositions = malloc(sizeMax * sizeof(int));

  int sizePosition = (NUMBER_POSITION_AX1E3 < 2 ? 2 : NUMBER_POSITION_AX1E3);
	for (int i = 0; i < sizeMax; i++) {
		path->patterns[i] = malloc(sizePosition * sizeof(Point_t*));
    for (int j = 0; j < sizePosition; j++) {
      path->patterns[i][j] = malloc(MAX_NB_ATOMS_PATTERN * sizeof(Point_t));
      for (int k = 0; k <  MAX_NB_ATOMS_PATTERN; k++) {
        path->patterns[i][j][k] = PT_init();
      }
    }
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
    int sizePosition = (NUMBER_POSITION_AX1E3 < 2 ? 2 : NUMBER_POSITION_AX1E3);

    for (int i = 0; i <= path->sizeMax; i++) {
		  for (int j = 0; j < sizePosition; j++) {
        free(path->patterns[i][j]);
      }
      free(path->patterns[i]);
	  }
	  free(path->patterns);
	  //free(path->positionsBuffer);
    free(path->positionNum);
    free(path->patternNum);
    free(path->maxPositions);
    free(path);
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


void PTH_addPath(Shell_t* moc, Path_t* path) {

  flag(atom(moc, path->idStart)) = CARBON_F;
  
  int idHydrogen, idNewStart;
  int idCenter = path->idStart;

  for (int i = 2; i <= path->size; i++) {

    idHydrogen = SHL_addAtom(moc, path->patterns[i][path->positionNum[i]][1], -1);
    flag(atom(moc, idHydrogen)) = HYDROGEN_F;
    SHL_addEdge(moc, idCenter, idHydrogen);
    idHydrogen = SHL_addAtom(moc, path->patterns[i][path->positionNum[i]][2], -1);
    flag(atom(moc, idHydrogen)) = HYDROGEN_F;
    SHL_addEdge(moc, idCenter, idHydrogen);

    if (path->patternNum[i] == CYCLE_PATTERN) {

      int idStartCycle, idNeighbor;

      idStartCycle = SHL_addAtom(moc, path->patterns[i][path->positionNum[i]][3], -1);
      flag(atom(moc, idStartCycle)) = CARBON_F;
      SHL_addEdge(moc, idCenter, idStartCycle);

      idNewStart = SHL_addAtom(moc, path->patterns[i][path->positionNum[i]][4], -1);
      flag(atom(moc, idNewStart)) = CARBON_F;
      SHL_addEdge(moc, idStartCycle, idNewStart);

      int j = 5;
      while (j < 11) {
        if (j != 9) {
          idHydrogen = SHL_addAtom(moc, path->patterns[i][path->positionNum[i]][j++], -1);
          flag(atom(moc, idHydrogen)) = HYDROGEN_F;
          SHL_addEdge(moc, idNewStart, idHydrogen);
        }
        else {
          idCenter = SHL_addAtom(moc, path->patterns[i][path->positionNum[i]][0], -1);
          flag(atom(moc, idCenter)) = CARBON_F;
          SHL_addEdge(moc, idNewStart, idCenter);
        }
        idNeighbor = SHL_addAtom(moc, path->patterns[i][path->positionNum[i]][j++], -1);
        flag(atom(moc, idNeighbor)) = CARBON_F;
        SHL_addEdge(moc, idNewStart, idNeighbor);

        idNewStart = idNeighbor;
      }

      idHydrogen = SHL_addAtom(moc, path->patterns[i][path->positionNum[i]][12], -1);
      flag(atom(moc, idHydrogen)) = HYDROGEN_F;
      SHL_addEdge(moc, idNewStart, idHydrogen);
      SHL_addEdge(moc, idNewStart, idStartCycle);
    }
    else {
      idNewStart = SHL_addAtom(moc, path->patterns[i][path->positionNum[i]][0], -1);
      flag(atom(moc, idNewStart)) = CARBON_F;
      SHL_addEdge(moc, idCenter, idNewStart);

      idCenter = idNewStart;
    }
  }
  SHL_addEdge(moc, idCenter, path->idEnd);

  Point_t hydrogen1, hydrogen2, hydroEnd1, hydroEnd2;
  Point_t center = path->patterns[path->size][path->positionNum[path->size]][0];
  Point_t firstNeighbor = path->patterns[path->size - 1][path->positionNum[path->size - 1]][0];

  hydrogen1 = AX2E2(center, firstNeighbor, coords(atom(moc, path->idEnd)), DIST_ATOM_H);
	hydrogen2 = AX3E1(center, firstNeighbor, coords(atom(moc, path->idEnd)), hydrogen1, DIST_ATOM_H);
	hydroEnd1 = AX2E2(coords(atom(moc, path->idEnd)), coordsNeighbor(moc, path->idEnd, 0), center, DIST_ATOM_H);
	hydroEnd2 = AX3E1(coords(atom(moc, path->idEnd)), coordsNeighbor(moc, path->idEnd, 0), center, hydroEnd1, DIST_ATOM_H);

  idHydrogen = SHL_addAtom(moc, hydrogen1, -1);
	flag(atom(moc, idHydrogen)) = HYDROGEN_F;
  SHL_addEdge(moc, idCenter, idHydrogen);
	idHydrogen = SHL_addAtom(moc, hydrogen2, -1);
	flag(atom(moc, idHydrogen)) = HYDROGEN_F;
	SHL_addEdge(moc, idCenter, idHydrogen);
	idHydrogen = SHL_addAtom(moc, hydroEnd1, -1);
	flag(atom(moc, idHydrogen)) = HYDROGEN_F;
	SHL_addEdge(moc, path->idEnd, idHydrogen);
	idHydrogen = SHL_addAtom(moc, hydroEnd2, -1);
	flag(atom(moc, idHydrogen)) = HYDROGEN_F;
	SHL_addEdge(moc, path->idEnd, idHydrogen);
	flag(atom(moc, path->idEnd)) = CARBON_F;
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
    for (int k = 0; k < MAX_NB_ATOMS_PATTERN; k++) {
      printf("|");
      for (int j = 0; j <= path->sizeMax; j++) {

        Point_t point = path->patterns[j][i][k];
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
    }
    PTH_printTopBottom(path->sizeMax + 1);
  }
}

void PTH_printPath(Path_t* path) {

  PTH_printPointers(path->positionNum, path->sizeMax + 1, "Position number");
  PTH_printPointers(path->patternNum, path->sizeMax + 1, "Pattern number");
  PTH_printPointers(path->maxPositions, path->sizeMax + 1, "Max index of positions");
  PTH_printPatterns(path);
  printf("Path size : %d (max : %d)\n", path->size, path->sizeMax);
}