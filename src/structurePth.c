#include "structure.h"
#include "util.h"

/*  
Example - for 2 kept positions and 4 different orientations of a cycle. 
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
max positions | -1  |  0  |  1  |  1  |  0  | ... |
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
  //path->positionsBuffer = malloc((360 / ROTATION_ANGLE_AX1E3) * sizeof(Point_t)); // Unused for the moment. TODO use it to store positions computed in projectionAX1E3.
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

  Point_t end = coords(atom(moc, path->idEnd));
  Point_t center, firstNeighbor, secondNeighbor;
  Point_t hydrogen1, hydrogen2, hydroEnd1, hydroEnd2;
  Point_t startCycle, neighborStart, start, newStart, newNeighbor; // Variables for adding a cycle.
  Point_t normal, direction; // // Variables for adding a cycle.

  flag(atom(moc, path->idStart)) = CARBON_F;
  
  int idHydrogen, idNewStart;
  int idFirstNeighbor = neighbor(atom(moc, path->idStart), 0);
  int idCenter = path->idStart;
  
  firstNeighbor = path->patterns[0][path->positionNum[0]];
  center = path->patterns[1][path->positionNum[1]];

  for (int i = 1; i <= path->size; i++) {

    secondNeighbor = (i < path->size) ? path->patterns[i + 1][path->positionNum[i + 1]] : end;

    if (path->patternNum[i] == CYCLE_PATTERN) {
      startCycle = addPoint(firstNeighbor, normalization(vector(firstNeighbor, center), DIST_SIMPLE));
      neighborStart = path->patterns[i - 2][path->positionNum[i - 2]];
	    start = firstNeighbor;

		  int idstartCycle, idNeighborStart, idNewNeighbor;
		  int idEndCycle = -1;
      
	    normal = planNormal(startCycle, start, neighborStart);
	    normal = rotation(vector(startCycle, start), (path->orientations[i] + 2) * 45, normal);
		
		  idstartCycle = SHL_addAtom(moc, startCycle, -1);
		  flag(atom(moc, idstartCycle)) = CARBON_F;
		  SHL_addEdge(moc, idstartCycle, idFirstNeighbor);

		  neighborStart = startCycle;
		  idNeighborStart = idstartCycle;

		  direction = vector(startCycle, start);
		  newStart = addPoint(startCycle, normalization(rotation(normal, -120, direction), SIMPLE_CYCLE));

		  idNewStart = SHL_addAtom(moc, newStart, -1);
		  flag(atom(moc, idNewStart)) = CARBON_F;
		  SHL_addEdge(moc, idNewStart, idNeighborStart);

		  for (int j = 0; j < 4; j++) {

			  newNeighbor = AX1E2(newStart, neighborStart, normal, SIMPLE_CYCLE);
        idNewNeighbor = SHL_addAtom(moc, newNeighbor, -1);
			  flag(atom(moc, idNewNeighbor)) = CARBON_F;
			  SHL_addEdge(moc, idNewStart, idNewNeighbor);

        if (j != 2) {
          hydrogen1 = AX2E1(newStart, neighborStart, newNeighbor, DIST_ATOM_H);
          idHydrogen = SHL_addAtom(moc, hydrogen1, -1);
			    flag(atom(moc, idHydrogen)) = HYDROGEN_F;
          SHL_addEdge(moc, idNewStart, idHydrogen);
        }
        else {
          idEndCycle = idNewStart;
        }
			  neighborStart = newStart;
			  newStart = newNeighbor;

			  idNeighborStart = idNewStart;
			  idNewStart = idNewNeighbor;
		  }

		  hydrogen1 = AX2E1(newStart, neighborStart, startCycle, DIST_ATOM_H);
		
		  idHydrogen = SHL_addAtom(moc, hydrogen1, -1);
		  flag(atom(moc, idHydrogen)) = HYDROGEN_F;
		  SHL_addEdge(moc, idNewStart, idstartCycle);
		  SHL_addEdge(moc, idNewStart, idHydrogen);
      SHL_addEdge(moc, idEndCycle, idCenter);
    }

    hydrogen1 = AX2E2(center, firstNeighbor, secondNeighbor, DIST_ATOM_H);
		hydrogen2 = AX3E1(center, firstNeighbor, secondNeighbor, hydrogen1, DIST_ATOM_H);

		idHydrogen = SHL_addAtom(moc, hydrogen1, -1);
		flag(atom(moc, idHydrogen)) = HYDROGEN_F;
		SHL_addEdge(moc, idCenter, idHydrogen);
		idHydrogen = SHL_addAtom(moc, hydrogen2, -1);
		flag(atom(moc, idHydrogen)) = HYDROGEN_F;
		SHL_addEdge(moc, idCenter, idHydrogen);

    idNewStart = (i < path->size) ? SHL_addAtom(moc, secondNeighbor, -1) : path->idEnd;
		flag(atom(moc, idNewStart)) = CARBON_F;

    if (path->patternNum[i + 1] != CYCLE_PATTERN) {
      SHL_addEdge(moc, idCenter, idNewStart);
    }

    if (i < path->size) {
      idFirstNeighbor = idCenter;
      firstNeighbor = center;
      idCenter = idNewStart;
      center = secondNeighbor;
    }
  }

  hydrogen1 = AX2E2(center, firstNeighbor, end, DIST_ATOM_H);
	hydrogen2 = AX3E1(center, firstNeighbor, end, hydrogen1, DIST_ATOM_H);
	hydroEnd1 = AX2E2(end, coordsNeighbor(moc, path->idEnd, 0), center, DIST_ATOM_H);
	hydroEnd2 = AX3E1(end, coordsNeighbor(moc, path->idEnd, 0), center, hydroEnd1, DIST_ATOM_H);

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
	flag(atom(moc, path->idEnd)) = CARBON_F; // Change end atom (arrival) flag.
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