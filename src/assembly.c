#include "assembly.h"
#include "util.h"
#include "output.h"
#include "constant.h"
#include "voxelization.h"
#include <math.h>
#include <omp.h>

/**
 * @file assembly.c
 * @brief This file contains functions for generating connected cages.
 */

/**
 * @brief Checks if a point is far enough away from the other atoms 
 * of the cage and those of the substrate.
 * 
 * @param moc Molecular cage being generated.
 * @param substrate Substrate molecule.
 * @param path Path the atom is added to.
 * @param p  Point (atom) tested.
 * @return (int) 1 if not far enough, 0 otherwise.
 */
int isHindered(Shell_t* moc, Molecule_t* sub, Path_t* path, Point_t p) {

	for (int i = 0; i < size(moc); i++) {
		Point_t A = coords(atom(moc, i));
		if (distInf(A, p, DIST_GAP_CAGE))
			return 1;
	}
	for (int i = 0; i < size(sub); i++) {
		Point_t A = coords(atom(sub, i));
		if (distInf(A, p, DIST_GAP_SUBSTRATE)) 
			return 1;
	}

	for (int i = 2; i < path->size; i++) {
		int size = (path->patternNum[i] == CYCLE_PATTERN) ? MAX_NB_ATOMS_PATTERN : 3;
		for (int j = 0; j < size; j++) {
			if (distInf(path->patterns[i][path->positionNum[i]][j], p, DIST_GAP_CAGE))
				return 1;
		}
	}
	return 0;
}

/**
 * @brief Checks if a point is far enough away from the other atoms 
 * of the substrate.
 * 
 * @param sub Substrate molecule.
 * @param p  Point (atom) tested.
 * @return (int) 1 if not far enough, 0 otherwise.
 */
int isHinderedSubstrate(Molecule_t* sub, Point_t p) {

	for (int i = 0; i < size(sub); i++) {
		Point_t A = coords(atom(sub, i));
		if (dist(A, p) < DIST_GAP_SUBSTRATE) 
			return 1;
	}
	return 0;
}

/**************************************/
/********* Patterns addition **********/
/**************************************/

/**
 * @brief Adds an aromatic ring in the path (cycle pattern).
 * 
 * @param processedMoc Molecular cage being generated.
 * @param path Path being generated in the cage.
 * @param substrate Substrate molecule.
 */
int addAromaticRing(Shell_t* processedMoc, Path_t* path, Molecule_t* substrate) {

	Point_t neighborStart = path->patterns[path->size - 2][path->positionNum[path->size - 2]][0];
	Point_t start = path->patterns[path->size - 1][path->positionNum[path->size - 1]][0];
	Point_t startCycle = path->patterns[path->size][path->positionNum[path->size]][0];

	path->patterns[path->size][path->positionNum[path->size]][3] = startCycle;

	// Look for a normal to position the ring.
	Point_t normal = planNormal(startCycle, start, neighborStart);
	normal = rotation(vector(startCycle, start),  45, normal);

	for (int i = 0; i < 4; i++) { // Test 4 different orientations of the ring.
		
		int k = 4;
		Point_t newStart, newNeighbor, hydrogen, next;

		normal = rotation(vector(startCycle, start),  45, normal);

		Point_t direction = vector(startCycle, start);
		newStart = addPoint(startCycle, normalization(rotation(normal, -120, direction), SIMPLE_CYCLE));
		neighborStart = startCycle;

		if (isHindered(processedMoc, substrate, path, newStart)) {
			continue;
		}

		path->patterns[path->size][path->positionNum[path->size]][k++] = newStart;

		int wasHindered = 0;
		for (int j = 0; j < 4; j++) {

			newNeighbor = AX1E2(newStart, neighborStart, normal, SIMPLE_CYCLE);
			if (j != 2) {
				hydrogen = AX2E1(newStart, neighborStart, newNeighbor, DIST_ATOM_H);
			}
			else {
				hydrogen = AX2E1(newStart, neighborStart, newNeighbor, DIST_SIMPLE);
				next = hydrogen;
			}
			if (isHindered(processedMoc, substrate, path, newNeighbor) || isHindered(processedMoc, substrate, path, hydrogen)) {
				wasHindered = 1;
				break;
			}
			if (j != 2) {
				path->patterns[path->size][path->positionNum[path->size]][k++] = hydrogen;
			}
			path->patterns[path->size][path->positionNum[path->size]][k++] = newNeighbor;

			neighborStart = newStart;
			newStart = newNeighbor;
		}
		if (wasHindered) {
			continue;
		}
		hydrogen = AX2E1(newStart, neighborStart, startCycle, DIST_ATOM_H);
		
		if (isHindered(processedMoc, substrate, path, hydrogen)) {
			continue;
		}
		path->patterns[path->size][path->positionNum[path->size]][k] = hydrogen;
		path->patterns[path->size][path->positionNum[path->size]][0] = next;
		return 1; // We arbitrarily choose to keep only the first valid cycle inserted. 
	}
	return 0;
}

/**************************************/
/******** Projections location ********/
/**************************************/

// Projection for an atom with one neighbor. Gives multiple projections.
void projectionAX1E3(Shell_t* processedMoc, Path_t* path, Molecule_t* substrate, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap) {
	
	Point_t start = path->patterns[path->size - 1][path->positionNum[path->size - 1]][0];
	Point_t startNeighbor = path->patterns[path->size - 2][path->positionNum[path->size - 2]][0];
	// Also works for a cycle because atoms on either sides of the cycle are in the same plane.

	Point_t direction = vector(start, startNeighbor);
	Point_t normal = {-(direction.y + direction.z) / direction.x, 1, 1}; // Orthogonal vector to direction.
	
	Point_t newStart = AX1E3(start, startNeighbor, normal, DIST_SIMPLE);
	Point_t hydrogen1 = AX2E2(start, startNeighbor, newStart, DIST_ATOM_H);
	Point_t hydrogen2 = AX3E1(start, startNeighbor, newStart, hydrogen1, DIST_ATOM_H);
	
	List_s* positions = LSTs_init(); //TODO change for a buffer.

	if (!isHindered(processedMoc, substrate, path, newStart) && !isHindered(processedMoc, substrate, path, hydrogen1) && !isHindered(processedMoc, substrate, path, hydrogen2)) {
		//LSTs_addElement(positions, newStart);
			LSTs_addElementInOrder(positions, newStartPos, endPos, voxelGrid, vMap, nodeHeap);
	}
	
	for (int i = 0; i < (360 / ROTATION_ANGLE_AX1E3) - 1; i++) { // 360° rotation.
		normal = rotation(direction, ROTATION_ANGLE_AX1E3, normal);
		newStart = AX1E3(start, startNeighbor, normal, DIST_SIMPLE);
		hydrogen1 = AX2E2(start, startNeighbor, newStart, DIST_ATOM_H);
		hydrogen2 = AX3E1(start, startNeighbor, newStart, hydrogen1, DIST_ATOM_H);
		
		if (!isHindered(processedMoc, substrate, path, newStart) && !isHindered(processedMoc, substrate, path, hydrogen1) && !isHindered(processedMoc, substrate, path, hydrogen2)) {
			LSTs_addElement(positions, newStart); //path->positionsBuffer[i] = newStart;
			LSTs_addElementInOrder(positions, newStartPos, endPos, voxelGrid, vMap, nodeHeap);
		}
	}
	int i;
	for (i = 0; i < NUMBER_POSITION_AX1E3 && positions->first; i++) { // Best placed position (min distance to the end)
		//newStart = minDist_obstacle(positions, end,sub);
		newStart = positions->first->position;
		path->patterns[path->size][i][0] = newStart;
		hydrogen1 = AX2E2(start, startNeighbor, newStart, DIST_ATOM_H);
		hydrogen2 = AX3E1(start, startNeighbor, newStart, hydrogen1, DIST_ATOM_H);
		path->patterns[path->size][i][1] = hydrogen1;
		path->patterns[path->size][i][2] = hydrogen2;
		LSTs_removeFirst(positions);
	}
	path->maxPositions[path->size] = (i) ? i - 1 : -1;
	LSTs_delete(positions);
}

// Projection for an atom with two neighbors. Gives two projections.
void projectionAX2E2(Shell_t* processedMoc, Path_t* path, Molecule_t* sub) {
	
	// int idFirstNeighbor = neighbor(atom(processedMoc, idStart), 0);
	// int idSecondNeighbor = neighbor(atom(processedMoc, idStart), 1);
	// Point_t newStart = AX2E2(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), DIST_SIMPLE);

	// if (!isHindered(processedMoc, sub, newStart)) {
	// 	LSTs_addElement(positions, newStart);
	// }
	
	// Point_t newStartPos2 = AX3E1(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), newStart, DIST_SIMPLE);
	
	// if (!isHindered(processedMoc, sub, newStartPos2)) {
	// 	LSTs_addElement(positions, newStartPos2);
	// }
}																																																																																																																																						

// Projection for an atom with 3 neighbors. Gives one projection.
void projectionAX3E1(Shell_t* processedMoc, Path_t* path, Molecule_t* sub) {

	// int idFirstNeighbor = neighbor(atom(processedMoc, idStart), 0);
	// int idSecondNeighbor = neighbor(atom(processedMoc, idStart), 1);
	// int idThirdNeighbor = neighbor(atom(processedMoc, idStart), 2);
	// Point_t newStart = AX3E1(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), coords(atom(processedMoc, idThirdNeighbor)), DIST_SIMPLE);

	// if (!isHindered(processedMoc, sub, newStart)) {
	// 	LSTs_addElement(positions, newStart);
	// }
}

/**************************************/
/********** Generate paths ************/
/**************************************/


void choosePositions(Shell_t* processedMoc, Path_t* path, Molecule_t* substrate, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap) {

	if (path->size > 2) {
		projectionAX1E3(processedMoc, path, substrate);
	}
	else {
		int startNbNeighbors = LST_nbElements(neighborhood(atom(processedMoc, path->idStart)));
		if (startNbNeighbors == 1) {
			projectionAX1E3(processedMoc, path, substrate, voxelGrid, vMap, nodeHeap);
		}
		else if (startNbNeighbors == 2) {
			projectionAX2E2(processedMoc, path, substrate); // TODO check if we need to know the atom flag.
		}
		else {
			projectionAX3E1(processedMoc, path, substrate);
		}
	}
}

/**
 * @brief Generates paths between two grouping of bonding patterns.
 * 
 * @param substrate Substrate molecule.
 * @param mocsInProgress Stack of cages in construction to be processed.
 * @param processedMoc Molecular cage being generated.
 * @param path Structure of buffers to store the patterns of the path.
 * @param forceCycle Forces the presence of a cycle in the path if true. 
 */
void generatePaths(Molecule_t* substrate, mStack_t* mocsInProgress, Shell_t* processedMoc, Path_t* path, int forceCycle, VOXELGRID voxelGrid,  VMap*** vMap, NodeHeap nodeHeap) {
	
	int existsPathInProgress = 1;
	Point_t end = coords(atom(processedMoc, path->idEnd));

	while (existsPathInProgress) {

		int existsNewStart = (path->maxPositions[path->size] >= 0);

		if (path->size >= 2 && path->patternNum[path->size] == CYCLE_PATTERN) {
			existsNewStart = addAromaticRing(processedMoc, path, substrate);
		}
		else if (path->size <= path->sizeMax && !existsNewStart) {
			choosePositions(processedMoc, path, substrate);
			existsNewStart = (path->maxPositions[path->size] >= 0);
		}
		int closeToTheEnd = 0;
		if (existsNewStart) {
			
			Point_t newStart = path->patterns[path->size][path->positionNum[path->size]][0];

			if (dist(newStart, end) < DIST_SIMPLE + DIST_ERROR) {
				closeToTheEnd = 1;
				Point_t neighborNewStart = path->patterns[path->size - 1][path->positionNum[path->size - 1]][0];
				float beforeLastAngle = angle(newStart, end, neighborNewStart);
				float lastAngle = angle(end, newStart, coordsNeighbor(processedMoc, path->idEnd, 0));
				Point_t hydrogen1 = AX2E2(newStart, neighborNewStart, end, DIST_ATOM_H);
				Point_t hydrogen2 = AX3E1(newStart, neighborNewStart, end, hydrogen1, DIST_ATOM_H);
				Point_t hydroEnd1 = AX2E2(end, coordsNeighbor(processedMoc, path->idEnd, 0), newStart, DIST_ATOM_H);
				Point_t hydroEnd2 = AX3E1(end, coordsNeighbor(processedMoc, path->idEnd, 0), newStart, hydroEnd1, DIST_ATOM_H);
				
				if (!isHindered(processedMoc, substrate, path, hydrogen1) && !isHindered(processedMoc, substrate, path, hydrogen2) && !isHindered(processedMoc, substrate, path, hydroEnd1) && !isHindered(processedMoc, substrate, path, hydroEnd2)) {
					if(beforeLastAngle >= (END_ANGLE - ANGLE_ERROR) && beforeLastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle >= (END_ANGLE - ANGLE_ERROR)) {
						if(!forceCycle || (forceCycle && PTH_countAroRings(path) > 0)) { // Only if there is a cycle in the path and we force the presence of a cycle.
							Shell_t* moc = SHL_copy(processedMoc);
							PTH_addPath(moc, path);
							mSTK_addElement(mocsInProgress, moc);// Add to the list to be processed.
						}
					}
				}
			}
		}
		if (path->size == path->sizeMax || closeToTheEnd || !existsNewStart) {
			while (((path->patternNum[path->size] == CYCLE_PATTERN && path->positionNum[path->size] == path->maxPositions[path->size])
							|| path->maxPositions[path->size] < 0) && existsPathInProgress) {
				path->patternNum[path->size] = 0;
				path->positionNum[path->size] = 0;
				path->maxPositions[path->size] = -1;
				(path->size)--;

				if (path->size <= 1) {
					existsPathInProgress = 0;
				}
			}
			if (existsPathInProgress) {
				if (path->patternNum[path->size] == SIMPLE_PATTERN) {
					path->patternNum[path->size] = CYCLE_PATTERN;
				}
				else {
					(path->positionNum[path->size])++;
					path->patternNum[path->size] = 0;
				}
			}
		}
		else {
			(path->size)++;
		}
	}
}

/*************************************************/
/******* Start and end atoms of the paths ********/
/*************************************************/

/**
 * @brief Recursively does a depth-first search.
 * 
 * @param s Cage without any added paths.
 * @param markedAtoms List of previously visited atoms.
 * @param index1 Starting atom index.
 * @param index2 Searched atom index.
 * @return (int) 1 if the searched atom is found, 0 otherwise. 
 */
int search(Shell_t* s, List_t* markedAtoms, int index1, int index2) {
	
	AtomShl_t* a = atom(s, index1);
	LST_addElement(markedAtoms, index1);
	
	if (neighborhoodSize(a) == 0) {
		return 0;
	}
	else {
		for (int i = 0; forEachNeighbor(a, i); i++) {
			
			if (neighbor(a, i) == index2) {
				return 1;
			}
			else {
				if (!LST_check(markedAtoms, neighbor(a, i))) {
					if (search(s, markedAtoms, neighbor(a, i), index2)) {
						return 1;
					}
				}
			}
		}
	}
	return 0;
}

/**
 * @brief Checks if two atoms in the cage are connected.
 * Uses a Depth-first search.
 * 
 * @param s Cage without any added paths.
 * @param index1 Starting atom index.
 * @param index2 Ending atom index.
 * @return (int) 1 if they are connected, 0 otherwise.
 */
int checkExistsPath(Shell_t* s, int index1, int index2) {
	
	List_t* markedAtoms = NULL;
	markedAtoms = LST_create();
	
	int exists = search(s, markedAtoms, index1, index2);
	
	LST_delete(markedAtoms);
	
	return exists;
}

/**
 * @brief Generates a list of every pairs of atoms that can be linked.
 * They must belong to different connected components.
 * 
 * @param s Cage without any added paths.
 * @param sub Substrate molecule.
 * @param voxelGrid Grid of voxel.
 * @return (Pair_t*) List of pair of atoms to be connected.
 */
Pair_t* chooseStartAndEndPairs(Shell_t* s, Molecule_t* sub, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap) {
	
	Pair_t* startEndAtoms = PR_init();
	int* isHindered = malloc(size(s) * sizeof(int));

	for (int i = 0; i < size(s); i++) {
		if(flag(atom(s, i)) == LINKABLE_F && isHinderedSubstrate(sub, coords(atom(s, i)))) {
			isHindered[i] = 1;
		}
		else {
			isHindered[i] = 0;
		}
	}
	for (int i = 0; i < size(s) - 1; i++) {
		if (flag(atom(s, i)) == LINKABLE_F && !isHindered[i]) {
			for (int j = i + 1; j < size(s); j++) {
				if (flag(atom(s, j)) == LINKABLE_F && !isHindered[j]) {
					if (!checkExistsPath(s, i, j)) {
						LST_pairs_addElementInOrder(s, &startEndAtoms, i, j, voxelGrid, vMap, nodeHeap);
						//LST_pairs_addElement(&startEndAtoms, i, j);
					}
				}
			}
		}
	}
	free(isHindered);
	return startEndAtoms;
}

/**
 * @brief Initializes the list of moc with the first pathless cage generated.
 * Deletes the list of mocs of the main structure.
 * 
 * @param m Grouping of the main structures (substrate and envelope).
 * @return (mStack_t*) Initialized list of mocs. 
 */
mStack_t* initMocsInProgress(Main_t* m){
	
	mStack_t* mocsInProgress = mSTK_init();
	
	for (int i = 0; i < mocSize(m); i++) {
		if (moc(m,i) != NULL) {
			if (i == 0) { // Only the first moc. TODO? change it
				mSTK_addElement(mocsInProgress, SHL_copyCageAtoms(moc(m, i)));
			}
			if (i != 0) {
				SHL_delete(moc(m,i)); // delete from final solutions list.
			}
		}
	}
	free(m->mocs);
	m->mocs = NULL;
	mocSize(m) = 0;
	
	return mocsInProgress;
}

/**************************************/
/*********** Main function ************/
/**************************************/

/**
 * @brief Generates connected cages and writes them to the results directory.
 * 
 * @param m Grouping of the main structures (substrate and envelope).
 * @param options Grouping of inputfile, alpha, sizeMax, maxResults.
 */
void generateWholeCages(Main_t* m, Options_t options) {
	printf("\n####### Start of paths generation #######\n");
	VOXELGRID voxelGrid = voxelization(substrat(m));
	VMap*** vMap = VMap_alloc();
	NodeHeap nodeHeap = NH_initAlloc();
	mStack_t* mocsInProgress = initMocsInProgress(m); // ! Take only the first moc.
	static int countResults = 0;
	int pathelessMocSize = SHL_nbAtom(mocsInProgress->first->moc); // Allows to recover the size before the addition of the paths, only if we keep one moc line (TODO modify otherwise).
	while (mocsInProgress->first) { // As long as there is a moc to process.	

		Pair_t* startEndAtoms = chooseStartAndEndPairs(mocsInProgress->first->moc, substrat(m), voxelGrid, vMap, nodeHeap);
		Pair_t* currentPair;
		
		if (!startEndAtoms) { // If there is only one grouping of patterns left (connected cage).
			if (countResults++ < options.maxResults) {
				writeShellOutput(options.input, mocsInProgress->first->moc, pathelessMocSize);
				mSTK_removeFirst(mocsInProgress);
			}
			else {
				free(startEndAtoms);
				free(mocsInProgress);
				return;
			}
		}
		else { // If there are at least 2 groupings of patterns. 	
			Shell_t* processedMoc = mocsInProgress->first->moc;
			mocsInProgress->first->moc = NULL;
			mSTK_removeFirst(mocsInProgress);
			//#pragma omp parallel
			{
				currentPair = startEndAtoms;
				//#pragma omp single
				{
					while (currentPair) { // For all pairs of atoms to connect.
					//#pragma omp task firstprivate(currentPair)
					{
						Path_t* path = PTH_init(options.sizeMax, currentPair);
						path->patterns[0][0][0] = coordsNeighbor(processedMoc, path->idStart, 0);
						(path->size)++;
						path->patterns[1][0][0] = coords(atom(processedMoc, path->idStart));
						path->maxPositions[1] = 0;

						int forceCycle = 0;
						float startEndDist = dist(coords(atom(processedMoc, path->idStart)),coords(atom(processedMoc, path->idEnd)));
						if (startEndDist <= DIST_SIMPLE_PATTERN * options.sizeMax + DIST_SIMPLE + DIST_ERROR) {
							if (startEndDist <= DIST_SIMPLE_PATTERN * (options.sizeMax - 1/*NB_ATOMS_IN_CYCLE*/) + DIST_CYCLE_PATTERN + DIST_SIMPLE + DIST_ERROR
									&& startEndDist > DIST_CYCLE_PATTERN + (1 * DIST_SIMPLE_PATTERN) + DIST_SIMPLE + DIST_ERROR) {
										forceCycle = 1;
							}
							generatePaths(substrat(m), mocsInProgress, processedMoc, path, forceCycle, voxelGrid, vMap, nodeHeap);
						}
						PTH_delete(path);
					}
						currentPair = currentPair->next;
					}
				}
			}
			SHL_delete(processedMoc);
			PR_delete(startEndAtoms);
		}
	}
	freeVoxelGrid(voxelGrid);
	VMap_free(vMap);
	NH_free(nodeHeap);
	free(mocsInProgress);
}