#include "assembly.h"
#include "util.h"
#include "output.h"
#include "constant.h"
#include <math.h>
#include <omp.h>

/**
 * @brief Checks if a point is far enough away from the other atoms 
 * of the cage and those of the substrate.
 * 
 * @param moc Molecular cage being generated.
 * @param sub Substrate molecule.
 * @param p  Point (atom) tested.
 * @return (int) 1 if not far enough, 0 otherwise.
 */
int isHindered(Shell_t* moc, Molecule_t* sub, Point_t p) {

	for (int i = 0; i < size(moc); i++) {
			Point_t A = coords(atom(moc, i));
			if (dist(A, p) < DIST_GAP_CAGE)
				return 1;
	}
	for (int i = 0; i < size(sub); i++) {
		Point_t A = coords(atom(sub, i));
		if (dist(A, p) < DIST_GAP_SUBSTRATE) 
			return 1;
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
 * @brief Adds an aromatic ring.
 * 
 * @param processedMoc Molecular cage being generated.
 * @param path Path being generated in the cage.
 * @param sub Substrate molecule.
 */
void addAromaticRing(Shell_t* processedMoc, Path_t* path, Molecule_t* sub) {

	int index = path->size;
	Point_t neighborStartPos = path->patterns[index - 2][path->positionNum[index - 2]];
	Point_t startPos = path->patterns[index - 1][path->positionNum[index - 1]];
	Point_t startCyclePos = path->patterns[index][path->positionNum[index]];

	// Look for a normal to position the ring.
	Point_t normal = planNormal(startCyclePos, startPos, neighborStartPos);
	normal = rotation(vector(startCyclePos, startPos),  45, normal);

	for (int i = 0; i < 4; i++) { // Test 4 different orientations of the ring.
		int idNext = -1;
	  int idAtomCycle = -1;
	  Point_t nextPos, hydroPos, hydro2Pos;
	  //Shell_t* moc = SHL_copy(processedMoc);
		Point_t newStartPos = startCyclePos;
		normal = rotation(vector(newStartPos, startPos),  45, normal);

		// Add the first ring atom of already known position.
		// int idnewStart = SHL_addAtom(moc, newStartPos, -1);
		// flag(atom(moc, idnewStart)) = CARBON_F;
		// SHL_addEdge(moc, idStart, idnewStart);

		// Position the other atoms of the cycle.
		neighborStartPos = AX1E2(newStartPos, startPos, normal, SIMPLE_CYCLE);
		newStartPos = AX2E1(newStartPos, startPos, neighborStartPos, SIMPLE_CYCLE); 
		if (isHindered(processedMoc, sub, newStartPos)) {
			// SHL_delete(moc);
			continue;
		}

		idAtomCycle = SHL_addAtom(moc, newStartPos, -1);
		flag(atom(moc, idAtomCycle)) = CARBON_F;
		SHL_addEdge(moc, idnewStart, idAtomCycle);

		int idNeighbor = idAtomCycle;
		int wasHindered = 0;

		for (int j = 0; j < 4; j++) {
			neighborStartPos = coords(atom(moc, neighbor(atom(moc, idNeighbor), 0)));
			newStartPos = AX1E2(newStartPos, neighborStartPos, normal, SIMPLE_CYCLE);

			if (j != 2) {
				hydroPos = AX2E1(coords(atom(moc, idNeighbor)), neighborStartPos, newStartPos, DIST_ATOM_H);
				if (isHindered(processedMoc, sub, hydroPos)) {
					// SHL_delete(moc);
					wasHindered = 1;
					break;
				}
			}
			if (j == 3) {
				hydro2Pos = AX2E1(newStartPos, coords(atom(moc, idNeighbor)), coords(atom(moc, idnewStart)), DIST_ATOM_H);
				if (isHindered(processedMoc, sub, hydro2Pos)) {
					// SHL_delete(moc);
					wasHindered = 1;
					break;
				}
			}
			if (isHindered(processedMoc, sub, newStartPos)) {
				// SHL_delete(moc);
				wasHindered = 1;
				break;
			}
			idAtomCycle = SHL_addAtom(moc, newStartPos, -1);
			flag(atom(moc, idAtomCycle)) = CARBON_F;
			SHL_addEdge(moc, idNeighbor, idAtomCycle);

			if (j != 2) { 
				int idHydrogen = SHL_addAtom(moc, hydroPos, -1);
				flag(atom(moc, idHydrogen)) = HYDROGEN_F;
				SHL_addEdge(moc, idNeighbor, idHydrogen);
			}
			if (j == 3) {
				int idHydrogen2 = SHL_addAtom(moc, hydro2Pos, -1);
				flag(atom(moc, idHydrogen2)) = HYDROGEN_F;
				SHL_addEdge(moc, idAtomCycle, idHydrogen2);
			}

			if (j == 1) {// Position the next starting atom to continue the path.
				nextPos = newStartPos;
				idNext = idAtomCycle;
			}
			idNeighbor = idAtomCycle;
		}
		if (wasHindered) {
			continue;
		}

		SHL_addEdge(moc, idnewStart, idAtomCycle);
			
		// Position atom after the cycle.
		neighborStartPos = coords(atom(moc, neighbor(atom(moc, idNext), 0)));
		Point_t v2 = coords(atom(moc, neighbor(atom(moc, idNext), 1)));
		newStartPos = AX2E1(nextPos, neighborStartPos, v2, DIST_SIMPLE); 
		if (isHindered(processedMoc, sub, newStartPos)) {
			SHL_delete(moc);
			continue;
		}
		int idSuiv2 = SHL_addAtom(moc, newStartPos, -1);
		flag(atom(moc, idSuiv2)) = CARBON_F;
		SHL_addEdge(moc, idNext, idSuiv2);
		

		LSTm_addElement(mocsInProgress, moc);
		mocsInProgress->first->nbPatterns = nbPatterns + 1;
		mocsInProgress->first->nbCycles = nbCycles + 1;
		LSTd_addElement(newStarts, idSuiv2);
		return; // We arbitrarily choose to keep only the first valid cycle inserted. 
	}
}

/**
 * @brief Adds the projected atom(s) to the cage being generated.
 * 
 * @param processedMoc Molecular cage being generated.
 * @param mocsInProgress List of cages in construction to be processed.
 * @param idStart Index of the first linkable atom in the path in construction.
 * @param newStarts Stack of atoms starting a path in construction.
 * @param numPattern Pattern number (0, 1) in the main loop. 
 * @param newStartPos Position (point) of the next atom added in the path.
 * @param sub Substrate molecule.
 */
void addProjection(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, int numPattern, Point_t newStartPos, Molecule_t* sub, int nbPatterns, int nbCycles) {
	
	if (numPattern == CYCLE_PATTERN) {
		addAromaticRing(processedMoc, mocsInProgress, idStart, newStarts, newStartPos, sub, nbPatterns, nbCycles);
		SHL_delete(processedMoc);
	}
	else {
		int idnewStart = SHL_addAtom(processedMoc, newStartPos, -1);
		flag(atom(processedMoc, idnewStart)) = CARBON_F;
		SHL_addEdge(processedMoc, idStart, idnewStart);
		
		LSTm_addElement(mocsInProgress, processedMoc);
		mocsInProgress->first->nbPatterns = nbPatterns + 1;
		mocsInProgress->first->nbCycles = nbCycles;
		LSTd_addElement(newStarts, idnewStart);
	}
}

/**************************************/
/******** Projections location ********/
/**************************************/

// Projection for an atom with one neighbor. Gives multiple projections.
void projectionAX1E3(Shell_t* processedMoc, Path_t* path, Molecule_t* substrate) {
	
	int index = path->size;
	Point_t startPos = path->patterns[index - 1][path->positionNum[index - 1]];
	Point_t startNeighborPos = path->patterns[index - 2][path->positionNum[index - 2]];
	// Also works for a cycle because atoms on either sides of the cycle are in the same plane.

	Point_t direction = vector(startPos, startNeighborPos);
	Point_t normal = {-(direction.y + direction.z) / direction.x, 1, 1}; // Orthogonal vector to direction.
	
	Point_t newStartPos = AX1E3(startPos, startNeighborPos, normal, DIST_SIMPLE);
	Point_t hydro1Pos = AX2E2(startPos, startNeighborPos, newStartPos, DIST_ATOM_H);
	Point_t hydro2Pos = AX3E1(startPos, startNeighborPos, newStartPos, hydro1Pos, DIST_ATOM_H);
	
	List_s* positions = LSTs_init(); //TODO change for a buffer.

	if (!isHindered(processedMoc, substrate, newStartPos) && !isHindered(processedMoc, substrate, hydro1Pos) && !isHindered(processedMoc, substrate, hydro2Pos)) {
		LSTs_addElement(positions, newStartPos);
	}
	
	for (int i = 0; i < (360 / ROTATION_ANGLE_AX1E3) - 1; i++) { // 360° rotation.
		normal = rotation(direction, ROTATION_ANGLE_AX1E3, normal);
		newStartPos = AX1E3(startPos, startNeighborPos, normal, DIST_SIMPLE);
		hydro1Pos = AX2E2(startPos, startNeighborPos, newStartPos, DIST_ATOM_H);
		hydro2Pos = AX3E1(startPos, startNeighborPos, newStartPos, hydro1Pos, DIST_ATOM_H);
		
		if (!isHindered(processedMoc, substrate, newStartPos) && !isHindered(processedMoc, substrate, hydro1Pos) && !isHindered(processedMoc, substrate, hydro2Pos)) {
			LSTs_addElement(positions, newStartPos); //path->positionsBuffer[i] = newStartPos;
		}
	}
	int i;
	for (i = 0; i < NUMBER_POSITION_AX1E3 && positions->first; i++) { // Best placed position (min distance to the end)
		//newStartPos = minDist_obstacle(positions, endPos,sub);
		newStartPos = minDist(positions, coords(atom(processedMoc, path->idEnd)));
		path->patterns[index][i] = newStartPos;
		LSTs_removeElement(positions, newStartPos);
		// Shell_t* moc = SHL_copy(processedMoc);
		// hydrogen1Pos = AX2E2(startPos, startNeighborPos, newStartPos, DIST_ATOM_H);
		// hydrogen2Pos = AX3E1(startPos, startNeighborPos, newStartPos, hydrogen1Pos, DIST_ATOM_H);
		// int idHydrogen = SHL_addAtom(moc, hydrogen1Pos, -1);
		// flag(atom(moc, idHydrogen)) = HYDROGEN_F;
		// SHL_addEdge(moc, idStart, idHydrogen);
		// idHydrogen = SHL_addAtom(moc, hydrogen2Pos, -1);
		// flag(atom(moc, idHydrogen)) = HYDROGEN_F;
		// SHL_addEdge(moc, idStart, idHydrogen);
	}
	path->maxPositions[index] = (i) ? i - 1 : -1;
	LSTs_delete(positions);
}

// Projection for an atom with two neighbors. Gives two projections.
void projectionAX2E2(Shell_t* processedMoc, Path_t* path, Molecule_t* sub) {
	
	// int idFirstNeighbor = neighbor(atom(processedMoc, idStart), 0);
	// int idSecondNeighbor = neighbor(atom(processedMoc, idStart), 1);
	// Point_t newStartPos = AX2E2(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), DIST_SIMPLE);

	// if (!isHindered(processedMoc, sub, newStartPos)) {
	// 	LSTs_addElement(positions, newStartPos);
	// }
	
	// Point_t newStartPos2 = AX3E1(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), newStartPos, DIST_SIMPLE);
	
	// if (!isHindered(processedMoc, sub, newStartPos2)) {
	// 	LSTs_addElement(positions, newStartPos2);
	// }
}																																																																																																																																						

// Projection for an atom with 3 neighbors. Gives one projection.
void projectionAX3E1(Shell_t* processedMoc, Path_t* path, Molecule_t* sub) {

	// int idFirstNeighbor = neighbor(atom(processedMoc, idStart), 0);
	// int idSecondNeighbor = neighbor(atom(processedMoc, idStart), 1);
	// int idThirdNeighbor = neighbor(atom(processedMoc, idStart), 2);
	// Point_t newStartPos = AX3E1(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), coords(atom(processedMoc, idThirdNeighbor)), DIST_SIMPLE);

	// if (!isHindered(processedMoc, sub, newStartPos)) {
	// 	LSTs_addElement(positions, newStartPos);
	// }
}

/**************************************/
/********** Generate paths ************/
/**************************************/


void choosePositions(Shell_t* processedMoc, Path_t* path, Molecule_t* substrate) {

	if (path->size > 2) {
		projectionAX1E3(processedMoc, path, substrate);
	}
	else {
		int startNbNeighbors = LST_nbElements(neighborhood(atom(processedMoc, path->idStart)));
		if (startNbNeighbors == 1) {
			projectionAX1E3(processedMoc, path, substrate);
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
 * @brief Iteratively generates paths between two grouping of bonding patterns.
 * 
 * @param substrate Substrate molecule.
 * @param mocsInProgress Stack of cages in construction to be processed.
 * @param processedMoc Molecular cage being generated.
 * @param path Structure of buffers to store the patterns.
 * @param forceCycle Forces the presence of a cycle in the path if true. 
 */
void generatePaths(Molecule_t* substrate, List_m* mocsInProgress, Shell_t* processedMoc, Path_t* path, int forceCycle) {
	
	int pathsInProgress = 1;
	Point_t endPos = coords(atom(processedMoc, path->idEnd));

	int i = 0;
	while (pathsInProgress && i<5) {
		i++;
		int existsNewStart = (path->maxPositions[path->size] >= 0);
		if (path->patternNum[path->size] == CYCLE_PATTERN) {
			addAromaticRing(processedMoc, path, substrate);
			existsNewStart = (path->orientations[path->size] >= 0);
		}
		else if (path->size <= path->sizeMax && !existsNewStart) {
			choosePositions(processedMoc, path, substrate);
			existsNewStart = (path->maxPositions[path->size] >= 0);
		}
		int closeToTheEnd = 0;
		if (existsNewStart) {
			Point_t newStartPos = path->patterns[path->size][path->positionNum[path->size]];
			if (dist(newStartPos, endPos) < DIST_SIMPLE + DIST_ERROR) {
				closeToTheEnd = 1;
				Point_t neighborNewStartPos = path->patterns[path->size - 1][path->positionNum[path->size - 1]];
				float beforeLastAngle = angle(newStartPos, endPos, neighborNewStartPos);
				float lastAngle = angle(endPos, newStartPos, coordsNeighbor(processedMoc, path->idEnd, 0));
				Point_t hydro1Pos = AX2E2(newStartPos, neighborNewStartPos, endPos, DIST_ATOM_H);
				Point_t hydro2Pos = AX3E1(newStartPos, neighborNewStartPos, endPos, hydro1Pos, DIST_ATOM_H);
				Point_t hydroEnd1Pos = AX2E2(endPos, coordsNeighbor(processedMoc, path->idEnd, 0), newStartPos, DIST_ATOM_H);
				Point_t hydroEnd2Pos = AX3E1(endPos, coordsNeighbor(processedMoc, path->idEnd, 0), newStartPos, hydroEnd1Pos, DIST_ATOM_H);
				// Check hinderance.
				// Add path to a moc copy.
				// Add copy to mocsInProgress.
				if (!isHindered(processedMoc, substrate, hydro1Pos) && !isHindered(processedMoc, substrate, hydro2Pos) && !isHindered(processedMoc, substrate, hydroEnd1Pos) && !isHindered(processedMoc, substrate, hydroEnd2Pos)) {
					if(beforeLastAngle >= (END_ANGLE - ANGLE_ERROR) && beforeLastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle >= (END_ANGLE - ANGLE_ERROR)) {
						if(!forceCycle || (forceCycle && PTH_countAroRings(path) > 0)) { // Only if there is a cycle in the path and we force the presence of a cycle.
							// int idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydro1Pos, -1);
							// flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
							// SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, idHydrogen);
							// idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydro2Pos, -1);
							// flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
							// SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, idHydrogen);
							// idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydroEnd1Pos, -1);
							// flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
							// SHL_addEdge(localMocsInProgress->first->moc, path->idEnd, idHydrogen);
							// idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydroEnd2Pos, -1);
							// flag(atom(processedMoc, idHydrogen)) = HYDROGEN_F;
							// SHL_addEdge(localMocsInProgress->first->moc, path->idEnd, idHydrogen);
							// flag(atom(processedMoc, path->idEnd)) = CARBON_F; // Change end atom (arrival) flag.
							// SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, path->idEnd); //Add a link between last atom of the path and arrival.
							// LSTm_addElement(mocsInProgress, SHL_copy(localMocsInProgress->first->moc));// Add to the list to be processed.
						}
					}
				}
			}
		}
		if (path->size == path->sizeMax || closeToTheEnd || !existsNewStart) {
			while (((path->patternNum[path->size] == CYCLE_PATTERN && path->positionNum[path->size] == path->maxPositions[path->size])
							|| path->maxPositions[path->size] < 0) && pathsInProgress) {
				path->patternNum[path->size] = 0;
				path->orientations[path->size] = -1;
				path->positionNum[path->size] = 0;
				path->maxPositions[path->size] = -1;
				(path->size)--;

				if (path->size == 1) {
					pathsInProgress = 0;
				}
			}
			if (pathsInProgress) {
				if (path->patternNum[path->size] == SIMPLE_PATTERN) {
					path->patternNum[path->size] = CYCLE_PATTERN;
				}
				else {
					(path->positionNum[path->size])++;
					path->patternNum[path->size] = 0;
					path->orientations[path->size] = -1;
				}
			}
		}
		else {
			(path->size)++;
		}
		PTH_printPath(path);
		printf("\n");
	}
	exit(0);
	


	// Count the number of patterns.
	// nbPatterns++;
	
	// for (int i = 0; i < NB_PATTERNS; i++) {
	// 	List_m* localMocsInProgress = LSTm_init();
	// 	List_d* newStarts = LSTd_init();
		
	// 	// Count the number of aromatic rings.
	// 	if (i == CYCLE_PATTERN) {
	// 		nbAroRings++;																																																		
	// 	}
	// 	//insertPattern(processedMoc, localMocsInProgress, idStart, newStarts, i, path->idEnd, substrate, 0, 0);

	// 	while (localMocsInProgress->first) {
	// 		if(sizeMax >= nbPatterns && nbAroRings <= 3) {
	// 			if (dist( coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), coords(atom(processedMoc, path->idEnd)) ) < DIST_SIMPLE + DIST_ERROR) {
	// 				float beforeLastAngle = angle(coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)),coords(atom(localMocsInProgress->first->moc, path->idEnd)),coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, newStarts->first->idAtom),0))));
	// 				float lastAngle = angle(coords(atom(localMocsInProgress->first->moc, path->idEnd)),coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)),coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, path->idEnd),0))));
	// 				Point_t hydrogen1 = AX2E2(coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, newStarts->first->idAtom),0))), coords(atom(processedMoc, path->idEnd)), DIST_ATOM_H);
	// 				Point_t hydrogen2 = AX3E1(coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, newStarts->first->idAtom),0))), coords(atom(processedMoc, path->idEnd)), hydrogen1, DIST_ATOM_H);
	// 				Point_t hydrogenEnd1 = AX2E2(coords(atom(processedMoc, path->idEnd)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, path->idEnd),0))), coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), DIST_ATOM_H);
	// 				Point_t hydrogenEnd2 = AX3E1(coords(atom(processedMoc, path->idEnd)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, path->idEnd),0))), coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), hydrogenEnd1, DIST_ATOM_H);
				
	// 				if (!isHindered(localMocsInProgress->first->moc, substrate, hydrogen1) && !isHindered(localMocsInProgress->first->moc, substrate, hydrogen2) && !isHindered(localMocsInProgress->first->moc, substrate, hydrogenEnd1) && !isHindered(localMocsInProgress->first->moc, substrate, hydrogenEnd2)) {
	// 					if(beforeLastAngle >= (END_ANGLE - ANGLE_ERROR) && beforeLastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle >= (END_ANGLE - ANGLE_ERROR)) {
	// 						if(!forceCycle || (forceCycle && nbAroRings > 0)) { // Only if there is a cycle in the path and we force the presence of a cycle.
	// 							int idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogen1, -1);
	// 							flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
	// 							SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, idHydrogen);
	// 							idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogen2, -1);
	// 							flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
	// 							SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, idHydrogen);
	// 							idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogenEnd1, -1);
	// 							flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
	// 							SHL_addEdge(localMocsInProgress->first->moc, path->idEnd, idHydrogen);
	// 							idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogenEnd2, -1);
	// 							flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
	// 							SHL_addEdge(localMocsInProgress->first->moc, path->idEnd, idHydrogen);
	// 							flag(atom(localMocsInProgress->first->moc, path->idEnd)) = CARBON_F; // Change end atom (arrival) flag.
	// 							SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, path->idEnd); //Add a link between last atom of the path and arrival.
	// 							LSTm_addElement(mocsInProgress, SHL_copy(localMocsInProgress->first->moc));// Add to the list to be processed.
	// 						}
	// 					}
	// 				}
	// 			}
	// 			else {
	// 				generatePaths(substrate, mocsInProgress, localMocsInProgress->first->moc, path, 
	// 											nbPatterns, nbAroRings, inputFile, sizeMax, forceCycle);
	// 			}
	// 		}
	// 		LSTm_removeFirst(localMocsInProgress); // TODO test whether the cage was added, and delete if not the case. Otherwise, loose the pointer.
	// 		LSTd_removeFirst(newStarts);
	// 	}
	// 	LSTm_delete(localMocsInProgress);
	// 	LSTd_delete(newStarts);
	// }
}
// void generatePaths(Molecule_t* substrate, List_m* mocsInProgress, Shell_t* processedMoc, Path_t* path, int nbPatterns, int nbAroRings, int forceCycle) {
	
// 	if (path->size < path->sizeMax - 1) {
// 		int startNbNeighbors = 1;
// 		Point_t** keptPositions = choosePositions(processedMoc, startNbNeighbors, path, substrate);
// 	}
// 	// Count the number of patterns.
// 	nbPatterns++;
	
// 	for (int i = 0; i < NB_PATTERNS; i++) {
// 		List_m* localMocsInProgress = LSTm_init();
// 		List_d* newStarts = LSTd_init();
		
// 		// Count the number of aromatic rings.
// 		if (i == CYCLE_PATTERN) {
// 			nbAroRings++;																																																		
// 		}
// 		insertPattern(processedMoc, localMocsInProgress, idStart, newStarts, i, path->idEnd, substrate, 0, 0);

// 		while (localMocsInProgress->first) {
// 			if(sizeMax >= nbPatterns && nbAroRings <= 3) {
// 				if (dist( coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), coords(atom(processedMoc, path->idEnd)) ) < DIST_SIMPLE + DIST_ERROR) {
// 					float beforeLastAngle = angle(coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)),coords(atom(localMocsInProgress->first->moc, path->idEnd)),coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, newStarts->first->idAtom),0))));
// 					float lastAngle = angle(coords(atom(localMocsInProgress->first->moc, path->idEnd)),coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)),coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, path->idEnd),0))));
// 					Point_t hydrogen1 = AX2E2(coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, newStarts->first->idAtom),0))), coords(atom(processedMoc, path->idEnd)), DIST_ATOM_H);
// 					Point_t hydrogen2 = AX3E1(coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, newStarts->first->idAtom),0))), coords(atom(processedMoc, path->idEnd)), hydrogen1, DIST_ATOM_H);
// 					Point_t hydrogenEnd1 = AX2E2(coords(atom(processedMoc, path->idEnd)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, path->idEnd),0))), coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), DIST_ATOM_H);
// 					Point_t hydrogenEnd2 = AX3E1(coords(atom(processedMoc, path->idEnd)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, path->idEnd),0))), coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), hydrogenEnd1, DIST_ATOM_H);
				
// 					if (!isHindered(localMocsInProgress->first->moc, substrate, hydrogen1) && !isHindered(localMocsInProgress->first->moc, substrate, hydrogen2) && !isHindered(localMocsInProgress->first->moc, substrate, hydrogenEnd1) && !isHindered(localMocsInProgress->first->moc, substrate, hydrogenEnd2)) {
// 						if(beforeLastAngle >= (END_ANGLE - ANGLE_ERROR) && beforeLastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle >= (END_ANGLE - ANGLE_ERROR)) {
// 							if(!forceCycle || (forceCycle && nbAroRings > 0)) { // Only if there is a cycle in the path and we force the presence of a cycle.
// 								int idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogen1, -1);
// 								flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
// 								SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, idHydrogen);
// 								idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogen2, -1);
// 								flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
// 								SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, idHydrogen);
// 								idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogenEnd1, -1);
// 								flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
// 								SHL_addEdge(localMocsInProgress->first->moc, path->idEnd, idHydrogen);
// 								idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogenEnd2, -1);
// 								flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
// 								SHL_addEdge(localMocsInProgress->first->moc, path->idEnd, idHydrogen);
// 								flag(atom(localMocsInProgress->first->moc, path->idEnd)) = CARBON_F; // Change end atom (arrival) flag.
// 								SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, path->idEnd); //Add a link between last atom of the path and arrival.
// 								LSTm_addElement(mocsInProgress, SHL_copy(localMocsInProgress->first->moc));// Add to the list to be processed.
// 							}
// 						}
// 					}
// 				}
// 				else {
// 					generatePaths(substrate, mocsInProgress, localMocsInProgress->first->moc, path, 
// 												nbPatterns, nbAroRings, inputFile, sizeMax, forceCycle);
// 				}
// 			}
// 			LSTm_removeFirst(localMocsInProgress); // TODO test whether the cage was added, and delete if not the case. Otherwise, loose the pointer.
// 			LSTd_removeFirst(newStarts);
// 		}
// 		LSTm_delete(localMocsInProgress);
// 		LSTd_delete(newStarts);
// 	}
// }

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
		printf("here\n");
		return 0;
	}
	else {
		for (int i = 0; i < neighborhoodSize(a); i++) {

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
 * @return (Pair_t*) List of pair of atoms to be connected.
 */
Pair_t* chooseStartAndEndPairs(Shell_t* s, Molecule_t* sub) {
	
	Pair_t* startEndAtoms = LST_pairs_init();
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
						LST_pairs_addElementInOrder(s, &startEndAtoms, i, j);
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
 * @return (List_m*) Initialized list of mocs. 
 */
List_m* initMocsInProgress(Main_t* m){
	
	List_m* mocsInProgress = LSTm_init();
	
	for (int i = 0; i < mocSize(m); i++) {
		if (moc(m,i) != NULL) {
			if (i == 0) { // Only the first moc. TODO? change it
				LSTm_addElement(mocsInProgress, SHL_copyCageAtoms(moc(m, i)));
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
	List_m* mocsInProgress = initMocsInProgress(m); // ! Take only the first moc.
	static int countResults = 0;

	int pathelessMocSize = SHL_nbAtom(mocsInProgress->first->moc); // Allows to recover the size before the addition of the paths, only if we keep one moc line (TODO modify otherwise).
	while (mocsInProgress->first) { // As long as there is a moc to process.	

		Pair_t* startEndAtoms = chooseStartAndEndPairs(mocsInProgress->first->moc, substrat(m));
		Pair_t* currentPair;
		
		if (!startEndAtoms) { // If there is only one grouping of patterns left (connected cage).
			if (countResults++ < options.maxResults) {
				writeShellOutput(options.input, mocsInProgress->first->moc, pathelessMocSize);
				LSTm_removeFirst(mocsInProgress);
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
			LSTm_removeFirst(mocsInProgress);
			//#pragma omp parallel
			{
				currentPair = startEndAtoms;
				//#pragma omp single
				{
					while (currentPair) { // For all pairs of atoms to connect.
					//#pragma omp task firstprivate(currentPair)
					{
						Path_t* path = PTH_init(options.sizeMax, currentPair);
						path->patterns[0][0] = coordsNeighbor(processedMoc, path->idStart, 0);
						(path->size)++;
						path->patterns[1][0] = coords(atom(processedMoc, path->idStart));
						path->maxPositions[1] = 0;

						int forceCycle = 0;
						float startEndDist = dist(coords(atom(processedMoc, path->idStart)),coords(atom(processedMoc, path->idEnd)));
						if (startEndDist <= DIST_SIMPLE_PATTERN * options.sizeMax + DIST_SIMPLE + DIST_ERROR) {
							if (startEndDist <= DIST_SIMPLE_PATTERN * (options.sizeMax - 1/*NB_ATOMS_IN_CYCLE*/) + DIST_CYCLE_PATTERN + DIST_SIMPLE + DIST_ERROR
									&& startEndDist > DIST_CYCLE_PATTERN + (1 * DIST_SIMPLE_PATTERN) + DIST_SIMPLE + DIST_ERROR) {
										forceCycle = 1;
							}
							Shell_t* appendedMoc = SHL_copy(processedMoc); // Create a new moc in the list to process.
							flag(atom(appendedMoc, path->idStart)) = CARBON_F;
							generatePaths(substrat(m), mocsInProgress, appendedMoc, path, forceCycle);
						}
						PTH_delete(path);
					}
						currentPair = currentPair->next;
					}
				}
			}
			SHL_delete(processedMoc);
			LST_pairs_delete(startEndAtoms);
		}
	}
	free(mocsInProgress);
}