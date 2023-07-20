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
 * @brief Adds an aromatic ring (pattern 4) perpendicular to the plane 
 * with its neighbor.
 * 
 * @param processedMoc Molecular cage being generated.
 * @param mocsInProgress List of cages in construction to be processed.
 * @param idStart Index of the first linkable atom in the path in construction.
 * @param newStarts Stack of atoms starting a path in construction.
 * @param newStartPos Position (point) of the next atom added in the path.
 * @param sub Substrate molecule.
 */
void addAromaticRing(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, Point_t newStartPos, Molecule_t* sub, int nbPatterns, int nbCycles) {
	
	for (int i = 0; forEachNeighbor((atom(processedMoc,idStart)), i); i++) { // For every possible plans with starting atom's neighbors
		int idNext = -1;
		int idAtomCycle = -1;
		Point_t nextPos;
		Point_t hydrogen;
		Point_t hydrogen2;
		Point_t copyNewStartPos = newStartPos;
		Shell_t* moc = SHL_copy(processedMoc);
			
		Point_t neighborStartPos = coords(atom(moc, neighbor(atom(moc, idStart), i)));
		Point_t startPos = coords(atom(moc, idStart));
			
		// Add the first ring atom of already known position.
		int idnewStart = SHL_addAtom(moc, copyNewStartPos, -1);
		flag(atom(moc, idnewStart)) = CARBON_F;
		SHL_addEdge(moc, idStart, idnewStart);

		// Look for the normal to position the ring.
		Point_t normal = planNormal(copyNewStartPos, startPos, neighborStartPos);
		normal = rotation(normalization(vector(copyNewStartPos, startPos), 1),  90, normal); // Perpendicular
			
		// Position the other atoms of the cycle.
		neighborStartPos = AX1E2(copyNewStartPos, coords(atom(moc, idStart)), normal, SIMPLE_CYCLE); // Neighbor
		copyNewStartPos = AX2E1(copyNewStartPos, coords(atom(moc, idStart)), neighborStartPos, SIMPLE_CYCLE); 
		if (isHindered(moc, sub, copyNewStartPos)) {
			SHL_delete(moc);
			continue;
		}

		idAtomCycle = SHL_addAtom(moc, copyNewStartPos, -1);
		flag(atom(moc, idAtomCycle)) = CARBON_F;
		SHL_addEdge(moc, idnewStart, idAtomCycle);

		int idNeighbor = idAtomCycle;
		int wasHindered = 0;
		for (int i = 0; i < 4; i++) {
			neighborStartPos = coords(atom(moc, neighbor(atom(moc, idNeighbor), 0)));
			copyNewStartPos = AX1E2(copyNewStartPos, neighborStartPos, normal, SIMPLE_CYCLE);
			if (i != 2) {
				hydrogen = AX2E1(coords(atom(moc, idNeighbor)), neighborStartPos, copyNewStartPos, DIST_ATOM_H);
				if (isHindered(moc, sub, hydrogen)) {
					SHL_delete(moc);
					wasHindered = 1;
					break;
				}
			}
			if (i == 3) {
				hydrogen2 = AX2E1(copyNewStartPos, coords(atom(moc, idNeighbor)), coords(atom(moc, idnewStart)), DIST_ATOM_H);
				if (isHindered(moc, sub, hydrogen2)) {
					SHL_delete(moc);
					wasHindered = 1;
					break;
				}
			}
			if (isHindered(moc, sub, copyNewStartPos)) {
				SHL_delete(moc);
				wasHindered = 1;
				break;
			}
			idAtomCycle = SHL_addAtom(moc, copyNewStartPos, -1);
			flag(atom(moc, idAtomCycle)) = CARBON_F;
			SHL_addEdge(moc, idNeighbor, idAtomCycle);
			if (i != 2) { 
				int idHydrogen = SHL_addAtom(moc, hydrogen, -1);
				flag(atom(moc, idHydrogen)) = HYDROGEN_F;
				SHL_addEdge(moc, idNeighbor, idHydrogen);
			}
			if (i == 3) {
				int idHydrogen2 = SHL_addAtom(moc, hydrogen2, -1);
				flag(atom(moc, idHydrogen2)) = HYDROGEN_F;
				SHL_addEdge(moc, idAtomCycle, idHydrogen2);
			}

			if (i == 1) {// Position the next starting atom to continue the path.
				nextPos = copyNewStartPos;
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
		copyNewStartPos = AX2E1(nextPos, neighborStartPos, v2, DIST_SIMPLE); 
		if (isHindered(moc, sub, copyNewStartPos)) {
			SHL_delete(moc);
			continue;
		}
		int idSuiv2 = SHL_addAtom(moc, copyNewStartPos, -1);
		flag(atom(moc, idSuiv2)) = CARBON_F;
		SHL_addEdge(moc, idNext, idSuiv2);
			
		LSTm_addElement(mocsInProgress, moc);
		mocsInProgress->first->nbPatterns = nbPatterns + 1;
		mocsInProgress->first->nbCycles = nbCycles + 1;
		LSTd_addElement(newStarts, idSuiv2);
		return; // Keep only the first valid cycle inserted.
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
void projectionAX1E3(Shell_t* processedMoc, Path_t* path, Point_t** keptPositions, Molecule_t* sub) {
	
	int idBeforeStart = neighbor(atom(processedMoc, idStart), 0);
	int idFirstNeighbor = neighbor(atom(processedMoc,idBeforeStart), 0);
	int idSecondNeighbor = neighbor(atom(processedMoc,idBeforeStart), 1);

	Point_t startPos = coords(atom(processedMoc, idStart));
	Point_t beforeStartPos = coords(atom(processedMoc, idBeforeStart));
	Point_t beforeStartNeighborPos = (idFirstNeighbor == idStart) ? coords(atom(processedMoc, idSecondNeighbor)) : 
																	                                coords(atom(processedMoc, idFirstNeighbor));

	Point_t normal = planNormal(startPos, beforeStartPos, beforeStartNeighborPos);
	
	Point_t newStartPos = AX1E3(startPos, beforeStartPos, normal, DIST_SIMPLE);
	Point_t hydrogen1Pos = AX2E2(startPos, beforeStartPos, newStartPos, DIST_ATOM_H);
	Point_t hydrogen2Pos = AX3E1(startPos, beforeStartPos, newStartPos, hydrogen1Pos, DIST_ATOM_H);
	
	List_s* positions = LSTs_init();

	if (!isHindered(processedMoc, sub, newStartPos) && !isHindered(processedMoc, sub, hydrogen1Pos) && !isHindered(processedMoc, sub, hydrogen2Pos)) {
			LSTs_addElement(positions, newStartPos);
	}
	
	for (int i = 0; i < (360 / ROTATION_ANGLE_AX1E3) - 1; i++) { // 360° rotation.
		normal = rotation(normalization(vector(startPos, beforeStartPos), 1),  ROTATION_ANGLE_AX1E3, normal);
		newStartPos = AX1E3(startPos, beforeStartPos, normal, DIST_SIMPLE);
		hydrogen1Pos = AX2E2(startPos, beforeStartPos, newStartPos, DIST_ATOM_H);
		hydrogen2Pos = AX3E1(startPos, beforeStartPos, newStartPos, hydrogen1Pos, DIST_ATOM_H);
		
		if (!isHindered(processedMoc, sub, newStartPos) && !isHindered(processedMoc, sub, hydrogen1Pos) && !isHindered(processedMoc, sub, hydrogen2Pos)) {
			LSTs_addElement(positions, newStartPos); //TODO changer
		}
	}

	for (int i = 0; i < NUMBER_POSITION_AX1E3 && positions->first; i++) { // Best placed position (min distance to the end)
		//newStartPos = minDist_obstacle(positions, endPos,sub);
		newStartPos = minDist(positions, coords(atom(processedMoc, path->idEnd)));
		LSTs_removeElement(positions, newStartPos);
		path->positionsBuffer[i] = newStartPos;

		Shell_t* moc = SHL_copy(processedMoc);
		hydrogen1Pos = AX2E2(startPos, beforeStartPos, newStartPos, DIST_ATOM_H);
		hydrogen2Pos = AX3E1(startPos, beforeStartPos, newStartPos, hydrogen1Pos, DIST_ATOM_H);
		int idHydrogen = SHL_addAtom(moc, hydrogen1Pos, -1);
		flag(atom(moc, idHydrogen)) = HYDROGEN_F;
		SHL_addEdge(moc, idStart, idHydrogen);
		idHydrogen = SHL_addAtom(moc, hydrogen2Pos, -1);
		flag(atom(moc, idHydrogen)) = HYDROGEN_F;
		SHL_addEdge(moc, idStart, idHydrogen);
	}
	LSTs_delete(positions);
}

// Projection for an atom with two neighbors. Gives two projections.
void projectionAX2E2(Shell_t* processedMoc, Path_t* path, Point_t** keptPositions, Molecule_t* sub) {
	
	int idFirstNeighbor = neighbor(atom(processedMoc, idStart), 0);
	int idSecondNeighbor = neighbor(atom(processedMoc, idStart), 1);
	Point_t newStartPos = AX2E2(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), DIST_SIMPLE);

	if (!isHindered(processedMoc, sub, newStartPos)) {
		LSTs_addElement(positions, newStartPos);
	}
	
	Point_t newStartPos2 = AX3E1(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), newStartPos, DIST_SIMPLE);
	
	if (!isHindered(processedMoc, sub, newStartPos2)) {
		LSTs_addElement(positions, newStartPos2);
	}
}																																																																																																																																						

// Projection for an atom with 3 neighbors. Gives one projection.
void projectionAX3E1(Shell_t* processedMoc, Path_t* path, Point_t** keptPositions, Molecule_t* sub) {

	int idFirstNeighbor = neighbor(atom(processedMoc, idStart), 0);
	int idSecondNeighbor = neighbor(atom(processedMoc, idStart), 1);
	int idThirdNeighbor = neighbor(atom(processedMoc, idStart), 2);
	Point_t newStartPos = AX3E1(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), coords(atom(processedMoc, idThirdNeighbor)), DIST_SIMPLE);

	if (!isHindered(processedMoc, sub, newStartPos)) {
		LSTs_addElement(positions, newStartPos);
	}
}

/**************************************/
/********** Generate paths ************/
/**************************************/

/**
 * @brief Inserts the inputted pattern in a path.
 * 
 * @param processedMoc Molecular cage being generated.
 * @param mocsInProgress List of cages in construction to be processed.
 * @param idStart Index of the first linkable atom in the path in construction.
 * @param newStarts Stack of atoms starting a path in construction.
 * @param numPattern Pattern number (0, 1, 2) in the main loop. 
 * @param idEnd Index of the atom the path in construction is to be connected to.
 * @param sub Substrate molecule.
 */
void insertPattern(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, int numPattern, int idEnd, Molecule_t* sub, int nbPatterns, int nbCycles) {
	
	int startNumberOfNeighbors = LST_nbElements(neighborhood(atom(processedMoc, idStart)));
	List_s* positions = LSTs_init();
	
	if (startNumberOfNeighbors == 1) {
		Point_t startPos = coords(atom(processedMoc, idStart));
		Point_t beforeStartPos = coords(atom(processedMoc, neighbor(atom(processedMoc, idStart), 0)));
		Point_t newStartPos, hydrogen1Pos, hydrogen2Pos;

		projectionAX1E3(processedMoc, positions, idStart, sub);

		for (int i = 0; i < NUMBER_POSITION_AX1E3 && positions->first; i++) { // Best placed position (min distance to the end)
			//newStartPos = minDist_obstacle(positions, endPos,sub);
			newStartPos = minDist(positions, coords(atom(processedMoc, idEnd)));
		
			LSTs_removeElement(positions, newStartPos);
			Shell_t* moc = SHL_copy(processedMoc);
			hydrogen1Pos = AX2E2(startPos, beforeStartPos, newStartPos, DIST_ATOM_H);
			hydrogen2Pos = AX3E1(startPos, beforeStartPos, newStartPos, hydrogen1Pos, DIST_ATOM_H);
			int idHydrogen = SHL_addAtom(moc, hydrogen1Pos, -1);
			flag(atom(moc, idHydrogen)) = HYDROGEN_F;
			SHL_addEdge(moc, idStart, idHydrogen);
			idHydrogen = SHL_addAtom(moc, hydrogen2Pos, -1);
			flag(atom(moc, idHydrogen)) = HYDROGEN_F;
			SHL_addEdge(moc, idStart, idHydrogen);
		
			
		}
		LSTs_delete(positions);
		addProjection(moc, mocsInProgress, idStart, newStarts, numPattern, newStartPos, sub, nbPatterns, nbCycles);
	}
	else if (flag(atom(processedMoc, idStart)) == CARBON_F) {
		if (startNumberOfNeighbors == 2) {
			projectionAX2E2(processedMoc, positions, idStart, sub);
		}
		else { 
			projectionAX3E1(processedMoc, positions, idStart, sub);
		}
		while (positions->first) {
			Shell_t* moc = SHL_copy(processedMoc);
			addProjection(moc, mocsInProgress, idStart, newStarts, numPattern, positions->first->position, sub, nbPatterns, nbCycles);
			LSTs_removeFirst(positions);
		}
	}
}

Point_t** choosePositions(Shell_t* processedMoc int startNbNeighbors, Path_t* path, Molecule_t* substrate) {

	Point_t** keptPositions = keptPositions_init();
	if (startNbNeighbors == 1) {
		projectionAX1E3(processedMoc, path, keptPositions, substrate);
	}
	else if (flag(atom(processedMoc, idStart)) == CARBON_F) {
		if (startNbNeighbors == 2) {
			projectionAX2E2(processedMoc, path, keptPositions, substrate);
		}
		else {
			projectionAX3E1(processedMoc, path, keptPositions, substrate);
		}
		
	}
}

/**
 * @brief Recursively generates paths between two grouping of bonding patterns.
 * 
 * @param substrate Substrate molecule.
 * @param mocsInProgress Stack of cages in construction to be processed.
 * @param processedMoc Molecular cage being generated.
 * @param idStart Index of the first linkable atom in the path in construction.
 * @param idEnd Index of the atom the path in construction is to be connected to.
 * @param nbPatterns Number of patterns in the path to compare to the autorized number.
 * @param nbAroRings Number of aromatic rings recquiried.
 * @param inputFile Name of the substrate's file.
 * @param sizeMax Maximale size (in patterns) of a path.
 * @param forceCycle Forces the presence of a cycle in the path if true. 
 */
void generatePaths(Molecule_t* substrate, List_m* mocsInProgress, Shell_t* processedMoc, Path_t* path, 
									int nbPatterns, int nbAroRings, char* inputFile, int sizeMax, int forceCycle) {
	
	if (path->index < sizeMax - 1) {
		int startNbNeighbors = 1;
		Point_t** keptPositions = choosePositions(processedMoc, startNbNeighbors, path, substrate);
	}
	// Count the number of patterns.
	nbPatterns++;
	
	for (int i = 0; i < NB_PATTERNS; i++) {
		List_m* localMocsInProgress = LSTm_init();
		List_d* newStarts = LSTd_init();
		
		// Count the number of aromatic rings.
		if (i == CYCLE_PATTERN) {
			nbAroRings++;																																																		
		}
		insertPattern(processedMoc, localMocsInProgress, idStart, newStarts, i, path->idEnd, substrate, 0, 0);

		while (localMocsInProgress->first) {
			if(sizeMax >= nbPatterns && nbAroRings <= 3) {
				if (dist( coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), coords(atom(processedMoc, path->idEnd)) ) < DIST_SIMPLE + DIST_ERROR) {
					float beforeLastAngle = angle(coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)),coords(atom(localMocsInProgress->first->moc, path->idEnd)),coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, newStarts->first->idAtom),0))));
					float lastAngle = angle(coords(atom(localMocsInProgress->first->moc, path->idEnd)),coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)),coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, path->idEnd),0))));
					Point_t hydrogen1 = AX2E2(coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, newStarts->first->idAtom),0))), coords(atom(processedMoc, path->idEnd)), DIST_ATOM_H);
					Point_t hydrogen2 = AX3E1(coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, newStarts->first->idAtom),0))), coords(atom(processedMoc, path->idEnd)), hydrogen1, DIST_ATOM_H);
					Point_t hydrogenEnd1 = AX2E2(coords(atom(processedMoc, path->idEnd)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, path->idEnd),0))), coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), DIST_ATOM_H);
					Point_t hydrogenEnd2 = AX3E1(coords(atom(processedMoc, path->idEnd)), coords(atom(localMocsInProgress->first->moc,neighbor(atom(localMocsInProgress->first->moc, path->idEnd),0))), coords(atom(localMocsInProgress->first->moc, newStarts->first->idAtom)), hydrogenEnd1, DIST_ATOM_H);
				
					if (!isHindered(localMocsInProgress->first->moc, substrate, hydrogen1) && !isHindered(localMocsInProgress->first->moc, substrate, hydrogen2) && !isHindered(localMocsInProgress->first->moc, substrate, hydrogenEnd1) && !isHindered(localMocsInProgress->first->moc, substrate, hydrogenEnd2)) {
						if(beforeLastAngle >= (END_ANGLE - ANGLE_ERROR) && beforeLastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle >= (END_ANGLE - ANGLE_ERROR)) {
							if(!forceCycle || (forceCycle && nbAroRings > 0)) { // Only if there is a cycle in the path and we force the presence of a cycle.
								int idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogen1, -1);
								flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
								SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, idHydrogen);
								idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogen2, -1);
								flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
								SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, idHydrogen);
								idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogenEnd1, -1);
								flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
								SHL_addEdge(localMocsInProgress->first->moc, path->idEnd, idHydrogen);
								idHydrogen = SHL_addAtom(localMocsInProgress->first->moc, hydrogenEnd2, -1);
								flag(atom(localMocsInProgress->first->moc, idHydrogen)) = HYDROGEN_F;
								SHL_addEdge(localMocsInProgress->first->moc, path->idEnd, idHydrogen);
								flag(atom(localMocsInProgress->first->moc, path->idEnd)) = CARBON_F; // Change end atom (arrival) flag.
								SHL_addEdge(localMocsInProgress->first->moc, newStarts->first->idAtom, path->idEnd); //Add a link between last atom of the path and arrival.
								LSTm_addElement(mocsInProgress, SHL_copy(localMocsInProgress->first->moc));// Add to the list to be processed.
							}
						}
					}
				}
				else {
					generatePaths(substrate, mocsInProgress, localMocsInProgress->first->moc, path, 
												nbPatterns, nbAroRings, inputFile, sizeMax, forceCycle);
				}
			}
			LSTm_removeFirst(localMocsInProgress); // TODO test whether the cage was added, and delete if not the case. Otherwise, loose the pointer.
			LSTd_removeFirst(newStarts);
		}
		LSTm_delete(localMocsInProgress);
		LSTd_delete(newStarts);
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
 * @return (Element*) List of pair of atoms to be connected.
 */
Element* chooseStartAndEndPairs(Shell_t* s, Molecule_t* sub) {
	
	Element* startEndAtoms = LST_pairs_init();
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

		Element* startEndAtoms = chooseStartAndEndPairs(mocsInProgress->first->moc, substrat(m));
		Element* currentPair;
		
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
			#pragma omp parallel
			{
				currentPair = startEndAtoms;
				#pragma omp single
				{
					while (currentPair) { // For all pairs of atoms to connect.
					#pragma omp task firstprivate(currentPair)
					{
						Path_t* path = PTH_init(options.sizeMax, currentPair);
						PTH_fillFirstColumn(path, processedMoc);
						int forceCycle = 0;
						float startEndDist = dist(coords(atom(processedMoc, path->idStart)),coords(atom(processedMoc, path->idEnd)));
						if (startEndDist <= DIST_SIMPLE_PATTERN * options.sizeMax + DIST_SIMPLE + DIST_ERROR) {
							if (startEndDist <= DIST_SIMPLE_PATTERN * (options.sizeMax - 1/*NB_ATOMS_IN_CYCLE*/) + DIST_CYCLE_PATTERN + DIST_SIMPLE + DIST_ERROR
									&& startEndDist > DIST_CYCLE_PATTERN + (1 * DIST_SIMPLE_PATTERN) + DIST_SIMPLE + DIST_ERROR) {
										forceCycle = 1;
							}
							Shell_t* appendedMoc = SHL_copy(processedMoc); // Create a new moc in the list to process.
							flag(atom(appendedMoc, path->idStart)) = CARBON_F;
							generatePaths(substrat(m), mocsInProgress, appendedMoc, path, 0, 0, options.input, options.sizeMax, forceCycle);
						}
						PTH_delete(path, options.sizeMax);
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