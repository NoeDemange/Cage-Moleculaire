#include "assembly.h"
#include "util.h"
#include "output.h"
#include "constant.h"
#include <math.h>

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

/**************************************/
/********* Patterns addition **********/
/**************************************/

/**
 * @brief Returns the type of the inserted atom.
 * 
 * @param numPattern Pattern number (0, 1, 2) in the main loop.
 * @return (int) Number corresponding to the pattern type.
 */
int insertType(int numPattern) {

	if (numPattern == 0) {
		return OXYGEN;
	}
	else if (numPattern == 1) {
		return NITROGEN;
	}
	else {
		return CARBON;
	}
}

/**
 * @brief Adds an aromatic ring (pattern 4) perpendicular to the plane 
 * with its neighbor in a path.
 * 
 * @param processedMoc Molecular cage being generated.
 * @param mocsInProgress List of cages in construction to be processed.
 * @param idStart Index of the first linkable atom in the path in construction.
 * @param newStarts Stack of atoms starting a path in construction.
 * @param newStartPos Position (point) of the next atom added in the path.
 * @param sub Substrate molecule.
 */
void addAromaticRing(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, Point_t newStartPos, Molecule_t* sub) {
	
	for (int i = 0; forEachNeighbor((atom(processedMoc,idStart)), i); i++) { // For every possible plans with starting atom's neighbors
		int idNext, idAtomCycle;
		Point_t nextPos;
		Point_t copyNewStartPos = newStartPos;
		Shell_t* moc = SHL_copy(processedMoc);
			
		Point_t neighborStartPos = coords(atom(moc, neighbor(atom(moc, idStart), i)));
		Point_t startPos = coords(atom(moc, idStart));
			
		// Add the first ring atom of already known position.
		int idnewStart = SHL_addAtom(moc, copyNewStartPos, -1);
		flag(atom(moc, idnewStart)) = CARBON;
		SHL_addEdge(moc, idStart, idnewStart);

		// Look for the normal to position the ring.
		Point_t normal = planNormal(copyNewStartPos, startPos, neighborStartPos);
		normal = rotation(normalization(vector(copyNewStartPos, startPos), 1),  90, normal); // Perpendicular
			
		// Position the other atoms of the cycle.
		//writeShellOutput("demos/substrates/ADENOS10.xyz",processedMoc,0);
		neighborStartPos = AX1E2(copyNewStartPos, coords(atom(moc, idStart)), normal, SIMPLE_CYCLE); // Neighbor
		//writeShellOutput("demos/substrates/ADENOS10.xyz",processedMoc,0);
		copyNewStartPos = AX2E1(copyNewStartPos, coords(atom(moc, idStart)), neighborStartPos, SIMPLE_CYCLE); 
		//writeShellOutput("demos/substrates/ADENOS10.xyz",processedMoc,0);
		if (isHindered(moc, sub, copyNewStartPos)) {
			SHL_delete(moc);
			return;
		}

		idAtomCycle = SHL_addAtom(moc, copyNewStartPos, -1);
		flag(atom(moc, idAtomCycle)) = CARBON;
		SHL_addEdge(moc, idnewStart, idAtomCycle);
			
		int idNeighbor = idAtomCycle;
		for (int i = 0; i < 4; i++) {
			neighborStartPos = coords(atom(moc, neighbor(atom(moc, idNeighbor), 0)));
			copyNewStartPos = AX1E2(copyNewStartPos, neighborStartPos, normal, SIMPLE_CYCLE);
			if (isHindered(moc, sub, copyNewStartPos)) {
				SHL_delete(moc);
				return;
			}
			idAtomCycle = SHL_addAtom(moc, copyNewStartPos, -1);
			flag(atom(moc, idAtomCycle)) = CARBON;
			SHL_addEdge(moc, idNeighbor, idAtomCycle);
				
			if (i == 1) {// Position the next starting atom to continue the path.
				nextPos = copyNewStartPos;
				idNext = idAtomCycle;
			}
			idNeighbor = idAtomCycle;
		}

		SHL_addEdge(moc, idnewStart, idAtomCycle);
			
		// Position atom after the cycle.
		neighborStartPos = coords(atom(moc, neighbor(atom(moc, idNext), 0)));
		Point_t v2 = coords(atom(moc, neighbor(atom(moc, idNext), 1)));
		copyNewStartPos = AX2E1(nextPos, neighborStartPos, v2, DIST_SIMPLE); 
		if (isHindered(moc, sub, copyNewStartPos) == 1) {
			SHL_delete(moc);
			return;
		}
		int idSuiv2 = SHL_addAtom(moc, copyNewStartPos, -1);
		flag(atom(moc, idSuiv2)) = CARBON;
		SHL_addEdge(moc, idNext, idSuiv2);
			
		LSTm_addElement(mocsInProgress, moc);
		LSTd_addElement(newStarts, idSuiv2);
	}
}

/**
 * @brief Add the oxygen atom of a carbonyl pattern (pattern 3).
 * 
 * @param processedMoc Molecular cage being generated.
 * @param idStart Index of the first linkable atom in the path in construction.
 * @param sub Substrate molecule.
 * @return (List_m*) List of cages with the added oxygen to be processed.
 */
List_m* addOxygenOfCarbonyl(Shell_t* processedMoc, int idStart, Molecule_t* sub) {
	
	List_m* mocsInProgress = LSTm_init();
	Point_t startPos = coords(atom(processedMoc, idStart));
	int idNeighbor1 = neighbor(atom(processedMoc, idStart), 0); // Neighbor
	Point_t neighbor1Pos = coords(atom(processedMoc, idNeighbor1));
	
	for (int i = 0; forEachNeighbor(atom(processedMoc,idNeighbor1), i); i++) { // For every possible plans with starting atom's neighbors.
		if (neighbor(atom(processedMoc, idNeighbor1), i) != idStart) {
			Shell_t* firstMoc = SHL_copy(processedMoc);
			Shell_t* secondMoc = SHL_copy(processedMoc);
			
			Point_t neighbor2Pos = coords(atom(firstMoc, neighbor(atom(firstMoc, idNeighbor1), i)));
						
			// Look for the normal to position the oxygen.
			Point_t normal = planNormal(startPos, neighbor1Pos, neighbor2Pos);
			
			// Oxygen
			// First position.
			Point_t oxygenPos = AX1E2(startPos, neighbor1Pos, normal, DIST_SIMPLE);
			int idOxygen;
			
			if (!isHindered(firstMoc, sub, oxygenPos)) {
				idOxygen = SHL_addAtom(firstMoc, oxygenPos, -1);
				flag(atom(firstMoc, idOxygen)) = OXYGEN;
				SHL_addEdge(firstMoc, idStart, idOxygen);
				
				LSTm_addElement(mocsInProgress, firstMoc);
			}
			else {
				SHL_delete(firstMoc);
			}
						
			// Second position
			oxygenPos = AX2E1(startPos, neighbor1Pos, oxygenPos, DIST_SIMPLE);
			
			if (!isHindered(secondMoc, sub, oxygenPos)) {
				int idOxygen = SHL_addAtom(secondMoc, oxygenPos, -1);
				flag(atom(secondMoc, idOxygen)) = OXYGEN;
				SHL_addEdge(secondMoc, idStart, idOxygen);
				
				LSTm_addElement(mocsInProgress, secondMoc);
			}
			else {
				SHL_delete(secondMoc);
			}
		}
	}
	return mocsInProgress;
}

/**
 * @brief Add a carbonyl pattern (c=0) with its neighbor.
 * 
 * @param processedMoc Molecular cage being generated.
 * @param mocsInProgress List of cages in construction to be processed.
 * @param idStart Index of the first linkable atom in the path in construction.
 * @param newStarts Stack of atoms starting a path in construction.
 * @param newStartPos Position (point) of the next atom added in the path.
 * @param sub Substrate molecule.
 */
void addCarbonyl(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, Point_t newStartPos, Molecule_t* sub) {
	
	for (int i = 0; forEachNeighbor(atom(processedMoc,idStart), i); i++) { // For every possible plans with starting atom's neighbors.
		Shell_t* firstMoc = SHL_copy(processedMoc);
		Shell_t* secondMoc = SHL_copy(processedMoc);
			
		Point_t neighbor1Pos = coords(atom(firstMoc, neighbor(atom(firstMoc, idStart), i)));
		Point_t startPos = coords(atom(firstMoc, idStart));
			
		// Carbon
		int idNewStartMoc1 = SHL_addAtom(firstMoc, newStartPos, -1);
		flag(atom(firstMoc, idNewStartMoc1)) = CARBON;
		SHL_addEdge(firstMoc, idStart, idNewStartMoc1);
		
		int idNewStartMoc2 = SHL_addAtom(secondMoc, newStartPos, -1);
		flag(atom(secondMoc, idNewStartMoc2)) = CARBON;
		SHL_addEdge(secondMoc, idStart, idNewStartMoc2);
			
		// Look for the normal to position the oxygen.
		Point_t normal = planNormal(newStartPos, startPos, neighbor1Pos);
		
		// Oxygen
		// First position.
		Point_t oxygenPos = AX1E2(newStartPos, startPos, normal, DIST_SIMPLE);
		int idOxygen;
			
		if (!isHindered(firstMoc, sub, oxygenPos)) {
			int idOxygen = SHL_addAtom(firstMoc, oxygenPos, -1);
			flag(atom(firstMoc, idOxygen)) = OXYGEN;
			SHL_addEdge(firstMoc, idNewStartMoc1, idOxygen);
				
			LSTm_addElement(mocsInProgress, firstMoc);
			LSTd_addElement(newStarts, idNewStartMoc1);
		}
		else {
			SHL_delete(firstMoc);
		}
			
		// Second position.
		oxygenPos = AX2E1(newStartPos, startPos, oxygenPos, DIST_SIMPLE);
			
		if (!isHindered(secondMoc, sub, oxygenPos)) {
			idOxygen = SHL_addAtom(secondMoc, oxygenPos, -1);
			flag(atom(secondMoc, idOxygen)) = OXYGEN;
			SHL_addEdge(secondMoc, idNewStartMoc2, idOxygen);
				
			LSTm_addElement(mocsInProgress, secondMoc);
			LSTd_addElement(newStarts, idNewStartMoc2);
		}
		else {
			SHL_delete(secondMoc);
		}		
	}
}

/**
 * @brief Add the projected atom(s) to the cage being generated.
 * 
 * @param processedMoc Molecular cage being generated.
 * @param mocsInProgress List of cages in construction to be processed.
 * @param idStart Index of the first linkable atom in the path in construction.
 * @param newStarts Stack of atoms starting a path in construction.
 * @param numPattern Pattern number (0, 1, 2) in the main loop. 
 * @param newStartPos Position (point) of the next atom added in the path.
 * @param sub Substrate molecule.
 */
void addProjection(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, int numPattern, Point_t newStartPos, Molecule_t* sub) {
	
	if (numPattern == 3) {
		addCarbonyl(processedMoc, mocsInProgress, idStart, newStarts, newStartPos, sub);
	}
	else if (numPattern == 4) {
		addAromaticRing(processedMoc, mocsInProgress, idStart, newStarts, newStartPos, sub);
	}
	else {
		Shell_t* moc = SHL_copy(processedMoc);
	
		int idnewStart = SHL_addAtom(moc, newStartPos, -1);
		flag(atom(moc, idnewStart)) = insertType(numPattern);
		SHL_addEdge(moc, idStart, idnewStart);
		
		LSTm_addElement(mocsInProgress, moc);
		LSTd_addElement(newStarts, idnewStart);
	}
}

/**************************************/
/********* Projection location ********/
/**************************************/

// Projection for an atom with one neighbor.
void projectionOCN_AX1E3(Shell_t* processedMoc, List_m* mocsInsProgress, int idStart, int idEnd, List_d* newStarts, int numPattern, Molecule_t* sub) {
	
	List_s* positions = LSTs_init();
	Point_t startPos = coords(atom(processedMoc, idStart));
	Point_t endPos = coords(atom(processedMoc, idEnd));
	int idFirstNeighborStart = neighbor(atom(processedMoc, idStart), 0);
	Point_t firstNeighborStartPos = coords(atom(processedMoc, idFirstNeighborStart));
	Point_t neighborOfFirstOnePos; // Neighbor of the the starting point's first neighbor.
	
	int idFirstNeighbor = neighbor(atom(processedMoc,idFirstNeighborStart),0);
	int idSecondNeighbor = neighbor(atom(processedMoc,idFirstNeighborStart),1);
	if (idFirstNeighbor == idStart)
		neighborOfFirstOnePos = coords(atom(processedMoc, idSecondNeighbor));
	else
		neighborOfFirstOnePos = coords(atom(processedMoc, idFirstNeighbor));

	Point_t normal = planNormal(startPos, firstNeighborStartPos, neighborOfFirstOnePos);
	
	Point_t newStartPos = AX1E3(startPos, firstNeighborStartPos, normal, DIST_SIMPLE);
	LSTs_addElement(positions, newStartPos);
	
	for (int i = 0; i < 11; i++) { // 360° rotation.
		normal = rotation(normalization(vector(startPos, firstNeighborStartPos), 1),  30, normal); // 30° rotation from normal.
		newStartPos = AX1E3(startPos, firstNeighborStartPos, normal, DIST_SIMPLE);
		
		if (!isHindered(processedMoc, sub, newStartPos)) {	
			LSTs_addElement(positions, newStartPos);
		}
	}
	
	for (int i = 0; i < 2 && positions->first; i++) { // Best placed position (min distance to the end)
		newStartPos = minDist(positions, endPos); 
		LSTs_removeElement(positions, newStartPos);
		addProjection(processedMoc, mocsInsProgress, idStart, newStarts, numPattern, newStartPos, sub);
	}
	LSTs_delete(positions);
}

// Projection for a nitrogen with two neighbors.
void projectionN_AX2E2(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, int numPattern, Molecule_t* sub) {
	
	int idFirstNeighbor = neighbor(atom(processedMoc, idStart), 0);
	int idSecondNeighbor = neighbor(atom(processedMoc, idStart), 1);
	Point_t newStartPos = AX2E2(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), DIST_SIMPLE);
	
	if (!isHindered(processedMoc, sub, newStartPos)) {
		addProjection(processedMoc, mocsInProgress, idStart, newStarts, numPattern, newStartPos, sub);
	}
}

// Projection for a carbon with two neighbors and one is an oxygen.
void projectionC_AX2E1(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, int numPattern, Molecule_t* sub) {
	
	int idFirstNeighbor = neighbor(atom(processedMoc, idStart), 0);
	int idSecondNeighbor = neighbor(atom(processedMoc, idStart), 1);
	Point_t newStartPos = AX2E1(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), DIST_SIMPLE);
	
	if (!isHindered(processedMoc, sub, newStartPos)) {
		addProjection(processedMoc, mocsInProgress, idStart, newStarts, numPattern, newStartPos, sub);
	}
}

// Projection for a carbon with two neighbors.
void projectionC_AX2E2(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, int numPattern, Molecule_t* sub) {
	
	int idFirstNeighbor = neighbor(atom(processedMoc, idStart), 0);
	int idSecondNeighbor = neighbor(atom(processedMoc, idStart), 1);
	Point_t newStartPos = AX2E2(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), DIST_SIMPLE);
	
	if (!isHindered(processedMoc, sub, newStartPos)) {
		addProjection(processedMoc, mocsInProgress, idStart, newStarts, numPattern, newStartPos, sub);
	}
	
	Point_t newStartPos2 = AX3E1(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), newStartPos, DIST_SIMPLE);
	
	if (!isHindered(processedMoc, sub, newStartPos2) == 0) {
		addProjection(processedMoc, mocsInProgress, idStart, newStarts, numPattern, newStartPos2, sub);
	}
}

// Projection for a carbone with 3 neighbors.
void projectionC_AX3E1(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, int numPattern,Molecule_t* sub) {
	
	int idFirstNeighbor = neighbor(atom(processedMoc, idStart), 0);
	int idSecondNeighbor = neighbor(atom(processedMoc, idStart), 1);
	int idThirdNeighbor = neighbor(atom(processedMoc, idStart), 2);
	Point_t newStartPos = AX3E1(coords(atom(processedMoc, idStart)), coords(atom(processedMoc, idFirstNeighbor)), coords(atom(processedMoc, idSecondNeighbor)), coords(atom(processedMoc, idThirdNeighbor)), DIST_SIMPLE);
	
	if (!isHindered(processedMoc, sub, newStartPos)) {
		addProjection(processedMoc, mocsInProgress, idStart, newStarts, numPattern, newStartPos, sub);
	}
}

/**************************************/
/*********** Generate path ************/
/**************************************/

/**
 * @brief Insert the pattern corresponding to number passed in argument.
 * 
 * @param processedMoc Molecular cage being generated.
 * @param mocsInProgress List of cages in construction to be processed.
 * @param idStart Index of the first linkable atom in the path in construction.
 * @param newStarts Stack of atoms starting a path in construction.
 * @param numPattern Pattern number (0, 1, 2) in the main loop. 
 * @param idEnd Index of the atom the path in construction is to be connected to.
 * @param sub Substrate molecule.
 */
void insertPattern(Shell_t* processedMoc, List_m* mocsInProgress, int idStart, List_d* newStarts, int numPattern, int idEnd, Molecule_t* sub){
	
	int numberOfNeighborsStart = LST_nbElements(neighborhood(atom(processedMoc, idStart)));
	if (numberOfNeighborsStart == 1) {
		//Projections
		//Diff rotations
		projectionOCN_AX1E3(processedMoc, mocsInProgress, idStart, idEnd, newStarts, numPattern, sub);
	}
	else if (flag(atom(processedMoc, idStart)) == NITROGEN && numberOfNeighborsStart == 2) {
		//Projection
		projectionN_AX2E2(processedMoc, mocsInProgress, idStart, newStarts, numPattern, sub);
	}
	else if (flag(atom(processedMoc, idStart)) == CARBON) {
		if (numberOfNeighborsStart == 2) {
			int idFirstNeighborStart = neighbor(atom(processedMoc, idStart), 0);
			int idSecondNeighborStart = neighbor(atom(processedMoc, idStart), 1);
			if (flag(atom(processedMoc, idFirstNeighborStart)) == OXYGEN || flag(atom(processedMoc, idSecondNeighborStart)) == OXYGEN) {
				// Projection
				projectionC_AX2E1(processedMoc, mocsInProgress, idStart, newStarts, numPattern, sub);
			}
			else {
				// 2 Projections
				projectionC_AX2E2(processedMoc, mocsInProgress, idStart, newStarts, numPattern, sub);
			}
		}
		else { 
			// 3 neighbors
			//Projection
			projectionC_AX3E1(processedMoc, mocsInProgress, idStart, newStarts, numPattern, sub);
		}
	}
}

/**
 * @brief Recursively generate paths between two grouping of bonding patterns.
 * 
 * @param m Grouping of the main structures (substrate and envelope).
 * @param mocsInProgress List of cages in construction to be processed.
 * @param processedMoc Molecular cage being generated.
 * @param idStart Index of the first linkable atom in the path in construction.
 * @param idEnd Index of the atom the path in construction is to be connected to.
 * @param nbCarbonyls Number of consecutive carbonyl patterns autorized.
 * @param nbAroRings Number of aromatic rings recquiried.
 * @param inputFile Name of the substrate's file.
 * @param sizeMax Maximale size (in atoms) of a path.
 * @param startingMocSize Size (in atoms) of the cage before adding the path.
 */
void generatePaths(Main_t* m, List_m* mocsInProgress, Shell_t* processedMoc, int idStart, int idEnd, int nbCarbonyls, int nbAroRings, char* inputFile, int sizeMax, int startingMocSize, int forceCycle) {
	/*************** Check distances bewteen atoms *****/
	Point_t B = coords(atom(processedMoc, idStart));
	for (int i = 0; i < size(processedMoc); i++) {
		if (i != idStart) {
			Point_t A = coords(atom(processedMoc, i));
			if (dist(A, B) < DIST_GAP_CAGE) {
				return;
			}
		}
	}
	for (int i = 0; i < size(substrat(m)); i++){
		Point_t A = coords(atom(substrat(m), i));
		if (dist(A, B) < DIST_GAP_SUBSTRATE) {
			return;
		}
	}
	/***************************************************/
	for (int i = 2; i < NB_PATTERNS; i++) {
		if(i == 3) i++; // TEMP exclude carbonyl pattern
		List_m* tempMocsInProg = LSTm_init();
		List_d* newStarts = LSTd_init();
		
		insertPattern(processedMoc, tempMocsInProg, idStart, newStarts, i, idEnd, substrat(m));
		
		while (tempMocsInProg->first) {
			// Count the number of consecutive carbonyls.
			if (i == 3) {
				nbCarbonyls++;
			}
			else {
				nbCarbonyls = 0;
			}
			
			// Count the number of aromatic rings.
			if (i == 4) {
				nbAroRings++;
			}
			
			if(sizeMax >= SHL_nbAtom(tempMocsInProg->first->moc) - startingMocSize) {
				if (dist( coords(atom(tempMocsInProg->first->moc, newStarts->first->idAtom)), coords(atom(processedMoc, idEnd)) ) < DIST_SIMPLE + DIST_ERROR) {
					float trA = dist( coords(atom(tempMocsInProg->first->moc, newStarts->first->idAtom)), coords(atom(tempMocsInProg->first->moc,neighbor(atom(tempMocsInProg->first->moc, newStarts->first->idAtom),0))));
					float trB = dist( coords(atom(tempMocsInProg->first->moc, newStarts->first->idAtom)), coords(atom(tempMocsInProg->first->moc, idEnd)));
					float trC = dist(coords(atom(tempMocsInProg->first->moc,neighbor(atom(tempMocsInProg->first->moc, newStarts->first->idAtom),0))), coords(atom(tempMocsInProg->first->moc, idEnd)));
					float trD = dist( coords(atom(tempMocsInProg->first->moc, idEnd)), coords(atom(tempMocsInProg->first->moc,neighbor(atom(tempMocsInProg->first->moc, idEnd),0))));
					float trE = dist(coords(atom(tempMocsInProg->first->moc,neighbor(atom(tempMocsInProg->first->moc, idEnd),0))), coords(atom(tempMocsInProg->first->moc, newStarts->first->idAtom)));
					float beforeLastAngle = radianToDegre(acosf(((trC * trC) - (trA * trA) - (trB * trB)) / (-2 * trA * trB)));
					float lastAngle = radianToDegre(acosf(((trE * trE) - (trD * trD) - (trB * trB)) / (-2 * trD * trB)));

					if(beforeLastAngle >= (END_ANGLE - ANGLE_ERROR) && beforeLastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle <= (END_ANGLE + ANGLE_ERROR) && lastAngle >= (END_ANGLE - ANGLE_ERROR)) {
						if(!forceCycle || (forceCycle && nbAroRings > 0)) { // Only if there is a cycle in the path and we force the presence of a cycle.
							flag(atom(tempMocsInProg->first->moc, idEnd)) = CARBON; // Change end atom (arrival) flag.
							SHL_addEdge(tempMocsInProg->first->moc, newStarts->first->idAtom, idEnd); //Add a link between last atom of the path and arrival.
							LSTm_addElement(mocsInProgress, SHL_copy(tempMocsInProg->first->moc));// Add to the list to be processed.
						}
					}
				}
				else if (nbCarbonyls < 5 && nbAroRings < 3) {
					generatePaths(m, mocsInProgress, tempMocsInProg->first->moc, newStarts->first->idAtom, idEnd, nbCarbonyls, nbAroRings, inputFile, sizeMax, startingMocSize, forceCycle);
				}
			}
			LSTm_removeFirst(tempMocsInProg);
			LSTd_removeFirst(newStarts);
		}
		LSTm_delete(tempMocsInProg);
		LSTd_delete(newStarts);
	}
}

/*************************************************/
/******* Start and end atoms of the paths ********/
/*************************************************/

/**
 * @brief Depth-first search.
 * 
 * @param s Cage without any added paths.
 * @param markedAtoms List of previously visited atoms.
 * @param index1 Starting atom index.
 * @param index2 Searched atom index.
 * @return int 
 */
int search(Shell_t* s, List_t* markedAtoms, int index1, int index2) {
	
	AtomShl_t* a = atom(s, index1);
	LST_addElement(markedAtoms, index1);
	
	if (neighborhoodSize(a) == 0) {
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
 * @brief Check if two atoms in the cage are connected.
 * Use a Depth-first search.
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
 * @brief Generate a list of every pairs of atoms that can be linked.
 * They must belong to different connected component.
 * 
 * @param s Cage without any added paths.
 * @return (List_p*) List of pair of atoms to be connected.
 */
List_p* chooseStartAndEndPairs(Shell_t* s) {
	
	List_p* startEndAtoms = LST2_init();
	
	for (int i = 0; i < size(s) - 1; i++) {
		if (flag(atom(s, i)) == LINKABLE_F) {
			for (int j = i + 1; j < size(s); j++) {
				if (flag(atom(s, j)) == LINKABLE_F) {
					if (!checkExistsPath(s, i, j)) {
						LST2_addElement(startEndAtoms, i, j);
					}
				}
			}
		}
	}
	return startEndAtoms;
}

/**
 * @brief Initialize the list of moc with the first pathless cage generated.
 * Delete the list of mocs of the main structure.
 * 
 * @param m Grouping of the main structures (substrate and envelope).
 * @return (List_m*) Empty list of mocs. 
 */
List_m* initMocsInProgress(Main_t* m){
	List_m* mocsInProgress = LSTm_init();
	
	for (int i = 0; i < mocSize(m); i++) {
		if (moc(m,i) != NULL) {
			if (i == 0) { // Only the first moc. TODO? change it
				LSTm_addElement(mocsInProgress, moc(m, i));
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
 * @brief Generate connected cages and write them to the results directory.
 * 
 * @param m Grouping of the main structures (substrate and envelope).
 * @param options Grouping of inputfile, alpha, sizeMax, maxResults.
 */
void generateWholeCages(Main_t* m, Options_t options) {
	
	printf("\n####### Start of paths generation #######\n");
	List_m* mocsInProgress = initMocsInProgress(m); // ! Take only the first moc.
	static int countResults = 0;

	// Remove the envelope's atoms.
	Shell_t* trimmedMoc;
	for (int j = 0; j < size(mocsInProgress->first->moc); j++) {
		if (flag(atom(mocsInProgress->first->moc,j)) == 0) {
			SHL_removeAtom(mocsInProgress->first->moc, j);
		}
	}
	trimmedMoc = SHL_copyCageAtoms(mocsInProgress->first->moc);
	LSTm_removeFirst(mocsInProgress);
	LSTm_addElement(mocsInProgress, trimmedMoc);

	int pathelessMocSize = SHL_nbAtom(mocsInProgress->first->moc); // Allows to recover the size before the addition of the paths, only if we keep one moc line (TODO modify otherwise).
	while (mocsInProgress->first) { // As long as there is a moc to process.	

		int startingMocSize = SHL_nbAtom(mocsInProgress->first->moc);
		List_p* startEndAtoms = chooseStartAndEndPairs(mocsInProgress->first->moc);
		
		if (!startEndAtoms->first) { // If there is only one grouping of patterns left (connected cage).
			if (countResults++ < options.maxResults) {
				writeShellOutput(options.input, mocsInProgress->first->moc, pathelessMocSize);
				LSTm_removeFirst(mocsInProgress);
			}
			else {
				LST2_delete(startEndAtoms);
				free(mocsInProgress);
				return;
			}
		}
		else { // If there are at least 2 groupings of patterns. 	
			Shell_t* processedMoc = mocsInProgress->first->moc;
			mocsInProgress->first->moc = NULL;
			LSTm_removeFirst(mocsInProgress);
			
			while (startEndAtoms->first) { // For all pairs of atoms to connect.
			
				int idStart = startEndAtoms->first->start;
				int idEnd = startEndAtoms->first->end;
				int forceCycle = 0;
				float startEndDist = dist(coords(atom(processedMoc,idStart)),coords(atom(processedMoc,idEnd)));
				
				if (startEndDist <= DIST_SIMPLE_PATTERN * options.sizeMax + DIST_SIMPLE + DIST_ERROR) {
					if (startEndDist <= DIST_SIMPLE_PATTERN * (options.sizeMax - NUMBER_ATOM_CYCLE_PATTERN) + DIST_CYCLE_PATTERN + DIST_SIMPLE + DIST_ERROR
					&& startEndDist > DIST_CYCLE_PATTERN + (1 * DIST_SIMPLE_PATTERN) + DIST_SIMPLE + DIST_ERROR) {
						forceCycle = 1;
					}
					Shell_t* appendedMoc = SHL_copy(processedMoc); // Create a new moc in the list to process.
					
		//#pragma omp parallel for
					for (int i = 2; i < 3/*4 with carbonyl*/; i++) { // Assignment of all types to the starting atom (atom on the edges).

						flag(atom(appendedMoc, idStart)) = insertType(i);
						if (i == 3) {
							if (LST_nbElements(neighborhood(atom(appendedMoc, idStart))) == 1) // Carbonyl possible only if the starting atom has only one neighbor.
							{
								List_m* mocsWithCarbonyl = addOxygenOfCarbonyl(appendedMoc, idStart,substrat(m));

								while (mocsWithCarbonyl->first) { // Process all mocs generated by this addition.
									generatePaths(m, mocsInProgress, mocsWithCarbonyl->first->moc, idStart, idEnd, 0, 0, options.input, options.sizeMax, startingMocSize, forceCycle);
									LSTm_removeFirst(mocsWithCarbonyl);
								}
								LSTm_delete(mocsWithCarbonyl);
							}
						}
						else {	
							generatePaths(m, mocsInProgress, appendedMoc, idStart, idEnd, 0, 0, options.input, options.sizeMax, startingMocSize, forceCycle);
						}
					}
					SHL_delete(appendedMoc);
				}
				LST2_removeFirst(startEndAtoms);
			}
			SHL_delete(processedMoc);
		}
		LST2_delete(startEndAtoms);
	}
	free(mocsInProgress);
}
