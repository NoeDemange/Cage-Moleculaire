#ifndef __CONSTANT_H
#define __CONSTANT_H

// Main
#define OPTSTR "i:a:s:r:h"
#define USAGE_FMT  "usage : [-i inputfile] [-a alpha (default : %1.f)] [-s sizemax (default : %d)] [-r maxresults (default : %d)] [-h]\n"
#define DEFLT_ALPHA 3.
#define DEFLT_SIZEMAX 5
#define DEFLT_MAX_RESULTS 10

#define PATHNAME "alphashape.R"

// Structure
#define REALLOCSIZE 4 // TODO could it be decreased? 

// Initialization of the substrate
#define EDGE_ERROR 20 // Acceptable error for the computation of the edges between atoms

// Distance
#define DIST_HYDRO 1.8 // Hydrogen bond size.
#define DIST_SIMPLE 1.5 // Simple covalent bond size.
#define DIST_ERROR 0.5 // Acceptable error for the covalent bond size. 
#define MINDIS 0.75 // Minimal distance bewteen two atoms (otherwise they are merged).
#define DIST_ATOM_H (DIST_SIMPLE/2)+(MINDIS/2) // Distance between one atom and an atom of hydrogen.
#define DIST_GAP_CAGE ((DIST_SIMPLE/2)+(MINDIS/2)-0.0001) // Distance between two cage's atoms (defined by trial). distance H-C
#define DIST_GAP_SUBSTRATE 1.8 // Distance between an atom of the cage and an atom of the substrate (at least a hydrogen bond size by trial).
#define DIST_SIMPLE_PATTERN 1.22 // AC/2, with ABC a triangle where each of its vertex is an atom of the path involved in a simple pattern, AB = BC = 1.5 and angle ABC = 109°C.
/*  B
 \ / \ /
  A   C
*/
#define DIST_CYCLE_PATTERN 5.8 // Distance bewteen two neighbors of a cycle ([AB] = 1.5 + 0.7 + 1.4 + 0.7 + 1.5).
/*   ___
A___/   \___B
    \___/
*/
#define NB_ATOMS_IN_CYCLE 7 // Number of atoms in an aromatic ring pattern.

// Distance in generateCycle (TODO document why these values were chosen)
#define SIMPLE_CYCLE 1.4 // Simple covalent bond size between an atom involved in a cycle and a neighboring atom outside of the cycle.
#define MINDIS_CYCLE 0.7 // Minimal distance bewteen two atoms when one of them belong to a cycle (otherwise they are merged).
#define MAXDIS_CYCLE 1.7 // Maximal distance between two atoms of a cycle (otherwise they can't be both in the same cycle).

// Angle
#define END_ANGLE 109.47
#define ANGLE_ERROR 10

// Path choice
#define NUMBER_POSITION_AX1E3 3 //Number of position keep to make path
#define ROTATION_ANGLE_AX1E3 15 //rotation angle to generate position (degree)

/********* not to be modified (the incremental order must be preserved) */
// Flags atoms in the envelope and cage
#define NOT_DEF_F -1 // Atom not used
#define SHELL_F 0 // Atom of the shell
#define LINKABLE_F 1 // Atom at the edge of a pattern that can still make another bond (unless it's a hydrogen).
#define CYCLE_F 2 // Atom in an aromatic ring that can't make another bond.
#define HYDRO_PATTERN_F 3 // Atom involved in a hydrogen pattern that can't make another bond.
/***********************************************************************/

// Path (in the cage)
#define NB_PATTERNS 2 // Number of patterns available to make a path
#define CYCLE_PATTERN 1 // Number for cycle pattern

// Flags atoms in paths, ! must be different of flags in the envelope
#define CARBON_F 6
#define NITROGEN_F 5
#define OXYGEN_F 4
#define HYDROGEN_F 7

// Voxelization
#define GRID_SIZE 201 //le nombre de voxel qu'on veut +1 car on place les points
#define START_GRID -30.0
#define LENGTH_GRID (-START_GRID*2)/GRID_SIZE

#endif