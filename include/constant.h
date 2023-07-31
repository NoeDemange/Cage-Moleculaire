#ifndef __CONSTANT_H
#define __CONSTANT_H

/** @file constant.h
 *  @brief Constants and options for the program.
 */

// Main : Options for execution

/** @def OPTSTR
 *  @brief Command-line options string for getopt.
 */
#define OPTSTR "i:a:s:r:h"

/** @def USAGE_FMT
 *  @brief Usage format string for displaying program options.
 */
#define USAGE_FMT  "usage : [-i inputfile] [-a alpha (default : %1.f)] [-s sizemax (default : %d)] [-r maxresults (default : %d)] [-h]\n"

/** @def DEFLT_ALPHA
 *  @brief Default alpha value for alphashape3d.
 */
#define DEFLT_ALPHA 3. 

/** @def DEFLT_SIZEMAX
 *  @brief Default maximum size of a path created in number of patterns.
 */
#define DEFLT_SIZEMAX 5 

/** @def DEFLT_MAX_RESULTS
 *  @brief Default number of results generated.
 */
#define DEFLT_MAX_RESULTS 10 

/** @def PATHNAME
 *  @brief File path for the R script.
 */
#define PATHNAME "alphashape.R"

// Structure

/** @def REALLOCSIZE
 *  @brief Size of reallocation for dynamic memory allocation.
 */
#define REALLOCSIZE 4 // TODO could it be decreased? 

// Initialization of the substrate

/** @def EDGE_ERROR
 *  @brief Acceptable error for the computation of the edges between atoms.
 */
#define EDGE_ERROR 20

// Distance

/** @def DIST_HYDRO
 *  @brief Hydrogen bond size.
 */
#define DIST_HYDRO 1.8

/** @def DIST_SIMPLE
 *  @brief Simple covalent bond size.
 */
#define DIST_SIMPLE 1.5

/** @def DIST_ERROR
 *  @brief Acceptable error for the covalent bond size.
 */
#define DIST_ERROR 0.5

/** @def MINDIS
 *  @brief Minimal distance between two atoms (otherwise they are merged).
 */
#define MINDIS 0.75 

/** @def DIST_ATOM_H
 *  @brief Distance between one atom and an atom of hydrogen.
 */
#define DIST_ATOM_H (DIST_SIMPLE/2)+(MINDIS/2)

/** @def DIST_GAP_CAGE
 *  @brief Distance between two cage's atoms (defined by trial). distance H-C.
 */
#define DIST_GAP_CAGE ((DIST_SIMPLE/2)+(MINDIS/2)-0.0001)

/** @def DIST_GAP_SUBSTRATE
 *  @brief Distance between an atom of the cage and an atom of the substrate (at least a hydrogen bond size by trial).
 */
#define DIST_GAP_SUBSTRATE 1.8

/** @def DIST_SIMPLE_PATTERN
 *  @brief AC/2, with ABC a triangle where each of its vertex is an atom of the path involved in a simple pattern,
 *         AB = BC = 1.5 and angle ABC = 109°C.
 *            B
 *         \ / \ /
 *          A   C
 */
#define DIST_SIMPLE_PATTERN 1.22


/** @def DIST_CYCLE_PATTERN
 *  @brief Distance between two neighbors of a cycle ([AB] = 1.5 + 0.7 + 1.4 + 0.7 + 1.5).
 *        ___
 *   A___/   \___B
 *       \___/
 */
#define DIST_CYCLE_PATTERN 5.8


/** @def NB_ATOMS_IN_CYCLE
 *  @brief Number of atoms in an aromatic ring pattern.
 */
#define NB_ATOMS_IN_CYCLE 7

// Distance in generateCycle (TODO document why these values were chosen)

/** @def SIMPLE_CYCLE
 *  @brief Simple covalent bond size between an atom involved in a cycle and a neighboring atom outside of the cycle.
 */
#define SIMPLE_CYCLE 1.4

/** @def MINDIS_CYCLE
 *  @brief Minimal distance between two atoms when one of them belongs to a cycle (otherwise they are merged).
 */
#define MINDIS_CYCLE 0.7

/** @def MAXDIS_CYCLE
 *  @brief Maximal distance between two atoms of a cycle (otherwise they can't be both in the same cycle).
 */
#define MAXDIS_CYCLE 1.7

// Angle

/** @def END_ANGLE
 *  @brief The end angle for path choice.
 */
#define END_ANGLE 109.47

/** @def ANGLE_ERROR
 *  @brief Acceptable error for angle calculations.
 */
#define ANGLE_ERROR 10

// Path choice

/** @def NUMBER_POSITION_AX1E3
 *  @brief Number of positions to keep to make a path.
 */
#define NUMBER_POSITION_AX1E3 2

/** @def ROTATION_ANGLE_AX1E3
 *  @brief Rotation angle to generate position (in degrees).
 */
#define ROTATION_ANGLE_AX1E3 30

/********* not to be modified (the incremental order must be preserved) */

// Flags atoms in the envelope and cage

/** @def NOT_DEF_F
 *  @brief Atom not used.
 */
#define NOT_DEF_F -1

/** @def SHELL_F
 *  @brief Atom of the shell.
 */
#define SHELL_F 0

/** @def LINKABLE_F
 *  @brief Atom at the edge of a pattern that can still make another bond (unless it's a hydrogen).
 */
#define LINKABLE_F 1

/** @def CYCLE_F
 *  @brief Atom in an aromatic ring that can't make another bond.
 */
#define CYCLE_F 2

/** @def HYDRO_PATTERN_F
 *  @brief Atom involved in a hydrogen pattern that can't make another bond.
 */
#define HYDRO_PATTERN_F 3

/***********************************************************************/

// Path (in the cage)

/** @def NB_PATTERNS
 *  @brief Number of patterns available to make a path.
 */
#define NB_PATTERNS 2

/** @def CYCLE_PATTERN
 *  @brief Number for the cycle pattern.
 */
#define CYCLE_PATTERN 1

// Flags atoms in paths, ! must be different of flags in the envelope

/** @def CARBON_F
 *  @brief Flag for carbon atom in paths.
 */
#define CARBON_F 6

/** @def NITROGEN_F
 *  @brief Flag for nitrogen atom in paths.
 */
#define NITROGEN_F 5

/** @def OXYGEN_F
 *  @brief Flag for oxygen atom in paths.
 */
#define OXYGEN_F 4

/** @def HYDROGEN_F
 *  @brief Flag for hydrogen atom in paths.
 */
#define HYDROGEN_F 7

// Voxelization

/** @def GRID_SIZE
 *  @brief Number of voxels (+1 because we are placing points).
 */
#define GRID_SIZE 51

/** @def START_GRID
 *  @brief Origins of the grids.
 */
#define START_GRID -30.0

/** @def LENGTH_GRID
 *  @brief Real length between two points.
 */
#define LENGTH_GRID (-START_GRID*2)/GRID_SIZE

#endif