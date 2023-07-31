#ifndef __STRUCTURE_H
#define __STRUCTURE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constant.h"

/**
 * @file structure.h
 * @brief Structure Header File
 *
 * This file contains data structures and macros used throughout the application.
 */

// MACROS FOR ACCESSING STRUCTURE ELEMENTS

/**
 * @brief Access the address of an atom in the shell structure.
 *
 * This macro returns the address of the atom at the specified index in the shell structure.
 *
 * @param o The shell structure.
 * @param i The index of the atom.
 * @return The address of the atom.
 */
#define atom(o,i) ((o)->atoms+(i))

/**
 * @brief Access the neighborhood of an atom.
 *
 * This macro returns the address of the neighborhood list of an atom in the shell structure.
 *
 * @param a The atom in the shell structure.
 * @return The address of the neighborhood list of the atom.
 */
#define neighborhood(a) (a)->neighborhood

/**
 * @brief Access the size of the shell structure.
 *
 * This macro returns the size of the shell structure.
 *
 * @param o The shell structure.
 * @return The size of the shell structure.
 */
#define size(o) (o)->size

/**
 * @brief Check if an atom is part of the cycle.
 *
 * This macro checks if the atom at the specified index is part of the cycle.
 *
 * @param o The shell structure.
 * @param i The index of the atom.
 * @return 1 if the atom is part of the cycle, 0 otherwise.
 */
#define cycle(o, i) LST_check((o)->cycle, (i))

/**
 * @brief Check if an atom is a vertex in the Graph structure.
 *
 * This macro checks if the atom at the specified index is a vertex in the Graph structure.
 *
 * @param o The Graph structure.
 * @param i The index of the atom.
 * @return 1 if the atom is a vertex, 0 otherwise.
 */
#define checkVertex(o, i) GPH_checkVertex((o)->bond, (i))

/**
 * @brief Check if there is a bond between two atoms in the Graph structure.
 *
 * This macro checks if there is a bond between the atoms at the specified indices in the Graph structure.
 *
 * @param o The Graph structure.
 * @param i The index of the first atom.
 * @param j The index of the second atom.
 * @return 1 if there is a bond between the atoms, 0 otherwise.
 */
#define checkBond(o, i, j)	GPH_checkBond((o)->bond, (i), (j))

/**
 * @brief Access the bond graph in the molecule structure.
 *
 * This macro returns the address of the bond graph in the molecule structure.
 *
 * @param o The molecule structure.
 * @return The address of the bond graph.
 */
#define bond(o)	(o)->bond

/**
 * @def elts(l, i)
 * @brief Access the i-th element of the list l.
 */
#define elts(l, i) (l)->elts[(i)]

/**
 * @def forEachElement(l, i)
 * @brief Macro to iterate through a list until it finds an unused value (-1).
 * Only used to scroll through lists where unused values are at the end.
 */
#define forEachElement(l,i) (i) < size((l)) && elts((l),(i)) != -1

/**
 * @def coords(a)
 * @brief Access the coordinates of an atom (a).
 */
#define coords(a) (a)->coords

/**
 * @def atomX(a)
 * @brief Access the coordinate x of an atom (a).
 */
#define atomX(a) (a)->coords.x

/**
 * @def atom(a)
 * @brief Access the coordinate y of an atom (a).
 */
#define atomY(a) (a)->coords.y

/**
 * @def atomZ(a)
 * @brief Access the coordinate z of an atom (a).
 */
#define atomZ(a) (a)->coords.z

/**
 * @def neighbor(a, i)
 * @brief Access the i-th neighbor of an atom (a).
 */
#define neighbor(a,i) (a)->neighborhood->elts[(i)]

/**
 * @def neighborhoodSize(a)
 * @brief Get the size of the neighborhood of an atom (a).
 */
#define neighborhoodSize(a) (a)->neighborhood->size

/**
 * @def forEachNeighbor(a, i)
 * @brief Macro to iterate through the neighbors of an atom until it finds an unused value (-1).
 * The unused values of neighbors must be at the end.
 */
#define forEachNeighbor(a,i) (i) < neighborhoodSize((a)) && neighbor((a),(i)) != -1

/**
 * @def symbol(a)
 * @brief Access the symbol of an atom (a) in a molecule.
 */
#define symbol(a) (a)->info.symbol

/**
 * @def radius(a)
 * @brief Access the radius of an atom (a) in a molecule.
 */
#define radius(a) (a)->info.radius

/**
 * @def ligands(a)
 * @brief Access the number of ligands of an atom (a) in a molecule.
 */
#define ligands(a) (a)->info.ligands

/**
 * @def lonePairs(a)
 * @brief Access the number of lone pairs of an atom (a) in a molecule.
 */
#define lonePairs(a) (a)->info.lonePairs

/**
 * @def steric(a)
 * @brief Calculate the steric hindrance of an atom (a) in a molecule (sum of ligands and lone pairs).
 */
#define steric(a) (a)->info.ligands + (a)->info.lonePairs


/**
 * @def flag(a)
 * @brief Access the flag of a shell atom (a).
 */
#define flag(a) (a)->flag

/**
 * @def parentAtom(a)
 * @brief Access the parent atom index of a shell atom (a).
 */
#define parentAtom(a) (a)->parentAtom

/**
 * @def vertex(g, i)
 * @brief Access the i-th vertex of a graph (g).
 */
#define vertex(g, i) ((g)->vertices+(i))

/**
 * @def id(v)
 * @brief Access the ID of a graph vertex (v).
 */
#define id(v)	(v)->id

/**
 * @def nbNeighbors(v)
 * @brief Access the number of neighbors of a graph vertex (v).
 */
#define nbNeighbors(v) (v)->nbNeighbors

//macro main
/**
 * @def substrat(m)
 * @brief Access the substrate of the main structure (m).
 */
#define substrat(m) (m)->substrat

/**
 * @def envelope(m)
 * @brief Access the envelope of the main structure (m).
 */
#define envelope(m) (m)->envelope

/**
 * @def envarom(m)
 * @brief Access the aromatic envelope of the main structure (m).
 */
#define envarom(m)	(m)->envarom

/**
 * @def moc(m, i)
 * @brief Access the i-th moc (Molecular Operating Cage) of the main structure (m).
 */
#define moc(m,i) (m)->mocs[i]

/**
 * @def mocSize(m)
 * @brief Access the size of the moc array of the main structure (m).
 */
#define mocSize(m) (m)->mocSize

/**
 * @def getName(var)
 * @brief Get the name of a variable as a string.
 */
#define getName(var)  #var

/**************************************/
/* POINT ******************************/
/**************************************/

/**
 * @struct Point3D
 * @brief Represents a 3D point with integer coordinates.
 */
typedef struct {
    int x; /**< The x-coordinate of the point. */
    int y; /**< The y-coordinate of the point. */
    int z; /**< The z-coordinate of the point. */
} Point3D;

/**
 * @struct Point_t
 * @brief Represents a 3D point with floating-point coordinates.
 */
typedef struct {
	float x; /**< The x-coordinate of the point. */
	float y; /**< The y-coordinate of the point. */
	float z; /**< The z-coordinate of the point. */
} Point_t;

/**************************************/
/* LISTE ******************************/
/**************************************/
/**
 * @struct List_t
 * @brief Represents a list of integers.
 */
typedef struct {
	int* elts; /**< An array of integers representing the elements in the list. */
	unsigned size; /**< The size of the list. */	
} List_t;

/**
 * @struct Element
 * @brief Represents an element in the list of vertices to be connected.
 */
typedef struct Element Element;
struct Element {
	int start; /**< The starting point of the edge. */
	int end; /**< The ending point of the edge. */
	float distance; /**< The distance between the starting and ending points. */
	Element *next; /**< Pointer to the next element in the list. */
};

/**
 * @struct Elem_s
 * @brief Represents an element in the list of intermediate vertices.
 */
typedef struct Elem_s Elem_s;
struct Elem_s {
	Point_t position; /**< The position of the intermediate vertex. */
	float distance; /**< The distance of the intermediate vertex from the starting point. */
	Elem_s *next; /**< Pointer to the next element in the list. */
};

/**
 * @struct List_s
 * @brief Represents a list of intermediate vertices.
 */
typedef struct {
	Elem_s *first; /**< Pointer to the first element in the list. */
} List_s;

/**
 * @struct Elem_d
 * @brief Represents an element in the list of integers.
 */
typedef struct Elem_d Elem_d;
struct Elem_d {
	int idAtom; /**< The ID of the atom in the list. */
	Elem_d *next; /**< Pointer to the next element in the list. */	

};

/**
 * @struct List_d
 * @brief Represents a list of integers.
 */
typedef struct {
	Elem_d *first; /**< Pointer to the first element in the list. */	
} List_d;

/**************************************/
/* GRAPHE *****************************/
/**************************************/
/**
 * @struct Vertex_t
 * @brief Represents a vertex in the graph.
 */
typedef struct {
	unsigned id; /**< The ID of the vertex. */
	// Voisinage
	List_t* neighborhood; /**< The list of neighbors of the vertex. */
	unsigned nbNeighbors; /**< The number of neighbors of the vertex. */
} Vertex_t;

/**
 * @struct Graph_t
 * @brief Represents a graph.
 */
typedef struct {
	Vertex_t* vertices; /**< The array of vertices in the graph. */
	unsigned size; /**< The size of the graph (number of vertices). */
} Graph_t;

/**
 * @struct Node
 * @brief Represents a node in the pathfinding algorithm.
 */
typedef struct {
    Point3D point; /**< The 3D point representing the node. */
    float g; /**< The g-value of the node. */
    float h; /**< The h-value of the node. */
    float f; /**< The f-value of the node. */
} Node;

/**************************************/
/* VOXEL ***************************/
/**************************************/

/**
 * @typedef VOXELGRID
 * @brief Represents a 3D voxel grid.
 */
typedef int*** VOXELGRID;


/**************************************/
/* MOLECULE ***************************/
/**************************************/
/**
 * @struct ChemicalInfo_t
 * @brief Represents information about a chemical element.
 */
typedef struct {
	char symbol[2]; /**< The chemical symbol of the element. */
	int radius; /**< The covalent radius of the element. */
	int ligands; /**< The number of ligands of the element. */
	int lonePairs; /**< The number of lone pairs of the element. */
} ChemicalInfo_t;

/**
 * @struct Atom_t
 * @brief Represents an atom in a molecule.
 */
typedef struct {
	ChemicalInfo_t info; /**< Information about the chemical element of the atom. */
	Point_t coords; /**< The 3D coordinates of the atom. */
	// Voisinage
	List_t* neighborhood; /**< The list of neighboring atoms of the atom. */
} Atom_t;

/**
 * @struct Molecule_t
 * @brief Represents a molecule.
 *
 * This structure holds information about a molecule, including its atoms, cycles, bonds, and size.
 */
typedef struct {
	Atom_t* atoms; /**< Array of atoms in the molecule. */
	List_t* cycle; /**< List of vertices belonging to a cycle. */
	Graph_t* bond; /**< Graph representing bonds between atoms. */
	unsigned size; /**< Size of the molecule (number of atoms). */
} Molecule_t;

/**************************************/
/* SHELL ******************************/
/**************************************/
/**
 * @struct AtomShl_t
 * @brief Represents an atom in a shell.
 *
 * This structure represents an atom within a shell, including its flag, coordinates,
 * parent atom index, and neighborhood.
 */
typedef struct {
	int flag; /**< Flag representing the state of the atom. */
	Point_t	coords; /**< Coordinates of the atom. */
	// Lien
	unsigned parentAtom; /**< Index of the parent atom in the shell. */
	// Voisinage
	List_t* neighborhood; /**< List of neighbors of the atom in the shell. */
} AtomShl_t;

/**
 * @struct Shell_t
 * @brief Represents a shell.
 *
 * This structure holds information about a shell, including its atoms, cycles, bonds, and size.
 */
typedef struct {
	AtomShl_t* atoms; /**< Array of atoms in the shell. */
	List_t* cycle; /**< List of vertices belonging to a cycle in the shell. */
	Graph_t* bond; /**< Graph representing bonds between atoms in the shell. */
	unsigned size; /**< Size of the shell (number of atoms). */
} Shell_t;

/**************************************/
/* MAIN *******************************/
/**************************************/
/**
 * @struct Main_t
 * @brief Represents the main structures of the application.
 *
 * This structure contains the main structures used in the application, including the substrate,
 * envelope, aromatic shell, molecular cages, and their sizes.
 */
typedef struct {
	Molecule_t* substrat; /**< Pointer to the substrate molecule. */
	Shell_t* envelope; /**< Pointer to the envelope shell. */
	Shell_t* envarom; /**< Pointer to the aromatic shell. */
	Shell_t** mocs; /**< Array of pointers to molecular cages. */
	unsigned mocSize; /**< Size of the array of molecular cages. */
} Main_t;

/**************************************/
/* ASHAPE3D ***************************/
/**************************************/
/**
 * @struct Ashape_t
 * @brief Represents the data for Ashape3D.
 *
 * This structure holds data for Ashape3D, including the number of triangles, edges, vertices,
 * and alpha values, along with arrays to store their values.
 */
typedef struct {
	int nb_triang; /**< Number of triangles. */
	int nb_edge; /**< Number of edges. */
	int nb_vertex; /**< Number of vertices. */
	int nb_x; /**< Number of x-values. */
	int nb_alpha; /**< Number of alpha values. */
	double* triang; /**< Array to store triangle values. */
	double* edge; /**< Array to store edge values. */
	double* vertex; /**< Array to store vertex values. */
	double* x; /**< Array to store x-values. */
	double* alpha; /**< Array to store alpha values. */
} Ashape_t;

/**************************************/
/* LISTE MOC **************************/
/**************************************/

/**
 * @struct Elem
 * @brief Represents an element in the list of molecular cages.
 *
 * This structure represents an element in the list of molecular cages to be processed,
 * containing information about the molecular cage, the number of patterns, and the number of cycles.
 */
typedef struct Elem Elem;
struct Elem {
	Shell_t* moc; /**< Pointer to the molecular cage. */
	int nbPatterns; /**< Number of patterns in the molecular cage. */
	int nbCycles; /**< Number of cycles in the molecular cage. */
	Elem* next; /**< Pointer to the next element in the list. */
};

/**
 * @struct List_m
 * @brief Represents a list of molecular cages.
 *
 * This structure represents a list of molecular cages to be processed.
 */
typedef struct {
	Elem* first; /**< Pointer to the first element in the list. */
} List_m;

/**************************************/
/* NODE HEAP **************************/
/**************************************/
/**
 * @struct VMap
 * @brief Represents a voxel map for pathfinding.
 *
 * This structure represents a voxel map for pathfinding, containing the distance and index heap information
 * for pathfinding algorithms.
 */
typedef struct {
	float dist; /**< Distance for pathfinding. */
	int indexHeap; /**< Index heap information for pathfinding. (-1 if already seen) */
} VMap;

/**
 * @struct NodeHeap
 * @brief Represents a node heap for pathfinding.
 *
 * This structure represents a node heap used in pathfinding algorithms, containing the size and array of nodes.
 */
typedef struct {
	int size; /**< Size of the node heap. */
	Node* node; /**< Array of nodes in the heap. */
} NodeHeap;

//NodeHeap
/**
 * @brief Initializes and allocates memory for a new NodeHeap.
 *
 * This function initializes a new NodeHeap structure, sets its size to 0, and allocates memory
 * for the node array based on the GRID_SIZE constant. The function returns the initialized NodeHeap.
 *
 * @return The initialized NodeHeap.
 */
NodeHeap NH_initAlloc();

/**
 * @brief Frees the memory allocated for a NodeHeap.
 *
 * This function frees the memory allocated for the node array in the NodeHeap.
 *
 * @param nodeHeap The NodeHeap to be freed.
 */
void NH_free(NodeHeap);

/**
 * @brief Inserts a new node into the NodeHeap.
 *
 * This function inserts a new node into the NodeHeap and adjusts the heap to maintain its binary property.
 *
 * @param nodeHeap Pointer to the NodeHeap structure.
 * @param p The new node to be inserted.
 * @param vMap Pointer to the 3D array of VMap structures.
 */
void NH_insert(NodeHeap*, Node, VMap***);

/**
 * @brief Extracts the node with the minimum f value from the NodeHeap.
 *
 * This function extracts the node with the minimum f value (the root) from the NodeHeap,
 * adjusts the heap, and returns the extracted node.
 *
 * @param nodeHeap Pointer to the NodeHeap structure.
 * @param vMap Pointer to the 3D array of VMap structures.
 * @return The node with the minimum f value.
 */
Node NH_extractMin(NodeHeap*, VMap***);

/**
 * @brief increases the priority of a node in the NodeHeap.
 *
 * This function increases the priority (decreases f value) of a node in the NodeHeap
 * and adjusts the heap to maintain its binary property.
 *
 * @param nodeHeap Pointer to the NodeHeap structure.
 * @param vMap Pointer to the 3D array of VMap structures.
 * @param point The 3D coordinates of the node to be updated.
 * @param distG The new g (distance from the start node) value for the node.
 */
void NH_increase_priority(NodeHeap* nodeHeap, VMap*** vMap, Point3D point, float distG);

/**
 * @brief Prints the priorities of nodes in the NodeHeap.
 *
 * This function prints the priorities (f values) of all nodes in the NodeHeap.
 * It is used for debugging and testing purposes.
 *
 * @param nodeHeap The NodeHeap to be printed.
 */
void NH_print(NodeHeap);

/**
 * @brief Allocates memory for the VMap data structure.
 *
 * This function allocates memory for the VMap 3D array, which is used for storing
 * node indices in the NodeHeap. The function returns the pointer to the allocated VMap.
 *
 * @return Pointer to the allocated VMap 3D array.
 */
VMap*** VMap_alloc();

/**
 * @brief Frees the memory allocated for the VMap data structure.
 *
 * This function frees the memory allocated for the VMap 3D array.
 *
 * @param vMap Pointer to the VMap 3D array to be freed.
 */
void VMap_free(VMap***);

//Point
/**
 * @brief Initializes a new Point_t with equal scalar values.
 *
 * This function initializes a new Point_t with the same scalar value for all coordinates (x, y, z).
 *
 * @param scal The scalar value to set for all coordinates of the Point_t.
 * @return The initialized Point_t.
 */
Point_t PT_init(float);

/**
 * @brief Adds two Point_t together and returns the result.
 *
 * This function adds two Point_t (A and B) together and returns the resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the addition operation as a new Point_t.
 */
Point_t PT_add(Point_t, Point_t);

/**
 * @brief Subtracts one Point_t from another and returns the result.
 *
 * This function subtracts Point_t B from Point_t A and returns the resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the subtraction operation as a new Point_t.
 */
Point_t PT_sub(Point_t, Point_t);

/**
 * @brief Multiplies a Point_t by a scalar value and returns the result.
 *
 * This function multiplies each coordinate of the Point_t A by the scalar value "scal" and returns the resulting Point_t.
 *
 * @param A The Point_t operand.
 * @param scal The scalar value to multiply each coordinate of the Point_t A.
 * @return The result of the multiplication operation as a new Point_t.
 */
Point_t PT_mul(Point_t, float);

/**
 * @brief Divides a Point_t by a scalar value and returns the result.
 *
 * This function divides each coordinate of the Point_t A by the scalar value "scal" and returns the resulting Point_t.
 * If the scalar value is 0, it returns a Point_t with all coordinates set to 0.
 *
 * @param A The Point_t operand.
 * @param scal The scalar value to divide each coordinate of the Point_t A.
 * @return The result of the division operation as a new Point_t.
 */
Point_t PT_div(Point_t, float);

/**
 * @brief Merges two Point_t by taking their average and returns the result.
 *
 * This function takes the average of Point_t A and Point_t B (by adding them and dividing by 2) and returns the resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the merging operation as a new Point_t.
 */
Point_t PT_merge(Point_t A, Point_t B);

/**
 * @brief Checks if two Point_t are equal.
 *
 * This function checks if two Point_t (A and B) are equal by comparing their x, y, and z coordinates.
 * If they are equal, the function returns 1; otherwise, it returns 0.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return 1 if Point_t A and Point_t B are equal; otherwise, 0.
 */
int PT_equal(Point_t A, Point_t B);

/**
 * @brief Creates a Point3D from a Point_t.
 *
 * This function creates a new Point3D from the given Point_t "p" by converting its float coordinates to integer indices
 * based on the START_GRID and LENGTH_GRID constants. The function returns the created Point3D.
 *
 * @param p The Point_t from which to create the Point3D.
 * @return The created Point3D with integer indices.
 */
Point3D createPoint3D(Point_t p);

/**
 * @brief Creates a Point_t from a Point3D.
 *
 * This function creates a new Point_t from the given Point3D "point" by converting its integer indices to float coordinates
 * based on the START_GRID and LENGTH_GRID constants. The function returns the created Point_t.
 *
 * @param point The Point3D from which to create the Point_t.
 * @return The created Point_t with float coordinates.
 */
Point_t createPoint_t(Point3D point);

//Liste
/**
 * @brief Initializes a new list.
 *
 * This function initializes a new list by setting its elements to NULL and size to 0.
 *
 * @param l Pointer to the list to initialize.
 */
void LST_init(List_t*);

/**
 * @brief Adds allocated memory to the list.
 *
 * This function adds allocated memory to the list to accommodate new elements.
 *
 * @param l Pointer to the list to which memory is added.
 */
void LST_addAlloc(List_t*);

/**
 * @brief Returns the number of elements in the list.
 *
 * This function returns the number of elements in the list.
 *
 * @param l Pointer to the list.
 * @return The number of elements in the list.
 */
unsigned LST_nbElements(List_t*);

/**
 * @brief Returns the index of the first available free element in the list.
 *
 * This function returns the index of the first available free element in the list.
 * If no free element is available, it adds and allocates memory for more elements.
 *
 * @param l Pointer to the list.
 * @return The index of the first available free element in the list.
 */
unsigned LST_getIndiceFree(List_t*);

/**
 * @brief Returns the index of the element with the specified identifier in the list.
 *
 * This function returns the index of the element with the specified identifier in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to search for.
 * @return The index of the element with the specified identifier, or -1 if not found.
 */
unsigned LST_getIndice(List_t*, unsigned);

/**
 * @brief Checks if an element with the specified identifier exists in the list.
 *
 * This function checks if an element with the specified identifier exists in the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to check for.
 * @return 1 if the element with the specified identifier exists, 0 otherwise.
 */
unsigned LST_check(List_t*, unsigned);

/**
 * @brief Adds an element to the list if it does not exist.
 *
 * This function adds an element to the list if it does not already exist.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be added.
 */
void LST_addElement(List_t*, unsigned);

/**
 * @brief Removes an element with the specified identifier from the list.
 *
 * This function removes an element with the specified identifier from the list.
 *
 * @param l Pointer to the list.
 * @param id Identifier of the element to be removed.
 */
void LST_removeElement(List_t* , unsigned);

/**
 * @brief Creates a new integer list.
 *
 * This function allocates memory for a new integer list and initializes its fields.
 *
 * @return Pointer to the newly allocated integer list.
 */
List_t* LST_create();

/**
 * @brief Creates a copy of an integer list.
 *
 * This function creates a copy of an integer list by copying its elements.
 *
 * @param l Pointer to the integer list to be copied.
 * @return Pointer to the newly created copy of the integer list.
 */
List_t* LST_copy(List_t*);

/**
 * @brief Copies the list with an offset of the numbering of elements 
 * according to the value of the array passed as an argument.
 * 
 * @param l List to copy.
 * @param shifts Array of offsets.
 * @return (List_t*) Copied list with shifts. 
 */
List_t* LST_copyWithShift(List_t* l, int* mod_pos_nei);

/**
 * @brief Merges two integer lists.
 *
 * This function merges two integer lists, adding elements from the second list to the first.
 * The function then deletes both input lists and returns the merged list.
 *
 * @param l1 First integer list.
 * @param l2 Second integer list.
 * @return Merged integer list.
 */
List_t* LST_addList(List_t*, List_t*);

/**
 * @brief Deletes an integer list and frees the associated memory.
 *
 * This function deletes an integer list and frees the memory associated with its elements.
 *
 * @param l Pointer to the integer list to be deleted.
 */
void LST_delete(List_t*);

/**
 * @brief Initializes a new Element list.
 *
 * This function initializes a new Element list by setting its first element to NULL.
 *
 * @return Pointer to the newly allocated Element list.
 */
Element* LST_pairs_init(void);

/**
 * @brief Adds an Element to the beginning of the list.
 *
 * This function adds an Element to the beginning of the list.
 *
 * @param list Pointer to the Element list.
 * @param start Start index of the element to add.
 * @param end End index of the element to add.
 */
void LST_pairs_addElement(Element** list, int start, int end);

/**
 * @brief Adds the pairs of atoms to connect by sorting them 
 * by ascending A* distance.
 * 
 * @param s Cage in progress.
 * @param list List of pairs of atoms to connect.
 * @param start Index of the starting atom of the path in the cage.
 * @param end Index of the ending atom of the path in the cage.
 * @param voxelGrid Grid of voxelization.
 */
void LST_pairs_addElementInOrder(Shell_t* s, Element** list, int start, int end, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap);

/**
 * @brief Removes the first Element from the list.
 *
 * This function removes the first Element from the list and frees the associated memory.
 *
 * @param list Pointer to the Element list.
 */
void LST_pairs_removeFirst(Element* list);

/**
 * @brief Deletes the Element list and frees the associated memory.
 *
 * This function deletes the Element list and frees the memory associated with its elements.
 *
 * @param list Pointer to the Element list to be deleted.
 */
void LST_pairs_delete(Element* list);

/**
 * @brief Initializes a new List_m (Shell_t) list.
 *
 * This function initializes a new List_m (used with Shell_t elements) by setting its first element to NULL.
 *
 * @return Pointer to the newly allocated List_m.
 */
List_m* LSTm_init();

/**
 * @brief Adds a Shell_t element at the beginning of the list.
 *
 * This function adds a Shell_t element to the list.
 *
 * @param list Pointer to the list.
 * @param moc Pointer to the Shell_t element to add.
 */
void LSTm_addElement(List_m* list, Shell_t* moc);

/**
 * @brief Removes the first Shell_t element from the list.
 *
 * This function removes the first Shell_t element from the list and frees the associated memory.
 *
 * @param list Pointer to the list.
 */
void LSTm_removeFirst(List_m* list);

/**
 * @brief Deletes the List_m list and frees the associated memory.
 *
 * This function deletes the List_m list and frees the memory associated with its elements.
 *
 * @param list Pointer to the List_m list to be deleted.
 */
void LSTm_delete(List_m* list);

/**
 * @brief Initializes a new List_s (Point_t) list.
 *
 * This function initializes a new List_s (used with Point_t elements) by setting its first element to NULL.
 *
 * @return Pointer to the newly allocated List_s.
 */
List_s* LSTs_init();

/**
 * @brief Adds a Point_t element at the beginning of the List_s (Point_t) list.
 *
 * This function adds a Point_t element to the List_s (Point_t) list.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param sommet Point_t element to add.
 */
void LSTs_addElement(List_s* list, Point_t sommet);

/**
 * @brief Adds a Point_t element to the List_s (Point_t) list in ascending order of A* distance.
 *
 * This function adds a Point_t element to the List_s (Point_t) list in ascending order of A* distance.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param startPos Position of the atom.
 * @param endPos Position of the objective atom.
 * @param voxelGrid Grid of voxelization.
 */
void LSTs_addElementInOrder(List_s* list, Point_t startPos, Point_t endPos, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap);

/**
 * @brief Removes the first Point_t element from the List_s (Point_t) list.
 *
 * This function removes the first Point_t element from the List_s (Point_t) list and frees the associated memory.
 *
 * @param list Pointer to the List_s (Point_t) list.
 */
void LSTs_removeFirst(List_s* list);

/**
 * @brief Deletes the List_s (Point_t) list and frees the associated memory.
 *
 * This function deletes the List_s (Point_t) list and frees the memory associated with its elements.
 *
 * @param list Pointer to the List_s (Point_t) list to be deleted.
 */
void LSTs_delete(List_s* list);

/**
 * @brief Removes a Point_t element from the List_s (Point_t) list.
 *
 * This function removes a Point_t element from the List_s (Point_t) list.
 *
 * @param list Pointer to the List_s (Point_t) list.
 * @param p Point_t element to be removed.
 */
void LSTs_removeElement(List_s* list, Point_t p);

//Point_t minDist(List_s* list, Point_t p) ;

/**
 * @brief Initializes a new List_d (integer) list.
 *
 * This function initializes a new List_d (used with integers) by setting its first element to NULL.
 *
 * @return Pointer to the newly allocated List_d.
 */
List_d* LSTd_init();

/**
 * @brief Adds an integer element at the beginning of the List_d (integer) list.
 *
 * This function adds an integer element to the List_d (integer) list.
 *
 * @param list Pointer to the List_d (integer) list.
 * @param sommet Integer element to add.
 */
void LSTd_addElement(List_d* list, int sommet);

/**
 * @brief Removes the first integer element from the List_d (integer) list.
 *
 * This function removes the first integer element from the List_d (integer) list and frees the associated memory.
 *
 * @param list Pointer to the List_d (integer) list.
 */
void LSTd_removeFirst(List_d* list);

/**
 * @brief Removes a specific integer element from the List_d (integer) list.
 *
 * This function removes a specific integer element from the List_d (integer) list.
 *
 * @param list Pointer to the List_d (integer) list.
 * @param sommet Integer element to be removed.
 */
void LSTd_removeSommet(List_d* list, int sommet);

/**
 * @brief Deletes the List_d (integer) list and frees the associated memory.
 *
 * This function deletes the List_d (integer) list and frees the memory associated with its elements.
 *
 * @param list Pointer to the List_d (integer) list to be deleted.
 */
void LSTd_delete(List_d* list);


//Graphe

/**
 * @brief Adds a neighbor to a vertex in the graph.
 *
 * This function adds a neighbor to a specified vertex in the graph.
 *
 * @param v Pointer to the vertex to which the neighbor is to be added.
 * @param id Identifier of the neighbor vertex to be added.
 */
void GPH_addNeighbor(Vertex_t*, unsigned);

/**
 * @brief Removes a neighbor from a vertex in the graph.
 *
 * This function removes a neighbor from a specified vertex in the graph.
 *
 * @param v Pointer to the vertex from which the neighbor is to be removed.
 * @param id Identifier of the neighbor vertex to be removed.
 */
void GPH_removeNeighbor(Vertex_t*, unsigned);

/**
 * @brief Returns the number of vertices in the graph.
 *
 * This function returns the number of vertices in the graph.
 *
 * @param g Pointer to the graph.
 * @return The number of vertices in the graph.
 */
int GPH_nbVertex(Graph_t*);

/**
 * @brief Adds and allocates memory for a specified number of vertices in the graph.
 *
 * This function adds and allocates memory for a specified number of vertices in the graph.
 *
 * @param g Pointer to the graph to which the vertices are to be added.
 * @param size Number of vertices to be added.
 */
void GPH_addAlloc(Graph_t*, unsigned);

/**
 * @brief Returns the index of the vertex with the specified identifier in the graph.
 *
 * This function returns the index of the vertex with the specified identifier in the graph.
 *
 * @param g Pointer to the graph.
 * @param id Identifier of the vertex.
 * @return The index of the vertex with the specified identifier, or -1 if not found.
 */
int GPH_getIndice(Graph_t* g, unsigned id);

/**
 * @brief Adds a vertex with the specified identifier to the graph if it does not exist.
 *
 * This function adds a vertex with the specified identifier to the graph if it does not exist.
 * If the vertex already exists, it returns the index of the existing vertex.
 *
 * @param g Pointer to the graph.
 * @param id Identifier of the vertex to be added.
 * @return The index of the newly added vertex or the index of the existing vertex.
 */
unsigned GPH_addVertex(Graph_t*, unsigned);

/**
 * @brief Removes a vertex with the specified identifier from the graph.
 *
 * This function removes a vertex with the specified identifier from the graph.
 *
 * @param g Pointer to the graph.
 * @param id Identifier of the vertex to be removed.
 */
void GPH_removeVertex(Graph_t*, unsigned);

/**
 * @brief Adds an edge between two vertices in the graph.
 *
 * This function adds an edge between two vertices in the graph.
 *
 * @param g Pointer to the graph.
 * @param id1 Identifier of the first vertex.
 * @param id2 Identifier of the second vertex.
 */
void GPH_addEdge(Graph_t*, unsigned, unsigned);

/**
 * @brief Removes an edge between two vertices in the graph.
 *
 * This function removes an edge between two vertices in the graph.
 *
 * @param g Pointer to the graph.
 * @param id1 Identifier of the first vertex.
 * @param id2 Identifier of the second vertex.
 */
void GPH_removeEdge(Graph_t*, unsigned, unsigned);

/**
 * @brief Determines the vertices in the graph that belong to a cycle.
 *
 * This function determines the vertices in the graph that belong to a cycle.
 * The cycle must not contain more than 6 vertices.
 *
 * @param g Pointer to the graph in which cycles are sought.
 * @return List of vertices belonging to a cycle in the graph.
 */
List_t* GPH_seekCycle(Graph_t*);

/**
 * @brief Checks if a vertex with the specified identifier exists in the graph.
 *
 * This function checks if a vertex with the specified identifier exists in the graph.
 *
 * @param g Pointer to the graph.
 * @param id Identifier of the vertex.
 * @return 1 if the vertex with the specified identifier exists, 0 otherwise.
 */
unsigned GPH_checkVertex(Graph_t*, unsigned);

/**
 * @brief Checks if an edge between two vertices exists in the graph.
 *
 * This function checks if an edge between two vertices exists in the graph.
 *
 * @param g Pointer to the graph.
 * @param id1 Identifier of the first vertex.
 * @param id2 Identifier of the second vertex.
 * @return 1 if the edge exists, 0 otherwise.
 */
unsigned GPH_checkBond(Graph_t*, unsigned, unsigned);

/**
 * @brief Creates a new graph.
 *
 * This function allocates memory for a new graph and initializes its fields.
 *
 * @return Pointer to the newly allocated graph.
 */
Graph_t* GPH_create();

/**
 * @brief Deletes a graph and frees the associated memory.
 *
 * This function deletes a graph and frees the memory associated with its vertices and edges.
 *
 * @param g Pointer to the graph to be deleted.
 */
void GPH_delete(Graph_t*);

/**
 * @brief Creates a copy of a graph.
 *
 * This function creates a copy of a graph by copying its vertices and edges.
 *
 * @param g Pointer to the graph to be copied.
 * @return Pointer to the newly created copy of the graph.
 */
Graph_t* GPH_copy(Graph_t*);

/**
 * @brief Helper function to create a new node with specified properties.
 *
 * This function creates a new Node with the given properties, including the 3D point,
 * the cost from the starting point (g), and the heuristic cost to the goal (h).
 * The function also computes the total cost (f) as the sum of g and h.
 *
 * @param point The 3D point for the node.
 * @param g The cost from the starting point to the node.
 * @param h The heuristic cost from the node to the goal.
 * @return The newly created Node with the specified properties.
 */
Node createNode(Point3D point, float g, float h);

//Molecule

/**
 * @brief Calculates the number of ligands (doublets liants) of an atom.
 *
 * The number of ligands of an atom is equal to the number of neighbors it has in the molecule.
 *
 * @param a Pointer to the Atom_t representing the atom.
 */
void MOL_nbLigands(Atom_t*);

/**
 * @brief Calculates the number of lone pairs (doublets non liants) of an atom.
 *
 * This function calculates the number of lone pairs of an atom based on its properties,
 * neighboring atoms' properties, and whether the atom belongs to a cycle.
 *
 * @param a Pointer to the Atom_t representing the atom.
 * @param alpha Average angle formed by the neighbors of the atom.
 * @param stericNeighbor Number of doublets of the neighbor of the atom (used when it has only one neighbor).
 * @param cycle Boolean indicating if the atom belongs to a cycle.
 */
void MOL_nbLonePairs(Atom_t*, float, int, unsigned);

/**
 * @brief Finds all vertices belonging to a cycle in a molecule.
 *
 * This function finds all vertices (atoms) belonging to a cycle in a molecule
 * and updates the molecule's cycle list accordingly.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 */
void MOL_seekCycle(Molecule_t*);

/**
 * @brief Calculates the number of edges (bonds) in a molecule.
 *
 * The number of edges (bonds) in a molecule is equal to half of the total number of ligands
 * (doublets liants) in all atoms of the molecule.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 * @return The number of edges (bonds) in the molecule.
 */
int MOL_nbEdges(Molecule_t*);

/**
 * @brief Adds an edge (bond) between two vertices (atoms) in a molecule.
 *
 * This function adds a bond between two vertices (atoms) in a molecule.
 * The edge is added by updating the neighborhood lists of the two atoms to include each other.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 * @param id1 Identifier of the first vertex (atom).
 * @param id2 Identifier of the second vertex (atom).
 */
void MOL_addEdge(Molecule_t*, unsigned, unsigned);

/**
 * @brief Removes an edge (bond) between two vertices (atoms) in a molecule.
 *
 * This function removes a bond between two vertices (atoms) in a molecule.
 * The edge is removed by updating the neighborhood lists of the two atoms to exclude each other.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 * @param id1 Identifier of the first vertex (atom).
 * @param id2 Identifier of the second vertex (atom).
 */
void MOL_removeEdge(Molecule_t*, unsigned, unsigned);

/**
 * @brief Finds the normal vector of a vertex (atom) in a molecule.
 *
 * This function finds the normal vector of a vertex (atom) in a molecule based on its neighbors.
 * It recursively calculates the normal vector for the vertex if it has only one neighbor.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 * @param ida Identifier of the vertex (atom) for which to find the normal vector.
 * @param dad Identifier of the parent vertex (atom) in the recursive search (used for atoms with one neighbor).
 * @return The normal vector of the vertex (atom) in the molecule.
 */
Point_t MOL_seekNormal(Molecule_t*, unsigned, unsigned);

/**
 * @brief Creates the dependency graph of a molecule.
 *
 * The dependency graph groups atoms in the molecule that can participate in hydrogen bonds.
 * An edge exists between two atoms if they cannot participate in hydrogen bonds simultaneously.
 *
 * @param m Pointer to the Molecule_t representing the molecule.
 */
void MOL_createBond(Molecule_t*);

/**
 * @brief Creates and initializes a new molecule.
 *
 * This function creates and initializes a new molecule with the specified number of vertices (atoms).
 * The atoms are created with default values, and the molecule's cycle and bond properties are set to NULL.
 *
 * @param size Number of vertices (atoms) in the molecule.
 * @return Pointer to the newly created molecule.
 */
Molecule_t* MOL_create(unsigned);

/**
 * @brief Deletes an atom from a molecule.
 *
 * Deleting an atom involves freeing the memory used by its neighborhood list.
 *
 * @param a Atom to be deleted.
 */
void MOL_deleteAtom(Atom_t*);

/**
 * @brief Deletes a molecule and its associated elements.
 *
 * Deleting a molecule involves freeing the memory used by its atoms, cycle list, and bond graph.
 *
 * @param m Pointer to the Molecule_t representing the molecule to be deleted.
 */
void MOL_delete(Molecule_t*);

//Shell

/**
 * @brief Gets the number of neighbors in the neighborhood of an AtomShl_t.
 *
 * This function returns the number of neighbors in the neighborhood of the given AtomShl_t.
 *
 * @param a Pointer to the AtomShl_t.
 * @return The number of neighbors in the AtomShl_t's neighborhood.
 */
int SHL_nbNeighborhood(AtomShl_t*);

/**
 * @brief Get the number of atoms in the Shell structure.
 *
 * This function returns the number of atoms in the specified Shell_t structure.
 *
 * @param s Pointer to the Shell_t structure.
 * @return The number of atoms in the Shell structure.
 */
int SHL_nbAtom(Shell_t*);

/**
 * @brief Get the number of edges in the Shell structure.
 *
 * This function returns the number of edges (connections between atoms) in the specified Shell_t structure.
 *
 * @param s Pointer to the Shell_t structure.
 * @return The number of edges in the Shell structure.
 */
int SHL_nbEdges(Shell_t*);

/**
 * @brief Adds an edge (bond) between two vertices with the specified IDs.
 *
 * This function adds an edge (bond) between two vertices (atoms) with the given IDs to the Shell data structure.
 * It also adds the reverse edge to create an undirected graph representation.
 *
 * @param s Pointer to the Shell data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void SHL_addEdge(Shell_t*, unsigned, unsigned);

/**
 * @brief Removes an edge (bond) between two vertices with the specified IDs.
 *
 * This function removes an edge (bond) between two vertices (atoms) with the given IDs from the Shell data structure.
 * It also removes the reverse edge to maintain the undirected graph representation.
 *
 * @param s Pointer to the Shell data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void SHL_removeEdge(Shell_t*, unsigned, unsigned);

/**
 * @brief Adds a new atom to the Shell with the specified coordinates and parent atom ID.
 *
 * This function adds a new atom to the Shell data structure with the specified coordinates and parent atom ID.
 * It also assigns a unique ID (indice) to the new atom.
 *
 * @param s Pointer to the Shell data structure.
 * @param coords The coordinates (Point_t) of the new atom.
 * @param parent The ID of the parent atom for the new atom.
 * @return The unique ID (indice) assigned to the newly added atom.
 */
unsigned SHL_addAtom(Shell_t*, Point_t, unsigned);

/**
 * @brief Removes an atom from the Shell with the specified ID.
 *
 * This function removes an atom with the given ID from the Shell data structure.
 * It also updates the neighborhood and other properties accordingly.
 *
 * @param s Pointer to the Shell data structure.
 * @param id The ID of the atom to be removed.
 */
void SHL_removeAtom(Shell_t*, unsigned);

/**
 * @brief Adds a vertex to the Shell with the specified ID.
 *
 * This function adds a vertex with the given ID to the Shell data structure.
 *
 * @param s Pointer to the Shell data structure.
 * @param id The ID of the vertex to be added.
 * @return The ID of the added vertex.
 */
unsigned SHL_addVertex(Shell_t*, unsigned);

/**
 * @brief Removes a vertex from the Shell with the specified ID.
 *
 * This function removes a vertex with the given ID from the Shell data structure.
 *
 * @param s Pointer to the Shell data structure.
 * @param id The ID of the vertex to be removed.
 */
void SHL_removeVertex(Shell_t*, unsigned);

/**
 * @brief Adds a bond between two vertices with the specified IDs.
 *
 * This function adds a bond between two vertices (atoms) with the given IDs to the Shell data structure.
 *
 * @param s Pointer to the Shell data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void SHL_addBond(Shell_t*, unsigned, unsigned);

/**
 * @brief Removes a bond between two vertices with the specified IDs.
 *
 * This function removes a bond between two vertices (atoms) with the given IDs from the Shell data structure.
 *
 * @param s Pointer to the Shell data structure.
 * @param id1 The ID of the first vertex (atom).
 * @param id2 The ID of the second vertex (atom).
 */
void SHL_removeBond(Shell_t*, unsigned, unsigned);
//void SHL_avoir2(Shell_t*, List_t*, List_t*);
//void SHL_linkBorder(Shell_t*, unsigned, List_t*);

/**
 * @brief Adds a cycle to the Shell with the specified atom ID.
 *
 * This function adds a cycle to the Shell data structure with the specified atom ID.
 *
 * @param s Pointer to the Shell data structure.
 * @param id The ID of the atom to be added to the cycle.
 */
void SHL_addCycle(Shell_t*, unsigned);

/**
 * @brief Merges two atoms in the Shell.
 *
 * This function merges (combines) two atoms in the Shell data structure.
 * It moves all the edges and properties from the eaten atom to the eater atom,
 * and then removes the eaten atom from the Shell.
 *
 * @param s Pointer to the Shell data structure.
 * @param eater The ID of the atom that will "eat" (absorb) the other atom.
 * @param eaten The ID of the atom that will be "eaten" (absorbed) by the other atom.
 */
void SHL_mergeAtom(Shell_t*, unsigned, unsigned);
//void SHL_testDis(Shell_t*);

/**
 * @brief Creates a new Shell data structure.
 *
 * This function creates and initializes a new Shell data structure.
 *
 * @return A pointer to the newly created Shell data structure.
 */
Shell_t* SHL_create();

/**
 * @brief Creates a deep copy of a Shell data structure.
 *
 * This function creates a deep copy of the given Shell data structure.
 *
 * @param s Pointer to the original Shell data structure to be copied.
 * @return A pointer to the newly created deep copy of the Shell data structure.
 */
Shell_t* SHL_copy(Shell_t*);

/**
 * @brief Creates a trimmed copy of a Shell data structure.
 *
 * This function creates a trimmed copy of the given Shell data structure by removing unused atoms or atoms belonging to the envelope.
 *
 * @param s Pointer to the original Shell data structure containing unused atoms.
 * @return A pointer to the newly created trimmed copy of the Shell data structure.
 */
Shell_t* SHL_copyCageAtoms(Shell_t* s);
//Shell_t* SHL_avoir(Shell_t*);

/**
 * @brief Deletes a Shell data structure.
 *
 * This function deletes the entire Shell data structure, including all atoms and the graph representation.
 *
 * @param s Pointer to the Shell data structure to be deleted.
 */
void SHL_delete(Shell_t*);

/**
 * @brief Deletes an atom in the Shell.
 *
 * This function deletes an atom (node) in the Shell data structure.
 *
 * @param a Pointer to the atom to be deleted.
 */
void SHL_deleteAtom(AtomShl_t* a);

/**
 * @brief Convert a molecule to a graph representation.
 *
 * This function takes a Molecule_t pointer and converts it to a Graph_t pointer.
 * The molecule's atoms and their neighboring atoms are used to create vertices and edges in the graph.
 *
 * @param m The Molecule_t pointer representing the molecule.
 * @return The Graph_t pointer representing the converted graph.
 */
Graph_t* MolToGph(Molecule_t*);

/**
 * @brief Convert a shell structure to a graph representation.
 *
 * This function takes a Shell_t pointer and converts it to a Graph_t pointer.
 * The shell's atoms and their neighboring atoms (with defined flags) are used to create
 * vertices and edges in the graph.
 *
 * @param s The Shell_t pointer representing the shell structure.
 * @return The Graph_t pointer representing the converted graph.
 */
Graph_t* ShlToGph(Shell_t*);

/**
 * @brief Allocates memory for a new Ashape_t structure.
 *
 * This function allocates memory for a new Ashape_t structure and initializes its fields.
 *
 * @return A pointer to the newly allocated Ashape_t structure.
 */
Ashape_t* ASP_create();

/**
 * @brief Frees memory used by the Ashape_t structure.
 *
 * This function frees the memory used by the Ashape_t structure and its associated arrays.
 *
 * @param as3d Pointer to the Ashape_t structure to be destroyed.
 */
void ASP_delete(Ashape_t*);

/**
 * @brief Gets the index of the next free moc in the Main_t.
 *
 * This function searches for the next free moc in the Main_t structure and
 * returns its index. If there are no free mocs, it allocates additional memory
 * to accommodate more mocs.
 *
 * @param m Pointer to the Main_t structure.
 * @return Index of the next free moc in the mocs array.
 */
unsigned MN_getIndiceFree(Main_t* m);

/**
 * @brief Gets the index of the next free moc in the Main_t (adding only one moc).
 *
 * This function searches for the next free moc in the Main_t structure and
 * returns its index. If there are no free mocs, it allocates memory for one more moc.
 *
 * @param m Pointer to the Main_t structure.
 * @return Index of the next free moc in the mocs array.
 */
unsigned MN_getIndiceFree2(Main_t* m);

/**
 * @brief Copies a Shell_t (moc) and adds it to the Main_t.
 *
 * This function copies a given Shell_t (moc) and adds it to the Main_t structure.
 * It returns the index of the added moc in the mocs array.
 *
 * @param m Pointer to the Main_t structure.
 * @param s Pointer to the Shell_t (moc) to be copied and added.
 * @return Index of the added moc in the mocs array.
 */
unsigned MN_copyMoc(Main_t*, Shell_t*);

/**
 * @brief Creates a new Main_t structure.
 *
 * This function allocates memory for a new Main_t structure and initializes its
 * elements (substrate, envelope, envarom, mocs) to NULL or zero.
 *
 * @return Pointer to the newly allocated Main_t structure.
 */
Main_t* MN_create();

/**
 * @brief Deletes a Main_t structure and its associated elements.
 *
 * This function deletes a Main_t structure and frees the memory associated with its
 * elements (substrate, envelope, envarom, mocs).
 *
 * @param m Pointer to the Main_t structure to be deleted.
 */
void MN_delete(Main_t*);

#endif
