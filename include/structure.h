#ifndef __STRUCTURE_H
#define __STRUCTURE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constant.h"

//retourne l'adresse
#define atom(o,i) ((o)->atoms+(i)) //adresse de l'atome
#define neighborhood(a) (a)->neighborhood //adresse du voisinage : list

#define size(o) (o)->size
#define cycle(o, i) LST_check((o)->cycle, (i))
#define checkVertex(o, i) GPH_checkVertex((o)->bond, (i))
#define checkBond(o, i, j)	GPH_checkBond((o)->bond, (i), (j))
#define bond(o)	(o)->bond
#define elts(l, i) (l)->elts[(i)]

// (in a loop) Browse a list of elements until it finds a -1 element (an unused value).
// Only used to scroll through list where unused values are at the end. 
#define forEachElement(l,i) (i) < size((l)) && elts((l),(i)) != -1

#define coords(a) (a)->coords
#define atomX(a) (a)->coords.x
#define atomY(a) (a)->coords.y
#define atomZ(a) (a)->coords.z
#define neighbor(a,i) (a)->neighborhood->elts[(i)]
#define neighborhoodSize(a) (a)->neighborhood->size

// (in a loop) Browse a list of neighbors until it finds a -1 element (an unused value).
// The unused values of neighbors must be at the end. 
#define forEachNeighbor(a,i) (i) < neighborhoodSize((a)) && neighbor((a),(i)) != -1

///macro molecule
#define symbol(a) (a)->info.symbol
#define radius(a) (a)->info.radius
#define ligands(a) (a)->info.ligands
#define lonePairs(a) (a)->info.lonePairs
#define steric(a) (a)->info.ligands + (a)->info.lonePairs


///macro shell
#define flag(a) (a)->flag
#define parentAtom(a) (a)->parentAtom

//macro graph
#define vertex(g, i) ((g)->vertices+(i))
#define id(v)	(v)->id
#define nbNeighbors(v) (v)->nbNeighbors

//macro main
#define substrat(m) (m)->substrat
#define envelope(m) (m)->envelope
#define envarom(m)	(m)->envarom
#define moc(m,i) (m)->mocs[i]
#define mocSize(m) (m)->mocSize

#define getName(var)  #var

typedef struct {
	float x;
	float y;
	float z;
} Point_t;

/**************************************/
/* LISTE ******************************/
/**************************************/
typedef struct {
	
	int* elts;
	unsigned size;	

} List_t;

// Liste des sommets a relier
typedef struct Element Element;
struct Element {
	
	int start;
	int end;
	float distance;
	Element *next;
};

// Liste des sommets intermediaires
typedef struct Elem_s Elem_s;
struct Elem_s {
	
	Point_t position;
	float distance;
	Elem_s *next;
};

typedef struct {
	
	Elem_s *first;
} List_s;

// Liste d'entiers
typedef struct Elem_d Elem_d;
struct Elem_d {
	
	int idAtom;
	Elem_d *next;	

};

typedef struct {
	
	Elem_d *first;	

} List_d;

/**************************************/
/* POINT ***************************/
/**************************************/

typedef struct {
    int x, y, z;
} Point3D;


/**************************************/
/* GRAPHE *****************************/
/**************************************/
typedef struct {

	unsigned id;

	// Voisinage
	List_t* neighborhood;
	unsigned nbNeighbors;
} Vertex_t;

typedef struct {
	
	Vertex_t* vertices;
	unsigned size;
} Graph_t;

typedef struct {
    Point3D point;
    float g, h, f;
} Node;

/**************************************/
/* VOXEL ***************************/
/**************************************/

typedef int*** VOXELGRID;


/**************************************/
/* MOLECULE ***************************/
/**************************************/
typedef struct {

	char symbol[2];
	int radius;
	int ligands;
	int lonePairs;

} ChemicalInfo_t;

typedef struct {

	ChemicalInfo_t info;
	Point_t coords;

	// Voisinage
	List_t* neighborhood;
} Atom_t;

typedef struct {
	
	Atom_t* atoms;
	List_t* cycle; // liste des sommets appartenant à un cycle
	Graph_t* bond;
	unsigned size;
} Molecule_t;

/**************************************/
/* SHELL ******************************/
/**************************************/
typedef struct {

	int flag;
	Point_t	coords;

	// Lien
	unsigned parentAtom;

	// Voisinage
	List_t* neighborhood;
} AtomShl_t;

typedef struct {
	
	AtomShl_t* atoms;
	List_t* cycle;
	Graph_t* bond;

	unsigned size;
} Shell_t;

/**************************************/
/* MAIN *******************************/
/**************************************/
typedef struct {

	Molecule_t* substrat;
	Shell_t* envelope;
	Shell_t* envarom;
	Shell_t** mocs;
	unsigned mocSize;

} Main_t;

/**************************************/
/* ASHAPE3D ***************************/
/**************************************/
typedef struct {
	int nb_triang;
	int nb_edge;
	int nb_vertex;
	int nb_x;
	int nb_alpha;

	double* triang;
	double* edge;
	double* vertex;
	double* x;
	double* alpha;
}Ashape_t;

/**************************************/
/* LISTE MOC **************************/
/**************************************/

// Liste des mocs a traiter
typedef struct Elem Elem;
struct Elem {
	
	Shell_t* moc;
	int nbPatterns;
	int nbCycles;
	Elem *next;	
};

typedef struct {
	
	Elem *first;	
} List_m;

/**************************************/
/* NODE HEAP **************************/
/**************************************/
typedef struct{
	float dist;
	int indexHeap; //-1 already seen
}VMap; //voxel map for pathfinding

typedef struct{
	int size;
	Node* node;
}NodeHeap;

//NodeHeap
NodeHeap NH_initAlloc();
void NH_free(NodeHeap);
void NH_insert(NodeHeap*, Node, VMap***);
Node NH_extractMin(NodeHeap*, VMap***);
void NH_decrease_priority(NodeHeap*, VMap***, Point3D, float distG);
void NH_print(NodeHeap);
VMap*** VMap_alloc();
void VMap_free(VMap***);

//Point
Point_t PT_init(float);
Point_t PT_add(Point_t, Point_t);
Point_t PT_sub(Point_t, Point_t);
Point_t PT_mul(Point_t, float);
Point_t PT_div(Point_t, float);
Point_t PT_merge(Point_t A, Point_t B);
int PT_equal(Point_t A, Point_t B);
Point3D createPoint3D(Point_t p);
Point_t createPoint_t(Point3D point);

//Liste

void LST_init(List_t*);
void LST_addAlloc(List_t*);
unsigned LST_nbElements(List_t*);
unsigned LST_getIndiceFree(List_t*);
unsigned LST_getIndice(List_t*, unsigned);
unsigned LST_check(List_t*, unsigned);
void LST_addElement(List_t*, unsigned);
void LST_removeElement(List_t* , unsigned);
List_t* LST_create();
List_t* LST_copy(List_t*);
List_t* LST_copyWithShift(List_t* l, int* mod_pos_nei);
List_t* LST_addList(List_t*, List_t*);
void LST_delete(List_t*);

Element* LST_pairs_init(void);
void LST_pairs_addElement(Element** list, int start, int end);
void LST_pairs_addElementInOrder(Shell_t* s, Element** list, int start, int end, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap);
void LST_pairs_removeFirst(Element* list);
void LST_pairs_delete(Element* list);

List_m* LSTm_init();
void LSTm_addElement(List_m* list, Shell_t* moc);
void LSTm_removeFirst(List_m* list);
void LSTm_delete(List_m* list);

List_s* LSTs_init();
void LSTs_addElement(List_s* list, Point_t sommet);
void LSTs_addElementInOrder(List_s* list, Point_t startPos, Point_t endPos, VOXELGRID voxelGrid, VMap*** vMap, NodeHeap nodeHeap);
void LSTs_removeFirst(List_s* list);
void LSTs_delete(List_s* list);
void LSTs_removeElement(List_s* list, Point_t p);
Point_t minDist(List_s* list, Point_t p) ;

List_d* LSTd_init();
void LSTd_addElement(List_d* list, int sommet);
void LSTd_removeFirst(List_d* list);
void LSTd_removeSommet(List_d* list, int sommet);
void LSTd_delete(List_d* list);


//Graphe

void GPH_addNeighbor(Vertex_t*, unsigned);
void GPH_removeNeighbor(Vertex_t*, unsigned);
int GPH_nbVertex(Graph_t*);
void GPH_addAlloc(Graph_t*, unsigned);
int GPH_getIndice(Graph_t* g, unsigned id);
unsigned GPH_addVertex(Graph_t*, unsigned);
void GPH_removeVertex(Graph_t*, unsigned);
void GPH_addEdge(Graph_t*, unsigned, unsigned);
void GPH_removeEdge(Graph_t*, unsigned, unsigned);
List_t* GPH_seekCycle(Graph_t*);
unsigned GPH_checkVertex(Graph_t*, unsigned);
unsigned GPH_checkBond(Graph_t*, unsigned, unsigned);
Graph_t* GPH_create();
void GPH_delete(Graph_t*);
Graph_t* GPH_copy(Graph_t*);
Node createNode(Point3D point, float g, float h);

//Molecule

void MOL_nbLigands(Atom_t*);
void MOL_nbLonePairs(Atom_t*, float, int, unsigned);
void MOL_seekCycle(Molecule_t*);
int MOL_nbEdges(Molecule_t*);
void MOL_addEdge(Molecule_t*, unsigned, unsigned);
void MOL_removeEdge(Molecule_t*, unsigned, unsigned);
Point_t MOL_seekNormal(Molecule_t*, unsigned, unsigned);
void MOL_createBond(Molecule_t*);
Molecule_t* MOL_create(unsigned);
void MOL_deleteAtom(Atom_t*);
void MOL_delete(Molecule_t*);

int SHL_nbNeighborhood(AtomShl_t*);
int SHL_nbAtom(Shell_t*);
int SHL_nbEdges(Shell_t*);
void SHL_addEdge(Shell_t*, unsigned, unsigned);
void SHL_removeEdge(Shell_t*, unsigned, unsigned);
unsigned SHL_addAtom(Shell_t*, Point_t, unsigned);
void SHL_removeAtom(Shell_t*, unsigned);
unsigned SHL_addVertex(Shell_t*, unsigned);
void SHL_removeVertex(Shell_t*, unsigned);
void SHL_addBond(Shell_t*, unsigned, unsigned);
void SHL_removeBond(Shell_t*, unsigned, unsigned);
//void SHL_avoir2(Shell_t*, List_t*, List_t*);
void SHL_linkBorder(Shell_t*, unsigned, List_t*);
void SHL_addCycle(Shell_t*, unsigned);
void SHL_mergeAtom(Shell_t*, unsigned, unsigned);
//void SHL_testDis(Shell_t*);
Shell_t* SHL_create();
Shell_t* SHL_copy(Shell_t*);
Shell_t* SHL_copyCageAtoms(Shell_t* s);
//Shell_t* SHL_avoir(Shell_t*);
void SHL_delete(Shell_t*);
void SHL_deleteAtom(AtomShl_t* a);

Graph_t* MolToGph(Molecule_t*);
Graph_t* ShlToGph(Shell_t*);

Ashape_t* ASP_create();
void ASP_delete(Ashape_t*);

unsigned MN_getIndiceFree(Main_t* m);
unsigned MN_getIndiceFree2(Main_t* m);
unsigned MN_copyMoc(Main_t*, Shell_t*);
Main_t* MN_create();
void MN_delete(Main_t*);

#endif
