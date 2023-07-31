#include "structure.h"

/**
 * @file structureGph.c
 * @brief Graph operations functions for allocating, manipulating, and checking a graph.
 *
 * This file contains functions for creating, manipulating, and checking a graph data structure.
 * A graph is represented by a collection of vertices and edges, and the functions here allow for
 * the addition, removal, and copying of vertices and edges. It also includes functions for checking
 * whether a vertex or edge exists in the graph and for determining cycles within the graph.
 * 
 */

/**
 * @brief Initializes a new vertex in the graph.
 *
 * This function allocates memory for the neighborhood of the new vertex in the graph and initializes its fields.
 *
 * @param v Pointer to the vertex to initialize.
 */
void GPH_initVertex(Vertex_t* v) {

	v->id = -1;

	neighborhood(v) = LST_create();
	nbNeighbors(v) = 0;
}

/**
 * @brief Adds a neighbor to a vertex in the graph.
 *
 * This function adds a neighbor to a specified vertex in the graph.
 *
 * @param v Pointer to the vertex to which the neighbor is to be added.
 * @param id Identifier of the neighbor vertex to be added.
 */
void GPH_addNeighbor(Vertex_t* v, unsigned id) {

	LST_addElement(neighborhood(v), id);
	nbNeighbors(v)++;
}

/**
 * @brief Removes a neighbor from a vertex in the graph.
 *
 * This function removes a neighbor from a specified vertex in the graph.
 *
 * @param v Pointer to the vertex from which the neighbor is to be removed.
 * @param id Identifier of the neighbor vertex to be removed.
 */
void GPH_removeNeighbor(Vertex_t* v, unsigned id) {

	LST_removeElement(neighborhood(v), id);
	nbNeighbors(v)--;
}

/**
 * @brief Deletes a vertex from the graph.
 *
 * This function frees the memory used by a vertex in the graph.
 *
 * @param v Pointer to the vertex to be deleted.
 */
void GPH_deleteVertex(Vertex_t* v) {

	LST_delete(neighborhood(v));
}

/**
 * @brief Returns the number of vertices in the graph.
 *
 * This function returns the number of vertices in the graph.
 *
 * @param g Pointer to the graph.
 * @return The number of vertices in the graph.
 */
int GPH_nbVertex(Graph_t* g) {

	int i, cpt = 0;
	for (i=0; i<size(g); i++)
		if (id(vertex(g,i)) != -1)
			cpt++;

	return cpt;
}

/**
 * @brief Adds and allocates memory for a specified number of vertices in the graph.
 *
 * This function adds and allocates memory for a specified number of vertices in the graph.
 *
 * @param g Pointer to the graph to which the vertices are to be added.
 * @param size Number of vertices to be added.
 */
void GPH_addAlloc(Graph_t* g, unsigned size) {

	int i;

	g->vertices = realloc(g->vertices, (size(g)+size)*sizeof(Vertex_t));

	for (i=0; i<size; i++)
		GPH_initVertex(vertex(g,size(g)+i));

	size(g) += size;
}

/**
 * @brief Returns the index of the first available free vertex in the graph.
 *
 * This function returns the index of the first available free vertex in the graph.
 * If no free vertex is available, it adds and allocates memory for more vertices.
 *
 * @param g Pointer to the graph.
 * @return The index of the first available free vertex in the graph.
 */
int GPH_getIndiceFree(Graph_t* g) {

	int i;

	for (i=0; i<size(g); i++)
		if (g->vertices[i].id == -1)
			return i;

	GPH_addAlloc(g, REALLOCSIZE);
	return i;
}

/**
 * @brief Returns the index of the vertex with the specified identifier in the graph.
 *
 * This function returns the index of the vertex with the specified identifier in the graph.
 *
 * @param g Pointer to the graph.
 * @param id Identifier of the vertex.
 * @return The index of the vertex with the specified identifier, or -1 if not found.
 */
int GPH_getIndice(Graph_t* g, unsigned id) {

	int i;

	for (i=0; i<size(g); i++)
		if (g->vertices[i].id == id)
			return i;

	return -1;
}

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
unsigned GPH_addVertex(Graph_t* g, unsigned id) {

	unsigned indice = GPH_getIndice(g, id);

	if (indice == -1) {

		indice = GPH_getIndiceFree(g);
		id(vertex(g,indice)) = id;
	}

	return indice;
}

/**
 * @brief Removes a vertex with the specified identifier from the graph.
 *
 * This function removes a vertex with the specified identifier from the graph.
 *
 * @param g Pointer to the graph.
 * @param id Identifier of the vertex to be removed.
 */
void GPH_removeVertex(Graph_t* g, unsigned id) {

	int i, indice = GPH_getIndice(g, id);

	if (indice != -1) {

		Vertex_t* v = vertex(g, indice);

		for (i=0; i<nbNeighbors(v); i++)
			GPH_removeNeighbor(vertex(g, GPH_getIndice(g, neighbor(v,i))), id);
	
		GPH_deleteVertex(v);
		GPH_initVertex(v);
	}
}

/**
 * @brief Adds an edge between two vertices in the graph.
 *
 * This function adds an edge between two vertices in the graph.
 *
 * @param g Pointer to the graph.
 * @param id1 Identifier of the first vertex.
 * @param id2 Identifier of the second vertex.
 */
void GPH_addEdge(Graph_t* g, unsigned id1, unsigned id2) {

	int indice1 = GPH_getIndice(g, id1), indice2 = GPH_getIndice(g, id2);
	if (indice1 != -1 && indice2 != -1) {

		GPH_addNeighbor(vertex(g, indice1), id2);
		GPH_addNeighbor(vertex(g, indice2), id1);
	}
}

/**
 * @brief Removes an edge between two vertices in the graph.
 *
 * This function removes an edge between two vertices in the graph.
 *
 * @param g Pointer to the graph.
 * @param id1 Identifier of the first vertex.
 * @param id2 Identifier of the second vertex.
 */
void GPH_removeEdge(Graph_t* g, unsigned id1, unsigned id2) {

	int indice1 = GPH_getIndice(g, id1), indice2 = GPH_getIndice(g, id2);
	if (indice1 != -1 && indice2 != -1) {

		GPH_removeNeighbor(vertex(g, indice1), id2);
		GPH_removeNeighbor(vertex(g, indice2), id1);
	}
}

/**
 * @brief Determines recursively if a given vertex belongs to a cycle in the graph.
 *
 * This function determines recursively if a given vertex belongs to a cycle in the graph.
 * The cycle must not contain more than 6 vertices.
 *
 * @param g Pointer to the graph to which the vertex belongs.
 * @param l Pointer to the list of visited vertices during the recursion.
 * @param id Identifier of the vertex being processed.
 * @param idP Identifier of the previously processed vertex.
 * @return 1 if the vertex belongs to a cycle, 0 otherwise.
 */
int GPH_cycle(Graph_t* g, List_t* l, unsigned id, unsigned idP) {

	int i, tmp = 0;
	Vertex_t* v;

	if (id == elts(l,0)) {
		return 1;
	}
	if (LST_nbElements(l) > 5 || LST_check(l, id)) {
		return 0;
	}
	LST_addElement(l, id);
	v = vertex(g, GPH_getIndice(g, id));

	for (i = 0; i < nbNeighbors(v) && tmp == 0; i++) {
		if (neighbor(v, i) != idP) {
			tmp = GPH_cycle(g, l, neighbor(v, i), id);
		}
	}
	LST_removeElement(l, id);
	return tmp;
}

/**
 * @brief Determines the vertices in the graph that belong to a cycle.
 *
 * This function determines the vertices in the graph that belong to a cycle.
 * The cycle must not contain more than 6 vertices.
 *
 * @param g Pointer to the graph in which cycles are sought.
 * @return List of vertices belonging to a cycle in the graph.
 */
List_t* GPH_seekCycle(Graph_t* g) {

	/********** Réduction du graphe ************/
	int i;
	List_t* out = LST_create();
	List_t* l = LST_create();
	Vertex_t* v, *n;
	for (i = 0; i < size(g); i++) {
		if (nbNeighbors(vertex(g, i)) == 1)
			LST_addElement(l, id(vertex(g,i)));
	}

	while (LST_nbElements(l) != 0) {

		v = vertex(g, GPH_getIndice(g, elts(l,0)));
		n = vertex(g, GPH_getIndice(g, neighbor(v,0)));

		if (nbNeighbors(n) == 2) {
			LST_addElement(l, id(n));
		}
		GPH_removeVertex(g, elts(l,0));
		LST_removeElement(l, elts(l,0));
	}
	/***** Recherche des sommets appartenant à un cycle *******/
	for (i = 0; i < size(g); i++) {
		if (id(vertex(g,i))!=-1 && GPH_cycle(g, l, id(vertex(g,i)), -1)) {
			LST_addElement(out, id(vertex(g,i)));
		}
	}
	LST_delete(l);

	return out;
}

/**
 * @brief Checks if a vertex with the specified identifier exists in the graph.
 *
 * This function checks if a vertex with the specified identifier exists in the graph.
 *
 * @param g Pointer to the graph.
 * @param id Identifier of the vertex.
 * @return 1 if the vertex with the specified identifier exists, 0 otherwise.
 */
unsigned GPH_checkVertex(Graph_t* g, unsigned id) {

	int i;

	for (i=0; i<size(g); i++)
		if (id(vertex(g,i)) == id)
			return 1;
	
	return 0;
}

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
unsigned GPH_checkBond(Graph_t* g, unsigned id1, unsigned id2) {

	return LST_check(neighborhood(vertex(g, GPH_getIndice(g, id1))), id2);
}

/**
 * @brief Creates a new graph.
 *
 * This function allocates memory for a new graph and initializes its fields.
 *
 * @return Pointer to the newly allocated graph.
 */
Graph_t* GPH_create() {

	Graph_t* g = malloc(sizeof(Graph_t));

	g->size = 0;
	g->vertices = NULL;

	return g;
}

/**
 * @brief Deletes a graph and frees the associated memory.
 *
 * This function deletes a graph and frees the memory associated with its vertices and edges.
 *
 * @param g Pointer to the graph to be deleted.
 */
void GPH_delete(Graph_t* g) {

	int i;

	for (i=0; i<size(g); i++)
		GPH_deleteVertex(vertex(g,i));

	free(g->vertices);
	free(g);
}

/**
 * @brief Creates a copy of a graph.
 *
 * This function creates a copy of a graph by copying its vertices and edges.
 *
 * @param g Pointer to the graph to be copied.
 * @return Pointer to the newly created copy of the graph.
 */
Graph_t* GPH_copy(Graph_t* g) {

	int i;
	Graph_t* copy = GPH_create();
	unsigned indice;

	GPH_addAlloc(copy, GPH_nbVertex(g));

	for (i=0; i<size(g); i++) {
		if (id(vertex(g,i)) != -1) {
			indice = GPH_addVertex(copy, id(vertex(g,i)));
			LST_delete(neighborhood(vertex(copy, indice)));
			neighborhood(vertex(copy, indice)) = LST_copy(neighborhood(vertex(g,i)));
			nbNeighbors(vertex(copy, indice)) = nbNeighbors(vertex(g,i));
		}
	}

	return copy;
}
