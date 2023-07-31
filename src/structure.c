#include "structure.h"

/**
 * @file structure.c
 * @brief Functions for converting molecular structures to graphs and helper functions.
 *
 * This file contains functions for converting a molecule and a shell structure to a graph.
 * It also includes a helper function to create a new node. The functions here facilitate the
 * conversion of molecular structures to graph representations for further processing and analysis.
 *
 */

/**
 * @brief Convert a molecule to a graph representation.
 *
 * This function takes a Molecule_t pointer and converts it to a Graph_t pointer.
 * The molecule's atoms and their neighboring atoms are used to create vertices and edges in the graph.
 *
 * @param m The Molecule_t pointer representing the molecule.
 * @return The Graph_t pointer representing the converted graph.
 */
Graph_t* MolToGph(Molecule_t* m) {

	int i, j;
	Graph_t* g = GPH_create();
	GPH_addAlloc(g, size(m));

	for (i=0; i<size(m); i++) {

		GPH_addVertex(g, i);
		for (j=0; j<neighborhoodSize(atom(m,i)); j++) {
			if (neighbor(atom(m,i),j) != -1)
				GPH_addNeighbor(vertex(g,i), neighbor(atom(m,i),j));
		}
	}

	return g;
}

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
Graph_t* ShlToGph(Shell_t* s) {

	int i, j;
	Graph_t* g = GPH_create();

	for (i=0; i<size(s); i++) {

		if (flag(atom(s,i)) != NOT_DEF_F) {

			GPH_addVertex(g, i);
			for (j=0; j<neighborhoodSize(atom(s,i)); j++) {
				if (neighbor(atom(s,i),j) != -1)
					GPH_addNeighbor(vertex(g,i), neighbor(atom(s,i),j));
			}
		}
	}

	return g;
}

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
Node createNode(Point3D point, float g, float h) {
    Node node;
    node.point = point;
    node.g = g;
    node.h = h;
    node.f = g + h;
    return node;
}
