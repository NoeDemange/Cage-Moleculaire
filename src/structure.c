#include "structure.h"

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

// Helper function to create a new node
Node createNode(Point3D point, int g, int h) {
    Node node;
    node.point = point;
    node.g = g;
    node.h = h;
    node.f = g + h;
    return node;
}

// Helper function to create a new Point3D from Point_t
Point3D createPoint3D(Point_t p) {
    Point3D point;
    point.x = (int)((abs(START_GRID_X) + (p.x))/(LENGTH_GRID_X));
    if(point.x<0){ 
		printf("grid too small\n");
		point.x = 0;
	} //Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    if(point.x>=GRID_SIZE_X) {
		printf("grid too small\n");
		point.x = GRID_SIZE_X-1;
	}//Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    point.y = (int)((abs(START_GRID_Y) + (p.y))/(LENGTH_GRID_Y));
    if(point.y<0) {
		printf("grid too small\n");
		point.y = 0;
	}//Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    if(point.y>=GRID_SIZE_Y) {
		printf("grid too small\n");
		point.y = GRID_SIZE_Y-1;
	}//Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    point.z = (int)((abs(START_GRID_Z) + (p.z))/(LENGTH_GRID_Z));
    if(point.z<0) {
		printf("grid too small\n");
		point.z = 0;
	}//Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    if(point.z>=GRID_SIZE_Z) {
		printf("grid too small\n");
		point.z = GRID_SIZE_Z-1;
	}//Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    return point;
}

// Helper function to create a new Point_t from Point3D
Point_t createPoint_t(Point3D point) {
    Point_t p;
    p.x = (START_GRID_X + (point.x*LENGTH_GRID_X));
    p.y = (START_GRID_Y + (point.y*LENGTH_GRID_Y));
    p.z = (START_GRID_Z + (point.z*LENGTH_GRID_Z));
    return p;
}