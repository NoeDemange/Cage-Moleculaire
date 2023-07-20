#include "structure.h"
#include <math.h>

Point_t PT_init(float scal) {
	Point_t _new;

	_new.x = scal;
	_new.y = scal;
	_new.z = scal;

	return _new;
}

Point_t PT_add(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = A.x + B.x;
	_new.y = A.y + B.y;
	_new.z = A.z + B.z;

	return _new;
}

Point_t PT_sub(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = A.x - B.x;
	_new.y = A.y - B.y;
	_new.z = A.z - B.z;

	return _new;
}

Point_t PT_mul(Point_t A, float scal) {
	Point_t _new;

	_new.x = scal * A.x;
	_new.y = scal * A.y;
	_new.z = scal * A.z;

	return _new;
}

Point_t PT_div(Point_t A, float scal) {
	Point_t _new;

	if (scal == 0)
		return PT_init(0);

	_new.x = A.x / scal;
	_new.y = A.y / scal;
	_new.z = A.z / scal;

	return _new;
}

Point_t PT_merge(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = (A.x + B.x)/2;
	_new.y = (A.y + B.y)/2;
	_new.z = (A.z + B.z)/2;

	return _new;
}

int PT_equal(Point_t A, Point_t B) { //1 if equal, if not 0
	if(A.x == B.x && A.y == B.y && A.z == B.z) return 1;
	return 0;
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