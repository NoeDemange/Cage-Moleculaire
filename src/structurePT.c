#include "structure.h"
#include <math.h>

/**
 * @file structurePT.c
 * @brief Point_t Data Structure Implementation and related operations
 *
 * This file contains the implementation of the Point_t data structure and related operations.
 * Point_t represents a 3D point with float coordinates (x, y, z).
 */

/**
 * @brief Initializes a new Point_t with equal scalar values.
 *
 * This function initializes a new Point_t with the same scalar value for all coordinates (x, y, z).
 *
 * @param scal The scalar value to set for all coordinates of the Point_t.
 * @return The initialized Point_t.
 */
Point_t PT_init(float scal) {
	Point_t _new;

	_new.x = scal;
	_new.y = scal;
	_new.z = scal;

	return _new;
}

/**
 * @brief Adds two Point_t together and returns the result.
 *
 * This function adds two Point_t (A and B) together and returns the resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the addition operation as a new Point_t.
 */
Point_t PT_add(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = A.x + B.x;
	_new.y = A.y + B.y;
	_new.z = A.z + B.z;

	return _new;
}

/**
 * @brief Subtracts one Point_t from another and returns the result.
 *
 * This function subtracts Point_t B from Point_t A and returns the resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the subtraction operation as a new Point_t.
 */
Point_t PT_sub(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = A.x - B.x;
	_new.y = A.y - B.y;
	_new.z = A.z - B.z;

	return _new;
}

/**
 * @brief Multiplies a Point_t by a scalar value and returns the result.
 *
 * This function multiplies each coordinate of the Point_t A by the scalar value "scal" and returns the resulting Point_t.
 *
 * @param A The Point_t operand.
 * @param scal The scalar value to multiply each coordinate of the Point_t A.
 * @return The result of the multiplication operation as a new Point_t.
 */
Point_t PT_mul(Point_t A, float scal) {
	Point_t _new;

	_new.x = scal * A.x;
	_new.y = scal * A.y;
	_new.z = scal * A.z;

	return _new;
}

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
Point_t PT_div(Point_t A, float scal) {
	Point_t _new;

	if (scal == 0)
		return PT_init(0);

	_new.x = A.x / scal;
	_new.y = A.y / scal;
	_new.z = A.z / scal;

	return _new;
}

int PT_compare(Point_t A, Point_t B) {
	return A.x == B.x && A.y == B.y && A.z == B.z;
}

/**
 * @brief Merges two Point_t by taking their average and returns the result.
 *
 * This function takes the average of Point_t A and Point_t B (by adding them and dividing by 2) and returns the resulting Point_t.
 *
 * @param A The first Point_t operand.
 * @param B The second Point_t operand.
 * @return The result of the merging operation as a new Point_t.
 */
Point_t PT_merge(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = (A.x + B.x)/2;
	_new.y = (A.y + B.y)/2;
	_new.z = (A.z + B.z)/2;

	return _new;
}

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
int PT_equal(Point_t A, Point_t B) { //1 if equal, if not 0
	if(A.x == B.x && A.y == B.y && A.z == B.z) return 1;
	return 0;
}

/**
 * @brief Creates a Point3D from a Point_t.
 *
 * This function creates a new Point3D from the given Point_t "p" by converting its float coordinates to integer indices
 * based on the START_GRID and LENGTH_GRID constants. The function returns the created Point3D.
 *
 * @param p The Point_t from which to create the Point3D.
 * @return The created Point3D with integer indices.
 */
Point3D createPoint3D(Point_t p) {
    Point3D point;
    point.x = (int)((abs(START_GRID) + (p.x))/(LENGTH_GRID));
    if(point.x<0){ 
		printf("grid too small\n");
		point.x = 0;
	} //Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    if(point.x>=GRID_SIZE) {
		printf("grid too small\n");
		point.x = GRID_SIZE-1;
	}//Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    point.y = (int)((abs(START_GRID) + (p.y))/(LENGTH_GRID));
    if(point.y<0) {
		printf("grid too small\n");
		point.y = 0;
	}//Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    if(point.y>=GRID_SIZE) {
		printf("grid too small\n");
		point.y = GRID_SIZE-1;
	}//Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    point.z = (int)((abs(START_GRID) + (p.z))/(LENGTH_GRID));
    if(point.z<0) {
		printf("grid too small\n");
		point.z = 0;
	}//Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    if(point.z>=GRID_SIZE) {
		printf("grid too small\n");
		point.z = GRID_SIZE-1;
	}//Pour éviter d'être en dehors de la grille peut-être à changer par autre méthode meilleur
    return point;
}

/**
 * @brief Creates a Point_t from a Point3D.
 *
 * This function creates a new Point_t from the given Point3D "point" by converting its integer indices to float coordinates
 * based on the START_GRID and LENGTH_GRID constants. The function returns the created Point_t.
 *
 * @param point The Point3D from which to create the Point_t.
 * @return The created Point_t with float coordinates.
 */
Point_t createPoint_t(Point3D point) {
    Point_t p;
    p.x = (START_GRID + (point.x*LENGTH_GRID));
    p.y = (START_GRID + (point.y*LENGTH_GRID));
    p.z = (START_GRID + (point.z*LENGTH_GRID));
    return p;
}