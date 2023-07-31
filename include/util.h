#ifndef __UTIL_H
#define __UTIL_H

#include "structure.h"

/**
 * @file util.h
 * @brief Utility Functions Header File
 *
 * This file contains function prototypes for various utility functions used in the project.
 */

/**
 * @brief Convert radians to degrees.
 *
 * This function converts an angle from radians to degrees.
 *
 * @param a The angle in radians.
 * @return The angle in degrees.
 */
float radianToDegre(float);

/**
 * @brief Convert degrees to radians.
 *
 * This function converts an angle from degrees to radians.
 *
 * @param a The angle in degrees.
 * @return The angle in radians.
 */
float degreToRadian(float);

/**
 * @brief Compute the Euclidean distance between two points.
 *
 * This function calculates the Euclidean distance between two 3D points A and B.
 *
 * @param A The first point.
 * @param B The second point.
 * @return The Euclidean distance between A and B.
 */
float dist(Point_t, Point_t);

/**
 * @brief Compute the Manhattan distance between two points.
 *
 * This function calculates the Manhattan distance between two 3D points A and B.
 * The Manhattan distance is the sum of the absolute differences in the x, y, and z coordinates.
 *
 * @param A The first point.
 * @param B The second point.
 * @return The Manhattan distance between A and B.
 */
float dist_manhattan(Point_t, Point_t);
//float dist_cercle(Point_t, Point_t);
//float dist_obstacle(Point_t A, Point_t B, Molecule_t* sub);

/**
 * @brief Normalize a vector to a specified length.
 *
 * This function normalizes a given vector to the specified length. The resulting vector will have
 * the same direction as the original vector but with a magnitude equal to the specified length.
 *
 * @param normal The vector to be normalized (Point_t structure).
 * @param length The desired length of the normalized vector.
 * @return The normalized vector (Point_t structure).
 */
Point_t normalization(Point_t, float);

/**
 * @brief Compute the angle at A of a ABC triangle. 
 * 
 * * This function calculates the angle at vertex A of a triangle ABC using the law of cosines.
 * 
 * @param A The first vertex of the triangle.
 * @param B The second vertex of the triangle.
 * @param C The third vertex of the triangle.
 * @return The angle at vertex A in degrees. 
 */
float angle(Point_t, Point_t, Point_t);

/**
 * @brief Compute the vector between two points.
 *
 * This function calculates the vector (direction) between two 3D points A and B.
 *
 * @param A The starting point.
 * @param B The ending point.
 * @return The vector from A to B.
 */
Point_t vector(Point_t, Point_t);

/**
 * @brief Add the third point to a triangular pattern around a point. 
 * 
 * This function adds the third point to a triangular pattern around a central point A. The pattern
 * is defined by two other points B and C, and the distance between A and the new point is specified by scal.
 *
 * @param A The point at the center of the pattern.
 * @param B The first point at the edge of the triangular pattern.
 * @param C The second point at the edge of the triangular pattern.
 * @param scal The distance between A and the new third point.
 * @return Point_t The third point at the edge to add.
 */
Point_t addThirdPoint(Point_t, Point_t, Point_t, float);

/**
 * @brief Compute the normal vector of a plane defined by three points.
 *
 * This function calculates the normal vector of a plane defined by three 3D points A, B, and C.
 *
 * @param A The first point of the plane.
 * @param B The second point of the plane.
 * @param C The third point of the plane.
 * @return The normal vector of the plane.
 */
Point_t planNormal(Point_t, Point_t, Point_t);

/**
 * @brief Perform rotation around a point.
 *
 * This function rotates a point A around a given vector vec by a specified angle alpha.
 *
 * @param vec The vector used for rotation.
 * @param alpha The angle of rotation in degrees.
 * @param A The point to be rotated.
 * @return The rotated point.
 */
Point_t rotation(Point_t, float, Point_t);


/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of vector x1.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized vector x1.
 *
 * @param a The starting point.
 * @param x1 The vector direction.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t AX1E1(Point_t, Point_t, float);

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of the sum of vectors x1 and x2.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized sum of vectors x1 and x2.
 *
 * @param a The starting point.
 * @param x1 The first vector direction.
 * @param x2 The second vector direction.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t AX2E1(Point_t, Point_t, Point_t, float);

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of the rotated vector x1.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized vector x1 rotated around the specified normal vector.
 *
 * @param a The starting point.
 * @param x1 The vector direction to rotate.
 * @param normal The normal vector around which to perform the rotation.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t AX1E2(Point_t, Point_t, Point_t, float);

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of the sum of vectors x1, x2, and x3.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized sum of vectors x1, x2, and x3.
 *
 * @param a The starting point.
 * @param x1 The first vector direction.
 * @param x2 The second vector direction.
 * @param x3 The third vector direction.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t AX3E1(Point_t, Point_t, Point_t, Point_t, float);

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of the sum of vectors x1 and x2 and performing a rotation.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized sum of vectors x1 and x2, and then performing a rotation around a calculated
 * normal vector.
 *
 * @param a The starting point.
 * @param x1 The first vector direction.
 * @param x2 The second vector direction.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t AX2E2(Point_t, Point_t, Point_t, float);

/**
 * @brief Compute the point obtained by moving a distance length from point a in the direction of the rotated vector x1.
 *
 * This function computes the point obtained by moving a distance length from point a in the direction
 * of the normalized vector x1 rotated around the specified normal vector.
 *
 * @param a The starting point.
 * @param x1 The vector direction to rotate.
 * @param normal The normal vector around which to perform the rotation.
 * @param length The distance to move.
 * @return The new point after the move.
 */
Point_t AX1E3(Point_t, Point_t, Point_t, float);

#endif
