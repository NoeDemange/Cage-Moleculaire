#include "util.h"
#include <math.h>
#include "output.h"

float radianToDegre(float a) {
	return a * 180 / M_PI;
}

float degreToRadian(float a) {
	return a * M_PI / 180;
}

Point_t initPoint(float scal) {
	Point_t _new;

	_new.x = scal;
	_new.y = scal;
	_new.z = scal;

	return _new;
}

Point_t addPoint(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = A.x + B.x;
	_new.y = A.y + B.y;
	_new.z = A.z + B.z;

	return _new;
}

Point_t subPoint(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = A.x - B.x;
	_new.y = A.y - B.y;
	_new.z = A.z - B.z;

	return _new;
}

Point_t mulPoint(Point_t A, float scal) {
	Point_t _new;

	_new.x = scal * A.x;
	_new.y = scal * A.y;
	_new.z = scal * A.z;

	return _new;
}

Point_t divPoint(Point_t A, float scal) {
	Point_t _new;

	_new.x = A.x / scal;
	_new.y = A.y / scal;
	_new.z = A.z / scal;

	return _new;
}

Point_t merPoint(Point_t A, Point_t B) {

	Point_t _new;

	_new.x = (A.x + B.x)/2;
	_new.y = (A.y + B.y)/2;
	_new.z = (A.z + B.z)/2;

	return _new;
}

int EqualPoint(Point_t A, Point_t B) { //1 if equal, if not 0
	if(A.x == B.x && A.y == B.y && A.z == B.z) return 1;
	return 0;
}

//Calcul de la distance euclidienne entre deux points.
float dist(Point_t A, Point_t B) {
	return sqrt(pow((A.x - B.x), 2) + pow((A.y - B.y), 2)	+ pow((A.z - B.z), 2)); //euclidian
}

//Calcul de la distance de manhattan entre deux points.
float dist_manhattan(Point_t A, Point_t B) {
	return (fabs(A.x - B.x)+fabs(A.y - B.y)+fabs(A.z - B.z));//manhattan
}

//Calcul de la distance entre deux points (demi-périmètre d'un cercle). TEST
float dist_cercle(Point_t A, Point_t B) {
	float r = dist(A,B)/2;
	return M_PI * r;
}

//intersection ligne et sphere
void interLineSphere(Point_t Start, Point_t End, Point_t P, Point_t* x1, Point_t* x2){
		Point_t vec = vector(Start, End);
		float a = pow(vec.x,2)+pow(vec.y,2)+pow(vec.z,2);
		float b = 2*(vec.x*(Start.x-P.x)+vec.y*(Start.y-P.y)+vec.z*(Start.z-P.z));
		float c = pow(Start.x-P.x,2)+pow(Start.y-P.y,2)+pow(Start.z-P.z,2)-pow(DIST_GAP_SUBSTRATE,2);
		float delta = pow(b,2)-4*a*c;
		if(delta>0){
			float d = (-b+sqrt(delta))/(2*a);
			Point_t dir;
			dir.x = d * vec.x;
			dir.y = d * vec.y;
			dir.z = d * vec.z;
			*x1 = addPoint(Start,dir);
			d = (-b-sqrt(delta))/(2*a);
			dir.x = d * vec.x;
			dir.y = d * vec.y;
			dir.z = d * vec.z;
			*x2 = addPoint(Start,dir);
		}
}

Point_t newPointOnSphere(Point_t x1, Point_t x2, Point_t P){
	Point_t vec = vector(x1,x2);
	Point_t dir, middle;
	float d = dist(x1,x2)/2;
	dir.x = d * vec.x;
	dir.y = d * vec.y;
	dir.z = d * vec.z;
	middle = addPoint(x1,dir);
	vec = vector(P,middle);
	dir = normalization(vec,DIST_GAP_SUBSTRATE);
	return addPoint(P,dir);
}


//Calcul de la distance entre deux points avec obstacle. TEST
float dist_obstacle(Point_t Start, Point_t End, Molecule_t* sub) {
	//printf("dist init : %f\n", dist(Start,End));
	int passage = 0;
	/*Shell_t* moc = SHL_create();
	int idStart = SHL_addAtom(moc,Start,-1);
	int idEnd = SHL_addAtom(moc,End,-1);
	flag(atom(moc,idEnd)) = CARBON_F;*/

	float minDistInter;
	float distTotal = 0.0;
	int minID;
	Point_t x1, x2, minPInter, secondPInter;

	while(!EqualPoint(Start,End)){
		minDistInter = __FLT_MAX__;
		minID = -1;
		x1 = initPoint(-1);
		x2 = initPoint(-1);
		for(int i = 0; i < size(sub); i++) {
			interLineSphere(Start,End,coords(atom(sub,i)),&x1,&x2);
			if(x1.x != -1 && x1.y != -1 && x1.z != -1 && x2.x != -1 && x2.y != -1 && x2.z != -1){
				float distx1 = dist(Start,x1);
				float distx2 = dist(Start,x2);
				if(distx1<distx2){
					if(distx1<minDistInter){
						minDistInter = distx1;
						minID = i;
						minPInter = x1;
						secondPInter = x2;
					}
				}else{
					if(distx2<minDistInter){
						minDistInter = distx2;
						minID = i;
						minPInter = x2;
						secondPInter = x1;
					}
				}
			}
		}
		float distInit = dist(Start,End);
		//printf("minID = %d ,distInit : %f et minDistInter : %f\n",minID,distInit,minDistInter);
		if(minID != -1 && distInit>minDistInter && minDistInter>0.001){
			/*int idx1 = SHL_addAtom(moc,minPInter,-1);
			int idx2 = SHL_addAtom(moc,secondPInter,-1);
			flag(atom(moc,idx1)) = CYCLE_F;
			flag(atom(moc,idx2)) = CYCLE_F;*/
			Point_t deriv = newPointOnSphere(minPInter, secondPInter, coords(atom(sub,minID)));
			distTotal += dist(Start,deriv);
			/*int idDeriv = SHL_addAtom(moc,deriv,-1);
			flag(atom(moc,idDeriv)) = LINKABLE_F;
			idStart = idDeriv;*/
			Start = deriv;
			passage++;
			/*if(passage >4){
			SHL_writeMol2("Test_dist_inter.mol2",moc);
			exit(EXIT_SUCCESS);
			}*/
		}
		else{
			distTotal += distInit;
			Start = End;
		}
	}
	//printf("dist: %f\n",distTotal);
	/*if(passage>3){
		SHL_writeMol2("Test_dist.mol2",moc);
		exit(EXIT_SUCCESS);}*/
	return distTotal;
}

//Normaliser un vecteur à la longueur length.
Point_t normalization(Point_t normal, float length) {
	///A vérifier
	Point_t a;
	float z;

	z = sqrt(pow(length,2) / (pow(normal.x,2) + pow(normal.y,2) + pow(normal.z,2)));
	a.x = z * normal.x;
	a.y = z * normal.y;
	a.z = z * normal.z;
	return a;
}

/**
 * @brief Compute the angle at A of a ABC triangle. 
 * 
 * @param A The first vertex of the triangle.
 * @param B The second vertex of the triangle.
 * @param C The third vertex of the triangle.
 * @return (float) Angle value in degree. 
 */
float angle(Point_t A, Point_t B, Point_t C) {
	float AB = dist(A,B), AC = dist(A,C), BC = dist(B,C);

	return acos( (pow(AC,2)+pow(AB,2)-pow(BC,2)) / (2*AC*AB) ) * 180 / M_PI;
}

Point_t vector(Point_t A, Point_t B) {
	Point_t _new;

	_new.x = B.x - A.x;
	_new.y = B.y - A.y;
	_new.z = B.z - A.z;

	return _new;
}



//Retourne la normal d'un plan à partir de 3 points du plan.
Point_t planNormal(Point_t A, Point_t B, Point_t C) {
	Point_t normal;

	normal.x = (B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y);
  normal.y = (B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z);
  normal.z = (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);

  return normalization(normal, 1);
}

//Rotation à partir d'un vecteur rotation, d'un angle et d'un point.
//Alpha doit être en degree
Point_t rotation(Point_t vec, float alpha, Point_t A) {
	Point_t rot;
	alpha = degreToRadian(alpha);
	vec = normalization(vec, 1);

	rot.x = (vec.x*vec.x + (1 - vec.x*vec.x)*cos(alpha)) * A.x +
			(vec.x*vec.y * (1 - cos(alpha)) - vec.z * sin(alpha)) * A.y +
			(vec.x*vec.z * (1 - cos(alpha)) + vec.y * sin(alpha)) * A.z;
	rot.y = (vec.y*vec.y + (1 - vec.y*vec.y)*cos(alpha)) * A.y +
			(vec.x*vec.y * (1 - cos(alpha)) + vec.z * sin(alpha)) * A.x +
			(vec.y*vec.z * (1 - cos(alpha)) - vec.x * sin(alpha)) * A.z;
	rot.z = (vec.z*vec.z + (1 - vec.z*vec.z)*cos(alpha)) * A.z +
			(vec.x*vec.z * (1 - cos(alpha)) - vec.y * sin(alpha)) * A.x +
			(vec.y*vec.z * (1 - cos(alpha)) + vec.x * sin(alpha)) * A.y;

	return rot;
}

/**
 * @brief Add the third point to a triangular pattern around a point. 
 * 
 * @param A The point at the center of the pattern.
 * @param B The first point at the edge of the triangular pattern.
 * @param C The second point at the edge of the triangular pattern.
 * @param scal The distance between A and the new third point.
 * @return Point_t The third point at the edge to add.
 */
Point_t addThirdPoint(Point_t A, Point_t B, Point_t C, float scal) {
	Point_t normal;

	normal.x = 2 * A.x - B.x - C.x;
 	normal.y = 2 * A.y - B.y - C.y;
 	normal.z = 2 * A.z - B.z - C.z;
 	normal = normalization(normal, scal);
  			
  return addPoint(A, normal);
}

Point_t AX1E1(Point_t a, Point_t x1, float length) {

	return addPoint(a, normalization(vector(x1, a), length));
}

Point_t AX2E1(Point_t a, Point_t x1, Point_t x2, float length) {

	Point_t v1, v2;

	v1 = normalization(vector(x1, a), 1);
	v2 = normalization(vector(x2, a), 1);

	return addPoint(a, normalization(addPoint(v1, v2), length));
}

Point_t AX1E2(Point_t a, Point_t x1, Point_t normal, float length) {

	Point_t v1;

	v1 = normalization(vector(a, x1), 1);

	return addPoint(a, normalization(rotation(normal, 120, v1), length));
}

Point_t AX3E1(Point_t a, Point_t x1, Point_t x2, Point_t x3, float length) {

	Point_t v1, v2, v3;

	v1 = normalization(vector(x1, a), 1);
	v2 = normalization(vector(x2, a), 1);
	v3 = normalization(vector(x3, a), 1);

	return addPoint(a, normalization(addPoint(v1, addPoint(v2,v3)), length));
}

Point_t AX2E2(Point_t a, Point_t x1, Point_t x2, float length) {

	Point_t v1, v2, other, normal, zero = {0, 0, 0};
	float angle = 180 - (109.47/2);

	v1 = normalization(vector(a,x1), 1);
	v2 = normalization(vector(a,x2), 1);

	other = normalization(addPoint(v1,v2), 1);
	normal = normalization(planNormal(zero, planNormal(a, x1, x2), other), 1);


	return addPoint(a, normalization(rotation(normal, angle, other), length));
}

Point_t AX1E3(Point_t a, Point_t x1, Point_t normal, float length) {

	Point_t v1;

	v1 = normalization(vector(a, x1), 1);

	return addPoint(a, normalization(rotation(normal, 109.47, v1), length));
}
