#include "transforms.h"

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

namespace CGL {

Vector2D operator*(const Matrix3x3 &m, const Vector2D &v) {
	Vector3D mv = m * Vector3D(v.x, v.y, 1);
	return Vector2D(mv.x / mv.z, mv.y / mv.z);
}

Matrix3x3 translate(float dx, float dy) {
	// Part 3: Fill this in.
    Matrix3x3 mat;
    for (int i = 0; i < 3; i++) {
        mat[i] = {0};
    }
    mat[0][0] += dx;
    mat[1][1] += dy;
    mat[2][2] += 1;
//	return Matrix3x3();
    return mat;
}


Matrix3x3 scale(float sx, float sy) {
	// Part 3: Fill this in.
    Matrix3x3 mat;
    for (int i = 0; i < 3; i++) {
        mat[i] = {0};
    }
    mat[0][0] = sx;
    mat[1][1] = sy;
    mat[2][2] += 1;
    return mat;
//	return Matrix3x3();
}

// The input argument is in degrees counterclockwise
Matrix3x3 rotate(float deg) {
	// Part 3: Fill this in.
    Matrix3x3 mat;
    for (int i = 0; i < 3; i++) {
        mat[i] = {0};
    }
    mat[0][0] = sin(deg);
    mat[0][1] = -1 * cos(deg);
    mat[1][0] = cos(deg);
    mat[1][1] = sin(deg);
    mat[2][2] += 1;
    return mat;

//	return Matrix3x3();
}

}
