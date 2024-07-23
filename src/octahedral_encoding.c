#include <stdio.h>
#include <math.h>

typedef struct {
	float x, y, z;
} v3;

typedef struct {
	float x, y;
} v2;

float v3_norm(v3 v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

v3 v3_unit(v3 v) {
	float norm = v3_norm(v);
	return (v3){v.x/norm, v.y/norm, v.z/norm};
}

float v3_dot(v3 v, v3 w) {
	return v.x*w.x + v.y*v.y + v.z*w.z;
}

v2 octahedral_encode(v3 v) {
	v2 w = {0};

	(void)v;

	// TODO: Map unit vector to octahedron.
	// TODO: Map resulting vector to 2D space.
	// TODO: Map 2D space to UV space.

	return w;
}

int main() {
	printf("Init\n");
	
	return 0;
}
