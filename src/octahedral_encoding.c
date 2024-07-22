#include <stdio.h>

typedef struct {
	float x, y, z;
} v3;

typedef struct {
	float x, y;
} v2;

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
