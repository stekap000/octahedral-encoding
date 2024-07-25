#include <stdio.h>
#include <math.h>

float float_abs(float x) {
	union {float f; unsigned int i; } bits = {x};
	bits.i &= 0x7fffffff;
	return bits.f;
}

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

v3 scale(v3 v, float a) {
	return (v3){v.x*a, v.y*a, v.z*a};
}

// v is assumed to be unit.
v2 octahedral_encode(v3 v) {
	/*
	  o Intersection of direction vector with octahedron

	    d		      - direction vector
		t		      - scalar
				      
		R		      - arbitrary point on the plane of interest
		n		      - plane normal
				      
		P		      - point on the plane and in the direction d
				      
		P = t*d       - equation of line in direction d
		(P - R)*n = 0 - equation of plane with normal n that contains point R

		Joining two equation together we get:
		    (t*d - R)*n = 0  =>  t*d*n = R*n  =>  t = (R*n)/(d*n)

		Since we want to find a point of intersection between direction vector
		and octahedron, we would need to test if direction vector intersects
		one of 8 triangles/planes of octahedron.
		
		To avoid this, we can recognize symmetry of octahedron and conclude that
		if we can find scalar t for first octant (x,y and z positive), thus
		obtaining intersection point t*d in this octant, we can also do it for
		direction vectors by representing them in first octant, finding t, and
		scaling them with given t in their original form.

		To find intersection point in first octant, we need to define parameters R,n
		for the plane that contains octahedron side located in first octant.
		Notice that it is enough to test intersection with the plane and not
		necessary to test intersection with a triangle since any direction vector
		in first octant must hit this plane/side of octahedron.
		We pick most obvious parameters:
		    R = (0, 0, 1) - Any vertex from first octant side.
			n = (1, 1, 1) - If this is not intuitive, you can prove it by forming
			                two vectors to represent two edges of the given octahedron
							side, taking their cross product and seeing that resulting
							vector has all the same coordinates (so we pick 1, 1, 1
							because it is easiest to compute).

		Using these parameters in formula for t, we get:
		    t = 1 / (dx + dy + dz) - dx, dy, dz are coordinates for direction vector

		This formula works for first octant. To make it usable for other direction
		vectors, we will map those vector to first octant. Notice that this can be
		done by changing all their negative coordinates to positive ie. we can just
		use absolute value of their coordinates.

		Formula that works for arbitrary direction vector is this one:
		    t = 1 / (|dx| + |dy| + |dz|)

		Why does this work?

		Because the scale of any direction vector for intersection point just depends
		on relative relation between that vector and it's octant octahedron side, not
		on the absolute positions in 3D space. This means that we can move that vector
		and it's corresponding octahedron side anywhere in 3D space and scalar t will
		not change, as long as they are in the same relative position.

		Finally, when we obtain this t, we just scale our vector with it.

		Conceptually, we have done this:
		    Take arbitrary unit vector -> Map it to first octant -> Find scalar for
			intersection -> Map vector back to it's octant -> Scale it with fould scalar

		Of course, we don't need to literally do all these steps. As an example,
		we don't need to map vector back, since we already have that vector.
	*/
	float t = 1 / (float_abs(v.x) + float_abs(v.y) + float_abs(v.z));
	v3 intersection_point = scale(v, t);

	// TODO: Map resulting vector to 2D space.
	// TODO: Map 2D space to UV space.

	return (v2){0,0};
}

int main() {
	printf("Init\n");

	
	
	return 0;
}
