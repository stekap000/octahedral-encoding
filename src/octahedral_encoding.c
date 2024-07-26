#include <stdio.h>
#include <math.h>

float float_abs(float x) {
	union {float f; unsigned int i; } bits = {x};
	bits.i &= 0x7fffffff;
	return bits.f;
}

float float_sign(float x) {
	return (x > 0) - (x < 0);
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

		Octahedron that we use is the one with vertices at:
		    (1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)
		
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
	float denom = float_abs(v.x) + float_abs(v.y) + float_abs(v.z);
	if(denom == 0) return (v2){0, 0};
	float t = 1 / denom;
	v3 intersection_point = scale(v, t);
	/*
	  o Mapping intersection points to 2D space

	    We will project upper octahedron portion to xy plane, and then patch 2D map
		into square by adjusting it with lower octahedron portion.

		First, we project intersection points on the upper sides to xy plane by setting
		z = 0. At this point, our 2D map looks like a unit diamond with vertices at
		(1, 0), (0, 1), (-1, 0), (0, -1).

		Next, we will complete this map into a square, by unfolding lower sides of the
		octahedron to their corresponding xy plane regions, where they are connected
		with their corresponding upper side. Then, we will just scale them, so that
		they complete map into a square that stretches from -1 to 1 on both axes.
		Computationally, we can do this by projecting intersection points from lower
		sides to xy plane, thus getting 4 triangles for those 4 sides, and then
		unfolding them by rotating them for 180 deg around their corresponding
		diamond map sides (or equivalently, making their reflection in mentioned sides).
		
		To do this in xy plane, with the reflection in side y = -x + 1, we will map
		(x, y) point to (1-y, 1-x) point.

		This mapping is completely valid, but it requires temp save when new x and y coordinates
		are calculated. This is because we need to store old x before we calculate new one since
		calculation of new y requires old x. Because of this, we can imagine additional mapping
		that we apply after mentioned rotation that rotates given triangle in first quadrant around
		y = x line. This way, point (x,y) gets mapped into (1-x, 1-y) which does not require temp
		save in calculation since new x coordinate only depends on previous x, and new y depends
		only on previous y.

		Of course, we use symmetry so that we can use the same calculation for all
		sides in the similar way as for finding intersection points. For this purpose,
		we again use absolute function. After calculation in the first quadrant, we just move
		result into specific quadrant corresponding to specific point by using signs of x and y.
	*/
	if(intersection_point.z >= 0) return (v2){intersection_point.x, intersection_point.y};

	return (v2){
		intersection_point.x = float_sign(intersection_point.x)*(1 - float_abs(intersection_point.x)),
		intersection_point.y = float_sign(intersection_point.y)*(1 - float_abs(intersection_point.y))
	};

	// TODO: Map 2D space to UV space.
}

int main() {
	printf("Init\n");

	
	
	return 0;
}
