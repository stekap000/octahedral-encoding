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

// In more realistic case, functions would accept one batch of vectors, and not just one vector.
// Also, they would produce one batch of mapped vectors. Number would depent on cache sizes.

v2 octahedral_encode(v3 v) {
	// v is assumed to be unit.
	
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
	v3 intersection_point = scale(v, 1 / denom);
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
}

v3 octahedral_decode(v2 v) {
	
}

v2 octahedral_map_to_uv(v2 v) {
	/*
	  o Mapping octahedral coordinates to UV space
	  
	    Here we map 2D octahedral coordinates, which are in range [-1, 1], to UV space
		coordinates, which are in range [0, 1].

		Reason for doing this is that we can easily access texture information based on
		stored coordinates.

		One example where this can be used is the following:
		    We have a complicated scene and we are focusing on rendering specific object. This object
			can have reflective and other properties. In such cases, it is usually not feasible to try
			to determine light transport in real time. For this reason, we can run simulation separately,
			determine color values for object in every direction, and there store those results. Later,
			when we are rendering this object in that particular scene, we can just query this map with
			direction vectors, get colors and render object. This kind of map from directions to colors
			is usually called environment mapping.
			
			One thing here where we can use octahedral encoding with UV coordinates is when we are decide
			to simulate object colors in certain number of directions. This way, we will get pairs of
			vectors and colors ie. ((x, y, z), (r, g, b)). We can map vectors to octahedral UV coordinates,
			thus getting ((u, v), (r, g, b)), where u and v are in range [0, 1]. Then, we can determine
			the size of texture where we will store these colors, and store them in such a way that some
			color (r1, g1, b1) with UV coordinates (u1, v1) will be placed on pixel with coordinates
			(u1*texture_width, v1*texture_height). We can fill this texture with more colors by doing
			simulations for more directions, and we can also try to interpolate between some color values.
			Regardless, we now have texture whose colors represent some object colors (in particular
			direction) in some scene.

			Now, when we are rendering this object in that scene, and we need it's color in specific
			direction, we can take that direction vector, map it to octahedral UV space, and index
			constructed texture by multiplying it's (u, v) coordinates with texture width and height.
			As a result, we get object color for requested direction.

		We don't need to store just colors in this way. We can store reflections, lighting and whatever
		else we need.
	*/
	return (v2){(v.x + 1)*0.5, (v.y + 1)*0.5};
}

int main() {

	return 0;
}
