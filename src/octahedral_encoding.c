#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float float_abs(float x) {
	union {float f; unsigned int i; } bits = {x};
	bits.i &= 0x7fffffff;
	return bits.f;
}

float rand_float() {
	return ((rand() / (float)RAND_MAX) - 0.5) * 2;
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

v3 v3_scale(v3 v, float a) {
	return (v3){v.x*a, v.y*a, v.z*a};
}

// In more realistic case, functions would accept one batch of vectors, and not just one vector.
// Also, they would produce one batch of mapped vectors. Number would depent on cache sizes.

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
	v = v3_unit(v);
	float denom = float_abs(v.x) + float_abs(v.y) + float_abs(v.z);
	if(denom == 0) return (v2){0, 0};
	v3 intersection_point = v3_scale(v, 1 / denom);
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
	if(v.z >= 0) return (v2){intersection_point.x, intersection_point.y};

	return (v2){
		intersection_point.x = float_sign(intersection_point.x)*(1 - float_abs(intersection_point.x)),
		intersection_point.y = float_sign(intersection_point.y)*(1 - float_abs(intersection_point.y))
	};
}

v3 octahedral_decode(v2 v) {
	/*
	  o Decoding octahedral encoded vector

	    Again, we will operate within first quadrant by using absolute values and adjust for other
		quadrants after that.
		
	    First we need to check whether given vector is within or outside the octahedral diamond in the map.
		We can do this by recognizing that we are inside triangle (first quadrant) as long as x + y < 1.

		Another thing that we need is a way to reconstruct z coordinate. For this, we use intersection of
		ray and plane, just like for encoding, with an adjustment to the ray form. Let's say we have point
		(x, y) within triangle and that it is result of projection of upper sides of octahedron. In order
		to reconstruct z, we can find intersection point of plane representing octahedron side and ray
		that starts at point (x, y) and goes upward or downward. Two equation are:
		    (P - R)*n = 0
			P = p + d*t, where p is mentioned (x, y)

		Whether we look at upward or downward ray depends on whether we are inside or outside triangle.
		If the point is inside, we know that we got it by projecting intersection points from upper side
		to xy plane. In this case, we will send ray upward.
		If the point is outside, then we know that it belonged to lower side of octahedron, so we will
		point ray downward. In this case, we also need to map outside point to the inside of the triangle,
		which will just reverse computation done during encoding. After doing this, we now have ray origin
		and direction.
			
		We can find t because we know everything else.
		Solution for upward ray is   : t = 1 - (x + y)
		Solution for downward ray is : t = (x + y) - 1
		
		Another subtle thing to notice is how we can calculate t for the outside case. As mentioned, in
		order to calculate this t, we need to map point to the inside of the triangle and then direct ray
		downward. If we imagine upper octahedron side (or rather that plane) extending further through
		xy plane, we can notice that if we take outside point (x, y) and attach upward ray to it, we
		will get negative t that corresponds to that ray intersecting this extended side. This negative t
		actually gives us z value for decoded point when point is on the outside. This means that we can
		just calculate one t from input point, and use it for decoding in both cases, without having to
		look for different t (based on mapping of outside point to inside point) for outside case.
	*/
	float temp = float_abs(v.x) + float_abs(v.y);
	if(temp > 1) {
		return v3_unit((v3){
				float_sign(v.x)*(1 - float_abs(v.x)),
				float_sign(v.y)*(1 - float_abs(v.y)),
				1 - temp
			});
	}
	
	return v3_unit((v3){v.x, v.y, 1 - temp});
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

v3 random_unit_vector() {
	v3 v = {rand_float(), rand_float(), rand_float()};
	return v3_unit(v);
}

void example1() {
	/*
	  This example just shows how encoding is done.
	*/
	v3 u = (v3){12, -21, -8};
	u = v3_unit(u);
	v2 v = octahedral_encode(u);
	printf("ORIGINAL : %f, %f %f\n", u.x, u.y, u.z);
	printf("ENCODED  : %f, %f\n", v.x, v.y);
	v2 uv = octahedral_map_to_uv(v);
	printf("UV       : %f, %f\n", uv.x, uv.y);
	v3 w = octahedral_decode(v);
	printf("DECODED  : %f, %f %f\n", w.x, w.y, w.z);
}

void example2(int n) {
	/*
	  This example assigns color to each unit vector and then encodes that vector to
	  uv coordinates, preserving the color value. The result of this example is a
	  PPM image that shows encoded vectors (shows their colors).
	  Picture contains more colors for larger number of mapped vectors n.
	  For n=10000000, the whole image if filled with colors.
	  Firstly, we can see the central diamond that represents upper portion of octahedron
	  mapped to underlying plane (these colors correspond to direction vectors in the
	  upper portion of 3D space).
	  Secondly, we see outer 4 triangles that correspond to lower portion of octahedron
	  in a way described earlier. To see that it is a previously described mapping,
	  notice, for example, the lower right triangle. Now imagine a line that goes from
	  center of the picture (diamond) to the lower right corner. Rotate lower right
	  triangle around this line by 180 deg. After that, colors that are joining this
	  triangle and a diamond will match, and there will be no visible edge (green will
	  connect with green and pink with pink). Same conclusion can be drawn for other
	  triangles and their rotation around respective lines.
	  The fact that there is a visible edge and the colors do not match is the result of
	  our choice of mapping that allowed us to avoid temporary save (described earlier).

	  Additionally, try to imagine how colored octahedron looks just based on an image.
	  Central diamond is a coloring for all upper sides of octahedron.
	  Outside triangles, after mentioned rotation, are lower sides of octahedron.
	*/
#define MATRIX_DIM 400
	FILE *file = fopen("coded.ppm", "wb");
	if(file == NULL) return;

	fprintf(file, "P6\n%d %d\n255\n", MATRIX_DIM, MATRIX_DIM);
	if(ferror(file)) {
		fclose(file);
		return;
	}

	unsigned int matrix[MATRIX_DIM*MATRIX_DIM] = {};

	v3 v;
	v2 encoded;
	int color = 0;
	for(int i = 0; i < n; ++i) {
		v = random_unit_vector();
		encoded = octahedral_map_to_uv(octahedral_encode(v));
		color = (unsigned char)(floor((v.z + 1)*0.5*255));
		color <<= 8;
		color |= (unsigned char)(floor((v.y + 1)*0.5*255));
		color <<= 8;
		color |= (unsigned char)(floor((v.x + 1)*0.5*255));
		matrix[(int)(encoded.y*MATRIX_DIM)*MATRIX_DIM + (int)(encoded.x*MATRIX_DIM)] = color;
	}
	
	for(int i = 0; i < MATRIX_DIM*MATRIX_DIM; ++i){
		fwrite(&matrix[i], sizeof(char), 3, file);
		if(ferror(file)) return;
	}
	
	fclose(file);
#undef MATRIX_DIM
}

int main(void) {
	//example1();
	example2(10000000);
	
	return 0;
}
