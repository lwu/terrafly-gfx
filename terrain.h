#ifndef TERRAIN_H
#define TERRAIN_H

#include <gfx/raster.h>
#include <gfx/geom3d.h>
#include "radiusImg.h"
#include <memory.h>

#include <queue>

class Plane
{
public:
	Plane() {
		plane_d = 0.0f;
	}

	void init(const Vec3f& normal, const Vec3f& p) {
		unit_normal = normal / norm(normal);
		plane_d = normal * p;
	}

	float proj_distance(const Vec3f& p) {
		return (p*unit_normal) - plane_d;
	}

private:
	Vec3f  unit_normal;
	float plane_d;
};

enum vis_type {
	VISIBLE_ALL, VISIBLE_MAYBE, VISIBLE_NONE
};

struct Vec2i {
	Vec2i() : x(0), y(0) {}
	Vec2i(int x_0, int y_0) : x(x_0), y(y_0) {}

	int x;
	int y;
};

struct VertexError {
	VertexError() {
	}
	VertexError(const Vec2i& v0, const Vec2i& v1, float v_error, int lvl, vis_type v)
		: parent(v0), child(v1), error(v_error), level(lvl), vis(v) {}
	Vec2i parent;
	Vec2i child;
	float error;
	short level;
	vis_type vis;
};

#if 0

template <class T, int max_size>
class Queue
{
public:
	Queue() {}

	void push(const T& elem) {
		if (index > max_size) {
			return;
		}
		data[index++] = elem;
	}

	VertexError& front() {
		return data[index-1];
	}

	VertexError& top() {
		return data[index-1];
	}

	void pop() {
		index--;
	}

	bool empty() const {
		return index <= 0;
	}

	int size() const {
		return index;
	}
private:
	T data[max_size];
	int index;
};

#endif // 0

struct IndexArray {
	static const unsigned int INVALID;	

	IndexArray() {
		data = 0;		
	}

	~IndexArray() {
		delete[] data;
	}

	void init(int x, int y) {
		dim_x = x;
		dim_y = y;		
		data = new unsigned int[dim_x*dim_y];
		clear();
	}

	void clear();
	
	unsigned int& operator()(int x, int y) {
		return data[dim_x*y + x];
	}
private:
	int dim_x;
	int dim_y;	
	unsigned int* data;	
};

struct TriRange
{
	TriRange() {
		begin = end = tri_limit = 0;
	}

	TriRange(int rbegin, int rend, int num_tris) {
		begin = rbegin;
		end = rend;
		tri_limit = num_tris;
	}

	int begin;
	int end;
	int tri_limit;

	static const TriRange ALL;
};

bool operator==(const TriRange& t1, const TriRange& t2);

class Terrain
{
public:
    // The actual image data that defines the terrain.  Note that the Z image
    // is grayscale and the texture image is RGB.
    ByteRaster *z;
    ByteRaster *texture;

	// Pre-calculated error and distance (radius) bounds
	ByteRaster *errorImg;
	//ByteRaster *GradiusImg;
	RadiusImg radiusImg;
	
    // Horizontal distance (in meters) between pixels
    int pixel_spacing;

    // Vertical distance (in meters) corresponding to a unit change in pixel
    // value (on the 0..255 scale) and the real-world height corresponding to
    // graylevel 0.
    float height_unit;
    float height_base;

    // Artificial stretch factor to accentuate terrain features.
    float height_stretch;
	float inv_w, inv_h;
	float height_factor;

	// Subdivision
	int subdiv_level;
	int max_subdiv_level;
	bool always_subdivide;

	// Culling
	bool do_view_frustum_culling;

	// Error control
	int tau;
	int lambda;
	int tri_limit;

	// Drawing and triangle limit parameters
	bool draw_style_fast;
	int total_tri_count;
	int tri_count;
	int double_cnt;
	int single_cnt;
public:
    Terrain() { 
		z = texture = NULL; 

		subdiv_level = 1;
		max_subdiv_level = 1;
		always_subdivide = false;
		do_view_frustum_culling = true;

		set_error_tolerance(320 /*lambda*/, 1/*tau*/);	

		tri_limit = 32000;
	}
    ~Terrain() { if(z) delete z; }

    Vec3f pixel_to_space(int i, int j);
    void emit_point(int i, int j);

	void draw();

	void draw_fast(const TriRange& range);
	void draw_correct(const TriRange& range);

	void set_eye(const Vec3f& pos, const Vec3f& dir, const Vec3f& up,
		float view_angle, float aspect, float dtheta, float phi);

	void set_error_tolerance(int l, int t) {
		tau = t;
		lambda = l;
		one_over_kappa = (float)lambda / (float)tau;
	}

	void init() {
		vertex_indices.init(z->width(), z->height());
	}

private:
	void draw_triangle(const Vec2i& v0, const Vec2i& v1, int lvl, vis_type vis);
	void render_triangle(const Vec2i& v0, const Vec2i& v1);
	vis_type calc_visibility(const Vec3f& v, float radius);

	float calc_vertex_error(const Vec2i& v1, vis_type vis, vis_type* new_vis);


	Vec3f terrain_viewpoint_pos;
	Vec3f terrain_viewpoint_dir;
	float one_over_kappa;

	// View frustum
	Plane view_frustum_planes[4];

	std::priority_queue<VertexError> errors;
	//std::stack<VertexError> errors;
	//std::queue<VertexError> errors;
	//Queue<VertexError, 65536> errors;

	void render_indices();

	IndexArray vertex_indices;
};

#endif