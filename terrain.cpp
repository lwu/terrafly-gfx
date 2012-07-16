#include <gfx/gui.h>
#include "terrain.h"
#include "globals.h"
#include "glwin.h"
#include <assert.h>

using namespace std;

extern bool use_vertex_array;
extern bool rot_compensate;
extern bool phi_compensate;

Vec3f rotate_direction(const Vec3f& moving, float angle, const Vec3f& axis);

struct VertexNode {
	VertexNode() {}
	VertexNode(float a, float b, float r, float s, float t) :
		u(a), v(b), x(r), y(s), z(t) {}

	float u, v;
	float x, y, z;
};

bool operator<(const VertexError& v1, const VertexError& v2)
{
	return v1.error < v2.error;
}


const unsigned int IndexArray::INVALID = -1;

#define MAX_TRIANGLES	(65536*5)
#define MAX_INDICES		(MAX_TRIANGLES*3)

struct geo_data {
	geo_data() {
		active_index = 0;
		num_indices  = 0;
	}

	VertexNode active_verts[MAX_TRIANGLES];
	int active_index;

	unsigned int indices[MAX_INDICES];
	int num_indices;

	unsigned short index_pairs[MAX_INDICES*2];
};

geo_data geometry[2];
geo_data* geo_updating = &geometry[0];
geo_data* geo_drawing = &geometry[1];

int updating_frame = 0;
const int max_update_frame = 3;

// Converts data for pixel (i,j) into world space (x,y,z) coordinates.
Vec3f Terrain::pixel_to_space(int i, int j)
{   
    return Vec3f(i * pixel_spacing,
		j * pixel_spacing,
		z->pixel(i,j)[0] * height_factor);		
}

void IndexArray::clear()
{	
	int n = dim_x*dim_y;
	for (int i=0; i < n; i++) {
		data[i] = -1;
	}
}

void Terrain::render_indices()
{
	int n = sizeof(geo_data);
	
	if (updating_frame == 0) {
		glInterleavedArrays(GL_T2F_V3F, 0, geo_drawing->active_verts);
		
		glLockArraysEXT(0, geo_drawing->active_index);
	}

	glDrawElements(GL_TRIANGLES, geo_drawing->num_indices, 
		GL_UNSIGNED_INT, geo_drawing->indices);

	if (updating_frame == max_update_frame) {
		glUnlockArraysEXT();
	}
}


void Terrain::emit_point(int i, int j)
{	
	if (use_vertex_array) {
		if (geo_updating->num_indices >= MAX_INDICES) {
			return;
		}
				
		unsigned int& new_index = vertex_indices(i,j);
		geo_updating->index_pairs[geo_updating->num_indices*2] = i;
		geo_updating->index_pairs[geo_updating->num_indices*2+1] = j;
		
		if (new_index == IndexArray::INVALID) {
			new_index = geo_updating->active_index;
			
			float u = i * inv_w;
			float v = j * inv_h;
			Vec3f p = pixel_to_space(i,j);
			
			geo_updating->active_verts[geo_updating->active_index++] = 
				VertexNode(u, v, p[0], p[1], p[2]);
		} 
		
		geo_updating->indices[geo_updating->num_indices++] = new_index;

	} else {		
		float u = i * inv_w;
		float v = j * inv_h;
		Vec3f p = pixel_to_space(i,j);
		
		glTexCoord2f(u, v);
		glVertex3fv(p);	
	}
}

void Terrain::set_eye(const Vec3f& pos, const Vec3f& dir, const Vec3f& up,
					  float y_view_angle, float aspect, float dtheta, float phi)
{
	if (!use_vertex_array || updating_frame == 0) {
		Vec3f left = cross(up,dir);	
		terrain_viewpoint_pos = pos;
		terrain_viewpoint_dir = dir;	 
		
		float left_extra = 1.0f;
		float right_extra = 1.0f;
		
		if (use_vertex_array && rot_compensate) {
			// increase view angle to compensate for buffering,
			//  proportional to viewer rotation (in one direction)
			if (dtheta < 0.0f) {
				dtheta = min(fabs(dtheta), 4.0f);
				right_extra = 1.0f + dtheta * 0.25f;
			} else if (dtheta > 0.0f) {
				dtheta = min(fabs(dtheta), 4.0f);
				left_extra = 1.0f + dtheta * 0.25f; 
			}
		}

		if (phi_compensate) {
			// increase view angle if user looks down 
			//  (which makes view frustum culling look ugly)
			if (phi > 0.0f) {
				 float dphi = phi / 90.0f;
				 left_extra *= 1.0f + dphi;
				 right_extra *= 1.0f + dphi;
			}
		}
		
		y_view_angle *= 0.5f;
		float x_view_angle = y_view_angle * aspect;
		
		Vec3f left_frustum   =  cross( rotate_direction(terrain_viewpoint_dir,  
			x_view_angle*left_extra, up), up );
		Vec3f right_frustum  = -cross( rotate_direction(terrain_viewpoint_dir, 
			-x_view_angle*right_extra, up), up );
		
		//Top + bottom frustum doesn't work properly with Lindstrom's algorithm	
		//Vec3f bottom_frustum =  cross( rotate_direction(dir, +y_view_angle, left), left );
		//Vec3f top_frustum    = -cross( rotate_direction(dir, -y_view_angle, left), left );	
		
		view_frustum_planes[0].init(left_frustum, pos);
		view_frustum_planes[1].init(right_frustum, pos);
		//view_frustum_planes[2].init(bottom_frustum, pos);
		//view_frustum_planes[3].init(top_frustum, pos);	
	}
}

vis_type Terrain::calc_visibility(const Vec3f& p, float radius)
{	
	bool visible_all = true;

	for (int i=0; i < 2; i++) {
		Plane& plane = view_frustum_planes[i];

		float dist = plane.proj_distance(p);
		if (dist > radius) {
			// totally visible
		} else if (dist < -radius) {
			// totally invisible
			return VISIBLE_NONE;
		} else {
			// maybe visible
			visible_all = false;
		}
	}

	return visible_all ? VISIBLE_ALL : VISIBLE_MAYBE;
}

float Terrain::calc_vertex_error(const Vec2i& v1, vis_type vis, vis_type *new_vis)
{
	Vec3f p = pixel_to_space(v1.x, v1.y);
	float radius = radiusImg(v1.x, v1.y) * pixel_spacing;

	if (vis == VISIBLE_ALL) {
		// totally visible, no need to recalc for subtris
		*new_vis = VISIBLE_ALL;
	} else if (vis == VISIBLE_NONE) {
		// invisible; return
		*new_vis = VISIBLE_NONE;

		return -1.0f;
	} else { // vis == VISIBLE_MAYBE
		*new_vis = calc_visibility(p, radius);

		if (*new_vis == VISIBLE_NONE) {
			return -2.0f;
		}
	}

	Vec3f dp = p - terrain_viewpoint_pos;
	float rhs = norm2(dp);

	// See equation (3) in Lindstrom's "Terrain Visualization" paper
	float obj_space_error = (errorImg->pixel(v1.x, v1.y)[0]) * height_factor;	

	float lhs = one_over_kappa * obj_space_error + radius;	
	lhs *= lhs;	
	
	if (display_error_lines) {
		if (!use_vertex_array) {
			glEnd();
		}
		
		glColor3f(1.0f, 0.0f, 0.0f);
		glBegin(GL_LINES);
		glVertex3f(p[0], p[1], 0.0f);		
		glVertex3f(p[0], p[1], radius * 0.5f);
		glEnd();
				
		glColor3f(1.0f, 1.0f, 1.0f);
		
		if (!use_vertex_array) {
			glBegin(GL_TRIANGLES);
		}
	}

	return (lhs - rhs);
}

typedef pair<Vec2i, Vec2i> PairVec2i;

inline PairVec2i get_children(const Vec2i& v0, const Vec2i& v1)
{
	short dx = -(v1.y - v0.y);
	short dy = (v1.x - v0.x);
	
	return make_pair(
		Vec2i((v0.x+ v1.x-dx)/2, (v0.y+ v1.y-dy)/2),   // left  child
		Vec2i((v0.x+ v1.x+dx)/2, (v0.y+ v1.y+dy)/2) ); // right child		
}

void Terrain::draw_triangle(const Vec2i& v0, const Vec2i& v1, int lvl, vis_type vis)
{
	int dx = -(v1.y - v0.y);
	int dy = v1.x - v0.x;

	vis_type new_vis;	
	float error = calc_vertex_error(v1, vis, &new_vis);

	bool subdivide = (error >= 0.0f);
	
	if (subdivide) {
		if (lvl > 1) {
			PairVec2i kids = get_children(v0, v1);

			draw_triangle(v1, kids.first,  lvl-1, new_vis);
			draw_triangle(v1, kids.second, lvl-1, new_vis);
		} else {
			tri_count += 2;
			double_cnt++;
			emit_point(v0.x, v0.y);
			emit_point(v1.x-dx, v1.y-dy);
			emit_point(v1.x, v1.y);

			emit_point(v0.x, v0.y);
			emit_point(v1.x, v1.y);
			emit_point(v1.x+dx, v1.y+dy);
		}
	} else {
		tri_count++;
		single_cnt++;
		emit_point(v0.x, v0.y);
		emit_point(v1.x-dx, v1.y-dy);
		emit_point(v1.x+dx, v1.y+dy);
	}
}


void Terrain::render_triangle(const Vec2i& v0, const Vec2i& v1)
{
	short dx = -(v1.y - v0.y);
	short dy = v1.x - v0.x;

	tri_count++;
	single_cnt++;
	emit_point(v0.x, v0.y);
	emit_point(v1.x-dx, v1.y-dy);
	emit_point(v1.x+dx, v1.y+dy);
}

void Terrain::draw()
{	
	TriRange range;

	if (!use_vertex_array) {
		total_tri_count = 0;
		tri_count = 0;
		single_cnt = 0;
		double_cnt = 0;

		range = TriRange::ALL;
	} else { // vertex_array -> partial update
		int range_tri_limit;
		if (draw_style_fast) {
			range_tri_limit = total_tri_count;
		} else { // correct
			range_tri_limit = tri_limit;
		}

		range = TriRange(updating_frame, max_update_frame, range_tri_limit);
	}

	double t0 = get_cpu_time();

	if (draw_style_fast) {
		draw_fast(range);
	} else {
		draw_correct(range);
	}
	
	if (use_vertex_array) {
		render_indices();

		if (updating_frame == max_update_frame) {
			// invalidate used part of vertex index arrays
			for (int i=0; i < geo_updating->num_indices; i++) {
				int x = geo_updating->index_pairs[i*2];
				int y = geo_updating->index_pairs[i*2+1];
				vertex_indices(x,y) = IndexArray::INVALID;
			}

			// swap geometry array pointers
			swap(geo_updating, geo_drawing);

			// reset vars
			geo_updating->num_indices = 0;
			geo_updating->active_index = 0;
			
			total_tri_count = tri_count;
			tri_count = 0;
			single_cnt = 0;
			double_cnt = 0;
			
			updating_frame = 0;
		} else {			
			updating_frame++;
		}
	} else {		
		total_tri_count = tri_count;
	}
}

const TriRange TriRange::ALL(0,0,0);

bool operator==(const TriRange& t1, const TriRange& t2) {
	return t1.begin == t2.begin && t1.end == t2.end && t1.tri_limit == t2.tri_limit;
}


void Terrain::draw_fast(const TriRange& range)
{
	int h = z->height();
	int w = z->width();

	int aspect_ratio = (w-1) / (h-1);
	int num_squares_x = aspect_ratio;

	int dx = (w-1) / num_squares_x;
	int center_x = dx / 2;
	int center_y = (h-1) / 2;	

	int lvl = subdiv_level;

	Vec2i v0(0,0);
	Vec2i v1(center_x, center_y);
	Vec2i v2(center_x*2, center_y*2);	

	vis_type vis;
	if (do_view_frustum_culling) {
		vis = VISIBLE_MAYBE;
	} else {
		vis = VISIBLE_ALL;
	}	
	
	if (range == TriRange::ALL) {
		glBegin(GL_TRIANGLES);	
		// draw as much as possible
		for (int i=0; i < num_squares_x; i++) {
			draw_triangle(v0, v1, lvl, vis);
			draw_triangle(v2, v1, lvl, vis);
			
			v0.x += dx;
			v1.x += dx;
			v2.x += dx;
		}
		glEnd();
	} else {
		// do a partial update, based on triangle range parameter
		int num_tri = 0;
		
		for (int i=0; i < num_squares_x; i++) {
			if (num_tri == range.begin || num_tri > range.end) {
				draw_triangle(v0, v1, lvl, vis); 
			}
			num_tri++;
			if (num_tri == range.begin || num_tri > range.end) {
				draw_triangle(v2, v1, lvl, vis); 
			}
			num_tri++;
						
			v0.x += dx;
			v1.x += dx;
			v2.x += dx;
		}
	}		
}

template <class Container>
void clear(Container& c) {
	// Faster than repeated c.pop() while not empty()?
	c = Container();
}

// Guarantee triangle limit
void Terrain::draw_correct(const TriRange& range)
{
	int h = z->height();
	int w = z->width();

	int aspect_ratio = (w-1) / (h-1);
	int num_squares_x = aspect_ratio;

	int dx = (w-1) / num_squares_x;
	int center_x = dx / 2;
	int center_y = (h-1) / 2;	

	int lvl = subdiv_level;
	
	Vec2i v0(0,0);
	Vec2i v1(center_x, center_y);
	Vec2i v2(center_x*2, center_y*2);	

	vis_type vis;
	if (do_view_frustum_culling) {
		vis = VISIBLE_MAYBE;
	} else {
		vis = VISIBLE_ALL;
	}
	
	glBegin(GL_TRIANGLES);
	
	bool partial_update = !(range == TriRange::ALL);
	int partial_tri_limit = 0;
	
	if ( !partial_update ) {
		if( !errors.empty() ) {
			clear(errors);
		}
	} 
	
	if ( !partial_update || updating_frame == 0 ) {
		// Add initial top-level vertex errors into priority queue		
		for (int i=0; i < num_squares_x; i++) {
			vis_type new_vis;
			float error = calc_vertex_error(v1, vis, &new_vis);
			
			errors.push( VertexError(v0, v1, error, lvl, new_vis) );
			errors.push( VertexError(v2, v1, error, lvl, new_vis) );
			
			v0.x += dx;
			v1.x += dx;
			v2.x += dx;
		}
	} 
	
	if (partial_update) {
		if (range.begin < range.end)  {
			partial_tri_limit = (float)(range.begin+1) / (float)(range.end+1) * range.tri_limit;
		} else {
			// might need to draw a few more triangles to take care of
			//  the corner case when we draw a few triangles after last
			partial_update = false;
		}
	}
	
	
	bool limit_hit = (tri_count + errors.size() >= tri_limit);
	
	glColor3f(1.0f, 1.0f, 1.0f);
	
	// While still errors left, subdivide if triangle quota not met,
	//  unless we're doing a partial update
	while ( !errors.empty() && (!partial_update || tri_count <= partial_tri_limit) ) {
		VertexError v( errors.top() );
		errors.pop();			
		
		// check tri limit
		glColor3f(1.0f, 1.0f, 1.0f);
		if (!limit_hit && tri_count + errors.size() >= tri_limit) {
			limit_hit = true;				
			
			// draw the couple of triangles that happen to have
			//  the same error-valued vertex, to avoid occasional cracks
			typedef std::vector<VertexError> error_vec;
			error_vec limit;
			limit.push_back( v );
			while (v.error == errors.top().error) {
				assert( !errors.empty() );
				limit.push_back( errors.top() );
				errors.pop();
			}
			
			int size = limit.size();				
			
			//bool uva = use_vertex_array;			
			//use_vertex_array = false;
						
			//glColor3f(0.0f, 1.0f, 0.0f);						
			for (int i=0; i < size; i++) {
				VertexError& lim = limit[i];
				
				if (lim.error < 0.0f) {
					bool stop = true;
				}
				
				if (lim.level > 1) {
					PairVec2i kids = get_children(lim.parent, lim.child);
					render_triangle(lim.child, kids.first);
					render_triangle(lim.child, kids.second);
				} else {
					render_triangle(lim.parent, lim.child);
				}
			}
			//glColor3f(0.0f, 1.0f, 0.0f);
						
			assert( !errors.empty() );
			
			v = errors.top();
			errors.pop();
		}
		
		// draw normally
		if (v.level > 1 && !limit_hit) {
			int new_level = v.level - 1;
			
			PairVec2i kids = get_children(v.parent, v.child);
			
			vis_type new_vis_first;
			vis_type new_vis_second;
			float error_first  = calc_vertex_error(kids.first, v.vis, &new_vis_first);
			float error_second = calc_vertex_error(kids.second, v.vis, &new_vis_second);
			
			if (error_first >= 0.0f) {
				errors.push( VertexError(v.child, kids.first, error_first, new_level, new_vis_first) );
			} else {
				render_triangle(v.child, kids.first);
			}
			
			if (error_second >= 0.0f) {
				errors.push( VertexError(v.child, kids.second, error_second, new_level, new_vis_second) );
			} else {
				render_triangle(v.child, kids.second);
			}
		} else {
			render_triangle(v.parent, v.child);
		}
		
		
	}
	
	glEnd();
}
