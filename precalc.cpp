
#include <FL/fl_ask.H>
#include <assert.h>
#include "precalc.h"

// VC6 hack
#define min std::_cpp_min
#define max std::_cpp_max

#include <set>
#include <vector>
#include <algorithm>

using namespace std;

typedef Raster<short> ShortRaster; 

void error_popup(const char* err_msg);

template <class T>
void fill_raster(Raster<T>& img, T val)
{
	T* head = img.head();
	int len = img.length();
	for (int i=0; i < len; i++) {
		*head++ = val;
	}
}

void write_tiff_shortraster(const char* fname, ShortRaster& radiusImg)
{
	int data_length = radiusImg.length();
	unsigned char* img = new unsigned char[data_length];

	short* p = radiusImg.head();	
	short max_val = 1;

	int i;

	// Find set of radii
	set<int> radii;
	p = radiusImg.head();
	for (i=0; i < data_length; i++) {		
		radii.insert(*p++);
	}

	cout << "[Lookup table: ";
	
	// Copy to byte vector
	typedef vector<int> int_vec;
	int_vec vradii;
	set<int>::iterator sit = radii.begin();
	for ( ; sit != radii.end(); ++sit) {
		static bool first = false;
		if (first) cout << ", ";

		cout << *sit;		
		vradii.push_back(*sit);

		first = true;
	}
	cout << "] ";

	int_vec::iterator vbegin = vradii.begin();
	int_vec::iterator vend = vradii.end();
	int_vec::iterator it;

	// Store lookup table indices
	p = radiusImg.head();
	unsigned char* cp = img;
	for (i=0; i < data_length; i++) {
		int val = *p++;
		it = lower_bound(vbegin, vend, val);
		assert(it != vend);

		int index = it - vbegin;
		*cp++ = index;
		assert(vradii[index] == val);
	}

	FILE *f = fopen("radii.dat", "wb");
	if (!f) {
		error_popup("Couldn't open radii datafile for writing");
	}

	// Pad to 256
	while (vradii.size() < 256) {
		vradii.push_back(0);
	}	

	// Write to disk: lookup table first, then radii indices
	int rv;
	rv = fwrite((void*)&vradii[0], sizeof(int), 256, f);
	assert(rv == 256);
	rv = fwrite(img, sizeof(unsigned char), data_length, f);	
	assert(rv == data_length);

	fclose(f);
	
	delete[] img;
}

// Sentinel value signifying yet-unknown values
#define UNDEFINED_VAL	0xFE

// Returns true iff a point is outside the valid area of an image
bool outside_image(const ByteRaster& img, int pt_x, int pt_y)
{
	return (pt_x < 0 || pt_x >= img.width() ||
			pt_y < 0 || pt_y >= img.height());
}

void check_terrain_size(int w, int h)
{
	if ((w & 0x1) && (h & 0x1)) {
		int new_w = w - 3;
		int new_h = h - 3;

		if (!(new_w & w) && !(new_h & h)) {
			return;
		}
	}

	fl_alert("I can't handle images that aren't a power of two, plus one");
}

// Calculate object space error by finding the difference between the
// real z-midpoint value and the average of two of the triangle endpoints
int calc_obj_space_error(const ByteRaster& terrain, ByteRaster& errorImg,
						 int pt_x,int dx, int pt_y,int dy)
{
	int x = pt_x+dx;
	int y = pt_y+dy;

	std::swap(dx, dy);
	dy = -dy;	
	
	if (outside_image(terrain, x+dx, y+dy) || (outside_image(terrain, x-dx, y-dy))) {
		return 0;
	}
	
	int v1_z = *terrain.pixel(x+dx,y+dy);
	int v2_z = *terrain.pixel(x-dx,y-dy);	
	
	return abs(*terrain.pixel(x,y) - (v1_z + v2_z)/2);
}

struct vdata {
	vdata() {}
	vdata(int e, int r) : error(e), radius(r) {}
	int error;
	int radius;
};


vdata calc_error(const ByteRaster& terrain, ByteRaster& errorImg, ShortRaster& radiusImg,
				 int pt_x,int dx, int pt_y,int dy, int edge_len, bool split_normal)
{
	int x = pt_x + dx;
	int y = pt_y + dy;
	if (outside_image(errorImg, x, y)) {
		return vdata(-1, -1);
	}

	int errorVal  = *errorImg.pixel(x,y);
	int radiusVal = *radiusImg.pixel(x,y);
	if (errorVal != UNDEFINED_VAL) {
		// already computed
		return vdata(errorVal, radiusVal);
	} else {
		// need to compute
		vdata c0, c1, c2, c3;
		int new_edge_len = edge_len;
		int ospace_error = calc_obj_space_error(terrain, errorImg, 
			pt_x,dx, pt_y,dy);
		short radius = 0;
	
		if (split_normal) {
			new_edge_len /= 2;
		}

		if (edge_len >= 1) {
			// recurse		
			if (split_normal) {
				// split n, e, s, w
				split_normal = !split_normal;
				c0 = calc_error(terrain, errorImg, radiusImg, x,0        , y,-edge_len, new_edge_len, split_normal);
				c1 = calc_error(terrain, errorImg, radiusImg, x,+edge_len, y,0,         new_edge_len, split_normal);
				c2 = calc_error(terrain, errorImg, radiusImg, x,0        , y,+edge_len, new_edge_len, split_normal);
				c3 = calc_error(terrain, errorImg, radiusImg, x,-edge_len, y,0,         new_edge_len, split_normal);

				radius = edge_len;
			} else {
				// split ne, se, sw, nw
				split_normal = !split_normal;
				c0 = calc_error(terrain, errorImg, radiusImg, x,+edge_len, y,-edge_len, new_edge_len, split_normal);
				c1 = calc_error(terrain, errorImg, radiusImg, x,+edge_len, y,+edge_len, new_edge_len, split_normal);
				c2 = calc_error(terrain, errorImg, radiusImg, x,-edge_len, y,+edge_len, new_edge_len, split_normal);
				c3 = calc_error(terrain, errorImg, radiusImg, x,-edge_len, y,-edge_len, new_edge_len, split_normal);				

				if (edge_len == 1) {
					radius = 2;
				} else {
					radius = (float)edge_len * sqrt(2) + 0.5f;
				}
			}

			int max_child_error = max( max(c0.error,c1.error), max(c2.error,c3.error) );
			ospace_error = max(ospace_error, max_child_error);

			radius += max( max(c0.radius,c1.radius), max(c2.radius,c3.radius) );
		}

		assert(ospace_error >= 0 && ospace_error < 255);
		assert(radius >= 0 && radius < 65535);

		// memoize calculated error and radius
		*errorImg.pixel(x,y)  = ospace_error;
		*radiusImg.pixel(x,y) = radius;

		return vdata(ospace_error, radius);
	}
}

void precalc_terrain(const ByteRaster& terrain)
{
	cout << "Pre-calculating error bounds for " << terrain.width() 
		 << "x" << terrain.height() << " terrain... ";

	// Per vertex, calculates delta^star_i -- upper bound
	// on object space error. See Lindstrom and Pascucci's
	// "Visualization of Large Terrains Made Easy" for details

	int nchan = 1;
	int w = terrain.width();
	int h = terrain.height();

	check_terrain_size(w,h);

	ByteRaster errorImg(w, h, nchan);
	ShortRaster radiusImg(w, h, nchan);

	fill_raster<unsigned char>(errorImg, UNDEFINED_VAL);
	fill_raster<short>(radiusImg, UNDEFINED_VAL);

	if (w >= h) {
		int aspect_ratio = (w-1) / (h-1);
		int edge_len = h / 2;

		int num_squares_x = aspect_ratio;

		int center_x = (w-1) / num_squares_x / 2;
		int center_y = (h-1) / 2;
		const bool split_normal = true;
		int dx = (w-1) / aspect_ratio;
		int dy = (h-1);
		int x  = 0;

		for (int i=0; i < num_squares_x; i++) {
			*errorImg.pixel(x,0)     = 0;
			*errorImg.pixel(x+dx,0)  = 0;
			*errorImg.pixel(x,dy)    = 0;
			*errorImg.pixel(x+dx,dy) = 0;
			
			*radiusImg.pixel(x,0)     = 0;
			*radiusImg.pixel(x+dx,0)  = 0;
			*radiusImg.pixel(x,dy)    = 0;
			*radiusImg.pixel(x+dx,dy) = 0;			

			vdata err0 = calc_error(terrain, errorImg, radiusImg,
				x,dx/2, 0,dy/2, edge_len, split_normal);
			assert(err0.error >= 0);

			x += dx;
		}
	} else {
		fl_alert("I can't handle images where width < height");
	}	

	bool rv = write_tiff_image("errorImg.tif", errorImg);

	cout << "done" << endl;

	cout << "Serializing compressed radius data and lookup table... ";
	write_tiff_shortraster("radiusImg.tif", radiusImg);	
	cout << "done" << endl;
}


