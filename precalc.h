
#ifndef PRECALC_H
#define PRECALC_H

#include <gfx/raster.h>

/*typedef unsigned char byte;
struct radii_header {
	byte magic_cookie;
	byte lookup_table[256];
	byte magic_constant;
	byte radii[0];
};*/

void precalc_terrain(const ByteRaster& terrainImg);

#endif // PRECALC_H