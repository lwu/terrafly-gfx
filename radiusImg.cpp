
#include <stdio.h>
#include <assert.h>
#include "radiusImg.h"
#include "error.h"

void RadiusImg::init(const char* fname, int w, int h)
{
	W = w;
	H = h;

	FILE *f = fopen(fname, "rb");
	if (!f) {
		error_popup("Couldn't read radii data file");
	}

	int rv;
	rv = fread( (void*)this->lookup_table, sizeof(int), 256, f);	
	assert(rv == 256);

	int num_bytes = w*h;
	this->data = new unsigned char[num_bytes];
	rv = fread( this->data, sizeof(unsigned char), num_bytes, f);
	assert(rv == num_bytes);

	fclose(f);
}
