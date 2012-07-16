
#ifndef RADIUS_IMG_H
#define RADIUS_IMG_H

// Radius data is stored in a compressed format.
//
// The stored indices point to a 256-element lookup table

struct RadiusImg
{
public:
	RadiusImg() {
		data = 0;
	}

	~RadiusImg() {
		delete[] data;
	}

	void init(const char* fname, int w, int h);

	int operator()(int i, int j) {		
		return lookup_table[*(data + j*W + i)];
	}

	int lookup_table[256];
	unsigned char* data;

	int W, H;
};

#endif
