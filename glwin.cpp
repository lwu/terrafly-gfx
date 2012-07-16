#include "glwin.h"
#include "error.h"

PFNGLOCKARRAYSEXTPROC    glLockArraysEXT   = NULL;
PFNGLUNLOCKARRAYSEXTPROC glUnlockArraysEXT = NULL;


void setup_gl_extensions()
{	
	const unsigned char* extensions = glGetString(GL_EXTENSIONS);

	//cout << extensions << endl;
	bool cva_support = ( strstr((const char*)extensions, "GL_EXT_compiled_vertex_array") ? true : false);

	if (!cva_support) {
		error_popup("Compiled vertex array extension not supported");
	}

	glLockArraysEXT = (PFNGLOCKARRAYSEXTPROC)wglGetProcAddress("glLockArraysEXT");
	glUnlockArraysEXT = (PFNGLUNLOCKARRAYSEXTPROC)wglGetProcAddress("glUnlockArraysEXT");

}
