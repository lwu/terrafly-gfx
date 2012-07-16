#ifndef GLWIN_H
#define GLWIN_H

#include <gfx/gui.h>

typedef void (APIENTRY *PFNGLOCKARRAYSEXTPROC)(GLint first, GLsizei count);
typedef void (APIENTRY *PFNGLUNLOCKARRAYSEXTPROC)(void);

extern PFNGLOCKARRAYSEXTPROC    glLockArraysEXT;
extern PFNGLUNLOCKARRAYSEXTPROC glUnlockArraysEXT;

void setup_gl_extensions();

#endif