/*!************************************************************************
 *
 * \file callPy.hxx
 * \brief File included into the physics code
 *
 * Pycall
 *
 **************************************************************************
 * Dmitry Meyerson added to use with BOUT++
 *
 **************************************************************************/

#ifndef __PYCALL_H__
#define __PYCALL_H__

#include <python2.6/Python.h>

//post processing with python
int callPy(int argc, char **argv);

#endif // __PYCALL_H__
