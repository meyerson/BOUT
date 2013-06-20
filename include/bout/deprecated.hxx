/*
 * Defines a DEPRECATED macro, to mark functions for future removal
 */

#ifndef __DEPRECATED_H__
#define __DEPRECATED_H__

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED(func) __declspec(deprecated) func
#else
//#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED(func) func
#endif

#endif // __DEPRECATED_H__

