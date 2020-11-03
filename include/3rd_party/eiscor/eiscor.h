/* define STDERR as 0 if undefined */
#ifndef STDERR 
#define STDERR (0)
#endif 

/* if DEBUG is defined, redefine as .TRUE. otherwise define as .FALSE. */
#ifdef DEBUG 
#undef DEBUG
#define DEBUG (.TRUE.)
#else 
#define DEBUG (.FALSE.)
#endif   

/* if VERBOSE is defined, redefine as .TRUE. otherwise define as .FALSE. */
#ifdef VERBOSE 
#undef VERBOSE
#define VERBOSE (.TRUE.)
#else 
#define VERBOSE (.FALSE.)
#endif 

/* define EISCOR_DBL_EPS as unit roundoff */
#ifdef EISCOR_DBL_EPS 
#undef EISCOR_DBL_EPS
#define EISCOR_DBL_EPS (epsilon(1d0))
#else 
#define EISCOR_DBL_EPS (epsilon(1d0))
#endif

/* define EISCOR_DBL_INF as largest supported positive double precision number */
#ifdef EISCOR_DBL_INF 
#undef EISCOR_DBL_INF
#define EISCOR_DBL_INF (huge(1d0))
#else 
#define EISCOR_DBL_INF (huge(1d0))
#endif

/* define EISCOR_DBL_PI as nearest double precision number to pi */
#ifdef EISCOR_DBL_PI 
#undef EISCOR_DBL_PI
#define EISCOR_DBL_PI (3.141592653589793239d0)
#else 
#define EISCOR_DBL_PI (3.141592653589793239d0)
#endif

