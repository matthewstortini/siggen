#ifndef Siggen_Utils
#define Siggen_Utils

namespace Siggen
{

//global variable for verbosity_level
enum verbosity_level {TERSE, NORMAL, CHATTY};
extern verbosity_level verbosity_setting;

#define TELL_NORMAL if (verbosity_setting >= NORMAL) tell
#define TELL_CHATTY if (verbosity_setting >= CHATTY) tell
void tell(const char *format, ...);
void error(const char *format, ...);

//matrix-mult for electron velocity model
void mat_mult( double a[3][3], double b[3][3], double out[3][3]);
void dot(double a[3][3], double b[3], double out[3]);
double dot(double a[3], double b[3]);
void transpose(double a[3][3], double out[3][3]);

#define SQ(x) ((x)*(x))

/* enum to identify cylindrical or cartesian coords */
#define CYL 0
#define CART 1

#define MAX_LINE 512


}

#endif
