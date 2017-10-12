#include "Utils.h"
#include <stdio.h>
#include <stdarg.h>

namespace Siggen
{
  verbosity_level verbosity_setting = NORMAL;

  /* tell
     write to stdout, provided that verb_level is above the threshold */
  void tell(const char *format, ...){
    va_list ap;

    va_start(ap, format);
    vprintf(format, ap);
    va_end(ap);
    return;
  }

  /*error
    report error messages to stderr */
  void error(const char *format, ...) {
    va_list ap;

    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    return;
  }

  void mat_mult(double a[3][3], double b[3][3], double out[3][3]){
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
        out[i][j] = 0;
        for(int k = 0; k < 3; ++k){
          out[i][j] += a[i][k] * b[k][j];
        }
      }
    }
  }

  void dot(double a[3][3], double b[3], double out[3]){

    for(int i = 0; i < 3; ++i){
      out[i] = 0;
      for(int k = 0; k < 3; ++k){
        out[i] += a[i][k] * b[k];
      }
    }
  }

  double dot(double a[3], double b[3]){
    double out = 0;

    for(int k = 0; k < 3; ++k){
      out += a[k] * b[k];
    }
    return out;
  }

  void transpose(double a[3][3], double out[3][3]){

    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
        out[i][j] = a[j][i];
      }
    }
  }



}
