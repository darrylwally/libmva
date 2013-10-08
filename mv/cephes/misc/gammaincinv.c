/*
Copyright (c) 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright (c) 2003-2012 SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of Enthought nor the names of the SciPy Developers
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

//#include <Python.h>
//#include <numpy/npy_math.h>

#include <stdio.h>
#include <math.h>

//#include "../cephes.h"

#include <cephes/cmath/mconf.h>
#undef fabs
#include "misc.h"

extern double igam(double, double);
extern double igami(double, double);

/* Limits after which to issue warnings about non-convergence */
#define ALLOWED_ATOL (1e-306)
#define ALLOWED_RTOL (1e-6)

//void scipy_special_raise_warning(char *fmt, ...);

/*
  Inverse of the (regularised) incomplete Gamma integral.

  Given a, find x such that igam(a, x) = y.
  For y not small, we just use igami(a, 1-y) (inverse of the complemented
  incomplete Gamma integral). For y small, however, 1-y is about 1, and we
  lose digits.

*/

extern double MACHEP, MAXNUM;
/* double precision NaN define */
static const long long _NaN = 0x7ff8000000000000;
#define GAMMAINCINV_NaN (*(double *)&_NaN)


static double
gammainc(double x, double params[2])
{
    return /*cephes_*/igam(params[0], x) - params[1];
}

double
gammaincinv(double a, double y)
{
    double lo = 0.0, hi;
    double flo = -y, fhi = 0.25 - y;
    double params[2];
    double best_x, best_f, errest;
    fsolve_result_t r;

    if (a <= 0.0 || y <= 0.0 || y >= 0.25) {
        return igami(a, 1-y);
    }

    /* Note: flo and fhi must have different signs (and be != 0),
     *       otherwise fsolve terminates with an error.
     */

    params[0] = a;
    params[1] = y;
    hi = igami(a, 0.75);
    /* I found Newton to be unreliable. Also, after we generate a small
       interval by bisection above, false position will do a large step
       from an interval of width ~1e-4 to ~1e-14 in one step (a=10, x=0.05,
       but similiar for other values).
     */

    r = false_position(&lo, &flo, &hi, &fhi,
                       (objective_function)gammainc, params,
                       2*MACHEP, 2*MACHEP, 1e-2*a,
                       &best_x, &best_f, &errest);
    if (!(r == FSOLVE_CONVERGED || r == FSOLVE_EXACT) &&
            errest > ALLOWED_ATOL + ALLOWED_RTOL*fabs(best_x)) {
        best_x = GAMMAINCINV_NaN;
    }
    return best_x;
}
