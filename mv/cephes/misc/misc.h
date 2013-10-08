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

#ifndef C_MISC_MISC_H
#define C_MISC_MISC_H

typedef enum {
  /* An exact solution was found, in which case the first point
     on the interval is the value */
  FSOLVE_EXACT,
  /* Interval width is less than the tolerance */
  FSOLVE_CONVERGED,
  /* Not a bracket */
  FSOLVE_NOT_BRACKET,
  /* Root-finding didn't converge in a set number of iterations. */
  FSOLVE_MAX_ITERATIONS
} fsolve_result_t;

typedef double (*objective_function)(double, void *);

fsolve_result_t false_position(double *a, double *fa, double *b, double *fb,
                       objective_function f, void *f_extra,
                       double abserr, double relerr, double bisect_til,
                       double *best_x, double *best_f, double *errest);

double besselpoly(double a, double lambda, double nu);
double gammaincinv(double a, double x);

#define gammaincinv_doc """gammaincinv(a, y) returns x such that gammainc(a, x) = y."""

#endif /* C_MISC_MISC_H */
