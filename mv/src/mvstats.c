/*! \file mvstats.c
  \brief Contains functions to compute necessary statistics for
  \author Darryl Wallace - wallacdj@gmail.com, darryl@wallynet.ca
  \date 2012
s
  THIS CODE IS COPYRIGHT DARRYL WALLACE (c) 2012.
  UNAUTHORIZED DISTRIBUTION, USE, OR VIEWING THE CONTENTS OF THIS FILE
  IS STRICTLY PROHIBITED.
  IF YOU HAVE RECEIVED THIS CODE WITHOUT AUTHORIZATION IT MUST BE DESTROYED
  IMMEDIATELY
  */

#include <cephes/cmath/mconf.h>
#include <cephes/cmath/protos.h>
#include "mvstats.h"
#include "mvmodel.h"

// this is in gammaincinv.c
extern double gammaincinv(double alpha, double df);

double t_ppf(double alpha, double df, double loc, double scale)
{
    // unused variables:
    (void)loc; (void) scale;
    return stdtri(df, alpha);
}

double F_ppf(double alpha, double N1, double N2, double loc, double scale)
{
    // unused variables:
    (void)loc; (void) scale;
    return fdtri(N1, N2, alpha);
}

double chi2_ppf(double alpha, double df, double loc, double scale)
{
    // unused variables:
    (void)loc; (void) scale;
    return chdtri(df, 1-alpha);
}

double gamma_ppf(double alpha, double df, double loc, double scale)
{
    // unused variables:
    (void)loc; (void) scale;
    return gammaincinv(df, alpha);
}

int mvHT2(mvMat *output, const mvMat *t, const mvMat *t_stddev,
          int first_component, int last_component)
{
    int i,j;
    first_component--;
    if (!(output->nrows == t->nrows && output->ncolumns == 1 &&
          first_component > -1 && last_component <= t->nrows &&
          t_stddev->nrows == 1 && t_stddev->ncolumns == t->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i < output->nrows; i++)
    {
        output->data[i][0] = 0.0;
        for (j = first_component; j<last_component; j++)
        {
            if (t->mask[i][j] == DATA_PRESENT)
            {
                output->data[i][0] += (t->data[i][j] * t->data[i][j]) /
                        (t_stddev->data[0][j] * t_stddev->data[0][j]);
            }
        }
    }
    return SUCCESS;
}

double mvHT2Limit(double alpha, int A, int N)
{
    double _A = A;
    double _N = N;
    return _A*(_N*_N-1.0)/ (_N*(_N-_A)) * F_ppf(alpha, _A, _N-_A, 0, 0);
}

int mvSPE(mvMat *output, const mvMat *residuals)
{
    return mvMatRowSS(output, residuals);
}

double mvSPELimit(double alpha, const mvMat *modelSPE_values)
{
    double m,v,g,h;
    mvMat *mean = mvAllocMat(1,1);
    mvMat *var = mvAllocMat(1,1);

    mvColumnMean(mean, modelSPE_values);
    // sample variance.
    mvColumnVar(var, modelSPE_values, 1);
    m = mean->data[0][0];
    v = var->data[0][0];
    mvFreeMat(&var);
    mvFreeMat(&mean);

    g = v / (2.0 * m);
    h = (2.0 * m * m) / v;

    return g * chi2_ppf(alpha, h, 0, 1);
}