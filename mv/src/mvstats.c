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

int mvHT2(MVMat *output, const MVMat *t, const MVMat *t_stddev,
          int first_component, int last_component)
{
    int i,j;
    first_component--;
    if (!(output->nrows == t->nrows && output->ncolumns == 1 &&
          first_component > -1 && last_component <= t->ncolumns &&
          t_stddev->nrows == 1 && t_stddev->ncolumns <= t->ncolumns))
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

int mvSPE(MVMat *output, const MVMat *residuals)
{
    return mvmat_row_ss(output, residuals);
}


int mvSPEXFromObs(MVMat *output, const MVModel *model, const MVMat *Xobs, const MVMat *tobs, int num_components)
{
    // todo: error checking
    MVMat *E = mvmat_alloc(Xobs->nrows, Xobs->ncolumns);
    if (tobs == NULL)
    {
        MVMat *t = mvmat_alloc(Xobs->nrows, num_components);
        mvmodel_new_obs_scores_t(t, E, Xobs, model, num_components, MV_NEW_SCORE_SCP);
        mvSPE(output, E);
        mvmat_free(&t);
    }
    else
    {
        MVMat *Xhat = E;
        mvmodel_compute_xpred(Xhat, model, tobs, num_components);
        mvmat_subtract(E, Xobs, Xhat);
        mvSPE(output, E);
    }

    mvmat_free(&E);

    return SUCCESS;
}

int mvSPEYFromObs(MVMat *output, const MVModel *model, const MVMat *Yobs, const MVMat *tobs, int num_components)
{
    // todo: error checking
    MVMat *F = mvmat_alloc(Yobs->nrows, Yobs->ncolumns);
    if (tobs == NULL)
    {
        MVMat *t = mvmat_alloc(Yobs->nrows, num_components);
        mvmodel_new_obs_scores_t(t, F, Yobs, model, num_components, MV_NEW_SCORE_SCP);
        mvSPE(output, F);
        mvmat_free(&t);
    }
    else
    {
        MVMat *Yhat = F;
        mvmodel_compute_ypred(Yhat, model, tobs, num_components);
        mvmat_subtract(F, Yobs, Yhat);
        mvSPE(output, F);
    }

    mvmat_free(&F);

    return SUCCESS;
}

double mvSPELimit(double alpha, const MVMat *modelSPE_values, int component)
{
    double m,v,g,h;
    MVMat *mean = mvmat_alloc(1, modelSPE_values->ncolumns);
    MVMat *var = mvmat_alloc(1, modelSPE_values->ncolumns);

    mvmat_column_mean(mean, modelSPE_values);
    // sample variance.
    mvmat_column_var(var, modelSPE_values, 1);
    m = mean->data[0][component-1];
    v = var->data[0][component-1];
    mvmat_free(&var);
    mvmat_free(&mean);

    g = v / (2.0 * m);
    h = (2.0 * m * m) / v;

    return g * chi2_ppf(alpha, h, 0, 1);
}
