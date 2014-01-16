/*! \file mvpreprocess.c
  \brief Implementation of preprocess information.
  \author Darryl Wallace - wallacdj@gmail.com, darryl@wallynet.ca
  \date 2013

  THIS CODE IS COPYRIGHT DARRYL WALLACE (c) 2013.
  UNAUTHORIZED DISTRIBUTION, USE, OR VIEWING THE CONTENTS OF THIS FILE
  IS STRICTLY PROHIBITED.
  IF YOU HAVE RECEIVED THIS CODE WITHOUT AUTHORIZATION IT MUST BE DESTROYED
  IMMEDIATELY
  */


#include "mvpreprocess.h"
#include "mvconstants.h"

#include <string.h>
#include <math.h>


MVPreprocessContext * mvpreprocess_alloc()
{
    MVPreprocessContext *out = malloc(sizeof(MVPreprocessContext));
    out->preprocess_info = NULL;
    out->centering = NULL;
    out->scaling = NULL;
    return out;
}

MVPreprocessContext * mvpreprocess_alloc_init_mat(int ncolumns)
{
    MVPreprocessContext *out = mvpreprocess_alloc();
    if (out)
    {
        mvpreprocess_init_mat(out, ncolumns);
    }
    return out;
}

MVPreprocessContext * mvpreprocess_alloc_mcuv(MVMat *matrix)
{
    MVPreprocessContext *out = mvpreprocess_alloc();
    if (out)
    {
        mvpreprocess_init_mcuv(out, matrix);
    }
    return out;
}

int mvpreprocess_init_mat(MVPreprocessContext *ctx, int ncolumns)
{
    if (!ctx)
    {
        return ERROR;
    }
    int size = ncolumns * sizeof(MVPreprocessColumnInfo);
    ctx->preprocess_info = (MVPreprocessColumnInfo *) malloc(size);
    memset(ctx->preprocess_info, 0, size);
    ctx->centering = mvmat_alloc(1, ncolumns);
    ctx->scaling = mvmat_alloc(1, ncolumns);
    return SUCCESS;
}

int mvpreprocess_init_mcuv(MVPreprocessContext *ctx, MVMat *matrix)
{
    if (!ctx)
    {
        return ERROR;
    }
    mvpreprocess_init_mat(ctx, matrix->ncolumns);
    mvmat_column_mean(ctx->centering, matrix);
    mvmat_column_stddev(ctx->scaling, matrix, 1);
    int i;
    MVPreprocessColumnInfo s = {MV_TRANSFORM_NONE,
                            {1.0, 0.0, 1.0},
                            MV_CENTERING_MEAN,
                            MV_SCALING_UV,
                            1.0, NULL, NULL};
    for (i = 0; i < matrix->ncolumns; i++)
    {
        ctx->preprocess_info[i] = s;
    }

    return SUCCESS;
}

int mvpreprocess_free(MVPreprocessContext **ctx)
{
    if (ctx)
    {
        MVPreprocessContext *c = *ctx;
        if (c)
        {
            mvmat_free(&c->centering);
            mvmat_free(&c->scaling);
            free(c->preprocess_info);
            *ctx = NULL;
            return SUCCESS;
        }
    }
    return ATTEMPT_TO_FREE_NULL_MATRIX;
}

int mvpreprocess_set_column(MVPreprocessContext *ctx, int column,
                            const MVPreprocessColumnInfo *col_info)
{
    ctx->preprocess_info[column] = *col_info;
    MVPreprocessColumnInfo *c = &ctx->preprocess_info[column];
    switch (ctx->preprocess_info[column].t)
    {
    case MV_TRANSFORM_NONE:
        c->t_func = NULL;
        c->inv_t_func = NULL;
        break;
    case MV_TRANSFORM_LINEAR:
        c->t_func = mvpreprocess_linear;
        c->inv_t_func = mvpreprocess_invlinear;
        break;
    case MV_TRANSFORM_LOGARITHMIC:
        c->t_func = mvpreprocess_log;
        c->inv_t_func = mvpreprocess_invlog;
        break;
    case MV_TRANSFORM_EXPONENTIAL:
        c->t_func = mvpreprocess_exp;
        c->inv_t_func = mvpreprocess_invexp;
        break;
    case MV_TRANSFORM_POWER:
        c->t_func = mvpreprocess_power;
        c->inv_t_func = mvpreprocess_invpower;
        break;
    }

    return SUCCESS;
}

int mvpreprocess_get_column(MVPreprocessContext *ctx,
                            MVPreprocessColumnInfo *col_info_out, int column)
{
    *col_info_out = ctx->preprocess_info[column];
    return SUCCESS;
}

int mvpreprocess_prep(MVPreprocessContext *ctx, MVMat *matrix)
{
    (void) ctx; (void) matrix;
    if (ctx->centering->ncolumns != matrix->ncolumns)
    {
        return INCORRECT_DIMENSIONS;
    }
    int j;
    MVMat *tmp = mvmat_alloc(matrix->nrows, matrix->ncolumns);
    MVMAT_FUNC_PTR * funcs = (MVMAT_FUNC_PTR *)malloc(matrix->ncolumns * sizeof(MVMAT_FUNC_PTR));
    MVPreprocessCoeffs **opaques = (MVPreprocessCoeffs **)malloc(matrix->ncolumns * sizeof(MVPreprocessCoeffs *));
    for (j = 0; j < matrix->ncolumns; j++)
    {
        funcs[j] = ctx->preprocess_info[j].t_func;
        opaques[j] = (MVPreprocessCoeffs *)&ctx->preprocess_info[j].t_coeffs;
    }

    /* Apply column functions */
    mvmat_column_func(tmp, matrix, funcs, opaques);

    /* Compute centering and scaling*/
    for (j = 0; j < matrix->ncolumns; j++)
    {
        double centering = 0.0;
        double scaling = 1.0;


        switch (ctx->preprocess_info[j].c)
        {
        case MV_CENTERING_MEAN:
            mvmat_colidx_mean(&centering, matrix, j);
            break;
        case MV_CENTERING_NONE: // delibarate fall through
        default:
            centering = 0.0;
            break;
        }

        switch (ctx->preprocess_info[j].s)
        {
        case MV_SCALING_UV:
            mvmat_colidx_stddev(&scaling, matrix, 1, j);
            break;
        case MV_SCALING_PARETO:
            mvmat_colidx_stddev(&scaling, matrix, 1, j);
            if (!MVISNAN_FUNC(scaling))
            {
                scaling = sqrt(scaling);
            }
            break;
        case MV_SCALING_NONE:
        default:
            scaling = 1.0;
            break;
        }

        scaling = ctx->preprocess_info[j].multiplier / scaling;

        mvmat_set_elem(ctx->centering, 0, j, centering);
        mvmat_set_elem(ctx->scaling, 0, j, scaling);
    }

    free(funcs);
    free(opaques);
    mvmat_free(&tmp);
    return SUCCESS;
}

int mvpreprocess_do(MVPreprocessContext *ctx, MVMat *preprocessed_out, MVMat *raw_in)
{
    (void) ctx; (void) preprocessed_out; (void) raw_in;
    return SUCCESS;
}

int mvpreprocess_undo(MVPreprocessContext *ctx, MVMat *raw_out, MVMat *preprocessed_in)
{
    (void) ctx; (void) raw_out; (void) preprocessed_in;
    return SUCCESS;
}

double mvpreprocess_linear(double x, void *opaque)
{
    MVPreprocessCoeffs *coeffs = (MVPreprocessCoeffs *) opaque;

    if (!coeffs)
    {
        return x;
    }
    return coeffs->A * x + coeffs->B;
}

double mvpreprocess_invlinear(double y, void *opaque)
{
    MVPreprocessCoeffs *coeffs = (MVPreprocessCoeffs *) opaque;

    if (!coeffs)
    {
        return y;
    }
    if (coeffs->A == 0)
    {
        return mv_NaN();
    }
    return ( y - coeffs->B) / coeffs->A;
}

double mvpreprocess_log(double x, void *opaque)
{
    MVPreprocessCoeffs *coeffs = (MVPreprocessCoeffs *) opaque;

    if (!coeffs)
    {
        return x;
    }
    double tmp = coeffs->A * x + coeffs->B;
    if (tmp == 0.0)
    {
        return mv_NaN();
    }
    return log10(coeffs->A * x + coeffs->B);
}

double mvpreprocess_invlog(double y, void *opaque)
{
    MVPreprocessCoeffs *coeffs = (MVPreprocessCoeffs *) opaque;

    if (!coeffs)
    {
        return y;
    }

    if (coeffs->A == 0)
    {
        return mv_NaN();
    }
    return (pow(10.0, y) - coeffs->B) / coeffs->A;
}

double mvpreprocess_exp(double x, void *opaque)
{
    MVPreprocessCoeffs *coeffs = (MVPreprocessCoeffs *) opaque;

    if (!coeffs)
    {
        return x;
    }
    return pow(M_E, (coeffs->A * x + coeffs->B));
}

double mvpreprocess_invexp(double y, void *opaque)
{
    MVPreprocessCoeffs *coeffs = (MVPreprocessCoeffs *) opaque;

    if (!coeffs)
    {
        return y;
    }
    if (coeffs->A == 0)
    {
        return mv_NaN();
    }
    if (y == 0.0)
    {
        return mv_NaN();
    }
    return (log(y) - coeffs->B) / coeffs->A;
}


double mvpreprocess_power(double x, void *opaque)
{
    MVPreprocessCoeffs *coeffs = (MVPreprocessCoeffs *) opaque;

    if (!coeffs)
    {
        return x;
    }
    if (coeffs->C == 0.0)
    {
        return 1.0;
    }
    return pow((coeffs->A * x + coeffs->B), coeffs->C);
}

double mvpreprocess_invpower(double y, void *opaque)
{
    MVPreprocessCoeffs *coeffs = (MVPreprocessCoeffs *) opaque;

    if (!coeffs)
    {
        return y;
    }
    if (coeffs->C == 0.0 || coeffs->A == 0.0)
    {
        return mv_NaN();
    }
    return (pow(y, 1.0 / coeffs->C) - coeffs->B) / coeffs->A;
}

