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
                            {1.0, 1.0, 1.0},
                            MV_CENTERING_MEAN,
                            MV_SCALING_UV,
                            1.0};
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
