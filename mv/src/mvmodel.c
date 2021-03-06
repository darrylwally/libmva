/*
* Filename: mvmodel.c
* Description: Implementation of structures and functions for  mv models.
* Author: Darryl Wallace <wallacdj@gmail.com
* Copyright (c) 2014 - Darryl Wallace
*
* License:
* libmva is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/.
*
* libmva is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "mvmodel.h"
#include "mvconstants.h"
#include "mvstats.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#define MAX_NIPALS_ITER 500
#define AUTOFIT_THRESHOLD 0.05

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

// Function prototypes
static int __mvmodel_add_PCA_component(MVModel *model, int perform_cv);
static int __mvmodel_add_PLS_component(MVModel *model, int perform_cv);

static int __free_cv_data(MVCrossValData **cv_data)
{
    int i;
    MVCrossValData * cvd = NULL;
    if (cv_data == NULL)
        return -1;
    cvd = *cv_data;
    if (cvd == NULL)
        return -1;

    mvmat_free (&cvd->PRESS);
    mvmat_free (&cvd->PRESSV);
    for (i=0; i< cvd->num_rounds; i++)
    {
        mvmodel_free((MVModel **)&cvd->models[i]);
    }
    free(cvd->models);
    *cv_data = NULL;
    return MV_SUCCESS;
}

static MVCrossValData * __alloc_cv_data(int num_rounds)
{
    int i;
    MVCrossValData *output = (MVCrossValData*) malloc(sizeof(MVCrossValData));
    if (!output)
        return NULL;
    output->models = malloc(num_rounds * sizeof(MVModel*));
    for (i=0; i< num_rounds; i++)
    {
        output->models[i] = NULL;
    }
    output->PRESSV = NULL;
    output->PRESS = NULL;
    output->num_rounds = num_rounds;
    return output;
}

/*! Performs column-wise regression of a vector on a matrix.

  When the regression vector(matrix) has the same number of columns (Kx1) as the
  matrix (NxK) and the result has the same number of rows as the matrix (Nx1)

  E.g., in the NIPALS algorithm t = X p / p'p

  Computed as per:
  Grunge, B., & Manne, R. (1998). Missing values in principal component analysis.
  Chemometrics and Intelligent Laboratory Systems, 42, 125-139.
  */
static int __mvmodel_regress_col(MVMat *output, const MVMat *X, const MVMat *v)
{

    /* TODO size checking.  Right now I'm avoiding it because it's an internal
       function that I think should have checked sizes already.. */
    double numerator, denominator;
    int i,j;

    for (i=0; i<X->nrows; i++)
    {
        numerator=denominator=0.0;
        for (j=0; j<X->ncolumns; j++)
        {
            if (X->mask[i][j] == DATA_PRESENT)
            {
                numerator += v->data[j][0] * X->data[i][j];
                denominator += v->data[j][0] * v->data[j][0];
            }
        }
        if (denominator==0.0)  //case where entire column is missing
        {
            output->data[i][0]=0.0;
        }
        else
        {
            output->data[i][0]=numerator/denominator;
        }
    }

    return 0;
}

/*! Performs row-wise regression of a vector on a matrix.

  When the regression vector(matrix) has the same number of rows (Nx1) as the
  matrix (NxK) and the result has the same number of columns as the matrix (Kx1)

  E.g., in the NIPALS algorithm p = t' X / t't

  Computed as per:
  Grunge, B., & Manne, R. (1998). Missing values in principal component analysis.
  Chemometrics and Intelligent Laboratory Systems, 42, 125-139.
  */
static int __mvmodel_regress_row(MVMat *output, const MVMat *X, const MVMat *v)
{
    /* TODO size checking.  Right now I'm avoiding it because it's an internal
       function that I think should have checked sizes already.. */
    double numerator, denominator;
    int i,j;
    for (j=0;j<X->ncolumns;j++)
    {
        numerator=denominator=0.0;

        for (i=0; i<X->nrows;i++)
        {
            if (X->mask[i][j] == DATA_PRESENT)
            {
                numerator += v->data[i][0] * X->data[i][j];
                denominator += v->data[i][0] * v->data[i][0];
            }
        }
        if (denominator==0.0)  //case where entire column is missing
        {
            output->data[j][0]=0.0;
        }
        else
        {
            output->data[j][0]=numerator/denominator;
        }
    }
    return 0;
}

/*! \internal WStar will be computed for the a'th component
  This is based on equation equation 30 in Dayal, B.S. & MacGregor, J.F. (1997)
  Improved PLS algorithms. J. Chemometrics 11:73-85
  \arg wstar_out is a pre-allocated vector of dim Kx1.
*/
static int __compute_wstar(MVMat *wstar_out, const MVMat * w, const MVMat *p,
                          const MVMat *cur_wstar)
{
    int i,j,a;
    int K = w->nrows;
    MVMat * wa = mvmat_alloc(K, 1); // Ath column of W (permanently)
    MVMat *p_i = mvmat_alloc(K, 1);
    MVMat *cur_wstar_i = mvmat_alloc(K,1);
    // start wStarOut as the the most recent vector of w
    a = cur_wstar->ncolumns;
    for (i=0; i< w->nrows; i++)
    {
        double val;
        mvmat_get_elem(w, &val, i, a);
        mvmat_set_elem(wstar_out, i, 0, val);
        mvmat_set_elem(wa, i, 0, val);
    }

    for (i=0; i<a; i++)
    {
        // extract the ith column of p and wStarSoFar
        double pTwa;
        for (j=0; j<K; j++)
        {
            double val;
            mvmat_get_elem(p, &val, j, i);
            mvmat_set_elem(p_i, j, 0, val);
            mvmat_get_elem(cur_wstar, &val, j, i);
            mvmat_set_elem(cur_wstar_i, j, 0, val);
        }
        pTwa = mvmat_dot_product(p_i, wa);
        mvmat_mult_scalar(cur_wstar_i, cur_wstar_i, pTwa);
        mvmat_subtract(wstar_out, wstar_out, cur_wstar_i);
    }

    mvmat_free(&wa);
    mvmat_free(&p_i);
    mvmat_free(&cur_wstar_i);

    return MV_SUCCESS;
}

static int __cross_validate_PCA_FAST(MVModel *pca)
{
    int round;
    int i;
    MVCrossValData *cv = pca->cvd;
    int numRounds = cv->num_rounds;
    MVMat * PRESS = mvmat_allocz(1,1);
    MVMat * PRESSV = mvmat_allocz(1, pca->X->ncolumns);
    MVMat * PRESSV_temp = mvmat_allocz(1, pca->X->ncolumns); // acts as PRESSV temporary addition location and (PRESSV / SSV)
    MVMat * component_slice = mvmat_alloc(1,1);
    double Q2;
    MVMat * Q2V = mvmat_alloc(1, pca->X->ncolumns);
    MVMat * X = NULL;

    // slice the correct matrix.
    if (pca->_A == 1)
    {
        X = pca->X;
    }
    else // pca->_A > 1
    {
        X = pca->E;
    }

    // Check for "Leave one out" CV mode.
    if (numRounds <= 0)
    {
        numRounds = X->nrows;
    }

    // allocate the models with the sliced rows of the residual matrix.
    for (round = 0; round < numRounds; round++)
    {
        MVModel *model_round;
        MVMat * X_model = NULL; // the matrix that the model is built
        MVMat * X_pred = NULL;  // the matrix of the data that is to be predicted
        MVMat * E_pred = NULL;  // will function as the predicted values and the residuals.
        MVMat * t_pred = NULL;
        MVMat * slice_pred = mvmat_range(round, pca->E->nrows, numRounds);

        mvmat_slice_rows_ref(&X_pred, X, slice_pred);
        mvmat_delete_rows_ref(&X_model, X, slice_pred);

        E_pred = mvmat_alloc(X_pred->nrows, X_pred->ncolumns);
        t_pred = mvmat_alloc(X_pred->nrows, 1);

        // Initialize model, add one component, compute new observation scores.
        model_round = mvmodel_alloc_init_pca(X_model);
        __mvmodel_add_PCA_component(model_round, 0);
        mvmodel_t_scores_from_obs(t_pred, E_pred, X_pred, model_round, 1, MV_NEW_SCORE_SCP);

        // Sum the columns of E_pred and then add those to PRESSV
        mvmat_column_ss(PRESSV_temp, E_pred);
        mvmat_add(PRESSV, PRESSV, PRESSV_temp);

        mvmat_free(&slice_pred);
        mvmat_free(&t_pred);
        mvmat_free(&E_pred);
        mvmat_free(&X_pred);
        mvmat_free(&X_model);
        mvmodel_free(&model_round);
    }
    mvmat_set_elem(PRESS, 0, 0, mvmat_sum(PRESSV));
    // Compute Q2 and Q2V
    {
        MVMat * SSXV_ref = NULL;
        component_slice->data[0][0] = pca->_A - 1;
        Q2 = 1.0 - PRESS->data[0][0] / pca->SSX->data[pca->_A-1][0];
        mvmat_slice_rows_ref(&SSXV_ref, pca->SSXV, component_slice);
        mvmat_column_div(PRESSV_temp, PRESSV, SSXV_ref);
        mvmat_mult_scalar(PRESSV_temp, PRESSV, -1.0);
        mvmat_add_scalar(Q2V, PRESSV_temp, 1.0);
        mvmat_free(&SSXV_ref);
    }

    if (pca->_A == 1)
    {
        // Q2's have not yet been allocated.
        pca->Q2 = mvmat_alloc(1,1);
        pca->Q2cum = mvmat_alloc(1,1);
        pca->Q2V = Q2V;
        pca->Q2Vcum = mvmat_alloc(1,X->ncolumns);

        cv->PRESS = PRESS;
        cv->PRESSV = PRESSV;

        // for the first component, Q2[V] = Q2[V]cum
        mvmat_copy(pca->Q2Vcum, pca->Q2V);
        mvmat_set_elem(pca->Q2, 0, 0, Q2);
        mvmat_set_elem(pca->Q2cum, 0, 0, Q2);
    }
    else
    {
        double Q2cum=1.0;
        MVMat *newPRESS, *newPRESSV, *PRESSV_SSV;
        MVMat *newQ2, *newQ2V, *newQ2cum, *newQ2Vcum, *newQ2Vcum_ref;
        newPRESS = mvmat_alloc(pca->_A, 1);
        newPRESSV = mvmat_alloc(pca->_A, pca->X->ncolumns);
        mvmat_concat_rows(newPRESS, cv->PRESS, PRESS);
        mvmat_concat_rows(newPRESSV, cv->PRESSV, PRESSV);
        mvmat_free(&cv->PRESS);
        mvmat_free(&cv->PRESSV);
        mvmat_free(&PRESS);
        mvmat_free(&PRESSV);
        cv->PRESS = newPRESS;
        cv->PRESSV = newPRESSV;

        // new Q2
        newQ2 = mvmat_alloc(pca->_A, 1);
        for (i=0; i<pca->_A-1; i++)
        {
            newQ2->data[i][0] = pca->Q2->data[i][0];
        }
        newQ2->data[pca->_A-1][0] = Q2;
        mvmat_free(&pca->Q2);
        pca->Q2 = newQ2;

        // new Q2V
        newQ2V = mvmat_alloc(pca->_A, X->ncolumns);
        mvmat_concat_rows(newQ2V, pca->Q2V, Q2V);
        mvmat_free(&pca->Q2V);
        mvmat_free(&Q2V);
        pca->Q2V = newQ2V;

        // new Q2cum
        newQ2cum = mvmat_alloc(pca->_A, 1);
        for(i=0; i<pca->_A-1; i++)
        {
            newQ2cum->data[i][0] = pca->Q2cum->data[i][0];
            Q2cum *= cv->PRESS->data[i][0] / pca->SSX->data[i][0];
        }
        Q2cum *= cv->PRESS->data[pca->_A-1][0] / pca->SSX->data[pca->_A-1][0];
        newQ2cum->data[pca->_A-1][0] = 1.0 - Q2cum;
        mvmat_free(&pca->Q2cum);
        pca->Q2cum = newQ2cum;

        // new Q2Vcum
        newQ2Vcum = mvmat_alloc(pca->_A, X->ncolumns);
        // first copy existing data into newQVcum
        for (i=0; i < pca->_A-1; i++)
        {
            int j;
            for (j=0; j<newQ2Vcum->ncolumns; j++)
            {
                newQ2Vcum->data[i][j] = pca->Q2Vcum->data[i][j];
            }
        }
        // set PRESSV_temp to be full of 1's
        mvmat_set(PRESSV_temp, 1.0);
        PRESSV_SSV = mvmat_alloc(1, X->ncolumns);
        for (i=0; i < pca->_A; i++)
        {
            MVMat *PRESSV_ref, *SSXV_ref;
            PRESSV_ref = NULL;
            SSXV_ref = NULL;
            component_slice->data[0][0] = i;
            mvmat_slice_rows_ref(&PRESSV_ref, cv->PRESSV, component_slice);
            mvmat_slice_rows_ref(&SSXV_ref, pca->SSXV, component_slice);
            mvmat_column_div(PRESSV_SSV, PRESSV_ref, SSXV_ref);
            mvmat_elem_mult(PRESSV_temp, PRESSV_temp, PRESSV_SSV);

            mvmat_free(&PRESSV_ref);
            mvmat_free(&SSXV_ref);
        }
        component_slice->data[0][0] = pca->_A - 1;
        newQ2Vcum_ref = NULL;
        mvmat_slice_rows_ref(&newQ2Vcum_ref, newQ2Vcum, component_slice);
        mvmat_mult_scalar(PRESSV_temp, PRESSV_SSV, -1.0);
        mvmat_add_scalar(newQ2Vcum_ref, PRESSV_temp, 1.0);
        mvmat_free(&PRESSV_SSV);
        mvmat_free(&newQ2Vcum_ref);
        mvmat_free(&pca->Q2Vcum);
        pca->Q2Vcum = newQ2Vcum;
    }
    mvmat_free(&component_slice);
    mvmat_free(&PRESSV_temp);
    return 0;
}

static int __cross_validate_PCA(MVModel *pca)
{
    switch (pca->cv_type)
    {
    case MV_CROSSVAL_TYPE_FAST:
        return __cross_validate_PCA_FAST(pca);
//    case FULL:
//        return __crossValidatePCA_FULL(pca);
    default:
        return MV_UNKNOWN_MODEL_TYPE;
    }

    return 0;
}

static int __cross_validate_PLS_FAST(MVModel *pls)
{

    int round;
    int i;
    MVCrossValData *cv = pls->cvd;
    int num_rounds = cv->num_rounds;
    MVMat * PRESS = mvmat_allocz(1,1);
    MVMat * PRESSV = mvmat_allocz(1, pls->Y->ncolumns);
    MVMat * PRESSV_temp = mvmat_allocz(1, pls->Y->ncolumns); // acts as PRESSV temporary addition location and (PRESSV / SSV)
    MVMat * component_slice = mvmat_alloc(1,1);
    double Q2;
    MVMat * Q2V = mvmat_alloc(1, pls->Y->ncolumns);
    MVMat * X = NULL;
    MVMat * Y = NULL;

    // slice the correct matrix.
    if (pls->_A == 1)
    {
        X = pls->X;
        Y = pls->Y;
    }
    else // pls->_A > 1
    {
        X = pls->E;
        Y = pls->F;
    }

    // Check for "Leave one out" CV mode.
    if (num_rounds <= 0)
    {
        num_rounds = X->nrows;
    }

    // allocate the models with the sliced rows of the residual matrix.
    for (round = 0; round < num_rounds; round++)
    {
        MVModel *model_round;
        MVMat * Y_model = NULL; // the matrix that the model is built
        MVMat * Y_pred = NULL;  // the matrix of the data that is to be predicted
        MVMat * X_model = NULL; // the matrix that the model is built
        MVMat * X_pred = NULL;  // the matrix of the data that is to be predicted
        MVMat * F_pred = NULL;  // will function as the predicted values and the residuals.
        MVMat * t_pred = NULL;
        MVMat * u_pred = NULL;
        MVMat * slice_pred = mvmat_range(round, pls->F->nrows, num_rounds);

        mvmat_slice_rows_ref(&Y_pred, Y, slice_pred);
        mvmat_delete_rows_ref(&Y_model, Y, slice_pred);
        mvmat_slice_rows_ref(&X_pred, X, slice_pred);
        mvmat_delete_rows_ref(&X_model, X, slice_pred);

        F_pred = mvmat_alloc(Y_pred->nrows, Y_pred->ncolumns);
        t_pred = mvmat_alloc(Y_pred->nrows, 1);
        u_pred = mvmat_alloc(Y_pred->nrows, 1);

        // Initialize model, add one component, compute new observation scores.
        model_round = mvmodel_alloc_init_pls(X_model, Y_model);
        __mvmodel_add_PLS_component(model_round, 0);
        mvmodel_t_scores_from_obs(t_pred, NULL, X_pred, model_round, 1, MV_NEW_SCORE_SCP);
        mvmodel_u_scores_from_obs(u_pred, F_pred, Y_pred, t_pred, model_round, 1, MV_NEW_SCORE_SCP);

        // Sum the columns of F_pred and then add those to PRESSV
        mvmat_column_ss(PRESSV_temp, F_pred);
        mvmat_add(PRESSV, PRESSV, PRESSV_temp);


        mvmat_free(&slice_pred);
        mvmat_free(&u_pred);
        mvmat_free(&t_pred);
        mvmat_free(&F_pred);
        mvmat_free(&X_pred);
        mvmat_free(&X_model);
        mvmat_free(&Y_pred);
        mvmat_free(&Y_model);
        mvmodel_free(&model_round);
    }
    mvmat_set_elem(PRESS, 0, 0, mvmat_sum(PRESSV));
    // Compute Q2 and Q2V
    {
        MVMat * SSYV_ref = NULL;
        component_slice->data[0][0] = pls->_A - 1;
        Q2 = 1.0 - PRESS->data[0][0] / pls->SSY->data[pls->_A-1][0];
        mvmat_slice_rows_ref(&SSYV_ref, pls->SSYV, component_slice);
        mvmat_column_div(PRESSV_temp, PRESSV, SSYV_ref);
        mvmat_mult_scalar(PRESSV_temp, PRESSV, -1.0);
        mvmat_add_scalar(Q2V, PRESSV_temp, 1.0);
        mvmat_free(&SSYV_ref);
    }

    if (pls->_A == 1)
    {
        // Q2's have not yet been allocated.
        pls->Q2 = mvmat_alloc(1,1);
        pls->Q2cum = mvmat_alloc(1,1);
        pls->Q2V = Q2V;
        pls->Q2Vcum = mvmat_alloc(1, Y->ncolumns);

        cv->PRESS = PRESS;
        cv->PRESSV = PRESSV;

        // for the first component, Q2[V] = Q2[V]cum
        mvmat_copy(pls->Q2Vcum, pls->Q2V);
        mvmat_set_elem(pls->Q2, 0, 0, Q2);
        mvmat_set_elem(pls->Q2cum, 0, 0, Q2);
    }
    else
    {
        double Q2cum=1.0;
        MVMat *newPRESS, *newPRESSV, *PRESSV_SSV;
        MVMat *newQ2, *newQ2V, *newQ2cum, *newQ2Vcum, *newQ2Vcum_ref;
        newPRESS = mvmat_alloc(pls->_A, 1);
        newPRESSV = mvmat_alloc(pls->_A, pls->Y->ncolumns);
        mvmat_concat_rows(newPRESS, cv->PRESS, PRESS);
        mvmat_concat_rows(newPRESSV, cv->PRESSV, PRESSV);
        mvmat_free(&cv->PRESS);
        mvmat_free(&cv->PRESSV);
        mvmat_free(&PRESS);
        mvmat_free(&PRESSV);
        cv->PRESS = newPRESS;
        cv->PRESSV = newPRESSV;

        // new Q2
        newQ2 = mvmat_alloc(pls->_A, 1);
        for (i=0; i<pls->_A-1; i++)
        {
            newQ2->data[i][0] = pls->Q2->data[i][0];
        }
        newQ2->data[pls->_A-1][0] = Q2;
        mvmat_free(&pls->Q2);
        pls->Q2 = newQ2;

        // new Q2V
        newQ2V = mvmat_alloc(pls->_A, Y->ncolumns);
        mvmat_concat_rows(newQ2V, pls->Q2V, Q2V);
        mvmat_free(&pls->Q2V);
        mvmat_free(&Q2V);
        pls->Q2V = newQ2V;

        // new Q2cum
        newQ2cum = mvmat_alloc(pls->_A, 1);
        for(i=0; i<pls->_A-1; i++)
        {
            newQ2cum->data[i][0] = pls->Q2cum->data[i][0];
            Q2cum *= cv->PRESS->data[i][0] / pls->SSY->data[i][0];
        }
        Q2cum *= cv->PRESS->data[pls->_A-1][0] / pls->SSY->data[pls->_A-1][0];
        newQ2cum->data[pls->_A-1][0] = 1.0 - Q2cum;
        mvmat_free(&pls->Q2cum);
        pls->Q2cum = newQ2cum;

        // new Q2Vcum
        newQ2Vcum = mvmat_alloc(pls->_A, Y->ncolumns);
        // first copy existing data into newQVcum
        for (i=0; i < pls->_A-1; i++)
        {
            int j;
            for (j=0; j<newQ2Vcum->ncolumns; j++)
            {
                newQ2Vcum->data[i][j] = pls->Q2Vcum->data[i][j];
            }
        }
        // set PRESSV_temp to be full of 1's
        mvmat_set(PRESSV_temp, 1.0);
        PRESSV_SSV = mvmat_alloc(1, Y->ncolumns);
        for (i=0; i < pls->_A; i++)
        {
            MVMat *PRESSV_ref, *SSYV_ref;
            PRESSV_ref = NULL;
            SSYV_ref = NULL;
            component_slice->data[0][0] = i;
            mvmat_slice_rows_ref(&PRESSV_ref, cv->PRESSV, component_slice);
            mvmat_slice_rows_ref(&SSYV_ref, pls->SSYV, component_slice);
            mvmat_column_div(PRESSV_SSV, PRESSV_ref, SSYV_ref);
            mvmat_elem_mult(PRESSV_temp, PRESSV_temp, PRESSV_SSV);

            mvmat_free(&PRESSV_ref);
            mvmat_free(&SSYV_ref);
        }
        component_slice->data[0][0] = pls->_A - 1;
        newQ2Vcum_ref = NULL;
        mvmat_slice_rows_ref(&newQ2Vcum_ref, newQ2Vcum, component_slice);
        mvmat_mult_scalar(PRESSV_temp, PRESSV_SSV, -1.0);
        mvmat_add_scalar(newQ2Vcum_ref, PRESSV_temp, 1.0);
        mvmat_free(&PRESSV_SSV);
        mvmat_free(&newQ2Vcum_ref);
        mvmat_free(&pls->Q2Vcum);
        pls->Q2Vcum = newQ2Vcum;
    }
    mvmat_free(&component_slice);
    mvmat_free(&PRESSV_temp);
    return 0;
}

static int __cross_validate_PLS(MVModel *pls)
{
    switch (pls->cv_type)
    {
    case MV_CROSSVAL_TYPE_FAST:
        return __cross_validate_PLS_FAST(pls);
//    case FULL:
//        return __crossValidatePLS_FULL(pls);
    default:
        return MV_UNKNOWN_MODEL_TYPE;
    }

    return 0;
}

MVModel * mvmodel_alloc_init_pca(MVMat *X)
{
    MVModel * output = (MVModel *) malloc(sizeof(MVModel));
    if (!output)
    {
        return NULL;
    }
    output->model_type = MV_MODEL_TYPE_PCA;
    output->X = X;
    output->E = mvmat_allocz(X->nrows, X->ncolumns);
    output->p = NULL;
    output->t = NULL;
    output->t_stddev = NULL;
    output->R2X = NULL;
    output->cvd = NULL;
    output->cv_type = MV_CROSSVAL_TYPE_FAST;
    output->num_cv_rounds = 7;
    output->A = output->_A = 0;
    output->SSX = mvmat_alloc(1,1);
    output->SSXV = mvmat_alloc(1, X->ncolumns);
    mvmat_column_ss(output->SSXV, X);
    output->SSX->data[0][0] = mvmat_sum(output->SSXV);
    output->Q2 = NULL;
    output->Q2V = NULL;
    output->Q2cum = NULL;
    output->Q2Vcum = NULL;
    output->SPEX = NULL;
    output->SPEY = NULL;
    /* Set PLS options to NULL just because */
    output->SSY = NULL;
    output->SSYV = NULL;
    output->Y = NULL;
    output->F = NULL;
    output->w = NULL;
    output->u = NULL;
    output->wstar = NULL;
    output->c = NULL;
    output->R2Y =NULL;
    /* TODO ADD Cross validation stuff */
    return output;
}

MVModel * mvmodel_alloc_init_pls(MVMat *X, MVMat *Y)
{
    MVModel *output;
    if (X->nrows != Y->nrows)
    {
        return NULL;
    }
    output = (MVModel *) malloc(sizeof(MVModel));
    if (!output)
    {
        return NULL;
    }
    output->model_type = MV_MODEL_TYPE_PLS;
    output->X = X;
    output->Y = Y;
    output->E = mvmat_allocz(X->nrows, X->ncolumns);
    output->F = mvmat_allocz(Y->nrows, Y->ncolumns);
    output->u = NULL;
    output->w = NULL;
    output->wstar = NULL;
    output->p = NULL;
    output->t = NULL;
    output->t_stddev = NULL;
    output->cvd = NULL;
    output->cv_type = MV_CROSSVAL_TYPE_FAST;
    output->num_cv_rounds = 7;
    output->SSX = mvmat_alloc(1,1);
    output->SSXV = mvmat_alloc(1, X->ncolumns);
    mvmat_column_ss(output->SSXV, X);
    output->SSX->data[0][0] = mvmat_sum(output->SSXV);
    output->SSY = mvmat_alloc(1,1);
    output->SSYV = mvmat_alloc(1, Y->ncolumns);
    mvmat_column_ss(output->SSYV, Y);
    output->SSY->data[0][0] = mvmat_sum(output->SSYV);
    output->SPEX = NULL;
    output->SPEY = NULL;
    output->R2X = NULL;
    output->R2Y = NULL;
    output->Q2 = NULL;
    output->Q2V = NULL;
    output->Q2cum = NULL;
    output->Q2Vcum = NULL;
    output->A = output->_A = 0;
    /* TODO ADD Cross validation stuff */
    return output;

}

static int __mv_free_PCA_model(MVModel **model)
{
    MVModel *m = NULL;
    if (!model)
        return -1;
    m = *model;
    if (!m)
        return -1;
    mvmat_free(&m->E);
    mvmat_free(&m->t);
    mvmat_free(&m->t_stddev);
    mvmat_free(&m->p);
    mvmat_free(&m->R2X);
    mvmat_free(&m->SSX);
    mvmat_free(&m->SSXV);
    mvmat_free(&m->R2X);
    mvmat_free(&m->Q2);
    mvmat_free(&m->Q2V);
    mvmat_free(&m->Q2cum);
    mvmat_free(&m->Q2Vcum);
    mvmat_free(&m->SPEX);
    __free_cv_data((MVCrossValData**) &m->cvd);
    free(m);
    *model = NULL;
    return MV_SUCCESS;
}

static int __mv_free_PLS_model(MVModel **model)
{
    MVModel *m = NULL;
    if (!model)
        return -1;
    m = *model;
    if (!m)
        return -1;
    mvmat_free(&m->E);
    mvmat_free(&m->F);
    mvmat_free(&m->u);
    mvmat_free(&m->w);
    mvmat_free(&m->wstar);
    mvmat_free(&m->t);
    mvmat_free(&m->t_stddev);
    mvmat_free(&m->p);
    mvmat_free(&m->R2X);
    mvmat_free(&m->R2Y);
    mvmat_free(&m->SSX);
    mvmat_free(&m->SSXV);
    mvmat_free(&m->SSY);
    mvmat_free(&m->SSYV);
    mvmat_free(&m->Q2);
    mvmat_free(&m->Q2V);
    mvmat_free(&m->Q2cum);
    mvmat_free(&m->Q2Vcum);
    mvmat_free(&m->SPEX);
    mvmat_free(&m->SPEY);
    __free_cv_data((MVCrossValData **) &m->cvd);
    free(m);
    *model = NULL;
    return MV_SUCCESS;
}

int mvmodel_free(MVModel **model)
{
    MVModel *m;
    m = *model;
    if (!m)
    {
        return -1;
    }
    switch (m->model_type)
    {
    case MV_MODEL_TYPE_PCA:
        return __mv_free_PCA_model(model);
    case MV_MODEL_TYPE_PLS:
        return __mv_free_PLS_model(model);
    default:
        return MV_UNKNOWN_MODEL_TYPE;
    }
    return MV_UNKNOWN_MODEL_TYPE;
}

static int __mvmodel_add_PCA_component(MVModel *model, int perform_cv)
{
    int num_iter = 0;
    int N = model->X->nrows;
    int K = model->X->ncolumns;
    MVMat *p = mvmat_alloc_setval(K, 1, 1.0);
    MVMat *p_last = mvmat_alloc(K, 1);
    MVMat *p_diff = mvmat_alloc(K, 1);
    MVMat *t = mvmat_alloc(N, 1);
    MVMat *R2 = mvmat_alloc(1,1);
    MVMat *iter = mvmat_allocz(1,1);
    MVMat *X;   // reference -> no clean up req'd.
    MVMat *pT, *tpT;
    double vector_norm;
    double SSE;
    MVMat *SSEV = mvmat_alloc(1, K);

    // Initialize p as a unit length vector.
    mvmat_mult_scalar(p, p, mvmat_vector_norm(p));

    if( model->_A==0)
    {
        X = model->X;
    }
    else
    {
        X = model->E;
    }
    // PCA NIPALS Algorithm
    do{
        num_iter++;
        __mvmodel_regress_col(t, X, p);    // t= Xp / (p'p);

        mvmat_copy(p_last, p);         // store contents of p into pOld

        __mvmodel_regress_row(p, X, t);    // p = t'X / (t't);

        vector_norm = mvmat_vector_norm(p);

        mvmat_mult_scalar(p, p, 1.0/vector_norm);  // normalize p

        mvmat_subtract(p_diff, p, p_last);          // get the difference of P

        vector_norm = mvmat_vector_norm(p_diff);
    } while (vector_norm > MV_SQRT_EPS && num_iter < MAX_NIPALS_ITER);

    // no longer need, pOld or pDiff
    mvmat_free(&p_last);
    mvmat_free(&p_diff);

    // Increment number of internal components
    model->_A++;

    // Set iter
    mvmat_set_elem(iter, 0, 0 , num_iter);

    // XXX: Cross validation must be performed before the new residual E is
    // computed and depends on _A == 1 for the first component (NOT ZERO!)
    if (perform_cv)
    {
        if (!model->cvd)
        {
            model->cvd = __alloc_cv_data(model->num_cv_rounds);
        }
        __cross_validate_PCA(model);
    }

    // compute residual
    tpT = mvmat_alloc(X->nrows, X->ncolumns);
    pT = mvmat_alloc(1, X->ncolumns);
    mvmat_transpose(pT, p);
    mvmat_mult(tpT, t, pT);
    mvmat_subtract(model->E, X, tpT);
    mvmat_free(&pT);
    mvmat_free(&tpT);

    // Compute R2;
    mvmat_column_ss(SSEV, model->E);
    SSE = mvmat_ss(model->E);
    mvmat_set_elem(R2, 0, 0, 1.0 - SSE/model->SSX->data[0][0]);

    // store new components.
    if (model->_A > 1)
    {
        MVMat *newIter = mvmat_alloc(model->_A, 1);
        MVMat *newR2 = mvmat_alloc(model->_A, 1);
        MVMat *newT = mvmat_alloc(model->X->nrows, model->_A);
        MVMat *newP = mvmat_alloc(model->X->ncolumns, model->_A);
        MVMat *newSPE = mvmat_alloc(model->X->nrows, model->_A);
        MVMat *SPE = mvmat_alloc(model->X->nrows, 1);

        // T and P
        mvmat_concat_columns(newT, model->t, t);
        mvmat_free(&model->t);
        model->t=newT;
        mvmat_concat_columns(newP, model->p, p);
        mvmat_free(&model->p);
        model->p=newP;
        mvmat_free(&t);
        mvmat_free(&p);
        mvmat_free(&model->t_stddev);
        model->t_stddev = mvmat_alloc(1, model->_A);
        mvmat_column_stddev(model->t_stddev, model->t, 1);

        // R2X
        mvmat_concat_rows(newR2, model->R2X, R2);
        mvmat_free(&model->R2X);
        model->R2X=newR2;
        mvmat_free(&R2);

        // SPE
        mvstats_spe(SPE, model->E);
        mvmat_concat_columns(newSPE, model->SPEX, SPE);
        mvmat_free(&model->SPEX);
        mvmat_free(&SPE);
        model->SPEX = newSPE;

        // Iter
        mvmat_concat_rows(newIter, model->iter, iter);
        mvmat_free(&model->iter);
        model->iter = newIter;
        mvmat_free(&iter);

    }
    else
    {
        model->t = t;
        model->t_stddev = mvmat_alloc(1,1);
        mvmat_column_stddev(model->t_stddev, model->t, 1);
        model->p = p;
        model->R2X = R2;
        model->iter = iter;
        model->SPEX = mvmat_alloc(model->E->nrows, 1);
        mvstats_spe(model->SPEX, model->E);
    }

    // Store new sum of squares values.
    {
        int i=0;
        MVMat *SSX = NULL;
        MVMat *SSXV = NULL;

        //SSX
        SSX = mvmat_alloc(model->SSX->nrows+1, 1);
        // Copy the data in so that we don't need another container.
        for (i=0; i<model->SSX->nrows; i++)
        {
            SSX->data[i][0] = model->SSX->data[i][0];
        }
        SSX->data[model->SSX->nrows][0] = SSE;
        mvmat_free(&model->SSX);
        model->SSX = SSX;

        // SSXV
        SSXV = mvmat_alloc(model->SSXV->nrows+1, 1);
        mvmat_concat_rows(SSXV, model->SSXV, SSEV);
        mvmat_free(&model->SSXV);
        mvmat_free(&SSEV);
        model->SSXV = SSXV;
    }

    model->A=model->_A;
    return 0;
}

static int __mvmodel_add_PLS_component(MVModel *model, int perform_cv)
{
    int num_iter = 0;
    int N = model->X->nrows;
    int K = model->X->ncolumns;
    int M = model->Y->ncolumns;
    MVMat *w = mvmat_alloc_setval(K, 1, 1.0);
    MVMat *w_last = mvmat_alloc(K, 1);
    MVMat *w_diff = mvmat_alloc(K, 1);
    MVMat *p = mvmat_alloc(K, 1);
    MVMat *c = mvmat_alloc(M, 1);
    MVMat *t = mvmat_alloc(N, 1);
    MVMat *u = mvmat_alloc(N, 1);
    MVMat *R2X = mvmat_alloc(1, 1);
    MVMat *R2Y = mvmat_alloc(1, 1);
    MVMat *iter = mvmat_allocz(1, 1);
    MVMat *X, *Y;   // reference -> no clean up req'd.
    MVMat *pT, *tpT, *cT, *tcT;
    double SSE, SSF;
    MVMat *SSEV = mvmat_alloc(1, K);
    MVMat *SSFV = mvmat_alloc(1, M);

    // Initialize w as a unit length vector.
    mvmat_mult_scalar(w, w, mvmat_vector_norm(w));

    if( model->_A==0)
    {
        X = model->X;
        Y = model->Y;
    }
    else
    {
        X = model->E;
        Y = model->F;
    }
    // PCA NIPALS Algorithm
    do{
        num_iter++;

        __mvmodel_regress_col(t, X, w);    // t = X w / w'w

        __mvmodel_regress_row(c, Y, t);    // c = t'Y / t't

        __mvmodel_regress_col(u, Y, c);    // u = Y c / c'c

        mvmat_copy(w_last, w);

        __mvmodel_regress_row(w, X, u);    // w = u'X / u'u

        mvmat_mult_scalar(w, w, 1.0/mvmat_vector_norm(w));  // normalize w

        mvmat_subtract(w_diff, w, w_last);          // get the difference of w

    } while (mvmat_vector_norm(w_diff)>MV_SQRT_EPS && num_iter < MAX_NIPALS_ITER);

    // compute loading p after loop
    __mvmodel_regress_row(p, X, t);

    // no longer need, wOld or wDiff
    mvmat_free(&w_last);
    mvmat_free(&w_diff);

    // Increment number of internal components
    model->_A++;

    // Set iter
    mvmat_set_elem(iter, 0, 0 , num_iter);

    // XXX: Cross validation must be performed before the new residual E is
    // computed and depends on _A == 1 for the first component (NOT ZERO!)
    if (perform_cv)
    {
        if (!model->cvd)
        {
            model->cvd = __alloc_cv_data(model->num_cv_rounds);
        }
        __cross_validate_PLS(model);
    }

    // compute residual E
    tpT = mvmat_alloc(X->nrows, X->ncolumns);
    pT = mvmat_alloc(1, X->ncolumns);
    mvmat_transpose(pT, p);
    mvmat_mult(tpT, t, pT);
    mvmat_subtract(model->E, X, tpT);
    mvmat_free(&pT);
    mvmat_free(&tpT);

    // Compute residual F
    tcT = mvmat_alloc (Y->nrows, Y->ncolumns);
    cT = mvmat_alloc(1, Y->ncolumns);
    mvmat_transpose(cT, c);
    mvmat_mult(tcT, t, cT);
    mvmat_subtract(model->F, Y, tcT);
    mvmat_free(&cT);
    mvmat_free(&tcT);

    // Compute R2;
    mvmat_column_ss(SSEV, model->E);
    mvmat_column_ss(SSFV, model->F);
    SSE = mvmat_sum(SSEV);
    SSF = mvmat_sum(SSFV);
    mvmat_set_elem(R2X, 0, 0, 1.0 - SSE / model->SSX->data[0][0]);
    mvmat_set_elem(R2Y, 0, 0, 1.0 - SSF / model->SSY->data[0][0]);

    // store new components
    if (model->_A > 1)
    {
        MVMat *newT, *newP, *newW, *newU, *newC, *wStar, *newWStar;
        MVMat *newR2X, *newR2Y, *newIter, *newSPE, *SPE;
        //t
        newT = mvmat_alloc(model->X->nrows, model->A+1);
        mvmat_concat_columns(newT, model->t, t);
        mvmat_free(&model->t);
        model->t=newT;
        mvmat_free(&t);
        mvmat_free(&model->t_stddev);
        model->t_stddev = mvmat_alloc(1, model->_A);
        mvmat_column_stddev(model->t_stddev, model->t, 1);

        //p
        newP = mvmat_alloc(model->X->ncolumns, model->A+1);
        mvmat_concat_columns(newP, model->p, p);
        mvmat_free(&model->p);
        model->p=newP;
        mvmat_free(&p);

        //w
        newW = mvmat_alloc(model->X->ncolumns, model->A+1);
        mvmat_concat_columns(newW, model->w, w);
        mvmat_free(&model->w);
        model->w=newW;
        mvmat_free(&w);

        //u
        newU = mvmat_alloc(model->Y->nrows, model->A+1);
        mvmat_concat_columns(newU, model->u, u);
        mvmat_free(&model->u);
        model->u=newU;
        mvmat_free(&u);

        //c
        newC = mvmat_alloc(model->Y->ncolumns, model->A+1);
        mvmat_concat_columns(newC, model->c, c);
        mvmat_free(&model->c);
        model->c=newC;
        mvmat_free(&c);

        //W*
        wStar = mvmat_alloc(K, 1);
        newWStar = mvmat_alloc(model->X->ncolumns, model->A+1);
        __compute_wstar(wStar, model->w, model->p, model->wstar);
        mvmat_concat_columns(newWStar, model->wstar, wStar);
        mvmat_free(&model->wstar);
        model->wstar = newWStar;
        mvmat_free(&wStar);

        //R2X
        newR2X = mvmat_alloc(model->A+1, 1);
        mvmat_concat_rows(newR2X, model->R2X, R2X);
        mvmat_free(&model->R2X);
        model->R2X = newR2X;
        mvmat_free(&R2X);

        //R2Y
        newR2Y = mvmat_alloc(model->A+1, 1);
        mvmat_concat_rows(newR2Y, model->R2Y, R2Y);
        mvmat_free(&model->R2Y);
        model->R2Y = newR2Y;
        mvmat_free(&R2Y);

        // SPE X and Y
        SPE = mvmat_alloc(model->X->nrows, 1);
        newSPE = mvmat_alloc(model->X->nrows, model->_A);
        mvstats_spe(SPE, model->E);
        mvmat_concat_columns(newSPE, model->SPEX, SPE);
        mvmat_free(&model->SPEX);
        model->SPEX = newSPE;

        newSPE = mvmat_alloc(model->Y->nrows, model->_A);
        mvstats_spe(SPE, model->F);
        mvmat_concat_columns(newSPE, model->SPEY, SPE);
        mvmat_free(&model->SPEY);
        model->SPEY = newSPE;
        mvmat_free(&SPE);


        // Iter
        newIter = mvmat_alloc(model->A+1, 1);
        mvmat_concat_rows(newIter, model->iter, iter);
        mvmat_free(&model->iter);
        model->iter = newIter;
        mvmat_free(&iter);

    }
    else
    {
        model->t=t;
        model->t_stddev = mvmat_alloc(1,1);
        mvmat_column_stddev(model->t_stddev, model->t, 1);
        model->p=p;
        model->c=c;
        model->u=u;
        model->w=w;
        model->wstar = mvmat_alloc_copy(w); // W* = W for the first component
        model->R2X = R2X;
        model->R2Y = R2Y;
        model->iter = iter;
        model->SPEX = mvmat_alloc(model->E->nrows, 1);
        mvstats_spe(model->SPEX, model->E);
        model->SPEY = mvmat_alloc(model->F->nrows, 1);
        mvstats_spe(model->SPEY, model->F);
    }

    // Store new sum of squares values.
    {
        int i=0;
        MVMat *SSX = NULL;
        MVMat *SSXV = NULL;
        MVMat *SSY = NULL;
        MVMat *SSYV = NULL;

        //SSX
        SSX = mvmat_alloc(model->SSX->nrows+1, 1);
        // Copy the data in so that we don't need another container.
        for (i=0; i<model->SSX->nrows; i++)
        {
            SSX->data[i][0] = model->SSX->data[i][0];
        }
        SSX->data[model->SSX->nrows][0] = SSE;
        mvmat_free(&model->SSX);
        model->SSX = SSX;

        // SSXV
        SSXV = mvmat_alloc(model->SSXV->nrows+1, 1);
        mvmat_concat_rows(SSXV, model->SSXV, SSEV);
        mvmat_free(&model->SSXV);
        mvmat_free(&SSEV);
        model->SSXV = SSXV;

        //SSY
        SSY = mvmat_alloc(model->SSY->nrows+1, 1);
        // Copy the data in so that we don't need another container.
        for (i=0; i<model->SSY->nrows; i++)
        {
            SSY->data[i][0] = model->SSY->data[i][0];
        }
        SSY->data[model->SSY->nrows][0] = SSF;
        mvmat_free(&model->SSY);
        model->SSY = SSY;

        // SSYV
        SSYV = mvmat_alloc(model->SSYV->nrows+1, 1);
        mvmat_concat_rows(SSYV, model->SSYV, SSFV);
        mvmat_free(&model->SSYV);
        mvmat_free(&SSFV);
        model->SSYV = SSYV;
    }

    model->A=model->_A;
    return 0;
}

int mvmodel_add_component(MVModel *model)
{
    switch(model->model_type)
    {
    case MV_MODEL_TYPE_PCA:
        return __mvmodel_add_PCA_component(model, 1);
    case MV_MODEL_TYPE_PLS:
        return __mvmodel_add_PLS_component(model, 1);
    default:
        return MV_UNKNOWN_MODEL_TYPE;
    }
    return MV_UNKNOWN_MODEL_TYPE;
}


static int __mv_new_obs_PCA_T(MVMat *t, MVMat *E, const MVMat *new_X, const MVModel *model,
                  int num_components, MVNewScoreCalcType method)
{
    int a, i, j;
    int freeE=0;
    double numerator, denominator;
    MVMat *p, *EHat;
    MVMat *_t, *_p, *_p_T, *_slice; // slices for calculation of E-hat
    if ( !(t->nrows == new_X->nrows && new_X->ncolumns == model->p->nrows &&
           model->p->ncolumns >= num_components && t->ncolumns >= num_components))
    {
        return MV_INCORRECT_DIMENSIONS;
    }

    if (E)
    {
        if (! (E->nrows == new_X->nrows && E->ncolumns == new_X->ncolumns))
        {
            return MV_INCORRECT_DIMENSIONS;
        }
        mvmat_copy(E, new_X);
    }
    else
    {
        E = mvmat_alloc_copy(new_X);
        freeE = 1;
    }

    // method is not yet used.
    (void)method;

    p = model->p;

    _t = mvmat_alloc(t->nrows, 1);
    _p = mvmat_alloc(p->nrows, 1);
    _p_T = mvmat_alloc(1, p->nrows);
    _slice = mvmat_alloc(1, 1);
    EHat = mvmat_alloc(E->nrows, E->ncolumns);

    for(a=0; a < num_components; a++)
    {
        for (i=0; i < E->nrows; i++)
        {
            numerator = denominator = 0.0;
            for (j=0; j < E->ncolumns; j++)
            {
                if (E->mask[i][j] == DATA_PRESENT)
                {
                    numerator += E->data[i][j] * p->data[j][a];
                    denominator += p->data[j][a] * p->data[j][a];
                }
            }
            if (denominator == 0.0)
            {
                t->data[i][a] = 0.0;
            }
            else
            {
                t->data[i][a] = numerator / denominator;
            }
        }
        // compute residual and start over.
        _slice->data[0][0]=(double)a;
        mvmat_slice_columns(_p, p, _slice);
        mvmat_slice_columns(_t, t, _slice);
        mvmat_transpose(_p_T, _p);
        mvmat_mult(EHat, _t, _p_T);  // XHat = tpT
        mvmat_subtract(E, E, EHat);
    }

    mvmat_free(&EHat);
    mvmat_free(&_slice);
    mvmat_free(&_p_T);
    mvmat_free(&_p);
    mvmat_free(&_t);
    if (freeE)
    {
        mvmat_free(&E);
    }

    return MV_SUCCESS;
}

static int __mv_new_obs_PLS_T(MVMat *t, MVMat *E, const MVMat *new_X, const MVModel *model,
                  int num_components, MVNewScoreCalcType method)
{
    int a, i, j;
    int freeE = 0;
    double numerator, denominator;
    MVMat *p, *w, *EHat;
    MVMat *_t, *_p, *_p_T, *_slice; // slices for calculation of E-hat
    if ( !(t->nrows == new_X->nrows && new_X->ncolumns == model->p->nrows &&
           model->p->ncolumns >= num_components &&
           t->ncolumns == num_components && model->w->ncolumns >= num_components))
    {
        return MV_INCORRECT_DIMENSIONS;
    }

    if (E)
    {
        if (! (E->nrows == new_X->nrows && E->ncolumns == new_X->ncolumns))
        {
            return MV_INCORRECT_DIMENSIONS;
        }
        mvmat_copy(E, new_X);
    }
    else
    {
        E = mvmat_alloc_copy(new_X);
        freeE = 1;
    }

    // method is not yet used
    (void) method;

    p = model->p;
    w = model->w;

    _t = mvmat_alloc(t->nrows, 1);
    _p = mvmat_alloc(p->nrows, 1);
    _p_T = mvmat_alloc(1, p->nrows);
    _slice = mvmat_alloc(1, 1);
    EHat = mvmat_alloc(E->nrows, E->ncolumns);

    for(a=0; a < num_components; a++)
    {
        for (i=0; i < E->nrows; i++)
        {
            numerator = denominator = 0.0;
            for (j=0; j < E->ncolumns; j++)
            {
                if (E->mask[i][j] == DATA_PRESENT)
                {
                    numerator += E->data[i][j] * w->data[j][a];
                    denominator += w->data[j][a] * w->data[j][a];
                }
            }
            if (denominator == 0.0)
            {
                t->data[i][a] = 0.0;
            }
            else
            {
                t->data[i][a] = numerator / denominator;
            }
        }

        // compute residual and start over.
        _slice->data[0][0]=(double)a;
        mvmat_slice_columns(_p, p, _slice);
        mvmat_slice_columns(_t, t, _slice);
        mvmat_transpose(_p_T, _p);
        mvmat_mult(EHat, _t, _p_T);  // XHat = tpT
        mvmat_subtract(E, E, EHat);
    }

    mvmat_free(&EHat);
    mvmat_free(&_slice);
    mvmat_free(&_p_T);
    mvmat_free(&_p);
    mvmat_free(&_t);
    if (freeE)
    {
        mvmat_free(&E);
    }

    return MV_SUCCESS;
}

int mvmodel_u_scores_from_obs(MVMat *u, MVMat *F, const MVMat *new_Y, const MVMat *new_T,
                  const MVModel *model, int num_components,
                  MVNewScoreCalcType method)
{
    int a, i, j;
    int freeF = 0;
    double numerator, denominator;
    MVMat *c, *FHat;
    MVMat *_t, *_c, *_c_T, *_slice; // slices for calculation of E-hat
    if ( !(u->nrows == new_Y->nrows && new_Y->ncolumns == model->c->nrows &&
           model->c->ncolumns >= num_components &&
           u->ncolumns == num_components && u->nrows == new_T->nrows &&
           u->ncolumns == new_T->ncolumns))
    {
        return MV_INCORRECT_DIMENSIONS;
    }
    else if (model->model_type != MV_MODEL_TYPE_PLS)
    {
        return MV_WRONG_MODEL_TYPE;
    }
    if (F)
    {
        if (! (F->nrows == new_Y->nrows && F->ncolumns == new_Y->ncolumns))
        {
            return MV_INCORRECT_DIMENSIONS;
        }
        mvmat_copy(F, new_Y);
    }
    else
    {
        F = mvmat_alloc_copy(new_Y);
        freeF = 1;
    }

    // method is not yet used
    (void) method;
    c = model->c;

    _t = mvmat_alloc(new_T->nrows, 1);
    _c = mvmat_alloc(c->nrows, 1);
    _c_T = mvmat_alloc(1, c->nrows);
    _slice = mvmat_alloc(1, 1);
    FHat = mvmat_alloc(F->nrows, F->ncolumns);

    for(a=0; a < num_components; a++)
    {
        for (i=0; i < F->nrows; i++)
        {
            numerator = denominator = 0.0;
            for (j=0; j < F->ncolumns; j++)
            {
                if (F->mask[i][j] == DATA_PRESENT)
                {
                    numerator += F->data[i][j] * c->data[j][a];
                    denominator += c->data[j][a] * c->data[j][a];
                }
            }
            if (denominator == 0.0)
            {
                u->data[i][a] = 0.0;
            }
            else
            {
                u->data[i][a] = numerator / denominator;
            }
        }
        // compute residual and start over.
        _slice->data[0][0]=(double)a;
        mvmat_slice_columns(_c, c, _slice);
        mvmat_slice_columns(_t, new_T, _slice);
        mvmat_transpose(_c_T, _c);
        mvmat_mult(FHat, _t, _c_T);  // FHat = tcT
        mvmat_subtract(F, F, FHat);
    }

    mvmat_free(&FHat);
    mvmat_free(&_slice);
    mvmat_free(&_c_T);
    mvmat_free(&_c);
    mvmat_free(&_t);
    if (freeF)
    {
        mvmat_free(&F);
    }

    return MV_SUCCESS;
}


int mvmodel_t_scores_from_obs(MVMat *t, MVMat * E, const MVMat *new_X, const MVModel *model,
              int num_components, MVNewScoreCalcType method)
{
    switch (model->model_type)
    {
    case MV_MODEL_TYPE_PCA:
        return __mv_new_obs_PCA_T(t, E, new_X, model, num_components, method);
    case MV_MODEL_TYPE_PLS:
        return __mv_new_obs_PLS_T(t, E, new_X, model, num_components, method);
    default:
        return MV_UNKNOWN_MODEL_TYPE;
    }
    return MV_UNKNOWN_MODEL_TYPE;
}

int mvmodel_autofit(MVModel *model)
{
    int return_val = MV_SUCCESS;
    int component_is_valid = 1;
    model->A = 0;
    while (component_is_valid)
    {
        double iter = 0.0;
        double Q2 = 0.0;
        component_is_valid = 0;
        mvmodel_add_component(model);
        mvmat_get_elem(model->iter, &iter, model->A-1, 0);      
        mvmat_get_elem(model->Q2, &Q2, model->A-1, 0);
        /* Rule 3 */
        if (model->A > MIN(model->X->nrows, model->X->ncolumns))
        {
            return_val = MV_CROSSVAL_RULE3;
        }
        /* RULE 4 */
        else if (iter >= MAX_NIPALS_ITER)
        {
            return_val = MV_CROSSVAL_RULE4;
        }
        /* Rule 1 */
        else if (Q2 >= AUTOFIT_THRESHOLD)
        {
            component_is_valid = 1;
        }
        else
        {
            int i;
            double QV;
            if (model->model_type == MV_MODEL_TYPE_PCA)
            {
                int n_QV_gt = 0;
                for (i = 0; i < model->X->ncolumns; i++)
                {
                    mvmat_get_elem(model->Q2V, &QV, i, model->A-1);
                    if (QV >= AUTOFIT_THRESHOLD)
                    {
                        n_QV_gt++;
                    }
                }
                if (sqrt(n_QV_gt) > model->X->ncolumns)
                {
                    component_is_valid = 1;
                }
            }
            else if (model->model_type == MV_MODEL_TYPE_PLS)
            {
                for (i = 0; i < model->Y->ncolumns; i++)
                {
                    mvmat_get_elem(model->Q2V, &QV, i, model->A-1);
                    if (QV >= AUTOFIT_THRESHOLD)
                    {
                        component_is_valid = 1;
                        break;
                    }
                }
            }
        }

    }
    /* Subtract one at the end because we've left the loop because the last
       component was invalid / insignificant */
    model->A--;
    return MV_SUCCESS;
}


static int __mv_compute_pred(MVMat *pred, const MVMat *scores, const MVMat *weights, int num_components)
{
    // todo: error checking
    MVReturnCode ret;
    MVMat _weightsT, _scores;
    MVMat *weightsT = mvmat_alloc(weights->ncolumns, weights->nrows);
    mvmat_transpose(weightsT, weights);
    _weightsT.nrows = num_components;
    _weightsT.ncolumns = weightsT->ncolumns;
    _weightsT.data = weightsT->data;
    _weightsT.mask = weightsT->mask;
    _weightsT.isReference = 1;

    _scores.nrows = scores->nrows;
    _scores.ncolumns = num_components;
    _scores.data = scores->data;
    _scores.mask = scores->mask;
    _scores.isReference = 1;

    ret = mvmat_mult(pred, &_scores, &_weightsT);

    mvmat_free(&weightsT);
    return ret;
}


int mvmodel_compute_xpred(MVMat *Xhat, const MVModel *model, const MVMat *t, int num_components)
{
    MVMat *weights = NULL;
    if (model->model_type == MV_MODEL_TYPE_PCA)
    {
        weights = model->p;
    }
    else if(model->model_type == MV_MODEL_TYPE_PLS)
    {
        weights = model->wstar;
    }
    return __mv_compute_pred(Xhat, t, weights, num_components);
}


int mvmodel_compute_ypred(MVMat *Yhat, const MVModel *model, const MVMat *t, int num_components)
{
    // todo: error checking
    return __mv_compute_pred(Yhat, t, model->c, num_components);
}
