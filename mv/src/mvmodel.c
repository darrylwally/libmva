/*! \file mvalgorithms.h
  \brief Contains structures and functions for initializing multivariate
  structures and performing pca and pls.
  \author Darryl Wallace - wallacdj@gmail.com, darryl@wallynet.ca
  \date 2012

  THIS CODE IS COPYRIGHT DARRYL WALLACE (c) 2012.
  UNAUTHORIZED DISTRIBUTION, USE, OR VIEWING THE CONTENTS OF THIS FILE
  IS STRICTLY PROHIBITED.
  IF YOU HAVE RECEIVED THIS CODE WITHOUT AUTHORIZATION IT MUST BE DESTROYED
  IMMEDIATELY
  */

#include "mvmodel.h"
#include "mvconstants.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#define MAX_NIPALS_ITER 500

// Function prototypes
static int __mvAddPCAComponent(mvModel *model, int performCrossValidation);
static int __mvAddPLSComponent(mvModel *model, int performCrossValidation);

static int __freeCrossValData(crossValData **cvData)
{
    int i;
    crossValData * cvd = NULL;
    if (cvData == NULL)
        return -1;
    cvd = *cvData;
    if (cvd == NULL)
        return -1;

    mvFreeMat (&cvd->PRESS);
    mvFreeMat (&cvd->PRESSV);
    for (i=0; i< cvd->numRounds; i++)
    {
        mvFreeModel((mvModel **)&cvd->models[i]);
    }
    free(cvd->models);
    *cvData = NULL;
    return SUCCESS;
}

static crossValData * __allocCrossValData(int numRounds)
{
    int i;
    crossValData *output = (crossValData*) malloc(sizeof(crossValData));
    if (!output)
        return NULL;
    output->models = malloc(numRounds * sizeof(mvModel*));
    for (i=0; i< numRounds; i++)
    {
        output->models[i] = NULL;
    }
    output->PRESSV = NULL;
    output->PRESS = NULL;
    output->numRounds = numRounds;
    return output;
}

/*! Performs regression of a vector on a matrix.

  When the regression vector has the same number of columns (Kx1) as the
  matrix (NxK) and the result has the same number of rows as the matrix (Nx1)

  E.g., in the NIPALS algorithm t = X p / p'p

  */
static int __mvRegressCol(mvMat *output, const mvMat *X, const mvMat *vector)
{

    /* TODO size checking.  Right now I'm avoiding it because it's an internal
       function that I think should have checked sizes already.. */
    double oNum, oDen;
    int i,j;

    for (i=0; i<X->nrows; i++)
    {
        oNum=oDen=0.0;
        for (j=0; j<X->ncolumns; j++)
        {
            if (X->mask[i][j] == DATA_PRESENT)
            {
                oNum += vector->data[j][0] * X->data[i][j];
                oDen += vector->data[j][0] * vector->data[j][0];
            }
        }
        if (oDen==0.0)  //case where entire column is missing
        {
            output->data[i][0]=0.0;
        }
        else
        {
            output->data[i][0]=oNum/oDen;
        }
    }

    return 0;
}

/*! Performs regression of a vector on a matrix.

  When the regression vector has the same number of rows (Nx1) as the
  matrix (NxK) and the result has the same number of columns as the matrix (Kx1)

  E.g., in the NIPALS algorithm p = t' X / t't

  */
static int __mvRegressRow(mvMat *output, const mvMat *X, const mvMat *vector)
{
    /* TODO size checking.  Right now I'm avoiding it because it's an internal
       function that I think should have checked sizes already.. */
    double oNum, oDen;
    int i,j;
    for (j=0;j<X->ncolumns;j++)
    {
        oNum=oDen=0.0;

        for (i=0; i<X->nrows;i++)
        {
            if (X->mask[i][j] == DATA_PRESENT)
            {
                oNum += vector->data[i][0] * X->data[i][j];
                oDen += vector->data[i][0] * vector->data[i][0];
            }
        }
        if (oDen==0.0)  //case where entire column is missing
        {
            output->data[j][0]=0.0;
        }
        else
        {
            output->data[j][0]=oNum/oDen;
        }
    }
    return 0;
}

/*! \internal WStar will be computed for the a'th component
  This is based on equation equation 30 in Dayal, B.S. & MacGregor, J.F. (1997)
  Improved PLS algorithms. J. Chemometrics 11:73-85
  \arg wStarOut is a pre-allocated vector of dim Kx1.
*/
static int __computeWStar(mvMat *wStarOut, const mvMat * w, const mvMat *p,
                          const mvMat *wStarSoFar)
{
    int i,j,a;
    int K = w->nrows;
    mvMat * wa = mvAllocMat(K, 1); // Ath column of W (permanently)
    mvMat *p_i = mvAllocMat(K, 1);
    mvMat *wStarSoFar_i = mvAllocMat(K,1);
    // start wStarOut as the the most recent vector of w
    a = wStarSoFar->ncolumns;
    for (i=0; i< w->nrows; i++)
    {
        double val;
        mvMatGetElem(w, &val, i, a);
        mvMatSetElem(wStarOut, i, 0, val);
        mvMatSetElem(wa, i, 0, val);
    }

    for (i=0; i<a; i++)
    {
        // extract the ith column of p and wStarSoFar
        double pTwa;
        for (j=0; j<K; j++)
        {
            double val;
            mvMatGetElem(p, &val, j, i);
            mvMatSetElem(p_i, j, 0, val);
            mvMatGetElem(wStarSoFar, &val, j, i);
            mvMatSetElem(wStarSoFar_i, j, 0, val);
        }
        pTwa = mvDotProduct(p_i, wa);
        mvMatMultS(wStarSoFar_i, wStarSoFar_i, pTwa);
        mvSubtractMat(wStarOut, wStarOut, wStarSoFar_i);
    }

    mvFreeMat(&wa);
    mvFreeMat(&p_i);
    mvFreeMat(&wStarSoFar_i);

    return SUCCESS;
}

static int __crossValidatePCA_FAST(mvModel *pca)
{
    int round;
    int i;
    crossValData *cv = pca->cvd;
    int numRounds = cv->numRounds;
    mvMat * PRESS = mvAllocMatZ(1,1);
    mvMat * PRESSV = mvAllocMatZ(1, pca->X->ncolumns);
    mvMat * PRESSV_temp = mvAllocMatZ(1, pca->X->ncolumns); // acts as PRESSV temporary addition location and (PRESSV / SSV)
    mvMat * componentSlice = mvAllocMat(1,1);
    double Q2;
    mvMat * Q2V = mvAllocMat(1, pca->X->ncolumns);
    mvMat * X = NULL;

    // slice the correct matrix.
    if (pca->_A == 1)
    {
        X = pca->X;
    }
    else // pca->_A > 1x
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
        mvModel *modelRound;
        mvMat * X_model = NULL; // the matrix that the model is built
        mvMat * X_pred = NULL;  // the matrix of the data that is to be predicted
        mvMat * E_pred = NULL;  // will function as the predicted values and the residuals.
        mvMat * t_pred = NULL;
        mvMat * slice_pred = mvRange(round, pca->E->nrows, numRounds);

        mvMatSliceRowsRef(&X_pred, X, slice_pred);
        mvMatDeleteRowsRef(&X_model, X, slice_pred);

        E_pred = mvAllocMat(X_pred->nrows, X_pred->ncolumns);
        t_pred = mvAllocMat(X_pred->nrows, 1);

        // Initialize model, add one component, compute new observation scores.
        modelRound = mvInitPCAModel(X_model);
        __mvAddPCAComponent(modelRound, 0);
        mvNewObsT(t_pred, E_pred, X_pred, modelRound, 1, SCP);

        // Sum the columns of E_pred and then add those to PRESSV
        mvMatColumnSS(PRESSV_temp, E_pred);
        mvAddMat(PRESSV, PRESSV, PRESSV_temp);

        mvFreeMat(&slice_pred);
        mvFreeMat(&t_pred);
        mvFreeMat(&E_pred);
        mvFreeMat(&X_pred);
        mvFreeMat(&X_model);
        mvFreeModel(&modelRound);
    }
    mvMatSetElem(PRESS, 0, 0, mvMatSum(PRESSV));
    // Compute Q2 and Q2V
    {
        mvMat * SSXV_ref = NULL;
        componentSlice->data[0][0] = pca->_A - 1;
        Q2 = 1.0 - PRESS->data[0][0] / pca->SSX->data[pca->_A-1][0];
        mvMatSliceRowsRef(&SSXV_ref, pca->SSXV, componentSlice);
        mvMatColumnDiv(PRESSV_temp, PRESSV, SSXV_ref);
        mvMatMultS(PRESSV_temp, PRESSV, -1.0);
        mvAddMatS(Q2V, PRESSV_temp, 1.0);
        mvFreeMat(&SSXV_ref);
    }

    if (pca->_A == 1)
    {
        // Q2's have not yet been allocated.
        pca->Q2 = mvAllocMat(1,1);
        pca->Q2cum = mvAllocMat(1,1);
        pca->Q2V = Q2V;
        pca->Q2Vcum = mvAllocMat(1,X->ncolumns);

        cv->PRESS = PRESS;
        cv->PRESSV = PRESSV;

        // for the first component, Q2[V] = Q2[V]cum
        mvMatCopy(pca->Q2Vcum, pca->Q2V);
        mvMatSetElem(pca->Q2, 0, 0, Q2);
        mvMatSetElem(pca->Q2cum, 0, 0, Q2);
    }
    else
    {
        double Q2cum=1.0;
        mvMat *newPRESS, *newPRESSV, *PRESSV_SSV;
        mvMat *newQ2, *newQ2V, *newQ2cum, *newQ2Vcum, *newQ2Vcum_ref;
        newPRESS = mvAllocMat(pca->_A, 1);
        newPRESSV = mvAllocMat(pca->_A, pca->X->ncolumns);
        mvConcatRows(newPRESS, cv->PRESS, PRESS);
        mvConcatRows(newPRESSV, cv->PRESSV, PRESSV);
        mvFreeMat(&cv->PRESS);
        mvFreeMat(&cv->PRESSV);
        mvFreeMat(&PRESS);
        mvFreeMat(&PRESSV);
        cv->PRESS = newPRESS;
        cv->PRESSV = newPRESSV;

        // new Q2
        newQ2 = mvAllocMat(pca->_A, 1);
        for (i=0; i<pca->_A-1; i++)
        {
            newQ2->data[i][0] = pca->Q2->data[i][0];
        }
        newQ2->data[pca->_A-1][0] = Q2;
        mvFreeMat(&pca->Q2);
        pca->Q2 = newQ2;

        // new Q2V
        newQ2V = mvAllocMat(pca->_A, X->ncolumns);
        mvConcatRows(newQ2V, pca->Q2V, Q2V);
        mvFreeMat(&pca->Q2V);
        mvFreeMat(&Q2V);
        pca->Q2V = newQ2V;

        // new Q2cum
        newQ2cum = mvAllocMat(pca->_A, 1);
        for(i=0; i<pca->_A-1; i++)
        {
            newQ2cum->data[i][0] = pca->Q2cum->data[i][0];
            Q2cum *= cv->PRESS->data[i][0] / pca->SSX->data[i][0];
        }
        Q2cum *= cv->PRESS->data[pca->_A-1][0] / pca->SSX->data[pca->_A-1][0];
        newQ2cum->data[pca->_A-1][0] = 1.0 - Q2cum;
        mvFreeMat(&pca->Q2cum);
        pca->Q2cum = newQ2cum;

        // new Q2Vcum
        newQ2Vcum = mvAllocMat(pca->_A, X->ncolumns);
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
        mvMatSet(PRESSV_temp, 1.0);
        PRESSV_SSV = mvAllocMat(1, X->ncolumns);
        for (i=0; i < pca->_A; i++)
        {
            mvMat *PRESSV_ref, *SSXV_ref;
            PRESSV_ref = NULL;
            SSXV_ref = NULL;
            componentSlice->data[0][0] = i;
            mvMatSliceRowsRef(&PRESSV_ref, cv->PRESSV, componentSlice);
            mvMatSliceRowsRef(&SSXV_ref, pca->SSXV, componentSlice);
            mvMatColumnDiv(PRESSV_SSV, PRESSV_ref, SSXV_ref);
            mvMatElemMult(PRESSV_temp, PRESSV_temp, PRESSV_SSV);

            mvFreeMat(&PRESSV_ref);
            mvFreeMat(&SSXV_ref);
        }
        componentSlice->data[0][0] = pca->_A - 1;
        newQ2Vcum_ref = NULL;
        mvMatSliceRowsRef(&newQ2Vcum_ref, newQ2Vcum, componentSlice);
        mvMatMultS(PRESSV_temp, PRESSV_SSV, -1.0);
        mvAddMatS(newQ2Vcum_ref, PRESSV_temp, 1.0);
        mvFreeMat(&PRESSV_SSV);
        mvFreeMat(&newQ2Vcum_ref);
        mvFreeMat(&pca->Q2Vcum);
        pca->Q2Vcum = newQ2Vcum;
    }
    mvFreeMat(&componentSlice);
    mvFreeMat(&PRESSV_temp);
    return 0;
}

static int __crossValidatePCA_FULL(mvModel *pca)
{
    int round;
    int i;
    crossValData *cv = pca->cvd;
    int numRounds = cv->numRounds;
    mvMat * PRESS = mvAllocMatZ(1,1);
    mvMat * PRESSV = mvAllocMatZ(1, pca->X->ncolumns);
    mvMat * PRESSV_temp = mvAllocMatZ(1, pca->X->ncolumns); // acts as PRESSV temporary addition location and (PRESSV / SSV)
    mvMat * componentSlice = mvAllocMat(1,1);
    double Q2;
    mvMat * Q2V = mvAllocMat(1, pca->X->ncolumns);
    mvMat * X = NULL;

    // slice the correct matrix.
    if (pca->_A == 1)
    {
        X = pca->X;
    }
    else // pca->_A > 1x
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
        mvModel *modelRound;
        mvMat * X_model = NULL; // the matrix that the model is built
        mvMat * X_pred = NULL;  // the matrix of the data that is to be predicted
        mvMat * E_pred = NULL;  // will function as the predicted values and the residuals.
        mvMat * t_pred = NULL;
        mvMat * slice_pred = mvRange(round, pca->E->nrows, numRounds);

        mvMatSliceRowsRef(&X_pred, X, slice_pred);
        mvMatDeleteRowsRef(&X_model, X, slice_pred);

        E_pred = mvAllocMat(X_pred->nrows, X_pred->ncolumns);
        t_pred = mvAllocMat(X_pred->nrows, 1);

        // Initialize model, add one component, compute new observation scores.
        modelRound = mvInitPCAModel(X_model);
        __mvAddPCAComponent(modelRound, 0);
        mvNewObsT(t_pred, E_pred, X_pred, modelRound, 1, SCP);

        // Sum the columns of E_pred and then add those to PRESSV
        mvMatColumnSS(PRESSV_temp, E_pred);
        mvAddMat(PRESSV, PRESSV, PRESSV_temp);

        mvFreeMat(&slice_pred);
        mvFreeMat(&t_pred);
        mvFreeMat(&E_pred);
        mvFreeMat(&X_pred);
        mvFreeMat(&X_model);
        mvFreeModel(&modelRound);
    }
    mvMatSetElem(PRESS, 0, 0, mvMatSum(PRESSV));
    // Compute Q2 and Q2V
    {
        mvMat * SSXV_ref = NULL;
        componentSlice->data[0][0] = pca->_A - 1;
        Q2 = 1.0 - PRESS->data[0][0] / pca->SSX->data[pca->_A-1][0];
        mvMatSliceRowsRef(&SSXV_ref, pca->SSXV, componentSlice);
        mvMatColumnDiv(PRESSV_temp, PRESSV, SSXV_ref);
        mvMatMultS(PRESSV_temp, PRESSV, -1.0);
        mvAddMatS(Q2V, PRESSV_temp, 1.0);
        mvFreeMat(&SSXV_ref);
    }

    if (pca->_A == 1)
    {
        // Q2's have not yet been allocated.
        pca->Q2 = mvAllocMat(1,1);
        pca->Q2cum = mvAllocMat(1,1);
        pca->Q2V = Q2V;
        pca->Q2Vcum = mvAllocMat(1,X->ncolumns);

        cv->PRESS = PRESS;
        cv->PRESSV = PRESSV;

        // for the first component, Q2[V] = Q2[V]cum
        mvMatCopy(pca->Q2Vcum, pca->Q2V);
        mvMatSetElem(pca->Q2, 0, 0, Q2);
        mvMatSetElem(pca->Q2cum, 0, 0, Q2);
    }
    else
    {
        double Q2cum=1.0;
        mvMat *newPRESS, *newPRESSV, *PRESSV_SSV;
        mvMat *newQ2, *newQ2V, *newQ2cum, *newQ2Vcum, *newQ2Vcum_ref;
        newPRESS = mvAllocMat(pca->_A, 1);
        newPRESSV = mvAllocMat(pca->_A, pca->X->ncolumns);
        mvConcatRows(newPRESS, cv->PRESS, PRESS);
        mvConcatRows(newPRESSV, cv->PRESSV, PRESSV);
        mvFreeMat(&cv->PRESS);
        mvFreeMat(&cv->PRESSV);
        mvFreeMat(&PRESS);
        mvFreeMat(&PRESSV);
        cv->PRESS = newPRESS;
        cv->PRESSV = newPRESSV;

        // new Q2
        newQ2 = mvAllocMat(pca->_A, 1);
        for (i=0; i<pca->_A-1; i++)
        {
            newQ2->data[i][0] = pca->Q2->data[i][0];
        }
        newQ2->data[pca->_A-1][0] = Q2;
        mvFreeMat(&pca->Q2);
        pca->Q2 = newQ2;

        // new Q2V
        newQ2V = mvAllocMat(pca->_A, X->ncolumns);
        mvConcatRows(newQ2V, pca->Q2V, Q2V);
        mvFreeMat(&pca->Q2V);
        mvFreeMat(&Q2V);
        pca->Q2V = newQ2V;

        // new Q2cum
        newQ2cum = mvAllocMat(pca->_A, 1);
        for(i=0; i<pca->_A-1; i++)
        {
            newQ2cum->data[i][0] = pca->Q2cum->data[i][0];
            Q2cum *= cv->PRESS->data[i][0] / pca->SSX->data[i][0];
        }
        Q2cum *= cv->PRESS->data[pca->_A-1][0] / pca->SSX->data[pca->_A-1][0];
        newQ2cum->data[pca->_A-1][0] = 1.0 - Q2cum;
        mvFreeMat(&pca->Q2cum);
        pca->Q2cum = newQ2cum;

        // new Q2Vcum
        newQ2Vcum = mvAllocMat(pca->_A, X->ncolumns);
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
        mvMatSet(PRESSV_temp, 1.0);
        PRESSV_SSV = mvAllocMat(1, X->ncolumns);
        for (i=0; i < pca->_A; i++)
        {
            mvMat *PRESSV_ref, *SSXV_ref;
            PRESSV_ref = NULL;
            SSXV_ref = NULL;
            componentSlice->data[0][0] = i;
            mvMatSliceRowsRef(&PRESSV_ref, cv->PRESSV, componentSlice);
            mvMatSliceRowsRef(&SSXV_ref, pca->SSXV, componentSlice);
            mvMatColumnDiv(PRESSV_SSV, PRESSV_ref, SSXV_ref);
            mvMatElemMult(PRESSV_temp, PRESSV_temp, PRESSV_SSV);

            mvFreeMat(&PRESSV_ref);
            mvFreeMat(&SSXV_ref);
        }
        componentSlice->data[0][0] = pca->_A - 1;
        newQ2Vcum_ref = NULL;
        mvMatSliceRowsRef(&newQ2Vcum_ref, newQ2Vcum, componentSlice);
        mvMatMultS(PRESSV_temp, PRESSV_SSV, -1.0);
        mvAddMatS(newQ2Vcum_ref, PRESSV_temp, 1.0);
        mvFreeMat(&PRESSV_SSV);
        mvFreeMat(&newQ2Vcum_ref);
        mvFreeMat(&pca->Q2Vcum);
        pca->Q2Vcum = newQ2Vcum;
    }
    mvFreeMat(&componentSlice);
    mvFreeMat(&PRESSV_temp);
    return 0;
}


static int __crossValidatePCA(mvModel *pca)
{
    switch (pca->crossValType)
    {
    case FAST:
        return __crossValidatePCA_FAST(pca);
    case FULL:
        return __crossValidatePCA_FULL(pca);
    default:
        return UNKNOWN_MODEL_TYPE;
    }

    return 0;
}

static int __crossValidatePLS_FAST(mvModel *pls)
{

    int round;
    int i;
    crossValData *cv = pls->cvd;
    int numRounds = cv->numRounds;
    mvMat * PRESS = mvAllocMatZ(1,1);
    mvMat * PRESSV = mvAllocMatZ(1, pls->Y->ncolumns);
    mvMat * PRESSV_temp = mvAllocMatZ(1, pls->Y->ncolumns); // acts as PRESSV temporary addition location and (PRESSV / SSV)
    mvMat * componentSlice = mvAllocMat(1,1);
    double Q2;
    mvMat * Q2V = mvAllocMat(1, pls->Y->ncolumns);
    mvMat * X = NULL;
    mvMat * Y = NULL;

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
    if (numRounds <= 0)
    {
        numRounds = X->nrows;
    }

    // allocate the models with the sliced rows of the residual matrix.
    for (round = 0; round < numRounds; round++)
    {
        mvModel *modelRound;
        mvMat * Y_model = NULL; // the matrix that the model is built
        mvMat * Y_pred = NULL;  // the matrix of the data that is to be predicted
        mvMat * X_model = NULL; // the matrix that the model is built
        mvMat * X_pred = NULL;  // the matrix of the data that is to be predicted
        mvMat * F_pred = NULL;  // will function as the predicted values and the residuals.
        mvMat * t_pred = NULL;
        mvMat * u_pred = NULL;
        mvMat * slice_pred = mvRange(round, pls->F->nrows, numRounds);

        mvMatSliceRowsRef(&Y_pred, Y, slice_pred);
        mvMatDeleteRowsRef(&Y_model, Y, slice_pred);
        mvMatSliceRowsRef(&X_pred, X, slice_pred);
        mvMatDeleteRowsRef(&X_model, X, slice_pred);

        F_pred = mvAllocMat(Y_pred->nrows, Y_pred->ncolumns);
        t_pred = mvAllocMat(Y_pred->nrows, 1);
        u_pred = mvAllocMat(Y_pred->nrows, 1);

        // Initialize model, add one component, compute new observation scores.
        modelRound = mvInitPLSModel(X_model, Y_model);
        __mvAddPLSComponent(modelRound, 0);
        mvNewObsT(t_pred, NULL, X_pred, modelRound, 1, SCP);
        mvNewObsU(u_pred, F_pred, Y_pred, t_pred, modelRound, 1, SCP);

        // Sum the columns of F_pred and then add those to PRESSV
        mvMatColumnSS(PRESSV_temp, F_pred);
        mvAddMat(PRESSV, PRESSV, PRESSV_temp);


        mvFreeMat(&slice_pred);
        mvFreeMat(&u_pred);
        mvFreeMat(&t_pred);
        mvFreeMat(&F_pred);
        mvFreeMat(&X_pred);
        mvFreeMat(&X_model);
        mvFreeMat(&Y_pred);
        mvFreeMat(&Y_model);
        mvFreeModel(&modelRound);
    }
    mvMatSetElem(PRESS, 0, 0, mvMatSum(PRESSV));
    // Compute Q2 and Q2V
    {
        mvMat * SSYV_ref = NULL;
        componentSlice->data[0][0] = pls->_A - 1;
        Q2 = 1.0 - PRESS->data[0][0] / pls->SSY->data[pls->_A-1][0];
        mvMatSliceRowsRef(&SSYV_ref, pls->SSYV, componentSlice);
        mvMatColumnDiv(PRESSV_temp, PRESSV, SSYV_ref);
        mvMatMultS(PRESSV_temp, PRESSV, -1.0);
        mvAddMatS(Q2V, PRESSV_temp, 1.0);
        mvFreeMat(&SSYV_ref);
    }

    if (pls->_A == 1)
    {
        // Q2's have not yet been allocated.
        pls->Q2 = mvAllocMat(1,1);
        pls->Q2cum = mvAllocMat(1,1);
        pls->Q2V = Q2V;
        pls->Q2Vcum = mvAllocMat(1, Y->ncolumns);

        cv->PRESS = PRESS;
        cv->PRESSV = PRESSV;

        // for the first component, Q2[V] = Q2[V]cum
        mvMatCopy(pls->Q2Vcum, pls->Q2V);
        mvMatSetElem(pls->Q2, 0, 0, Q2);
        mvMatSetElem(pls->Q2cum, 0, 0, Q2);
    }
    else
    {
        double Q2cum=1.0;
        mvMat *newPRESS, *newPRESSV, *PRESSV_SSV;
        mvMat *newQ2, *newQ2V, *newQ2cum, *newQ2Vcum, *newQ2Vcum_ref;
        newPRESS = mvAllocMat(pls->_A, 1);
        newPRESSV = mvAllocMat(pls->_A, pls->Y->ncolumns);
        mvConcatRows(newPRESS, cv->PRESS, PRESS);
        mvConcatRows(newPRESSV, cv->PRESSV, PRESSV);
        mvFreeMat(&cv->PRESS);
        mvFreeMat(&cv->PRESSV);
        mvFreeMat(&PRESS);
        mvFreeMat(&PRESSV);
        cv->PRESS = newPRESS;
        cv->PRESSV = newPRESSV;

        // new Q2
        newQ2 = mvAllocMat(pls->_A, 1);
        for (i=0; i<pls->_A-1; i++)
        {
            newQ2->data[i][0] = pls->Q2->data[i][0];
        }
        newQ2->data[pls->_A-1][0] = Q2;
        mvFreeMat(&pls->Q2);
        pls->Q2 = newQ2;

        // new Q2V
        newQ2V = mvAllocMat(pls->_A, Y->ncolumns);
        mvConcatRows(newQ2V, pls->Q2V, Q2V);
        mvFreeMat(&pls->Q2V);
        mvFreeMat(&Q2V);
        pls->Q2V = newQ2V;

        // new Q2cum
        newQ2cum = mvAllocMat(pls->_A, 1);
        for(i=0; i<pls->_A-1; i++)
        {
            newQ2cum->data[i][0] = pls->Q2cum->data[i][0];
            Q2cum *= cv->PRESS->data[i][0] / pls->SSY->data[i][0];
        }
        Q2cum *= cv->PRESS->data[pls->_A-1][0] / pls->SSY->data[pls->_A-1][0];
        newQ2cum->data[pls->_A-1][0] = 1.0 - Q2cum;
        mvFreeMat(&pls->Q2cum);
        pls->Q2cum = newQ2cum;

        // new Q2Vcum
        newQ2Vcum = mvAllocMat(pls->_A, Y->ncolumns);
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
        mvMatSet(PRESSV_temp, 1.0);
        PRESSV_SSV = mvAllocMat(1, Y->ncolumns);
        for (i=0; i < pls->_A; i++)
        {
            mvMat *PRESSV_ref, *SSYV_ref;
            PRESSV_ref = NULL;
            SSYV_ref = NULL;
            componentSlice->data[0][0] = i;
            mvMatSliceRowsRef(&PRESSV_ref, cv->PRESSV, componentSlice);
            mvMatSliceRowsRef(&SSYV_ref, pls->SSYV, componentSlice);
            mvMatColumnDiv(PRESSV_SSV, PRESSV_ref, SSYV_ref);
            mvMatElemMult(PRESSV_temp, PRESSV_temp, PRESSV_SSV);

            mvFreeMat(&PRESSV_ref);
            mvFreeMat(&SSYV_ref);
        }
        componentSlice->data[0][0] = pls->_A - 1;
        newQ2Vcum_ref = NULL;
        mvMatSliceRowsRef(&newQ2Vcum_ref, newQ2Vcum, componentSlice);
        mvMatMultS(PRESSV_temp, PRESSV_SSV, -1.0);
        mvAddMatS(newQ2Vcum_ref, PRESSV_temp, 1.0);
        mvFreeMat(&PRESSV_SSV);
        mvFreeMat(&newQ2Vcum_ref);
        mvFreeMat(&pls->Q2Vcum);
        pls->Q2Vcum = newQ2Vcum;
    }
    mvFreeMat(&componentSlice);
    mvFreeMat(&PRESSV_temp);
    return 0;
}

static int __crossValidatePLS(mvModel *pls)
{
    switch (pls->crossValType)
    {
    case FAST:
        return __crossValidatePLS_FAST(pls);
//    case FULL:
//        return __crossValidatePLS_FULL(pls);
    default:
        return UNKNOWN_MODEL_TYPE;
    }

    return 0;
}

mvModel * mvInitPCAModel(mvMat *X)
{
    mvModel * output = (mvModel *) malloc(sizeof(mvModel));
    if (!output)
    {
        return NULL;
    }
    output->modelType = PCA;
    output->X = X;
    output->E = mvAllocMatZ(X->nrows, X->ncolumns);
    output->p = NULL;
    output->t = NULL;
    output->R2X = NULL;
    output->cvd = NULL;
    output->crossValType = FAST;
    output->numCrossValRounds = 7;
    output->A = output->_A = 0;
    output->SSX = mvAllocMat(1,1);
    output->SSXV = mvAllocMat(1, X->ncolumns);
    mvMatColumnSS(output->SSXV, X);
    output->SSX->data[0][0] = mvMatSum(output->SSXV);
    output->Q2 = NULL;
    output->Q2V = NULL;
    output->Q2cum = NULL;
    output->Q2Vcum = NULL;
    /* Set PLS options to NULL just because */
    output->SSY = NULL;
    output->SSYV = NULL;
    output->Y = NULL;
    output->F = NULL;
    output->w = NULL;
    output->u = NULL;
    output->wStar = NULL;
    output->c = NULL;
    output->R2Y =NULL;
    /* TODO ADD Cross validation stuff */
    return output;
}

mvModel * mvInitPLSModel(mvMat *X, mvMat *Y)
{
    mvModel *output;
    if (X->nrows != Y->nrows)
    {
        return NULL;
    }
    output = (mvModel *) malloc(sizeof(mvModel));
    if (!output)
    {
        return NULL;
    }
    output->modelType = PLS;
    output->X = X;
    output->Y = Y;
    output->E = mvAllocMatZ(X->nrows, X->ncolumns);
    output->F = mvAllocMatZ(Y->nrows, Y->ncolumns);
    output->u = NULL;
    output->w = NULL;
    output->wStar = NULL;
    output->p = NULL;
    output->t = NULL;
    output->cvd = NULL;
    output->crossValType = FAST;
    output->numCrossValRounds = 7;
    output->SSX = mvAllocMat(1,1);
    output->SSXV = mvAllocMat(1, X->ncolumns);
    mvMatColumnSS(output->SSXV, X);
    output->SSX->data[0][0] = mvMatSum(output->SSXV);
    output->SSY = mvAllocMat(1,1);
    output->SSYV = mvAllocMat(1, Y->ncolumns);
    mvMatColumnSS(output->SSYV, Y);
    output->SSY->data[0][0] = mvMatSum(output->SSYV);
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

static int __mvFreePCAModel(mvModel **model)
{
    mvModel *m = NULL;
    if (!model)
        return -1;
    m = *model;
    if (!m)
        return -1;
    mvFreeMat(&m->E);
    mvFreeMat(&m->t);
    mvFreeMat(&m->p);
    mvFreeMat(&m->R2X);
    mvFreeMat(&m->SSX);
    mvFreeMat(&m->SSXV);
    mvFreeMat(&m->R2X);
    mvFreeMat(&m->Q2);
    mvFreeMat(&m->Q2V);
    mvFreeMat(&m->Q2cum);
    mvFreeMat(&m->Q2Vcum);
    __freeCrossValData((crossValData**) &m->cvd);
    free(m);
    *model = NULL;
    return SUCCESS;
}

static int __mvFreePLSModel(mvModel **model)
{
    mvModel *m = NULL;
    if (!model)
        return -1;
    m = *model;
    if (!m)
        return -1;
    mvFreeMat(&m->E);
    mvFreeMat(&m->F);
    mvFreeMat(&m->u);
    mvFreeMat(&m->w);
    mvFreeMat(&m->wStar);
    mvFreeMat(&m->t);
    mvFreeMat(&m->p);
    mvFreeMat(&m->R2X);
    mvFreeMat(&m->R2Y);
    mvFreeMat(&m->SSX);
    mvFreeMat(&m->SSXV);
    mvFreeMat(&m->SSY);
    mvFreeMat(&m->SSYV);
    mvFreeMat(&m->Q2);
    mvFreeMat(&m->Q2V);
    mvFreeMat(&m->Q2cum);
    mvFreeMat(&m->Q2Vcum);
    __freeCrossValData((crossValData **) &m->cvd);
    free(m);
    *model = NULL;
    return SUCCESS;
}

int mvFreeModel(mvModel **model)
{
    mvModel *m;
    m = *model;
    if (!m)
    {
        return -1;
    }
    switch (m->modelType)
    {
    case PCA:
        return __mvFreePCAModel(model);
    case PLS:
        return __mvFreePLSModel(model);
    default:
        return UNKNOWN_MODEL_TYPE;
    }
    return UNKNOWN_MODEL_TYPE;
}

static int __mvAddPCAComponent(mvModel *model, int performCrossValidation)
{
    int iter = 0;
    int N = model->X->nrows;
    int K = model->X->ncolumns;
    mvMat *p = mvAllocMatVal(K, 1, 1.0);
    mvMat *pOld = mvAllocMat(K, 1);
    mvMat *pDiff = mvAllocMat(K, 1);
    mvMat *t = mvAllocMat(N, 1);
    mvMat *R2 = mvAllocMat(1,1);
    mvMat *X;   // reference -> no clean up req'd.
    mvMat *pT, *tpT;
    double vectorNorm;
    double SSE;
    mvMat *SSEV = mvAllocMat(1, K);

    // Initialize p as a unit length vector.
    mvMatMultS(p, p, mvVectorNorm(p));

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
        iter++;
        __mvRegressCol(t, X, p);    // t= Xp / (p'p);

        mvMatCopy(pOld, p);         // store contents of p into pOld

        __mvRegressRow(p, X, t);    // p = t'X / (t't);

        vectorNorm = mvVectorNorm(p);

        mvMatMultS(p, p, 1.0/vectorNorm);  // normalize p

        mvSubtractMat(pDiff, p, pOld);          // get the difference of P

        vectorNorm = mvVectorNorm(pDiff);
    } while (vectorNorm > MV_DBL_SQRT_EPS && iter < MAX_NIPALS_ITER);

    // no longer need, pOld or pDiff
    mvFreeMat(&pOld);
    mvFreeMat(&pDiff);

    // Increment number of internal components
    model->_A++;

    // XXX: Cross validation must be performed before the new residual E is
    // computed and depends on _A == 1 for the first component (NOT ZERO!)
    if (performCrossValidation)
    {
        if (!model->cvd)
        {
            model->cvd = __allocCrossValData(model->numCrossValRounds);
        }
        __crossValidatePCA(model);
    }

    // compute residual
    tpT = mvAllocMat(X->nrows, X->ncolumns);
    pT = mvAllocMat(1, X->ncolumns);
    mvTransposeMat(pT, p);
    mvMatMult(tpT, t, pT);
    mvSubtractMat(model->E, X, tpT);
    mvFreeMat(&pT);
    mvFreeMat(&tpT);

    // Compute R2;
    mvMatColumnSS(SSEV, model->E);
    SSE = mvMatSS(model->E);
    mvMatSetElem(R2, 0, 0, 1.0 - SSE/model->SSX->data[0][0]);

    // store new components.
    if (model->_A > 1)
    {
        mvMat *newR2 = mvAllocMat(model->_A, 1);
        mvMat *newT = mvAllocMat(model->X->nrows, model->_A);
        mvMat *newP = mvAllocMat(model->X->ncolumns, model->_A);
        // T and P
        mvConcatColumns(newT, model->t, t);
        mvFreeMat(&model->t);
        model->t=newT;
        mvConcatColumns(newP, model->p, p);
        mvFreeMat(&model->p);
        model->p=newP;
        mvFreeMat(&t);
        mvFreeMat(&p);

        // R2X
        mvConcatRows(newR2, model->R2X, R2);
        mvFreeMat(&model->R2X);
        model->R2X=newR2;
        mvFreeMat(&R2);
    }
    else
    {
        model->t = t;
        model->p = p;
        model->R2X = R2;
    }

    // Store new sum of squares values.
    {
        int i=0;
        mvMat *SSX = NULL;
        mvMat *SSXV = NULL;

        //SSX
        SSX = mvAllocMat(model->SSX->nrows+1, 1);
        // Copy the data in so that we don't need another container.
        for (i=0; i<model->SSX->nrows; i++)
        {
            SSX->data[i][0] = model->SSX->data[i][0];
        }
        SSX->data[model->SSX->nrows][0] = SSE;
        mvFreeMat(&model->SSX);
        model->SSX = SSX;

        // SSXV
        SSXV = mvAllocMat(model->SSXV->nrows+1, 1);
        mvConcatRows(SSXV, model->SSXV, SSEV);
        mvFreeMat(&model->SSXV);
        mvFreeMat(&SSEV);
        model->SSXV = SSXV;
    }

    model->A=model->_A;
    return 0;
}

static int __mvAddPLSComponent(mvModel *model, int performCrossValidation)
{
    int iter = 0;
    int N = model->X->nrows;
    int K = model->X->ncolumns;
    int M = model->Y->ncolumns;
    mvMat *w = mvAllocMatVal(K, 1, 1.0);
    mvMat *wOld = mvAllocMat(K, 1);
    mvMat *wDiff = mvAllocMat(K, 1);
    mvMat *p = mvAllocMat(K, 1);
    mvMat *c = mvAllocMat(M, 1);
    mvMat *t = mvAllocMat(N, 1);
    mvMat *u = mvAllocMat(N, 1);
    mvMat *R2X = mvAllocMat(1, 1);
    mvMat *R2Y = mvAllocMat(1, 1);
    mvMat *X, *Y;   // reference -> no clean up req'd.
    mvMat *pT, *tpT, *cT, *tcT;
    double SSE, SSF;
    mvMat *SSEV = mvAllocMat(1, K);
    mvMat *SSFV = mvAllocMat(1, M);

    // Initialize w as a unit length vector.
    mvMatMultS(w, w, mvVectorNorm(w));

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
        iter++;

        __mvRegressCol(t, X, w);    // t = X w / w'w

        __mvRegressRow(c, Y, t);    // c = t'Y / t't

        __mvRegressCol(u, Y, c);    // u = Y c / c'c

        mvMatCopy(wOld, w);

        __mvRegressRow(w, X, u);    // w = u'X / u'u

        mvMatMultS(w, w, 1.0/mvVectorNorm(w));  // normalize w

        mvSubtractMat(wDiff, w, wOld);          // get the difference of w

    } while (mvVectorNorm(wDiff)>MV_DBL_SQRT_EPS && iter < MAX_NIPALS_ITER);

    // compute loading p after loop
    __mvRegressRow(p, X, t);

    // no longer need, wOld or wDiff
    mvFreeMat(&wOld);
    mvFreeMat(&wDiff);

    // Increment number of internal components
    model->_A++;

    // XXX: Cross validation must be performed before the new residual E is
    // computed and depends on _A == 1 for the first component (NOT ZERO!)
    if (performCrossValidation)
    {
        if (!model->cvd)
        {
            model->cvd = __allocCrossValData(model->numCrossValRounds);
        }
        __crossValidatePLS(model);
    }

    // compute residual E
    tpT = mvAllocMat(X->nrows, X->ncolumns);
    pT = mvAllocMat(1, X->ncolumns);
    mvTransposeMat(pT, p);
    mvMatMult(tpT, t, pT);
    mvSubtractMat(model->E, X, tpT);
    mvFreeMat(&pT);
    mvFreeMat(&tpT);

    // Compute residual F
    tcT = mvAllocMat (Y->nrows, Y->ncolumns);
    cT = mvAllocMat(1, Y->ncolumns);
    mvTransposeMat(cT, c);
    mvMatMult(tcT, t, cT);
    mvSubtractMat(model->F, Y, tcT);
    mvFreeMat(&cT);
    mvFreeMat(&tcT);

    // Compute R2;
    mvMatColumnSS(SSEV, model->E);
    mvMatColumnSS(SSFV, model->F);
    SSE = mvMatSum(SSEV);
    SSF = mvMatSum(SSFV);
    mvMatSetElem(R2X, 0, 0, 1.0 - SSE / model->SSX->data[0][0]);
    mvMatSetElem(R2Y, 0, 0, 1.0 - SSF / model->SSY->data[0][0]);

    // store new components
    if (model->_A > 1)
    {
        mvMat *newT, *newP, *newW, *newU, *newC, *wStar, *newWStar;
        mvMat *newR2X, *newR2Y;
        //t
        newT = mvAllocMat(model->X->nrows, model->A+1);
        mvConcatColumns(newT, model->t, t);
        mvFreeMat(&model->t);
        model->t=newT;
        mvFreeMat(&t);

        //p
        newP = mvAllocMat(model->X->ncolumns, model->A+1);
        mvConcatColumns(newP, model->p, p);
        mvFreeMat(&model->p);
        model->p=newP;
        mvFreeMat(&p);

        //w
        newW = mvAllocMat(model->X->ncolumns, model->A+1);
        mvConcatColumns(newW, model->w, w);
        mvFreeMat(&model->w);
        model->w=newW;
        mvFreeMat(&w);

        //u
        newU = mvAllocMat(model->Y->nrows, model->A+1);
        mvConcatColumns(newU, model->u, u);
        mvFreeMat(&model->u);
        model->u=newU;
        mvFreeMat(&u);

        //c
        newC = mvAllocMat(model->Y->ncolumns, model->A+1);
        mvConcatColumns(newC, model->c, c);
        mvFreeMat(&model->c);
        model->c=newC;
        mvFreeMat(&c);

        //W*
        wStar = mvAllocMat(K, 1);
        newWStar = mvAllocMat(model->X->ncolumns, model->A+1);
        __computeWStar(wStar, model->w, model->p, model->wStar);
        mvConcatColumns(newWStar, model->wStar, wStar);
        mvFreeMat(&model->wStar);
        model->wStar = newWStar;
        mvFreeMat(&wStar);

        //R2X
        newR2X = mvAllocMat(model->A+1, 1);
        mvConcatRows(newR2X, model->R2X, R2X);
        mvFreeMat(&model->R2X);
        model->R2X = newR2X;
        mvFreeMat(&R2X);

        //R2Y
        newR2Y = mvAllocMat(model->A+1, 1);
        mvConcatRows(newR2Y, model->R2Y, R2Y);
        mvFreeMat(&model->R2Y);
        model->R2Y = newR2Y;
        mvFreeMat(&R2Y);

    }
    else
    {
        model->t=t;
        model->p=p;
        model->c=c;
        model->u=u;
        model->w=w;
        model->wStar = mvAllocMatCopy(w); // W* = W for the first component
        model->R2X = R2X;
        model->R2Y = R2Y;
    }

    // Store new sum of squares values.
    {
        int i=0;
        mvMat *SSX = NULL;
        mvMat *SSXV = NULL;
        mvMat *SSY = NULL;
        mvMat *SSYV = NULL;

        //SSX
        SSX = mvAllocMat(model->SSX->nrows+1, 1);
        // Copy the data in so that we don't need another container.
        for (i=0; i<model->SSX->nrows; i++)
        {
            SSX->data[i][0] = model->SSX->data[i][0];
        }
        SSX->data[model->SSX->nrows][0] = SSE;
        mvFreeMat(&model->SSX);
        model->SSX = SSX;

        // SSXV
        SSXV = mvAllocMat(model->SSXV->nrows+1, 1);
        mvConcatRows(SSXV, model->SSXV, SSEV);
        mvFreeMat(&model->SSXV);
        mvFreeMat(&SSEV);
        model->SSXV = SSXV;

        //SSY
        SSY = mvAllocMat(model->SSY->nrows+1, 1);
        // Copy the data in so that we don't need another container.
        for (i=0; i<model->SSY->nrows; i++)
        {
            SSY->data[i][0] = model->SSY->data[i][0];
        }
        SSY->data[model->SSY->nrows][0] = SSF;
        mvFreeMat(&model->SSY);
        model->SSY = SSY;

        // SSYV
        SSYV = mvAllocMat(model->SSYV->nrows+1, 1);
        mvConcatRows(SSYV, model->SSYV, SSFV);
        mvFreeMat(&model->SSYV);
        mvFreeMat(&SSFV);
        model->SSYV = SSYV;
    }

    model->A=model->_A;
    return 0;
}

int mvModelAddComponent(mvModel *model)
{
    switch(model->modelType)
    {
    case PCA:
        return __mvAddPCAComponent(model, 1);
    case PLS:
        return __mvAddPLSComponent(model, 1);
    default:
        return UNKNOWN_MODEL_TYPE;
    }
    return UNKNOWN_MODEL_TYPE;
}


static int __mvNewObsPCA_T(mvMat *t, mvMat *E, const mvMat *newX, const mvModel *model,
                  int num_components, MVNewScoreCalcType method)
{
    int a, i, j;
    int freeE=0;
    double numerator, denominator;
    mvMat *p, *EHat;
    mvMat *_t, *_p, *_p_T, *_slice; // slices for calculation of E-hat
    if ( !(t->nrows == newX->nrows && newX->ncolumns == model->p->nrows &&
           model->p->ncolumns >= num_components && t->ncolumns == num_components))
    {
        return INCORRECT_DIMENSIONS;
    }

    if (E)
    {
        if (! (E->nrows == newX->nrows && E->ncolumns == newX->ncolumns))
        {
            return INCORRECT_DIMENSIONS;
        }
        mvMatCopy(E, newX);
    }
    else
    {
        E = mvAllocMatCopy(newX);
        freeE = 1;
    }

    // method is not yet used.
    (void)method;

    p = model->p;

    _t = mvAllocMat(t->nrows, 1);
    _p = mvAllocMat(p->nrows, 1);
    _p_T = mvAllocMat(1, p->nrows);
    _slice = mvAllocMat(1, 1);
    EHat = mvAllocMat(E->nrows, E->ncolumns);

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
        mvMatSliceColumns(_p, p, _slice);
        mvMatSliceColumns(_t, t, _slice);
        mvTransposeMat(_p_T, _p);
        mvMatMult(EHat, _t, _p_T);  // XHat = tpT
        mvSubtractMat(E, E, EHat);
    }

    mvFreeMat(&EHat);
    mvFreeMat(&_slice);
    mvFreeMat(&_p_T);
    mvFreeMat(&_p);
    mvFreeMat(&_t);
    if (freeE)
    {
        mvFreeMat(&E);
    }

    return SUCCESS;
}

static int __mvNewObsPLS_T(mvMat *t, mvMat *E, const mvMat *newX, const mvModel *model,
                  int num_components, MVNewScoreCalcType method)
{
    int a, i, j;
    int freeE = 0;
    double numerator, denominator;
    mvMat *p, *w, *EHat;
    mvMat *_t, *_p, *_p_T, *_slice; // slices for calculation of E-hat
    if ( !(t->nrows == newX->nrows && newX->ncolumns == model->p->nrows &&
           model->p->ncolumns >= num_components &&
           t->ncolumns == num_components && model->w->ncolumns >= num_components))
    {
        return INCORRECT_DIMENSIONS;
    }

    if (E)
    {
        if (! (E->nrows == newX->nrows && E->ncolumns == newX->ncolumns))
        {
            return INCORRECT_DIMENSIONS;
        }
        mvMatCopy(E, newX);
    }
    else
    {
        E = mvAllocMatCopy(newX);
        freeE = 1;
    }

    // method is not yet used
    (void) method;

    p = model->p;
    w = model->w;

    _t = mvAllocMat(t->nrows, 1);
    _p = mvAllocMat(p->nrows, 1);
    _p_T = mvAllocMat(1, p->nrows);
    _slice = mvAllocMat(1, 1);
    EHat = mvAllocMat(E->nrows, E->ncolumns);

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
        mvMatSliceColumns(_p, p, _slice);
        mvMatSliceColumns(_t, t, _slice);
        mvTransposeMat(_p_T, _p);
        mvMatMult(EHat, _t, _p_T);  // XHat = tpT
        mvSubtractMat(E, E, EHat);
    }

    mvFreeMat(&EHat);
    mvFreeMat(&_slice);
    mvFreeMat(&_p_T);
    mvFreeMat(&_p);
    mvFreeMat(&_t);
    if (freeE)
    {
        mvFreeMat(&E);
    }

    return SUCCESS;
}

int mvNewObsU(mvMat *u, mvMat *F, const mvMat *newY, const mvMat *newT,
                  const mvModel *model, int num_components,
                  MVNewScoreCalcType method)
{
    int a, i, j;
    int freeF = 0;
    double numerator, denominator;
    mvMat *c, *FHat;
    mvMat *_t, *_c, *_c_T, *_slice; // slices for calculation of E-hat
    if ( !(u->nrows == newY->nrows && newY->ncolumns == model->c->nrows &&
           model->c->ncolumns >= num_components &&
           u->ncolumns == num_components && u->nrows == newT->nrows &&
           u->ncolumns == newT->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }
    else if (model->modelType != PLS)
    {
        return WRONG_MODEL_TYPE;
    }
    if (F)
    {
        if (! (F->nrows == newY->nrows && F->ncolumns == newY->ncolumns))
        {
            return INCORRECT_DIMENSIONS;
        }
        mvMatCopy(F, newY);
    }
    else
    {
        F = mvAllocMatCopy(newY);
        freeF = 1;
    }

    // method is not yet used
    (void) method;
    c = model->c;

    _t = mvAllocMat(newT->nrows, 1);
    _c = mvAllocMat(c->nrows, 1);
    _c_T = mvAllocMat(1, c->nrows);
    _slice = mvAllocMat(1, 1);
    FHat = mvAllocMat(F->nrows, F->ncolumns);

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
        mvMatSliceColumns(_c, c, _slice);
        mvMatSliceColumns(_t, newT, _slice);
        mvTransposeMat(_c_T, _c);
        mvMatMult(FHat, _t, _c_T);  // FHat = tcT
        mvSubtractMat(F, F, FHat);
    }

    mvFreeMat(&FHat);
    mvFreeMat(&_slice);
    mvFreeMat(&_c_T);
    mvFreeMat(&_c);
    mvFreeMat(&_t);
    if (freeF)
    {
        mvFreeMat(&F);
    }

    return SUCCESS;
}


int mvNewObsT(mvMat *t, mvMat * E, const mvMat *newX, const mvModel *model,
              int num_components, MVNewScoreCalcType method)
{
    switch (model->modelType)
    {
    case PCA:
        return __mvNewObsPCA_T(t, E, newX, model, num_components, method);
    case PLS:
        return __mvNewObsPLS_T(t, E, newX, model, num_components, method);
    default:
        return UNKNOWN_MODEL_TYPE;
    }
    return UNKNOWN_MODEL_TYPE;
}
