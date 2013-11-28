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

#ifndef MVALGORITHMS_H
#define MVALGORITHMS_H


#include "mvmatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TODO Define error code enums */

/*! New score calculation methods:
    References: Arteaga, F and Ferrer, A. Journal of Chemometrics, 16, pp408-418, 2002.
*/
typedef enum MVNewScoreCalcType_enum {
    MV_NEW_SCORE_SCP//, /*! Single component projection */
    //PMP, /*! Projection to model plane */
    //TSR  /*! Trimmed score regression */
} MVNewScoreCalcType;

typedef enum MVModelType_enum {
    MV_MODEL_TYPE_PCA,
    MV_MODEL_TYPE_PLS
} MVModelType;

typedef enum MVCrossValType_enum {
    MV_CROSSVAL_TYPE_FAST,
    MV_CROSSVAL_TYPE_FULL /* Not yet implemented */
} MVCrossValType;

typedef struct MVCrossValData_s{
    MVMat * PRESS;      /*! _Ax1 */
    MVMat * PRESSV;     /*! _AxK for PCA or _AxM for PLS */
    void ** models;     /*! array of models of size NUM_ROUNDS */
    int numRounds;      /*! Number of rounds in this cross val data */
} MVCrossValData;

/*! Structure for both PCA and PLS models
  */
typedef struct MVModel_s{
    MVModelType modelType;
    MVMat * X;           /*! Reference to X-matrix of size(NxK)*/
    MVMat * Y;           /*! Reference to Y-matrix of size(NxM) - PLS only */
    MVMat * E;           /*! X-Residual matrix */
    MVMat * F;           /*! Y-Residual matrix */
    MVMat * t;           /*! T-scores matrix of size (NxA) */
    MVMat * t_stddev;    /*! Standard deviation of T-scores of size (1xA); */
    MVMat * p;           /*! P-loadings matrix of size (KxA) */
    MVMat * u;           /*! U-scores matrix of size (NxA) */
    MVMat * w;           /*! W-weightings matrix of size (KxA) */
    MVMat * wStar;       /*! W*-weightings matrix of size (KxA) */
    MVMat * c;           /*! C-loadings matrix of size (MxA) */
    MVMat * R2X;         /*! Column vector MVMat (_Ax1)*/
    MVMat * R2Y;         /*! Column vector MVMat (_Ax1)*/
    MVMat * Q2cum;       /*! Column vector MVMat (_Ax1)*/
    MVMat * Q2Vcum;      /*! Column vector MVMat (_Ax1)*/
    MVMat * Q2;          /*! Column vector MVMat (_Ax1)*/
    MVMat * Q2V;         /*! Q2 per Variable matrix of size KxA (PCA) or MxA (PLS) */
    MVMat * SPEX;        /*! Squared Prediction Error of X (size NxA) */
    MVMat * SPEY;        /*! Squared Prediction Error of Y - PLS only (size NxA) */
    MVCrossValData * cvd;  /*! Data required for crossValidation */
    MVCrossValType crossValType; /*! Type of cross validation */
    int numCrossValRounds; /*! Number of rounds to for cross validation (default = 7) */
    MVMat * SSX;         /*! Sum of squares of X-matrix for each component size (_A+1x1)*/
    MVMat * SSY;         /*! Sum of squares of Y-matrix for each component size (_A+1x1)*/
    MVMat * SSXV;        /*! Sum of squares of columns of X-matrix for each component size (_A+1xK)*/
    MVMat * SSYV;        /*! Sum of squares of columns of Y-matrix for each component size (_A+1xY)*/
    int A;      /*! Number of active components */
    int _A;     /*! Total number of computed components */
    MVMat * iter;        /*! Column vector MVMat (_Ax1) Number of iterations on NIPALS algorithm */

} MVModel;

/*! Initializes a PCA model

  This model needs to be freed with mvFreeModel.  This model does not take
  ownsership of X and it will not free X automatically.

  \arg X The X-Matrix of the PCA model
  \return mvModel of type PCA

  \sa mvFreeModel, MVModelType, MVCrossValType
  */
MVModel * mvmodel_alloc_init_pca(MVMat *X);

/*! Initializes a PLS model

  This model needs to be freed with mvFreeModel.  This model does not take
  ownsership of X or Y and it will not free X or Y automatically.

  \arg X The X-Matrix of the PLS model
  \arg Y The Y-Matrix of the PLS Model
  \return mvModel of type PLS
  \sa mvFreeModel, MVModelType
  */
MVModel * mvmodel_alloc_init_pls(MVMat *X, MVMat *Y);

/*! Frees the mvModel
  \arg model pointer to the mvModel pointer.
  \return 0 on success or MVModelErrorCode
  */
int mvmodel_free(MVModel **model);

/*! Add a component to the model

  Components are added via NIPALs algorithm.
  All calculations for relevant parameters, including cross-validation, are
  performed.

  References:
   - Summary of PCA and PLS algorithms: Westerhuis, J.A., Kourti, T., &
     MacGregor, J.F., ANALYSIS OF MULTIBLOCK AND HIERARCHICAL PCA AND PLS MODELS,
     J Chemometrics 12, 301–321 (1998)
   - Calculation of W*: Dayal, B.S. & MacGregor, J.F. (1997) Improved PLS
     algorithms. J. Chemometrics 11:73-85
   - "Fast" cross-validation is computed as in: Eriksson, L., Johansson, E.,
     Kettaneh-Wold, N., Trigg, J., C., W., & Wold.S. (2006). Multi- and
     Megavariate Data Analysis: Part 1 - Basic Principles and Applications. Umea,
     Sweden: Umetrics AB.

  \arg model The mvModel
  \return 0 on success or MVModelErrorCode;
  */
int mvmodel_add_component(MVModel *model);

/*! Performs auto-fit on PCA or PLS model

  Cross validation rules:
   - RULE 1: New components are added if Q2 > threshold
   - RULE 2: New components are added if any variable in Q2V > threshold (PLS)
     or at least sqrt(K) variables have Q2V > threshold.
   - RULE 3: New components stop if A > min(X.rows, X.cols)
   - RULE 4: Cross validation will stop of the NIPALS iteration threshold has
     been reached implying convergence cannot be reached which means the model
     may have no more data to extract.  The max number of NIPALS iterations is
     500.

   Rules 3 and 4 are fail rules.  If Rules 3 and 4 do not fail, Rule 1 or 2 must
   pass.

   References:
    - Cross validation rules as in: Eriksson, L., Johansson, E.,
     Kettaneh-Wold, N., Trigg, J., C., W., & Wold.S. (2006). Multi- and
     Megavariate Data Analysis: Part 1 - Basic Principles and Applications. Umea,
     Sweden: Umetrics AB.

   \arg model mvModel to auto-fit
   \return MVModelReturnCode - The rule that failed.
  */
int mvmodel_autofit(MVModel *model);

/*! Computes T-scores for new observations of X of a mvModel

  \arg t output preallocated vector of scores.  Must be size N' x A' where
       A' is num_components
  \arg E output preallocated vector of residuals. Must be size N'xK
  \arg newX matrix of new observations of size N'xK.
  \arg pca_model mvModel object containing model weights.
  \arg num_components Number of components of the new observation to compute.
  \arg method Method with which to compute new scores.

  \sa MVNewScoreCalcType
  */
int mvmodel_new_obs_scores_t(MVMat *t, MVMat * E, const MVMat *newX, const MVModel *pca_model,
                  int num_components, MVNewScoreCalcType method);

/*! Computes U-scores for a new observation of Y of a PLS Model

  \arg u output preallocated vector of scores.  Must be size N' x A' where
       A' is num_components
  \arg F output preallocated vector of residuals. Must be size N'xM
  \arg newY matrix of new observations of size N'xM.
  \arg pls_model mvModel object containing model weights.
  \arg num_components Number of components of the new observation to compute.
  \arg method Method with which to compute new scores.

  \sa MVNewScoreCalcType
  */
int mvmodel_new_obs_scores_u(MVMat *u, MVMat *F, const MVMat *newY, const MVMat *newT,
                  const MVModel *pls_model, int num_components,
                  MVNewScoreCalcType method);

/*! Compute the predicted X values based on scores and loadings provided.

    Predictions for X: Xhat = TP'

    \arg pred preallocated MVMat, must be size NxK
    \arg model mvModel.
    \arg t scores with which to compute Xhat.
    \arg num_components number of components to compute Xhat
*/
int mvmodel_compute_xpred(MVMat *Xhat, const MVModel *model, const MVMat *t, int num_components);

/*! Compute the predicted X values based on scores and loadings provided.

    Predictions for Y: Yhat = TC'

    \arg pred preallocated MVMat, must be size NxK
    \arg model mvModel must be PLS model.
    \arg t scores with which to compute Yhat
    \arg num_components number of components to compute Xhat
*/
int mvmodel_compute_ypred(MVMat *Yhat, const MVModel *model, const MVMat *t, int num_components);



#ifdef __cplusplus
}
#endif

#endif // MVALGORITHMS_H
