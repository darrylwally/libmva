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

/*!
    Arteaga and Ferrer
*/
typedef enum MVNewScoreCalcType_enum {
    SCP//, /*! Single component projection */
    //PMP, /*! Projection to model plane */
    //TSR  /*! Trimmed score regression */
} MVNewScoreCalcType;

typedef enum MVModelType_enum {
    PCA,
    PLS
} MVModelType;

typedef enum MVCrossValType_enum {
    FAST,
    FULL /* Not yet implemented */
} MVCrossValType;

typedef struct crossValData_s{
    mvMat * PRESS;      /*! _Ax1 */
    mvMat * PRESSV;     /*! _AxK for PCA or _AxM for PLS */
    void ** models;     /*! array of models of size NUM_ROUNDS */
    int numRounds;      /*! Number of rounds in this cross val data */
} crossValData;

/*! Structure for both PCA and PLS models
  */
typedef struct mvModel_s{
    MVModelType modelType;
    mvMat * X;           /*! Reference to X-matrix of size(NxK)*/
    mvMat * Y;           /*! Reference to Y-matrix of size(NxM) - PLS only */
    mvMat * E;           /*! X-Residual matrix */
    mvMat * F;           /*! Y-Residual matrix */
    mvMat * t;           /*! T-scores matrix of size (NxA) */
    mvMat * p;           /*! P-loadings matrix of size (KxA) */
    mvMat * u;           /*! U-scores matrix of size (NxA) */
    mvMat * w;           /*! W-weightings matrix of size (KxA) */
    mvMat * wStar;       /*! W*-weightings matrix of size (KxA) */
    mvMat * c;           /*! C-loadings matrix of size (MxA) */
    mvMat * R2X;         /*! Column vector mvMat (_Ax1)*/
    mvMat * R2Y;         /*! Column vector mvMat (_Ax1)*/
    mvMat * Q2cum;       /*! Column vector mvMat (_Ax1)*/
    mvMat * Q2Vcum;      /*! Column vector mvMat (_Ax1)*/
    mvMat * Q2;          /*! Column vector mvMat (_Ax1)*/
    mvMat * Q2V;         /*! Q2 per Variable matrix of size KxA (PCA) or MxA (PLS) */
    crossValData * cvd;  /*! Data required for crossValidation */
    MVCrossValType crossValType; /*! Type of cross validation */
    int numCrossValRounds; /*! Number of rounds to for cross validation (default = 7) */
    mvMat * SSX;         /*! Sum of squares of X-matrix for each component size (_A+1x1)*/
    mvMat * SSY;         /*! Sum of squares of Y-matrix for each component size (_A+1x1)*/
    mvMat * SSXV;        /*! Sum of squares of columns of X-matrix for each component size (_A+1xK)*/
    mvMat * SSYV;        /*! Sum of squares of columns of Y-matrix for each component size (_A+1xY)*/
    int A;      /*! Number of active components */
    int _A;     /*! Total number of computed components */

} mvModel;

/*! Initializes a PCA model

  This model needs to be freed with mvFreeModel.  This model does not take
  ownsership of X and it will not free X automatically.

  \arg X The X-Matrix of the PCA model
  \return mvModel of type PCA

  \sa mvFreeModel, MVModelType, MVCrossValType
  */
mvModel * mvInitPCAModel(mvMat *X);

/*! Initializes a PLS model

  This model needs to be freed with mvFreeModel.  This model does not take
  ownsership of X or Y and it will not free X or Y automatically.

  \arg X The X-Matrix of the PLS model
  \arg Y The Y-Matrix of the PLS Model
  \return mvModel of type PLS
  \sa mvFreeModel, MVModelType
  */
mvModel * mvInitPLSModel(mvMat *X, mvMat *Y);

/*! Frees the mvModel
  \arg model pointer to the mvModel pointer.
  \return 0 on success or MVModelErrorCode
  */
int mvFreeModel(mvModel **model);

/*! Add a component to the model

  Components will not be added if A >min(X->nrows, X->ncolumns)
  \arg model The mvModel
  \return 0 on success or MVModelErrorCode;
  */
int mvModelAddComponent(mvModel *model);

/*! Computes T-scores for new observations of X of a mvModel

  */
int mvNewObsT(mvMat *t, mvMat * E, const mvMat *newX, const mvModel *pca_model,
                  int num_components, MVNewScoreCalcType method);

/*! Computes U-scores for a new observation of Y of a PLS Model

  */
int mvNewObsU(mvMat *u, mvMat *F, const mvMat *newY, const mvMat *newT,
                  const mvModel *pls_model, int num_components,
                  MVNewScoreCalcType method);


#ifdef __cplusplus
}
#endif

#endif // MVALGORITHMS_H
