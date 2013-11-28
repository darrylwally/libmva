/*! \file mvstats.h
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

#ifndef MVSTATS_H
#define MVSTATS_H

#include "mvmatrix.h"
#include "mvmodel.h"

double t_ppf(double alpha, double df, double loc, double scale);

double F_ppf(double alpha, double N1, double N2, double loc, double scale);

double gamma_ppf(double alpha, double a, double loc, double scale);

double chi2_ppf(double alpha, double a, double loc, double scale);

int mvHT2(MVMat *output, const MVMat *scores, const MVMat *t_stddev,
          int first_component, int last_component);

double mvHT2Limit(double alpha, int A, int N);

int mvSPE(MVMat *output, const MVMat *residuals);
int mvSPEXFromObs(MVMat *output, const MVModel *model, const MVMat *Xobs, const MVMat *tobs, int num_components);
int mvSPEYFromObs(MVMat *output, const MVModel *model, const MVMat *Yobs, const MVMat *tobs, int num_components);

double mvSPELimit(double alpha, const MVMat *modelSPE_values, int component);


#endif // MVSTATS_H
