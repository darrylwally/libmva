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

double mvstats_t_ppf(double alpha, double df, double loc, double scale);

double mvstats_F_ppf(double alpha, double N1, double N2, double loc, double scale);

double mvstats_gamma_ppf(double alpha, double a, double loc, double scale);

double mvstats_chi2_ppf(double alpha, double a, double loc, double scale);

int mvstats_ht2(MVMat *output, const MVMat *scores, const MVMat *t_stddev,
          int first_component, int last_component);

double mvstats_ht2_limit(double alpha, int A, int N);

int mvstats_spe(MVMat *output, const MVMat *residuals);
int mvstats_spex_from_obs(MVMat *output, const MVModel *model, const MVMat *Xobs, const MVMat *tobs, int num_components);
int mvstats_spey_from_obs(MVMat *output, const MVModel *model, const MVMat *Yobs, const MVMat *tobs, int num_components);

double mvstats_spe_limit(double alpha, const MVMat *modelSPE_values, int component);


#endif // MVSTATS_H
