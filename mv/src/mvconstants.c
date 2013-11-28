/*! \file mvconstants.c
  \brief Contains various constants used in mvlib.
  \author Darryl Wallace - wallacdj@gmail.com, darryl@wallynet.ca
  \date 2012

  THIS CODE IS COPYRIGHT DARRYL WALLACE (c) 2012.
  UNAUTHORIZED DISTRIBUTION, USE, OR VIEWING THE CONTENTS OF THIS FILE
  IS STRICTLY PROHIBITED.
  IF YOU HAVE RECEIVED THIS CODE WITHOUT AUTHORIZATION IT MUST BE DESTROYED
  IMMEDIATELY
  */

#include "mvconstants.h"

const double MV_EPS = 2.2204460492503131e-016;
const double MV_SQRT_EPS = 1.4901161193847656e-08;
const float MV_EPS_F = 1.192092896e-07;
const float MV_SQRT_EPS_F = 3.452669831e-04;

static const long long __mv_nan = 0x7ff8000000000000;
static const int __mv_nanf = 0x7fc00000;
static const unsigned long long __mv_inf = 0x7ff0000000000000;
static const unsigned long long __mv_neg_inf = 0xfff0000000000000;
static const unsigned int __mv_inff = 0x7f800000;
static const unsigned int __mv_neg_inff = 0xff800000;

double mv_NaN()
{
    return *(double *) &__mv_nan;
}

float mv_NaNf()
{
    return *(float *)&__mv_nanf;
}

double mv_inf()
{

    return *(double *)&__mv_inf;
}

double mv_neg_inf()
{

    return *(double *)&__mv_neg_inf;
}

float mv_inff()
{
    return *(float *)&__mv_inff;
}

float mv_neg_inff()
{
    return *(float *)&__mv_neg_inff;
}
