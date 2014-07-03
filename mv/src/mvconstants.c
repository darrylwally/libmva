/*
* Filename: mvconstants.c
* Description: Implementation of various constants used in mvlib.
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

#include "mvconstants.h"

const double MV_EPS = 2.2204460492503131e-016;
const double MV_SQRT_EPS = 1.4901161193847656e-08;
const float MV_EPS_F = 1.192092896e-07f;
const float MV_SQRT_EPS_F = 3.452669831e-04f;

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
