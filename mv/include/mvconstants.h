/*
* Filename: mvconstants.h
* Description: Contains various constants used in mvlib.
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

#ifndef MVCONSTANTS_H
#define MVCONSTANTS_H

#ifdef __cplusplus
extern "C"{
#endif

#include <math.h>

#define DATA_MISSING 0
#define DATA_PRESENT 1

typedef enum MVReturnCode_enum {
    ERROR = -1,
    INCORRECT_DIMENSIONS = -2,
    INDEX_OUT_OF_BOUNDS = -3,
    ATTEMPT_TO_FREE_NULL_MATRIX = -4,
    REFERENCE_MAT_REQUIRED = -5,
    WRONG_MODEL_TYPE = -6,
    UNKNOWN_MODEL_TYPE = -7,
    UNKOWN_CROSSVALIDATION_TYPE = -8,
    SUCCESS = 0,
    CROSSVAL_RULE1 = 1,
    CROSSVAL_RULE2 = 2,
    CROSSVAL_RULE3 = 3,
    CROSSVAL_RULE4 = 4
} MVReturnCode;

#ifdef WIN32
#include <float.h>
#define MVISNAN_FUNC _isnan
#else
#define MVISNAN_FUNC isnan
#endif

// EPS
extern const double MV_EPS;
extern const double MV_SQRT_EPS;
extern const float MV_EPS_F;
extern const float MV_SQRT_EPS_F;


// NANS
double mv_NaN();
float mv_NaNf();

// Infinite
double mv_inf();
double mv_neg_inf();
float mv_inff();
float mv_neg_inff();

#ifdef __cplusplus
}
#endif

#endif // MVCONSTANTS_H
