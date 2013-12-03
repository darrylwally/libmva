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


MVPreprocessContext * mvpreprocess_alloc()
{
    MVPreprocessContext *out = malloc(sizeof(MVPreprocessContext));
    out->matrix = NULL;
    out->preprocess_info = NULL;
    return out;
}

MVPreprocessContext * mvpreprocess_alloc_mat(MVMat *matrix)
{
    MVPreprocessContext *out = mvpreprocess_alloc();
    if (out)
    {
        out->matrix = matrix;
    }
}
