/*! \file mvfile.c
  \brief Contains implementation of storing and retrieving mvmodels from files.
  \author Darryl Wallace - wallacdj@gmail.com, darryl@wallynet.ca
  \date 2012

  THIS CODE IS COPYRIGHT DARRYL WALLACE (c) 2012.
  UNAUTHORIZED DISTRIBUTION, USE, OR VIEWING THE CONTENTS OF THIS FILE
  IS STRICTLY PROHIBITED.
  IF YOU HAVE RECEIVED THIS CODE WITHOUT AUTHORIZATION IT MUST BE DESTROYED
  IMMEDIATELY
  */

#include "mvfile.h"
#include <stdint.h>
#include <stdio.h>

static const uint64_t VERSION = 0x000000001;


/* File header consists of the following information
  uint64_t - version
  uint64_t - number_of_models
  uint64_t * number of models - size in bytes of each model.
  MVModel * number of models - the mv model data is serialized.

  Models are written sequentially based on their size as tightly packed as
  possible.  Start position of each model could be inferred by taking the size
  of each model past the header.
  */
int  mvfile_write_header(FILE *file, MVModel *models, int num_models)
{
    int i;
    size_t total = fwrite(&VERSION, sizeof(uint64_t), 1, file);
    uint64_t num = (uint64_t) num_models;
    total += fwrite(&num, sizeof(uint64_t), 1, file);

    MVModel * model = models;
    for (i = 0; i < num_models; i++)
    {
        num = __mvfile_compute_model_size(model);
        total += fwrite(&num, sizeof(uint64_t, 1, file));
        model++;
    }
}
