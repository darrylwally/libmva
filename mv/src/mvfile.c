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

  Models are written sequentially based on their size as tightly packed as
  possible.  Start position of each model could be inferred by taking the size
  of each model past the header.
  */
int  mvfile_write_header(FILE *file, const MVModel *models, int num_models)
{
    int i;
    size_t total = fwrite(&VERSION, sizeof(uint64_t), 1, file);
    uint64_t num = (uint64_t) num_models;
    total += fwrite(&num, sizeof(uint64_t), 1, file);

    const MVModel * model = models;
    for (i = 0; i < num_models; i++)
    {
        num = __mvfile_compute_model_size(model);
        total += fwrite(&num, sizeof(uint64_t), 1, file);
        model++;
    }
    return 0;
}

static int __mvfile_serialize_mat(FILE *file, const MVMat *mat)
{
    fwrite(&mat->nrows, sizeof(mat->nrows), 1, file);
    fwrite(&mat->ncolumns, sizeof(mat->ncolumns), 1, file);
    fwrite(mat->data[0], sizeof(double), sizeof(double)*mat->nrows*mat->ncolumns, file);
    fwrite(mat->mask[0], sizeof(size_t), sizeof(size_t)*mat->nrows*mat->ncolumns, file);
    return 0;
}

static int __mvfile_deserialize_mat(MVMat **mat, FILE *file)
{
    MVMat *m = NULL;
    int nrows, ncolumns;
    fread(&nrows,sizeof(int),1, file);
    fread(&ncolumns,sizeof(int),1, file);

    if (nrows > 0 && ncolumns > 0)
    {
        m = mvmat_alloc(nrows, ncolumns);
        fread(m->data[0], sizeof(double), sizeof(double)*nrows*ncolumns, file);
        fread(m->mask[0], sizeof(size_t), sizeof(size_t)*nrows*ncolumns, file);
    }
    *mat = m;
    return 0;
}
