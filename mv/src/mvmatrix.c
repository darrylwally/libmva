/*! \file mvmatrix.c
  \brief Contains structures and functions for initializing matrix
  structures and performing standard matrix algebra funcitons.
  \author Darryl Wallace - wallacdj@gmail.com, darryl@wallynet.ca
  \date 2012

  THIS CODE IS COPYRIGHT DARRYL WALLACE (c) 2012.
  UNAUTHORIZED DISTRIBUTION, USE, OR VIEWING THE CONTENTS OF THIS FILE
  IS STRICTLY PROHIBITED.
  IF YOU HAVE RECEIVED THIS CODE WITHOUT AUTHORIZATION IT MUST BE DESTROYED
  IMMEDIATELY
  */

#include "mvmatrix.h"
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#ifdef WIN32
#define MVISNAN_FUNC _isnan
#else
#define MVISNAN_FUNC isnan
#endif

mvMat * mvAllocMat(int rows, int columns)
{
    int i;
    double *data;
    size_t *mask, *mask_end;
    mvMat *output;
    data = NULL;
    mask = NULL;
    output = (mvMat *) malloc(sizeof(mvMat));

    if (!output)
    {
        return NULL;
    }
    data = (double *) malloc(rows * columns * sizeof(double));
    if (!data)
    {
        free(output);
        return NULL;
    }
    mask = (size_t *) malloc(rows * columns * sizeof(size_t));
    if (!mask)
    {
        free(output);
        free(data);
        return NULL;
    }

    output->data = (double **) malloc(rows * sizeof(double *));
    if (!output->data)
    {
        free(output);
        free(data);
        free(mask);
        return NULL;
    }
    output->mask = (size_t **) malloc(rows * sizeof(size_t *));
    if (!output->mask)
    {
        free(output);
        free(data);
        free(mask);
        free(output->data);
        return NULL;
    }
    output->nrows = rows;
    output->ncolumns = columns;
    output->isReference = 0;

    // ensure all data is allocated contiguously.
    for(i=0; i<rows; i++)
    {
        output->data[i] = data + i*columns;
        output->mask[i] = mask + i*columns;
    }

    // need to loop through each mask value to set the mask to DATA_PRESENT
    mask_end = &mask[rows*columns - 1];
    *mask = DATA_PRESENT;
    while (mask != mask_end)
    {
        mask++;
        *mask=DATA_PRESENT;
    }
    return output;
}

mvMat * mvAllocMatZ(int rows, int columns)
{
    mvMat *output = mvAllocMat(rows, columns);
    memset(output->data[0], 0, rows * columns * sizeof(double));
    return output;
}

mvMat *mvAllocMatVal(int rows, int columns, double val)
{
    mvMat *output = mvAllocMat(rows, columns);
    mvMatSet(output, val);
    return output;
}

mvMat *mvAllocMatCopy(const mvMat *other)
{
    mvMat *output = mvAllocMat(other->nrows, other->ncolumns);
    mvMatCopy(output, other);
    return output;
}

mvMat * mvAllocIdentity(int dim)
{
    int i;
    mvMat *output;
    if (dim == 0)
    {
        return NULL;
    }
    output = mvAllocMatZ(dim, dim);
    if (!output)
        return NULL;

    for (i=0; i<dim; i++)
    {
        output->data[i][i]=1.0;
    }
    return output;
}

mvMat * mvRange(int start, int stop, int step)
{
    int range;
    int num_elements;
    int i;
    mvMat *output;
    double *data;
    if ( (start < stop && step<0) || (start > stop && step > 0) ||
         step==0.0 || fabs(step) <= MV_DBL_SQRT_EPS)
    {
        return NULL;
    }

    // Figure out size of range
    range = stop - start;
    num_elements = (int)ceil(fabs(range)/fabs(step));

    output = mvAllocMat(num_elements, 1);
    data = output->data[0];


    for (i=0; i<num_elements; i++)
    {
        data[i] = start + (step * i);
    }
    return output;
}

mvMat * mvLinspace(double start, double stop, int num_values, int end_point)
{
    double range, step, *data;
    mvMat *output;
    int i;
    if (num_values<=0)
    {
        return NULL;
    }

    range = stop - start;
    if (!end_point)
    {
        step = range/(double)num_values;
    }
    else
    {
        step = range/((double)(num_values - 1.0));
    }

    output = mvAllocMat(num_values, 1);
    data = output->data[0];

    for (i=0; i<num_values; i++)
    {
        data[i]=start + step * i;
    }

    return output;

}

mvMat * mvAllocMatRef(int rows)
{
    mvMat * output = (mvMat *)malloc(sizeof(mvMat));
    if (!output)
    {
        return NULL;
    }
    output->isReference = 1;
    output->data = (double **) malloc(rows * sizeof(double));
    output->mask = (size_t **) malloc(rows * sizeof(size_t));
    output->ncolumns = 0;
    output->nrows = rows;

    return output;
}

int mvFreeMat(mvMat **matrix)
{
    mvMat *mat = *matrix;
    if (!mat)
        return ATTEMPT_TO_FREE_NULL_MATRIX;
    if (!mat->isReference)
    {
        free(mat->data[0]);
        free(mat->mask[0]);
    }
    free(mat->data);
    free(mat->mask);
    free(mat);
    *matrix = NULL;
    return SUCCESS;
}

int mvMatCopy(mvMat *output, const mvMat *other)
{
    if (!(other->nrows==output->nrows && other->ncolumns==other->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }
    if (output->isReference || other->isReference)
    {
        int i,j;
        for (i=0; i<output->nrows; i++)
        {
            for (j=0; j<output->ncolumns; j++)
            {
                output->data[i][j] = other->data[i][j];
                output->mask[i][j] = other->mask[i][j];
            }
        }
    }
    else
    {
        int num_elements = other->nrows * other->ncolumns;
        int i;
        // This works because the contents are malloc'd in a contiguous array.
        double *output_data = output->data[0];
        size_t *output_mask = output->mask[0];
        double *other_data = other->data[0];
        size_t *other_mask = other->mask[0];
        for (i=0; i<num_elements; i++)
        {
            output_data[i]=other_data[i];
            output_mask[i]=other_mask[i];
        }
    }
    return SUCCESS;
}

int mvConcatColumns(mvMat *output, const mvMat *A, const mvMat *B)
{
    int i,j;
    if (!(output->ncolumns==(A->ncolumns+B->ncolumns) &&
          output->nrows==A->nrows && A->nrows==B->nrows))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<A->ncolumns; j++)
        {
            output->data[i][j] = A->data[i][j];
            output->mask[i][j] = A->mask[i][j];
        }
        for (j=A->ncolumns; j<output->ncolumns; j++)
        {
            output->data[i][j] = B->data[i][j-A->ncolumns];
            output->mask[i][j] = B->mask[i][j-A->ncolumns];
        }
    }
    return 0;
}

int mvConcatRows(mvMat *output, const mvMat *A, const mvMat *B)
{
    int i,j;
    if (!(output->nrows==(A->nrows+B->nrows) &&
          output->ncolumns==A->ncolumns && A->ncolumns==B->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<A->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            output->data[i][j] = A->data[i][j];
            output->mask[i][j] = A->mask[i][j];
        }
    }
    for (i=A->nrows; i<output->nrows; i++)
    {
        int bRow = i-A->nrows;
        for (j=0; j<output->ncolumns; j++)
        {
            output->data[i][j] = B->data[bRow][j];
            output->mask[i][j] = B->mask[bRow][j];
        }
    }
    return SUCCESS;
}

int mvMatSliceRows(mvMat *output, const mvMat *A, const mvMat *rows)
{
    int i,j;
    if ( !(output->nrows == rows->nrows && output->ncolumns == A->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<rows->nrows; i++)
    {
        int row = (int)rows->data[i][0];
        if (row >= A->nrows || row < 0)
        {
            return INDEX_OUT_OF_BOUNDS;
        }

        for (j=0; j<output->ncolumns; j++)
        {
            output->data[i][j] = A->data[row][j];
            output->mask[i][j] = A->mask[row][j];
        }
    }
    return SUCCESS;
}

int mvMatDeleteRows(mvMat *output, const mvMat *A, const mvMat *rows)
{
    int i;
    int num_rows = A->nrows - rows->nrows;
    int slice_index = 0;
    int out_index = 0;
    if ( !(output->nrows == num_rows && output->ncolumns == A->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }
    for (i=0; i<A->nrows; i++)
    {
        int row = -1;
        if (slice_index < (int)rows->nrows)
        {
            row = (int) rows->data[slice_index][0];
            if (row >= (int)A->nrows || row < 0)
            {
                return INDEX_OUT_OF_BOUNDS;
            }
        }
        if (row == i)
        {
            slice_index++;
        }
        else
        {
            output->data[out_index] = A->data[i];
            output->mask[out_index] = A->mask[i];
            out_index++;
        }
    }
    return SUCCESS;
}

int mvMatSliceRowsRef(mvMat **output, const mvMat *A, const mvMat *rows)
{
    int i;
    mvMat *o = *output;
    if (o && !o->isReference)
    {
        return REFERENCE_MAT_REQUIRED;
    }
    if (o && o->nrows != rows->nrows)
    {
        return INCORRECT_DIMENSIONS;
    }
    if (!o)
    {
        o = mvAllocMatRef(rows->nrows);
    }
    for (i=0; i<rows->nrows; i++)
    {
        int row = (int) rows->data[i][0];
        if (row >= A->nrows || row < 0)
        {
            mvFreeMat (&o);
            return INDEX_OUT_OF_BOUNDS;
        }
        o->data[i] = A->data[row];
        o->mask[i] = A->mask[row];
    }
    o->ncolumns = A->ncolumns;
    *output = o;
    return SUCCESS;
}

int mvMatDeleteRowsRef(mvMat **output, const mvMat *A, const mvMat *rows)
{
    int i;
    mvMat *o = *output;
    int num_rows = A->nrows - rows->nrows;
    int slice_index = 0;
    int out_index = 0;
    if (o && !o->isReference)
    {
        return REFERENCE_MAT_REQUIRED;
    }
    if (o && o->nrows != num_rows)
    {
        return INCORRECT_DIMENSIONS;
    }
    if (!o)
    {
        o = mvAllocMatRef(num_rows);
    }
    for (i=0; i<A->nrows; i++)
    {
        int row = -1;
        if (slice_index < (int)rows->nrows)
        {
            row = (int) rows->data[slice_index][0];
            if (row >= (int)A->nrows || row < 0)
            {
                mvFreeMat (&o);
                return INDEX_OUT_OF_BOUNDS;
            }
        }

        if (row == i)
        {
            slice_index++;
        }
        else
        {
            o->data[out_index] = A->data[i];
            o->mask[out_index] = A->mask[i];
            out_index++;
        }
    }
    o->ncolumns = A->ncolumns;
    *output = o;
    return SUCCESS;
}


int mvMatSliceColumns(mvMat *output, const mvMat *A, const mvMat *columns)
{
    int i,j;
    if ( !(output->ncolumns == columns->nrows && output->nrows == A->nrows))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (j=0; j<columns->nrows; j++)
    {
        int column = (int)columns->data[j][0];
        if (column >= A->ncolumns || column < 0)
        {
            return INDEX_OUT_OF_BOUNDS;
        }

        for (i=0; i<output->nrows; i++)
        {
            output->data[i][j] = A->data[i][column];
            output->mask[i][j] = A->mask[i][column];
        }
    }
    return SUCCESS;
}

int mvTransposeMat(mvMat *output, const mvMat *mat)
{
    int i,j;
    if (! (mat->ncolumns==output->nrows &&
           mat->nrows==output->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<mat->nrows; i++)
    {
        for (j=0; j<mat->ncolumns; j++)
        {
            output->data[j][i]=mat->data[i][j];
            output->mask[j][i]=mat->mask[i][j];
        }
    }
    return SUCCESS;
}

int mvMatSet(mvMat *mat, double val)
{
    if(!mat->isReference)
    {
        double *pointer = mat->data[0];
        double *end = &mat->data[mat->nrows-1][mat->ncolumns-1];
        /* This method works because we always allocate the data
           in a contiguous block */
        *pointer=val;
        while (pointer != end)
        {
            pointer++;
            *pointer=val;
        }
        if (MVISNAN_FUNC(val))
        {
            memset(mat->mask[0], DATA_MISSING,
                   mat->nrows*mat->ncolumns*sizeof(size_t));
        }
    }
    else
    {
        int i,j;
        size_t data_missing_val = MVISNAN_FUNC(val) ? DATA_MISSING : DATA_PRESENT;
        for (i=0; i<mat->nrows; i++)
        {
            for (j=0; j<mat->ncolumns; j++)
            {
                mat->data[i][j] = val;
                mat->mask[i][j] = data_missing_val;
            }
        }
    }
    return SUCCESS;
}

int mvMatSetElem(mvMat *mat, int row, int column, double value)
{
    if ( row> (mat->nrows-1)  || column>(mat->ncolumns-1) )
    {
        return INDEX_OUT_OF_BOUNDS;
    }

    mat->data[row][column] = value;
    if (MVISNAN_FUNC(value))
    {
        mat->mask[row][column] = DATA_MISSING;
    }
    else
    {
        mat->mask[row][column] = DATA_PRESENT;
    }
    return SUCCESS;
}

int mvMatGetElem(const mvMat *mat, double *value, int row, int column)
{
    if ( row> (mat->nrows-1)  || column>(mat->ncolumns-1) )
    {
        return INDEX_OUT_OF_BOUNDS;
    }
    if (mat->mask[row][column] == DATA_PRESENT)
    {
        *value = mat->data[row][column];
    }
    else
    {
        *value = mvNaN();
    }
    return SUCCESS;
}

int mvAddMat(mvMat *output, const mvMat *A, const mvMat *B)
{
    int i,j;
    if ( !(A->nrows == B->nrows && A->nrows == output->nrows &&
          A->ncolumns == B->ncolumns && A->ncolumns == output->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            double result = A->data[i][j] + B->data[i][j];
            output->data[i][j] = result;
            output->mask[i][j] = (MVISNAN_FUNC(result) ? DATA_MISSING : DATA_PRESENT);
        }
    }
    return SUCCESS;
}

int mvAddMatS(mvMat *output, const mvMat *A, double value)
{
    int i,j;
    if ( !(A->nrows == output->nrows && A->ncolumns == output->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            double result = A->data[i][j] + value;
            output->data[i][j] = result;
            output->mask[i][j] = (MVISNAN_FUNC(result) ? DATA_MISSING : DATA_PRESENT);
        }
    }
    return SUCCESS;
}

int mvSubtractMat(mvMat *output, const mvMat *A, const mvMat *B)
{
    int i,j;
    if ( !(A->nrows == B->nrows && A->nrows == output->nrows &&
          A->ncolumns == B->ncolumns && A->ncolumns == output->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            if (A->mask[i][j] == DATA_MISSING || B->mask[i][j] == DATA_MISSING)
            {
                output->data[i][j] = mvNaN();
                output->mask[i][j] = DATA_MISSING;
            }
            else
            {
                output->data[i][j] = A->data[i][j] - B->data[i][j];
                output->mask[i][j] = DATA_PRESENT;
            }
        }
    }
    return SUCCESS;
}


int mvMatMult(mvMat *output, const mvMat *A, const mvMat *B)
{
    int i,j,k;
    if ( !(A->ncolumns == B->nrows && A->nrows == output->nrows &&
          B->ncolumns == output->ncolumns ))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            double sum = 0.0;
            for (k=0; k<A->ncolumns; k++)
            {
                sum += A->data[i][k]*B->data[k][j];
            }
            output->data[i][j] = sum;
        }
    }
    return SUCCESS;
}

int mvMatMultS(mvMat *output, const mvMat *A, double scalar)
{
    int i,j;
    if (!(A->nrows == output->nrows && A->ncolumns == output->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            output->data[i][j] = A->data[i][j] * scalar;
            output->mask[i][j] = MVISNAN_FUNC(output->data[i][j]) ? DATA_MISSING : DATA_PRESENT;
        }
    }
    return SUCCESS;
}

int mvMatElemMult(mvMat *output, const mvMat *A, const mvMat *B)
{
    int i,j;
    if ( !(A->nrows == B->nrows && A->nrows == output->nrows &&
          A->ncolumns == B->ncolumns && A->ncolumns == output->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            output->data[i][j] = A->data[i][j] * B->data[i][j];
            output->mask[i][j] = MVISNAN_FUNC(output->data[i][j]) ? DATA_MISSING : DATA_PRESENT;
        }
    }
    return SUCCESS;
}

int mvMatElemDiv(mvMat *output, const mvMat *A, const mvMat *B)
{
    int i,j;
    if ( !(A->nrows == B->nrows && A->nrows == output->nrows &&
          A->ncolumns == B->ncolumns && A->ncolumns == output->ncolumns))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            if (B->data[i][j]==0.0)
            {
                output->data[i][j] = mvNaN();
                output->mask[i][j] = 1;
            }
            else
            {
                output->data[i][j] = A->data[i][j] / B->data[i][j];
                output->mask[i][j] = MVISNAN_FUNC(output->data[i][j]) ? DATA_MISSING : DATA_PRESENT;
            }
        }
    }
    return SUCCESS;
}

double mvVectorNorm(const mvMat *A)
{
    int i;
    double norm = 0.0;
    assert(A->nrows==1 || A->ncolumns==1);

    if (A->nrows>A->ncolumns)
    {
        for (i=0;i<A->nrows;i++)
        {
            if (A->mask[i][0]==1)
            {
                norm += A->data[i][0] * A->data[i][0];
            }
        }
    }
    else
    {
        for (i=0;i<A->ncolumns;i++)
        {
            if (A->mask[0][i]==1)
            {
                norm += A->data[0][i] * A->data[0][i];
            }
        }
    }
    return sqrt(norm);
}

int mvMatColumnSum(mvMat *output, const mvMat *A)
{
    int i,j, num_missing;
    if (output->ncolumns!=A->ncolumns)
    {
        return INCORRECT_DIMENSIONS;
    }

    for (j=0; j<A->ncolumns; j++)
    {
        num_missing=0;
        output->data[0][j]=0.0;
        for (i=0; i<A->nrows; i++)
        {
            if (A->mask[i][j]==DATA_PRESENT)
            {
                output->data[0][j] += A->data[i][j];
            }
        }
    }
    return SUCCESS;
}

int mvMatColumnSS(mvMat *output, const mvMat *A)
{
    int i,j;
    if (!(output->ncolumns == A->ncolumns && output->nrows == 1))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (j=0; j<A->ncolumns; j++)
    {
        output->data[0][j]=0.0;
        for (i=0; i<A->nrows; i++)
        {
            if (A->mask[i][j]==DATA_PRESENT)
            {
                output->data[0][j] += A->data[i][j] * A->data[i][j];
            }
        }
    }
    return SUCCESS;
}

int mvColumnMean(mvMat *output, const mvMat *A)
{
    int i,j, num_missing;
    if (output->ncolumns!=A->ncolumns)
    {
        return INCORRECT_DIMENSIONS;
    }

    for (j=0; j<A->ncolumns; j++)
    {
        num_missing=0;
        output->data[0][j]=0.0;
        for (i=0; i<A->nrows; i++)
        {
            if (A->mask[i][j]==DATA_MISSING)
            {
                num_missing++;
            }
            else
            {
                output->data[0][j] += A->data[i][j];
            }
        }
        if (num_missing==A->nrows)
        {
            output->data[0][j] = mvNaN();
            output->mask[0][j] = DATA_MISSING;
        }
        else
        {
            output->data[0][j] = output->data[0][j] / (A->nrows - num_missing);
        }
    }
    return SUCCESS;
}

int mvColumnVar(mvMat *output, const mvMat *A, int ddof)
{
    int i,j;
    int N, num_missing;
    double temp;
    mvMat *colmean;
    if (output->ncolumns!=A->ncolumns)
    {
        return INCORRECT_DIMENSIONS;
    }

    N = A->nrows - ddof;

    colmean = mvAllocMat(1, A->ncolumns);
    mvColumnMean(colmean, A);

    for (j=0; j<A->ncolumns; j++)
    {
        num_missing=0;
        output->data[0][j]=0.0;
        output->mask[0][j]=DATA_PRESENT;
        if ( MVISNAN_FUNC(colmean->data[0][j]) )
        {
            output->data[0][j] = mvNaN();
            output->mask[0][j] = DATA_MISSING;
            continue;
        }
        for (i=0; i<A->nrows; i++)
        {
            if (A->mask[i][j]==DATA_MISSING)
            {
                num_missing++;
            }
            else
            {
                temp = (A->data[i][j] - colmean->data[0][j]);
                temp = temp * temp;
                output->data[0][j] += temp;
            }
        }
        if ((N-num_missing) <= 0)
        {
            output->data[0][j] = mvNaN();
            output->mask[0][j] = DATA_MISSING;
        }
        else
        {
            output->data[0][j] = output->data[0][j] / ((double)(N-num_missing));
        }
    }

    mvFreeMat(&colmean);
    return SUCCESS;
}

int mvColumnStdDev(mvMat *output, const mvMat *A, int ddof)
{
    int j;
    int ret_val = mvColumnVar(output, A, ddof);
    if (ret_val)
        return ret_val;

    for (j=0; j<output->ncolumns; j++)
    {
        if (output->mask[0][j]==DATA_PRESENT)
        {

            output->data[0][j] = sqrt(output->data[0][j]);
        }
    }

    return SUCCESS;
}

int mvMatRowSS(mvMat *output, const mvMat *A)
{
    int i,j;
    if (!(output->nrows == A->nrows && output->ncolumns == 1))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<A->nrows; i++)
    {
        output->data[i][0]=0.0;
        for (j=0; j<A->ncolumns; j++)
        {
            if (A->mask[i][j]==DATA_PRESENT)
            {
                output->data[i][0] += A->data[i][j] * A->data[i][j];
            }
        }
    }
    return SUCCESS;
}

int mvMatColumnAdd(mvMat *output, const mvMat *A, const mvMat *columnValues)
{
    int i,j;
    if (!(output->nrows == A->nrows && output->ncolumns && A->ncolumns
            && columnValues->ncolumns == A->ncolumns && columnValues->nrows==1))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            if (columnValues->mask[0][j]==DATA_MISSING)
            {
                output->mask[i][j]=DATA_MISSING;
                output->data[i][j]=mvNaN();
            }
            else
            {
                output->data[i][j] = A->data[i][j] + columnValues->data[0][j];
                output->mask[i][j] = A->mask[i][j];
            }
        }
    }
    return SUCCESS;
}

int mvMatColumnSubtract(mvMat *output, const mvMat *A, const mvMat *columnValues)
{
    int i,j;
    if (!(output->nrows == A->nrows && output->ncolumns && A->ncolumns
            && columnValues->ncolumns == A->ncolumns && columnValues->nrows==1))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            if (columnValues->mask[0][j]==DATA_MISSING)
            {
                output->mask[i][j]=DATA_MISSING;
                output->data[i][j]=mvNaN();
            }
            else
            {
                output->data[i][j] = A->data[i][j] - columnValues->data[0][j];
                output->mask[i][j] = A->mask[i][j];
            }
        }
    }
    return SUCCESS;
}

int mvMatColumnMult(mvMat *output, const mvMat *A, const mvMat *columnValues)
{
    int i,j;
    if (!(output->nrows == A->nrows && output->ncolumns && A->ncolumns
            && columnValues->ncolumns == A->ncolumns && columnValues->nrows==1))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            if (columnValues->mask[0][j]==DATA_MISSING)
            {
                output->mask[i][j]=DATA_MISSING;
                output->data[i][j]=mvNaN();
            }
            else
            {
                output->data[i][j] = A->data[i][j] * columnValues->data[0][j];
                output->mask[i][j] = A->mask[i][j];
            }
        }
    }
    return SUCCESS;
}

int mvMatColumnDiv(mvMat *output, const mvMat *A, const mvMat *columnValues)
{
    int i,j;
    if (!(output->nrows == A->nrows && output->ncolumns == A->ncolumns
            && columnValues->ncolumns == A->ncolumns && columnValues->nrows==1))
    {
        return INCORRECT_DIMENSIONS;
    }

    for (i=0; i<output->nrows; i++)
    {
        for (j=0; j<output->ncolumns; j++)
        {
            if (columnValues->mask[0][j]==DATA_MISSING ||
                    columnValues->data[0][j]==0.0)
            {
                output->mask[i][j]=DATA_MISSING;
                output->data[i][j]=mvNaN();
            }
            else
            {
                output->data[i][j] = A->data[i][j] / columnValues->data[0][j];
                output->mask[i][j] = A->mask[i][j];
            }
        }
    }
    return SUCCESS;
}

int mvNumMissing(const mvMat *mat)
{
    int num_missing;
    int i,j;
    num_missing = 0;
    for (i=0; i<mat->nrows; i++)
    {
        for (j=0; j<mat->ncolumns; j++)
        {
            if (mat->mask[i][j]==DATA_MISSING)
                num_missing++;
        }
    }
    return num_missing;
}

double mvPctMissing(const mvMat *mat)
{
    double num_vals = (double) (mat->ncolumns * mat->nrows);
    double num_missing = mvNumMissing(mat);
    return num_missing/num_vals;

}

double mvDotProduct(const mvMat *a, const mvMat *b)
{
    int i, N;
    double result = 0.0;
    double *a_data, *b_data;
    size_t *a_mask, *b_mask;
    if ( !((a->nrows*a->ncolumns == b->nrows*b->ncolumns) &&
           (a->nrows==1 || a->ncolumns==1) && (b->nrows==1 || b->ncolumns) ))
    {
        return mvNaN();
    }
    N = a->nrows * a->ncolumns;
    a_data = a->data[0];
    b_data = b->data[0];
    a_mask = a->mask[0];
    b_mask = b->mask[0];

    // by now it should be a 1xN or Nx1 for each vector.
    for (i=0; i<N; i++)
    {
        if (a_mask[i]==DATA_PRESENT && b_mask[i]==DATA_PRESENT)
        {
            result += a_data[i] * b_data[i];
        }
    }
    return result;
}

double mvMatSS(const mvMat *A)
{
    double result = 0.0;
    if (!A->isReference)
    {
        double *a_data;
        size_t *a_mask;
        int i, N;
        N = A->nrows * A->ncolumns;

        a_data = A->data[0];
        a_mask = A->mask[0];

        for (i=0; i<N; i++)
        {
            if (a_mask[i]==DATA_PRESENT)
            {
                result += a_data[i] * a_data[i];
            }
        }
        return result;
    }
    else
    {
        int i,j;
        for (i=0; i<A->nrows; i++)
        {
            for (j=0; j<A->ncolumns; j++)
            {
                if(A->mask[i][j] == DATA_PRESENT)
                {
                    result += A->data[i][j] * A->data[i][j];
                }
            }
        }
        return result;
    }

}

double mvMatSum(const mvMat *A)
{
    double result = 0.0;
    if (!A->isReference)
    {
        double *a_data;
        size_t *a_mask;
        int i, N;
        N = A->nrows * A->ncolumns;

        a_data = A->data[0];
        a_mask = A->mask[0];

        for (i=0; i<N; i++)
        {
            if (a_mask[i]==DATA_PRESENT)
            {
                result += a_data[i];
            }
        }
        return result;
    }
    else
    {
        int i,j;
        for (i=0; i<A->nrows; i++)
        {
            for (j=0; j<A->ncolumns; j++)
            {
                if(A->mask[i][j] == DATA_PRESENT)
                {
                    result += A->data[i][j];
                }
            }
        }
        return result;
    }
}
