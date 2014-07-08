/*
* Filename: main.c
* Description: Test program for the libmva functions.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#define _USE_MATH_DEFINES
#include <windows.h>
#include <Psapi.h>
#endif

#include <mvconstants.h>
#include <mvmatrix.h>
#include <mvmodel.h>
#include <mvstats.h>
#include <mvpreprocess.h>

#include "data.h"

#ifdef WIN32
#define MVISNAN_FUNC _isnan
#else
#define MVISNAN_FUNC isnan
#endif

//int min(int a, int b) { return (a < b) ? a : b; }
//int max(int a, int b) { return (a > b) ? a : b; }

size_t get_memory_usage()
{
    size_t output = 0;
#ifdef WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    HANDLE hProcess = GetCurrentProcess();
    pmc.cb = sizeof(PROCESS_MEMORY_COUNTERS_EX);
    GetProcessMemoryInfo(hProcess, (PROCESS_MEMORY_COUNTERS*)&pmc, pmc.cb);
    output = pmc.PrivateUsage;
#else
    // nothin yet
#endif
    return output;
}

/* Two test functions that square and log the data.  For use in mvmat_row_func
  and mvmat_column_func */
double square(double x, void *opaque)
{
    (void) opaque;
    return x*x;
}

double log_missing(double x, void *opaque)
{
    (void) opaque;
    if (x <= 0.0)
    {
        return mv_NaN();
    }
    return log(x);
}

void mvmat_dump(MVMat *A)
{
    int i,j;
    for (i=0; i < A->nrows; i++)
    {
        if (i==0)
        {
            printf("[");
        }
        else
        {
            printf(" ");
        }
        for (j = 0; j < A->ncolumns; j++)
        {
            printf("%lf, ",A->data[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    int i,j,k;
    size_t memusage_before;
    size_t memusage_during;
    size_t memusage_after;
    MVMat *matrix1, *matrix2, *matrix3, *matrix4;
    char *data;
    (void)data;
    (void)matrix1;
    (void)matrix2;  // UNUSED
    (void)matrix3;
    (void)i; (void)j; (void)k;
    (void)argc; (void)argv;

    printf("****mvlib test suite****\n\n");

    printf("\n**Testing memory allocation and freeing**\n");
    memusage_before = get_memory_usage();
    matrix1 = mvmat_alloc(1024,1024);
    matrix2 = mvmat_allocz(1024,1024);
    matrix3 = mvmat_alloc_setval(1024, 1024, 3.1415962);
    matrix4 = mvmat_alloc_copy(matrix3);
    memusage_during = get_memory_usage();
    mvmat_free(&matrix1);
    mvmat_free(&matrix2);
    mvmat_free(&matrix3);
    mvmat_free(&matrix4);
    memusage_after = get_memory_usage();

    printf("\nMemory before allocation: %ld", memusage_before);
    printf("\nMemory after acllocation matrix: %ld.  Difference = %ld", memusage_during, (memusage_during-memusage_before));
    printf("\nMemory usage after mvFree: %ld.  Difference: %ld", memusage_after, (memusage_after-memusage_before));
    printf("\n");

    printf("\n**Testing setting, copying and getting data**\n");
    {
        matrix1 = mvmat_alloc(2,3);
        matrix2 = mvmat_alloc(2,3);
        for (i=0; i<2; i++)
        {
            for (j=0; j<3; j++)
            {
                mvmat_set_elem(matrix1,i,j, (double)i*j);
            }
        }
        mvmat_copy(matrix2, matrix1);
        for (i=0; i<2; i++)
        {
            for (j=0; j<3; j++)
            {
                double output;
                mvmat_get_elem(matrix2, &output, i,j);
                assert(output==(double)(i*j));
            }
        }
        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
    }

    printf("\n**Testing setting out of bounds**\n");
    {
        matrix1 = mvmat_alloc(2,3);
        assert(mvmat_set_elem(matrix1,3,0, M_PI)==INDEX_OUT_OF_BOUNDS);
        assert(mvmat_set_elem(matrix1,0,4, M_E)==INDEX_OUT_OF_BOUNDS);
        mvmat_free(&matrix1);
    }

    printf("\n**Testing getting out of bounds**\n");
    {
        double val;
        matrix1 = mvmat_alloc(4,5);
        assert(mvmat_get_elem(matrix1,&val, 0,6)==INDEX_OUT_OF_BOUNDS);
        assert(mvmat_get_elem(matrix1,&val,5,0)==INDEX_OUT_OF_BOUNDS);
        mvmat_free(&matrix1);

    }
    printf("\n**Testing mvAllocMatVal (and mvMatSet)\n");
    {
        int i,j;
        matrix1 = mvmat_alloc_setval(2,4, 1.00000);
        for (i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double value;
                mvmat_get_elem(matrix1, &value, i,j);
                assert(value == 1.000000);
            }
        }
        mvmat_free(&matrix1);
    }

    printf("\n**Testing concat rows with unmatching columns**\n");
    {
        matrix1 = mvmat_alloc_setval(2,3, M_PI);
        matrix2 = mvmat_alloc_setval(2,4, M_E);
        matrix3 = mvmat_alloc(4,3);
        assert(mvmat_concat_rows(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing concat rows**\n");
    {
        matrix1 = mvmat_alloc_setval(2,3, M_PI);
        matrix2 = mvmat_alloc_setval(2,3, M_E);
        matrix3 = mvmat_alloc(4,3);
        assert(mvmat_concat_rows(matrix3, matrix1, matrix2)==SUCCESS);
        for (i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double value;
                mvmat_get_elem(matrix3, &value, i,j);
                assert(value == M_PI);
            }
        }
        for (i=2; i<4; i++)
        {
            for (j=0; j<4; j++)
            {
                double value;
                mvmat_get_elem(matrix3, &value, i,j);
                assert(value == M_E);
            }
        }
        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
        mvmat_free(&matrix3);

    }

    printf("\n**Testing concat columns with unmatching rows**\n");
    {
        matrix1 = mvmat_alloc_setval(2,3, M_PI);
        matrix2 = mvmat_alloc_setval(3,4, M_E);
        matrix3 = mvmat_alloc(2,7);
        assert(mvmat_concat_columns(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing concat columns**\n");
    {
        matrix1 = mvmat_alloc_setval(2,3, M_PI);
        matrix2 = mvmat_alloc_setval(2,4, M_E);
        matrix3 = mvmat_alloc(2,7);
        assert(mvmat_concat_columns(matrix3, matrix1, matrix2)==SUCCESS);
        for (i=0; i<2; i++)
        {
            for (j=0; j<3; j++)
            {
                double value;
                mvmat_get_elem(matrix3, &value, i,j);
                assert(value == M_PI);
            }
        }
        for (i=0; i<2; i++)
        {
            for (j=3; j<7; j++)
            {
                double value;
                mvmat_get_elem(matrix3, &value, i,j);
                assert(value == M_E);
            }
        }
        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing mvTransposeMat with unmatching dimensions**\n");
    {
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc(4,3);

        assert(mvmat_transpose(matrix2, matrix1)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(3,2);
        assert(mvmat_transpose(matrix2, matrix1)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
    }

    printf("\n**Testing mvTransposeMat**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        for (j=0; j<4; j++)
        {
            mvmat_set_elem(matrix1, 0, j, M_E);
        }
        matrix2 = mvmat_alloc(4,2);
        assert(mvmat_transpose(matrix2, matrix1)==SUCCESS);
        for (i=0; i<4; i++)
        {
            for (j=0; j<2; j++)
            {
                double val;
                mvmat_get_elem(matrix2, &val, i, j);
                if (j==0)
                {
                    assert(val==M_E);
                }
                else if(j==1)
                {
                    assert(val==M_PI);
                }
            }
        }

        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
    }

    printf("\n**Testing mvAddMat with unmatching dimensions**\n");
    {
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc_setval(2, 5, 1.0);
        matrix3 = mvmat_alloc(2,4);

        assert(mvmat_add(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(3,4);

        assert(mvmat_add(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix3);
        matrix3 = mvmat_alloc(3,4);
        assert(mvmat_add(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing mvAddMat**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc_setval(2, 4, 1.0);
        matrix3 = mvmat_alloc(2,4);

        assert(mvmat_add(matrix3, matrix1, matrix2)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvmat_get_elem(matrix3,&val, i,j);
                assert(val == (M_PI + 1.0));
            }
        }
    }

    printf("\n**Testing mvSubtractMat with unmatching dimensions**\n");
    {
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc_setval(2, 5, 1.0);
        matrix3 = mvmat_alloc(2,4);

        assert(mvmat_subtract(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(3,4);

        assert(mvmat_subtract(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix3);
        matrix3 = mvmat_alloc(3,4);
        assert(mvmat_subtract(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing mvSubtractMat**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc_setval(2, 4, 1.0);
        matrix3 = mvmat_alloc(2,4);

        assert(mvmat_subtract(matrix3, matrix1, matrix2)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvmat_get_elem(matrix3,&val, i,j);
                assert(val == (M_PI - 1.0));
            }
        }
    }

    printf("\n**Testing mvAddMatS with unmatching dimensions**\n");
    {
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc(2,3);

        assert(mvmat_add_scalar(matrix2, matrix1, 1.0)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(3,4);
        assert(mvmat_add_scalar(matrix2, matrix1, 1.0)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
    }

    printf("\n**Testing mvAddMatS**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc(2,4);

        assert(mvmat_add_scalar(matrix2, matrix1, -1.0)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvmat_get_elem(matrix2,&val, i,j);
                assert(val == (M_PI - 1.0));
            }
        }
        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
    }

    printf("\n**Testing mvMatMult with incorrect dimensions.**\n");
    {
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc_setval(3, 3, 2.0);
        matrix3 = mvmat_alloc(2,3);

        assert(mvmat_mult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc_setval(4,3, 2.0);
        mvmat_free(&matrix3);
        matrix3 = mvmat_alloc(2,4);

        assert(mvmat_mult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
        mvmat_free(&matrix3);
    }
    printf("\n**Testing mvMatMult**\n");
    {
        int i;
        matrix1 = mvmat_alloc_setval(3, 3, 1.0);
        matrix2 = mvmat_alloc_setval(3, 1, 2.0);
        matrix3 = mvmat_alloc(3,1);

        assert(mvmat_mult(matrix3, matrix1, matrix2)==SUCCESS);

        for (i=0; i<3; i++)
        {
            double val;
            mvmat_get_elem(matrix3, &val, i, 0);
            assert(val == 6.0);
        }
    }
    printf("\n**Testing mvMultMatS with unmatching dimensions**\n");
    {
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc(2,3);

        assert(mvmat_mult_scalar(matrix2, matrix1, 2.0)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(3,4);
        assert(mvmat_mult_scalar(matrix2, matrix1, 2.0)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
    }

    printf("\n**Testing mvMultMatS**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc(2,4);

        assert(mvmat_mult_scalar(matrix2, matrix1, 2.0)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvmat_get_elem(matrix2,&val, i,j);
                assert(val == (M_PI * 2.0));
            }
        }
        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
    }

    printf("\n**Testing mvMatElemMult with unmatching dimensions**\n");
    {
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc_setval(2, 5, 1.0);
        matrix3 = mvmat_alloc(2,4);

        assert(mvmat_elem_mult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(3,4);

        assert(mvmat_elem_mult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix3);
        matrix3 = mvmat_alloc(3,4);
        assert(mvmat_elem_mult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing mvMatElemMult**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc_setval(2, 4, 2.0);
        matrix3 = mvmat_alloc(2,4);

        assert(mvmat_elem_mult(matrix3, matrix1, matrix2)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvmat_get_elem(matrix3,&val, i,j);
                assert(val == (M_PI * 2.0));
            }
        }
    }

    printf("\n**Testing mvMatElemDiv with unmatching dimensions**\n");
    {
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc_setval(2, 5, 1.0);
        matrix3 = mvmat_alloc(2,4);

        assert(mvmat_elem_div(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(3,4);

        assert(mvmat_elem_div(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix3);
        matrix3 = mvmat_alloc(3,4);
        assert(mvmat_elem_div(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing mvMatElemDiv**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc_setval(2, 4, M_PI);
        matrix2 = mvmat_alloc_setval(2, 4, 2.0);
        matrix3 = mvmat_alloc(2,4);

        assert(mvmat_elem_div(matrix3, matrix1, matrix2)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvmat_get_elem(matrix3,&val, i,j);
                assert(val == (M_PI / 2.0));
            }
        }
        mvmat_set(matrix2, 0.0);
        assert(mvmat_elem_div(matrix3,matrix1, matrix2)==SUCCESS);
        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvmat_get_elem(matrix3,&val, i,j);
                assert(MVISNAN_FUNC(val));
            }
        }
        mvmat_free(&matrix2);
        mvmat_free(&matrix1);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing vector norm**\n");
    {
        int i=0;
        double val;
        matrix1 = mvmat_alloc(4,1);
        for (i=0; i< 4; i++)
        {
            mvmat_set_elem(matrix1, i, 0, 1.0 + i);
        }
        val = mvmat_vector_norm(matrix1);
        assert(val == sqrt((1.0*1.0+2.0*2.0+3.0*3.0+4.0*4.0)));
        mvmat_set_elem(matrix1, 0,0, mv_NaN());
        val = mvmat_vector_norm(matrix1);
        assert(val == sqrt((2.0*2.0+3.0*3.0+4.0*4.0)));
        mvmat_set_elem(matrix1, 1,0, mv_NaN());
        val = mvmat_vector_norm(matrix1);
        assert(val == sqrt((3.0*3.0+4.0*4.0)));
        mvmat_set_elem(matrix1, 2,0, mv_NaN());
        val = mvmat_vector_norm(matrix1);
        assert(val == sqrt((4.0*4.0)));
        mvmat_set_elem(matrix1, 3,0, mv_NaN());
        val = mvmat_vector_norm(matrix1);
        assert(val == 0.0);
        mvmat_free(&matrix1);

        matrix1 = mvmat_alloc(1,4);
        for (i=0; i< 4; i++)
        {
            mvmat_set_elem(matrix1, 0, i, 1.0 + i);
        }
        val = mvmat_vector_norm(matrix1);
        assert(val == sqrt((1.0*1.0+2.0*2.0+3.0*3.0+4.0*4.0)));
        mvmat_set_elem(matrix1, 0,0, mv_NaN());
        val = mvmat_vector_norm(matrix1);
        assert(val == sqrt((2.0*2.0+3.0*3.0+4.0*4.0)));
        mvmat_set_elem(matrix1, 0,1, mv_NaN());
        val = mvmat_vector_norm(matrix1);
        assert(val == sqrt((3.0*3.0+4.0*4.0)));
        mvmat_set_elem(matrix1, 0,2, mv_NaN());
        val = mvmat_vector_norm(matrix1);
        assert(val == sqrt((4.0*4.0)));
        mvmat_set_elem(matrix1, 0,3, mv_NaN());
        val = mvmat_vector_norm(matrix1);
        assert(val == 0.0);
        mvmat_free(&matrix1);

    }

    printf("\n**Testing mvColumnMean**\n");
    {
        int i,j;
        matrix1=mvmat_alloc(3,4);
        matrix3=mvmat_alloc(1,4);    // known result;
        mvmat_set_elem(matrix3, 0, 0, 2.0);
        mvmat_set_elem(matrix3, 0, 1, 1.5);
        mvmat_set_elem(matrix3, 0, 2, 1.0);
        mvmat_set_elem(matrix3, 0, 3, mv_NaN());
        for(i=0;i<3; i++)
        {
            for(j=0;j<4;j++)
            {
                mvmat_set_elem(matrix1, i, j, (double)(i+1));
            }
        }

        mvmat_set_elem(matrix1, 2, 1, mv_NaN());
        mvmat_set_elem(matrix1, 1, 2, mv_NaN());
        mvmat_set_elem(matrix1, 2, 2, mv_NaN());
        mvmat_set_elem(matrix1, 0, 3, mv_NaN());
        mvmat_set_elem(matrix1, 1, 3, mv_NaN());
        mvmat_set_elem(matrix1, 2, 3, mv_NaN());

        matrix2=mvmat_alloc(1,5);
        assert(mvmat_column_mean(matrix2, matrix1)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        matrix2=mvmat_alloc(1,4);

        assert(mvmat_column_mean(matrix2, matrix1)==SUCCESS);
        for (j=0; j<4; j++)
        {
            double val1, val2;
            mvmat_get_elem(matrix2, &val1, 0, j);
            mvmat_get_elem(matrix3, &val2, 0, j);
            if (MVISNAN_FUNC(val2))
            {
                assert(MVISNAN_FUNC(val1));
            }
            else{
                assert(val1==val2);
            }

        }

        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing mvColumnVar**\n");
    {
        int i,j;
        matrix1=mvmat_alloc(3,4);
        matrix3=mvmat_alloc(1,4);    // known result;
        mvmat_set_elem(matrix3, 0, 0, 0.66666666666666663);
        mvmat_set_elem(matrix3, 0, 1, 0.25);
        mvmat_set_elem(matrix3, 0, 2, 0.0);
        mvmat_set_elem(matrix3, 0, 3, mv_NaN());
        for(i=0;i<3; i++)
        {
            for(j=0;j<4;j++)
            {
                mvmat_set_elem(matrix1, i, j, (double)(i+1));
            }
        }

        mvmat_set_elem(matrix1, 2, 1, mv_NaN());
        mvmat_set_elem(matrix1, 1, 2, mv_NaN());
        mvmat_set_elem(matrix1, 2, 2, mv_NaN());
        mvmat_set_elem(matrix1, 0, 3, mv_NaN());
        mvmat_set_elem(matrix1, 1, 3, mv_NaN());
        mvmat_set_elem(matrix1, 2, 3, mv_NaN());

        matrix2=mvmat_alloc(1,5);
        assert(mvmat_column_var(matrix2, matrix1, 0)==INCORRECT_DIMENSIONS);
        mvmat_free(&matrix2);
        matrix2=mvmat_alloc(1,4);

        assert(mvmat_column_var(matrix2, matrix1, 0)==SUCCESS);
        for (j=0; j<4; j++)
        {
            double val1, val2;
            mvmat_get_elem(matrix2, &val1, 0, j);
            mvmat_get_elem(matrix3, &val2, 0, j);
            if (MVISNAN_FUNC(val2))
            {
                assert(MVISNAN_FUNC(val1));
            }
            else{
                assert(val1==val2);
            }

        }

        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing mvMatColumnAdd**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc(4,2);
        matrix2 = mvmat_alloc(1,3);
        matrix3 = mvmat_alloc(4,2);
        for (i=0; i< 4; i++)
        {
            for (j=0; j<2; j++)
            {
                mvmat_set_elem(matrix1,i,j,(double)(i+1));
            }
        }
        assert(mvmat_column_add(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(1,2);
        mvmat_set_elem(matrix2, 0, 0, 2.0);
        mvmat_set_elem(matrix2, 0, 1, 4.0);
        assert(mvmat_column_add(matrix3, matrix1, matrix2)==SUCCESS);

        for (i=0; i<4; i++)
        {
            double val, mat1val;
            mvmat_get_elem(matrix3, &val, i, 0);
            mvmat_get_elem(matrix1, &mat1val, i, 0);
            assert(val==(mat1val + 2.0));
            mvmat_get_elem(matrix3, &val, i, 1);
            mvmat_get_elem(matrix1, &mat1val, i, 1);
            assert(val==(mat1val + 4.0));
        }
        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing mvMatColumnSubtract**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc(4,2);
        matrix2 = mvmat_alloc(1,3);
        matrix3 = mvmat_alloc(4,2);
        for (i=0; i< 4; i++)
        {
            for (j=0; j<2; j++)
            {
                mvmat_set_elem(matrix1,i,j,(double)(i+1));
            }
        }
        assert(mvmat_column_subtract(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(1,2);
        mvmat_set_elem(matrix2, 0, 0, 2.0);
        mvmat_set_elem(matrix2, 0, 1, 4.0);
        assert(mvmat_column_subtract(matrix3, matrix1, matrix2)==SUCCESS);

        for (i=0; i<4; i++)
        {
            double val, mat1val;
            mvmat_get_elem(matrix3, &val, i, 0);
            mvmat_get_elem(matrix1, &mat1val, i, 0);
            assert(val==(mat1val - 2.0));
            mvmat_get_elem(matrix3, &val, i, 1);
            mvmat_get_elem(matrix1, &mat1val, i, 1);
            assert(val==(mat1val - 4.0));
        }
        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing mvMatColumnMult**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc(4,2);
        matrix2 = mvmat_alloc(1,3);
        matrix3 = mvmat_alloc(4,2);
        for (i=0; i< 4; i++)
        {
            for (j=0; j<2; j++)
            {
                mvmat_set_elem(matrix1,i,j,(double)(i+1));
            }
        }
        assert(mvmat_column_mult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(1,2);
        mvmat_set_elem(matrix2, 0, 0, 2.0);
        mvmat_set_elem(matrix2, 0, 1, 4.0);
        assert(mvmat_column_mult(matrix3, matrix1, matrix2)==SUCCESS);

        for (i=0; i<4; i++)
        {
            double val, mat1val;
            mvmat_get_elem(matrix3, &val, i, 0);
            mvmat_get_elem(matrix1, &mat1val, i, 0);
            assert(val==(mat1val * 2.0));
            mvmat_get_elem(matrix3, &val, i, 1);
            mvmat_get_elem(matrix1, &mat1val, i, 1);
            assert(val==(mat1val * 4.0));
        }
        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
        mvmat_free(&matrix3);
    }
    printf("\n**Testing mvMatColumnDiv**\n");
    {
        int i,j;
        matrix1 = mvmat_alloc(4,2);
        matrix2 = mvmat_alloc(1,3);
        matrix3 = mvmat_alloc(4,2);
        for (i=0; i< 4; i++)
        {
            for (j=0; j<2; j++)
            {
                mvmat_set_elem(matrix1,i,j,(double)(i+1));
            }
        }
        assert(mvmat_column_div(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvmat_free(&matrix2);
        matrix2 = mvmat_alloc(1,2);
        mvmat_set_elem(matrix2, 0, 0, 2.0);
        mvmat_set_elem(matrix2, 0, 1, 0.0);
        assert(mvmat_column_div(matrix3, matrix1, matrix2)==SUCCESS);

        for (i=0; i<4; i++)
        {
            double val, mat1val;
            mvmat_get_elem(matrix3, &val, i, 0);
            mvmat_get_elem(matrix1, &mat1val, i, 0);
            assert(val==(mat1val / 2.0));
            mvmat_get_elem(matrix3, &val, i, 1);
            mvmat_get_elem(matrix1, &mat1val, i, 1);
            assert(MVISNAN_FUNC(val));
        }
        mvmat_free(&matrix1);
        mvmat_free(&matrix2);
        mvmat_free(&matrix3);
    }

    printf("\n**Testing mvNumMissing & mvPctMissing**\n");
    {
        matrix1 = mvmat_alloc(3,4);
        for(i=0;i<3; i++)
        {
            for(j=0;j<4;j++)
            {
                mvmat_set_elem(matrix1, i, j, (double)(i+1));
            }
        }

        mvmat_set_elem(matrix1, 2, 1, mv_NaN());
        mvmat_set_elem(matrix1, 1, 2, mv_NaN());
        mvmat_set_elem(matrix1, 2, 2, mv_NaN());
        mvmat_set_elem(matrix1, 0, 3, mv_NaN());
        mvmat_set_elem(matrix1, 1, 3, mv_NaN());
        mvmat_set_elem(matrix1, 2, 3, mv_NaN());

        assert(mvmat_pct_missing(matrix1)==6.0);
        assert(mvPctMissing(matrix1)==0.5);
        mvmat_free(&matrix1);
    }

    printf("\n**Testing mvRange **\n");
    {
        int i;
        MVMat * range;

        printf("\n\nmvRange(1, 5, -1);");
        range = mvmat_range(1, 5, -1);
        printf("\nPointer should be NULL - range = %p", range);
        assert(range==NULL);

        printf("\nmvRange(1, 5, 1);");
        range = mvmat_range(1, 5, 1);
        assert(range->nrows == 4 && range->ncolumns==1);
        for (i=0; i<range->nrows; i++)
        {
            printf("\nrange[%d] = %lf", i, range->data[i][0]);
        }
        mvmat_free(&range);

        printf("\n\nmvRange(1, 5, 3);");
        range = mvmat_range(1, 5, 3);
        assert(range->nrows == 2 && range->ncolumns == 1);
        for (i=0; i<range->nrows; i++)
        {
            printf("\nrange[%d] = %lf", i, range->data[i][0]);
        }
        mvmat_free(&range);

        printf("\n\nmvRange(5, -5, -2);");
        range = mvmat_range(5, -5, -2);
        assert(range->nrows == 5 && range->ncolumns == 1);
        for (i=0; i<range->nrows; i++)
        {
            printf("\nrange[%d] = %lf", i, range->data[i][0]);
        }
        mvmat_free(&range);
    }

    printf("\n\n**Testing LinSpace**\n");
    {
        int i;
        MVMat *linspace;

        printf("\nmvLinspace(1.0, 5.0, 0, 1);");
        linspace = mvmat_linspace(1.0, 5.0, 0, 1);
        printf("\nPointer should be NULL - linspace = %p", linspace);
        assert(linspace==NULL);

        printf("\n\nmvLinspace(1.0, 5.0, 4, 1)");
        linspace = mvmat_linspace(1.0, 5.0, 4, 1);
        assert(linspace->nrows == 4 && linspace->ncolumns == 1);
        for (i=0; i < linspace->nrows; i++)
        {
            printf("\nlinspace[%d] = %lf", i, linspace->data[i][0]);
        }
        mvmat_free(&linspace);

        printf("\n\nmvLinspace(1.0, 5.0, 10, 0)");
        linspace = mvmat_linspace(1.0, 5.0, 10, 0);
        assert(linspace->nrows == 10 && linspace->ncolumns == 1);
        for (i=0; i < linspace->nrows; i++)
        {
            printf("\nlinspace[%d] = %lf", i, linspace->data[i][0]);
        }
        mvmat_free(&linspace);

        printf("\n\nmvLinspace(1.0, -5.0, 4, 1)");
        linspace = mvmat_linspace(1.0, -5.0, 4, 1);
        assert(linspace->nrows == 4 && linspace->ncolumns == 1);
        for (i=0; i < linspace->nrows; i++)
        {
            printf("\nlinspace[%d] = %lf", i, linspace->data[i][0]);
        }
        mvmat_free(&linspace);
    }

    printf("\n\n**Testing mvMatSliceRowsRef && mvMatSliceRows**\n");
    {
        int i, j;
        const double x_data [FOODS_DATA_ROWS][FOODS_DATA_COLUMNS]=FOODS_DATA;

        MVMat *slice = mvmat_range(0, FOODS_DATA_ROWS, 7);
        MVMat *XSliceRef = NULL;
        MVMat *XSlice = mvmat_alloc(slice->nrows, FOODS_DATA_COLUMNS);
        MVMat *X = mvmat_alloc(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);

        for (i=0; i<FOODS_DATA_ROWS; i++)
        {
            for (j=0; j<FOODS_DATA_COLUMNS; j++)
            {
                mvmat_set_elem(X, i,j, x_data[i][j]);
            }
        }

        mvmat_slice_rows(XSlice, X, slice);
        mvmat_slice_rows_ref(&XSliceRef, X, slice);

        for (i=0; i<slice->nrows;i++)
        {
            int row = (int) slice->data[i][0];
            for (j=0; j<FOODS_DATA_COLUMNS; j++)
            {
                assert(X->data[row][j] == XSliceRef->data[i][j]);
                assert(X->data[row][j] == XSlice->data[i][j]);
                assert(XSlice->data[i][j] == XSliceRef->data[i][j]);
            }
        }

        assert(mvmat_ss(XSliceRef)==mvmat_ss(XSlice));
        mvmat_free(&XSliceRef);
        mvmat_free(&XSlice);
        mvmat_free(&slice);
        mvmat_free(&X);
    }

    printf("\n\n**Testing mvMatDeleteRowsRef**\n");
    {
        int i, j, round, num_rounds;
        const double x_data [FOODS_DATA_ROWS][FOODS_DATA_COLUMNS]=FOODS_DATA;

        MVMat *X = mvmat_alloc(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);



        for (i=0; i<FOODS_DATA_ROWS; i++)
        {
            for (j=0; j<FOODS_DATA_COLUMNS; j++)
            {
                mvmat_set_elem(X, i,j, x_data[i][j]);
            }
        }
        for (num_rounds = 0; num_rounds<FOODS_DATA_ROWS; num_rounds++)
        {
            for (round = 0; round<num_rounds; round++)
            {
                MVMat *slice = mvmat_range(round, FOODS_DATA_ROWS, num_rounds);
                MVMat *XSliceRef = NULL;


                mvmat_delete_rows_ref(&XSliceRef, X, slice);
                assert(XSliceRef->nrows == X->nrows - slice->nrows);

                for (i=0; i<XSliceRef->nrows; i++)
                {
                    double *dataXSlice = XSliceRef->data[i];
                    int row=0;
                    for (row = 0; row<slice->nrows; row++)
                    {
                        double *dataX = X->data[(int)slice->data[row][0]];
                        assert(dataX != dataXSlice);
                    }
                }
                mvmat_free(&XSliceRef);
                mvmat_free(&slice);
            }
        }
        mvmat_free(&X);
    }

    printf("\n**Testing mvmat_column_func**\n");
    {
        MVMAT_FUNC_PTR funcs[3] = {square, NULL, log_missing} ;
        void *opaques[3] = {NULL, NULL, NULL};

        printf("\nCol 1 func = square, Col 2 func = NULL, Col 3 func = log");

        MVMat * A = mvmat_alloc(3,3);
        mvmat_set_elem(A, 0, 0, 1.0);
        mvmat_set_elem(A, 0, 1, -2.0);
        mvmat_set_elem(A, 0, 2, 3.0);
        mvmat_set_elem(A, 1, 0, -4.0);
        mvmat_set_elem(A, 1, 1, 5.0);
        mvmat_set_elem(A, 1, 2, -6.0);
        mvmat_set_elem(A, 2, 0, 7.0);
        mvmat_set_elem(A, 2, 1, -8.0);
        mvmat_set_elem(A, 2, 2, 9.0);

        printf("\nIncoming A:\n");
        mvmat_dump(A);

        mvmat_column_func(A, A, funcs, opaques);

        printf("\nOutgoing A:\n");
        mvmat_dump(A);

        mvmat_free(&A);
    }

    printf("\n**Testing mvmat_row_func**\n");
    {
        MVMAT_FUNC_PTR funcs[3] = {square, NULL, log_missing} ;
        void *opaques[3] = {NULL, NULL, NULL};

        printf("\nRow 1 func = square, Row 2 func = NULL, Row 3 func = log");

        MVMat * A = mvmat_alloc(3,3);
        mvmat_set_elem(A, 0, 0, 1.0);
        mvmat_set_elem(A, 0, 1, -2.0);
        mvmat_set_elem(A, 0, 2, 3.0);
        mvmat_set_elem(A, 1, 0, -4.0);
        mvmat_set_elem(A, 1, 1, 5.0);
        mvmat_set_elem(A, 1, 2, -6.0);
        mvmat_set_elem(A, 2, 0, 7.0);
        mvmat_set_elem(A, 2, 1, -8.0);
        mvmat_set_elem(A, 2, 2, 9.0);

        printf("\nIncoming A:\n");
        mvmat_dump(A);

        mvmat_column_func(A, A, funcs, opaques);

        printf("\nOutgoing A:\n");
        mvmat_dump(A);

        mvmat_free(&A);
    }

    printf("\n**Testing preprocess do and undo **\n");
    {
        const double data [FOODS_DATA_ROWS][FOODS_DATA_COLUMNS]=FOODS_DATA;
        MVMat * X = mvmat_alloc(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);
        MVMat * X_pp = mvmat_alloc(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);
        MVMat * X_unpp = mvmat_alloc(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);
        MVMat * diff = mvmat_alloc(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);
        MVPreprocessContext *prepro = mvpreprocess_alloc_init_mat(FOODS_DATA_COLUMNS);
        MVPreprocessColumnInfo info;
        info.t_coeffs.A = 1.0;
        info.t_coeffs.B = 2.0;
        info.t_coeffs.C = 1.1;
        info.c = MV_CENTERING_MEAN;
        info.s = MV_SCALING_UV;
        info.multiplier = 1.0;
        int i, j;
        for (i=0; i<FOODS_DATA_ROWS; i++)
        {
            for (j=0; j<FOODS_DATA_COLUMNS; j++)
            {

                if (data[i][j]==FOODS_DATA_MASK)
                    mvmat_set_elem(X, i,j, mv_NaN());
                else
                    mvmat_set_elem(X, i,j, data[i][j]);
            }
        }

        // Set preprocess info for each column
        for (j = 0; j < FOODS_DATA_COLUMNS; j++)
        {
            info.t = j % 5;
            mvpreprocess_set_column(prepro, j, &info);
        }

        mvpreprocess_prep(prepro, X);
        mvpreprocess_do(prepro, X_pp, X);
        mvpreprocess_undo(prepro, X_unpp, X_pp);

        mvmat_subtract(diff, X, X_unpp);

        double diff_ss = mvmat_ss(diff);
        double ss_X = mvmat_ss(X);
        double ss_X_pp = mvmat_ss(X_pp);
        printf("\nX_pp[0][0] = %lf", X_pp->data[0][0]);

        printf("\n Diff result after preprocessing and undoing = %lf (%lf -> %lf)\n", diff_ss, ss_X, ss_X_pp);


        mvmat_free(&X);
        mvmat_free(&X_pp);
        mvmat_free(&X_unpp);
        mvpreprocess_free(&prepro);
    }

    printf("\n**Testing PCA**\n");
    {
        const double data [FOODS_DATA_ROWS][FOODS_DATA_COLUMNS]=FOODS_DATA;
        MVMat * X = mvmat_alloc(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);
        MVMat * X_mcuv = mvmat_alloc(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);
        MVMat *X_mean = mvmat_alloc(1,FOODS_DATA_COLUMNS);
        MVMat *X_std = mvmat_alloc(1,FOODS_DATA_COLUMNS);
        MVMat *new_t = NULL;
        MVMat *X_hat = mvmat_alloc(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);
        MVMat *HT2 = mvmat_allocz(X->nrows, 1);
        MVMat *SPE = mvmat_allocz(X->nrows, 1);
        MVModel *pca_model;
        int i, j;
        for (i=0; i<FOODS_DATA_ROWS; i++)
        {
            for (j=0; j<FOODS_DATA_COLUMNS; j++)
            {
                if (data[i][j]==FOODS_DATA_MASK)
                    mvmat_set_elem(X, i,j, mv_NaN());
                else
                    mvmat_set_elem(X, i,j, data[i][j]);
            }
        }


        mvmat_column_mean(X_mean, X);
        mvmat_column_stddev(X_std, X, 1);
        // mean center
        mvmat_column_subtract(X_mcuv, X, X_mean);
        // scale to unit variance.
        mvmat_column_div(X_mcuv, X_mcuv, X_std);

        pca_model = mvmodel_alloc_init_pca(X_mcuv);
        mvmodel_autofit(pca_model);

        printf("\nNumber of components computed: %d\nNumber of signification components: %d\n",
               pca_model->_A, pca_model->A);
        for(i=0; i< (int) pca_model->A; i++)
        {
            printf("PRESS[%2d] = %lf\n", i+1, pca_model->cvd->PRESS->data[i][0]);
        }

        for (i=0; i<FOODS_DATA_ROWS; i++)
        {
            for (j=0; j < pca_model->A; j++)
            {
                printf("T%d[%2d]=%1.8lf\t", j+1, i+1, pca_model->t->data[i][j]);
            }
            printf("\n");
        }

        printf("\nFoods R2X and Q2Cum values:\n");
        for (i = 0; i < pca_model->A; i++)
        {
            printf("R2X[%d]=%1.8lf\tQ2Cum[%d]=%1.8lf\n",
                  i+1, pca_model->R2X->data[i][0], i+1, pca_model->Q2cum->data[i][0]);
        }
        new_t = mvmat_allocz(X->nrows, pca_model->A);

        mvmodel_t_scores_from_obs(new_t, NULL, X_mcuv, pca_model, pca_model->A, MV_NEW_SCORE_SCP);

        printf("\nNew foods score:\n");
        for (i = 0; i < pca_model->A; i++)
        {
            printf("T%d[1] = %1.8lf\t", i+1, new_t->data[0][i]);
        }

        mvstats_ht2(HT2, pca_model->t, pca_model->t_stddev, 1, pca_model->A);

        printf("\nNew foods HT2 for [1-%d].", pca_model->A);

        for (i=0; i<FOODS_DATA_ROWS; i++)
        {
            printf("\nHT2[%d] = %1.8lf", i+1, HT2->data[i][0]);
        }

        int comp = 2;
        mvstats_spex_from_obs(SPE, pca_model, X_mcuv, pca_model->t, comp);
        printf("\nSPE Limits for Foods at A=%d. 0.95 = %1.8lf, 0.99 = %1.8lf",
               comp, mvstats_spe_limit(0.95, pca_model->SPEX, comp), mvstats_spe_limit(0.99, pca_model->SPEX, comp));
        for(i=0; i <FOODS_DATA_ROWS; i++)
        {
            printf("\nSPEFromObs[%d] = %1.8lf", i+1, SPE->data[i][0]);
        }

        for(i=0; i <FOODS_DATA_ROWS; i++)
        {
            printf("\nSPE[%d] = %1.8lf", i+1, pca_model->SPEX->data[i][comp-1]);
        }
        printf("\n");
        mvmat_free(&SPE);
        mvmat_free(&HT2);
        mvmat_free(&new_t);
        mvmat_free(&X);
        mvmat_free(&X_hat);
        mvmat_free(&X_mcuv);
        mvmat_free(&X_mean);
        mvmat_free(&X_std);
        mvmodel_free(&pca_model);
    }
    printf("\n**Testing PLS**\n");
    {
        int ioi = 207; // index of interest
        const double x_data [KAMYR_ROWS][KAMYR_X_COLUMNS]=KAMYR_X;
        const double y_data [KAMYR_ROWS][KAMYR_Y_COLUMNS]=KAMYR_Y;
        MVMat * X = mvmat_alloc(KAMYR_ROWS, KAMYR_X_COLUMNS);
        MVMat * X_mcuv = mvmat_alloc(KAMYR_ROWS, KAMYR_X_COLUMNS);
        MVMat *X_mean = mvmat_alloc(1,KAMYR_X_COLUMNS);
        MVMat *X_std = mvmat_alloc(1,KAMYR_X_COLUMNS);
        MVMat * Y = mvmat_alloc(KAMYR_ROWS, KAMYR_Y_COLUMNS);
        MVMat * Y_mcuv = mvmat_alloc(KAMYR_ROWS, KAMYR_Y_COLUMNS);
        MVMat *Y_mean = mvmat_alloc(1,KAMYR_Y_COLUMNS);
        MVMat *Y_std = mvmat_alloc(1,KAMYR_Y_COLUMNS);
        MVMat *E_pred = mvmat_allocz(KAMYR_ROWS, KAMYR_X_COLUMNS);
        MVMat *F_pred = mvmat_allocz(KAMYR_ROWS, KAMYR_Y_COLUMNS);
        MVMat *new_t = mvmat_allocz(X->nrows, 2);
        MVMat *new_u = mvmat_allocz(X->nrows, 2);
        MVModel *pls_model;
        int i, j;
        for (i=0; i<KAMYR_ROWS; i++)
        {
            for (j=0; j<KAMYR_X_COLUMNS; j++)
            {
                if (x_data[i][j]==KAMYR_MASK)
                    mvmat_set_elem(X, i,j, mv_NaN());
                else
                    mvmat_set_elem(X, i,j, x_data[i][j]);
            }
            for (j=0; j<KAMYR_Y_COLUMNS; j++)
            {
                if (y_data[i][j]==KAMYR_MASK)
                    mvmat_set_elem(Y, i,j, mv_NaN());
                else
                    mvmat_set_elem(Y, i,j, y_data[i][j]);
            }
        }


        mvmat_column_mean(X_mean, X);
        mvmat_column_stddev(X_std, X, 1);
        mvmat_column_mean(Y_mean, Y);
        mvmat_column_stddev(Y_std, Y, 1);
        // mean center
        mvmat_column_subtract(X_mcuv, X, X_mean);
        mvmat_column_subtract(Y_mcuv, Y, Y_mean);
        // scale to unit variance.
        mvmat_column_div(X_mcuv, X_mcuv, X_std);
        mvmat_column_div(Y_mcuv, Y_mcuv, Y_std);

        pls_model = mvmodel_alloc_init_pls(X_mcuv, Y_mcuv);
        mvmodel_add_component(pls_model);
        mvmodel_add_component(pls_model);
        mvmodel_add_component(pls_model);
        mvmodel_add_component(pls_model);

        printf("Kamyr Score for Obs. 1207: T1[%2d]=%1.8lf\tT2[%2d]=%1.8lf\n", ioi, pls_model->t->data[ioi][0],
                ioi, pls_model->t->data[ioi][1]);
        printf("Kamyr Score for Obs. 1207: U1[%2d]=%1.8lf\tU2[%2d]=%1.8lf\n", ioi, pls_model->u->data[ioi][0],
                ioi, pls_model->u->data[ioi][1]);
        ioi=6;
        printf("Kamyr W* for Variable Uczaa-3[6]: \nW*1[%d]=%1.8lf\nW*2[%d]=%1.8lf\nW*3[%d]=%1.8lf\nW*4[%d]=%1.8lf\n",
               ioi, pls_model->wstar->data[ioi][0],
               ioi, pls_model->wstar->data[ioi][1],
               ioi, pls_model->wstar->data[ioi][2],
               ioi, pls_model->wstar->data[ioi][3]);

        printf("\nKamyr R2X values:\nR2X[%d]=%1.8lf\nR2X[%d]=%1.8lf\nR2X[%d]=%1.8lf\nR2X[%d]=%1.8lf\n",
               1, pls_model->R2X->data[0][0],
               2, pls_model->R2X->data[1][0],
               3, pls_model->R2X->data[2][0],
               4, pls_model->R2X->data[3][0]);

        printf("\nKamyr R2Y values:\nR2Y[%d]=%1.8lf\nR2Y[%d]=%1.8lf\nR2Y[%d]=%1.8lf\nR2Y[%d]=%1.8lf\n",
               1, pls_model->R2Y->data[0][0],
               2, pls_model->R2Y->data[1][0],
               3, pls_model->R2Y->data[2][0],
               4, pls_model->R2Y->data[3][0]);

        printf("\nKamyr Q2cum values:\nQ2cum[%d]=%1.8lf\nQ2cum[%d]=%1.8lf\nQ2cum[%d]=%1.8lf\nQ2cum[%d]=%1.8lf\n",
               1, pls_model->Q2cum->data[0][0],
               2, pls_model->Q2cum->data[1][0],
               3, pls_model->Q2cum->data[2][0],
               4, pls_model->Q2cum->data[3][0]);

        mvmodel_t_scores_from_obs(new_t, E_pred, X_mcuv, pls_model, 2, MV_NEW_SCORE_SCP);


        ioi=207;
        printf("\nNew kamyr score T1[%2d] = %1.8lf, T2[%2d] = %1.8lf",
               ioi, new_t->data[ioi][0], ioi, new_t->data[ioi][1]);
        assert(fabs(new_t->data[ioi][0]-pls_model->t->data[ioi][0]) < MV_SQRT_EPS);
        assert(fabs(new_t->data[ioi][1]-pls_model->t->data[ioi][1]) < MV_SQRT_EPS);

        printf("\n\nR2X test %1.8lf, SSX = %lf, SSE0 = %lf, SSE1 = %lf, SSE2 = %lf, SSE_pred2 = %lf",
               1.0 - mvmat_ss(E_pred) / mvmat_ss(X_mcuv), mvmat_ss(X_mcuv),
               pls_model->SSX->data[0][0], pls_model->SSX->data[1][0], pls_model->SSX->data[2][0], mvmat_ss(E_pred));


        mvmodel_u_scores_from_obs(new_u, F_pred, Y_mcuv, new_t, pls_model, 2, MV_NEW_SCORE_SCP);
        printf("\nNew kamyr score U1[%2d] = %1.8lf, U2[%2d] = %1.8lf",
               ioi, new_u->data[ioi][0], ioi, new_u->data[ioi][1]);
        printf("\n\nR2Y test %1.8lf, SSY = %lf, SSF0 = %lf, SSF1 = %lf, SSF2 = %lf, SSF_pred2 = %lf",
               1.0 - mvmat_ss(F_pred) / mvmat_ss(Y_mcuv), mvmat_ss(Y_mcuv),
               pls_model->SSY->data[0][0], pls_model->SSY->data[1][0], pls_model->SSY->data[2][0], mvmat_ss(F_pred));

        mvmat_free(&F_pred);
        mvmat_free(&E_pred);
        mvmat_free(&new_u);
        mvmat_free(&new_t);
        mvmat_free(&X);
        mvmat_free(&X_mcuv);
        mvmat_free(&X_mean);
        mvmat_free(&X_std);
        mvmat_free(&Y);
        mvmat_free(&Y_mcuv);
        mvmat_free(&Y_mean);
        mvmat_free(&Y_std);
        mvmodel_free(&pls_model);
    }

    printf("\n**Testing F_ppf**\n");
    {
        double N = 16.0;
        double A = 2.0;
        double alpha = 0.95;
        //float(numComponents*(model.N-1)*(model.N+1) )/float((model.N*(model.N-numComponents)))
        printf("\nF_ppf(N1=2, N2=14, alpha=0.95) = %2.16lf", mvstats_F_ppf(0.95, 2, 14,0,0));
        printf("\nHT2(N=%d, A=%d, alpha=%lf) = %2.16lf", (int)N, (int)A, alpha,
               (A*(N-1)*(N+1) / (N*(N-A)))*mvstats_F_ppf(alpha, A, (N-A), 0, 0));
        printf("\nHT2(N=%d, A=%d, alpha=%lf) = %2.16lf", (int)N, (int)A, alpha,
               mvstats_ht2_limit(alpha, (int)A, (int)N));
    }

    printf("\n**Testing Chi2_ppf**\n");
    {
        double df = 10.0;
        double alpha = 0.95;
        printf("\nChi2_ppf(df = %lf, alpha=%lf) = %2.16lf", df, alpha, mvstats_chi2_ppf(alpha, df, 0, 0));
    }

    printf("\n**Testing Student's T ppf**\n");
    {
        double df = 10.2;
        double alpha =0.95;
        printf("\nStudent's T_ppf (df=%lf, alpha = %lf) = %2.16lf", df, alpha, mvstats_t_ppf(alpha,df,0,0));
    }

    printf("\n**Testing Gamma distribution**\n");
    {
        double df = 10.2;
        double alpha =0.95;
        printf("\nGamma (df=%lf, alpha = %lf) = %2.16lf", df, alpha, mvstats_gamma_ppf(alpha,df,0,0));
    }


    printf("\n\nDone!\n\n");
    return 0;
}

