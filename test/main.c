
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

int main(int argc, char *argv[])
{
    int i,j,k;
    size_t memusage_before;
    size_t memusage_during;
    size_t memusage_after;
    mvMat *matrix1, *matrix2, *matrix3, *matrix4;
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
    matrix1 = mvAllocMat(1024,1024);
    matrix2 = mvAllocMatZ(1024,1024);
    matrix3 = mvAllocMatVal(1024, 1024, 3.1415962);
    matrix4 = mvAllocMatCopy(matrix3);
    memusage_during = get_memory_usage();
    mvFreeMat(&matrix1);
    mvFreeMat(&matrix2);
    mvFreeMat(&matrix3);
    mvFreeMat(&matrix4);
    memusage_after = get_memory_usage();

    printf("\nMemory before allocation: %ld", memusage_before);
    printf("\nMemory after acllocation matrix: %ld.  Difference = %ld", memusage_during, (memusage_during-memusage_before));
    printf("\nMemory usage after mvFree: %ld.  Difference: %ld", memusage_after, (memusage_after-memusage_before));
    printf("\n");

    printf("\n**Testing setting, copying and getting data**\n");
    {
        matrix1 = mvAllocMat(2,3);
        matrix2 = mvAllocMat(2,3);
        for (i=0; i<2; i++)
        {
            for (j=0; j<3; j++)
            {
                mvMatSetElem(matrix1,i,j, (double)i*j);
            }
        }
        mvMatCopy(matrix2, matrix1);
        for (i=0; i<2; i++)
        {
            for (j=0; j<3; j++)
            {
                double output;
                mvMatGetElem(matrix2, &output, i,j);
                assert(output==(double)(i*j));
            }
        }
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
    }

    printf("\n**Testing setting out of bounds**\n");
    {
        matrix1 = mvAllocMat(2,3);
        assert(mvMatSetElem(matrix1,3,0, M_PI)==INDEX_OUT_OF_BOUNDS);
        assert(mvMatSetElem(matrix1,0,4, M_E)==INDEX_OUT_OF_BOUNDS);
        mvFreeMat(&matrix1);
    }

    printf("\n**Testing getting out of bounds**\n");
    {
        double val;
        matrix1 = mvAllocMat(4,5);
        assert(mvMatGetElem(matrix1,&val, 0,6)==INDEX_OUT_OF_BOUNDS);
        assert(mvMatGetElem(matrix1,&val,5,0)==INDEX_OUT_OF_BOUNDS);
        mvFreeMat(&matrix1);

    }
    printf("\n**Testing mvAllocMatVal (and mvMatSet)\n");
    {
        int i,j;
        matrix1 = mvAllocMatVal(2,4, 1.00000);
        for (i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double value;
                mvMatGetElem(matrix1, &value, i,j);
                assert(value == 1.000000);
            }
        }
        mvFreeMat(&matrix1);
    }

    printf("\n**Testing concat rows with unmatching columns**\n");
    {
        matrix1 = mvAllocMatVal(2,3, M_PI);
        matrix2 = mvAllocMatVal(2,4, M_E);
        matrix3 = mvAllocMat(4,3);
        assert(mvConcatRows(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing concat rows**\n");
    {
        matrix1 = mvAllocMatVal(2,3, M_PI);
        matrix2 = mvAllocMatVal(2,3, M_E);
        matrix3 = mvAllocMat(4,3);
        assert(mvConcatRows(matrix3, matrix1, matrix2)==SUCCESS);
        for (i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double value;
                mvMatGetElem(matrix3, &value, i,j);
                assert(value == M_PI);
            }
        }
        for (i=2; i<4; i++)
        {
            for (j=0; j<4; j++)
            {
                double value;
                mvMatGetElem(matrix3, &value, i,j);
                assert(value == M_E);
            }
        }
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix3);

    }

    printf("\n**Testing concat columns with unmatching rows**\n");
    {
        matrix1 = mvAllocMatVal(2,3, M_PI);
        matrix2 = mvAllocMatVal(3,4, M_E);
        matrix3 = mvAllocMat(2,7);
        assert(mvConcatColumns(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing concat columns**\n");
    {
        matrix1 = mvAllocMatVal(2,3, M_PI);
        matrix2 = mvAllocMatVal(2,4, M_E);
        matrix3 = mvAllocMat(2,7);
        assert(mvConcatColumns(matrix3, matrix1, matrix2)==SUCCESS);
        for (i=0; i<2; i++)
        {
            for (j=0; j<3; j++)
            {
                double value;
                mvMatGetElem(matrix3, &value, i,j);
                assert(value == M_PI);
            }
        }
        for (i=0; i<2; i++)
        {
            for (j=3; j<7; j++)
            {
                double value;
                mvMatGetElem(matrix3, &value, i,j);
                assert(value == M_E);
            }
        }
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing mvTransposeMat with unmatching dimensions**\n");
    {
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMat(4,3);

        assert(mvTransposeMat(matrix2, matrix1)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(3,2);
        assert(mvTransposeMat(matrix2, matrix1)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
    }

    printf("\n**Testing mvTransposeMat**\n");
    {
        int i,j;
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        for (j=0; j<4; j++)
        {
            mvMatSetElem(matrix1, 0, j, M_E);
        }
        matrix2 = mvAllocMat(4,2);
        assert(mvTransposeMat(matrix2, matrix1)==SUCCESS);
        for (i=0; i<4; i++)
        {
            for (j=0; j<2; j++)
            {
                double val;
                mvMatGetElem(matrix2, &val, i, j);
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

        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
    }

    printf("\n**Testing mvAddMat with unmatching dimensions**\n");
    {
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMatVal(2, 5, 1.0);
        matrix3 = mvAllocMat(2,4);

        assert(mvAddMat(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(3,4);

        assert(mvAddMat(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix3);
        matrix3 = mvAllocMat(3,4);
        assert(mvAddMat(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing mvAddMat**\n");
    {
        int i,j;
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMatVal(2, 4, 1.0);
        matrix3 = mvAllocMat(2,4);

        assert(mvAddMat(matrix3, matrix1, matrix2)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvMatGetElem(matrix3,&val, i,j);
                assert(val == (M_PI + 1.0));
            }
        }
    }

    printf("\n**Testing mvSubtractMat with unmatching dimensions**\n");
    {
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMatVal(2, 5, 1.0);
        matrix3 = mvAllocMat(2,4);

        assert(mvSubtractMat(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(3,4);

        assert(mvSubtractMat(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix3);
        matrix3 = mvAllocMat(3,4);
        assert(mvSubtractMat(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing mvSubtractMat**\n");
    {
        int i,j;
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMatVal(2, 4, 1.0);
        matrix3 = mvAllocMat(2,4);

        assert(mvSubtractMat(matrix3, matrix1, matrix2)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvMatGetElem(matrix3,&val, i,j);
                assert(val == (M_PI - 1.0));
            }
        }
    }

    printf("\n**Testing mvAddMatS with unmatching dimensions**\n");
    {
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMat(2,3);

        assert(mvAddMatS(matrix2, matrix1, 1.0)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(3,4);
        assert(mvAddMatS(matrix2, matrix1, 1.0)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
    }

    printf("\n**Testing mvAddMatS**\n");
    {
        int i,j;
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMat(2,4);

        assert(mvAddMatS(matrix2, matrix1, -1.0)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvMatGetElem(matrix2,&val, i,j);
                assert(val == (M_PI - 1.0));
            }
        }
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
    }

    printf("\n**Testing mvMatMult with incorrect dimensions.**\n");
    {
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMatVal(3, 3, 2.0);
        matrix3 = mvAllocMat(2,3);

        assert(mvMatMult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        matrix2 = mvAllocMatVal(4,3, 2.0);
        mvFreeMat(&matrix3);
        matrix3 = mvAllocMat(2,4);

        assert(mvMatMult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix3);
    }
    printf("\n**Testing mvMatMult**\n");
    {
        int i;
        matrix1 = mvAllocMatVal(3, 3, 1.0);
        matrix2 = mvAllocMatVal(3, 1, 2.0);
        matrix3 = mvAllocMat(3,1);

        assert(mvMatMult(matrix3, matrix1, matrix2)==SUCCESS);

        for (i=0; i<3; i++)
        {
            double val;
            mvMatGetElem(matrix3, &val, i, 0);
            assert(val == 6.0);
        }
    }
    printf("\n**Testing mvMultMatS with unmatching dimensions**\n");
    {
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMat(2,3);

        assert(mvMatMultS(matrix2, matrix1, 2.0)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(3,4);
        assert(mvMatMultS(matrix2, matrix1, 2.0)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
    }

    printf("\n**Testing mvMultMatS**\n");
    {
        int i,j;
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMat(2,4);

        assert(mvMatMultS(matrix2, matrix1, 2.0)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvMatGetElem(matrix2,&val, i,j);
                assert(val == (M_PI * 2.0));
            }
        }
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
    }

    printf("\n**Testing mvMatElemMult with unmatching dimensions**\n");
    {
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMatVal(2, 5, 1.0);
        matrix3 = mvAllocMat(2,4);

        assert(mvMatElemMult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(3,4);

        assert(mvMatElemMult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix3);
        matrix3 = mvAllocMat(3,4);
        assert(mvMatElemMult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing mvMatElemMult**\n");
    {
        int i,j;
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMatVal(2, 4, 2.0);
        matrix3 = mvAllocMat(2,4);

        assert(mvMatElemMult(matrix3, matrix1, matrix2)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvMatGetElem(matrix3,&val, i,j);
                assert(val == (M_PI * 2.0));
            }
        }
    }

    printf("\n**Testing mvMatElemDiv with unmatching dimensions**\n");
    {
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMatVal(2, 5, 1.0);
        matrix3 = mvAllocMat(2,4);

        assert(mvMatElemDiv(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(3,4);

        assert(mvMatElemDiv(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix3);
        matrix3 = mvAllocMat(3,4);
        assert(mvMatElemDiv(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing mvMatElemDiv**\n");
    {
        int i,j;
        matrix1 = mvAllocMatVal(2, 4, M_PI);
        matrix2 = mvAllocMatVal(2, 4, 2.0);
        matrix3 = mvAllocMat(2,4);

        assert(mvMatElemDiv(matrix3, matrix1, matrix2)==SUCCESS);

        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvMatGetElem(matrix3,&val, i,j);
                assert(val == (M_PI / 2.0));
            }
        }
        mvMatSet(matrix2, 0.0);
        assert(mvMatElemDiv(matrix3,matrix1, matrix2)==SUCCESS);
        for(i=0; i<2; i++)
        {
            for (j=0; j<4; j++)
            {
                double val;
                mvMatGetElem(matrix3,&val, i,j);
                assert(MVISNAN_FUNC(val));
            }
        }
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing vector norm**\n");
    {
        int i=0;
        double val;
        matrix1 = mvAllocMat(4,1);
        for (i=0; i< 4; i++)
        {
            mvMatSetElem(matrix1, i, 0, 1.0 + i);
        }
        val = mvVectorNorm(matrix1);
        assert(val == sqrt((1.0*1.0+2.0*2.0+3.0*3.0+4.0*4.0)));
        mvMatSetElem(matrix1, 0,0, mvNaN());
        val = mvVectorNorm(matrix1);
        assert(val == sqrt((2.0*2.0+3.0*3.0+4.0*4.0)));
        mvMatSetElem(matrix1, 1,0, mvNaN());
        val = mvVectorNorm(matrix1);
        assert(val == sqrt((3.0*3.0+4.0*4.0)));
        mvMatSetElem(matrix1, 2,0, mvNaN());
        val = mvVectorNorm(matrix1);
        assert(val == sqrt((4.0*4.0)));
        mvMatSetElem(matrix1, 3,0, mvNaN());
        val = mvVectorNorm(matrix1);
        assert(val == 0.0);
        mvFreeMat(&matrix1);

        matrix1 = mvAllocMat(1,4);
        for (i=0; i< 4; i++)
        {
            mvMatSetElem(matrix1, 0, i, 1.0 + i);
        }
        val = mvVectorNorm(matrix1);
        assert(val == sqrt((1.0*1.0+2.0*2.0+3.0*3.0+4.0*4.0)));
        mvMatSetElem(matrix1, 0,0, mvNaN());
        val = mvVectorNorm(matrix1);
        assert(val == sqrt((2.0*2.0+3.0*3.0+4.0*4.0)));
        mvMatSetElem(matrix1, 0,1, mvNaN());
        val = mvVectorNorm(matrix1);
        assert(val == sqrt((3.0*3.0+4.0*4.0)));
        mvMatSetElem(matrix1, 0,2, mvNaN());
        val = mvVectorNorm(matrix1);
        assert(val == sqrt((4.0*4.0)));
        mvMatSetElem(matrix1, 0,3, mvNaN());
        val = mvVectorNorm(matrix1);
        assert(val == 0.0);
        mvFreeMat(&matrix1);

    }

    printf("\n**Testing mvColumnMean**\n");
    {
        int i,j;
        matrix1=mvAllocMat(3,4);
        matrix3=mvAllocMat(1,4);    // known result;
        mvMatSetElem(matrix3, 0, 0, 2.0);
        mvMatSetElem(matrix3, 0, 1, 1.5);
        mvMatSetElem(matrix3, 0, 2, 1.0);
        mvMatSetElem(matrix3, 0, 3, mvNaN());
        for(i=0;i<3; i++)
        {
            for(j=0;j<4;j++)
            {
                mvMatSetElem(matrix1, i, j, (double)(i+1));
            }
        }

        mvMatSetElem(matrix1, 2, 1, mvNaN());
        mvMatSetElem(matrix1, 1, 2, mvNaN());
        mvMatSetElem(matrix1, 2, 2, mvNaN());
        mvMatSetElem(matrix1, 0, 3, mvNaN());
        mvMatSetElem(matrix1, 1, 3, mvNaN());
        mvMatSetElem(matrix1, 2, 3, mvNaN());

        matrix2=mvAllocMat(1,5);
        assert(mvColumnMean(matrix2, matrix1)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        matrix2=mvAllocMat(1,4);

        assert(mvColumnMean(matrix2, matrix1)==SUCCESS);
        for (j=0; j<4; j++)
        {
            double val1, val2;
            mvMatGetElem(matrix2, &val1, 0, j);
            mvMatGetElem(matrix3, &val2, 0, j);
            if (MVISNAN_FUNC(val2))
            {
                assert(MVISNAN_FUNC(val1));
            }
            else{
                assert(val1==val2);
            }

        }

        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing mvColumnVar**\n");
    {
        int i,j;
        matrix1=mvAllocMat(3,4);
        matrix3=mvAllocMat(1,4);    // known result;
        mvMatSetElem(matrix3, 0, 0, 0.66666666666666663);
        mvMatSetElem(matrix3, 0, 1, 0.25);
        mvMatSetElem(matrix3, 0, 2, 0.0);
        mvMatSetElem(matrix3, 0, 3, mvNaN());
        for(i=0;i<3; i++)
        {
            for(j=0;j<4;j++)
            {
                mvMatSetElem(matrix1, i, j, (double)(i+1));
            }
        }

        mvMatSetElem(matrix1, 2, 1, mvNaN());
        mvMatSetElem(matrix1, 1, 2, mvNaN());
        mvMatSetElem(matrix1, 2, 2, mvNaN());
        mvMatSetElem(matrix1, 0, 3, mvNaN());
        mvMatSetElem(matrix1, 1, 3, mvNaN());
        mvMatSetElem(matrix1, 2, 3, mvNaN());

        matrix2=mvAllocMat(1,5);
        assert(mvColumnVar(matrix2, matrix1, 0)==INCORRECT_DIMENSIONS);
        mvFreeMat(&matrix2);
        matrix2=mvAllocMat(1,4);

        assert(mvColumnVar(matrix2, matrix1, 0)==SUCCESS);
        for (j=0; j<4; j++)
        {
            double val1, val2;
            mvMatGetElem(matrix2, &val1, 0, j);
            mvMatGetElem(matrix3, &val2, 0, j);
            if (MVISNAN_FUNC(val2))
            {
                assert(MVISNAN_FUNC(val1));
            }
            else{
                assert(val1==val2);
            }

        }

        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing mvMatColumnAdd**\n");
    {
        int i,j;
        matrix1 = mvAllocMat(4,2);
        matrix2 = mvAllocMat(1,3);
        matrix3 = mvAllocMat(4,2);
        for (i=0; i< 4; i++)
        {
            for (j=0; j<2; j++)
            {
                mvMatSetElem(matrix1,i,j,(double)(i+1));
            }
        }
        assert(mvMatColumnAdd(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(1,2);
        mvMatSetElem(matrix2, 0, 0, 2.0);
        mvMatSetElem(matrix2, 0, 1, 4.0);
        assert(mvMatColumnAdd(matrix3, matrix1, matrix2)==SUCCESS);

        for (i=0; i<4; i++)
        {
            double val, mat1val;
            mvMatGetElem(matrix3, &val, i, 0);
            mvMatGetElem(matrix1, &mat1val, i, 0);
            assert(val==(mat1val + 2.0));
            mvMatGetElem(matrix3, &val, i, 1);
            mvMatGetElem(matrix1, &mat1val, i, 1);
            assert(val==(mat1val + 4.0));
        }
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing mvMatColumnSubtract**\n");
    {
        int i,j;
        matrix1 = mvAllocMat(4,2);
        matrix2 = mvAllocMat(1,3);
        matrix3 = mvAllocMat(4,2);
        for (i=0; i< 4; i++)
        {
            for (j=0; j<2; j++)
            {
                mvMatSetElem(matrix1,i,j,(double)(i+1));
            }
        }
        assert(mvMatColumnSubtract(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(1,2);
        mvMatSetElem(matrix2, 0, 0, 2.0);
        mvMatSetElem(matrix2, 0, 1, 4.0);
        assert(mvMatColumnSubtract(matrix3, matrix1, matrix2)==SUCCESS);

        for (i=0; i<4; i++)
        {
            double val, mat1val;
            mvMatGetElem(matrix3, &val, i, 0);
            mvMatGetElem(matrix1, &mat1val, i, 0);
            assert(val==(mat1val - 2.0));
            mvMatGetElem(matrix3, &val, i, 1);
            mvMatGetElem(matrix1, &mat1val, i, 1);
            assert(val==(mat1val - 4.0));
        }
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing mvMatColumnMult**\n");
    {
        int i,j;
        matrix1 = mvAllocMat(4,2);
        matrix2 = mvAllocMat(1,3);
        matrix3 = mvAllocMat(4,2);
        for (i=0; i< 4; i++)
        {
            for (j=0; j<2; j++)
            {
                mvMatSetElem(matrix1,i,j,(double)(i+1));
            }
        }
        assert(mvMatColumnMult(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(1,2);
        mvMatSetElem(matrix2, 0, 0, 2.0);
        mvMatSetElem(matrix2, 0, 1, 4.0);
        assert(mvMatColumnMult(matrix3, matrix1, matrix2)==SUCCESS);

        for (i=0; i<4; i++)
        {
            double val, mat1val;
            mvMatGetElem(matrix3, &val, i, 0);
            mvMatGetElem(matrix1, &mat1val, i, 0);
            assert(val==(mat1val * 2.0));
            mvMatGetElem(matrix3, &val, i, 1);
            mvMatGetElem(matrix1, &mat1val, i, 1);
            assert(val==(mat1val * 4.0));
        }
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix3);
    }
    printf("\n**Testing mvMatColumnDiv**\n");
    {
        int i,j;
        matrix1 = mvAllocMat(4,2);
        matrix2 = mvAllocMat(1,3);
        matrix3 = mvAllocMat(4,2);
        for (i=0; i< 4; i++)
        {
            for (j=0; j<2; j++)
            {
                mvMatSetElem(matrix1,i,j,(double)(i+1));
            }
        }
        assert(mvMatColumnDiv(matrix3, matrix1, matrix2)==INCORRECT_DIMENSIONS);

        mvFreeMat(&matrix2);
        matrix2 = mvAllocMat(1,2);
        mvMatSetElem(matrix2, 0, 0, 2.0);
        mvMatSetElem(matrix2, 0, 1, 0.0);
        assert(mvMatColumnDiv(matrix3, matrix1, matrix2)==SUCCESS);

        for (i=0; i<4; i++)
        {
            double val, mat1val;
            mvMatGetElem(matrix3, &val, i, 0);
            mvMatGetElem(matrix1, &mat1val, i, 0);
            assert(val==(mat1val / 2.0));
            mvMatGetElem(matrix3, &val, i, 1);
            mvMatGetElem(matrix1, &mat1val, i, 1);
            assert(MVISNAN_FUNC(val));
        }
        mvFreeMat(&matrix1);
        mvFreeMat(&matrix2);
        mvFreeMat(&matrix3);
    }

    printf("\n**Testing mvNumMissing & mvPctMissing**\n");
    {
        matrix1 = mvAllocMat(3,4);
        for(i=0;i<3; i++)
        {
            for(j=0;j<4;j++)
            {
                mvMatSetElem(matrix1, i, j, (double)(i+1));
            }
        }

        mvMatSetElem(matrix1, 2, 1, mvNaN());
        mvMatSetElem(matrix1, 1, 2, mvNaN());
        mvMatSetElem(matrix1, 2, 2, mvNaN());
        mvMatSetElem(matrix1, 0, 3, mvNaN());
        mvMatSetElem(matrix1, 1, 3, mvNaN());
        mvMatSetElem(matrix1, 2, 3, mvNaN());

        assert(mvNumMissing(matrix1)==6.0);
        assert(mvPctMissing(matrix1)==0.5);
        mvFreeMat(&matrix1);
    }

    printf("\n**Testing mvRange **\n");
    {
        int i;
        mvMat * range;

        printf("\n\nmvRange(1, 5, -1);");
        range = mvRange(1, 5, -1);
        printf("\nPointer should be NULL - range = %p", range);
        assert(range==NULL);

        printf("\nmvRange(1, 5, 1);");
        range = mvRange(1, 5, 1);
        assert(range->nrows == 4 && range->ncolumns==1);
        for (i=0; i<range->nrows; i++)
        {
            printf("\nrange[%d] = %lf", i, range->data[i][0]);
        }
        mvFreeMat(&range);

        printf("\n\nmvRange(1, 5, 3);");
        range = mvRange(1, 5, 3);
        assert(range->nrows == 2 && range->ncolumns == 1);
        for (i=0; i<range->nrows; i++)
        {
            printf("\nrange[%d] = %lf", i, range->data[i][0]);
        }
        mvFreeMat(&range);

        printf("\n\nmvRange(5, -5, -2);");
        range = mvRange(5, -5, -2);
        assert(range->nrows == 5 && range->ncolumns == 1);
        for (i=0; i<range->nrows; i++)
        {
            printf("\nrange[%d] = %lf", i, range->data[i][0]);
        }
        mvFreeMat(&range);
    }

    printf("\n\n**Testing LinSpace**\n");
    {
        int i;
        mvMat *linspace;

        printf("\nmvLinspace(1.0, 5.0, 0, 1);");
        linspace = mvLinspace(1.0, 5.0, 0, 1);
        printf("\nPointer should be NULL - linspace = %p", linspace);
        assert(linspace==NULL);

        printf("\n\nmvLinspace(1.0, 5.0, 4, 1)");
        linspace = mvLinspace(1.0, 5.0, 4, 1);
        assert(linspace->nrows == 4 && linspace->ncolumns == 1);
        for (i=0; i < linspace->nrows; i++)
        {
            printf("\nlinspace[%d] = %lf", i, linspace->data[i][0]);
        }
        mvFreeMat(&linspace);

        printf("\n\nmvLinspace(1.0, 5.0, 10, 0)");
        linspace = mvLinspace(1.0, 5.0, 10, 0);
        assert(linspace->nrows == 10 && linspace->ncolumns == 1);
        for (i=0; i < linspace->nrows; i++)
        {
            printf("\nlinspace[%d] = %lf", i, linspace->data[i][0]);
        }
        mvFreeMat(&linspace);

        printf("\n\nmvLinspace(1.0, -5.0, 4, 1)");
        linspace = mvLinspace(1.0, -5.0, 4, 1);
        assert(linspace->nrows == 4 && linspace->ncolumns == 1);
        for (i=0; i < linspace->nrows; i++)
        {
            printf("\nlinspace[%d] = %lf", i, linspace->data[i][0]);
        }
        mvFreeMat(&linspace);
    }

    printf("\n\n**Testing mvMatSliceRowsRef && mvMatSliceRows**\n");
    {
        int i, j;
        const double x_data [FOODS_DATA_ROWS][FOODS_DATA_COLUMNS]=FOODS_DATA;

        mvMat *slice = mvRange(0, FOODS_DATA_ROWS, 7);
        mvMat *XSliceRef = NULL;
        mvMat *XSlice = mvAllocMat(slice->nrows, FOODS_DATA_COLUMNS);
        mvMat *X = mvAllocMat(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);

        for (i=0; i<FOODS_DATA_ROWS; i++)
        {
            for (j=0; j<FOODS_DATA_COLUMNS; j++)
            {
                mvMatSetElem(X, i,j, x_data[i][j]);
            }
        }

        mvMatSliceRows(XSlice, X, slice);
        mvMatSliceRowsRef(&XSliceRef, X, slice);

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

        assert(mvMatSS(XSliceRef)==mvMatSS(XSlice));
        mvFreeMat(&XSliceRef);
        mvFreeMat(&XSlice);
        mvFreeMat(&slice);
        mvFreeMat(&X);
    }

    printf("\n\n**Testing mvMatDeleteRowsRef**\n");
    {
        int i, j, round, num_rounds;
        const double x_data [FOODS_DATA_ROWS][FOODS_DATA_COLUMNS]=FOODS_DATA;

        mvMat *X = mvAllocMat(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);



        for (i=0; i<FOODS_DATA_ROWS; i++)
        {
            for (j=0; j<FOODS_DATA_COLUMNS; j++)
            {
                mvMatSetElem(X, i,j, x_data[i][j]);
            }
        }
        for (num_rounds = 0; num_rounds<FOODS_DATA_ROWS; num_rounds++)
        {
            for (round = 0; round<num_rounds; round++)
            {
                mvMat *slice = mvRange(round, FOODS_DATA_ROWS, num_rounds);
                mvMat *XSliceRef = NULL;


                mvMatDeleteRowsRef(&XSliceRef, X, slice);
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
                mvFreeMat(&XSliceRef);
                mvFreeMat(&slice);
            }
        }
        mvFreeMat(&X);
    }

    printf("\n**Testing PCA**\n");
    {
        const double data [FOODS_DATA_ROWS][FOODS_DATA_COLUMNS]=FOODS_DATA;
        mvMat * X = mvAllocMat(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);
        mvMat * X_mcuv = mvAllocMat(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);
        mvMat *X_mean = mvAllocMat(1,FOODS_DATA_COLUMNS);
        mvMat *X_std = mvAllocMat(1,FOODS_DATA_COLUMNS);
        mvMat *new_t = NULL;
        mvMat *X_hat = mvAllocMat(FOODS_DATA_ROWS, FOODS_DATA_COLUMNS);
        mvMat *HT2 = mvAllocMatZ(X->nrows, 1);
        mvMat *SPE = mvAllocMatZ(X->nrows, 1);
        mvModel *pca_model;
        int i, j;
        for (i=0; i<FOODS_DATA_ROWS; i++)
        {
            for (j=0; j<FOODS_DATA_COLUMNS; j++)
            {
                if (data[i][j]==FOODS_DATA_MASK)
                    mvMatSetElem(X, i,j, mvNaN());
                else
                    mvMatSetElem(X, i,j, data[i][j]);
            }
        }


        mvColumnMean(X_mean, X);
        mvColumnStdDev(X_std, X, 1);
        // mean center
        mvMatColumnSubtract(X_mcuv, X, X_mean);
        // scale to unit variance.
        mvMatColumnDiv(X_mcuv, X_mcuv, X_std);

        pca_model = mvInitPCAModel(X_mcuv);
        mvAutoFit(pca_model);

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
        new_t = mvAllocMatZ(X->nrows, pca_model->A);

        mvNewObsT(new_t, NULL, X_mcuv, pca_model, pca_model->A, SCP);

        printf("\nNew foods score:\n");
        for (i = 0; i < pca_model->A; i++)
        {
            printf("T%d[1] = %1.8lf\t", i+1, new_t->data[0][i]);
        }

        mvHT2(HT2, pca_model->t, pca_model->t_stddev, 1, pca_model->A);

        printf("\nNew foods HT2 for [1-%d].", pca_model->A);

        for (i=0; i<FOODS_DATA_ROWS; i++)
        {
            printf("\nHT2[%d] = %1.8lf", i+1, HT2->data[i][0]);
        }
        mvSPE(SPE, pca_model->E);
        printf("\nSPE Limits for Foods at A=. 0.95 = %1.8lf, 0.99 = %1.8lf",
               mvSPELimit(0.95, SPE), mvSPELimit(0.99, SPE));
        for(i=0; i <FOODS_DATA_ROWS; i++)
        {
            printf("\nSPE[%d] = %1.8lf", i+1, SPE->data[i][0]);
        }
        printf("\n");
        mvFreeMat(&SPE);
        mvFreeMat(&HT2);
        mvFreeMat(&new_t);
        mvFreeMat(&X);
        mvFreeMat(&X_hat);
        mvFreeMat(&X_mcuv);
        mvFreeMat(&X_mean);
        mvFreeMat(&X_std);
        mvFreeModel(&pca_model);
    }
    printf("\n**Testing PLS**\n");
    {
        int ioi = 207; // index of interest
        const double x_data [KAMYR_ROWS][KAMYR_X_COLUMNS]=KAMYR_X;
        const double y_data [KAMYR_ROWS][KAMYR_Y_COLUMNS]=KAMYR_Y;
        mvMat * X = mvAllocMat(KAMYR_ROWS, KAMYR_X_COLUMNS);
        mvMat * X_mcuv = mvAllocMat(KAMYR_ROWS, KAMYR_X_COLUMNS);
        mvMat *X_mean = mvAllocMat(1,KAMYR_X_COLUMNS);
        mvMat *X_std = mvAllocMat(1,KAMYR_X_COLUMNS);
        mvMat * Y = mvAllocMat(KAMYR_ROWS, KAMYR_Y_COLUMNS);
        mvMat * Y_mcuv = mvAllocMat(KAMYR_ROWS, KAMYR_Y_COLUMNS);
        mvMat *Y_mean = mvAllocMat(1,KAMYR_Y_COLUMNS);
        mvMat *Y_std = mvAllocMat(1,KAMYR_Y_COLUMNS);
        mvMat *E_pred = mvAllocMatZ(KAMYR_ROWS, KAMYR_X_COLUMNS);
        mvMat *F_pred = mvAllocMatZ(KAMYR_ROWS, KAMYR_Y_COLUMNS);
        mvMat *new_t = mvAllocMatZ(X->nrows, 2);
        mvMat *new_u = mvAllocMatZ(X->nrows, 2);
        mvModel *pls_model;
        int i, j;
        for (i=0; i<KAMYR_ROWS; i++)
        {
            for (j=0; j<KAMYR_X_COLUMNS; j++)
            {
                if (x_data[i][j]==KAMYR_MASK)
                    mvMatSetElem(X, i,j, mvNaN());
                else
                    mvMatSetElem(X, i,j, x_data[i][j]);
            }
            for (j=0; j<KAMYR_Y_COLUMNS; j++)
            {
                if (y_data[i][j]==KAMYR_MASK)
                    mvMatSetElem(Y, i,j, mvNaN());
                else
                    mvMatSetElem(Y, i,j, y_data[i][j]);
            }
        }


        mvColumnMean(X_mean, X);
        mvColumnStdDev(X_std, X, 1);
        mvColumnMean(Y_mean, Y);
        mvColumnStdDev(Y_std, Y, 1);
        // mean center
        mvMatColumnSubtract(X_mcuv, X, X_mean);
        mvMatColumnSubtract(Y_mcuv, Y, Y_mean);
        // scale to unit variance.
        mvMatColumnDiv(X_mcuv, X_mcuv, X_std);
        mvMatColumnDiv(Y_mcuv, Y_mcuv, Y_std);

        pls_model = mvInitPLSModel(X_mcuv, Y_mcuv);
        mvModelAddComponent(pls_model);
        mvModelAddComponent(pls_model);
        mvModelAddComponent(pls_model);
        mvModelAddComponent(pls_model);

        printf("Kamyr Score for Obs. 1207: T1[%2d]=%1.8lf\tT2[%2d]=%1.8lf\n", ioi, pls_model->t->data[ioi][0],
                ioi, pls_model->t->data[ioi][1]);
        printf("Kamyr Score for Obs. 1207: U1[%2d]=%1.8lf\tU2[%2d]=%1.8lf\n", ioi, pls_model->u->data[ioi][0],
                ioi, pls_model->u->data[ioi][1]);
        ioi=6;
        printf("Kamyr W* for Variable Uczaa-3[6]: \nW*1[%d]=%1.8lf\nW*2[%d]=%1.8lf\nW*3[%d]=%1.8lf\nW*4[%d]=%1.8lf\n",
               ioi, pls_model->wStar->data[ioi][0],
               ioi, pls_model->wStar->data[ioi][1],
               ioi, pls_model->wStar->data[ioi][2],
               ioi, pls_model->wStar->data[ioi][3]);

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

        mvNewObsT(new_t, E_pred, X_mcuv, pls_model, 2, SCP);


        ioi=207;
        printf("\nNew kamyr score T1[%2d] = %1.8lf, T2[%2d] = %1.8lf",
               ioi, new_t->data[ioi][0], ioi, new_t->data[ioi][1]);
        assert(fabs(new_t->data[ioi][0]-pls_model->t->data[ioi][0]) < MV_DBL_SQRT_EPS);
        assert(fabs(new_t->data[ioi][1]-pls_model->t->data[ioi][1]) < MV_DBL_SQRT_EPS);

        printf("\n\nR2X test %1.8lf, SSX = %lf, SSE0 = %lf, SSE1 = %lf, SSE2 = %lf, SSE_pred2 = %lf",
               1.0 - mvMatSS(E_pred) / mvMatSS(X_mcuv), mvMatSS(X_mcuv),
               pls_model->SSX->data[0][0], pls_model->SSX->data[1][0], pls_model->SSX->data[2][0], mvMatSS(E_pred));


        mvNewObsU(new_u, F_pred, Y_mcuv, new_t, pls_model, 2, SCP);
        printf("\nNew kamyr score U1[%2d] = %1.8lf, U2[%2d] = %1.8lf",
               ioi, new_u->data[ioi][0], ioi, new_u->data[ioi][1]);
        printf("\n\nR2Y test %1.8lf, SSY = %lf, SSF0 = %lf, SSF1 = %lf, SSF2 = %lf, SSF_pred2 = %lf",
               1.0 - mvMatSS(F_pred) / mvMatSS(Y_mcuv), mvMatSS(Y_mcuv),
               pls_model->SSY->data[0][0], pls_model->SSY->data[1][0], pls_model->SSY->data[2][0], mvMatSS(F_pred));

        mvFreeMat(&F_pred);
        mvFreeMat(&E_pred);
        mvFreeMat(&new_u);
        mvFreeMat(&new_t);
        mvFreeMat(&X);
        mvFreeMat(&X_mcuv);
        mvFreeMat(&X_mean);
        mvFreeMat(&X_std);
        mvFreeMat(&Y);
        mvFreeMat(&Y_mcuv);
        mvFreeMat(&Y_mean);
        mvFreeMat(&Y_std);
        mvFreeModel(&pls_model);
    }

    printf("\n**Testing F_ppf**\n");
    {
        double N = 16.0;
        double A = 2.0;
        double alpha = 0.95;
        //float(numComponents*(model.N-1)*(model.N+1) )/float((model.N*(model.N-numComponents)))
        printf("\nF_ppf(N1=2, N2=14, alpha=0.95) = %2.16lf", F_ppf(0.95, 2, 14,0,0));
        printf("\nHT2(N=%d, A=%d, alpha=%lf) = %2.16lf", (int)N, (int)A, alpha,
               (A*(N-1)*(N+1) / (N*(N-A)))*F_ppf(alpha, A, (N-A), 0, 0));
        printf("\nHT2(N=%d, A=%d, alpha=%lf) = %2.16lf", (int)N, (int)A, alpha,
               mvHT2Limit(alpha, (int)A, (int)N));
    }

    printf("\n**Testing Chi2_ppf**\n");
    {
        double df = 10.0;
        double alpha = 0.95;
        printf("\nChi2_ppf(df = %lf, alpha=%lf) = %2.16lf", df, alpha, chi2_ppf(alpha, df, 0, 0));
    }

    printf("\n**Testing Student's T ppf**\n");
    {
        double df = 10.2;
        double alpha =0.95;
        printf("\nStudent's T_ppf (df=%lf, alpha = %lf) = %2.16lf", df, alpha, t_ppf(alpha,df,0,0));
    }

    printf("\n**Testing Gamma distribution**\n");
    {
        double df = 10.2;
        double alpha =0.95;
        printf("\nGamma (df=%lf, alpha = %lf) = %2.16lf", df, alpha, gamma_ppf(alpha,df,0,0));
    }


    printf("\n\nDone!\n\n");
    return 0;
}

