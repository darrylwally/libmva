/*! \file mvmatrix.h
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

#ifndef MVMATRIX_H
#define MVMATRIX_H

#include <stdlib.h>
#include "mvconstants.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Basic matrix data structure.
  This data structure can contain a mask of present/missing data which is the
  mask [optional].

  The mask, if set will be hold DATA_PRESENT for values that are present and
  DATA_MISSING for values that are missing.  Any computation on missing values will
  result in a missing value in the return location.

  rows and columns parameters determine the size of the matrix
  */
typedef struct MVMat_s {
    ///Array of data.  Will be allocated as 2d-array (matrix)
    double **data;

    /*! Array of missing data. */
    size_t **mask;

    /// Number of rows
    int nrows;

    /// Number of columns
    int ncolumns;

    /// Identifies if the data within this matrix is a reference.
    /*! If it is a reference, calling mvFreeMat will not delete the data.
      */
    int isReference;

} MVMat;

/*! Allocate an MVMat matrix structure.

  This structure can be freed with mvFreeMat;

  \arg rows the rows of the matrix.
  \arg columns the columns of the matrix.
  \return MVMat on success, NULL on failure

  \sa MVMat, mvFreeMat
  */
MVMat * mvAllocMat(int rows, int columns);

/*! Allocate an MVMat matrix structure and initialize all values to zero.

  \arg rows the rows of the matrix.
  \arg columns the columns of the matrix.
  \return MVMat on success, NULL on failure

  \sa MVMat, mvFreeMat
  */
MVMat * mvAllocMatZ(int rows, int columns);

/*! Allocate an MVMat matrix sturcture and initalize all values to val.

  \arg rows the rows of the matrix.
  \arg columns the columns of the matrix.
  \arg val the value of every element in the matrix.
  \return MVMat on success, NULL on failure

  \sa MVMat, mvFreeMat;
  */
MVMat *mvAllocMatVal(int rows, int columns, double val);

/*! Allocate an MVMat structure with the dimensions and contents of another
    effectively creating a copy.

  \arg other
  \return MVMat on success, NULL on failure

  \sa MVMat, mvFreeMat
  */
MVMat *mvAllocMatCopy(const MVMat *other);

/*! Generates a identity matrix of dimension NxN

  \arg dim The size of the matrix NxN where N > 0.
  \return MVMat * on success or NULL on fail
  \sa mvFreeMat
  */
MVMat *mvAllocIdentity(int dim);

/*! Generates a column vector with the range and step

  This function follows the Python style of the "range" function.

  \example To generate a vector from 1->5, set the following
  * start = 1
  * stop = 6
  * step = 1
  * result = [1,2,3,4,5]

  \example To generate a vector from 5->-5, i.e. count backwards, set the following
  * start = 5
  * stop = -6
  * step = -1
  * result = [5,4,3,2,1,0,-1,-2,-3,-4,-5]

  \note Depending on the step size, you may not stop at the "stop" value.  For example,
  a start=1, stop=5 and step = 3 would yield a vector with two values: [1, 4].

  \note This function will return NULL with ill-conditioned arguments.  For example,
  if start > stop and step > 0.  Will return NULL if step==0;

  \note Because MVMat data type is double, the contents will be double and if
  used for indexing they will need casting to integer types.

  \arg start starting value of the range.
  \arg stop ending value of the range.
  \arg step step size of the range.
  \return MVMat with Nx1 or NULL if error.  N is computed based on the step size and the range.
  \sa mvFreeMat
  */
MVMat *mvRange(int start, int stop, int step);

/*! Generates a column vector with evenly spaced values between the start and stop

  \example To generate a vector from 1->5 with 5 values including the end point
  set the following:
  * start = 1
  * stop = 5
  * num_values = 5
  * end_point = 1
  * result = [1,2,3,4,5]

  \example To generate a vector from 1->5 with 10 values excluding the end point
  set the following:
  * start = 1
  * stop = 5
  * num_values = 10
  * end_point = 0
  * result = [1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6]


  \example To generate a vector from 1->5 with 10 values including the end point
  set the following:
  * start = 1
  * stop = 5
  * num_values = 10
  * end_point = 1
  * result = [1.0, 1.44444, 1.88889, 2.33333, 2.77778, 3.22222, 3.66667, 4.11111,
              4.55556, 5.0]

  \arg start starting value
  \arg stop ending value
  \arg num_values number of values within the range.
  \arg end_point Include end point in the linear space (1=include, 0=don't include)
  \return MVMat with size (num_values x 1) or NULL if error.
  \sa mvFreeMat
 */
MVMat *mvLinspace(double start, double stop, int num_values, int end_point);

/*! Free an MVMat matrix structure.
  \arg mat pointer to the MvMat data structure.
  \return 0 on success, -1 on failure.
  */
int mvFreeMat(MVMat **mat);

/*! Copy the contents of other into output.

  \arg output The preallocated output matrix
  \arg other The matrix that will be copied into output.
  */
int mvMatCopy(MVMat *output, const MVMat *other);

/*! Concatenates, row-wise, two matrices together and put them in output.

  \arg output The pre-allocated output matrix of correct dimensions (N+M)xK
  \arg A The first matrix of size NxK
  \arg B The second matrix of size MxK
  \return 0 on success or MVErrorCode
  */
int mvConcatRows(MVMat *output, const MVMat *A, const MVMat *B);

/*! Concatenates, column-wise, two matrices together and put them in output.

  \arg output The pre-allocated output matrix of correct dimensions Nx(K+P)
  \arg A The first matrix of size NxK
  \arg B The second matrix of size NxP
  \return 0 on success or MVErrorCode
  */
int mvConcatColumns(MVMat *output, const MVMat *A, const MVMat *B);

/*! Performs row-wise slicing of A and stores the result in output.
  The "rows" vector which contains the rows from A that will be put into output.
  The rows from A will be copied into output in the order that they are listed
  in the "rows" vector

  See mvRange for conveniently creating a slice vector.

  \arg output The MVMat matrix of size(rows->nrows, A->ncolumns)
  \arg A The matrix that is being sliced.
  \arg rows MVMat that should be a column vector (Nx1) of row indexes.
  \return 0 on success or MVErrorCode
  \sa mvRange
  */
int mvMatSliceRows(MVMat *output, const MVMat *A, const MVMat *rows);

/*! Performs row-wise delete of A and stores the result in output.
  The "rows" vector which contains the rows that will be "deleted" from A.  The
  result is stored in output.
  The rows from A will be copied into output in the order that they are listed
  in the "rows" vector

  See mvRange for conveniently creating a slice vector.

  \arg output The MVMat matrix of size(rows->nrows, A->ncolumns)
  \arg A The matrix that is being sliced.
  \arg rows MVMat that should be a column vector (Nx1) of row indexes.
  \return 0 on SUCCESS or MVErrorCode
  \sa mvRange
  */
int mvMatDeleteRows(MVMat *output, const MVMat *A, const MVMat *rows);

/*! Performs row-wise slicing of A and stores the result in output.
  The "rows" vector which contains the rows from A that will be put into output.
  The rows from A will be referenced into output in the order that they are listed
  in the "rows" vector

  See mvRange for conveniently creating a slice vector.

  \arg output The double pointer MVMat matrix that will become an MVMat with
       is_reference = 1.  Can be NULL;
  \arg A The matrix that is being sliced.
  \arg rows MVMat that should be a column vector (Nx1) of row indexes.
  \return 0 on SUCCESS or MVErrorCode
  \sa mvRange
  */
int mvMatSliceRowsRef(MVMat **output, const MVMat *A, const MVMat *rows);

/*! Performs row-wise delete of A and stores the result in output.
  The "rows" vector which contains the rows that will be "deleted" from A.  The
  result is stored in output.
  The rows from A will be referenced into output in the order that they are listed
  in the "rows" vector

  See mvRange for conveniently creating a slice vector.

  \arg output The double pointer MVMat matrix that will become an MVMat with
       is_reference = 1.  Can be NULL
  \arg A The matrix that is being sliced.
  \arg rows MVMat that should be a column vector (Nx1) of row indexes.
  \return 0 on SUCCESS or MVErrorCode
  \sa mvRange
  */
int mvMatDeleteRowsRef(MVMat **output, const MVMat *A, const MVMat *rows);


/*! Performs column-wise slicing of A and stores the result in output.
  The "columns" vector which contains the columns from A that will be put into output.
  The columns from A will be copied into output in the order that they are listed
  in the "columns" vector

  See mvRange for conveniently creating a slice vector.

  \arg output The MVMat matrix of size(rows->nrows, A->ncolumns)
  \arg A The matrix that is being sliced.
  \arg columns MVMat that should be a column vector (Nx1) of column indexes.
  \return 0 on success or MVErrorCode
  \sa mvRange
  */
int mvMatSliceColumns(MVMat *output, const MVMat *A, const MVMat *columns);

/*! Convenience setter function
  Sets the value to every element in the matrix.

  \note To set every element in a matrix to "missing" pass in mvNaN() to the
  value.

  \arg mat The matrix that will be set.
  \arg value the value of each element.
  \return 0 on success or MVErrorCode
  \sa mvNaN
  */
int mvMatSet(MVMat *mat, double value);

/*! Convenience element setter function
  Sets the element in (row, column) of mat to value.  If the location of
  (row,column) is currently "missing" then it will be set to "present"; i.e.,
  its mask value will be set to TRUE if FALSE, unless value is NaN.

  \arg mat The MVMat structure.
  \arg row The row of the value we wish to set.
  \arg column The column of the value we wish to set.
  \arg value The value that will be set.
  \return 0 on success or a MVErrorCode.
  */
int mvMatSetElem(MVMat *mat, int row, int column, double value);

/*! Convenience element getter function
  Gets the element in (row, column) of mat to value.
  \arg mat The MVMat structure.
  \arg value Pointer to the where the value should be stored.  Will be NaN if
       value is missing or if index is out of bounds.
  \arg row The row of the value we wish to set.
  \arg column The column of the value we wish to set.
  \return 0 on success or MVErrorCode.
  */
int mvMatGetElem(const MVMat *mat, double *value, int row, int column);

/*! Performs the transponse of a matrix.
  This function takes the input matrix, mat, and transposes it into the output.

  \arg output Preallocated output matrix of correct dimensions.
  \arg mat The matrix that is to be transposed.
  \return 0 on success or a MVErrorCode.
  */

int mvTransposeMat(MVMat *output, const MVMat *mat);

/*! Performs addition of two matrices.
  This function will add two matrices, A and B, and store the result in output.

  \f[
    OUTPUT = A + B
  \f]

  \arg output The preallocated output matrix
  \arg A the first matrix
  \arg B the second matrix
  \return 0 on success or a MVErrorCode.
  */
int mvAddMat(MVMat *output, const MVMat *A, const MVMat *B);

/*! Performs subtraction of two matrices.
  This function will subtract two matrices, A and B, and store the result in
  output.

  \f[
    OUTPUT = A - B
  \f]

  \arg output The preallocated output matrix
  \arg A the first matrix
  \arg B the second matrix
  \return 0 on success or a MVErrorCode.
  */
int mvSubtractMat(MVMat *output, const MVMat *A, const MVMat *B);

/*! Performs addition a matrix with a scalar.
  This function will add matrix A with scalar and store the result in output.
  \f[
    OUTPUT = A + scalar
  \f]

  \note This function can operate on A in-place by specifying the output
  argument to be the same as A.

  \arg output The preallocated output matrix
  \arg A the first matrix
  \arg scalar the scalar value that will be subtracted from A
  \return 0 on success or a MVErrorCode.
  */
int mvAddMatS(MVMat *output, const MVMat *A, double scalar);

/*! Performs matrix multiplication
  This function performs matrix multiplication on A and B and stores the
  result in output.

  \f[
    OUTPUT = AB
  \f]

  The matrix dimensions must be appropriate; i.e, if A(m x p) then B(p x n)
  and the result must be Output(m x n)

  \arg output The preallocated output matrix
  \arg A the first matrix
  \arg B the second matrix
  \return 0 on success or a MVErrorCode.
  */
int mvMatMult(MVMat *output, const MVMat *A, const MVMat *B);

/*! Performs multiplication of a matrix with a scalar
    This function multiplies each element in matrix A by a scalar
    and stores the result in output.

    \f[
        OUTPUT = sA
    \f]
   If you wish to divide the matrix by a value, simply specify the scalar as
   1/scalar.

  \note This function can operate on A in-place by specifying the output
  argument to be the same as A.
  \arg output The preallocated output matrix
  \arg A the input matrix
  \arg scalar the scalar
  \return 0 on success or a MVErrorCode.
  */
int mvMatMultS(MVMat *output, const MVMat *A, double scalar);

/*! Performs element-wise matrix multiplication on matrices of similar dimension

  \note This function can operate on A in-place by specifying the output
  argument to be the same as A.
  \arg output The preallocated output matrix
  \arg A the first matrix
  \arg B the second matrix
  \return 0 on success or a MVErrorCode.
  */
int mvMatElemMult(MVMat *output, const MVMat *A, const MVMat *B);

/*! Performs element-wise matrix division on matrices of similar dimension

  Division by zero will result in the element becoming NaN and
  its mask element will be set to false.

  \note This function can operate on A in-place by specifying the output
  argument to be the same as A.

  \arg output The preallocated output matrix
  \arg A the first matrix
  \arg B the second matrix
  \return 0 on success or a MVErrorCode.
  */
int mvMatElemDiv(MVMat *output, const MVMat *A, const MVMat *B);

/*! Computes the norm of a vector
  A vector is supplied as an MVMat.  The vector must be Nx1 or 1xN in dimension.

  \note This function asserts that one of the dimensions is 1.

  \arg A MVMat that is Nx1 or 1xN in dimension.
  \return The vector's Euclidiean norm. An entirely missing vector will return 0.
  */
double mvVectorNorm(const MVMat * A);

/*! Computes the sum of each column in a matrix.
  The sum of each column is computed and stored in output.

  \note If A is NxM in size, then output must be 1xM in size.
  \arg output The output MVMat that is 1xM in dimension.
  \arg A the input matrix.
  \return 0 on success or a MVErrorCode.
  */
int mvMatColumnSum(MVMat *output, const MVMat * A);

/*! Computes the sum of squares for each column in a matrix.
  The sum of squares for each column is computed and stored in output.

  \note If A is NxM in size, then output must be 1xM in size.
  \arg output The output MVMat that is 1xM in dimension.
  \arg A the input matrix.
  \return 0 on success or a MVErrorCode.
  */
int mvMatColumnSS(MVMat *output, const MVMat *A);

/*! Computes the mean of each column in a matrix.
  The mean of each column is computed and stored in output.

  \note If A is NxM in size, then output must be 1xM in size.
  \arg output The output MVMat that is 1xM in dimension.
  \arg A the input matrix.
  \return 0 on success or a MVErrorCode.
  */
int mvColumnMean(MVMat * output, const MVMat * A);

/*! Computes the variance of each column in a matrix.
  The variance of each column is computed and stored in output.

  \note If A is NxM in size, then output must be 1xM in size.
  \arg output The output MVMat that is 1xM in dimension.
  \arg A the input matrix.
  \arg ddof N-ddof in the denominator.  Default should be zero or one.
  \return 0 on success or a MVErrorCode.
  */
int mvColumnVar(MVMat * output, const MVMat * A, int ddof);

/*! Computes the standard deviation of each column in a matrix.
  The standard deviation of each column is computed and stored in output.

  \note If A is NxM in size, then output must be 1xM in size.
  \arg output The output MVMat that is 1xM in dimension.
  \arg A the input matrix.
  \arg ddof N-ddof in the denominator.  Default should be zero or one.
  \return 0 on success or a MVErrorCode.
  */
int mvColumnStdDev(MVMat *output, const MVMat *A, int ddof);

/*! Computes the sum of squares for each row in a matrix.
  The sum of squares for each row is computed and stored in output.

  \note If A is NxM in size, then output must be Nx1 in size.
  \arg output The output MVMat that is Nx1 in dimension.
  \arg A the input matrix.
  \return 0 on success or a MVErrorCode.
  */
int mvMatRowSS(MVMat *output, const MVMat *A);

/*! Performs a column-wise addition of elements in column with elements of A.

  \note This function can operate on A in-place by specifying the output
  argument to be the same as A.
  \arg output The output matrix
  \arg A the input matrix.
  \arg columnValues The values that will be added to each column of A.
  \return 0 on success or a MVErrorCode.
  */
int mvMatColumnAdd(MVMat *output, const MVMat *A, const MVMat *columnValues);

/*! Performs a column-wise subtraction of elements in column from elements of A.

  \note This function can operate on A in-place by specifying the output
  argument to be the same as A.
  \arg output The output matrix
  \arg A the input matrix.
  \arg columnValues The values that will be subtracted from each column of A.
  \return 0 on success or a MVErrorCode.
  */
int mvMatColumnSubtract(MVMat *output, const MVMat *A, const MVMat *columnValues);

/*! Performs a column-wise multiplication of elements in column with elements of A.

  \note This function can operate on A in-place by specifying the output
  argument to be the same as A.
  \arg output The output matrix
  \arg A the input matrix.
  \arg columnValues The values that will be multiplied with each column of A.
  \return 0 on success or a MVErrorCode.
  */
int mvMatColumnMult(MVMat *output, const MVMat *A, const MVMat *columnValues);


/*! Performs a column-wise division of elements in column with elements of A.

  \note This function can operate on A in-place by specifying the output
  argument to be the same as A.
  \arg output The output matrix
  \arg A the input matrix.
  \arg columnValues The values that will be divided with each column of A.
  \return 0 on success or a MVErrorCode.
  */
int mvMatColumnDiv(MVMat *output, const MVMat *A, const MVMat *columnValues);

/*! Returns the number of missing values in the matrix
  \arg mat The MVMat matrix.
  \return number of missing values
  */
int mvNumMissing(const MVMat *mat);

/*! Returns the percentage of missing values in the matrix
  \arg mat The MVMat matrix.
  \return percentage of missing values (between 0 -> 1.0)
  */
double mvPctMissing(const MVMat *mat);

/*! Returns the dot product of vector's A and B
  Vectors must be 1 dimensional.  As long as 1 dimension is size(1) and the
  other dimension is size(N) in both vectors the dot product will be computed.

  \note if all values are missing then the result will be zero.
  \arg a the first vector of size 1xN
  \arg b the second vector of size Nx1.
  \arg the dot product or NaN if error.
  */
double mvDotProduct(const MVMat *a, const MVMat *b);

/*! Returns the SUM of SQUARES of an MVMat

  \note if all values are missing then the result will be zero.
  \arg A the matrix
  \return sum of squares;
  */
double mvMatSS(const MVMat *A);

/*! Returns the sum of all of the elements of an MVMat

  This function properly accounts for missing data.

  \note if all values are missing then the result will be zero.
  \arg A the matrix.
  \return sum
  */
double mvMatSum(const MVMat *A);

#ifdef __cplusplus
}
#endif

#endif // MVMATRIX_H
