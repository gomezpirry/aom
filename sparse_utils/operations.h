#ifndef OPERATIONS_H_
#define OPERATIONS_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!\brief Get max absolute value in array
 *
 * \param[in]    x    array
 * \param[in]    size Array size
 *
 */
int  maxabs(double x[], int size);

/*!\brief assign values in array based on array of index
 *
 * \param[in]    y    array to assign
 * \param[in]    x    array with values to assign
 * \param[in]    ind  array of index
 * \param[in]    k    size of array
 */
void vec_assign(double y[], double x[], int ind[], int k);

/*!\brief Back substitution Lower Lx=b
 *
 * \param[in]    A    A matrix
 * \param[in]    b    b vector
 * \param[in]    x    x vector
 * \param[in]    n    row of matrix
 * \param[in]    k    columns of matrix
 *
 */
void backsubst_L(double A[], double b[], double x[], int n, int k);

/*!\brief Back substitution Lower Transpose Lt
 *
 * \param[in]    A    A matrix
 * \param[in]    b    b vector
 * \param[in]    x    x vector
 * \param[in]    n    row of matrix
 * \param[in]    k    columns of matrix
 *
 */
void backsubst_Lt(double A[], double b[], double x[], int n, int k);

/*!\brief Cholesky factorization A=LLt
 *
* \param[in]     A    A matrix
 * \param[in]    b    b vector
 * \param[in]    x    x vector
 * \param[in]    n    row of A matrix
 * \param[in]    k    columns of A matrix
 *
 */
void cholsolve(double A[], double b[], double x[], int n, int k);

/*!\brief (alpha)A*b = y
 *
 * \param[in]    alpha    scalar value for multiplication
 * \param[in]    A        A matrix
 * \param[in]    x        x vector
 * \param[in]    y        y vector
 * \param[in]    n        rows of A matrix
 * \param[in]    m        columns of A matrix
 *
 */
void mat_vec(double alpha, double A[], double x[], double y[], int n, int m);

/*!\brief (alpha)x = y
 *
 * \param[in]    alpha    scalar value for multiplication
 * \param[in]    x        x vector
 * \param[in]    y        y vector
 * \param[in]    n        size of x vector
 *
 */
void vec_sum(double alpha, double x[], double y[], int n);

/*!\brief (alpha)At*b = y
 *
 * \param[in]    alpha    scalar value for multiplication
 * \param[in]    A        A matrix transpose
 * \param[in]    x        x vector
 * \param[in]    y        y vector
 * \param[in]    n        rows of A matrix
 * \param[in]    m        columns of A matrix
 *
 */
void matT_vec(double alpha, double A[], double x[], double y[], int n, int m);

/*!\brief (alpha)A*B = X
 *
 * \param[in]    alpha    scalar value for multiplication
 * \param[in]    A        A matrix
 * \param[in]    B        B matrix
 * \param[in]    X        X matrix
 * \param[in]    n        rows of A matrix
 * \param[in]    m        columns of A matrix
 * \param[in]    k        columns of B matrix
 *
 */
void mat_mat(double alpha, double A[], double B[], double X[], int n, int m, int k);

/*!\brief (alpha)At*B = X
 *
 * \param[in]    alpha    scalar value for multiplication
 * \param[in]    A        A matrix transpose
 * \param[in]    B        B matrix
 * \param[in]    X        X matrix
 * \param[in]    n        rows of A matrix
 * \param[in]    m        columns of A matrix
 * \param[in]    k        columns of B matrix
 *
 */
void matT_mat(double alpha, double A[], double B[], double X[], int n, int m, int k);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // OPERATIONS_H_
