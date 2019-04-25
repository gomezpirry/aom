#include "operations.h"

int maxabs(double x[], int size) {

  int max_id = 0;
  double abs_value, max_value = pow(*x, 2);   /* use square which is quicker than absolute value */

  for (int k = 1; k < size; ++k) {
    abs_value = pow(x[k], 2);
    if (abs_value > max_value) {
      max_value = abs_value;
      max_id = k;
    }
  }
  return max_id;
}

void vec_assign(double y[], double x[], int ind[], int k) {

  for (int i = 0; i < k; ++i)
    y[i] = x[ind[i]];
}

void backsubst_L(double A[], double b[], double x[], int n, int k) {
  double rhs;

  for (int i = 0; i < k; ++i) {
    rhs = b[i];
    for (int j = 0; j<i; ++j) {
      rhs -= A[j*n + i] * x[j];
    }
    x[i] = rhs / A[i*n + i];
  }
}

void backsubst_Lt(double A[], double b[], double x[], int n, int k) {

  double rhs;

  for (int i = k; i >= 1; --i) {
    rhs = b[i - 1];
    for (int j = i; j<k; ++j) {
      rhs -= A[(i - 1)*n + j] * x[j];
    }
    x[i - 1] = rhs / A[(i - 1)*n + i - 1];
  }
}

void cholsolve(double A[], double b[], double x[], int n, int k) {
  double *tmp = (double*)malloc(k * sizeof(double));

  backsubst_L(A, b, tmp, n, k);
  backsubst_Lt(A, tmp, x, n, k);

  free(tmp);
}

void mat_vec(double alpha, double A[], double x[], double y[], int n, int m) {
  int i_n;
  double* Ax = (double*) calloc(n, sizeof(double));

  for (int i = 0; i < m; ++i) {
    i_n = i * n;
    for (int j = 0; j < n; ++j) {
      Ax[j] += A[i_n + j] * x[i];
    }
  }

  for (int j = 0; j<n; ++j) {
    y[j] = alpha*Ax[j];
  }

  free(Ax);
}

void vec_sum(double alpha, double x[], double y[], int n) {

  for (int i = 0; i<n; ++i) {
    y[i] += alpha * x[i];
  }
}

void matT_vec(double alpha, double A[], double x[], double y[], int n, int m) {
  int j, n_i;
  double sum0, sum1, sum2, sum3;

  for (j = 0; j<m; ++j) {
    y[j] = 0;
  }

  /* use loop unrolling to accelerate computation */

  for (int i = 0; i<m; ++i) {
    n_i = n*i;
    sum0 = sum1 = sum2 = sum3 = 0;
    for (j = 0; j + 4<n; j += 4) {
      sum0 += A[n_i + j] * x[j];
      sum1 += A[n_i + j + 1] * x[j + 1];
      sum2 += A[n_i + j + 2] * x[j + 2];
      sum3 += A[n_i + j + 3] * x[j + 3];
    }
    y[i] += alpha * ((sum0 + sum1) + (sum2 + sum3));
    while (j<n) {
      y[i] += alpha*A[n_i + j] * x[j];
      j++;
    }
  }
}

void mat_mat(double alpha, double A[], double B[], double X[], int n, int m, int k) {
  int i1, i2, i3, iX, iA, i2_n;
  double b;

  for (i1 = 0; i1<n*k; i1++) {
    X[i1] = 0;
  }

  for (i2 = 0; i2<m; ++i2) {
    i2_n = i2*n;
    iX = 0;
    for (i3 = 0; i3<k; ++i3) {
      iA = i2_n;
      b = B[i2 + i3*m];
      for (i1 = 0; i1<n; ++i1) {
        X[iX++] += A[iA++] * b;
      }
    }
  }

  for (i1 = 0; i1<n*k; i1++) {
    X[i1] *= alpha;
  }
}

void matT_mat(double alpha, double A[], double B[], double X[], int n, int m, int k) {
  int i1, i2, i3, iX, iA, i2_n;
  double *x, sum0, sum1, sum2, sum3;

  for (i2 = 0; i2<m; ++i2) {
    for (i3 = 0; i3<k; ++i3) {
      sum0 = sum1 = sum2 = sum3 = 0;
      for (i1 = 0; i1 + 4<n; i1 += 4) {
        sum0 += A[i1 + 0 + i2*n] * B[i1 + 0 + i3*n];
        sum1 += A[i1 + 1 + i2*n] * B[i1 + 1 + i3*n];
        sum2 += A[i1 + 2 + i2*n] * B[i1 + 2 + i3*n];
        sum3 += A[i1 + 3 + i2*n] * B[i1 + 3 + i3*n];
      }
      X[i2 + i3*m] = (sum0 + sum1) + (sum2 + sum3);
      while (i1<n) {
        X[i2 + i3*m] += A[i1 + i2*n] * B[i1 + i3*n];
        i1++;
      }
    }
  }

  for (i1 = 0; i1<m*k; i1++) {
    X[i1] *= alpha;
  }
}
