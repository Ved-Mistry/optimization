#include <omp.h>
#include <x86intrin.h>

#include "compute.h"

// Computes the dot product of vec1 and vec2, both of size n
int dot(uint32_t n, int32_t *vec1, int32_t *vec2) {
  // TODO: implement dot product of vec1 and vec2, both of size n
  __m256i vector1; 
  __m256i vector2; 
  __m256i prod_v = _mm256_setzero_si256();
  for (int i = 0; i < n/8 * 8; i += 8) { // Vectorized loop
        vector1 = _mm256_loadu_si256((__m256i *) (vec1 + i));
        vector2 = _mm256_loadu_si256((__m256i *) (vec2 + i)); 
        vector1 = _mm256_mullo_epi32(vector1, vector2); 
        prod_v = _mm256_add_epi32(prod_v, vector1); 
  }
  int result[8];
  _mm256_storeu_si256((__m256i *) result, prod_v);
  int sum = 0;
  sum += result[0] + result[1] + result[2] + result[3] + result[4] + result[5] + result[6] + result[7];
  for (int i = n/8*8; i < n; i++) { // Handle tail case
    sum += vec1[i]*vec2[i];
  }

  return sum; 
}

// Computes the convolution of two matrices
int convolve(matrix_t *a_matrix, matrix_t *b_matrix, matrix_t **output_matrix) {
  // TODO: convolve matrix a and matrix b, and store the resulting matrix in
  // output_matrix
      int32_t arow = a_matrix -> rows;
      int32_t acol = a_matrix -> cols;
      int32_t brow = b_matrix -> rows;
      int32_t bcol = b_matrix -> cols;
      int32_t* a = a_matrix -> data;
      int32_t* b = b_matrix -> data;
      int32_t n_b = brow * bcol;
 
      for (int32_t i = 0; i < n_b/2; i++) {
          int32_t temp = b[i];
          b[i] = b[n_b - i - 1];
          b[n_b - i - 1] = temp;
      }
 
      int32_t out_cols = acol - bcol + 1;
      int32_t out_rows = arow - brow + 1;
 
      *output_matrix = malloc(sizeof(matrix_t));
      int32_t* out = calloc(out_cols * out_rows, sizeof(int32_t));
 
      (*output_matrix) -> cols = out_cols;
      (*output_matrix) -> rows = out_rows;
      (*output_matrix) -> data = out;
      #pragma omp parallel for
      for (int i = 0; i < out_rows; i++) {
         for (int j = 0; j < out_cols; j++) {
              int x = 0;
             for (int k = 0; k < brow; k++) {
                  x += dot(bcol, a + (i + k)*acol + j, b + k*bcol);
              }
              out[i*out_cols+j] = x; 
          }
      }
     return 0;
}

// Executes a task
int execute_task(task_t *task) {
  matrix_t *a_matrix, *b_matrix, *output_matrix;

  if (read_matrix(get_a_matrix_path(task), &a_matrix))
    return -1;
  if (read_matrix(get_b_matrix_path(task), &b_matrix))
    return -1;

  if (convolve(a_matrix, b_matrix, &output_matrix))
    return -1;

  if (write_matrix(get_output_matrix_path(task), output_matrix))
    return -1;

  free(a_matrix->data);
  free(b_matrix->data);
  free(output_matrix->data);
  free(a_matrix);
  free(b_matrix);
  free(output_matrix);
  return 0;
}
