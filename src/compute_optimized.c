#include <omp.h>
#include <x86intrin.h>

#include "compute.h"

// Computes the dot product of vec1 and vec2, both of size n
int dot(uint32_t n, int32_t *vec1, int32_t *vec2) {
  // TODO: implement dot product of vec1 and vec2, both of size n
  int result[8];
  __m256i prod_v = __mm_set1_epi32(0);
  for (int i = 0; i < n/8 * 8; i += 8) { // Vectorized loop
        __m256i temp = __mm_mullo_epi32(__mm_loadu_si256((__m256i *) (vec1 + i)), __mm_loadu_si256((__m256i *) (vec2 + i)));
        prod_v = _mm_add_epi32(prod_v, temp); 
  }
  mm_storeu_si128((__m256i *) result, prod_v);
  for (int i = n/4 * 4; i < n; i++) { // Handle tail case
  result[0] *= vec1[i]*vec2[i];
  }
  return result[0] + result[1] + result[2] + result[3] + result[4] + result[5] + result[6] + result[7];
}

// Computes the convolution of two matrices
int convolve(matrix_t *a_matrix, matrix_t *b_matrix, matrix_t **output_matrix) {
  // TODO: convolve matrix a and matrix b, and store the resulting matrix in
  // output_matrix

  return -1;
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
