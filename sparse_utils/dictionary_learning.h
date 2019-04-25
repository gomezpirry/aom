#ifndef DICTIONARY_LEARNING_H_
#define DICTIONARY_LEARNING_H_

#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#include "operations.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum Training_Alg {
  RANDOM,
  K_SVD
} Training;

struct SparseConfig {
  int sparse_coding;
  int k_codewords;
  int train_dict;
  int sparse_depth;
  int iterations;
  Training training;
};

typedef struct dictionary_t {
  double* data;
  unsigned int sample_size;
  unsigned int k;
} dictionary_t;

typedef struct sample_t {
  double* data;
  unsigned int sample_size;
  unsigned int number_blocks;
} sample_t;

typedef struct sparse_coeff_t {
  double* data;
  unsigned int k;
  unsigned int number_blocks;
} sparse_coeff_t;

/*!\brief Allocating storage for a dictionary
 *
 * \param[in]    dictionary    Pointer to storage blocks
 * \param[in]    sampleSize    block size of dictionary
 * \param[in]    k             number of codewords of dictionary
 *
 */
dictionary_t* dictionary_alloc(dictionary_t* dictionary, unsigned int sampleSize, unsigned int k);

/*!\brief Allocating storage for a samples
 *
 * \param[in]    samples          Pointer to storage sample blocks
 * \param[in]    block_size       block size of samples
 * \param[in]    width            image width
 * \param[in]    height           image height
 */
sample_t* samples_alloc(sample_t* samples, unsigned int sample_size, unsigned int number_samples);

/*!\brief Allocating storage for a sparse coefficients
 *
 * \param[in]    sparse_coeff     Pointer to storage sparse coefficients
 * \param[in]    k                block size of samples
 * \param[in]    number_samples   number of codewords of dictionary
 *
 */
sparse_coeff_t* sparse_coeff_alloc(sparse_coeff_t* sparse_coeff, unsigned int k, unsigned int number_blocks);

/*!\brief Select blocks of sample randomly
 *
 * \param[in]    samples    Pointer to samples blocks
 * \param[in]    dictionary Pointer to allocate dictionary
 * \param[in]    block_size Size of blocks
 *
 */
void init_random_dictionary(sample_t* samples, dictionary_t* dictionary, int block_size, int max_rand);

/*!\brief Batch OMP algorithm to calculate l0-norm
 *
 * \param[in]    iterations    number of iterations for algorithm
 *
 */
void  batch_OMP(int iterations, dictionary_t* dictionary, sample_t* samples, size_t number_samples, size_t sample_block, sparse_coeff_t* sparse_coefficients);

/*!\brief Select blocks of sample randomly
 *
 * \param[in]    array    Pointer to array to calculate norm
 * \param[in]    begin    first index
 * \param[in]    end      last index
 *
 */
double calculate_norm(double* array, int begin, int end);

void print_dictionary_csv(dictionary_t* dictionary, char* file_name);

void print_sparse_csv(sparse_coeff_t* sparse_coeff, char* file_name);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // DICTIONARY_LEARNING_H_
