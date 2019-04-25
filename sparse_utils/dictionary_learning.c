#include "dictionary_learning.h"

dictionary_t* dictionary_alloc(dictionary_t* dictionary, unsigned int sample_size, unsigned int k) {

  memset(dictionary, 0, sizeof(dictionary_t));

  dictionary->data = (double*)malloc(sizeof(double) * (sample_size * k));

  dictionary->sample_size = sample_size;
  dictionary->k = k;

  return dictionary;
}

sample_t* samples_alloc(sample_t* samples, unsigned int sample_size, unsigned int number_samples) {
  memset(samples, 0, sizeof(sample_t));

  samples->sample_size = sample_size;
  samples->number_blocks = number_samples;
  samples->data = (double*)malloc(sample_size * number_samples * sizeof(double));

  return samples;
}

sparse_coeff_t* sparse_coeff_alloc(sparse_coeff_t* sparse_coeff, unsigned int k, unsigned int number_blocks) {
  memset(sparse_coeff, 0, sizeof(sparse_coeff_t));

  sparse_coeff->data = (double*)calloc(k * number_blocks, sizeof(double));
  sparse_coeff->k = k;
  sparse_coeff->number_blocks = number_blocks;

  return sparse_coeff;
}

void init_random_dictionary(sample_t* samples, dictionary_t* dictionary, int all_block_size, int max_rand) {

  srand(time(NULL));

  if(all_block_size != (int)dictionary->sample_size) {
    fprintf(stderr, "Sample size and sample size of dictionary are not equal\n");
    exit(EXIT_FAILURE);
  }

  for(int i = 0; i < (int)dictionary->k; i++){
    int index = rand() %  max_rand;
    int index_offset = index * all_block_size;
    double norm = calculate_norm(samples->data, index_offset, index_offset + all_block_size);
    for(int j = index_offset, k = 0; j < index_offset + all_block_size; j++, k++) {
      dictionary->data[i * all_block_size + k] = samples->data[j] / norm;
    }
  }
}

void batch_OMP(int iterations, dictionary_t* dictionary, sample_t* samples, size_t number_samples, size_t sample_block, sparse_coeff_t* sparse_coefficients) {

  assert(dictionary->sample_size == sample_block);

  // allocate aux arrays
  int k_codewords   = dictionary->k;

  double* G         = (double*)calloc(k_codewords * k_codewords, sizeof(double));
  double* DtX       = (double*)calloc(k_codewords * number_samples, sizeof(double));
  double* XtX       = (double*)calloc(number_samples * number_samples, sizeof(double));
  double* Lchol     = (double*)malloc(sizeof(double) * (sample_block * k_codewords));
  double* alpha_0   = (double*)malloc(sizeof(double) * k_codewords);
  double* Gsub      = (double*)malloc(sizeof(double) * (k_codewords * k_codewords));
  double* tempvec1  = (double*)malloc(sizeof(double) * k_codewords);
  double* tempvec2  = (double*)malloc(sizeof(double) * k_codewords);
  double* c         = (double*)malloc(sizeof(double) * sample_block);

  int* selected_atoms = (int*)malloc(sizeof(int) * k_codewords);
  int* ind            = (int*) malloc(sizeof(int) * sample_block);

  // calculate G, DtX and XtX
  matT_mat(1, dictionary->data, dictionary->data, G, sample_block, k_codewords, k_codewords);
  matT_mat(1, samples->data, samples->data, XtX, sample_block, number_samples, number_samples);
  matT_mat(1, dictionary->data, samples->data, DtX, sample_block, k_codewords, number_samples);

  double eps2, deltaprev = 0;
  double resnorm = 1;
  int pos;

  for (int signum = 0; signum < number_samples; signum++) {

    if (resnorm>eps2 && iterations > 0) {

      matT_vec(1, dictionary->data, samples->data + sample_block * signum, DtX, sample_block, k_codewords);


      memcpy(alpha_0, DtX + k_codewords * signum, k_codewords * sizeof(double));

      for (int i = 0; i < k_codewords; i++) {
        selected_atoms[i] = 0;
      }
    }

    int iteration = 0;
    while (resnorm>eps2 && iteration < iterations) {
      pos = maxabs(alpha_0, k_codewords);

      if (selected_atoms[pos] || alpha_0[pos] * alpha_0[pos]<1E-14) {
        break;
      }

      ind[iteration]    = pos;
      selected_atoms[pos] = 1;


      memcpy(Gsub + iteration * k_codewords, G + pos* k_codewords, k_codewords * sizeof(double));

      if (iteration == 0) {
        *Lchol = 1;
      }
      else {
        vec_assign(tempvec1, Gsub + iteration * k_codewords, ind, iteration);
        backsubst_L(Lchol, tempvec1, tempvec2, sample_block, iteration);
        for (int j = 0; j < iteration; ++j) {
          Lchol[j * sample_block + iteration] = tempvec2[j];
        }

        double sum = 0;
        for (int j = 0; j < iteration; ++j) {
          sum += pow(Lchol[j * sample_block + iteration], 2);
        }
        if ((1 - sum) <= 1e-14) {
          break;
        }
        Lchol[iteration * sample_block + iteration] = sqrt(1 - sum);
      }

      iteration++;

      vec_assign(tempvec1, DtX + k_codewords * signum, ind, iteration);
      cholsolve(Lchol, tempvec1, c, sample_block, iteration);

      mat_vec(1, Gsub, c, tempvec1, k_codewords, iteration);
      memcpy(alpha_0, DtX + k_codewords * signum, k_codewords * sizeof(double));
      vec_sum(-1, tempvec1, alpha_0, k_codewords);

      for (int j = 0; j < iteration; ++j) {
        sparse_coefficients->data[k_codewords * signum + ind[j]] = c[j];
      }
    }
  }

  free(G);
  free(DtX);
  free(XtX);
  free(Lchol);
  free(alpha_0);
  free(Gsub);
  free(tempvec1);
  free(tempvec2);
  free(c);
  free(selected_atoms);
  free(ind);
}

double calculate_norm(double* array, int begin, int end) {

  double acc = 0;
   for(int i = begin; i < end; i++){
     acc += (array[i] * array[i]);
   }
   return sqrt(acc);
}

void print_dictionary_csv(dictionary_t* dictionary, char* file_name) {

  FILE *data;
  data = fopen(file_name, "w+");

  for (size_t y = 0; y < dictionary->sample_size; y++) {
    for (size_t x = 0; x < dictionary->k; x++) {
      fprintf(data, "%f,", dictionary->data[(x * dictionary->sample_size) + y]);
    }

    fprintf(data, "\n");
  }

  fclose(data);
}

void print_sparse_csv(sparse_coeff_t* sparse_coeff, char* file_name) {

  FILE *data;
  data = fopen(file_name, "w+");

  for (size_t y = 0; y < sparse_coeff->k; y++) {
    for (size_t x = 0; x < sparse_coeff->number_blocks; x++) {
      fprintf(data, "%f,", sparse_coeff->data[(x * sparse_coeff->k) + y]);
    }

    fprintf(data, "\n");
  }

  fclose(data);
}
