#ifndef IMAGE_HANDLER_H_
#define IMAGE_HANDLER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "dictionary_learning.h"
#include "../aom/aom_image.h"

#define HISTOGRAM_BIN 256
#define ALPHA_THRES 5.5
#define BLOCK_COMPARE 128

typedef struct image_t {
  unsigned char* image_data;
  uint32_t width;
  uint32_t height;
} image_t;

typedef struct block_t {
  unsigned char* data;
  uint8_t width;
  uint8_t height;
} block_t;

typedef struct motion_t {
  uint16_t x;
  uint16_t y;
} motion_t;

/*!\brief Split an image in blocks
 *
 * \param[in]    samples    Pointer to storage blocks
 * \param[in]    image      Image to split
 * \param[in]    blockSize  Block size to split image
 * \param[in]    rowMajor   Indicates if image is split based on row Major or column Major
 */
void split_in_samples(sample_t* samples, image_t* image, bool rowMajor);

void extract_block(aom_image_t* image, block_t* block, int x, int y, bool rowMajor);

/*!\brief Write an image using OpenCV
 *
 * \param[in]    filename   File path to write image
 * \param[in]    image      Image to write
 */
void write_image(const char* filename, image_t* image);

/*!\brief Open a descriptor, allocating storage for the underlying image
 *
 * Returns a descriptor for storing an image of the given format. The
 * storage for the descriptor is allocated on the heap.
 *
 * \param[in]    img       Pointer to storage for descriptor. If this parameter
 *                         is NULL, the storage for the descriptor will be
 *                         allocated on the heap.
 * \param[in]    d_w       Width of the image
 * \param[in]    d_h       Height of the image
 *
 * \return Returns a pointer to the initialized image descriptor. If the img
 *         parameter is non-null, the value of the img parameter will be
 *         returned.
 */
image_t* image_alloc(image_t *img, unsigned int d_w, unsigned int d_h);

block_t* block_alloc(block_t *img, uint8_t b_w, uint8_t b_h);

void read_image (image_t* image, char* file_name);

void copy_image(unsigned char* source_data, unsigned int source_width, image_t* dest);

void write_blocks(double* blocks, unsigned int sample_block, unsigned int number_blocks, image_t* image, bool rowMajor);

void calculate_histogram64(unsigned char* img, unsigned int width, unsigned int height, int* histogram);

void calculate_histogram(unsigned char* img, unsigned int width, unsigned int height, float* histogram);

int histogram_difference(float* histogram1, float* histogram2, float* difference, int size);

float  normalized_histogram_difference(float* histogram1, float* histogram2, float* difference, int pixels, int size);

float histogram_frame_difference(float* histogram1, float* histogram2, int size);

float histogram_mean(float* histogram, float total, int size);

float histogram_deviation(float* histogram, float mean, float total, int size);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // IMAGE_HANDLER_H_
