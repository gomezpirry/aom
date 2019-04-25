#include "image_handler.h"

void split_in_samples(sample_t* samples, image_t* image, bool rowMajor) {
  int block_it = 0;
  int block_size = (int)(sqrt(samples->sample_size));

  if (rowMajor) {
    for (uint32_t y = 0; y < image->height; y += block_size)
      for (uint32_t x = 0; x < image->width; x += block_size) {
        for (int row = 0; row < block_size; row++)
          for (int column = 0; column < block_size; column++) {
            samples->data[(block_it + (row * block_size * samples->number_blocks) + (column * samples->number_blocks))] = image->image_data[(x + column) + ((y + row) * image->width)] / 255.0f;
          }
        block_it++;
      }
  } else {
    for (uint32_t y = 0; y < image->height; y += block_size)
      for (uint32_t x = 0; x < image->width; x += block_size) {
        for (int row = 0; row < block_size; row++)
          for (int column = 0; column < block_size; column++) {
            samples->data[(block_it * samples->sample_size) + ((row * block_size) + column)] = image->image_data[(x + column) + ((y + row) * image->width)] / 255.0f;
          }
        block_it++;
      }
  }
}

void extract_block(aom_image_t* image, block_t* block, int x, int y, bool rowMajor) {

  if(rowMajor){
    for (int row = 0; row < block->height; row++)
      for (int column = 0; column < block->width; column++) {
        block->data[(row * block->width)  + (column)] = image->img_data[(x + column) + ((y + row) * image->d_w)];
      }
  } else {
    for (int row = 0; row < block->height; row++)
      for (int column = 0; column < block->width; column++) {
        block->data[(column * block->height) + row] = image->img_data[(x + column) + ((y + row) * image->d_w)];
      }
  }
}

void write_image(const char* filename, image_t* image) {

  IplImage* img = cvCreateImage(cvSize(image->width, image->height), IPL_DEPTH_8U, 1);
  img->imageData = (char*)image->image_data;
  cvSaveImage(filename, img);
}

void write_blocks(double* blocks, unsigned int sample_block, unsigned int number_blocks, image_t* image, bool rowMajor){
  int block_it = 0;
  int block_size = (int)(sqrt(sample_block));

    if (rowMajor) {
      for (uint32_t y = 0; y < image->height; y += block_size)
        for (uint32_t x = 0; x < image->width; x += block_size) {
          for (int row = 0; row < block_size; row++)
            for (int column = 0; column < block_size; column++) {
              image->image_data[(x + column) + ((y + row) * image->width)]  = (unsigned char)blocks[(block_it + (row * block_size * number_blocks) + (column * number_blocks))] * 255;
            }
          block_it++;
        }
    } else {
      for (uint32_t y = 0; y < image->height; y += block_size)
        for (uint32_t x = 0; x < image->width; x += block_size) {
          for (int row = 0; row < block_size; row++)
            for (int column = 0; column < block_size; column++) {
              image->image_data[(x + column) + ((y + row) * image->width)] = (unsigned char)blocks[(block_it * sample_block) + ((row * block_size) + column)] * 255;
            }
          block_it++;
        }
    }
}


image_t* image_alloc(image_t *img, unsigned int d_w, unsigned int d_h) {

  memset(img, 0, sizeof(image_t));

  img->width = d_w;
  img->height = d_h;

  img->image_data = (unsigned char*)malloc(sizeof(unsigned char) * d_w* d_h);

  return img;
}

block_t* block_alloc(block_t *block, uint8_t b_w, uint8_t b_h) {

  memset(block, 0, sizeof(block_t));

  block->width = b_w;
  block->height = b_h;

  block->data = (unsigned char*)malloc(sizeof(unsigned char) * b_w* b_h);

  return block;
}

void read_image(image_t* image, char* file_name) {
  memset(image, 0, sizeof(image_t));

  IplImage* img  = cvLoadImage(file_name, 0);

  image->width = img->width;
  image->height = img->height;

  image->image_data = (unsigned char*)malloc(sizeof(unsigned char) * img->width * img->height);
  memcpy(image->image_data, img->imageData, img->width * img->height);
}

void copy_image(unsigned char* source_data, unsigned int source_width, image_t* dest){

  for(unsigned int y = 0; y < dest->height; y++)
    for(unsigned int x = 0; x < dest->width; x++){
      dest->image_data[x + (y * dest->width)] = source_data[x + (y * source_width)];
    }
}

void calculate_histogram64(unsigned char* img, unsigned int width, unsigned int height, int* histogram){
  int size = width * height;

  for(int i = 0; i < size; i++) {
    int index = (int)(((float)img[i]/256) * HISTOGRAM_BIN);
    histogram[index]++;
  }
}

void calculate_histogram(unsigned char* img, unsigned int width, unsigned int height, float* histogram){
  int size = width * height;

  for(int i = 0; i < size; i++) {
    histogram[(int)img[i]]++;
  }
}

int histogram_difference(float* histogram1, float* histogram2, float* difference, int size) {

  int sum = 0;
  for(int i = 0; i < size; i++) {
    int value = abs(histogram1[i] - histogram2[i]);
    difference[i] = value;
    sum += value;
  }

  return sum;
}

float normalized_histogram_difference(float* histogram1, float* histogram2, float* difference, int pixels, int size) {

  float sum = 0;
  for(int i = 0; i < size; i++) {
    float value = (float)abs(histogram1[i] - histogram2[i]);
    difference[i] += value;
    sum += value;
  }

  return (sum * size) / pixels;
}

float histogram_frame_difference(float* histogram1, float* histogram2, int size){

  float difference = 0;
  for(int i = 0; i < size; i++) {
    difference += abs(histogram1[i] - histogram2[i]);
  }
  return difference/size;
}

float histogram_mean(float* histogram, float total, int size) {

  float sum = 0;

  for(int i = 0; i < size; i++){
    sum += histogram[i] * i;
  }

  return sum/total;
}

float histogram_deviation(float* histogram, float mean, float total, int size) {

  float sum = 0;

  for(int i = 0; i < size; i++){
    sum += (float)pow(histogram[i] - mean, 2);
  }

  return sqrt(sum/(total - 1));
}
