#ifndef ITLIB_H_
#define ITLIB_H_

#include "../libit/types.h"
#include "../filters/filters.h"
#include "../utils/utils.h"
#include "../hdr/hdr.h"
#include "../transforms/transforms.h"
#include "../segmentation/segmentation.h"
/**
    \brief The states of image transform
 */
typedef struct {
    void *pic[2];   /**<  Pointer to current image */
    int w;          /**<  Image widht */
    int h;          /**<  Image height */
    int bpp;        /**<  Bytes per pixel */
    int colort;     /**<  Color types */
    int bg;         /**<  Bayer grids pattern */
} TransState;

#endif /*ITLIB_H_*/
