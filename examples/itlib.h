#ifndef ITLIB_H_
#define ITLIB_H_

#include "../libit/types.h"

/**
    \brief The states of image transform
 */
typedef struct {
    int w;          /**<  Image widht */
    int h;          /**<  Image height */
    int bpp;        /**<  Bytes per pixel */
    int colort;     /**<  Color types */
} TransState;

#endif /*ITLIB_H_*/
