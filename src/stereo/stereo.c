#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../libit/types.h"
#include "./stereo.h"
#include "../transforms/transforms.h"

/** \brief Calulate disparity
    \param in	The input 16 bits first image.
    \param in1 	The input 16 bits second image.
    \param out	The output 16 bits disparity.
    \param buff	The temporary buffer.
    \param w	The image width.
    \param h	The imahe height.
*/
void stereo_disparity(const int16 *in, const int16 *in1, int16 *out, int16 *buff, const int w, const int h)
{
    int i, size = w*h;
    int16 *Y[2];
    Y[0] = buff;
    Y[1] = &buff[w*h];

    trans_rgb_to_yuv444(in , Y[0], &buff[w*h<<1], &buff[w*h<<1], w, h);
    trans_rgb_to_yuv444(in1, Y[1], &buff[w*h<<1], &buff[w*h<<1], w, h);

    for(i=0; i < size; i++){
        //out[i] = abs(Y[0][i] - Y[1][i]);
        out[i] = Y[0][i];
    }

}
