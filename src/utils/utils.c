#include "../libit/types.h"
#include "./utils.h"

/** \brief Copy image from the buffer
    \param in       The input buffer.
    \param out		The output image.
    \param w		The image width.
    \param h		The image height.
    \param bpp		The bits per pixel.
*/
void utils_image_copy_n(uint8 *in, int16 *out, uint32 w, uint32 h, uint32 bpp)
{
    uint32 i, size = w*h;
    //uint8 *in1 = &in[1];

    if(bpp > 8){
        for(i=0; i<size; i++) {
            //For Aptina sensor
            out[i] = ((in[(i<<1)+1]<<8) | in[(i<<1)]);
            //For Sony sensor
            //out[i] = ((in[(i<<1)]<<8) | in[(i<<1)+1]);
            //img[i] = ((buff[(i<<1)]<<8) | buff[(i<<1)+1]);
            //printf("MSB = %d LSB = %d img = %d shift = %d\n", buff[(i<<1)], buff[(i<<1)+1], ((buff[(i<<1)]) | buff[(i<<1)+1]<<8), shift);
        }
    } else
        for(i=0; i<size; i++) out[i] = in[i];
}
