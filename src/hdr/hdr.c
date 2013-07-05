#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../libit/types.h"
#include "./hdr.h"
#include "../filters/filters.h"


/**	\brief Automatic Color Enhancement algorithm.
    \param in       The input 16 bits image.
    \param out      The output 16 bits image.
    \param buff     The temporary buffer.
    \param bpp      The input image bits per pixel.
    \param bpp1     The output image bits per pixel.
    \param w        The image width.
    \param h        The image height.
    \param cs      The type of input image ColorSpace.
*/
void hdr_ace(const int16 *in, int16 *out, int *buff, const int w, const int h, const int bpp, const int bpp1, const int cs)
{
    int x, y, yx, yw, hs = 1<<bpp, size;
    int *hi;
    int64_t b, sum, shb = 40 - bpp1;

    if(cs == BAYER || cs == GREY){
        size = w*h;
    } else if(cs == RGB || cs == YUV444){
        size = w*h*3;
    } else if (cs == YUV420){
        size = w*h*3>>1;
    } else {
        fprintf(stderr, "hdr_ace Error! Can't support image format %d\n", cs);
        return;
    }

    b = (1LL<<40)/size;
    //printf("b = %d w*h = %d size = %d\n", b, w*h, size);

    hi = buff; //hl = &hi[hs]; hr = &hl[hs]; lt = &hr[hs];

    //Fill historgam
    memset(hi, 0, sizeof(int)*hs);
    for(x=0; x < size; x++) hi[in[x]]++;

    //Make lookup table
    sum = 0;
    for(x=0; x < hs; x++) {
        sum += hi[x];
        hi[x] = sum*b>>shb;
        //printf("%d %d\n", x, hi[x]);
        //hi[x] = (sum<<bpp1)/size;
    }

    for(x=0; x < size; x++) out[x] = hi[in[x]];
    /*
    for(y=0; y < h; y++){
        yw = y*w;
        for(x=0; x < w; x++){
            yx = yw + x;
            //For sony A55 sensor check
            //if(in[yx] > hs-1) in[yx] = hs-1;
            out[yx] = hi[in[yx]];
            //if(out[yx] < 0 || out[yx] > 255) printf("yx = %d out = %d in = %d lt = %d\n", yx, out[yx], in[yx], lt[in[yx]]);
        }
    }
    */
}

/**	\brief Get the differnce between 16 bits  original and averaged bayer image.
    \param in       The input 16 bits bayer image.
    \param avr      The input 16 bits averaged bayer image.
    \param dif      The output image.
    \param w        The image width.
    \param h        The image height.
    \param bpp      The input images bits per pixel.
*/
void hdr_diff(const int16 *in, const int16 *avr, int16 *dif, const int w, const int h, const int bpp)
{
    int x, y, yw, yx;
    int hs = 2, ws = 2, w2 = w<<1, whs = w*hs, bs;
    int sh = 1<<(bpp-1);

    //bs = (ws+1)*(hs+1); //
    //bs = ((ws<<1)+1)*((hs<<1)+1)>>2;
    //bs = 9;

    //printf("bs = %d\n", bs);

    //utils_integral_bayer(in, buff, w, h);

    for(y=0; y < h; y++){
        yw = y*w;
        for(x=0; x < w; x++){
            yx = yw + x;
            //out[yx] = (buff[yx+ws+whs] + buff[yx-ws-whs-w2-2] - buff[yx+ws-whs-w2] - buff[yx-ws+whs-2])/bs;
            dif[yx] = (in[yx]*15>>4) - (avr[yx]*17>>4) + sh;
            //dif[yx] = lb(dif[yx]);
            //if(dif[yx]< 0) dif[yx] = 0;
            //av[0][yxr] = avr/bs;
            //ws = 2; whs = ws*w; //bs = ((ws<<1)+1)*((ws<<1)+1);
            //av[1][yxr] = ing[yxr+ws+whs] + ing[yxr-ws-whs-w2-2] - ing[yxr+ws-whs-w2] - ing[yxr-ws+whs-2] - av[0][yxr];

            //av[1][yxr] = avr/bs;
            //ws = 6; whs = ws*w; bs = ((ws<<1)+1)*((ws<<1)+1)>>2;
            //avr  = ing[yxr+ws+whs] + ing[yxr-ws-whs-w2-2] - ing[yxr+ws-whs-w2] - ing[yxr-ws+whs-2];
            //av[2][yxr] = avr/bs;
        }
    }
}

/**	\brief Local Automatic Color Enhancement algorithm.
    \param in       The input 16 bits bayer image.
    \param out      The output 16 bits bayer image.
    \param buff     The temporary buffer.
    \param bpp      The input image bits per pixel.
    \param w        The image width.
    \param h        The image height.
*/
void hdr_ace_local(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bpp)
{
    int i, size = w*h;
    int16 *dif = &buff[w*h];

    filters_median_bayer(in, out, buff, w, h, 0);
    hdr_diff(in, out, dif, w, h, bpp);
    hdr_ace(in, out, (int*)buff, w, h, bpp, 8, GREY);
    hdr_ace(dif, in, (int*)buff, w, h, bpp, 6, GREY);

    for(i=0; i < size; i++) out[i] = in[i] + out[i] - 32;
}

