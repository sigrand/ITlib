#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../libit/types.h"
#include "./transforms.h"

//#define lb(x) (((x) < 0) ? 0 : (((x) > 255) ? 255 : (x)))

/**	\brief Bilinear algorithm for bayer to rgb interpolation.
    \param in	 	The input bayer image.
    \param rgb		The output rgb image.
    \param buff		The temporary 3 rows buffer
    \param w		The image width.
    \param h		The image height.
    \param bay		The Bayer grids pattern.
    \param bpp      The bits per pixel.
    \retval			Output rgb image..
*/
uint8* utils_bay16_to_rgb8_bi(const int16 *in, uint8 *rgb, int16 *buff, const int w, const int h, const BayerGrid bay, const int bpp)
{
    int x, x1, x2, xs, ys, y = 0, wy, xwy3, w2 = w<<1, yw = 0, h1, w1, h2, sh = bpp - 8;
    int16 *l0, *l1, *l2, *tm;
    l0 = buff; l1 = &buff[w+2]; l2 = &buff[(w+2)<<1];
    //printf("bpp = %d shift = %d\n", bpp, shift);

    switch(bay){
        case(BGGR):{ xs = 1; ys = 1; w1 = w+1; h1 = h+1; break; }
        case(GRBG):{ xs = 1; ys = 0; w1 = w+1; h1 = h;   break; }
        case(GBRG):{ xs = 0; ys = 1; w1 = w;   h1 = h+1; break; }
        case(RGGB):{ xs = 0; ys = 0; w1 = w;   h1 = h;   break; }
    }
    h2 = h1-1;
    //Create 3 rows buffer for transform
    l0[0] = in[w+1]; for(x=0; x < w; x++) l0[x+1] = in[w+x];  l0[w+1] = l0[w-1];
    l1[0] = in[1];   for(x=0; x < w; x++) l1[x+1] = in[x];    l1[w+1] = l1[w-1];

    for(y=ys, yw=0; y < h1; y++, yw+=w){
        wy = (y == h2) ? yw - w : yw + w;
        l2[0] = in[wy+1]; for(x=0; x < w; x++) l2[x+1] = in[wy + x];  l2[w+1] = l2[w-1];

        for(x=xs, x1=0; x < w1; x++, x1++){
            wy 	= x1 + yw;
            x2 = x1 + 1;
            xwy3 = wy*3;

            if(!(y&1) && !(x&1)){
                rgb[xwy3] 	= 	(l1[x2])>>sh;
                rgb[xwy3+1] = 	(((l0[x2] + l2[x2] + l1[x2-1] + l1[x2+1])>>2))>>sh;
                rgb[xwy3+2] = 	(((l0[x2+1] + l2[x2-1] + l0[x2-1] + l2[x2+1])>>2))>>sh;
            }else if (!(y&1) && (x&1)){
                rgb[xwy3] = 	(((l1[x2-1] + l1[x2+1])>>1))>>sh;
                rgb[xwy3+1] = 	(l1[x2])>>sh;
                rgb[xwy3+2] =	(((l0[x2] + l2[x2])>>1))>>sh;
            }else if ((y&1) && !(x&1)){
                rgb[xwy3] = 	(((l0[x2] + l2[x2])>>1))>>sh;
                rgb[xwy3+1] = 	(l1[x2])>>sh;
                rgb[xwy3+2] =	(((l1[x2-1] + l1[x2+1])>>1))>>sh;
            }else {
                rgb[xwy3] = 	(((l0[x2+1] + l2[x2-1] + l0[x2-1] + l2[x2+1])>>2))>>sh;
                rgb[xwy3+1] = 	(((l0[x2] + l2[x2] + l1[x2-1] + l1[x2+1])>>2))>>sh;
                rgb[xwy3+2] = 	(l1[x2])>>sh;
            }
        }
        tm = l0; l0 = l1; l1 = l2; l2 = tm;
    }
    return rgb;
}


/**	\brief Automatic Color Enhancement algorithm.
    \param in       The input 16 bits rgb image.
    \param out      The output 16 bits rgb image.
    \param buff     The temporary buffer.
    \param bpp      The input image bits per pixel.
    \param bpp1     The output image bits per pixel.
    \param w        The image width.
    \param h        The image height.
*/
void utils_ace(const int16 *in, int16 *out, int *buff, const int w, const int h, const int bpp, const int bpp1)
{
    int x, y, yx, yw, hs = 1<<bpp, size = w*h, sum;
    int *hi;
    int b = (1<<30)/size, shb = 30 - bpp1;

    hi = buff; //hl = &hi[hs]; hr = &hl[hs]; lt = &hr[hs];

    //Fill historgam
    memset(hi, 0, sizeof(int)*hs);
    for(x=0; x < size; x++) hi[in[x]]++;

    //Make lookup table
    sum = 0;
    for(x=0; x < hs; x++) {
        sum += hi[x];
        hi[x] = sum*b>>shb;
    }

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
}

