#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libit/types.h"
#include "./filters.h"

static inline void sort_16(int16 *s, int16 x1, int16 x2, int16 x3)
{
        if(x1 > x2){
                if(x2 > x3) { s[0] = x3; s[1] = x2; s[2] = x1; }
                else {
                        if(x1 > x3) { s[0] = x2; s[1] = x3; s[2] = x1; }
                        else 		{ s[0] = x2; s[1] = x1; s[2] = x3; }
                }
        } else {
                if(x1 > x3) { s[0] = x3; s[1] = x1; s[2] = x2; }
                else {
                        if(x3 > x2) { s[0] = x1; s[1] = x2; s[2] = x3; }
                        else 		{ s[0] = x1; s[1] = x3; s[2] = x2; }
                }
        }
}

static inline int16* sort_3_16(int16 *s)
{
    int16 tmp;
    if(s[0] > s[1]) { tmp = s[1]; s[1] = s[0]; s[0] = tmp; }
    if(s[1] > s[2]) { tmp = s[2]; s[2] = s[1]; s[1] = tmp; }
    return s;
}

static inline int16  max_3_16(int16 s0, int16 s1, int16 s2)
{
    return (s0 > s1) ? (s0 > s2 ? s0 : s2) : (s1 > s2 ? s1 : s2);
}

static inline int16  min_3_16(int16 s0, int16 s1, int16 s2)
{
    return (s0 > s1) ? (s1 > s2 ? s2 : s1) : (s0 > s2 ? s2 : s0);
}

static inline int16  median_3_16(int16 s0, int16 s1, int16 s2)
{
    return (s2 > s1) ? (s1 > s0 ? s1 : (s2 > s0 ? s0 : s2))
                     : (s2 > s0 ? s2 : (s1 > s0 ? s0 : s1));
}

void inline static cp_line(int16 *in, int16 *l, uint32 w, uint32 sh)
{
    uint32 i;
    for(i=0; i < sh; i++) l[i] = in[sh-i];
    for(i=0; i < w ; i++) l[i+sh] = in[i];
    for(i=0; i < sh; i++) l[i+sh+w] = in[w-sh-i];
}


/**	\brief	3x3 Median filter for 16 bits grey image.
    \param	in		The input image.
    \param	out     The output filtered image.
    \param	buff    The 3 lines + 18(int16) buffer.
    \param	w       The image width.
    \param  h       The image height.
    \param  type    The type of filter 0 - non adaptive, 1 - adaptive.
*/
void filters_median(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int type)
{
    // s[0]   s[1]  s[2]
    //|-----|-----|-----|
    //|     |     |     | l[0]
    //|-----|-----|-----|
    //|     | yx  |     | l[1]
    //|-----|-----|-----|
    //|     |     |     | l[2]
    //|-----|-----|-----|
    uint32 i, y, x, x2, xs, yx, yw, yw1, h1 = h-1, sh = 1, w2 = w + (sh<<1), sh1 = sh+1;
    int16 *s[3], *l[3], *tm, max, min, med;

    s[0] = buff, s[1] = &s[0][3], s[2] = &s[1][3];
    l[0] = &s[2][3];
    for(i=1; i < 3; i++) l[i] = &l[i-1][w2];

    //Prepare buffer
    cp_line(&in[w], l[0], w, sh);
    cp_line(&in[0], l[1], w, sh);

    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + w : yw - w;
        cp_line(&in[yw1], l[2], w, sh);

        sort_16(s[0], l[0][0], l[1][0], l[2][0]);
        sort_16(s[1], l[0][1], l[1][1], l[2][1]);

        for(x=0; x < w; x++){
            yx = yw + x;
            x2 = x + sh1;
            sort_16(s[2], l[0][x2], l[1][x2], l[2][x2]);
            max = max_3_16(s[0][2], s[1][2], s[2][2]);
            min = min_3_16(s[0][0], s[1][0], s[2][0]);
            med = median_3_16(  max_3_16    (s[0][0], s[1][0], s[2][0]),
                                median_3_16 (s[0][1], s[1][1], s[2][1]),
                                min_3_16    (s[0][2], s[1][2], s[2][2]));
            xs = x+sh;
            if(type){
                if(l[1][xs] == max || l[1][xs] == min)  out[yx] = med;
                else out[yx] = l[1][xs];
            } else {
                out[yx] = med;
            }

            tm = s[0]; s[0] = s[1]; s[1] = s[2]; s[2] = tm;
        }
        tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = tm;
    }
}


/**	\brief	3x3 Median filter for 16 bits bayer image.
    \param	in		The input image.
    \param	out     The output filtered image.
    \param	buff    The 5 lines + 18(int16) buffer.
    \param	w       The image width.
    \param  h       The image height.
    \param  type    The type of filter 0 - non adaptive, 1 - adaptive.
*/
void filters_median_bayer(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int type)
{
    // s00    s10   s01   s11   s02
    //|-----|-----|-----|-----|-----|
    //|     |     |     |     |     | l[0]
    //|-----|-----|-----|-----|-----|
    //|     |     |     |     |     | l[1]
    //|-----|-----|-----|-----|-----|
    //|     |     |  yx |     |     | l[2]
    //|-----|-----|-----|-----|-----|
    //|     |     |     |     |     | l[3]
    //|-----|-----|-----|-----|-----|
    //|     |     |     |     |     | l[4]
    //|-----|-----|-----|-----|-----|

    uint32 i, y, x, x2, xs, yx, yw, yw1, h1 = h-2, sh = 2, w2 = w + (sh<<1), ws = w<<1;
    int16 *s[2][3], *l[5], *tm, max, min, med;

    s[0][0] = buff, s[0][1] = &s[0][0][3], s[0][2] = &s[0][1][3];
    s[1][0] = &s[0][2][3], s[1][1] = &s[1][0][3], s[1][2] = &s[1][1][3];
    l[0] = &s[1][2][3];

    for(i=1; i < 5; i++) l[i] = &l[i-1][w2];

    //Prepare buffer
    cp_line(&in[w*2], l[0], w, sh);
    cp_line(&in[w  ], l[1], w, sh);
    cp_line(&in[0  ], l[2], w, sh);
    cp_line(&in[w  ], l[3], w, sh);

    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + ws : yw;
        cp_line(&in[yw1], l[4], w, sh);

        sort_16(s[0][0], l[0][0], l[2][0], l[4][0]);
        sort_16(s[0][1], l[0][2], l[2][2], l[4][2]);

        sort_16(s[1][0], l[0][1], l[2][1], l[4][1]);
        sort_16(s[1][1], l[0][3], l[2][3], l[4][3]);

        for(x=0; x < w; x++){
            yx = yw + x;
            x2 = x + 4;
            i = x&1 ? 1 : 0;
            //if(x&1) i = 1;
            //else i = 0;

            sort_16(s[i][2], l[0][x2], l[2][x2], l[4][x2]);
            max = max_3_16(s[i][0][2], s[i][1][2], s[i][2][2]);
            min = min_3_16(s[i][0][0], s[i][1][0], s[i][2][0]);
            med = median_3_16(  max_3_16    (s[i][0][0], s[i][1][0], s[i][2][0]),
                                median_3_16 (s[i][0][1], s[i][1][1], s[i][2][1]),
                                min_3_16    (s[i][0][2], s[i][1][2], s[i][2][2]));
            xs = x+sh;

            if(type){
                if(l[2][xs] == max || l[2][xs] == min)  out[yx] = med;
                else out[yx] = l[2][xs];
            } else {
                out[yx] = med;
            }
            //if(max-min > 100) img1[yx] = 0;

            //img1[yx] = max;
            //img2[yx] = max-min;
            //if(img1[yx] == med)
            //printf(" yx = %d min = %d med = %d max = %d in[x] = %d img1[x] = %d\n",
            //       yx, min+shift, med+shift, max+shift, l[2][xs]+shift, img1[yx]+shift);

            tm = s[i][0]; s[i][0] = s[i][1]; s[i][1] = s[i][2]; s[i][2] = tm;
        }
        tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = l[3]; l[3] = l[4]; l[4] = tm;
    }
}

/** \brief The simple Non-Local Means denoise  algorithm.
    \param in	The input 16 bits bayer image.
    \param avr	The input 16 bits average image.
    \param out	The output 16 bits bayer image.
    \param buff	The temporary buffer.
    \param br   The radius around the pixel.
    \param sg   The noise standard deviation.
    \param w    The image width.
    \param h 	The image height.
*/
void filters_NLM_denoise_bayer(int16 *in, int16 *avr, int16 *out, int16 *buff, const int br, const int sg, const int w, const int h)
{
    int i, x, xb, y, yb, yw, yx, yxr, yxb, ybw;
    int hg = sg*sg;
    int ex[256], smi, smi1, blm, cf;

    int sh, ns = (br<<1) + 1, bs = ((br>>1)+1)*((br>>1)+1);
    int w1 = w + (br<<1), h1 = h + (br<<1), h2 = h + br - 1, w2 = w + br - 1, br2 = br<<1;
    int16 *l[ns], *m[ns];
    int64_t b;

    //Finding b coefficient to change / to *
    for(i=1; bs>>i; i++);
    sh = 63 - i - 16;
    b = (1LL<<sh)/bs;

    l[0] = buff;
    for(i=1; i < ns; i++) l[i] = &l[i-1][w1];

    //Make lut table to remove exp
    for(i=0; i < 256; i++){
        ex[i] = (int)(exp(-(double)i*i/(double)hg)*512);
        //printf("%3d exp = %d\n", i, ex[i]);
    }

    //Prepare firbr lines
    for(i=0; i < ns - 1; i++){
        if(i < br) cp_line(&in[w*(br-i)], l[i], w, br);
        else cp_line(&in[w*(i-br)], l[i], w, br);
    }

    for(yb = br; yb <= h2; yb++){
        ybw = yb*w1;
        if(y > h2) cp_line(&in[w*(((h-1)<<1)+br-y)], l[ns-1], w, br);
        else       cp_line(&in[w*(y-br)], l[ns-1], w, br);

        for(xb = br; xb <= w2; xb++){
            yxb = ybw + xb;
            smi1 = smi = 0;

            for(y=0; y < ns; y+=2){
                yw = y*w;
                for(x=xb-br; x <= xb+br; x+=2){
                    yx = yw + x;
                    yxr = yx;

                    blm  = abs(l[y][x] - l[br][xb]);
                    //printf("blm = %d\n", blm);

                    cf = blm > 255 ? 0 : ex[blm];
                    smi1 += cf;
                    smi += m[y][x]*cf;
                }
            }


            for(y=yb-br; y <= yb+br; y+=2){
                yw = y*w;
                for(x=xb-br; x <= xb+br; x+=2){
                    yx = yw + x;
                    yxr = yx;

                    //blm = block_matching_xy(in, w, h, ws, hs, yxr, yxb)*b>>sh;

                    //blm  = (ing[yxr+ws+whs] + ing[yxr-ws-whs-w-1] - ing[yxr+ws-whs-w] - ing[yxr-ws+whs-1] - avr)/bs;
                    //blm  = (av[0][yxr] - av[0][yxb] + av[1][yxr] - av[1][yxb])/bs;

                    blm  = (abs(avr[yxr] - avr[yxb]) + abs(in[yxr] - in[yxb]));
                    //printf("blm = %d\n", blm);

                    cf = blm > 255 ? 0 : ex[blm];
                    smi1 += cf;
                    smi += in[yxr]*cf;

                }
            }
            out[yxb - br -br*w1] = smi/smi1;
            //printf("Start x = %d y = %d avr = %d in = %d out = %d blm = %d\n", xb, yb, avr[yxb], in[yxb], out[yxb], blm);
        }
    }
}


