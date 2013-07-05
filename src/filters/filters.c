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

void inline static clear_line_32(int *in, uint32 w, uint32 sh)
{
    uint32 i, w1 = w + (sh<<1);
    for(i=0; i < w1; i++) in[i] = 0;
}

void inline static cp_line_32(int16 *in, int *l, uint32 w, uint32 sh)
{
    uint32 i;
    for(i=0; i < sh; i++) l[i] = in[sh-i];
    for(i=0; i < w ; i++) l[i+sh] = in[i];
    for(i=0; i < sh; i++) l[i+sh+w] = in[w-sh-i];
}

void inline static cp_line_16(int16 *in, int16 *l, uint32 w, uint32 sh)
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
    cp_line_16(&in[w], l[0], w, sh);
    cp_line_16(&in[0], l[1], w, sh);

    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + w : yw - w;
        cp_line_16(&in[yw1], l[2], w, sh);

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
    cp_line_16(&in[w*2], l[0], w, sh);
    cp_line_16(&in[w  ], l[1], w, sh);
    cp_line_16(&in[0  ], l[2], w, sh);
    cp_line_16(&in[w  ], l[3], w, sh);

    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + ws : yw;
        cp_line_16(&in[yw1], l[4], w, sh);

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
            //printf(" yx = %d min = %d med = %d max = %d ing[x] = %d img1[x] = %d\n",
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
    int i, x, xi, y, yi, yw, yx;
    int hg = sg*sg;
    int ex[256], smi, smi1, blm, cf;

    int sh, ns = (br<<1) + 1; //, bs = (br+1)*((br>>1)+1);
    int w1 = w + (br<<1);
    int16 *l[ns], *m[ns], *tm;
    int64_t b;

    //Finding b coefficient to change / to *
    //for(i=1; bs>>i; i++);
    //sh = 63 - i - 16;
    //b = (1LL<<sh)/bs;

    //Make lut table to remove exp
    for(i=0; i < 256; i++){
        ex[i] = (int)(exp(-(double)i*i/(double)hg)*512);
        printf("%3d exp = %d\n", i, ex[i]);
    }

    //Rows buffer for input image
    l[0] = buff;
    for(i=1; i < ns; i++) l[i] = &l[i-1][w1];

    //Rows buffer for averasge image
    m[0] = &l[ns-1][w1];
    for(i=1; i < ns; i++) m[i] = &m[i-1][w1];

    //Prepare first raws
    for(i=0; i < ns - 1; i++){
        if(i < br) cp_line_16(&in[w*(br-i)], l[i], w, br);
        else cp_line_16(&in[w*(i-br)], l[i], w, br);
    }

    for(i=0; i < ns - 1; i++){
        if(i < br) cp_line_16(&avr[w*(br-i)], m[i], w, br);
        else cp_line_16(&avr[w*(i-br)], m[i], w, br);
    }

    for(y = 0; y < h; y++){
        yw = y*w;
        if(y+br > h-1) {
            cp_line_16(&in[w*(((h-1)<<1)-y)], l[ns-1], w, br);
            cp_line_16(&avr[w*(((h-1)<<1)-y)], m[ns-1], w, br);
        } else {
            cp_line_16(&in[w*(y+br)], l[ns-1], w, br);
            cp_line_16(&avr[w*(y+br)], m[ns-1], w, br);
        }

        for(x = 0; x < w; x++){
            yx = yw + x;
            smi1 = smi = 0;

            for(yi=0; yi < ns; yi+=2){
                for(xi=x; xi <= x+ns; xi+=2){

                    blm  = abs(m[yi][xi] - m[br][xi]);
                    //printf("blm = %d\n", blm);

                    cf = blm > 255 ? 0 : ex[blm];
                    smi1 += cf;
                    smi += l[yi][xi]*cf;
                }
            }
            //tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = tm;
            out[yx] = smi/smi1;

            //printf("Start x = %d y = %d avr = %d in = %d out = %d blm = %d\n", xb, yb, avr[yxb], ing[yxb], out[yxb], blm);
        }
        tm = l[0];
        for(i=1; i < ns; i++) l[i-1] = l[i];
        l[ns-1] =  tm;
        tm = m[0];
        for(i=1; i < ns; i++) m[i-1] = m[i];
        m[ns-1] =  tm;
    }
}


/** \brief Calculate the determinant of Hessian of grey image.
    \param in 	The input 16 bits image.
    \param out	The output 16 bits image.
    \param buff	The temporary buffer.
    \param w	The image width.
    \param h	The imahe height.
*/
void filters_hessian(int16 *in, int16 *out, uint32 *buff, const int w, const int h)
{
    //  Dxx                          Dyy                          Dxy
    //  0  0  0  0  0  0  0  0  0    0  0  1  1  1  1  1  0  0    0  0  0  0  0  0  0  0  0
    //  0  0  0  0  0  0  0  0  0    0  0  1  1  1  1  1  0  0    0  1  1  1  0 -1 -1 -1  0
    //  1  1  1 -2 -2 -2  1  1  1    0  0  1  1  1  1  1  0  0    0  1  1  1  0 -1 -1 -1  0
    //  1  1  1 -2 -2 -2  1  1  1    0  0 -2 -2 -2 -2 -2  0  0    0  1  1  1  0 -1 -1 -1  0
    //  1  1  1 -2 -2 -2  1  1  1    0  0 -2 -2 -2 -2 -2  0  0    0  0  0  0  0  0  0  0  0
    //  1  1  1 -2 -2 -2  1  1  1    0  0 -2 -2 -2 -2 -2  0  0    0 -1 -1 -1  0  1  1  1  0
    //  1  1  1 -2 -2 -2  1  1  1    0  0  1  1  1  1  1  0  0    0 -1 -1 -1  0  1  1  1  0
    //  0  0  0  0  0  0  0  0  0    0  0  1  1  1  1  1  0  0    0 -1 -1 -1  0  1  1  1  0
    //  0  0  0  0  0  0  0  0  0    0  0  1  1  1  1  1  0  0    0  0  0  0  0  0  0  0  0
    //
    //  det = Dxx*Dyy - Dxy*Dxy

    int x, y, x1, y1, yw, yw1, yx, yx1;
    int dxx1, dxx2, dxx3, dyy1, dyy2, dyy3, dxy1, dxy2, dxy3, dxy4;
    int dxx, dyy, dxy;
    int br = 5, h1 = h + (br<<1), w1 = w + (br<<1);
    int w2 = w1*2, w3 = w1*3, w4 = w1*4, w5 = w1*5;

    uint32 *ing;

    ing = (uint32*)malloc(w1*h1*sizeof(uint32));
    if (ing == NULL) {
        fprintf(stderr, "Error! utils_average: Can't allocate memory\n");
        return;
    }

    utils_integral(in, ing, buff, w, h, br);

    for(y=br, y1 = 0; y1 < h; y++, y1++){
        yw = y*w1;
        yw1 = y1*w;
        for(x=br, x1 = 0; x1 < w; x++, x1++){
            yx = yw + x;
            yx1 = yw1 + x1;

            dxx1 = ing[yx+w2-2] + ing[yx-w3-5] - ing[yx-w3-2] - ing[yx+w2-5];
            dxx2 = ing[yx+w2+1] + ing[yx-w3-2] - ing[yx-w3+1] - ing[yx+w2-2];
            dxx3 = ing[yx+w2+4] + ing[yx-w3+1] - ing[yx-w3+4] - ing[yx+w2+1];
            dxx = dxx1 + dxx3 - dxx2*2;
            //printf("dxx1 = %d dxx2 = %d dxx3 =%d dxx = %d\n", dxx1, dxx2, dxx3, dxx);

            dyy1 = ing[yx-w2+2] + ing[yx-w5-3] - ing[yx-w5+2] - ing[yx-w2-3];
            dyy2 = ing[yx+w +2] + ing[yx-w2-3] - ing[yx-w2+2] - ing[yx+w -3];
            dyy3 = ing[yx+w4+2] + ing[yx+w1-3] - ing[yx+w1+2] - ing[yx+w4-3];
            dyy = dyy1 + dyy3 - dyy2*2;
            //printf("dyy1 = %d dyy2 = %d dyy3 =%d dyy = %d\n", dyy1, dyy2, dyy3, dyy);

            dxy1 = ing[yx-w1-1] + ing[yx-w4-4] - ing[yx-w4-1] - ing[yx-w1-4];
            dxy2 = ing[yx-w1+3] + ing[yx-w4  ] - ing[yx-w4+3] - ing[yx-w1  ];
            dxy3 = ing[yx+w3-1] + ing[yx   -4] - ing[yx   -1] - ing[yx+w3-4];
            dxy4 = ing[yx+w3+3] + ing[yx     ] - ing[yx   +3] - ing[yx+w3  ];
            dxy = dxy1 - dxy2 - dxy3 + dxy4;
            //printf("dxy1 = %d dxy2 = %d dxy3 =%d dxy4 =%d dxy = %d\n", dxy1, dxy2, dxy3, dxy4, dxy);
            //hs[yx] = abs(dxx*dyy - dxy*dxy);
            out[yx1] = abs(dxy);

            //det = dxx*dyy - dxy*dxy;
            //img[yx] = abs(det) > 255 ? 255 : abs(det);
            //img[yx] = abs(det)>> 14;
            //printf("yx = %d img = %d\n", yx, det);
        }
    }
}

/** \brief The mean square error (MSE) regression of plane denoise filter
    \param in	The input 16 bits bayer image.
    \param out	The output 16 bits bayer image.
    \param buff	The temporary buffer.
    \param br   The radius around the pixel.
    \param w    The image width.
    \param h 	The image height.
*/
void filters_denoise_regression_bayer(int16 *in, int16 *out, int *buff, const int br, const int w, const int h)
{
    int i, x, xi, y, yi, yw, yx;
    //int hg = sg*sg;
    //int ex[256], smi, smi1, blm, cf;

    int sh, ns = (br<<1) + 1, bs = (br+1)*(br+1);
    int w1 = w + (br<<1);
    int xa, ya, xa2, ya2, avr, av, bw = w*br;
    int *l[ns], *m[ns], *tm;
    int a, b, c;
    int64_t d;


    for(i=1; bs>>i; i++);
    sh = 63 - i - 16;
    d = (1LL<<sh)/bs;

    //Precalculate sum(xi**2)
    xa2 = 0; ya2 = 0;
    for(y=-br; y <= br ; y+=2){
        for(x=-br; x <= br ; x+=2){
            printf("x = %d y = %d\n", x, y);
            xa2 += x*x;
            ya2 += y*y;
        }
    }
    printf(" xa2 = %d ya2 = %d bs = %d\n", xa2, ya2, bs);

    //Rows buffer for input image
    l[0] = buff;
    for(i=1; i < ns; i++) l[i] = &l[i-1][w1];

    //Rows buffer for averasge image
    m[0] = &l[ns-1][w1];
    for(i=1; i < ns; i++) m[i] = &m[i-1][w1];

    //Prepare first raws
    for(i=0; i < ns - 1; i++){
        if(i < br) cp_line_32(&in[w*(br-i)], l[i], w, br);
        else cp_line_32(&in[w*(i-br)], l[i], w, br);
    }

    for(i=0; i < ns - 1; i++){
        clear_line_32(m[i], w, br);
    }

    for(y = 0; y < h; y++){
        yw = y*w;
        if(y+br > h-1) {
            cp_line_32(&in[w*(((h-1)<<1)-y)], l[ns-1], w, br);
            clear_line_32(m[ns-1], w, br);
        } else {
            cp_line_32(&in[w*(y+br)], l[ns-1], w, br);
            clear_line_32(m[ns-1], w, br);
        }

        for(x = 0; x < w; x++){
            yx = yw + x;
            avr = 0;
            xa = 0; ya = 0;

            //Calculate plane z = ax + by + c coefficients
            for(yi=-br; yi <= br; yi+=2){
                for(xi=-br; xi <= br; xi+=2){
                    av = l[yi + br][x + xi + br];
                    avr += av;
                    xa += av*xi;
                    ya += av*yi;
                    //printf("blm = %d\n", blm);
                }
            }
            a = (xa<<10)/xa2; b = (ya<<10)/ya2; c = (avr<<10)/bs;

            for(yi=-br; yi <= br; yi+=2){
                for(xi=-br; xi <= br; xi+=2){
                    m[yi + br][x + xi + br] += (a*xi + b*yi + c)>>10;
                }
            }

            if(y >= br) out[yx - bw] = m[0][x+br]/bs;
            //yx - bwif(y >= br) printf("a = %d b = %d c = %d in = %d out = %d pl = %d\n", a>>8, b>>8, c>>8, in[yx - bw], out[yx - bw], (c)>>8);

            //printf("Start x = %d y = %d avr = %d in = %d out = %d blm = %d\n", xb, yb, avr[yxb], ing[yxb], out[yxb], blm);
        }
        tm = l[0];
        for(i=1; i < ns; i++) l[i-1] = l[i];
        l[ns-1] =  tm;
        tm = m[0];
        for(i=1; i < ns; i++) m[i-1] = m[i];
        m[ns-1] =  tm;
    }
}


/** \brief The simple MSE calculation in each pixel with radius br.
    \param in	The input 16 bits bayer image.
    \param avr	The input 16 bits average image.
    \param out	The output 16 bits bayer image.
    \param buff	The temporary buffer.
    \param br   The radius around the pixel.
    \param bpp  The bits per pixel.
    \param w    The image width.
    \param h 	The image height.
*/
void filters_MSE_bayer(int16 *in, int16 *avr, int16 *out, int16 *buff, const int br, const int bpp, const int w, const int h)
{
    int i, x, xi, y, yi, yw, yx;
    int  blm, tmp;

    int sh, ns = (br<<1) + 1, bs = (br+1)*(br+1);
    int w1 = w + (br<<1), s[4], v[4], per = 0;
    int min, min1, min2;
    int16 *l[ns], *m[ns], *tm;
    int64_t b;

    //Finding b coefficient to change / to *
    for(i=1; bs>>i; i++);
    sh = 63 - i - 16;
    b = (1LL<<sh)/bs;
    //b = (1LL<<sh);

    //Rows buffer for input image
    l[0] = buff;
    for(i=1; i < ns; i++) l[i] = &l[i-1][w1];

    //Rows buffer for average image
    m[0] = &l[ns-1][w1];
    for(i=1; i < ns; i++) m[i] = &m[i-1][w1];

    //Prepare first raws
    for(i=0; i < ns - 1; i++){
        if(i < br) cp_line_16(&in[w*(br-i)], l[i], w, br);
        else cp_line_16(&in[w*(i-br)], l[i], w, br);
    }

    for(i=0; i < ns - 1; i++){
        if(i < br) cp_line_16(&avr[w*(br-i)], m[i], w, br);
        else cp_line_16(&avr[w*(i-br)], m[i], w, br);
    }

    for(y = 0; y < h; y++){
        yw = y*w;
        if(y+br > h-1) {
            cp_line_16(&in[w*(((h-1)<<1)-(y+br))], l[ns-1], w, br);
            cp_line_16(&avr[w*(((h-1)<<1)-(y+br))], m[ns-1], w, br);
        } else {
            cp_line_16(&in[w*(y+br)], l[ns-1], w, br);
            cp_line_16(&avr[w*(y+br)], m[ns-1], w, br);
        }

        for(x = 0; x < w; x++){
            yx = yw + x;
            blm = 0;
            s[0] = 0; s[1] = 0; s[2] = 0; s[3] = 0;
            v[0] = 0; v[1] = 0; v[2] = 0; v[3] = 0;

            //if(l[br][x+br] > 100) printf("m = %d l = %d\n", m[br][x+br], l[br][x+br]);
            for(yi=0; yi < ns; yi+=2){
                for(xi=x; xi <= x+ns; xi+=2){
                    if(yi == br) {
                        s[0] += abs(l[yi][xi] - l[br][x+br]);
                        v[0] += l[yi][xi];
                    } else if (xi-x == yi){
                        s[1] += abs(l[yi][xi] - l[br][x+br]);
                        v[1] += l[yi][xi];
                    } else if (xi == x+br){
                        s[2] += abs(l[yi][xi] - l[br][x+br]);
                        v[2] += l[yi][xi];
                    } else if (ns-xi+x-1 == yi){
                        s[3] += abs(l[yi][xi] - l[br][x+br]);
                        v[3] += l[yi][xi];
                    }

                    blm  += abs(l[yi][xi] - m[br][x+br]);
                    //blm  += l[yi][xi];
                    //if(l[br][x+br] > 100) printf("%6d ", l[yi][xi]);
                }
                //if(l[br][x+br] > 100) printf("\n");
            }
            v[0] = v[0] - l[br][x+br];

            /*
            if(l[br][x+br]) {
                tmp = ((s[0]+s[1]+s[2]+s[3])>>4)*100/l[br][x+br];
                //printf("tmp = %d\n", tmp);
                if(tmp > per){

                    min1 = (s[0] > s[1]) ? 1 : 0;
                    min2 = (s[2] > s[3]) ? 3 : 2;
                    min = (s[min2] > s[min1]) ? min1 : min2;
                    out[yx] = (v[min] + l[br][x+br])/(br+1);
                } else {
                    out[yx] = blm/bs;
                }
            } else {
                out[yx] = blm/bs;
            }
            */
            out[yx] = blm/bs;

            //printf("mse = %d val = %d per = %d \n", blm/bs, l[br][x+br], blm*100/(bs*l[br][x+br]));

            //if(l[br][x+br] > 100) printf("blm = %d per = %d out = %d max = %d s = %d v = %d bs+1 = %d\n",
            //                             blm*b>>sh, (blm*b>>sh)*100/l[br][x+br], out[yx], min, s[min], v[min], br+1);
            //if(l[br][x+br] > 100) printf("s0 = %d s1 = %d s2 = %d s3 = %d\n", s[0], s[1], s[2], s[3]);
            //if(l[br][x+br] > 100) printf("s0 = %d s1 = %d s2 = %d s3 = %d\n", v[0], v[1], v[2], v[3]);

            //tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = tm;
            //tmp = blm*b>>sh;
            //out[yx] = tmp > max ? max : tmp;
            //printf("mean = %d\n", out[yx]);
            //out[yx] = blm*b>>sh;

            //printf("Start x = %d y = %d avr = %d in = %d out = %d blm = %d\n", xb, yb, avr[yxb], ing[yxb], out[yxb], blm);
        }
        tm = l[0];
        for(i=1; i < ns; i++) l[i-1] = l[i];
        l[ns-1] =  tm;
        tm = m[0];
        for(i=1; i < ns; i++) m[i-1] = m[i];
        m[ns-1] =  tm;
    }
}



/** \brief Culculate basis fuction for bicubic B-spline.
    \param n3	The output basis vector.
    \param t	The t parametr from 1 to 2.
*/
static void get_basis(int *n, const float t)
{
    float  t1, t2, t3, t4, t5, t6;

    t1 = 1.-t;
    t2 = 2.-t;
    t3 = 3.-t;
    t4 = 4.-t;
    t5 = 5.-t;
    t6 = 6.-t;

    n[0] = (t4*t4*t4)*8.;
    n[1] = (-t1*t4*t4 - t2*t4*t5 - t3*t5*t5)*8.;
    n[2] = (t2*t2*t4 + t2*t3*t5 + t3*t3*t6)*8.;
    n[3] = (-t3*t3*t3)*8.;

    printf("get_basis: n0 = %d n1 = %d n2 = %d n3 = %d\n", n[0], n[1], n[2], n[3]);
}

/** \brief Interpolation points in the central square for red and blue color
    \param p0	The 1 row of pixels.
    \param p1	The 2 row of pixels.
    \param p2	The 3 row of pixels.
    \param p3	The 4 row of pixels.
    \param nx   Basis vector for x direction.
    \param ny   Basis vector for y direction.
    \retval     The interpolation point.
*/
static inline int16 b_spline_rb(const int16 *p0, const int16 *p1, const int16 *p2, const int16 *p3, const int *nx, const int *ny)
{
    // p0[0]----p0[2]----p0[4]----p0[6]
    //  |        |        |        |
    //  |        |        |        |
    // p1[0]----p1[2]----p1[4]----p1[6]
    //  |        |  x     |        |
    //  |        |     x  |        |
    // p2[0]----p2[2]----p2[4]----p2[6]
    //  |        |        |        |
    //  |        |        |        |
    // p3[0]----p3[2]----p3[4]----p3[6]

    int tm[4];
    int64_t tmp;
    int b = (1<<30)/36;


    tm[0] = p0[0]*nx[0] + p0[2]*nx[1] + p0[4]*nx[2] + p0[6]*nx[3];
    tm[1] = p1[0]*nx[0] + p1[2]*nx[1] + p1[4]*nx[2] + p1[6]*nx[3];
    tm[2] = p2[0]*nx[0] + p2[2]*nx[1] + p2[4]*nx[2] + p2[6]*nx[3];
    tm[3] = p3[0]*nx[0] + p3[2]*nx[1] + p3[4]*nx[2] + p3[6]*nx[3];

    tmp = tm[0]*ny[0] + tm[1]*ny[1] + tm[2]*ny[2] + tm[3]*ny[3];

    return (tmp>>6)/36;
    //return (int)(tmp*b>>36);
}

/** \brief Interpolation points in the central square for green color
    \param p0	The 1 row of pixels.
    \param p1	The 2 row of pixels.
    \param p2	The 3 row of pixels.
    \param p3	The 4 row of pixels.
    \param p4	The 4 row of pixels.
    \param p5	The 5 row of pixels.
    \param p6	The 6 row of pixels.
    \param nx   Basis vector for x direction.
    \param ny   Basis vector for y direction.
    \retval     The interpolation point.
*/
static inline int16 b_spline_g(const int16 *p0, const int16 *p1, const int16 *p2, const int16 *p3,
                                const int16 *p4, const int16 *p5, const int16 *p6, const int *nx, const int *ny)
{
    //                p0[3]
    //              /      \
    //           p1[2]----p1[4]
    //          /     \ /     \
    //      p2[1]----p2[3]----p2[5]
    //     /     \  /     \  /    \
    // p3[0]----p3[2]----p3[4]----p3[6]
    //     \     /  \   /    \    /
    //      p4[1]----p4[3]----p4[5]
    //          \    /    \    /
    //           p5[2]----p5[4]
    //              \      /
    //                p6[3]

    int tm[4];
    int64_t tmp;
    int b = (1<<30)/36;


    tm[0] = p0[3]*nx[0] + p1[4]*nx[1] + p2[5]*nx[2] + p3[6]*nx[3];
    tm[1] = p1[2]*nx[0] + p2[3]*nx[1] + p3[4]*nx[2] + p4[5]*nx[3];
    tm[2] = p2[1]*nx[0] + p3[2]*nx[1] + p4[3]*nx[2] + p5[4]*nx[3];
    tm[3] = p3[0]*nx[0] + p4[1]*nx[1] + p5[2]*nx[2] + p6[3]*nx[3];

    tmp = tm[0]*ny[0] + tm[1]*ny[1] + tm[2]*ny[2] + tm[3]*ny[3];

    return (tmp>>6)/36;
    //return tmp*b>>36;
}

/** \brief Bicubic B-spline interpolation of 4x4 pixels patch.
    \param in       The input 16 bits bayer image.
    \param out      The output 16 bits RGB interpolated image.
    \param buff     The temporary buffer.
    \param w        The image width.
    \param h        The image height.
    \param bay		The Bayer grids pattern.
*/
void  b_spline_interpolation(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bay)
{
    int i,  x, y, yw, yx, yx3, st = 2, w3 = w*3;

    int br = 2, br1 = br-1, br2 = br<<1, br3 = br+1, ns = 8, w1 = w + (br2<<1);
    int n1[4], n2[4], n3[4];
    int16 *l[ns], *tm;
    /*
    switch(bay){
        case(BGGR):{ xs = 1; ys = 1; w1 = w+1; h1 = h+1; break; }
        case(GRBG):{ xs = 1; ys = 0; w1 = w+1; h1 = h;   break; }
        case(GBRG):{ xs = 0; ys = 1; w1 = w;   h1 = h+1; break; }
        case(RGGB):{ xs = 0; ys = 0; w1 = w;   h1 = h;   break; }
    }
*/
    //Prepare two basis vectors
    get_basis(n1, 3.);
    get_basis(n2, 3.5);
    get_basis(n3, 4.);

    //Rows buffer for input image
    l[0] = buff;
    for(i=1; i < ns; i++) l[i] = &l[i-1][w1];

    //Prepare first raws
    for(i=0; i < ns - 1; i++){
        if(i < br3) cp_line_16(&in[w*(br3-i)], l[i], w, br2);
        else cp_line_16(&in[w*(i-br3)], l[i], w, br2);
    }

    //for(y = 0; y < 3; y++){
    for(y = 0; y < h; y++){
        yw = y*w;
        if(y+br2 > h-1) {
            cp_line_16(&in[w*(((h-1)<<1)-(y+br2))], l[ns-1], w, br2);
        } else {
            cp_line_16(&in[w*(y+br2)], l[ns-1], w, br2);
        }
        //for(x = 0; x < 3; x++){
        for(x = 0; x < w; x++){
            yx = yw + x;
            yx3 = yx*3;
            if(!(x&1) && !(y&1)) {
                //Red
                out[yx3]        = b_spline_rb(&l[1][x+br], &l[3][x+br], &l[5][x+br], &l[7][x+br], n1, n1);
                out[yx3+3]      = b_spline_rb(&l[1][x+br], &l[3][x+br], &l[5][x+br], &l[7][x+br], n2, n1);
                out[yx3+w3]     = b_spline_rb(&l[1][x+br], &l[3][x+br], &l[5][x+br], &l[7][x+br], n1, n2);
                out[yx3+w3+3]   = b_spline_rb(&l[1][x+br], &l[3][x+br], &l[5][x+br], &l[7][x+br], n2, n2);
                //out[yx3] = 0; out[yx3+3] = 0; out[yx3+w3] = 0; out[yx3+w3+3]  = 0;
                //Green
                //out[yx3+1] = b_spline_g(&l[0][x+br1], &l[1][x+br1], &l[2][x+br1], &l[3][x+br1], &l[4][x+br1], &l[5][x+br1], &l[6][x+br1], n2, n2);
                //out[yx3+4] = b_spline_g(&l[0][x+br1], &l[1][x+br1], &l[2][x+br1], &l[3][x+br1], &l[4][x+br1], &l[5][x+br1], &l[6][x+br1], n3, n1);
                out[yx3+1] = 0; out[yx3+4] = 0;
            } else if ((x&1) && (y&1)){
                //Blue
                //out[yx3+2]      = b_spline_rb(&l[1][x+br], &l[3][x+br], &l[5][x+br], &l[7][x+br], n1, n1);
                //out[yx3+5]      = b_spline_rb(&l[1][x+br], &l[3][x+br], &l[5][x+br], &l[7][x+br], n2, n1);
                //out[yx3+w3+2]   = b_spline_rb(&l[1][x+br], &l[3][x+br], &l[5][x+br], &l[7][x+br], n1, n2);
                //out[yx3+w3+5]   = b_spline_rb(&l[1][x+br], &l[3][x+br], &l[5][x+br], &l[7][x+br], n2, n2);
                out[yx3+2] = 0; out[yx3+5] = 0; out[yx3+w3+2] = 0; out[yx3+w3+5] = 0;
                //Green
                //out[yx3+1] = b_spline_g(&l[0][x+br1], &l[1][x+br1], &l[2][x+br1], &l[3][x+br1], &l[4][x+br1], &l[5][x+br1], &l[6][x+br1], n2, n2);
                //out[yx3+4] = b_spline_g(&l[0][x+br1], &l[1][x+br1], &l[2][x+br1], &l[3][x+br1], &l[4][x+br1], &l[5][x+br1], &l[6][x+br1], n3, n1);
                out[yx3+1] = 0; out[yx3+4] = 0;
            }
            /*
            //out[yx] = in[yx]; out[yx+1] = in[yx]; out[yx+w] =in[yx]; out[yx+w+1] = in[yx];
            out[yx] = b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n1, n1);
            out[yx+1] = b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n2, n1);
            out[yx+w] = b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n1, n2);
            out[yx+w+1] = b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n2, n2);


            printf("x = %d y = %d\n", x, y);
            printf(" %4d %4d\n %4d %4d\n\n", in[yx], in[yx+st], in[yx+w*st], in[yx+st+w*st]);
            printf(" %4d %4d %4d\n %4d %4d %4d\n %4d %4d %4d\n\n",
                   b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n1, n1),
                    b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n2, n1),
                    b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n3, n1),
                    b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n1, n2),
                    b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n2, n2),
                    b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n3, n2),
                    b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n1, n3),
                    b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n2, n3),
                    b_spline_4x4(&l[0][x + br], &l[2][x + br], &l[4][x + br], &l[6][x + br], n3, n3));
            */
        }
        tm = l[0];
        for(i=1; i < ns; i++) l[i-1] = l[i];
        l[ns-1] =  tm;
    }

}
