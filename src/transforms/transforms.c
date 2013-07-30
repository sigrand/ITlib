#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../libit/types.h"
#include "./transforms.h"

//#define lb(x) (((x) < 0) ? 0 : (((x) > 255) ? 255 : (x)))

void inline static cp_line_16(int16 *in, int16 *l, uint32 w, uint32 sh)
{
    uint32 i;
    for(i=0; i < sh; i++) l[i] = in[sh-i];
    for(i=0; i < w ; i++) l[i+sh] = in[i];
    for(i=0; i < sh; i++) l[i+sh+w] = in[w-sh-i];
}

/**	\brief Bilinear algorithm from bayer to rgb interpolation.
    \param in	 	The input 16 bits bayer image.
    \param rgb		The output 16 bits rgb image.
    \param buff		The temporary 3 rows buffer
    \param w		The image width.
    \param h		The image height.
    \param bay		The Bayer grids pattern.
    \retval			Output rgb image..
*/
void trans_bay_to_rgb_bi(const int16 *in, int16 *rgb, int16 *buff, const int w, const int h, const BayerGrid bay)
{
    int x, x1, x2, xs, ys, y = 0, wy, xwy3, w2 = w<<1, yw = 0, h1, w1, h2;
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
                rgb[xwy3] 	= 	l1[x2];
                rgb[xwy3+1] = 	(l0[x2] + l2[x2] + l1[x2-1] + l1[x2+1])>>2;
                rgb[xwy3+2] = 	(l0[x2+1] + l2[x2-1] + l0[x2-1] + l2[x2+1])>>2;
            }else if (!(y&1) && (x&1)){
                rgb[xwy3] = 	(l1[x2-1] + l1[x2+1])>>1;
                rgb[xwy3+1] = 	l1[x2];
                rgb[xwy3+2] =	(l0[x2] + l2[x2])>>1;
            }else if ((y&1) && !(x&1)){
                rgb[xwy3] = 	(l0[x2] + l2[x2])>>1;
                rgb[xwy3+1] = 	l1[x2];
                rgb[xwy3+2] =	(l1[x2-1] + l1[x2+1])>>1;
            }else {
                rgb[xwy3] = 	(l0[x2+1] + l2[x2-1] + l0[x2-1] + l2[x2+1])>>2;
                rgb[xwy3+1] = 	(l0[x2] + l2[x2] + l1[x2-1] + l1[x2+1])>>2;
                rgb[xwy3+2] = 	l1[x2];
            }
        }
        tm = l0; l0 = l1; l1 = l2; l2 = tm;
    }
}

/**	\brief Bilinear algorithm for bayer to Y interpolation.
    \param in	 	The input 16 bits bayer image.
    \param out		The output 16 bits grey Y image.
    \param buff		The temporary 3 rows buffer
    \param w		The image width.
    \param h		The image height.
    \param bay		The Bayer grids pattern.
*/
void trans_bay_to_grey_bi(const int16 *in, int16 *out, int16 *buff, const int w, const int h, const BayerGrid bay)
{
    int x, x1, x2, xs, ys, y = 0, wy, yw = 0, h1, w1, h2;
    int r, g, b;
    int16 *l0, *l1, *l2, *tm;
    l0 = buff; l1 = &buff[w+2]; l2 = &buff[(w+2)<<1];

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

            if(!(y&1) && !(x&1)){
                r = l1[x2];
                g = (l0[x2] + l2[x2] + l1[x2-1] + l1[x2+1])>>2;
                b = (l0[x2+1] + l2[x2-1] + l0[x2-1] + l2[x2+1])>>2;
            }else if (!(y&1) && (x&1)){
                r = (l1[x2-1] + l1[x2+1])>>1;
                g = l1[x2];
                b =	(l0[x2] + l2[x2])>>1;
            }else if ((y&1) && !(x&1)){
                r = (l0[x2] + l2[x2])>>1;
                g = l1[x2];
                b =	(l1[x2-1] + l1[x2+1])>>1;
            }else {
                r = (l0[x2+1] + l2[x2-1] + l0[x2-1] + l2[x2+1])>>2;
                g = (l0[x2] + l2[x2] + l1[x2-1] + l1[x2+1])>>2;
                b = l1[x2];
            }
            out[wy] = ((306*(r - g) + 117*(b - g))>>10) + g;

        }
        tm = l0; l0 = l1; l1 = l2; l2 = tm;
    }
}

/** \brief Convert 16 bit RGB image to 16 bits YUV444.
    \param rgb 	The input RGB image.
    \param y	The output Y image.
    \param u	The output U image.
    \param v	The output V image.
    \param w	The image width.
    \param h	The image height.
*/
void trans_rgb_to_yuv444(const int16 *rgb, int16 *Y, int16 *U, int16 *V, const uint32 w, const uint32 h)
{
    int i, i3, sz = w*h;

    for(i=0; i < sz; i++){
        i3 = i*3;
        Y[i] = ((306*(rgb[i3]-rgb[i3 + 1]) + 117*(rgb[i3 + 2]-rgb[i3 + 1]))>>10) + rgb[i3 + 1];
        U[i] = 578*(rgb[i3 + 2] - Y[i])>>10;
        V[i] = 730*(rgb[i3]    -  Y[i])>>10;
    }
}

/** \brief Convert 16 bit RGB image to 16 bits YUV420.
    \param rgb 	The input RGB image.
    \param y	The output Y image.
    \param u	The output U image.
    \param v	The output V image.
    \param w	The image width.
    \param h	The image height.
*/
void trans_rgb_to_yuv420(const int16 *rgb, int16 *Y, int16 *U, int16 *V, const uint32 w, const uint32 h)
{
    int y, x, yw, yx, yx2, yx3, w2 = w>>1;

    for(y=0; y < h; y++){
        yw = y*w;
        for(x=0; x < w; x++){
            yx = yw + x;
            yx3 = yx*3;
            yx2 = (x>>1) + (y>>1)*w2;

            Y[yx] = ((306*(rgb[yx3]-rgb[yx3 + 1]) + 117*(rgb[yx3 + 2]-rgb[yx3 + 1]))>>10) + rgb[yx3 + 1];
            U[yx2] += 578*(rgb[yx3 + 2] - Y[yx])>>10;
            V[yx2] += 730*(rgb[yx3]     - Y[yx])>>10;
            if(x&1 && y&1) {
                U[yx2]>>=2;
                V[yx2]>>=2;
            }
        }
    }
}

/** \brief Culculate basis fuction for bicubic B-spline.
    \param n3	The output basis vector.
    \param t	The t parametr from 3 to 4.
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

/** \brief Bicubic B-spline interpolation of 16 bits bayer image.
    \param in       The input 16 bits bayer image.
    \param out      The output 16 bits RGB interpolated image.
    \param buff     The temporary buffer.
    \param w        The image width.
    \param h        The image height.
    \param bay		The Bayer grids pattern.
*/
void  trans_bay_to_rgb_b_spline(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bay)
{
    int i,  x, x1, xs, xt, y, ys, yt, yw, yx, yx3, dis = 0, size = w*h;

    int br = 2, br1 = br-1, br2 = br<<1, br3 = br+1, ns = 8, w1 = w + (br2<<1);
    int n1[4], n2[4], n3[4];
    int16 *l[ns], *tm;

    switch(bay){
        case(BGGR):{ xs = 1; ys = 1; break; }
        case(GRBG):{ xs = 1; ys = 0; break; }
        case(GBRG):{ xs = 0; ys = 1; break; }
        case(RGGB):{ xs = 0; ys = 0; break; }
    }

    //Prepare two basis vectors
    get_basis(n1, 3.);
    get_basis(n2, 3.5);
    get_basis(n3, 4.);

    //Rows buffer for input image
    l[0] = buff;
    for(i=1; i < ns; i++) l[i] = &l[i-1][w1];

    //Prepare first rows
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
            x1 = x-1;
            xt = x + xs; yt = y + ys;

            if(!(xt&1) && !(yt&1)) {
                //Red
                out[yx3]   = b_spline_rb(&l[1][x+br], &l[3][x+br], &l[5][x+br], &l[7][x+br], n1, n1);
                //Green
                out[yx3+1] = b_spline_g(&l[0][x+br1], &l[1][x+br1], &l[2][x+br1], &l[3][x+br1], &l[4][x+br1], &l[5][x+br1], &l[6][x+br1], n2, n2);
                //Blue
                out[yx3+2] = b_spline_rb(&l[0][x1+br], &l[2][x1+br], &l[4][x1+br], &l[6][x1+br], n2, n2);
                dis += abs(in[yx] - out[yx3]);

            } else if ((xt&1) && !(yt&1)){
                //Red
                out[yx3]   = b_spline_rb(&l[1][x1+br], &l[3][x1+br], &l[5][x1+br], &l[7][x1+br], n2, n1);
                //Green
                out[yx3+1] = b_spline_g(&l[0][x1+br1], &l[1][x1+br1], &l[2][x1+br1], &l[3][x1+br1], &l[4][x1+br1], &l[5][x1+br1], &l[6][x1+br1], n3, n1);
                //Blue
                out[yx3+2] = b_spline_rb(&l[0][x+br], &l[2][x+br], &l[4][x+br], &l[6][x+br], n1, n2);
                dis += abs(in[yx] - out[yx3+1]);

            } else if (!(xt&1) && (yt&1)){
                //Red
                out[yx3]   = b_spline_rb(&l[0][x+br], &l[2][x+br], &l[4][x+br], &l[6][x+br], n1, n2);
                //Green
                out[yx3+1] = b_spline_g(&l[0][x1+br1], &l[1][x1+br1], &l[2][x1+br1], &l[3][x1+br1], &l[4][x1+br1], &l[5][x1+br1], &l[6][x1+br1], n3, n1);
                //Blue
                out[yx3+2] = b_spline_rb(&l[1][x1+br], &l[3][x1+br], &l[5][x1+br], &l[7][x1+br], n2, n1);
                dis += abs(in[yx] - out[yx3+1]);

            } else if ((xt&1) && (yt&1)){
                //Red
                x1 = x-1;
                out[yx3]   = b_spline_rb(&l[0][x1+br], &l[2][x1+br], &l[4][x1+br], &l[6][x1+br], n2, n2);
                //Green
                out[yx3+1] = b_spline_g(&l[0][x+br1], &l[1][x+br1], &l[2][x+br1], &l[3][x+br1], &l[4][x+br1], &l[5][x+br1], &l[6][x+br1], n2, n2);
                //Blue
                out[yx3+2] = b_spline_rb(&l[1][x+br], &l[3][x+br], &l[5][x+br], &l[7][x+br], n1, n1);
                dis += abs(in[yx] - out[yx3+2]);
            }
            //out[yx3] = 0;
            //out[yx3+1] = 0;
            //out[yx3+2] = 0;

            /*
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
    printf("dist = %d\n", dis/size);

    //for(i=0; i< 20; i++) printf("r = %d g = %d b = %d\n", out[i*3], out[i*3+1], out[i*3+2]);
}

