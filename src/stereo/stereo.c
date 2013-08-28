#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../libit/types.h"
#include "./stereo.h"
#include "../transforms/transforms.h"
#include "../segmentation/segmentation.h"


void inline static cp_line_16(const int16 *in, int16 *l, uint32 w, uint32 sh)
{
    uint32 i;
    for(i=0; i < sh; i++) l[i] = in[sh-i];
    for(i=0; i < w ; i++) l[i+sh] = in[i];
    for(i=0; i < sh; i++) l[i+sh+w] = in[w-sh-i];
}

/**	\brief	Check if a pixel the end of edge.
    \param	l0		The pointer to 1 line.
    \param	l1		The pointer to 2 line.
    \param	l2		The pointer to 3 line.
    \param  x		Coordinate pixel in x direction.
    \retval			0 - if not
*/
static inline uint32 block_maching(int16 *l, int16 *r, const int xl, const int xr, const int w)
{
    int x, y, ywl, ywr, yxl, yxr;
    int sad = 0;
    for(y=0; y < 5; y++){
        ywl = y*w + xl;
        ywr = y*w + xr;
        for(x=-2; x < 3; x++){
            yxl = ywl + x;
            yxr = ywr + x;
            sad += abs(l[yxl] - r[yxr]);
            if(xl == xr) printf(" %3d %3d  ", l[yxl], r[yxr]);
        }
        if(xl == xr) printf("\n");
    }
    if(xl == xr) printf("\n");
    return sad;
}

/**	\brief	Check if a pixel.
    \param	l0		The pointer to 1 line.
    \param	l1		The pointer to 2 line.
    \param	l2		The pointer to 3 line.
    \param  x		Coordinate pixel in x direction.
    \retval			0 - if not
*/
static inline uint32 check_pixel(int16 *l1, int16 *l2, int16 *l3, const int x, int *j)
{
    int i = 0;
    *j = 0;

    if(l2[x-1]) { i++; (*j) += 1;   }
    if(l1[x-1]) { i++; (*j) += 2;   }
    if(l1[x  ]) { i++; (*j) += 4;   }
    if(l1[x+1]) { i++; (*j) += 8;   }
    if(l2[x+1]) { i++; (*j) += 16;  }
    if(l3[x+1]) { i++; (*j) += 32;  }
    if(l3[x  ]) { i++; (*j) += 64;  }
    if(l3[x-1]) { i++; (*j) += 128; }

    return i;
}

/** \brief Calulate disparity
    \param limg     The input 16 bits grey left image.
    \param rimg 	The input 16 bits grey right image.
    \param ledg     The input 16 bits left edges.
    \param redg 	The input 16 bits right edges.
    \param out      The output 16 bits disparity.
    \param buff     The temporary buffer.
    \param w        The image width.
    \param h        The imahe height.
*/
void stereo_maching(const int16 *limg, const int16 *rimg, const int16 *ledg, const int16 *redg, int16 *out, int16 *buff, const int w, const int h)
{
    int i, j, y, x, y1, x1, yx, yxl, yxr, yw, yw1; //, sh1 = sh+1;
    int d = 100, f = 100, ths = 200;
    int sad, sd[4], sadt, z, bl, gp = 0, tp = 0, nv, nvt, dir, dirt;

    int sh = 2, ls = (sh<<1)+1, h1 = h-sh, w2 = w + (sh<<1), w1 = w*sh;
    int16 *l[ls], *r[ls], *le[ls], *re[ls], *tm;

    //For left image
    l[0] = buff;
    for(i=1; i < ls; i++) l[i] = &l[i-1][w2];
    //Prepare buffer
    for(i=0; i < ls-1; i++) cp_line_16(&limg[w*abs(sh-i)], l[i], w, sh);

    //For right image
    r[0] =  &l[ls-1][w2];
    for(i=1; i < ls; i++) r[i] = &r[i-1][w2];
    //Prepare buffer
    for(i=0; i < ls-1; i++) cp_line_16(&rimg[w*abs(sh-i)], r[i], w, sh);

    //For left edges
    le[0] = &r[ls-1][w2];
    for(i=1; i < ls; i++) le[i] = &le[i-1][w2];
    //Prepare buffer
    for(i=0; i < ls-1; i++) cp_line_16(&ledg[w*abs(sh-i)], le[i], w, sh);

    //For right edges
    re[0] =  &le[ls-1][w2];
    for(i=1; i < ls; i++) re[i] = &re[i-1][w2];
    //Prepare buffer
    for(i=0; i < ls-1; i++) cp_line_16(&redg[w*abs(sh-i)], re[i], w, sh);

    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + w1 : w*(((h-1)<<1)-y);
        cp_line_16(&limg[yw1], l[ls-1], w, sh);
        cp_line_16(&rimg[yw1], r[ls-1], w, sh);
        cp_line_16(&ledg[yw1], le[ls-1], w, sh);
        cp_line_16(&redg[yw1], re[ls-1], w, sh);

        for(x=0; x < w-d; x++){
            yx = yw + x;
            out[yx] = 0;
            if(re[2][x+sh]){
                nv = check_pixel(re[1], re[2], re[3], x+sh, &dir);
                tp++;
                printf("nv = %d dir = %d\n", nv, dir);
                sadt = 0;
                for(i=0; i < d; i++){
                    if(le[2][x+d+sh+i]){
                        sad = 0;
                        for(j=0; j < 4; j++) sd[j] = 0;
                        //Block matching
                        nv = check_pixel(le[1], le[2], le[3], x+d+sh+i, &dir);
                        for(y1=0; y1 < 5; y1++){
                            for(x1=0; x1 < 5; x1++){
                                yxl = x+d+i+x1;
                                yxr = x + x1;
                                //sad += abs(l[y1][yxl] - r[y1][yxr]);
                                // 0  0  0  1  1
                                // 0  0  0  1  1
                                // 3  3     1  1
                                // 3  3  2  2  2
                                // 3  3  2  2  2
                                if(y1 < 2 && x1 < 3) sd[0] += abs(l[y1][yxl] - r[y1][yxr]);
                                if(y1 < 3 && x1 > 2) sd[1] += abs(l[y1][yxl] - r[y1][yxr]);
                                if(y1 > 2 && x1 > 1) sd[2] += abs(l[y1][yxl] - r[y1][yxr]);
                                if(y1 > 1 && x1 < 2) sd[3] += abs(l[y1][yxl] - r[y1][yxr]);
                                //if(!i) printf(" %3d %3d  ", l[y1][yxl], r[y1][yxr]);
                            }
                            //if(!i) printf("\n");
                        }
                        sad = sd[0] + sd[1] + sd[2] + sd[3] + abs(l[2][x+d+i+2] - r[2][x + 2]);
                        //if(!i) printf("\n");
                        printf("sad = %d sd0 = %d sd1 = %d sd2 = %d sd3 = %d\n", sad/25, sd[0]/6, sd[1]/6, sd[2]/6, sd[3]/6);

                        if(!sadt) { sadt = sad; bl = x+d+i; dirt = dir; nvt = nv; }
                        if(sad < sadt) { sadt = sad; bl = x+d+i; dirt = dir; nvt = nv; }
                        //printf("sadt = %d\n", sadt);
                    }
                }
                //if(sadt) {
                if(sadt < ths && sadt ) {
                    out[yx] = d*f/(bl - x);
                    gp++;
                } else {
                    printf("Bag point\n");
                }
                printf("rx = %d lx = %d sad = %d nv = %d dir = %d\n", x, bl, sadt/25, nvt, dirt);

                //if(bl - x) out[yx] = d*f/(bl - x);
                //else out[yx] = 255;
                //printf("\nsad = %d xl = %d xr = %d dis = %d z = %d\n", sadt, bl, x, bl - x, z);
            }
        }

        tm = l[0];  for(i=0; i < ls-1; i++) l[i]  = l[i+1];  l[ls-1]  = tm;
        tm = r[0];  for(i=0; i < ls-1; i++) r[i]  = r[i+1];  r[ls-1]  = tm;
        tm = le[0]; for(i=0; i < ls-1; i++) le[i] = le[i+1]; le[ls-1] = tm;
        tm = re[0]; for(i=0; i < ls-1; i++) re[i] = re[i+1]; re[ls-1] = tm;
        //tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = l[3]; l[3] = l[4]; l[4] = l[5]; l[5] = l[6]; l[6] = tm;
    }
    printf("Total pixels = %d good pixele = %d  %d \n", tp, gp, gp*100/tp);
}

/** \brief Calulate disparity
    \param in	The input 16 bits first image.
    \param in1 	The input 16 bits second image.
    \param out	The output 16 bits disparity.
    \param buff	The temporary buffer.
    \param w	The image width.
    \param h	The imahe height.
*/
void stereo_disparity(const int16 *left, const int16 *right, int16 *out, int16 *buff, const int w, const int h)
{
    int i, size = w*h, th = 10;
    int16 *Y[2][2], *buf;

    Y[0][0] = buff;
    Y[0][1] = &buff[size];
    Y[1][0] = &buff[size*2];
    Y[1][1] = &buff[size*3];
    buf = &buff[size*4];

    trans_rgb_to_grey(left , Y[0][0], w, h);
    trans_rgb_to_grey(right, Y[1][0], w, h);

    seg_canny_edge(Y[0][0], Y[0][1], buf, w, h, th);
    seg_canny_edge(Y[1][0], Y[1][1], buf, w, h, th);

    seg_end_of_edges(Y[0][1], Y[0][1], buf, w, h);
    seg_end_of_edges(Y[1][1], Y[1][1], buf, w, h);

    stereo_maching(Y[0][0], Y[1][0], Y[0][1], Y[1][1], out, buf, w, h);

}
