#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../libit/types.h"
#include "./stereo.h"
#include "../transforms/transforms.h"
#include "../segmentation/segmentation.h"


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
        }
    }
    return sad;
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
    int i, j, y, x, yx, yw, yw1; //, sh1 = sh+1;
    int d = 100, sad; //Cameras distance in pixels

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
            if(re[2][x+sh]){
                for(i=0; i < d; i++){
                    if(le[2][x+d+sh+i]){
                        sad = block_maching(l[0], r[0], x+d+sh+i, x+sh, w2);
                    }
                }
            }
        }

        tm = l[0];  for(i=0; i < ls-1; i++) l[i]  = l[i+1];  l[ls-1]  = tm;
        tm = r[0];  for(i=0; i < ls-1; i++) r[i]  = r[i+1];  r[ls-1]  = tm;
        tm = le[0]; for(i=0; i < ls-1; i++) le[i] = le[i+1]; le[ls-1] = tm;
        tm = re[0]; for(i=0; i < ls-1; i++) re[i] = re[i+1]; re[ls-1] = tm;
        //tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = l[3]; l[3] = l[4]; l[4] = l[5]; l[5] = l[6]; l[6] = tm;
    }
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

    for(i=0; i < size; i++){
        out[i] = Y[0][1][i] + Y[1][1][i];
       // out[i] = Y[0][i];
        //out[i] = in[i];
    }

}
