#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./segmentation.h"
#include "../libit/types.h"
#include "../filters/filters.h"

void inline static cp_line_16(int16 *in, int16 *l, uint32 w, uint32 sh)
{
    uint32 i;
    for(i=0; i < sh; i++) l[i] = in[sh-i];
    for(i=0; i < w ; i++) l[i+sh] = in[i];
    for(i=0; i < sh; i++) l[i+sh+w] = in[w-sh-i];
}

/**	\brief Get the differnce between 16 bits  original and averaged bayer image.
    \param in       The input 16 bits bayer image.
    \param avr      The input 16 bits averaged bayer image.
    \param dif      The output image.
    \param w        The image width.
    \param h        The image height.
    \param bpp      The input images bits per pixel.
*/
void seg_grad3(uint8 *img, uint8 *img1, uint8 *con, uint8 *di, uint32 w, uint32 h, int th)
{
    /// | |x| |      | | | |      |x| | |      | | |x|
    /// | |x| |      |x|x|x|      | |x| |      | |x| |
    /// | |x| |      | | | |      | | |x|      |x| | |
    ///  g[2]         g[0]         g[1]         g[3]
    /// Direction
    ///   n=0          n=2         n=3          n=1
    /// | | | |      | ||| |      | | |/|      |\| | |
    /// |-|-|-|      | ||| |      | |/| |      | |\| |
    /// | | | |      | ||| |      |/| | |      | | |\|
    /// 255 - intersection
    /// 254 - vertex
    /// 253 - direction curve
    /// 252 - local max

    uint32 y, x, yx, yw, w1 = w-1, h1 = h-1, w2 = w-2, h2 = h-2;
    uint8 in, col = 253;
    uint32 g[8];
    int max;

    //The up border
    img1[w + 1] = 255; for(x=2; x < w2; x++) img1[w + x] = col; img1[w + x] = 255;
    //con[w+1] = 255; for(x=2; x < w2; x++) con[w + x] = 64; con[w + x] = 255;
    //for(x=1; x < w2; x++) di[w + x] = 8; di[w + x] = 2;
    con[w + 1] = 80; for(x=2; x < w2; x++) con[w + x] = 17; con[w + x] = 65;

    for(y=2; y < h2; y++){
        yw = y*w;
        img1[yw + 1] = col; con[yw + 1] = 68;  //di[yw + 1] = 4; //di[yw + 1] = 32;
        for(x=2; x < w2; x++){
            yx = yw + x;

            g[0] = abs(img[yx-1  ] - img[yx+1  ]);
            g[1] = abs(img[yx-1-w] - img[yx+1+w]);
            g[2] = abs(img[yx-w  ] - img[yx+w  ]);
            g[3] = abs(img[yx+1-w] - img[yx-1+w]);

            /*
                        g[0] = abs(img[yx-2  ] - img[yx+2  ]);
                        g[1] = abs(img[yx-2-(w<<1)] - img[yx+2+(w<<1)]);
                        g[2] = abs(img[yx-(w<<1)  ] - img[yx+(w<<1)  ]);
                        g[3] = abs(img[yx+2-(w<<1)] - img[yx-2+(w<<1)]);
                        */

            /*
                        g[4] = abs(img[yx-2-w  ] - img[yx+2+w  ]);
                        g[5] = abs(img[yx-1-(w<<1)] - img[yx+1+(w<<1)]);
                        g[6] = abs(img[yx+1-(w<<1)] - img[yx-1+(w<<1)]);
                        g[7] = abs(img[yx+2-w] - img[yx-2+w]);
                        */
            max = (g[0] + g[1] + g[2] + g[3])>>1;
            //max = (g[0] + g[1] + g[2] + g[3] + g[4] + g[5] + g[6] + g[7])>>2;
            img1[yx] = (max-th) > 0 ? (max > 251 ? 251 : max) : 0;
            //printf("img = %d max = %d th = %d max-th = %d\n", img1[yx], max, th, max-th);
            //img1[yx] = max>>th;
            //max = (((g[0] + g[1] + g[2] + g[3])>>2)>>th)<<th;
            //img1[yx] = max > 252 ? 252 : max;
        }
        img1[yx + 1] = col; con[yx + 1] = 68; //di[yx + 1] = 64; //di[yx + 1] = 2;
        //img1[yx + 2] = col1;
    }
    //The bottom border
    yw = y*w;
    img1[yw + 1] = 255; for(x=2; x < w2; x++) img1[yw + x] = col; img1[yw + x] = 255;
    con[yw + 1] = 20; for(x=2; x < w2; x++) con[yw + x] = 17; con[yw + x] = 5;
    //di[yw+1] = 32; for(x=2; x < w1; x++) di[yw + x] = 128;
    //di[yw+1] = 4; for(x=2; x < w1; x++) di[yw + x] = 1;
    printf("Finish gradient\n");
}

/**	\brief	Gradient of 16 bits grey image.
    \param	in		The input 16 bit image.
    \param	out     The output 16 bit gradient image.
    \param	buff    The 3 lines  buffer.
    \param	w       The image width.
    \param  h       The image height.
    \param  th      The gradient threshould (if < th = 0).
*/
void seg_grad(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int th)
{
    uint32 i, y, x, x2, yx, yw, yw1, h1 = h-1, sh = 1, w2 = w + (sh<<1); //, sh1 = sh+1;
    int16  *tm, max, *l[3];
    int g[4];

    l[0] = buff;
    for(i=1; i < 3; i++) l[i] = &l[i-1][w2];

    //Prepare buffer
    cp_line_16(&in[w], l[0], w, sh);
    cp_line_16(&in[0], l[1], w, sh);

    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + w : yw - w;
        cp_line_16(&in[yw1], l[2], w, sh);

        for(x=0; x < w; x++){
            yx = yw + x;
            //x2 = x + sh1;
            max = 0;
            g[0] = abs(l[1][x] - l[1][x+2]);    max = max < g[0] ? g[0] : max;
            g[1] = abs(l[0][x] - l[2][x+2]);    max = max < g[1] ? g[1] : max;
            g[2] = abs(l[0][x+1] - l[2][x+1]);  max = max < g[2] ? g[2] : max;
            g[3] = abs(l[0][x+2] - l[2][x]);    max = max < g[3] ? g[3] : max;

            max = (g[0] + g[1] + g[2] + g[3]);
            out[yx] = (max-th) > 0 ? (max > 251 ? 251 : max) : 0;
        }
        tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = tm;
    }
}

/**	\brief	Check is a pixel the local maximum.
    \param	img		The pointer to gradient image.
    \param	dr		The pointer to search direction.
    \param	yx		The pixel coordinate (yx = y*w + x)
    \param  w		The image width.
    \retval			0 - if not
*/
static inline uint32 loc_max(const int16 *img, const int *dr, const int yx, const int w)
{
    uint32 i = 0;
    for(i=0; i < 24; i++) if(img[yx+dr[i]] > img[yx]) return 0;
    return img[yx];
}

/**	\brief	Gradient of 16 bits grey image.
    \param	in		The input 16 bit image.
    \param	out     The output 16 bit gradient image.
    \param	buff    The 3 lines  buffer.
    \param	w       The image width.
    \param  h       The image height.
    \param  th      The local max threshould (if < th = 0).
*/
uint32 seg_local_max(int16 *in, int16 *out, int16 *buff, uint32 w, uint32 h, uint32 th)
{
    //Direction dr[24]
    //
    //  8  9 10 11 12
    // 13  0  1  2 14
    // 15  7     3 16
    // 17  6  5  4 18
    // 19 20 21 22 23
    //
    //int w2 = w<<1, max;
    /*
    int dr[24] = {   -1-w, -w, +1-w, 1, 1+w, w, -1+w, -1,
                    -2-w2, -1-w2, -w2, 1-w2, 2-w2,
                    -2-w, 2-w,
                    -2, 2,
                    -2+w, 2+w,
                    -2+w2, -1+w2, w2, 1+w2, 2+w2
                 };
                 */
    uint32 i, y, x, x2, yx, yw, yw1,  sh = 2, h1 = h-sh, w2 = w + (sh<<1), w1 = w<<1, max; //, sh1 = sh+1;

    int16 *l[5], *tm;

    l[0] = buff;

    for(i=1; i < 5; i++) l[i] = &l[i-1][w2];

    //Prepare buffer
    cp_line_16(&in[w*2], l[0], w, sh);
    cp_line_16(&in[w  ], l[1], w, sh);
    cp_line_16(&in[0  ], l[2], w, sh);
    cp_line_16(&in[w  ], l[3], w, sh);

    i = 0;
    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + w1 : w*(((h-1)<<1)-y);
        cp_line_16(&in[yw1], l[2], w, sh);

        for(x=0; x < w; x++){
            yx = yw + x;
            //x2 = x + sh1;
            //max = loc_max(in, dr, yx, w);
            //if(l[2][x+2]


            if(max) {
                out[yx] = 252;
                i++;
            } else out[yx] = in[yx];
        }
        tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = l[3]; l[3] = l[4]; l[4] = tm;
    }

    printf("Local maxs = %d\n",i);
    return i;
}
