#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    // s[0]  s[1]  s[2]
    //|-----|-----|-----|
    //|     |     |     |
    //|-----|-----|-----|
    //|     | yx  |     |
    //|-----|-----|-----|
    //|     |     |     |
    //|-----|-----|-----|
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
            if(x&1) i = 1;
            else i = 0;

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

