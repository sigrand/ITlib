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


/**	\brief	Gradient of 16 bits grey image.
    \param	in		The input 16 bit image.
    \param	out     The output 16 bit gradient image.
    \param	buff    The 3 lines  buffer.
    \param	w       The image width.
    \param  h       The image height.
    \param  th      The gradient threshould (if < th = 0).
*/
void seg_gradient(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int th)
{
    uint32 i, y, x, x2, yx, yw, yw1, h1 = h-1, sh = 1, w2 = w + (sh<<1); //, sh1 = sh+1;
    int16  *tm, *l[3];
    int g[4], max, min[2];

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
            /*
            g[0] = abs(l[1][x] - l[1][x+2]);    max = max < g[0] ? g[0] : max;
            g[1] = abs(l[0][x] - l[2][x+2]);    max = max < g[1] ? g[1] : max;
            g[2] = abs(l[0][x+1] - l[2][x+1]);  max = max < g[2] ? g[2] : max;
            g[3] = abs(l[0][x+2] - l[2][x]);    max = max < g[3] ? g[3] : max;

            max = (g[0] + g[1] + g[2] + g[3]);
            out[yx] = (max-th) > 0 ? (max > 251 ? 251 : max) : 0;
            */
            /*
            g[0] = abs(l[1][x] - l[1][x+2]);    min[0] = g[0];
            g[2] = abs(l[0][x+1] - l[2][x+1]);  min[0] = min[0] > g[2] ? g[2] : min[0];

            g[1] = abs(l[0][x] - l[2][x+2]);    min[1] = g[1];
            g[3] = abs(l[0][x+2] - l[2][x]);    min[1] = min[1] > g[3] ? g[3] : min[1];

            if(min[0]>th || min[1]>th) {
                max = min[0] > min[1] ? min[0] : min[1];
                out[yx] = max > 251 ? 251 : max;
                //printf("max = %d out = %d\n", max, out[yx]);
            } else out[yx] = 0;
            */
            max = 0;
            g[0] = abs(l[1][x]   - l[1][x+2]);  if(g[0] > th && max < g[0]) { max = g[0]; i = 1; }
            g[1] = abs(l[0][x]   - l[2][x+2]);  if(g[1] > th && max < g[1]) { max = g[1]; i = 2; }
            g[2] = abs(l[0][x+1] - l[2][x+1]);  if(g[2] > th && max < g[2]) { max = g[2]; i = 3; }
            g[3] = abs(l[0][x+2] - l[2][x]  );  if(g[3] > th && max < g[3]) { max = g[3]; i = 4; }

            out[yx] = max ? ((g[0]+g[1]+g[2]+g[3])>>2) : 0;


        }
        tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = tm;
    }
}

/**	\brief	Average direction filter.
    \param	l0		The pointer to 1 line.
    \param	l1		The pointer to 2 line.
    \param	l2		The pointer to 3 line.
    \param  x		Coordinate pixel in x direction.
    \retval			0 - if not
*/
static inline uint32  avr_dir(uint8 *l0, uint8 *l1, uint8 *l2, const int x)
{
    int tmp;
    tmp =   l0[x-1] + l0[x] + l0[x+1] +
            l1[x-1] + l1[x] + l1[x+1] +
            l2[x-1] + l2[x] + l2[x+1];
    return tmp/9;
}

/**	\brief	Canny edge detection.
    \param	in		The input 16 bit image.
    \param	out     The output 16 bit gradient image.
    \param	buff    The 3 lines  buffer.
    \param	w       The image width.
    \param  h       The image height.
    \param  th      The gradient threshould (if < th = 0).
*/
void seg_canny_edge(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int th)
{
    uint32 i, y, x, x2, yx, yx1, yw, yw1, h1 = h-1, sh = 1, w2 = w + (sh<<1); //, sh1 = sh+1;
    int16  *tm, *l[3], *gr[3];
    int g[4], max, min[2];
    uint8 *dr[2], *tm1, drt;

    //Buffer for image
    l[0] = buff;
    for(i=1; i < 3; i++) l[i] = &l[i-1][w2];
    //Buffer for gradient
    gr[0] = &l[2][w2];
    for(i=1; i < 3; i++) gr[i] = &gr[i-1][w2];
    //Direction line
    dr[0] = (uint8*)&gr[2][w2];
    //for(i=1; i < 3; i++) dr[i] = &dr[i-1][w2];
    dr[1] = &dr[0][w2];
    dr[2] = &dr[1][w2];

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
            //printf("x = %d y = %d\n", x, y);
            /*
            g = abs(l[1][x]   - l[1][x+2]);  if(g > th && max < g) { max = g; i = 1; }
            g = abs(l[0][x]   - l[2][x+2]);  if(g > th && max < g) { max = g; i = 2; }
            g = abs(l[0][x+1] - l[2][x+1]);  if(g > th && max < g) { max = g; i = 3; }
            g = abs(l[0][x+2] - l[2][x]  );  if(g > th && max < g) { max = g; i = 4; }
            gr[2][x] = max;
            */
            g[0] = abs(l[1][x]   - l[1][x+2]);  if(g[0] > th && max < g[0]) { max = g[0]; i = 1; }
            g[1] = abs(l[0][x]   - l[2][x+2]);  if(g[1] > th && max < g[1]) { max = g[1]; i = 2; }
            g[2] = abs(l[0][x+1] - l[2][x+1]);  if(g[2] > th && max < g[2]) { max = g[2]; i = 3; }
            g[3] = abs(l[0][x+2] - l[2][x]  );  if(g[3] > th && max < g[3]) { max = g[3]; i = 4; }

            gr[2][x] = max ? ((g[0]+g[1]+g[2]+g[3])>>2) : 0;
            dr[2][x] = i;


            //printf("x = %d y = %d gr = %d dir = %d \n", x, y, max, i);

            yx1 = yx-w-1;

            if(x > 1 && y > 1) {
                drt = dr[1][x-1];
                //Avarage direction filter
                //drt = avr_dir(dr[0], dr[1], dr[2], x-1);

                out[yx1] = 0;
                //out[yx-w-1] = gr[1][x-1];
                if       (drt == 1){
                    if(gr[1][x-2] < gr[1][x-1] && gr[1][x] < gr[1][x-1]) out[yx1] = gr[1][x-1];
                } else if(drt == 2){
                    if(gr[0][x-2] < gr[1][x-1] && gr[2][x] < gr[1][x-1]) out[yx1] = gr[1][x-1];
                } else if(drt == 3){
                    if(gr[0][x-1] < gr[1][x-1] && gr[2][x-1] < gr[1][x-1]) out[yx1] = gr[1][x-1];
                } else if(drt == 4){
                    if(gr[0][x] < gr[1][x-1] && gr[2][x-2] < gr[1][x-1]) out[yx1] = gr[1][x-1];
                }

                //Check for neighborhood
                // x    x x   x x    x
                // x x    x   x    x x
                /*
                if(out[yx1]) {
                    if     (out[yx1-1] && out[yx1-w-1]) out[yx1-1] = gr[1][x-1]; //out[yx1-1] = 0;
                    else if(out[yx1-w] && out[yx1-w-1]) out[yx1-w] = gr[1][x-1]; //out[yx1-1] = 0;
                    else if(out[yx1-w] && out[yx1-w+1]) out[yx1-w] = gr[1][x-1]; //out[yx1-1] = 0;
                    else if(out[yx1-1] && out[yx1-w  ]) out[yx1]   = gr[1][x-1]; //out[yx1]   = 0;
                }
                */

                if(out[yx1]) {
                    if     (out[yx1-1] && out[yx1-w-1]) out[yx1-1] = 0;
                    else if(out[yx1-w] && out[yx1-w-1]) out[yx1-w] = 0;
                    else if(out[yx1-w] && out[yx1-w+1]) out[yx1-w] = 0;
                    else if(out[yx1-1] && out[yx1-w  ]) out[yx1]   = 0;
                }

                //if(out[yx1]) out[yx1] = 255;
                //else out[yx1] = gr[1][x-1];

                /*
                printf("out[%d] = %d \n", yx-w-1, out[yx-w-1]);
                printf(" %3d  %3d  %3d\n %3d  %3d  %3d\n %3d  %3d  %3d\n\n",
                       gr[0][x-2], gr[0][x-1], gr[0][x],  gr[1][x-2], gr[1][x-1], gr[1][x],  gr[2][x-2], gr[2][x-1], gr[2][x]);
                printf(" %3d  %3d  %3d\n %3d  %3d  %3d\n %3d  %3d  %3d\n\n",
                       dr[0][x-2], dr[0][x-1], dr[0][x],  dr[1][x-2], dr[1][x-1], dr[1][x],  dr[2][x-2], dr[2][x-1], dr[2][x]);
                */
            }
        }
        tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = tm;
        tm =  gr[0]; gr[0] = gr[1]; gr[1] = gr[2]; gr[2] = tm;
        tm1 = dr[0]; dr[0] = dr[1]; dr[1] = dr[2]; dr[2] = tm1;
    }
}


/**	\brief	Corner detector of 16 bits grey image.
    \param	in		The input 16 bit image.
    \param	out     The output 16 bit gradient image.
    \param	buff    The 5 lines  buffer.
    \param	w       The image width.
    \param  h       The image height.
    \param  th      The gradient threshould (if < th = 0).
*/
void seg_corners(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int th)
{
    int i, j, y, x, x5, x2, yx, sum, yw, yw1,  sh = 2, h1 = h-sh, w2 = w + (sh<<1), w3 = w2<<1, w1 = w<<1, max, sz = w*h; //, sh1 = sh+1;
    int16 *l[5], *tm;
    int more, less, dif, mo;

    l[0] = buff;

    for(i=1; i < 5; i++) l[i] = &l[i-1][w2];

    //Prepare buffer
    cp_line_16(&in[w*2], l[0], w, sh);
    cp_line_16(&in[w  ], l[1], w, sh);
    cp_line_16(&in[0  ], l[2], w, sh);
    cp_line_16(&in[w  ], l[3], w, sh);

    //Copy image
    for(j=0; j < sz; j++) out[j] = in[j];

    i = 0;
    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + w1 : w*(((h-1)<<1)-y);
        cp_line_16(&in[yw1], l[4], w, sh);

        for(x=0; x < w; x++){
            yx = yw + x;
            //yx1 =  x + 2 + w3;
            x5 = x+5;
            x2 = x+2;
            sum = 0; more = 0; less = 0;
            /*
            for(i=x; i < x5; i++) sum += l[0][i] - l[2][x+2];
            for(i=x; i < x5; i++) sum += l[1][i] - l[2][x+2];
            for(i=x; i < x5; i++) sum += l[2][i] - l[2][x+2];
            for(i=x; i < x5; i++) sum += l[3][i] - l[2][x+2];
            for(i=x; i < x5; i++) sum += l[4][i] - l[2][x+2];
            */
            for(j=0; j < 5; j++){
                for(i=x; i < x5; i++){
                    dif = l[j][i] - l[2][x+2];
                    if(abs(dif) > th){
                        if(dif > 0) more++;
                        else less++;
                    }
                }
            }
            dif = l[0][0] - l[2][x+2];
            if(abs(dif) > th){
                if(dif > 0) { more++; mo = 1; }
                else { less++; mo = 0; }
            }
            dif = l[0][1] - l[2][x+2];
            if(abs(dif) > th){
                if(dif > 0) { more++; mo = 1; }
                else { less++; mo = 0; }
            }
            dif = l[1][0] - l[2][x+2];
            if(abs(dif) > th){
                if(dif > 0) { more++; mo = 1; }
                else { less++; mo = 0; }
            }


            //printf("more = %d less = %d\n", more, less);
            out[yx] = abs(more-less)*10;
        }
        tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = l[3]; l[3] = l[4]; l[4] = tm;
    }

    //printf("Local maxs = %d\n",i);
    //return i;
}

/**	\brief	Check if a pixel the local maximum.
    \param	l0		The pointer to 1 line.
    \param	l1		The pointer to 2 line.
    \param	l2		The pointer to 3 line.
    \param	l3		The pointer to 4 line.
    \param	l4		The pointer to 5 line.
    \param  x		Coordinate pixel in x direction.
    \retval			0 - if not
*/
static inline uint32 loc_max(int16 *l0, int16 *l1, int16 *l2, int16 *l3, int16 *l4, const int x)
{
    int i = 0;
    for(i=-2; i < 3; i++) if(l2[x-i] > l2[x]) return 0;
    for(i=-2; i < 3; i++) if(l0[x-i] > l2[x]) return 0;
    for(i=-2; i < 3; i++) if(l1[x-i] > l2[x]) return 0;
    for(i=-2; i < 3; i++) if(l3[x-i] > l2[x]) return 0;
    for(i=-2; i < 3; i++) if(l4[x-i] > l2[x]) return 0;
    if(l2[x]) l2[x] = 252;
    return l2[x];
}

/**	\brief	Find local maxumum in 5x5 window of 16 bits grey image.
    \param	in		The input 16 bit image.
    \param	out     The output 16 bit gradient image.
    \param	buff    The 5 lines  buffer.
    \param	w       The image width.
    \param  h       The image height.
*/
uint32 seg_local_max(int16 *in, int16 *out, int16 *buff, const int w, const int h)
{
    int i, j, y, x, yx, yx1, yw, yw1,  sh = 2, h1 = h-sh, w2 = w + (sh<<1), w3 = w2<<1, w1 = w<<1, max, sz = w*h; //, sh1 = sh+1;
    int16 *l[5], *tm;

    l[0] = buff;

    for(i=1; i < 5; i++) l[i] = &l[i-1][w2];

    //Prepare buffer
    cp_line_16(&in[w*2], l[0], w, sh);
    cp_line_16(&in[w  ], l[1], w, sh);
    cp_line_16(&in[0  ], l[2], w, sh);
    cp_line_16(&in[w  ], l[3], w, sh);

    //Copy image
    for(j=0; j < sz; j++) out[j] = in[j];

    i = 0;
    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + w1 : w*(((h-1)<<1)-y);
        cp_line_16(&in[yw1], l[4], w, sh);

        for(x=0; x < w; x++){
            yx = yw + x;


            //max = loc_max(l[0], dr, yx1, w2);
            max = loc_max(l[0], l[1], l[2], l[3], l[4], x+sh);
            if(max) {
                out[yx] = 252;
                i++;
                //x+=2;
            }// else out[yx] = in[yx];
        }
        tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = l[3]; l[3] = l[4]; l[4] = tm;
    }

    printf("Local maxs = %d\n",i);
    return i;
}

/**	\brief	Check if a pixel the end of edge.
    \param	l0		The pointer to 1 line.
    \param	l1		The pointer to 2 line.
    \param	l2		The pointer to 3 line.
    \param  x		Coordinate pixel in x direction.
    \retval			0 - if not
*/
static inline uint32 end_of_edge(int16 *l0, int16 *l1, int16 *l2, const int x)
{
    int i = 0;
    /*
    if(l1[x-1]) i++;
    if(i<2 && l0[x-1]) i++;
    else return 0;
    if(i<2 && l0[x]) i++;
    else return 0;
    if(i<2 && l0[x+1]) i++;
    else return 0;
    if(i<2 && l1[x+1]) i++;
    else return 0;
    if(i<2 && l2[x+1]) i++;
    else return 0;
    if(i<2 && l2[x]) i++;
    else return 0;
    if(i<2 && l2[x-1]) i++;
    else return 0;
    if(i>1) return 0;
    */
    if(l1[x-1]) i++;
    if(l0[x-1]) i++;
    if(l0[x  ]) i++;
    if(l0[x+1]) i++;
    if(l1[x+1]) i++;
    if(l2[x+1]) i++;
    if(l2[x  ]) i++;
    if(l2[x-1]) i++;

    if(i == 1) return 1;
    else return 0;
}


/**	\brief	Find ends of edges of 16 bits grey image.
    \param	in		The input 16 bit image.
    \param	out     The output 16 bit gradient image.
    \param	buff    The 5 lines  buffer.
    \param	w       The image width.
    \param  h       The image height.
*/
uint32 seg_end_of_edges(int16 *in, int16 *out, int16 *buff, const int w, const int h)
{
    int i, j, y, x, yx, yx1, yw, yw1,  sh = 1, h1 = h-sh, w2 = w + (sh<<1), w1 = w<<(sh-1), max, sz = w*h; //, sh1 = sh+1;
    int16 *l[3], *tm;

    l[0] = buff;

    for(i=1; i < 3; i++) l[i] = &l[i-1][w2];

    //Prepare buffer
    //cp_line_16(&in[w*2], l[0], w, sh);
    cp_line_16(&in[w  ], l[0], w, sh);
    cp_line_16(&in[0  ], l[1], w, sh);
    //cp_line_16(&in[w  ], l[3], w, sh);

    //Copy image
    for(j=0; j < sz; j++) out[j] = in[j];

    i = 0;
    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + w1 : w*(((h-1)<<1)-y);
        cp_line_16(&in[yw1], l[2], w, sh);

        for(x=0; x < w; x++){
            yx = yw + x;

            if(l[1][x+sh]) {
                if(end_of_edge(l[0], l[1], l[2], x+sh)) { out[yx] = 252; i++; }
            }

        }
        tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = tm;
    }

    printf("End of edges = %d\n",i);
    return i;
}

static inline uint32 check_max(int16 *img, int *dr, int yx, uint8 val)
{
    if(img[yx+dr[0]] == val) return 1;
    if(img[yx+dr[1]] == val) return 1;
    if(img[yx+dr[2]] == val) return 1;
    if(img[yx+dr[3]] == val) return 1;
    if(img[yx+dr[4]] == val) return 1;
    if(img[yx+dr[5]] == val) return 1;
    if(img[yx+dr[6]] == val) return 1;
    if(img[yx+dr[7]] == val) return 1;
    return 0;
}

static inline uint32 check_diagonal(uint8 *con, int *dr, int yx, uint8 di)
{
    if      (di == 1 && con[yx+dr[0]]&8   && con[yx+dr[2]]&128) return 1;
    else if (di == 3 && con[yx+dr[2]]&32  && con[yx+dr[4]]&2  ) return 1;
    else if (di == 5 && con[yx+dr[4]]&128 && con[yx+dr[6]]&8  ) return 1;
    else if (di == 7 && con[yx+dr[6]]&2   && con[yx+dr[0]]&32 ) return 1;
    return 0;
}

static inline void add_dir(uint8 *con, int *dr, int yx, int d, int w)
{
    //vx->n++;
    if      (d == dr[0]) { con[yx] |= 1;     con[yx + d] |= 16; }
    else if (d == dr[1]) { con[yx] |= 2;     con[yx + d] |= 32; }
    else if (d == dr[2]) { con[yx] |= 4;     con[yx + d] |= 64; }
    else if (d == dr[3]) { con[yx] |= 8;     con[yx + d] |= 128;}
    else if (d == dr[4]) { con[yx] |= 16;    con[yx + d] |= 1;  }
    else if (d == dr[5]) { con[yx] |= 32;    con[yx + d] |= 2;  }
    else if (d == dr[6]) { con[yx] |= 64;    con[yx + d] |= 4;  }
    else if (d == dr[7]) { con[yx] |= 128;   con[yx + d] |= 8;  }
    //printf("finish d = %d dir = %o\n", d, dir);
}

static inline uint32 find_dir(int *dr, int in)
{
    if(dr[0] == in) return 0;
    if(dr[1] == in) return 1;
    if(dr[2] == in) return 2;
    if(dr[3] == in) return 3;
    if(dr[4] == in) return 4;
    if(dr[5] == in) return 5;
    if(dr[6] == in) return 6;
    if(dr[7] == in) return 7;
}

/**	\brief	Find the maximum around the pixel.
    \param	img		The pointer to gradient image.
    \param	yx		The pixel coordinate (yx = y*w + x)
    \param	w		The image width.
    \param  in1		The previous maximum direction.
    \retval			The direction of of local max.
*/
static inline int dir(uint8 *con, int16 *grad, int *dr, int yx, int in1, int w)
{
    uint8 max = 0,  i;
    int in = 0, tmp[3];
    /*
    if(yx == 812-1 + (542+1)*w) {
        printf("befor dir3: x = %d y = %d in1 = %d\n", yx%w, yx/w, in1);
        print_around(con, yx, w);
        print_around(grad, yx, w);

        print_con_grad(grad, con, yx, w);
    }*/
    //x = 1077 y = 1174
    if(in1 == 0){
        if(grad[yx+dr[0]] > max) { max = grad[yx+dr[0]]; in = dr[0]; }
        if(grad[yx+dr[1]] > max) { max = grad[yx+dr[1]]; in = dr[1]; }
        if(grad[yx+dr[2]] > max) { max = grad[yx+dr[2]]; in = dr[2]; }
        if(grad[yx+dr[3]] > max) { max = grad[yx+dr[3]]; in = dr[3]; }
        if(grad[yx+dr[4]] > max) { max = grad[yx+dr[4]]; in = dr[4]; }
        if(grad[yx+dr[5]] > max) { max = grad[yx+dr[5]]; in = dr[5]; }
        if(grad[yx+dr[6]] > max) { max = grad[yx+dr[6]]; in = dr[6]; }
        if(grad[yx+dr[7]] > max) { max = grad[yx+dr[7]]; in = dr[7]; }
        return in;
    }

    if      (in1 == dr[0]){ tmp[0] = dr[3]; tmp[1] = dr[4]; tmp[2] = dr[5]; }
    else if (in1 == dr[1]){ tmp[0] = dr[4]; tmp[1] = dr[5]; tmp[2] = dr[6]; }
    else if (in1 == dr[2]){ tmp[0] = dr[5]; tmp[1] = dr[6]; tmp[2] = dr[7]; }
    else if (in1 == dr[3]){ tmp[0] = dr[6]; tmp[1] = dr[7]; tmp[2] = dr[0]; }
    else if (in1 == dr[4]){ tmp[0] = dr[7]; tmp[1] = dr[0]; tmp[2] = dr[1]; }
    else if (in1 == dr[5]){ tmp[0] = dr[0]; tmp[1] = dr[1]; tmp[2] = dr[2]; }
    else if (in1 == dr[6]){ tmp[0] = dr[1]; tmp[1] = dr[2]; tmp[2] = dr[3]; }
    else if (in1 == dr[7]){ tmp[0] = dr[2]; tmp[1] = dr[3]; tmp[2] = dr[4]; }
    //printf("tmp0 = %d tmp1 = %d tmp2 = %d\n", tmp[0], tmp[1], tmp[2]);
    max = 0;
    if(grad[yx+tmp[0]] > max) { max = grad[yx+tmp[0]]; in = tmp[0]; }
    if(grad[yx+tmp[1]] > max) { max = grad[yx+tmp[1]]; in = tmp[1]; }
    if(grad[yx+tmp[2]] > max) { max = grad[yx+tmp[2]]; in = tmp[2]; }
    //printf("max = %d in = %d\n", max, in);

    if(con[yx+in]){
        if(check_diagonal(con, dr, yx, find_dir(dr, in)) ){
            if(in == tmp[0]) tmp[0] = tmp[2];
            else if(in == tmp[1]) tmp[1] = tmp[2];
            //printf("Remove diagonal 1 \n ");
            //print_con_grad(grad, con, yx, w);
            //printf("First ");
        } else if(check_max(grad, dr, yx+in, 255) ){
            if(in == tmp[0]) tmp[0] = tmp[2];
            else if(in == tmp[1]) tmp[1] = tmp[2];
        } else return in;

    } else {
        return in;
    }

    if(grad[yx+tmp[0]] > grad[yx+tmp[1]])   in = tmp[0];
    else                                    in = tmp[1];

    if(!check_diagonal(con, dr, yx, find_dir(dr, in))) {
        if(con[yx+in]){
            //if(max == 255) return in;
            if(check_max(grad, dr, yx+in, 255)){
                if(in == tmp[0]) in = tmp[1];
                else in = tmp[0];
                //printf("Second ");
            } else return in;
        } else {
            return in;
        }
    } else {
        if(in == tmp[0]) in = tmp[1];
        else in = tmp[0];
    }

    if(!check_diagonal(con, dr, yx, find_dir(dr, in))) {
        if(con[yx+in]){
            //if(max == 255) return in;
            if(check_max(grad, dr, yx+in, 255)){
                return 0;
                //printf("Second ");
            } else return in;
        } else {
            return in;
        }
    } else return 0;
}

/**	\brief	Edge detection.
    \param	in		The input 16 bit gradient image with local maximum.
    \param	out     The output 16 bit gradient image.
    \param	w       The image width.
    \param  h       The image height.
*/
void seg_edge_detection(int16 *in, uint8 *con, const int w, const int h)
{
    uint32 y, y1, x, yx, yw, yx1, yx2, yx3, yxt, i, h1 = h-1, w1 = w-1, max;//, yxl = 0;
    int d1, d2, d3, npix = 0, cr;
    int cn, ck = 0;
    int dr[8] = { -1, -1-w, -w, +1-w, 1, 1+w, w, -1+w };

    for(y=1; y < h1; y++){
        yw = y*w;
        for(x=1; x < w1; x++){
            yx = yw + x;
            if(in[yx] == 252){
                //printf("x = %d y = %d\n", x, y);
                //yx = hist[i];
                yx1 = yx; cr = 0; cn = 0;
                //con[yx1] = 255;
                //Store first search direction
                d1 = dir(con, in, dr, yx1, 0, w);
                //Store second search direction
                //d2 = dir3(con, in, dr, yx1, d1, w);
                d2 = d1;
                if(!con[yx1+d1]){
                    in[yx1] = 253;
                    //Start first search direction
                    while(!cr && cn < 2){
                        while(1){
                            yx1 = yx1 + d1;
                            if(con[yx1]) {
                                if(d1) add_dir(con, dr, yx1, -d1, w);
                                in[yx1] = 255; //con[yx1] = 255;
                                npix++;

#ifdef DEBUG
                                /*
                                if(check_vertex(in, yx1, 255, w)) {
                                    printf("seg_find_intersect7: x = %d y = %d yx = %d d1 = %d\n", yx1%w, yx1/w, yx1, d1);
                                    if((yx1%w) > 1 && (yx1/w > 1)){
                                        print_around(con, yx1, w);
                                        print_around(in, yx1, w);
                                        print_con_in(in, con, yx1, w);
                                    }
                                }
                                */
#endif
                                //Check for circle
                                if(yx == yx1) cr = 1;
                                break;
                            }
                            //add_dir1(&di[yxt], d1, w);
                            add_dir(con, dr, yx1, -d1, w);
                            in[yx1] = 253; //con[yx1] = 255;
                            //d3 = d1; //yxc = 0;
                            d1 = dir(con, in, dr, yx1, -d1, w);

                            if(!d1) {
                                in[yx1] = 255; //con[yx1] = 255;
                                npix++;
                                break;
                            }
                        }
                        cn++;
                        //d1 = d2; yx1 = yx;
                        yx1 = yx;
                        d1 = dir(con, in, dr, yx1, d2, w);
                        //printf("Change dir\n");
                    }
                }
            }
        }
    }
    printf("Numbers of intersection  = %d\n", npix);

}
