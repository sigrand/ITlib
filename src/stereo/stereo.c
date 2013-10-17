#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/**	\brief	Find abs differnce with two block 7x7 pixels size.
    \param	l		The pointer to left block.
    \param	r		The pointer to right block.
    \param	sd		The pointer to direction array.
    \param  xl		Coordinate of the center point of left block.
    \param  xr		Coordinate of the center point of right block.
    \retval			The summ of differnce.
*/
static inline uint32 block_maching(int16 **l, int16 **r, int *sd, const int xl, const int xr)
{
    int i, y, x, yxl, yxr, sad = 0, sad1 = 0;
    int ia[7][7] = {{0, 0, 0, 1, 1, 1, 2},
                    {7, 0, 0, 1, 1, 2, 2},
                    {7, 7, 0, 1, 2, 2, 2},
                    {7, 7, 7, 8, 3, 3, 3},
                    {6, 6, 6, 5, 4, 3, 3},
                    {6, 6, 5, 5, 4, 4, 3},
                    {6, 5, 5, 5, 4, 4, 4}};

    for(i=0; i < 9; i++) sd[i] = 0;
    for(y=0; y < 7; y++){
        for(x=0; x < 7; x++){
            yxl = x + xl;
            yxr = x + xr;
            sad1 += abs(l[y][yxl] - r[y][yxr]);
            sd[ia[y][x]] += abs(l[y][yxl] - r[y][yxr]);
            //if(!i) printf(" %3d %3d  ", l[y1][yxl], r[y1][yxr]);
        }
        //if(!i) printf("\n");
    }
    for(i=0; i < 9; i++) sad += sd[i];
    if(sad != sad1) printf(" sad = $d sad1 = %d\n", sad, sad1);
    return sad;
}

/**	\brief	Check if not only horizontal edge direction
    \param	l		The pointer to edge left block.
    \retval			The 1 if ok.
*/
static inline uint32 check_hdir(int16 **l, const int xb)
{
    int i, y, x, s[7], ed = xb+7;

    for(y=0; y < 7; y++){
        s[y] = 0;
        for(x=xb; x < ed; x++){
            s[y] += l[y][x];
        }
    }

    return 1;

    if(s[0] || s[1] || s[2] || s[4] || s[5] || s[6]) return 1;
    //if(s[0] || s[1] || s[2] || s[5] || s[6]) return 1;
    /*
    for(y=0; y < 7; y++){
        for(x=xb; x < ed; x++){
            printf("%3d ", l[y][x]);
        }
        printf("\n");
    }
    printf("\n");
    */
    return 0;
}

/**	\brief	Find best maching direction.
    \param	sd		The pointer to direction array.
    \param	ths		The maching threshould.
    \retval			The minimal of sam.
*/
static inline uint32 find_direction(int *sd)
{
    int i, j, sum, min;
    sum = sd[4] + sd[5] + sd[6] + sd[7] + sd[8]; //
    min = sum;
    for(i=0; i < 8; i++){
        j = i - 4;
        j = j < 0 ? 8 + j : j;
        //printf("sum = %d j = %d i = %d\n", sum, j, i);
        sum  = sum - sd[j] + sd[i];
        if(sum < min)  min = sum;
    }
    return min;
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

/**	\brief	Get range for left image matching
    \param	st		The left image start pixel.
    \param	en		The left image end pixel.
    \param	R		The maximum setereo distance.
    \param	r		The maximum stereo distance.
    \param  f		The focal length.
    \param  d       The distance between sensors.
    \param  w       The sensor width.
    \param  pw      The pixel width.
    \retval			The start pixel for the left image for matching.
*/
static uint32 get_left_array(int16 *st, int16 *en, double R, double r, double f, double d, double w, double pw)
{
    double a, b, a1, b1, c1;
    double x1, y1, x2, y2, xf, xn, yf, yn;
    int x, xr, xl;

    x1 = -d/2.;         y1 = 0;
    x2 = (w - d)/2.;  y2 = f;

    a = (y1 - y2)/(x1 - x2);
    b = y1 - a*x1;

    a1 = (a*a + 1);
    b1 = 2.*a*b;
    c1 = b*b - R*R;

    x1 = (-b1 + sqrt(b1*b1 - 4.*a1*c1))/(2.*a1);
    y1 = sqrt(R*R - x1*x1);
    x2 = d/2.; y2 = 0;

    a = (y1 - y2)/(x1 - x2);
    b = y1 - a*x1;

    xr = (int)(((f - b)/a - (d-w)/2.)/pw);
    xl = (int)((-(f - b)/a + (d+w)/2.)/pw);
    printf("xr = %d xl = %d\n", xr, xl);

    for(x=0; x < xr; x++){
        x1 = (d - w)/2. + (double)x*pw;
        y1 = f;
        x2 = d/2; y2 = 0;

        a = (y1 - y2)/(x1 - x2);
        b = y1 - a*x1;

        a1 = (a*a + 1);
        b1 = 2.*a*b;
        c1 = b*b - R*R;
        //Max range pixel
        xf = (-b1 - sqrt(b1*b1 - 4.*a1*c1))/(2.*a1);
        yf = sqrt(R*R - xf*xf);
        //The near pixel
        c1 = b*b - r*r;
        xn = (-b1 - sqrt(b1*b1 - 4.*a1*c1))/(2.*a1);
        yn = sqrt(r*r - xn*xn);

        x2 = -d/2.; y2 = 0;
        a = (yf - y2)/(xf - x2);
        b = yf - a*xf;

        st[x] = (int)(((f - b)/a + (d+w)/2.)/pw);

        a = (yn - y2)/(xn - x2);
        b = yn - a*xn;

        en[x] = (int)(((f - b)/a + (d+w)/2.)/pw);
        printf("x = %d sl = %d el = %d dif = %d  %d\n", x, st[x], en[x], en[x] - st[x], st[x] - x);
    }

    return xl;
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
    int i, j, y, x, xd, y1, y2, x1, yx, yxl, yxr, yw, yw1; //, sh1 = sh+1;
    int d = 100, ds = d>>2, f = 100, th = 5, ths = 49*th, thb = 6*4*th + th, size = w*h, f1;
    int sad, sad1, sd[9], sadt, z, bl, gp = 0, rp = 0, tp = 0, nv, nvt, dir, dirt;


    int sh = 3, ls = (sh<<1)+1, h1 = h-sh, w2 = w + (sh<<1), w1 = w*sh;
    int16 *l[ls], *r[ls], *le[ls], *re[ls], *tm;

    for(i=0; i < size; i++) out[i] = 0;

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
        //for(x=w; x > d ; x--){
            yx = yw + x;
            //check_hdir
            //if(le[sh][x+sh] && check_hdir(&le[0], x)){
            if(re[sh][x+sh] && check_hdir(&re[0], x)){
                //nv = check_pixel(re[sh-1], re[sh], re[sh+1], x+sh, &dir);
                tp++;
                //printf("nv = %d dir = %d\n", nv, dir);
                sadt = 0;
                //for(i=-ds; i < d; i++){
                    //if(le[sh][x+d+sh+i]){
                    //xd = x+d+i;
                f1 = (int)sqrt((double)((x - (w>>1))*(x - (w>>1)) + f*f));
                printf("x = %d x1 = %d x2 = %d\n", x, x+((d*f1)>>9), x+((d*f1)>>4));
                for(i=((d*f1)>>9); i < ((d*f1)>>4); i++){
                    //xd = x-d+sh-i;
                    //if(re[sh][x-d+sh-i] && check_hdir(&re[0], x-d-i)){
                    if(le[sh][x+sh+i] && check_hdir(&le[0], x+i)){
                        sad = 0; sad1 = 0;
                        //Block matching
                        j=0; {
                        //for(j=-1; j < 2; j++){
                            //xd = x-d-i+j;
                            //sad = block_maching7( &l[0], &r[0], sd, x, xd);
                            xd = x+i+j;
                            sad = block_maching( &l[0], &r[0], sd, xd, x);
                            sad = find_direction(sd);
                            //if(!i) printf("\n");
                            //if(y==1) printf("sad = %d sd0 = %d sd1 = %d sd2 = %d sd3 = %d\n", sad/25, sd[0]/6, sd[1]/6, sd[2]/6, sd[3]/6);
                            //if(y==1) printf("sad = %d sad1 = %d diff = %d\n", sad, sad1, sad-sad1);

                            if(!sadt) { sadt = sad; bl = xd; }
                            if(sad < sadt) { sadt = sad; bl = xd; }
                        }
                    }
                }
                /*
                for(i=-d+1; i < d-1; i++){
                    //xd = x-d+sh-i;
                    //if(re[sh][x-d+sh-i] && check_hdir(&re[0], x-d-i)){
                    if(le[sh][x+d+sh+i] && check_hdir(&le[0], x+d+i)){
                        sad = 0; sad1 = 0;
                        //Block matching
                        j=0; {
                        //for(j=-1; j < 2; j++){
                            //xd = x-d-i+j;
                            //sad = block_maching7( &l[0], &r[0], sd, x, xd);
                            xd = x+d+i+j;
                            sad = block_maching7( &l[0], &r[0], sd, xd, x);
                            sad = find_direction(sd);
                            //if(!i) printf("\n");
                            //if(y==1) printf("sad = %d sd0 = %d sd1 = %d sd2 = %d sd3 = %d\n", sad/25, sd[0]/6, sd[1]/6, sd[2]/6, sd[3]/6);
                            //if(y==1) printf("sad = %d sad1 = %d diff = %d\n", sad, sad1, sad-sad1);

                            if(!sadt) { sadt = sad; bl = xd; }
                            if(sad < sadt) { sadt = sad; bl = xd; }
                        }
                    }
                }
                */
                if(sadt) {
                    //Check vertical direction
                    if(sadt < thb){
                        out[yw + ((x + bl)>>1)] = d*f/abs(bl - x);
                        printf("bl = %d x = %d d = %d out = %d sadt = %d\n", bl, x, d, out[yw + ((x + bl)>>1)], sadt);
                        gp++;
                    /*} else if(find_direction(sd) < thb){
                        out[yx] = d*f/abs(bl - x);
                        rp++;
                        */
                    } /*else {
                        //Check around

                        sad = block_maching7( &l[0], &r[0], sd, x, xd-1);
                        if(sad < sadt) {
                            sadt = sad;
                            printf("left  sadt = %d sad = %d\n", sadt, sad);
                        }
                        sad = block_maching7( &l[0], &r[0], sd, x, xd+1);
                        if(sad < sadt) {
                            sadt = sad;
                            printf("Right sadt = %d sad = %d\n", sadt, sad);
                        }
                        //One more time vheck
                        if(sadt < ths){
                            out[yx] = d*f/abs(bl - x);
                            gp++;
                        } else {
                            if(find_direction(sd) < thb){
                                out[yx] = d*f/abs(bl - x); rp++;
                            }
                        }

                    }*/
                } else {
                    //printf("Don't have any pixeles\n");
                    //out[yx] = 255;
                }
                //printf("rx = %d lx = %d sad = %d nv = %d dir = %d\n", x, bl, sadt, nvt, dirt);

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


    printf("Total pixels = %d good pixele = %d  refiment  = %d %d \n", tp, gp, rp, gp*100/tp);

    tp = 0;
    for(i=0; i < size; i++) if(out[i]) tp++;
    printf("stereo_maching7: Cloud point number = %d\n", tp);
}


/** \brief Filtering outlier pixels.
    \param in	The input 16 bits disparity image.
    \param out	The output 16 bits filtered disparity image.
    \param buff	The temporary buffer.
    \param w	The image width.
    \param h	The imahe height.
*/
void stereo_filter(const int16 *in, int16 *out, int16 *buff, const int w, const int h)
{
    int i, j, y, x, y1, x1, yx, yw, yw1;
    int sum, rp = 0, nz;

    int sh = 3, ls = (sh<<1)+1, h1 = h-sh, w2 = w + (sh<<1), w1 = w*sh, sq = ls*ls;
    int16 *d[ls], *tm;

    //For disparity image
    d[0] = buff;
    for(i=1; i < ls; i++) d[i] = &d[i-1][w2];
    //Prepare buffer
    for(i=0; i < ls-1; i++) cp_line_16(&in[w*abs(sh-i)], d[i], w, sh);

    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y < h1 ? yw + w1 : w*(((h-1)<<1)-y);
        cp_line_16(&in[yw1], d[ls-1], w, sh);

        for(x=0; x < w; x++){
            yx = yw + x;
            out[yx] = in[yx];

            if(d[sh][x+sh]){
                sum = 0; nz = 0;
                for(y1=0; y1 < ls; y1++){
                    for(x1=x; x1 < ls+x; x1++){
                        if(d[y1][x1]){
                            sum += abs(d[y1][x1] - d[sh][x+sh]);
                            nz++;
                        }
                        //if(!i) printf(" %3d %3d  ", l[y1][yxl], r[y1][yxr]);
                    }
                    //if(!i) printf("\n");
                }

                sum = sum/nz;
                if(sum > 5) {
                    //printf("d = %d sum = %d\n", d[sh][x+sh], sum);
                    out[yx] = 0;
                    rp++;
                }
            }

        }
        tm = d[0];  for(i=0; i < ls-1; i++) d[i]  = d[i+1];  d[ls-1]  = tm;
    }
    printf("Remove pixels = %d\n", rp);
}


/** \brief Calculate disparity
    \param in	The input 16 bits first image.
    \param in1 	The input 16 bits second image.
    \param out	The output 16 bits disparity.
    \param buff	The temporary buffer.
    \param w	The image width.
    \param h	The imahe height.
*/
void stereo_disparity(const int16 *left, const int16 *right, int16 *out, int16 *buff, const int w, const int h)
{
    int i, size = w*h, th = 6;
    int16 *Y[2][2], *buf;
    int16 st[w], en[w];

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

    get_left_array(st, en, 5000., 500., 4., 80., 3.54, 3.54/1312.);

    //stereo_maching(Y[0][0], Y[1][0], Y[0][1], Y[1][1], out, buf, w, h);

    //stereo_filter(out, out, buf, w, h);

}

