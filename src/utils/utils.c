#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libit/types.h"
#include "./utils.h"

/** \brief Chage image bytes order.
    \param in       The input 16 bits image.
    \param w		The image width.
    \param h		The image height.
*/
void utils_cnange_bytes(int16 *in, const int w, const int h)
{
    uint32 i, size = w*h;
    for(i=0; i < size; i++) in[i] = ((in[i]&0x00FF)<<8) | ((in[i]&0xFF00)>>8);
}

/** \brief Get 16 bits grey image statistics.
    \param in       The input 16 bits image.
    \param w		The image width.
    \param h		The image height.
    \param bpp      The bits per pixel.
    \param min      The minimum image value.
    \param max      The maximum image value.
*/
void utils_get_stat(int16 *in, const int w, const int h, int *bpp, int *min, int *max)
{
    uint32 i, size = w*h;
    *max = in[0]; *min = in[0];
    for(i=1; i < size; i++){
        if      (in[i] > *max) *max = in[i];
        else if (in[i] < *min) *min = in[i];
    }
    for(*bpp=1; (*max)>>(*bpp); (*bpp)++);
}

/** \brief Make exponential lut table.
    \param ex       The input array.
    \param sd		The standard deviation
    \param sz		The array size
*/
void utils_lut_exp(int *ex, const int sd, const int sz)
{
    int i;
    double hg = sd*sd;
    for(i=0; i < sz; i++){
        ex[i] = (int)(exp(-(double)i*i/(double)hg)*sz);
        //printf("%3d exp = %d\n", i, ex[i]);
    }
}


/**	\brief Transform 16 bits image to 8 bits image.
    \param in	 		The input 16 bits rgb image.
    \param out	 		The output 8 bits rgb image.
    \param w 			The image width.
    \param h 			The image height.
    \param bpp          The bits per pixel.
    \param par          If  0 - direct transform grey to grey image,
                            1 - scale transform grey to grey image,
                            2 - direct transform grb to rgb image,
                            3 - scale transform grb to rgb image,
    \retval             The output 8 bits image.
*/
uint8* utils_16_to_8(const int16 *in, uint8 *out, const int w, const int h,  int bpp, const int par)
{
    int i, j, size, sh, tmp;

    if(bpp < 8) bpp = 8;

    if(par == 0) {
        size = w*h; sh = 0;
    } else if (par == 1) {
        size = w*h; sh = bpp - 8;
    } else if (par == 2) {
        size = w*h*3; sh = 0;
    } else if (par == 3) {
        size = w*h*3; sh = bpp - 8;
    }

    //printf("w = %d h = %d size = %d\n", w, h, size);
    for(i = 0; i < size; i++) {
        tmp = in[i] >> sh;
        out[i] = tmp > 255 ? 255 : tmp;
    }

    return out;
}

/**	\brief Zoom out the rgb 16 bits image in integer times.
    \param in	 		The input rgb image.
    \param out	 		The output rgb image.
    \param buff	 		The temporary buffer, should include 1 row of rgb image size = sizeof(uint32)*w*3).
    \param w 			The image width.
    \param h 			The image height.
    \param zoom 		The zoom parameter 1 - no zoom, 2 - twice, 3 - three times ...
*/
void utils_zoom_out_rgb16_to_rgb16(const int16 *in, int16 *out, uint32 *buff, const int w, const int h, const int zoom)
{
    int i, j, x, x1, y, y1, yw, yx, yxi, sq = zoom*zoom, w1 = w/zoom;
    uint32 max = 1<<31, sh = 0;
    uint32 *buff1, *buff2;

    buff1 = &buff[w1]; buff2 = &buff1[w1];

    memset(buff, 0, sizeof(uint32)*w1*3);

    //Find zoom value when / can changed to >>
    for(i=2; i < max; i<<=1) if((i|zoom) == i) sh++;

    for(y=0, y1=0; y < h; y+=zoom, y1++){

        for(j=0; j < zoom; j++){
            yw = (y+j)*w;
            for(x=0, x1=0; x < w; x+=zoom, x1++){
                yx = yw + x;
                for(i=0; i < zoom; i++) {
                    yxi = (yx+i)*3;
                    buff [x1] += in[yxi  ];
                    buff1[x1] += in[yxi+1];
                    buff2[x1] += in[yxi+2];
                }
                if(j == zoom-1) {
                    yxi = (y1*w1+x1)*3;
                    if(sh) {
                        out[yxi  ] =  buff [x1]>>zoom;
                        out[yxi+1] =  buff1[x1]>>zoom;
                        out[yxi+2] =  buff2[x1]>>zoom;
                    } else {
                        out[yxi  ] =  buff [x1]/sq;
                        out[yxi+1] =  buff1[x1]/sq;
                        out[yxi+2] =  buff2[x1]/sq;
                    }
                    buff[x1] = buff1[x1] = buff2[x1] = 0;
                }
            }
        }
    }
}

/**	\brief Zoom out of bayer image and convert to rgb format.
    \param in	 	The input 16 bits bayer image.
    \param out	 	The output 16 bits rgb image.
    \param buff	 	The temporary buffer, should include 2 row of bayer image size = sizeof(uint32)*w*2.
    \param w 		The image width.
    \param h 		The image height.
    \param zoom 	The zoom parameter 1 - 2x, 2 - 4x, 3 - 6x times ...
    \param bay		The Bayer grids pattern.
*/
void utils_zoom_out_bayer16_to_rgb16(const uint16 *in, uint16 *out, uint32 *buff, const int w, const int h, const int zoom, const BayerGrid bay)
{
    int i, j, x, x1, y, y1, yw, yx, sq = zoom*zoom, zoom2 = zoom<<1, w1 = w/zoom2, zm;
    uint32 max = 1<<31, sh = 0;
    uint32 *c[4];   //Three color buffer

    c[0] = buff; c[1] = &c[0][w1]; c[2] = &c[1][w1]; c[3] = &c[2][w1];

    memset(c[0], 0, sizeof(uint32)*w1<<2);

    //Find zoom value when / can changed to >>
    for(i=1; i < max; i<<=1) if((i|zoom) == i) sh++;

    zm = (zoom == 1) ? 0 : zoom-2;
    //printf("zoom = %d sh = %d\n", zoom, sh);

    for(y=0, y1=0; y < h; y+=zoom2, y1++){

        for(j=0; j < zoom2; j+=2){
            yw = (y+j)*w;
            for(x=0, x1=0; x < w; x+=zoom2, x1++){
                yx = yw + x;
                for(i=0; i < zoom2; i+=2) {
                    //printf("%d %d %d %d\n", in[yx+i], in[yx+i+1], in[yx+i+w], in[yx+i+w+1]);
                    c[0][x1] += in[yx+i];
                    c[1][x1] += in[yx+i+1];
                    c[2][x1] += in[yx+i+w];
                    c[3][x1] += in[yx+i+w+1];
                }
                if(j == zoom2-2) {
                    switch(bay){
                    case(BGGR):{
                        out[(y1*w1+x1)*3]   = sh ? c[3][x1]>>zm : c[3][x1]/sq;
                        out[(y1*w1+x1)*3+1] = sh ? (c[1][x1]+c[2][x1])>>(zm+1) : (c[1][x1]+c[2][x1])/(sq<<1);
                        out[(y1*w1+x1)*3+2] = sh ? c[0][x1]>>zm : c[0][x1]/sq;
                        break;
                    }
                    case(GRBG):{
                        out[(y1*w1+x1)*3]   = sh ? c[1][x1]>>zm : c[01][x1]/sq;
                        out[(y1*w1+x1)*3+1] = sh ? (c[0][x1]+c[3][x1])>>(zm+1) : (c[0][x1]+c[3][x1])/(sq<<1);
                        out[(y1*w1+x1)*3+2] = sh ? c[2][x1]>>zm : c[2][x1]/sq;
                        break;
                    }
                    case(GBRG):{
                        out[(y1*w1+x1)*3]   = sh ? c[2][x1]>>zm : c[2][x1]/sq;
                        out[(y1*w1+x1)*3+1] = sh ? (c[0][x1]+c[3][x1])>>(zm+1) : (c[0][x1]+c[3][x1])/(sq<<1);
                        out[(y1*w1+x1)*3+2] = sh ? c[1][x1]>>zm : c[1][x1]/sq;
                        break;
                    }
                    case(RGGB):{
                        out[(y1*w1+x1)*3]   = sh ? c[0][x1]>>zm : c[0][x1]/sq;
                        out[(y1*w1+x1)*3+1] = sh ? (c[1][x1]+c[2][x1])>>(zm+1) : (c[1][x1]+c[2][x1])/(sq<<1);
                        out[(y1*w1+x1)*3+2] = sh ? c[3][x1]>>zm : c[3][x1]/sq;
                        break;
                    }
                    }
                    c[0][x1] = c[1][x1] = c[2][x1] = c[3][x1] = 0;
                }
            }
        }
    }
}

/**	\brief Fill histogram from bayer image.
    \param in	 	The input 16 bits bayer image.
    \param R	 	The output red histogram.
    \param G	 	The output green histogram.
    \param B	 	The output blue histogram.
    \param Y	 	The output Y histogram.
    \param buff	 	The half line buffer for calculation Y histogram.
    \param w 		The image width.
    \param h 		The image height.
    \param bg		The Bayer grids pattern.
    \param bpp		The pits per pixel.
*/
void utils_fill_hist_bayer(const uint16 *in, int *R, int *G, int *B, int *Y, int *buff, const int w, const int h, const BayerGrid bg, const int bpp)
{
    int i, x, x1, y, yw, yx;
    int *l, ns = 1<<bpp;

    //Clear histogram
    for(i=0; i < ns; i++) R[i] = 0;
    for(i=0; i < ns; i++) G[i] = 0;
    for(i=0; i < ns; i++) B[i] = 0;
    for(i=0; i < ns; i++) Y[i] = 0;

    //Rows buffer for Y
    l = buff;

    for(y = 0; y < h; y++){
        yw = y*w;
        for(x = 0; x < w; x++){
            x1 = x>>1;
            yx = yw + x;
            switch(bg){
            case(BGGR):{
                if(y&1){
                    if(x&1) { R[in[yx]]++; Y[(612*in[yx] + l[x1])>>11]++; }
                    else    { G[in[yx]]++; l[x1] += 601*in[yx]; }
                } else {
                    if(x&1) { G[in[yx]]++; l[x1] += 601*in[yx]; }
                    else    { B[in[yx]]++; l[x1]  = 234*in[yx]; }
                }
                break;
            }
            case(GRBG):{
                if(y&1){
                    if(x&1) { G[in[yx]]++; Y[(601*in[yx] + l[x1])>>11]++; }
                    else    { B[in[yx]]++; l[x1] += 234*in[yx]; }

                } else {
                    if(x&1) { R[in[yx]]++; l[x1] += 612*in[yx]; }
                    else    { G[in[yx]]++; l[x1]  = 601*in[yx]; }
                }
                break;
            }
            case(GBRG):{
                if(y&1){
                    if(x&1) { G[in[yx]]++; Y[(601*in[yx] + l[x1])>>11]++; }
                    else    { R[in[yx]]++; l[x1] += 612*in[yx]; }

                } else {
                    if(x&1) { B[in[yx]]++; l[x1] += 234*in[yx]; }
                    else    { G[in[yx]]++; l[x1]  = 601*in[yx]; }
                }
                break;
            }
            case(RGGB):{
                if(y&1){
                    if(x&1) { B[in[yx]]++; Y[(234*in[yx] + l[x1])>>11]++; }
                    else    { G[in[yx]]++; l[x1] += 601*in[yx]; }
                } else {
                    if(x&1) { G[in[yx]]++; l[x1] += 601*in[yx]; }
                    else    { R[in[yx]]++; l[x1]  = 612*in[yx]; }
                }
                break;
            }
            //printf("Start in = %d out = %d blm = %d\n", in[yx], out[yx], blm);
            }
        }
    }
}

/**	\brief Fill histogram for white balancing.
    \param in	 	The input 16 bits bayer image.
    \param hi	 	The output bayer histogram.
    \param r	 	The avarage red.
    \param g	 	The avarage green.
    \param b	 	The avarage blue.
    \param buff	 	The half line buffer for calculation Y histogram.
    \param w 		The image width.
    \param h 		The image height.
    \param bg		The Bayer grids pattern.
    \param bpp		The pits per pixel.
*/
void utils_hist_wb(const uint16 *in, int *hi, int *r, int *g, int *b, const int w, const int h, const BayerGrid bg, const int bpp)
{
    int i, x,  y, yw, yx;
    int ns = 1<<bpp, sz = w*h;

    //Clear histogram
    for(i=0; i < ns; i++) hi[i] = 0;
    (*r) = 0; (*g) = 0; (*b) = 0;


    for(y = 0; y < h; y++){
        yw = y*w;
        for(x = 0; x < w; x++){
            yx = yw + x;
            hi[in[yx]]++;
            switch(bg){
            case(BGGR):{
                if(y&1){
                    if(x&1) (*r)+=in[yx];
                    else    (*g)+=in[yx];
                } else {
                    if(x&1) (*g)+=in[yx];
                    else    (*b)+=in[yx];
                }
                break;
            }
            case(GRBG):{
                if(y&1){
                    if(x&1) (*g)+=in[yx];
                    else    (*b)+=in[yx];
                } else {
                    if(x&1) (*r)+=in[yx];
                    else    (*g)+=in[yx];
                }
                break;
            }
            case(GBRG):{
                if(y&1){
                    if(x&1) (*g)+=in[yx];
                    else    (*r)+=in[yx];
                } else {
                    if(x&1) (*b)+=in[yx];
                    else    (*g)+=in[yx];
                }
                break;
            }
            case(RGGB):{
                if(y&1){
                    if(x&1) (*b)+=in[yx];
                    else    (*g)+=in[yx];
                } else {
                    if(x&1) (*g)+=in[yx];
                    else    (*r)+=in[yx];
                }
                break;
            }
            //printf("Start in = %d out = %d blm = %d\n", in[yx], out[yx], blm);
            }
        }
    }
    *r = (*r)/(sz>>2); *g = (*g)/(sz>>1); *b = (*b)/(sz>>2);
    //printf("Start r = %d g = %d b = %d\n", *r, *g, *b);
}

/**	\brief Calculate the white balance multiplier for red and blue color of 16 bits rgb image.
    \param in	The input 16 bits rgb image.
    \param rm   The pointer to the red multiplier.
    \param bm   The pointer to the blue multiplier.
    \param mp   The value of pedestal.
    \param mx   The value of for removing dark pixeles.
    \param buff The histogram buffer size = sizeof(uint32)*(1<<bpp).
    \param w    The image width.
    \param h 	The image height.
    \param sh   The shift for integer multiplier.
    \param bpp  The input image bits per pixel.
*/
void utils_wb(int16 *in, int *rm, int *bm, int mp, int mx,  int r, int g, int b, const int w, const int h, const int sh)
{
    int i, j, size = w*h, size3 = size*3;
    uint32 d, d1;
    //float s = 0.01,  th = 0.5;
    int m, mr,  mb, s = 10;
    /*
    uint32 *hi = buff;
    memset(hi, 0, sizeof(uint32)*hs);
    //Gray world algorithm the first step of iteration
    //Make Y histogram
    for(i = 0; i < size3; i+=3) {
        r += in[i  ];
        g += in[i+1];
        b += in[i+2];
        //Y = (306*in[i] + 601*in[i+1] + 117*in[i+2])>>10;
        //if(!Y) printf("r = %d g = %d b = %d Y = %d\n", in[i  ], in[i+1], in[i+2], Y);
        //hi[Y]++;
        hi[in[i  ]]++; hi[in[i+1]]++; hi[in[i+2]]++;
    }
    //for(i=0; i< 500; i++) printf("%d hi = %d\n", i, hi[i]);


    //Find pedestal
    //printf("tp = %d hs = %d\n", tp, hs);
    sum = 0;
    for(i=0; i < hs; i++) {
        //printf("%d sum = %d hi = %d\n", i, sum, hi[i]);
        sum += hi[i];
        if(sum > tp) break;
    }
    *mp = i;

    //Find threshold for removing dark pixels.
    sum = 0;
    for(i=0; i < hs; i++) {
        sum += hi[i];
        if(sum > ts) break;
    }
    mx = i;
    */
    //Remove pedestal and dark pixels
    for(i = 0; i < size3; i+=3) {
        //Y = (306*in[i] + 601*in[i+1] + 117*in[i+2])>>10;
        //if(Y < mx) {
        if(in[i] < mx || in[i+1] < mx || in[i+2] < mx) {
            in[i] = 0; in[i+1] = 0; in[i+2] = 0;
        } else {
            in[i  ] = (in[i  ] - mp) < 0 ? 0 : (in[i  ] - mp);
            in[i+1] = (in[i+1] - mp) < 0 ? 0 : (in[i+1] - mp);
            in[i+2] = (in[i+2] - mp) < 0 ? 0 : (in[i+2] - mp);
        }
    }

    //Gray world multiplier for first step approximation
    //printf("size = %d r = %d g = %d b = %d \n", size, r, g, b);
    //printf("r = %d g = %d b = %d \n", r, g, b);
    r = r-mp; g = g-mp; b = b-mp;
    //printf("r = %d g = %d b = %d \n", r, g, b);
    //mr = (double)g/(double)r;
    //mb = (double)g/(double)b;
    mr = (g<<sh)/r;
    mb = (g<<sh)/b;
    //printf("mr = %d mb = %d mr = %f mb = %f\n",mr, mb,(double)mr/(double)(1<<sh), (double)mb/(double)(1<<sh));
    //printf("cn = %d sz = %d p = %f\n", cn, size, (double)(size - cn)/(double)size);

    //Red color (g - r) minimization
    d = 0; m = mr;
    for(i = 0; i < size3; i+=3) d += abs(in[i+1] - (in[i]*m>>sh));
    for(j=0; ;j++){
        m = m + s;
        d1 = 0;
        for(i = 0; i < size3; i+=3) d1 += abs(in[i+1] - (in[i]*m>>sh));
        //printf("j = %d d = %d d1 = %d m = %f\n", j, d, d1, (double)m/(double)(1<<sh));
        if(!j && d1 > d) s = -s;
        if( j && d1 > d) break;
        d = d1;
    }
    *rm = m;

    //Blue color (g - b) minimization
    d = 0; m = mb;
    for(i = 0; i < size3; i+=3) d += abs(in[i+1] - (in[i+2]*m>>sh));
    for(j=0; ;j++){
        m = m + s;
        d1 = 0;
        for(i = 0; i < size3; i+=3) d1 += abs(in[i+1] - (in[i+2]*m>>sh));
        //printf("j = %d d = %d d1 = %d m = %f\n", j, d, d1, (double)m/(double)(1<<sh));
        if(!j && d1 > d) s = -s;
        if( j && d1 > d) break;
        d = d1;
    }
    *bm = m;
}

/**	\brief White balance 16 bits rgb image.
    \param in	The input 16 bits rgb image.
    \param out	The output 16 bits rgb image.
    \param buff	The the buffer size should be 1000000 bytes.
    \param w    The image width.
    \param h 	The image height.
    \param bpp  The image bits per pixel.
*/
void utils_wb_rgb(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bpp)
{
    int i, j, size3 = h*w*3, sh = 10, z, zoom, w1, h1, max = (1<<bpp)-1;
    int rm, bm, mp, mx, r, g, b;

    //Check zoom image size
    for(i=0; (w*h>>i) > 500000; i+=2);
    z = i>>1; zoom = 1<<z; w1 = w>>z; h1 = h>>z;

    utils_zoom_out_rgb16_to_rgb16(in, buff, (uint32*)&buff[w1*h1*3], w, h, zoom);

    utils_wb(buff, &rm, &bm, mp, mx, r, g, b, w1, h1, sh);
    printf("rm = %d bm = %d rm = %f bm = %f \n", rm, bm, (double)rm/(double)(1<<sh), (double)bm/(double)(1<<sh));

    for(i = 0; i < size3; i+=3) {
        out[i]   = (in[i]-mp)*rm>>sh;    out[i] = out[i] > max ? max : out[i]; out[i] = out[i] < 0 ? 0 : out[i];
        out[i+1] = in[i+1]-mp; out[i+1] = out[i+1] < 0 ? 0 : out[i+1];
        out[i+2] = (in[i+2]-mp)*bm>>sh;  out[i+2] = out[i+2] > max ? max : out[i+2]; out[i+2] = out[i+2] < 0 ? 0 : out[i+2];
    }
}

/**	\brief White balance 16 bits bayer image.
    \param in	The input 16 bits bayer image.
    \param out	The output 16 bits rgb image.
    \param buff	The the buffer size should be 2000000 bytes.
    \param w    The image width.
    \param h 	The image height.
    \param bpp  The image bits per pixel.
    \param bg   The Bayer grid pattern
*/
void utils_wb_bayer(const int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bpp, const int bg)
{
    int sh = 10, z, zoom, w1, h1, max = (1<<bpp)-1, sz = w*h;
    int i, j, x, y, yx, yw, rm, bm, mp, mx, sum, ts = sz>>5, tp = sz>>10;
    int *hi, hs = 1<<bpp;
    int r, g, b;
    hi = (int*)buff;

    //Check zoom image size
    for(i=0; (w*h>>i) > 1000000; i+=2);
    z = i>>1; zoom = 1<<z; w1 = w>>(z+1); h1 = h>>(z+1);
    //printf("z = %d zoom = %d  size = %d newsize = %d\n", z, zoom, w*h, w1*h1);

    utils_hist_wb(in, hi, &r, &g, &b, w, h, bg, bpp);

    //Calculate pedestal

    sum = 0; j = 0;
    for(i=0; i < hs; i++) {
        //printf("%d sum = %d hi = %d\n", i, sum, hi[i]);
        sum += hi[i];
        if(sum > tp && !j) j = i;
        if(sum > ts) break;
    }
    //Pedestal
    mp = j;
    //Threshold for removing dark pixels.
    mx = i;
    //printf("Pedestal = %d dark pixels = %d\n", mp, mx);

    //printf("zoom = %d , w1 = %d h1 = %d\n", zoom, w1, h1);
    utils_zoom_out_bayer16_to_rgb16(in, buff, (uint32*)&buff[w1*h1*3], w, h, zoom, bg);
    //printf("zoom = %d , w1 = %d h1 = %d\n", zoom, w1, h1);

    utils_wb(buff, &rm, &bm, mp, mx, r, g, b, w1, h1, sh);
    printf("rm = %d bm = %d rm = %f bm = %f mp = %d\n", rm, bm, (double)rm/(double)(1<<sh), (double)bm/(double)(1<<sh), mp);
    //rm = 5000;
    switch(bg){
    case(BGGR):{
        for(y=0; y < h; y++){
            yw = y*w;
            for(x=0; x < w; x++){
                yx = yw + x;
                if(y&1){
                    if(x&1) { out[yx] = (in[yx]-mp)*rm>>sh;  out[yx] = out[yx] > max ? max : out[yx]; out[yx] = out[yx] < 0 ? 0 : out[yx];}
                    else    { out[yx] = in[yx]-mp; out[yx] = out[yx] < 0 ? 0 : out[yx]; }
                } else {
                    if(x&1) { out[yx] = in[yx]-mp; out[yx] = out[yx] < 0 ? 0 : out[yx]; }
                    else    { out[yx] = (in[yx]-mp)*bm>>sh;  out[yx] = out[yx] > max ? max : out[yx]; out[yx] = out[yx] < 0 ? 0 : out[yx];}
                }
            }
        }
        break;
    }
    case(GRBG):{
        for(y=0; y < h; y++){
            yw = y*w;
            for(x=0; x < w; x++){
                yx = yw + x;
                if(y&1){
                    if(x&1) { out[yx] = in[yx]-mp; out[yx] = out[yx] < 0 ? 0 : out[yx]; }
                    else    { out[yx] = (in[yx]-mp)*bm>>sh;  out[yx] = out[yx] > max ? max : out[yx]; out[yx] = out[yx] < 0 ? 0 : out[yx];}
                } else {
                    if(x&1) { out[yx] = (in[yx]-mp)*rm>>sh;  out[yx] = out[yx] > max ? max : out[yx]; out[yx] = out[yx] < 0 ? 0 : out[yx];}
                    else    { out[yx] = in[yx]-mp; out[yx] = out[yx] < 0 ? 0 : out[yx]; }
                }
            }
        }
        break;
    }
    case(GBRG):{
        for(y=0; y < h; y++){
            yw = y*w;
            for(x=0; x < w; x++){
                yx = yw + x;
                if(y&1){
                    if(x&1) { out[yx] = in[yx]-mp; out[yx] = out[yx] < 0 ? 0 : out[yx];}
                    else    { out[yx] = (in[yx]-mp)*rm>>sh;  out[yx] = out[yx] > max ? max : out[yx]; out[yx] = out[yx] < 0 ? 0 : out[yx];}
                } else {
                    if(x&1) { out[yx] = (in[yx]-mp)*bm>>sh;  out[yx] = out[yx] > max ? max : out[yx]; out[yx] = out[yx] < 0 ? 0 : out[yx];}
                    else    { out[yx] = in[yx]-mp; out[yx] = out[yx] < 0 ? 0 : out[yx]; }
                }
            }
        }
        break;
    }
    case(RGGB):{
        for(y=0; y < h; y++){
            yw = y*w;
            for(x=0; x < w; x++){
                yx = yw + x;
                if(y&1){
                    if(x&1) { out[yx] = (in[yx]-mp)*bm>>sh;  out[yx] = out[yx] > max ? max : out[yx]; out[yx] = out[yx] < 0 ? 0 : out[yx];}
                    else    { out[yx] = in[yx]-mp; out[yx] = out[yx] < 0 ? 0 : out[yx]; }
                } else {
                    if(x&1) { out[yx] = in[yx]-mp; out[yx] = out[yx] < 0 ? 0 : out[yx]; }
                    else    { out[yx] = (in[yx]-mp)*rm>>sh;  out[yx] = out[yx] > max ? max : out[yx]; out[yx] = out[yx] < 0 ? 0 : out[yx];}
                }
            }
        }
        break;
    }
    }
}

void inline static cp_line(int16 *in, uint32 *l, uint32 w, uint32 br)
{
    uint32 i;
    for(i=0; i < br; i++) l[i] = in[br-i];
    for(i=0; i < w ; i++) l[i+br] = in[i];
    for(i=0; i < br; i++) l[i+br+w] = in[w-br-i];
}

/** \brief Make the integral matrix of 16 bits grey image with border.
    \param in	The pointer to a input image.
    \param ing 	The pointer to a integral matrix, the size = (w + 2*br)*(h + 2*br)*(sizeof(uint32)).
    \param buff	The two lines buffer, size should be (w + 2*br)*2*(sizeof(int16)) bytes.
    \param w	The image width.
    \param h	The imahe height.
    \param br	The border on the image ing, new w = w + 2*br and new h = h + 2*br;
*/
void utils_integral(int16 *in, uint32 *ing, uint32 *buff, const int w, const int h, const int br)
{
    uint32 x, y=0, yw, yx, w1 = w + (br<<1), h1 = h + (br<<1), h2 = h + br - 1;
    uint32 *l[2], *tm; //Two lines buffer

    l[0] = buff; l[1] = &buff[w1];

    //First line
    cp_line(&in[w*br], l[0], w, br);
    ing[0] = l[0][0];
    for(x=1; x < w1; x++) {
        l[0][x] = l[0][x-1] + l[0][x];
        ing[x] = l[0][x];
    }

    for(y=1; y < h1; y++){
        yw = y*w1;

        if(y < br)      cp_line(&in[w*(br-y)], l[1], w, br);
        else if(y > h2) cp_line(&in[w*(((h-1)<<1)+br-y)], l[1], w, br);
        else            cp_line(&in[w*(y-br)], l[1], w, br);

        l[1][0] = l[0][0] + l[1][0];
        ing[yw] = l[1][0];

        for(x=1; x < w1; x++){
            yx = yw + x;
            l[1][x] = l[1][x-1] + l[0][x] - l[0][x-1] + l[1][x];
            ing[yx] = l[1][x];
        }

        tm = l[1]; l[1] = l[0]; l[0] = tm;
    }
    /*
    uint32 sum = 0;
    for(x=0; x < w1*h1; x++){
        sum += in[x];
    }
    printf("Check the integral matrix = %d sum = %d test = %d\n", ing[w1*h1-1], sum, sum - ing[w1*h1-1]);
    */
}

/** \brief Make the integral matrix of 16 bits bayer image with border.
    \param in	The pointer to a input image.
    \param ing 	The pointer to a integral matrix, the size = (w + 2*br)*(h + 2*br)*(sizeof(uint32)).
    \param buff	The three lines buffer, size should be (w + 2*br)*5*(sizeof(int16)) bytes.
    \param w	The image width.
    \param h	The imahe height.
    \param br	The border on the image ing, br should be more then 0;
*/
void utils_integral_bayer(int16 *in, uint32 *ing, uint32 *buff, const int w, const int h, const int br)
{
    uint32 x, y=0, yw, yx, br2 = br, w1 = w + (br2<<1), h1 = h + (br2<<1), h2 = h + br2 - 1;
    uint32 *l[3], *tm; //Three lines buffer

    l[0] = buff; l[1] = &l[0][w1]; l[2] = &l[1][w1];

    //First two lines
    cp_line(&in[w*br2], l[0], w, br2);
    cp_line(&in[w*(br2-1)], l[1], w, br2);

    ing[0] = l[0][0];
    ing[1] = l[0][1];
    for(x=2; x < w1; x++) {
        l[0][x] = l[0][x-2] + l[0][x];
        ing[x] = l[0][x];
    }

    ing[w1] = l[1][0];
    ing[w1+1] = l[1][1];
    for(x=2; x < w1; x++) {
        l[1][x] = l[1][x-2] + l[1][x];
        ing[w1+x] = l[1][x];
    }

    for(y=2; y < h1; y++){
        yw = y*w1;

        if(y < br2)     cp_line(&in[w*(br2-y)], l[2], w, br2);
        else if(y > h2) cp_line(&in[w*(((h-1)<<1)+br2-y)], l[2], w, br2);
        else            cp_line(&in[w*(y-br2)], l[2], w, br2);

        l[2][0] = l[0][0] + l[2][0];
        ing[yw] = l[2][0];
        l[2][1] = l[0][1] + l[2][1];
        ing[yw+1] = l[2][1];

        for(x=2; x < w1; x++){
            yx = yw + x;
            l[2][x] = l[2][x-2] + l[0][x] - l[0][x-2] + l[2][x];
            ing[yx] = l[2][x];
        }

        tm = l[0]; l[0] = l[1]; l[1] = l[2]; l[2] = tm;
    }
    /*
    //For check only
    uint32 sum;
    sum = 0;
    for(y=0; y < h1; y+=2){ yw = y*w1; for(x=0; x < w1; x+=2) {sum += in[yw + x];}}
    printf("Check the integral sum = %d\n", sum);
    sum = 0;
    for(y=0; y < h1; y+=2){ yw = y*w1; for(x=1; x < w1; x+=2) {sum += in[yw + x];}}
    printf("Check the integral sum = %d\n", sum);
    sum = 0;
    for(y=1; y < h1; y+=2){ yw = y*w1; for(x=0; x < w1; x+=2) {sum += in[yw + x];}}
    printf("Check the integral sum = %d\n", sum);
    sum = 0;
    for(y=1; y < h1; y+=2){ yw = y*w1; for(x=1; x < w1; x+=2) {sum += in[yw + x];}}
    printf("Check the integral sum = %d\n", sum);
    printf("mat1 = %d mat2 = %d mat3 = %d mat4 = %d\n", ing[w1*h1-w1-2], ing[w1*h1-w1-1], ing[w1*h1-2], ing[w1*h1-1]);
    */
}

/** \brief Average each pixel in 16 bits image witn br border out of pixel
    \param in	The input 16 bits image.
    \param out 	The output 16 bits image.
    \param buff	The two lines buffer, size should be (w + 2*br)*2*(sizeof(int16)) bytes.
    \param w	The image width.
    \param h	The imahe height.
    \param br	The border of averaging if br = 1  3x3 = 9 pixel.
*/
void utils_average(int16 *in, int16 *out, uint32 *buff, const int w, const int h, const int br)
{
    int i,x, y, yw, yw1, yx, yx1;
    int w1 = w + ((br+1)<<1), h1 = h + ((br+1)<<1), bs = ((br<<1)+1)*((br<<1)+1), br2 = (br<<1) + 1;
    uint32 *ing, sh;
    int64_t b;

    //Finding b coefficient to change / to *
    for(i=1; bs>>i; i++);
    sh = 63 - i - 16;
    b = (1LL<<sh)/bs;
    //printf("i = %d b = %lld bs = %d sh = %d\n", i, b, bs, sh);

    ing = (uint32*)malloc(w1*h1*sizeof(uint32));
    if (ing == NULL) {
        fprintf(stderr, "Error! utils_average: Can't allocate memory\n");
        return;
    }

    utils_integral(in, ing, buff, w, h, br+1);

    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y*w1;
        for(x=0; x < w; x++){
            yx = yw + x;
            yx1 = yw1 + x;
            out[yx] = (ing[yx1 + br2 + br2*w1] + ing[yx1] - ing[yx1 + br2] - ing[yx1 + br2*w1])*b>>sh;
            //out[yx] = (ing[yx1 + br2 + br2*w1] + ing[yx1] - ing[yx1 + br2] - ing[yx1 + br2*w1])/bs;
            //yx1 = yw1 + x + br + 1;
            //out[yx] = (ing[yx1 + br + br*w1] + ing[yx1 - (br+1) - (br+1)*w1] - ing[yx1 + br - (br+1)*w1] - ing[yx1 + br*w1 - (br+1)])/bs;
            //yx1 = yw1 + x + br2;
            //out[yx] = (ing[yx1 + br + br*w1] + ing[yx1 - br2 - br2*w1] - ing[yx1 + br - br2*w1] - ing[yx1 + br*w1 - br2])/bs;
        }
    }

    if(ing) free(ing);
}

/** \brief Average each pixel in 16 bits bayer image witn br border out of pixel
    \param in	The input 16 bits image.
    \param out 	The output 16 bits image.
    \param buff	The two lines buffer, size should be (w + 2*br)*2*(sizeof(int16)) bytes.
    \param w	The image width.
    \param h	The imahe height.
    \param br	The border of averaging if br = 1  3x3 = 9 pixel.
*/
void utils_average_bayer(int16 *in, int16 *out, uint32 *buff, const int w, const int h, const int br)
{
    int i,x, y, yw, yw1, yx, yx1, br2 = (br+2);
    int w1 = w + (br2<<1), h1 = h + (br2<<1), bs = (br+1)*(br+1), br4 = (br<<1) + 2;
    uint32 *ing, sh;
    int64_t b;

    //Finding b coefficient to change / to *
    for(i=1; bs>>i; i++);
    sh = 63 - i - 16;
    b = (1LL<<sh)/bs;
    //printf("i = %d b = %lld bs = %d sh = %d\n", i, b, bs, sh);

    ing = (uint32*)malloc(w1*h1*sizeof(uint32));
    if (ing == NULL) {
        fprintf(stderr, "Error! utils_average: Can't allocate memory\n");
        return;
    }

    //printf("start utils_integral_bayer\n");
    utils_integral_bayer(in, ing, buff, w, h, br2);

    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = y*w1;
        for(x=0; x < w; x++){
            yx = yw + x;
            yx1 = yw1 + x;
            out[yx] = (ing[yx1 + br4 + br4*w1] + ing[yx1] - ing[yx1 + br4] - ing[yx1 + br4*w1])*b>>sh;
            //out[yx] = (ing[yx1 + br2 + br2*w1] + ing[yx1] - ing[yx1 + br2] - ing[yx1 + br2*w1])/bs;
            //yx1 = yw1 + x + br + 1;
            //out[yx] = (ing[yx1 + br + br*w1] + ing[yx1 - (br+1) - (br+1)*w1] - ing[yx1 + br - (br+1)*w1] - ing[yx1 + br*w1 - (br+1)])/bs;
            //yx1 = yw1 + x + br2;
            //out[yx] = (ing[yx1 + br + br*w1] + ing[yx1 - br2 - br2*w1] - ing[yx1 + br - br2*w1] - ing[yx1 + br*w1 - br2])/bs;
        }
    }
    if(ing) free(ing);
}

/** \brief Subtraction one grey image from another.
    \param in	The input 16 bits image.
    \param in1 	The input 16 bits image.
    \param out	The output 16 bits difference.
    \param w	The image width.
    \param h	The imahe height.
    \param bpp  The image bits per pixel.
*/
void utils_subtract(const int16 *in, const int16 *in1, int16 *out, const int w, const int h, const int bpp)
{
    int i, j, size = w*h, sh = 1<<(bpp-1), tmp, max = (1<<bpp)-1;
    uint32 sum = 0;

    for(i = 0; i < size; i++) {
        tmp = in[i] - in1[i] + sh;
        out[i] = tmp < 0 ? 0 : (tmp > max ? max : tmp);
        sum += abs(in[i] - in1[i]);
    }
    sum = sum/size;
    printf("diff = %d\n", sum);
}

/** \brief Add two grey images.
    \param in	The input 16 bits image.
    \param in1 	The input 16 bits image.
    \param out	The output 16 bits imafe.
    \param w	The image width.
    \param h	The imahe height.
    \param bpp  The image bits per pixel.
*/
void utils_add(const int16 *in, const int16 *in1, int16 *out, const int w, const int h, const int bpp, const int bpp1)
{
    int i, j, size = w*h, sh = 1<<(bpp-1), tmp, max = (1<<bpp1)-1;

    for(i = 0; i < size; i++) {
        //tmp = in[i] + in1[i] - sh;
        //out[i] = tmp < 0 ? 0 : (tmp > max ? max : tmp);
        out[i] = in1[i] ? in1[i] : in[i];
        //printf("in = %d in1 = %d sh = %d tmp = %d out = %d\n", in[i], in1[i], sh, tmp, out[i]);
    }
}

/** \brief Two times downsampling image.
    \param in	The input 16 bits image.
    \param out 	The output 16 bits image.
    \param buff	One line buffer
    \param w	The image width.
    \param h	The imahe height.
*/
void utils_resize_down_2(const int16 *in, int16 *out, int16 *buff, const int w, const int h)
{
    int x, y, yw, yx, yw1, w1 = w>>1, h1 = h>>1, x2 = x<<1;
    int16 *l = buff;


    for(y=0; y < h1; y++){
        yw1 = (y<<1)*w;
        for(x=0; x < w1; x++) l[x]  = in[yw1 + (x<<1)] + in[yw1 + (x<<1)+1];
        yw1 = yw1 + w;
        for(x=0; x < w1; x++) l[x] += in[yw1 + (x<<1)] + in[yw1 + (x<<1)+1];
        yw = y*w1;

        for(x=0; x < w1; x++){
            yx = yw + x;
            out[yx] = l[x]>>2;
        }
    }
}

