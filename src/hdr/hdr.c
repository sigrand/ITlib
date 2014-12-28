#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libit/types.h"
#include "./hdr.h"
#include "../filters/filters.h"
#include "../utils/utils.h"


/**	\brief Automatic Color Enhancement algorithm.
    \param in       The input 16 bits image.
    \param out      The output 16 bits image.
    \param buff     The temporary buffer.
    \param bpp      The input image bits per pixel.
    \param bpp1     The output image bits per pixel.
    \param w        The image width.
    \param h        The image height.
    \param cs       The type of input image ColorSpace.
*/
void hdr_ace(const int16 *in, int16 *out, int *buff, const int w, const int h, const int bpp, const int bpp1, const int cs)
{
    int x, y, yx, yw, hs = 1<<bpp, size;
    int *hi;
    int64_t b, sum, shb = 40 - bpp1;

    if(cs == BAYER || cs == GREY){
        size = w*h;
    } else if(cs == RGB || cs == YUV444){
        size = w*h*3;
    } else if (cs == YUV420){
        size = w*h*3>>1;
    } else {
        fprintf(stderr, "hdr_ace Error! Can't support image format %d\n", cs);
        return;
    }

    b = (1LL<<40)/size;
    //printf("b = %d w*h = %d size = %d\n", b, w*h, size);

    hi = buff; //hl = &hi[hs]; hr = &hl[hs]; lt = &hr[hs];

    //Fill historgam
    memset(hi, 0, sizeof(int)*hs);
    for(x=0; x < size; x++) hi[in[x]]++;

    //Make lookup table
    sum = 0;
    for(x=0; x < hs; x++) {
        sum += hi[x];
        //if(x<500) printf("%d %d ", x, hi[x]);
        hi[x] = sum*b>>shb;
        //if(x<500) printf("%d \n", hi[x]);
        //hi[x] = (sum<<bpp1)/size;
    }
    //printf("hs = %d sum = %d\n", hs, sum);
    for(x=0; x < size; x++) out[x] = hi[in[x]];
    /*
    for(y=0; y < h; y++){
        yw = y*w;
        for(x=0; x < w; x++){
            yx = yw + x;
            //For sony A55 sensor check
            //if(in[yx] > hs-1) in[yx] = hs-1;
            out[yx] = hi[in[yx]];
            //if(out[yx] < 0 || out[yx] > 255) printf("yx = %d out = %d in = %d lt = %d\n", yx, out[yx], in[yx], lt[in[yx]]);
        }
    }
    */
}


/**	\brief Get the differnce between 16 bits  original and averaged bayer image.
    \param in       The input 16 bits bayer image.
    \param avr      The input 16 bits averaged bayer image.
    \param dif      The output image.
    \param w        The image width.
    \param h        The image height.
    \param bpp      The input images bits per pixel.
*/
void hdr_diff(const int16 *in, const int16 *avr, int16 *dif, const int w, const int h, const int bpp)
{
    int x, y, yw, yx;
    int hs = 2, ws = 2, w2 = w<<1, whs = w*hs, bs;
    int sh = 1<<(bpp-1);

    //bs = (ws+1)*(hs+1); //
    //bs = ((ws<<1)+1)*((hs<<1)+1)>>2;
    //bs = 9;

    //printf("bs = %d\n", bs);

    //utils_integral_bayer(in, buff, w, h);

    for(y=0; y < h; y++){
        yw = y*w;
        for(x=0; x < w; x++){
            yx = yw + x;
            //out[yx] = (buff[yx+ws+whs] + buff[yx-ws-whs-w2-2] - buff[yx+ws-whs-w2] - buff[yx-ws+whs-2])/bs;
            dif[yx] = (in[yx]*15>>4) - (avr[yx]*17>>4) + sh;
            //dif[yx] = lb(dif[yx]);
            //if(dif[yx]< 0) dif[yx] = 0;
            //av[0][yxr] = avr/bs;
            //ws = 2; whs = ws*w; //bs = ((ws<<1)+1)*((ws<<1)+1);
            //av[1][yxr] = ing[yxr+ws+whs] + ing[yxr-ws-whs-w2-2] - ing[yxr+ws-whs-w2] - ing[yxr-ws+whs-2] - av[0][yxr];

            //av[1][yxr] = avr/bs;
            //ws = 6; whs = ws*w; bs = ((ws<<1)+1)*((ws<<1)+1)>>2;
            //avr  = ing[yxr+ws+whs] + ing[yxr-ws-whs-w2-2] - ing[yxr+ws-whs-w2] - ing[yxr-ws+whs-2];
            //av[2][yxr] = avr/bs;
        }
    }
}

/**	\brief Local Automatic Color Enhancement algorithm.
    \param in       The input 16 bits bayer image.
    \param out      The output 16 bits bayer image.
    \param buff     The temporary buffer.
    \param bpp      The input image bits per pixel.
    \param w        The image width.
    \param h        The image height.
*/
void hdr_ace_local(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bpp)
{
    int i, size = w*h;
    int16 *dif = &buff[w*h];

    filters_median_bayer(in, out, buff, w, h, 0);
    hdr_diff(in, out, dif, w, h, bpp);
    hdr_ace(in, out, (int*)buff, w, h, bpp, 8, GREY);
    hdr_ace(dif, in, (int*)buff, w, h, bpp, 6, GREY);

    for(i=0; i < size; i++) out[i] = in[i] + out[i] - 32;
}

/**	\brief Array convolution with exp.
    \param in       The input array.
    \param out      The output array.
    \param ex       The exp array.
    \param ns       The size of array.
    \param sz       The size of exp.
*/
static inline void convolution_exp(int *in, int *out, int *ex, const int ns, const int sz)
{
    int i, j;
    long long int ij, sum, sum1;

    for(i=0; i < ns; i++){
        if(in[i]){
            sum = 0; sum1 = 0;
            for(j=-sz+1; j < sz; j++){
                ij = i+j;
                if(ij >= 0 && ij < ns){
                    sum += ex[abs(j)]*in[ij];
                    sum1 += ex[abs(j)]*in[ij]*ij;
                    //printf("ij = %d sum1 = %i sum = %i\n", ij, sum1, sum);
                }
            }
            out[i] = sum1/sum;
            //printf("i = %d in = %d out = %d\n",i, in[i], out[i]);
        }
    }
}

/**	\brief Bayer image tone mapping
    \param in       The input 16 bits bayer image.
    \param out      The output 16 bits bayer image.
    \param buff     The temporary buffer.
    \param w        The image width.
    \param h        The image height.
    \param bay		The Bayer grids pattern.
    \param bpp		The pits per pixel.
    \param br       The radius around the pixel.
    \param sd		The standard deviation.
*/
void hdr_tone_bayer(int16 *in, int16 *out, int16 *buff, const int w, const int h, const BayerGrid bay, const int bpp, const int br, const int sd)
{
    int x, y, yw, yx, sz = w*h;
    int i, ns = 1<<bpp, sz1 = 128;
    int *R, *G, *B, *Y; //Histograms
    int *rl, *gl, *bl;  //LUTs
    int ex[sz];
    int16 *iml = buff, *imh = &iml[sz], *buf = &imh[sz];
    int min = 0, max = 0;

    //Prepare histograms and LUTs
    /*
    R = buff; G = &R[ns]; B = &G[ns]; Y = &B[ns];
    rl = &Y[ns]; gl = &rl[ns]; bl = &gl[ns];

    utils_fill_hist_bayer(in, R, G, B, Y, &bl[ns], w, h, bay, bpp);

    for(i=0; i < 1000; i++) printf("%4d  %d\n", i, Y[i]);
    for(i=0; i < ns; i++) tmp += R[i];
    printf("red w*h = %d  hist = %d\n", w*h>>2, tmp);
    tmp = 0;
    for(i=0; i < ns; i++) tmp += G[i];
    printf("green w*h = %d  hist = %d\n", w*h>>1, tmp);
    tmp = 0;
    for(i=0; i < ns; i++) tmp += B[i];
    printf("blue w*h = %d  hist = %d\n", w*h>>2, tmp);
    tmp = 0;
    for(i=0; i < ns; i++) tmp += Y[i];
    printf("blue w*h = %d  hist = %d\n", w*h>>2, tmp);

    //for(i=0; i< 500; i++) printf("%d Y = %d\n", i, Y[i]);

    utils_lut_exp(ex, sd, sz);
    //for(i=0; i < sz; i++) printf("ex = %d\n", ex[i]);

    convolution_exp(R, rl, ex, ns, sz1);
    convolution_exp(G, gl, ex, ns, sz1);
    convolution_exp(B, bl, ex, ns, sz1);
    */
    //for(i=0; i < 500; i++) printf("%4d %d\n", i, rl[i]);

    // Low frequency iml
    filters_bilateral_denoise_bayer(in, iml, buf, br, sd, w, h);

    // High frequency imh
    for(i=0; i < sz; i++) imh[i] = in[i] - iml[i];

    //Reduce contrast
    hdr_ace(iml, iml, (int*)buf, w, h, bpp, 8, BAYER);

    //for(i=0; i < sz; i++) out[i] = imh[i] + 128;


    //Making HDR image
    for(i=0; i < sz; i++) {
        out[i] = imh[i] + iml[i];
        if     (out[i] > max) max = out[i];
        else if(out[i] < min) min = out[i];
    }
    printf("min = %d max = %d\n", min, max);

    //rang = (max - min);

    for(i=0; i < sz; i++) out[i] = (out[i] - min)*255/(max - min);

    /*
    for(y = 0; y < h; y++){
        yw = y*w;
        for(x = 0; x < w; x++){
            yx = yw + x;

            switch(bay){
            case(BGGR):{
                if(y&1){
                    if(x&1) out[yx] = rl[in[yx]];
                    else    out[yx] = gl[in[yx]];
                } else {
                    if(x&1) out[yx] = gl[in[yx]];
                    else    out[yx] = bl[in[yx]];
                }
                break;
            }
            case(GRBG):{
                if(y&1){
                    if(x&1) out[yx] = gl[in[yx]];
                    else    out[yx] = bl[in[yx]];

                } else {
                    if(x&1) out[yx] = rl[in[yx]];
                    else    out[yx] = gl[in[yx]];
                }
                break;
            }
            case(GBRG):{
                if(y&1){
                    if(x&1) out[yx] = gl[in[yx]];
                    else    out[yx] = rl[in[yx]];

                } else {
                    if(x&1) out[yx] = bl[in[yx]];
                    else    out[yx] = gl[in[yx]];
                }
                break;
            }
            case(RGGB):{
                if(y&1){
                    if(x&1) out[yx] = bl[in[yx]];
                    else    out[yx] = gl[in[yx]];
                } else {
                    if(x&1) out[yx] = gl[in[yx]];
                    else    out[yx] = rl[in[yx]];
                }
                break;
            }
            //printf("Start in = %d out = %d blm = %d\n", in[yx], out[yx], blm);
            }
        }
    }
    */
}

void gamma_table(int in, int out)
{
    double gam[512], gam1[512], gam2[512], ia, a = 0.02, a1 = 20, b;
    uint32 i, j, gamma[512], vl0, vl1, hdr = 0;
    in = 511; out = 1023;

    if(hdr == 0){
        for(i=0; i < 512; i++){
            //ia = in*a;
            gam[i] = out*((log(i + in*a) - log(in*a))/(log(in + in*a) - log(in*a)));
        }

        vl0 = 0;
        for(i=1; i < 512; i++){
            vl1 = gam[i];
            gamma[i-1] = (vl0<<10) | (vl1 - vl0);
            //gamma[i-1] = vl0;
            //printf("%3d  %3d\n", i, vl0);
            printf("%7d, ", gamma[i-1]);
            if(!(i&15)) printf("\n");
            vl0 = vl1;
        }
        gamma[511] = vl0<<10;
        printf("%7d, ", gamma[511]);
        //printf("%3d  %3d,", vl0, (vl1 - vl0));
    } else if (hdr == 1){
        for(i=0; i < 512; i++){
            //ia = in*a;
            gam1[i] = out*((log(i + in*a) - log(in*a))/(log(in + in*a) - log(in*a)));
            gam2[i] = (out - out*((log((511-i) + in*a1) - log(in*a1))/(log(in + in*a1) - log(in*a1))));
            gam[i] = (gam1[i]*(511-i) +  gam2[i]*i)/511;
        }

        vl0 = 0;
        for(i=1; i < 512; i++){
            vl1 = gam[i];
            gamma[i-1] = (vl0<<10) | (vl1 - vl0);
            //rintf("%3d   %3d  %3d\n", i-1, vl0, (vl1 - vl0));
            printf("%7d, ", gamma[i-1]);
            if(!(i&15)) printf("\n");
            vl0 = vl1;
        }
        gamma[511] = vl0<<10;
        printf("%7d, \n", gamma[511]);
        //printf("%3d   %3d  %3d\n", 511, vl0, (vl1 - vl0));

    } else if(hdr == 2){
        float px[6], py[6];
        px[0] = 0; px[1] = 20;  px[2] = 40;  px[3] = 512;  //px[3] = px[5] - 20; px[4] = px[5] - 10;
        py[0] = 0; py[1] = 200; py[2] = 300; py[3] = 1023; //py[3] = py[5] - 40; py[4] = py[5] - 20;

        for(j=0; j < 3; j++){
            a = (py[j] - py[j+1])/(px[j] - px[j+1]);
            b = (px[j]*py[j+1] - px[j+1]*py[j])/(px[j] - px[j+1]);
            for(i=px[j]; i < px[j+1]; i++){
                gam[i] = a*i + b;
                printf("i = %d gam = %f a = %f b = %f j = %d\n", i, gam[i], a, b, j);
            }
        }

        vl0 = 0;
        for(i=1; i < 512; i++){
            vl1 = gam[i];
            gamma[i-1] = (vl0<<10) | (vl1 - vl0);
            //printf("%3d  %3d\n", i, vl0);
            printf("%7d, ", gamma[i-1]);
            if(!(i&15)) printf("\n");
            vl0 = vl1;
        }
        gamma[511] = vl0<<10;
        printf("%7d, ", gamma[511]);
        //printf("%3d  %3d,", vl0, (vl1 - vl0));

    }
}
