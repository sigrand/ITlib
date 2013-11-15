#ifndef UTILS_H_
#define UTILS_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void utils_cnange_bytes(int16 *in, const int w, const int h);
void utils_get_stat(int16 *in, const int w, const int h, int *bpp, int *min, int *max);
uint8* utils_16_to_8(const int16 *in, uint8 *out, const int w, const int h, int bpp, const int par);

void utils_zoom_out_rgb16_to_rgb16(const int16 *in, int16 *out, uint32 *buff, const int w, const int h, const int zoom);
void utils_zoom_out_bayer16_to_rgb16(const uint16 *in, uint16 *out, uint32 *buff, const int w, const int h, const int zoom, const BayerGrid bay);

void utils_wb(int16 *in, int *rm, int *bm, int *mp, uint32 *buff, const int w, const int h, const int sh, const int bpp);
void utils_wb_rgb(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bpp);
void utils_wb_bayer(const int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bpp, const int bg);

void utils_integral(int16 *in, uint32 *ing, uint32 *buff, const int w, const int h, const int br);
void utils_integral_bayer(int16 *in, uint32 *ing, uint32 *buff, const int w, const int h, const int br);

void utils_average(int16 *in, int16 *out, uint32 *buff, const int w, const int h, const int br);
void utils_average_bayer(int16 *in, int16 *out, uint32 *buff, const int w, const int h, const int br);

void utils_subtract(const int16 *in, const int16 *in1, int16 *out, const int w, const int h, const int bpp);
void utils_add(const int16 *in, const int16 *in1, int16 *out, const int w, const int h, const int bpp, const int bpp1);

void utils_resize_down_2(const int16 *in, int16 *out, int16 *buff, const int w, const int h);

void utils_fill_hist_bayer(const uint16 *in, int *R, int *G, int *B, int *Y, int *buff, const int w, const int h, const BayerGrid bay, const int bpp);
void utils_lut_exp(int *ex, const int sd, const int sz);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /*UTILS_HH_*/
