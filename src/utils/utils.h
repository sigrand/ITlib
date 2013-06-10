#ifndef UTILS_H_
#define UTILS_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void utils_cnange_bytes(int16 *in, const int w, const int h);
void utils_get_stat(int16 *in, const int w, const int h, int *bpp, int *min, int *max);
uint8* utils_16_to_8(const int16 *in, uint8 *out, const int w, const int h, const int bpp, const int par);

void utils_zoom_out_rgb16_to_rgb16(const int16 *in, int16 *out, uint32 *buff, const int w, const int h, const int zoom);
void utils_zoom_out_bayer16_to_rgb16(const uint16 *in, uint16 *out, uint32 *buff, const int w, const int h, const int zoom, const BayerGrid bay);

void utils_wb(int16 *in, int *rm, int *bm, uint32 *buff, const int w, const int h, const int sh, const int bpp);
void utils_wb_rgb(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bpp);
void utils_wb_bayer(const int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bpp, const int bg);

void utils_integral_grey(const int16 *in, int *ing, const int w, const int h);
void utils_integral_bayer(const int16 *in, uint32 *ing, const int w, const int h);
void utils_integral_grey_br(int16 *in, uint32 *ing, uint32 *buff, const int w, const int h, const int br);

void utils_average(int16 *in, int16 *out, uint32 *buff, const int w, const int h, const int br);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /*UTILS_HH_*/
