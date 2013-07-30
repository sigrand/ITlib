#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

void utils_bay_to_rgb_bi(const int16 *in, int16 *rgb, int16 *buff, const int w, const int h, const BayerGrid bay);
void utils_bay_to_grey_bi(const int16 *in, int16 *out, int16 *buff, const int w, const int h, const BayerGrid bay);
void trans_rgb_to_yuv444(const int16 *rgb, int16 *Y, int16 *U, int16 *V, const uint32 w, const uint32 h);
void trans_rgb_to_yuv420(const int16 *rgb, int16 *Y, int16 *U, int16 *V, const uint32 w, const uint32 h);
void  trans_bay_to_rgb_b_spline(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bay);

#ifdef __cplusplus
}
#endif /*__cplusplus*/
#endif /*TRANSFORMS_H_*/
