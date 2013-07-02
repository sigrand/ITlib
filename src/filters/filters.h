#ifndef FILTERS_H_
#define FILTERS_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void filters_median(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int type);
void filters_median_bayer(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int type);

void filters_NLM_denoise_bayer(int16 *in, int16 *avr, int16 *out, int16 *buff, const int st, const int sg, const int w, const int h);

void filters_hessian(int16 *in, int16 *out, uint32 *buff, const int w, const int h);

void filters_denoise_regression_bayer(int16 *in, int16 *out, int *buff, const int br, const int w, const int h);

void filters_MSE_bayer(int16 *in, int16 *avr, int16 *out, int16 *buff, const int br, const int bpp, const int w, const int h);

void  b_spline_aproximation(int16 *in, int16 *out, int16 *buff, const int w, const int h);
#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /*FILTERS_HH_*/
