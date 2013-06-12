#ifndef FILTERS_H_
#define FILTERS_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void filters_median(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int type);
void filters_median_bayer(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int type);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /*FILTERS_HH_*/
