#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

uint8* utils_bay16_to_rgb8_bi(const int16 *in, uint8 *rgb, int16 *buff, const int w, const int h, const BayerGrid bay, const int bpp);

#ifdef __cplusplus
}
#endif /*__cplusplus*/
#endif /*TRANSFORMS_H_*/
