#ifndef UTILS_H_
#define UTILS_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void utils_image_copy_n(uint8 *in, const int w, const int h);
void utils_get_stat(int16 *in, const int w, const int h, int *bpp, int *min, int *max);
uint8* utils_16_to_8(const int16 *in, uint8 *out, const int w, const int h, const int bpp, const int par);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /*UTILS_HH_*/
