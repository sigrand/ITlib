#ifndef UTILS_H_
#define UTILS_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void utils_image_copy_n(const uint8 *in, int16 *out, uint32 w, uint32 h, uint32 bpp);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /*UTILS_HH_*/
