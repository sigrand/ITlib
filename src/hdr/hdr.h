#ifndef HDR_H_
#define HDR_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */
void hdr_ace(const int16 *in, int16 *out, int *buff, const int w, const int h, const int bpp, const int bpp1, const int cs);
void hdr_diff(const int16 *in, const int16 *avr, int16 *dif, const int w, const int h, const int bpp);
void hdr_ace_local(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int bpp);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /*HDR_HH_*/
