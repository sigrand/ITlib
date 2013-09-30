#ifndef STEREO_H_
#define STEREO_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void stereo_disparity(const int16 *in, const int16 *in1, int16 *out, int16 *buff, const int w, const int h);
void stereo_filter(const int16 *in, int16 *out, int16 *buff, const int w, const int h);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /*STEREO_HH_*/
