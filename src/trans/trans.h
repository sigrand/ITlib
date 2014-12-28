#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

void trans(const int16 *in, int16 *rgb, int16 *buff, const int w, const int h, const BayerGrid bay);

#ifdef __cplusplus
}
#endif /*__cplusplus*/
#endif /*TRANSFORMS_H_*/
