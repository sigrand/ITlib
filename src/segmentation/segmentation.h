#ifndef SEG_H_
#define SEG_H_

#include "../libit/types.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void seg_gradient(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int th);
uint32 seg_local_max(int16 *in, int16 *out, int16 *buff, const int w, const int h);
void seg_edge_detection(int16 *in, uint8 *con, const int w, const int h);
void seg_canny_edge(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int th);
uint32 seg_end_of_edges(int16 *in, int16 *out, int16 *buff, const int w, const int h);
void seg_corners(int16 *in, int16 *out, int16 *buff, const int w, const int h, const int th);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /*SEG_H_*/
