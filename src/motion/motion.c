#include "motion.h"
#include "../libit/types.h"
#include "../utils/utils.h"

/** \brief Find the objects in integral matrix
    \param ing 	The pointer to a integral matrix, the size = (w + 2*br)*(h + 2*br)*(sizeof(uint32)).
    \param w	The image width.
    \param h	The imahe height.
    \param br	The border on the image ing, new w = w + 2*br and new h = h + 2*br;
*/
void find_object(uint32 *ing, uint16 *out, const int w, const int h, const int br, Object *obj)
{
    int x, y, yw, yx, wb = w + br, hb = h + br;
    int w1 = w + br*2, ox = (obj->sz.x-1)>>1, oy = (obj->sz.y-1)>>1;


    for(y=br; y < hb; y++){
        yw = y*w1;
        for(x=br; x < wb; x++){
            yx = yw + x;
            out[yx] = ing[yx + ox + oy*w1] + ing[yx - (ox+1) - (oy+1)*w1] - ing[yx + ox - (oy+1)*w1] - ing[yx - (ox+1) + oy*w1];
        }
    }
}
