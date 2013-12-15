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
    int x, y, yw, yw1, yx, yx1, w1 = w + br*2, h1 = h + br*2;
    int ox = obj->sz.x, oy = obj->sz.y, oxh = (obj->sz.x>>1)+1, oyh = (obj->sz.y>>1)+1;


    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = (y + br - oyh)*w1;
        for(x=0; x < w; x++){
            yx = yw + x;
            yx1 = yw1 + (x + br - oxh);
            out[yx] = (ing[yx1 + ox + oy*w1] + ing[yx1] - ing[yx1 + ox] - ing[yx1 + oy*w1]);
        }
    }
}
