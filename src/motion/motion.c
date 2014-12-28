#include <stdio.h>
#include <string.h>

#include "motion.h"
#include "../libit/types.h"
#include "../utils/utils.h"

/** \brief Find the objects in integral matrix
    \param ing 	The pointer to a integral matrix, the size = (w + 2*br)*(h + 2*br)*(sizeof(uint32)).
    \param w	The image width.
    \param h	The imahe height.
    \param br	The border on the image ing, new w = w + 2*br and new h = h + 2*br;
*/
void find_object(uint32 *ing, uint16 *out, uint32 *hist, const int w, const int h, const int br, Object *obj)
{
    int i, x, y, yw, yw1, yx, yx1, w1 = w + br*2, cn;
    int ox = obj->sz.x, oy = obj->sz.y, oxh = (obj->sz.x>>1)+1, oyh = (obj->sz.y>>1)+1;
    int dr[8] = {-1, -1+w, -w, 1-w, 1, 1+w, w, -1+w};
    //uint16 lm[w*h];     //Local max array
    //uint32 hist[2048];  //Histogram for adaptive threshould
    Object oba[64];     //Objects array
    Object *obp[w*h];

    memset(hist, 0, sizeof(hist));

    for(y=0; y < h; y++){
        yw = y*w;
        yw1 = (y + br - oyh)*w1;
        for(x=0; x < w; x++){
            yx = yw + x;
            yx1 = yw1 + (x + br - oxh);
            out[yx] = (ing[yx1 + ox + oy*w1] + ing[yx1] - ing[yx1 + ox] - ing[yx1 + oy*w1]);
            //Update histogram
            hist[out[yx] ]++;
        }
    }

    //Find local max
    for(y=1; y < h-1; y++){
        yw = y*w;
        for(x=1; x < w; x++){
            yx = yw + x;
            cn = 0;
            for(i=0; i < 8; i++) if(out[yx] >= out[yx + dr[i]]) cn++;
            //Local max
            if(cn == 8) {

            }
        }
    }
}


