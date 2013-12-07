#include "alg_motion.h"
//#include <df/log.h>

//#define DFTRACE() dflog(LOG_INFO, "%s():%u", __func__, __LINE__)



void *ALG_motionDetectCreate(ALG_MotionDetectCreate *create, ALG_MotionDetectCreateStatus *status)
{
        ALG_MotionObj  *pObj = NULL;

        pObj = OSA_memAlloc(sizeof(ALG_MotionObj));

        if(pObj==NULL)
            return NULL;

        memset(pObj, 0, sizeof(*pObj));
        if( create )
            pObj->createPrm = *create;

        if( status )
            pObj->createStatus = *status;

        /*Start motion detection function after delay start_cnt frames*/
        pObj->frame_count 	= 0;
        pObj->start_cnt 	= 200;

        return pObj;
}


int ALG_motionDetectGetThres(ALG_MotionObj  *pObj)
{
    int block_size;
    int th_per_pix = pObj->runPrm.motioncvalue; // was: 30
    if( pObj == NULL )
        return ALG_MOTION_S_FAIL;

    pObj->frame_width   = (pObj->runPrm.ImageWidth  >> 4); // Number of macroblocks in frame width
    pObj->frame_height  = (pObj->runPrm.ImageHeight >> 4); // Number of macroblocks in frame height

    /* Set the motion block size base on image size */
    pObj->win_width 	= pObj->runPrm.windowWidth >> 4; //Window width in macroblocks
    pObj->win_height 	= pObj->runPrm.windowHeight >> 4; //Window height in macroblocks

    block_size = pObj->win_width * pObj->win_height;

#if 0 //FIXME: debug only
    dflog(LOG_INFO, "%s():%u Image: %u x %u, frame: %u x %u, window: %u x %u, win: %u x %u", __func__, __LINE__,
          pObj->runPrm.ImageWidth,
          pObj->runPrm.ImageHeight,

          pObj->frame_width,
          pObj->frame_height,

          pObj->runPrm.windowWidth,
          pObj->runPrm.windowHeight,

          pObj->win_width,
          pObj->win_height
         );
#endif

    //+ Windows Motion Enabled and Detected Maps
    size_t frame_square = pObj->frame_height * pObj->frame_width;

    if (pObj->Enabled == NULL)
        pObj->Enabled = calloc(1, frame_square);

    if (pObj->Detected == NULL)
        pObj->Detected = malloc(frame_square);
    memset(pObj->Detected, 0, frame_square);
    //- Windows Motion Enabled and Detected Maps

    pObj->SAD_THRESHOLD = th_per_pix<<8;

    return ALG_MOTION_S_OK;
}

int ALG_motionDetectCalc(ALG_MotionObj  *pObj)
{
    ALG_MotionDetectMbMvInfo *mbMV_data;
    int i, j, status;

    mbMV_data = pObj->runPrm.mbMvInfo + pObj->MvDataOffset;

    pObj->warning_count = 0;

    for (i = 0; i < pObj->win_height; i++)
    {
        for(j = 0; j < pObj->win_width; j++)
        {
            if(mbMV_data->SAD > pObj->SAD_THRESHOLD)
            {
                pObj->warning_count++;
            }
            mbMV_data ++;
        }
        mbMV_data = mbMV_data + (pObj->frame_width - pObj->win_width);
     }

    /* If the pObj->warning_count is bigger than pObj->threshold,
    the function will return a alarm signal*/

    status = ( pObj->warning_count >= pObj->threshold ) ? ALG_MOTION_S_DETECT : ALG_MOTION_S_NO_DETECT;

    return status;
}

int ALG_motionDetectStart(ALG_MotionObj  *pObj)
{
    int detect_cnt = 0, avr = 0, max = 0, min = 0xfffffff, Sad_Threshold;
    size_t x, y, yx, w = pObj->frame_width, h = pObj->frame_height, sz = w*h;
    char *dc = pObj->Detected;
    char *en = pObj->Enabled;
    ALG_MotionDetectMbMvInfo *mbMV_data;
    mbMV_data = pObj->runPrm.mbMvInfo;

    //dflog(LOG_INFO, "%s():%u %u x %u", __func__, __LINE__, w, h);
    pObj->MvDataOffset = 0;

    //Get noise statistics and calculate threshold
    for(y=0; y < h; y++) {
        for(x=0; x < w; x++) {
            yx = y*w + x;
            avr += mbMV_data->SAD;
            if      (max < mbMV_data->SAD) max = mbMV_data->SAD;
            else if (min > mbMV_data->SAD) min = mbMV_data->SAD;
        }
    }
    min = (min>>8);
    max = (max>>8);
    avr = (avr>>8)/sz;
    Sad_Threshold = avr*4;
    //dflog(LOG_INFO, "%s():%u  min = %d avr = %d max = %d Sad_Threshold = %d", __func__, __LINE__, min, avr, max, Sad_Threshold);

    //Find active cells
    for(y=0; y < h; y++) {
        for(x=0; x < w; x++) {
            yx = y*w + x;
            //FXIME: check pObj->Enabled[yx] when it will be supplied
            if (1 /* en[yx] */)  {
                if(mbMV_data->SAD > Sad_Threshold){
                    dc[yx] = 1; //Motion Detected
                    detect_cnt++;
                } else {
                    dc[yx] = 0; //Motion Not Detected
                }
                mbMV_data++;
            }
        }
    }

    //Check if MB have two or more neighborhood
    if(mbMV_data){
        for(y=1; y < h-1; ++y) {
            for(x=1; x < w-1; ++x) {
                size_t cn = 0;
                yx = y*w + x;
                if (1 /* en[yx] */)  {
                    if (dc[yx]) {
                        if (dc[yx-1  ]) ++cn;
                        if (dc[yx-1-w]) ++cn;
                        if (dc[yx-w  ]) ++cn;
                        if (dc[yx+1-w]) ++cn;
                        if (dc[yx+1  ]) ++cn;
                        if (dc[yx+1+w]) ++cn;
                        if (dc[yx+w  ]) ++cn;
                        if (dc[yx-1+w]) ++cn;
                        if (cn > 2) return ALG_MOTION_S_DETECT;
                    }
                }
            }
        }

        //Find clusters
        for(y=1; y < h-1; ++y) {
            for(x=1; x < w-1; ++x) {
                size_t cn = 0;
                yx = y*w + x;
                if (1 /* en[yx] */)  {
                    if (dc[yx]) {
                        if (dc[yx-1  ]) ++cn;
                        if (dc[yx-1-w]) ++cn;
                        if (dc[yx-w  ]) ++cn;
                        if (dc[yx+1-w]) ++cn;
                        if (dc[yx+1  ]) ++cn;
                        if (dc[yx+1+w]) ++cn;
                        if (dc[yx+w  ]) ++cn;
                        if (dc[yx-1+w]) ++cn;
                        if (cn > 2) return ALG_MOTION_S_DETECT;
                    }
                }
            }
        }
    }

#if 0 // FIXME: debug only
    for(y=0; y < h; y++) {
        for(x=0; x < w; x++) {
            yx = y*w + x;
            dflog_(LOG_INFO, "%1d ", pObj->runStatus.windowMotionDetected[yx]);
        }
        dflog_flush(LOG_INFO);
    }
#endif

    return ALG_MOTION_S_NO_DETECT;
}

int ALG_motionDetectRun(void *hndl, ALG_MotionDetectRunPrm *prm, ALG_MotionDetectRunStatus *status)
{
    ALG_MotionObj  *pObj = (ALG_MotionObj*)hndl;

    if( pObj == NULL )
        return ALG_MOTION_S_FAIL;

    /*Parm tranfer*/
    pObj->runPrm 	= *prm;
    pObj->runStatus = *status;
    pObj->frame_count++;

    if((pObj->runPrm.isKeyFrame == TRUE) || (pObj->frame_count < pObj->start_cnt))
    {
        return ALG_MOTION_S_NO_DETECT;
    }

    ALG_motionDetectGetThres(pObj);

    return ALG_motionDetectStart(pObj);
}

int ALG_motionDetectDelete(void *hndl)
{
    ALG_MotionObj *pObj=(ALG_MotionObj *)hndl;

    if(pObj==NULL)
        return ALG_MOTION_S_FAIL;

    //Free Windows Motion Enable and Detected Maps
    OSA_memFree(pObj->Detected);
    OSA_memFree(pObj->Enabled);

    OSA_memFree(pObj);

    return ALG_MOTION_S_OK;
}
