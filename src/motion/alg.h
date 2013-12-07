#ifndef _ALG_H_
#define _ALG_H_

#include <osa.h>
#include <drv.h>

//#define ALG_DEBUG

#define ALG_VID_DATA_FORMAT_YUV422    (DRV_DATA_FORMAT_YUV422)
#define ALG_VID_DATA_FORMAT_YUV420    (DRV_DATA_FORMAT_YUV420)

#define ALG_VID_DATA_FORMAT_PROGRESSIVE  (0<<4)
#define ALG_VID_DATA_FORMAT_INTERLACED   (1<<4)

typedef enum{
    ALG_VID_CODEC_H264=0,
    ALG_VID_CODEC_MPEG4,
    ALG_VID_CODEC_MJPEG,
    ALG_VID_CODEC_VNF,
    ALG_VID_CODEC_NONE=-1
} ALG_VID_CODEC_TYPE;

typedef enum{
    ALG_AUD_CODEC_G711=0,
    ALG_AUD_CODEC_AAC,
    ALG_AUD_CODEC_SWAAC,
    ALG_AUD_CODEC_NONE=-1
} ALG_AUD_CODEC_TYPE;

int ALG_sysInit(void);
int ALG_sysExit(void);

#endif
