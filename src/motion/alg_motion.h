#ifndef _ALG_MOTION_DETECT_H_
#define _ALG_MOTION_DETECT_H_

#include "alg.h"

#define ALG_MOTION_DETECT_HORZ_REGIONS 	(4)
#define ALG_MOTION_DETECT_VERT_REGIONS 	(3)
#define ALG_MOTION_DETECT_REGIONS   	(ALG_MOTION_DETECT_HORZ_REGIONS*ALG_MOTION_DETECT_VERT_REGIONS)

//Following values are Sensitivity levels that can come from UI
#define ALG_MOTION_DETECT_SENSITIVITY_LOW       	0
#define ALG_MOTION_DETECT_SENSITIVITY_MEDIUM    	1
#define ALG_MOTION_DETECT_SENSITIVITY_HIGH         	2

enum{
    ALG_MOTION_S_DETECT = 100,
    ALG_MOTION_S_NO_DETECT ,

    ALG_MOTION_S_FAIL 	= OSA_EFAIL,
    ALG_MOTION_S_OK 	= OSA_SOK
};

typedef struct
{
    short MVx;
    short MVy;
    int   SAD;
}ALG_MotionDetectMbMvInfo;

typedef struct {

  Uint16 codec;
  Uint16 width;
  Uint16 height;
  Uint16 numRegions;

} ALG_MotionDetectCreate;

typedef struct {

  Uint16 numMbH;
  Uint16 numMbV;

} ALG_MotionDetectCreateStatus;

typedef struct {

    ALG_MotionDetectMbMvInfo *mbMvInfo;

    int ImageWidth;
    int ImageHeight;
    int windowWidth;
    int windowHeight;
    int isKeyFrame;
    int codecType;
    int	isDateTimeDraw;
    int motionenable;
    int motioncenable;
    int motioncvalue;
    int motionlevel;
    int motionsens;
    int blockNum;
    char windowEnableMap[ALG_MOTION_DETECT_HORZ_REGIONS][ALG_MOTION_DETECT_VERT_REGIONS];

} ALG_MotionDetectRunPrm;

//The following structure will store the information for motion detected for a particular selected window.
typedef struct {
    char windowMotionDetected[ALG_MOTION_DETECT_HORZ_REGIONS][ALG_MOTION_DETECT_VERT_REGIONS];
}ALG_MotionDetectRunStatus;

typedef struct
{
    int frame_width;
    int	frame_height;
    int frame_count;
    int start_cnt;
    int SAD_THRESHOLD;
    int threshold;
    int win_width;
    int win_height;
    int MvDataOffset;
    int warning_count;
    ALG_MotionDetectCreate 				createPrm;
    ALG_MotionDetectCreateStatus	createStatus;
    ALG_MotionDetectRunPrm				runPrm;
    char *Detected;
    char *Enabled;
    ALG_MotionDetectRunStatus runStatus;

}ALG_MotionObj;

typedef struct
{
    int x;
    int y;
}Vector;

typedef struct
{
    Vector top_left;
    Vector down_right;
    Vector CenterOfGravity;
    int NCells; //Number of Cells in the object
}Object;

void *ALG_motionDetectCreate(ALG_MotionDetectCreate *create, ALG_MotionDetectCreateStatus *status);
int ALG_motionDetectRun(void *hndl, ALG_MotionDetectRunPrm *prm, ALG_MotionDetectRunStatus* status);
int ALG_motionDetectDelete(void *hndl);

#endif
