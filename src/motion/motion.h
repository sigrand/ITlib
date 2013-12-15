#ifndef _ALG_MOTION_DETECT_H_
#define _ALG_MOTION_DETECT_H_

typedef struct
{
    int x;
    int y;
} Vector;

typedef struct
{
    Vector sz;  //Size w and h
    Vector cog; //Center Of Gravity
    int nc;     //Number of Cells in the object

}   Object;


#endif
