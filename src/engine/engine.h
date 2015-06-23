#ifndef ENGINE_H_
#define ENGINE_H_

#include "../libit/types.h"

#define PI 3.1415926535
//#define E 2.71828182
#define MU 1.256637E-6


typedef struct COIL {
    double R;   //Resistence
    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double D;   //Density kg/gm**3
    double T;   //The thickness of the insulation mm
    double R1, R2; //Internal and external radius
    double L;   //Width
    double FF;  //Filing factor
    double C;   //Conductivity
    double N;   //Number for turn
    double S;   //square
} COIL;

typedef struct CORE {
    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double R1, R2; //Internal and external radius
    double D;   //Density kg/gm**3
    double L;   //Width
    double w;   //Teeth width
    double p;   //Period lenght
    double a;   //Step

} CORE;

typedef struct ROTOR {
    CORE CR;

    double R1, R2; //Internal and external radius
    double M;   //Mass kg
    double V;   //Volume dm**3
    double D;   //Density kg/gm**3
    double L;   //Length
    double N;  //Numberes of teeth

} ROTOR;

typedef struct STATOR {
    CORE CR;
    COIL CL;

    double R1, R2; //Internal and external radius
    double M;   //Mass kg
    double V;   //Volume dm**3
    double D;   //Density kg/gm**3
    double L;   //Length
    double N;   //Numberes of teeth

} STATOR;

typedef struct ENGINE {
    ROTOR  RT;
    STATOR ST;

    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double PW;  //Power in kW
    double h;   //The height of the coil with spiral
    double L;   //The length of engine
    double I;   //Current
    double D;   //Density kg/gm**3
    double EX;  //Design increase in persent
    double R;   //External radius
    double Ph;  //Number of pahses
    double Pn;  //Groups of pahses

} ENGINE;


#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

void trans(void);

#ifdef __cplusplus
}
#endif /*__cplusplus*/
#endif /*TRANSFORMS_H_*/
