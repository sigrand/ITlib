#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

#include "../libit/types.h"

typedef struct COIL {
    double U;   //Voltage V
    double I;   //Current A
    double T;   //The thickness of the wire insulation mm
    double R;   //The external radius mm
    double s;   //The cross-sectional area of the wire mm**2
    double s1;   //The cross-sectional area of the wire mm**2 without isolation
    double L;   //The lenght of coil wire m
    double L1;   //The lenght of coil wire m
    double N;   //The number of turns
    double N1, N2;   //The number of turns
    double Nr;   //The number of turns in radial derection
    double Nh;   //The number of turns vertical direction
    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double Rz;   //Resistance
    double Rl;  //Induction Resistance
    double h;   //The wire height
    double w;   //The wire width with spiral
    double w1;   //The wire width
    double C;   //Conductivity
    double D;   //Density aluminum kg/gm**3
    double S;   //square м**2
    double Sw;   //The surface area of the wire m**2
    double a,b; //Elips parameters
    double sp[2]; //The cross-sectional area of the wire range
    double Np[2]; //The number of turns range and step
    double Rp[2];   //Radius search paraments
} COIL;

typedef struct MCORE {
    double R;   //The external radius mm
    double R1;  //The tor radius mm
    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double D;   //Density of transformer steel
    double Sfc; //Core Stacking Factor
    double Lc;  //Loss in magnetic core W/kg
    double Mu;  //Magnetic permeability of transformer steel
    double B;   //The max magnetic field T
    double B1;   //The max magnetic field T
    double H;   //The magnetic industion
    double S;   //square м**2
    double a,b; //Elips parameters
    double Rp[2];   //Radius search paraments
    double R1p[2];   //Tor Radius search paraments
} MCORE;

typedef struct INS {
    double R;   //The external radius mm
    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double D;   //Density of insulation kg/gm**3
    double T;   //The thickness of the insulation mm
    double S;   //square м**2
} INS;

typedef struct TRANS {
    MCORE m;     //Magnetic core
    INS   i[2];  //Insulation
    COIL  c[2];  //First coil

    double PW;  //Power in kW
    double NN;  //Transformation coefficient
    double M;   //Mass kg
    double V;   //Volume dm**3
    double S;   //square м**2
    double P;   //The loss power W
    double H;   //The height of the coil with spiral
    double I;   //No-load current
    double W;   //The distance between the coils
    double T;   //Temperature of transformer
    double Hp[2];   //The height of the coilrange mm
    double Ks[2];   //Skin effect paramrters
} TRANS;

typedef double (*fp)(double, double);
typedef double (*fp1)(double , double ,  double , double );

enum Trans_types
{
    Round_3phase        = 0,
    Round_120_3phase	= 1,
    Tore_1phase         = 2
};


#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

void trans(void);

#ifdef __cplusplus
}
#endif /*__cplusplus*/
#endif /*TRANSFORMS_H_*/
