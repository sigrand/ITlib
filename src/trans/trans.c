#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libit/types.h"
#include "./trans.h"

#define PI 3.1415926535

struct save {
    double h;
    double s[2];
    double M[5];
    double P[5];
    double R[5];
    double MT;
    double PT;
};

//Calculate square of the coil,
//s - square of one wire in mm**2
//n - number of coils
inline double sq(double s, double n)
{
    return n*s;
}

//Calculate size of wire,
//s - square of one wire in mm**2
inline double ln(double s)
{
    return sqrt(s);
}

//Calculate number of ring in vertical direction
//h - vertical size in mm
//s - square of one wire in mm**2
inline double Nh(double h, double s)
{
    return h/ln(s);
}

//Calculate number of ring radial direction
//h - vertical size in mm
//s - square of one wire in mm**2
//n - number of coils
inline double Nr(double h, double s, double n)
{
    return sq(s,n)/(h*ln(s));
}

//Calculate lenght of coil in m
//h - vertical size in mm
//s - square of one wire in mm**2
//n - number of coils
//r - inside radius of coil
inline double Ln(double h, double s, double n, double r)
{
    return 2.*PI*Nh(h,s)*Nr(h,s,n)*(r + (Nr(h,s,n)-1)*ln(s)/2.)/1000.;
}

//Calculate radial size of the coil
//h - vertical size in mm
//s - square of one wire in mm**2
//n - number of coils
inline double Rr(double h, double s, double n)
{
    return sq(s,n)/h;
}

//Calculate loss power in coil in W
//h - vertical size in mm
//s - square of one wire in mm**2
//n - number of coils
//r - inside radius of coil
//C - Conductivity
//I - Current
inline double Pw(double h, double s, double n, double r, double C, double I)
{
    return C*Ln(h, s, n, r)*I*I/s;
}

//Calculate valume
//h - vertical size in mm
//R1 - internal radius in mm
//R2 - External radius
inline double Vol(double h, double R1, double R2)
{
    return PI*h*(R2*R2 - R1*R1)/1000.;
}

//Calculate valume outside
//R1 - internal radius in mm
//R2 - External radius
inline double Vou(double R1, double R2)
{
    return PI*PI*R1*R1*R2/1000.;
}

//Calculate mass
//h - vertical size in mm
//R1 - internal radius in mm
//R2 - External radius im mm
//D - Density im kg/dm**3
inline double Ms(double h, double R1, double R2, double D)
{
    return Vol(h, R1, R2)*D/1000.;
}

void trans(void)
{
    double Cal = 0.0282; //Conductivity aluminum OM*m/mm**2
    double Dal = 2.6989; //Density aluminum kg/gm**3
    double Dfe = 7.874; //Density of transformer steel kg/gm**3
    double R[5], V[5], M[5], N[2], I[2], s[2], P[2];
    double L1, Pw1, h = 300., s1 = 5.;
    double MT, PT, Mmin=100000;
    struct save min;


    R[0] = 60.; //Radius of magnetic core mm
    N[0] = 2000.; //the number of coils in 10Kv coil
    N[1] = 80.; //the number of coils in 400v coil
    I[0] = 10./sqrt(3.); //the max current of 10Kv coil
    I[1] = 250./sqrt(3.); //the max current of 400v coil

    /*
    L1 = Ln(h, s1, N[0], R[0]);
    Pw1 = Pw(h, s1, N[0], R[0], Cal, I[0]);
    V[1] = Vol(h, R[0], R[0]+Rr(h, s1, N[0]));
    M[1] = Ms(h, R[0], R[0]+Rr(h, s1, N[0]), Dal);

    printf("The length of coil = %f M Nh = %f Nr = %f ln(s) = %f Power = %f Vol = %f M = %f\n", L1, Nh(h,s1), Nr(h,s1, N[0]), ln(s1), Pw1, V[1], M[1]);

    h = 300.; s[0] = 5; s[1] = 50;
    */

    for(h=200.; h < 1000.; h+=10.){
        for(s[0]=3.; s[0] < 30.; s[0]+=0.2){
            for(s[1]=50.; s[1] < 600.; s[1]+=1.){

                R[1] = R[0]+Rr(h, s[0], N[0]);
                R[2] = R[1]+Rr(h, s[1], N[1]);
                M[0] = (Vol(h, 0, R[0]) + Vou(R[0],R[2]))*Dfe/1000.;            //Mass steel
                M[2] = Vol(h, R[0], R[1])*Dal/1000.;                            //Mass 10kV coil
                M[4] = Vol(h, R[1], R[2])*Dal/1000.;                            //Mass 400V coil
                MT = 3.*(M[0] + M[2] + M[4]);

                P[0] = M[0];                                //Loss in steel
                P[2] = Pw(h, s[0], N[0], R[0], Cal, I[0]);  //Loss in 10kV coil
                P[4] = Pw(h, s[1], N[1], R[1], Cal, I[1]);  //Loss in 400V coil
                PT = 3.*(P[0] + P[2] + P[4]);

                if(PT <= 1000) {
                    if(MT < Mmin) {
                        Mmin = MT;
                        min.h = h;
                        min.s[0] = s[0];
                        min.s[1] = s[1];
                        min.M[0] = M[0];
                        min.M[2] = M[2];
                        min.M[4] = M[4];
                        min.P[0] = P[0];
                        min.P[2] = P[2];
                        min.P[4] = P[4];
                        min.R[1] = R[1];
                        min.R[2] = R[2];
                        min.MT = MT;
                        min.PT = PT;
                    }

                }
                /*
                printf("h = %f s1 = %f s2 = %f PT = %f MT = %f\n",h, s[0], s[1], PT, MT);
                printf("R1 = %f R2 = %f M0 = %f M2 = %f M4 = %f P0 = %f P2 = %f P4 = %f\n",
                       R[1], R[2], M[0], M[2], M[4], P[0], P[2], P[4]);
                */

            }
        }
    }
    printf("h = %f s1 = %f s2 = %f PT = %f MT = %f\n",min.h, min.s[0], min.s[1], min.PT, min.MT);
    printf("I0 = %f I1 = %f R1 = %f R2 = %f M0 = %f M2 = %f M4 = %f P0 = %f P2 = %f P4 = %f\n",
           I[0], I[1], min.R[1], min.R[2], min.M[0], min.M[2], min.M[4], min.P[0], min.P[2], min.P[4]);

}

