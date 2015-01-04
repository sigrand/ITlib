#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libit/types.h"
#include "./trans.h"

#define PI 3.1415926535
#define MU 1.256637E-6

struct save {
    double h;
    double s[2];
    double M[5];
    double P[5];
    double R[5];
    double L[2];
    double N[2];
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

//Calculate number of ring radial direction
//h - vertical size in mm
//s - square of one wire in mm**2
//n - number of coils
inline double Nr1(double h, double s, double n)
{
    return round(sq(s,n)/(h*ln(s)));
}

//Calculate size of wire
//s - square of one wire in mm**2
inline double lr(double h, double s, double n, double Nr)
{
    return sqrt(s);
}

//Calculate lenght of coil in m
//h - vertical size in mm
//s - square of one wire in mm**2
//n - number of coils
//r - inside radius of coil
inline double Ln(double h, double s, double n, double r)
{
    return 2.*PI*Nh(h,s)*Nr(h,s,n)*(r + ln(s)/2. + (Nr(h,s,n)-1)*ln(s)/2.)/1000.;
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
//Ti - The thickness of the insulation
inline double Pw(double h, double s, double n, double r, double C, double I, double Ti)
{
    return C*Ln(h, s, n, r)*I*I/((sqrt(s) - 2.*Ti)*(sqrt(s) - 2.*Ti));
}

//Calculate volume
//h - vertical size in mm
//R1 - internal radius in mm
//R2 - External radius
inline double Vol(double h, double R1, double R2)
{
    return PI*h*(R2*R2 - R1*R1)/1000.;
}

//Calculate vоlume outside
//R1 - internal radius in mm
//R2 - External radius
inline double Vou(double R1, double R2)
{
    return 4.*(PI*R1*R1*PI*R1*5./(4*6.) + (R2-R1)*PI*R1*R1/2.)/1000.;
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

//Magnetic lenght in m
//h - vertical size in mm
//R1 - internal radius in mm
//R2 - External radius im mm
inline double Lm(double h, double R1, double R2)
{
    return 2*h + R1*2*PI + 4*(R2-R1);
}

//Idle current
//h - vertical size in mm
//R1 - internal radius in mm
//R2 - External radius im mm
//U  - Input Voltage
//mu - magnetic permeability
//N - number of coils
//F - frequency
inline double Iid(double h, double R1, double R2, double U, double mu, double N, double F)
{
    return U*Lm(h, R1, R2)*1000./(2.*PI*PI*F*R1*R1*mu*MU*N*N);
}

//Idle current
//h - vertical size in mm
//R1 - internal radius in mm
//R2 - External radius im mm
//B  - Max Magnetic field
//mu - magnetic permeability
//N - number of coils
inline double Iid1(double h, double R1, double R2, double B, double mu, double N)
{
    return (B/sqrt(2.))*Lm(h, R1, R2)/(mu*MU*N*1000.);
}

//The magnetic core radius
//B - Max Magnetic field
//U  - Input Voltage
//N - number of coils
//F - frequency
inline double Rc(double B, double U, double N, double F)
{
    return sqrt(U*2./(sqrt(2.)*3.*PI*PI*F*B*N));
}

void trans(void)
{
    double Cal = 0.0282; //Conductivity aluminum OM*m/mm**2
    double Dal = 2.6989; //Density aluminum kg/gm**3
    double Dfe = 7.65; //Density of transformer steel kg/gm**3
    double Din = 1.5;  //Density of insulation kg/gm**3
    double R[5], U[2], M[5], N[2], I[4], s[2], P[2], L[2], S[2], Rz[2], B, T[2];
    double LM, mu, h, In, Rm, Sfc, Lc, PW, NN;
    double MT, PT, Mmin=100000, Pmin = 100000;
    struct save min;
    int i;

    PW = 100; //Power in kW
    U[0] = 315.;      //10000kV
    U[1] = 400./sqrt(3.); //400V
    I[0] = PW*1000./U[0]/3.; //the max current of 10Kv coil
    I[1] = PW*1000/U[1]/3.; //the max current of 400v coil
    mu = 20000.;     //Magnetic permeability of transformer steel
    B = 1.85;    //The max magnetic field
    Lc = 1; //Loss in magnetic core W/kg
    In = 2.;    //The thickness of the insulation
    Sfc = 0.955; //Core Stacking Factor
    T[0] = 0.06; //The thickness of the insulation of first coil
    T[1] = 0.06; //The thickness of the insulation of second coil
    NN = U[0]/U[1]/sqrt(3.); // Transformation coefficient

    /*
    U[0] = 240.;
    U[1] = 40.;
    I[0] = 6.25;
    I[1] = 37.5;
    R[0] = 24.;
    N[0] = 300.;
    N[1] = 49.;
    h = 2.14*14.;
    */
    i = 1;
    for(N[!i]=2; N[!i] <= 150; N[!i]++){
        for(h=100.; h < 400.; h+=5.){
            for(s[i]=50; s[i] < 200.; s[i]+=1){
            //for(s[0]=0.1; s[0] < 5.; s[0]+=0.1){
                for(s[!i]=50.; s[!i] < 200.; s[!i]+=1.){

                    N[i] = N[!i]*NN;

                    R[0] = 1000.*Rc(B, U[0], N[i], 50)*sqrt(2.-Sfc);
                    R[1] = R[0] + In;
                    R[2] = R[1]+Rr(h, s[0], N[0]);
                    R[3] = R[2] + In;
                    R[4] = R[3]+Rr(h, s[1], N[1]);

                    L[0] = Ln(h, s[0], N[0], R[1]);
                    L[1] = Ln(h, s[1], N[1], R[3]);

                    M[0] = Vol(h, 0, R[0])*Dfe/1000. + Vou(R[0],R[4])*Dfe/1000.;           //Mass steel
                    M[1] = Vol(h, R[0], R[1])*Din/1000.;        //Mass of insulation
                    M[2] = Vol(h, R[1], R[2])*Dal/1000.;        //Mass 10kV coil
                    M[3] = Vol(h, R[2], R[3])*Din/1000.;        //Mass of insulation
                    M[4] = Vol(h, R[3], R[4])*Dal/1000.;        //Mass 400V coil
                    MT = 3.*(M[0] + M[1] + M[2] + M[3] + M[4]);

                    P[0] = M[0]*Lc;                         //Loss in steel
                    P[2] = Pw(h, s[0], N[0], R[1], Cal, I[0], T[0]);  //Loss in 10kV coil
                    P[4] = Pw(h, s[1], N[1], R[3], Cal, I[1], T[1]);  //Loss in 400V coil
                    PT = 3.*(P[0] + P[2] + P[4]);

                    if(PT <= 1000.) {
                         if(MT < Mmin) {
                         //if(PT < Pmin){
                            Mmin = MT;
                            Pmin = PT;
                            min.h = h;
                            min.s[0] = s[0];
                            min.s[1] = s[1];
                            min.M[0] = M[0];
                            min.M[1] = M[1];
                            min.M[2] = M[2];
                            min.M[3] = M[3];
                            min.M[4] = M[4];
                            min.P[0] = P[0];
                            min.P[2] = P[2];
                            min.P[4] = P[4];
                            min.R[0] = R[0];
                            min.R[1] = R[1];
                            min.R[2] = R[2];
                            min.R[3] = R[3];
                            min.R[4] = R[4];
                            min.L[0] = L[0];
                            min.L[1] = L[1];
                            min.N[0] = N[0];
                            min.N[1] = N[1];
                            min.MT = MT;
                            min.PT = PT;
                        }
                    }
                //printf("h = %f s1 = %f s2 = %f PT = %f MT = %f\n",h, s[0], s[1], PT, MT);
                //printf("R1 = %f R2 = %f M0 = %f M2 = %f M4 = %f P0 = %f P2 = %f P4 = %f\n",
                //       R[1], R[2], M[0], M[2], M[4], P[0], P[2], P[4]);
                }
            }
        }
    }
    S[0] = 4.*min.L[0]*sqrt(min.s[0])/1000.;
    S[1] = 4.*min.L[1]*sqrt(min.s[1])/1000.;
    LM = Lm(min.h, min.R[0], min.R[4])/1000.;
    I[2] = Iid(min.h, min.R[0], min.R[4], U[0], mu, N[0], 50);
    I[3] = Iid1(min.h, min.R[0], min.R[4], B, mu, N[0]);
    Rz[0] = min.L[0]*Cal/min.s[0];
    Rz[1] = min.L[1]*Cal/min.s[1];
    //B = I[2]*mu*MU*N[0]/LM;
    Rm = 1000.*Rc(1.8, U[0], N[0], 50);

    printf("h = %f s1 = %f s2 = %f PT = %f MT = %f\n",min.h, min.s[0], min.s[1], min.PT, min.MT);
    printf("I0 = %f I1 = %f L1 = %f L2 = %f S1 = %f S2 = %f R1 = %f R2 = %f Lm = %f Iid = %f \n",
           I[0], I[1], min.L[0], min.L[1], S[0], S[1], Rz[0], Rz[1], LM, I[2]);
    printf("N0 = %f N1 = %f R0 = %f R1 = %f R2 = %f R3 = %f R4 = %f M0 = %f M1 = %f M2 = %f M3 = %f M4 = %f P0 = %f P2 = %f P4 = %f\n",
           min.N[0], min.N[1], min.R[0], min.R[1], min.R[2], min.R[3], min.R[4], min.M[0], min.M[1], min.M[2], min.M[3], min.M[4], min.P[0], min.P[2], min.P[4]);
    printf("Входные параметры\n");
    printf("Схема подключения треугольник - звезда\n");
    printf("Мощность %f кВ\n", PW);
    printf("Напряжение первичной обмотки %f В\n", U[0]);
    printf("Напряжение вторичной обмотки %f В\n", U[1]);
    printf("Максимальный ток первичной обмотки %f А\n", I[0]);
    printf("Максимальный ток вторичной обмотки %f А\n", I[1]);
    printf("Максимальное поле в магнитном сердечнике %f Тл\n", B);
    printf("Потери в магнитном сердечнике %f Вт/кг\n", Lc);
    printf("Толщина изоляции между катушками %f мм\n", In);
    printf("Коэффициент заполнения магнитопровода %f %%\n",Sfc*100.);
    printf("Толщина изоляции первичной обмотки %f мм\n",T[0]);
    printf("Толщина изоляции вторичной обмотки %f мм\n",T[1]);
    printf("Коэффициент трансформации обмоток %f \n", NN);

    printf("\nВыходные параметры\n");
    printf("Суммарные потери %f Вт\n", min.PT);
    printf("Масса %f кг\n", min.MT);
    printf("Ток холостого хода %f А %f A\n", I[2], I[3]);

    printf("\nМагнитный сердечник\n");
    printf("Масса  %f кг\n", min.M[0]);
    printf("Потери %f Вт\n", min.P[0]);
    printf("Радиус %f мм\n", min.R[0]);

    printf("\nПервичная катушка\n");
    printf("Количество витков %f\n", min.N[0]);
    printf("Масса  %f кг\n", min.M[2]);
    printf("Потери %f Вт\n", min.P[2]);
    printf("Сопротивление %f Ом\n", Rz[0]);
    printf("Длина %f м\n", min.L[0]);
    printf("Сечение %f мм**2\n", min.s[0]);
    printf("Площадь поверхности провода %f м**2\n", S[0]);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", min.R[1], min.R[2], min.h);
    printf("Число витков по высоте %f по радиусу %f \n", min.h/sqrt(min.s[0]), (min.R[2] - min.R[1])/sqrt(min.s[0]));

    printf("\nВторичная катушка\n");
    printf("Количество витков %f\n", min.N[1]);
    printf("Масса  %f кг\n", min.M[4]);
    printf("Потери %f Вт\n", min.P[4]);
    printf("Сопротивление %f Ом\n", Rz[1]);
    printf("Длина %f м\n", min.L[1]);
    printf("Сечение %f мм**2\n", min.s[1]);
    printf("Площадь поверхности провода %f м**2\n", S[1]);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", min.R[3], min.R[4], min.h);
    printf("Число витков по высоте %f по радиусу %f \n", min.h/sqrt(min.s[1]), (min.R[4] - min.R[3])/sqrt(min.s[1]));


}

