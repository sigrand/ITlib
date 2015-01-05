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

struct COIL {
    double U;   //Voltage V
    double I;   //Current A
    double T;   //The thickness of the wire insulation mm
    double In;  //The coil insulation mm
    double H;   //The height of the coil mm
    double R1;  //The internal radius mm
    double R2;  //The external radius mm
    double s;   //The cross-sectional area of the wire mm**2
    double L;   //The lenght of coil wire m
    double N;   //The number of turns
    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double R;   //Resistance
    double h;   //The wire height
    double w;   //The wire width
    double C;   //Conductivity
    double D;   //Density aluminum kg/gm**3
};

struct MCORE {
    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double D;   //Density of transformer steel
    double Sfc; //Core Stacking Factor
    double Lc;  //Loss in magnetic core W/kg
    double Mu;  //Magnetic permeability of transformer steel
    double B;   //The max magnetic field T
};

struct INS {
    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double D;   //Density of insulation kg/gm**3
};

struct TRANS {
    struct MCORE m;  //Magnetic core
    struct INS  i[2];  //Insulation
    struct COIL c[2];  //First coil

    double PW;  //Power in kW
    double NN;  //Transformation coefficient
    double M;   //Mass kg
    double V;   //Volume dm**3
    double S;   //square м**2
    double P;   //The loss power W
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

//Calculate lenght of coil in m
//Nh - number of turns in vertical
//Nr - number of turns in radial
//H - the height of the wire
//R1 - inside radius of coil
inline double Ln1(double Nh, double Nr, double H, double R1)
{
    return 2.*PI*Nh*Nr*(R1 + H/2. + (Nr-1)*H/2.)/1000.;
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

//Calculate loss power in coil in W
//h - vertical size in mm
//s - square of one wire in mm**2
//N - number of coils
//R1 - inside radius of coil
//R2 - outside radius of coil
//C - Conductivity
//I - Current
//Ti - The thickness of the insulation
inline double Pw1(double h, double s, double N, double R1, double R2, double C, double I, double Ti)
{
    double sd = ((sqrt(s) - 2.*Ti)*(sqrt(s) - 2.*Ti))/s;
    return PI*C*I*I*N*N*(R2 + R1)/(1000.*h*sd*(R2 - R1));
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
    double R[5], U[2], M[5], N[2], Nr[2], Nh[2], I[4], s[2], P[2], L[4], S[2], Rz[4], B, T[2], In[2], H[2], W[2], P1[2];
    double LM, mu, h, Rm, Sfc, Lc, PW, NN;
    double MT, PT, Mmin=100000, Pmin = 100000;
    struct save min;
    int i=1;

    PW = 100; //Power in kW
    U[i]  = 315.;      //10000kV
    U[!i] = 400./sqrt(3.); //400V
    I[i]  = PW*1000./U[i]/3.; //the max current of 10Kv coil
    I[!i] = PW*1000/U[!i]/3.; //the max current of 400v coil
    T[i]  = 0.06; //The thickness of the insulation of first coil
    T[!i] = 0.06; //The thickness of the insulation of second coil
    mu = 20000.;     //Magnetic permeability of transformer steel
    B = 1.85;    //The max magnetic field
    Lc = 1; //Loss in magnetic core W/kg
    In[0] = 0.;    //The thickness of the insulation
    In[1] = 2.;    //The thickness of the insulation
    Sfc = 0.955; //Core Stacking Factor
    NN = U[i]/U[!i]/sqrt(3.); // Transformation coefficient

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

    for(N[!i]=10; N[!i] <= 100; N[!i]++){
        for(h=200.; h < 300.; h+=1.){
            //{ s[i] = 90.;
            for(s[i]=50; s[i] < 200.; s[i]+=1){
            //for(s[i]=0.1; s[i] < 5.; s[i]+=0.1){
                //{ s[!i] = 90.;
                for(s[!i]=50.; s[!i] < 200.; s[!i]+=1.){

                    N[i] = N[!i]*NN;

                    R[0] = 1000.*Rc(B, U[i], N[i], 50)*sqrt(2.-Sfc);
                    R[1] = R[0] + In[0];
                    R[2] = R[1] + Rr(h, s[0], N[0]);
                    R[3] = R[2] + In[1];
                    R[4] = R[3] + Rr(h, s[1], N[1]);

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
                    //printf("h = %f s1 = %f s2 = %f N1 = %f N2 = %f PT = %f MT = %f\n",h, s[0], s[1], N[0], N[1], PT, MT);
                    //printf("N0 = %f N1 = %f R0 = %f R1 = %f R2 = %f R3 = %f R4 = %f M0 = %f M1 = %f M2 = %f M3 = %f M4 = %f P0 = %f P2 = %f P4 = %f\n",
                    //       N[0], N[1], R[0], R[1], R[2], R[3], R[4], M[0], M[1], M[2], M[3], M[4], P[0], P[2], P[4]);

                }
            }
        }
    }
    Nr[0] = ceil((min.R[2] - min.R[1])/sqrt(min.s[0]));
    Nh[0] = round(min.N[0]/Nr[0]);
    N[0] = Nr[0]*Nh[0];
    H[0] = (min.R[2] - min.R[1])/Nr[0];
    W[0] = min.h/Nh[0];
    L[2] = Ln1(Nh[0], Nr[0], H[0], min.R[1]);
    Rz[2] = L[2]*Cal/((H[0] - 2.*T[0])*(W[0] - 2.*T[0]));
    P1[0] = Rz[2]*I[0]*I[0];
    printf("Nr[0] = %f Nh[0] = %f N0 = %f H0 = %f W = %f L = %f R = %f P = %f\n", Nr[0], Nh[0], N[0], H[0], W[0], L[2], Rz[2], P1[0]);

    Nr[1] = ceil((min.R[4] - min.R[3])/sqrt(min.s[1]));
    Nh[1] = round(min.N[1]/Nr[1]);
    N[1] = Nr[1]*Nh[1];
    H[1] = (min.R[4] - min.R[3])/Nr[1];
    W[1] = min.h/Nh[1];
    L[3] = Ln1(Nh[1], Nr[1], H[1], min.R[3]);
    Rz[3] = L[3]*Cal/((H[1] - 2.*T[1])*(W[1] - 2.*T[1]));
    P1[1] = Rz[3]*I[1]*I[1];
    printf("Nr[1] = %f Nh[1] = %f N1 = %f H1 = %f W1 = %f L = %f R = %f P = %f\n", Nr[1], Nh[1], N[1], H[1], W[1], L[3], Rz[3], P1[1]);

    printf("Transformation = %f U = %f\n", N[i]/N[!i], U[i]*N[i]/N[!i]);

    S[0] = 4.*min.L[0]*sqrt(min.s[0])/1000.;
    S[1] = 4.*min.L[1]*sqrt(min.s[1])/1000.;
    LM = Lm(min.h, min.R[0], min.R[4])/1000.;
    I[2] = Iid(min.h, min.R[0], min.R[4], U[i], mu, N[i], 50);
    I[3] = Iid1(min.h, min.R[0], min.R[4], B, mu, N[i]);
    Rz[0] = min.L[0]*Cal/min.s[0];
    Rz[1] = min.L[1]*Cal/min.s[1];
    //B = I[2]*mu*MU*N[0]/LM;
    Rm = 1000.*Rc(1.8, U[i], N[0], 50);

    //printf("Power = %f\n", Pw1(min.h, min.s[0], min.N[0], min.R[1], min.R[2], Cal, I[0], T[0]));
    //printf("Power = %f\n", Pw1(min.h, min.s[1], min.N[1], min.R[3], min.R[4], Cal, I[1], T[1]));

    //printf("h = %f s1 = %f s2 = %f PT = %f MT = %f\n",min.h, min.s[0], min.s[1], min.PT, min.MT);
    //printf("I0 = %f I1 = %f L1 = %f L2 = %f S1 = %f S2 = %f R1 = %f R2 = %f Lm = %f Iid = %f \n",
    //       I[0], I[1], min.L[0], min.L[1], S[0], S[1], Rz[0], Rz[1], LM, I[2]);
    //printf("N0 = %f N1 = %f R0 = %f R1 = %f R2 = %f R3 = %f R4 = %f M0 = %f M1 = %f M2 = %f M3 = %f M4 = %f P0 = %f P2 = %f P4 = %f\n",
    //       min.N[0], min.N[1], min.R[0], min.R[1], min.R[2], min.R[3], min.R[4], min.M[0], min.M[1], min.M[2], min.M[3], min.M[4], min.P[0], min.P[2], min.P[4]);

    printf("Входные параметры\n");
    printf("Схема подключения треугольник - звезда\n");
    printf("Мощность %f кВ\n", PW);
    printf("Максимальное поле в магнитном сердечнике %f Тл\n", B);
    printf("Потери в магнитном сердечнике %f Вт/кг\n", Lc);
    printf("Толщина изоляции между катушками %f мм\n", In);
    printf("Коэффициент заполнения магнитопровода %f %%\n",Sfc*100.);
    printf("Коэффициент трансформации обмоток %f \n", NN);

    printf("\nВыходные параметры\n");
    printf("Суммарные потери %f Вт\n", min.PT);
    printf("Масса %f кг\n", min.MT);
    printf("Ток холостого хода %f А %f A\n", I[2], I[3]);

    printf("\nМагнитный сердечник\n");
    printf("Масса  %f кг\n", min.M[0]);
    printf("Потери %f Вт\n", min.P[0]);
    printf("Радиус %f мм\n", min.R[0]);

    printf("\nПервая катушка\n");
    printf("Напряжение %f В\n", U[0]);
    printf("Максимальный ток %f А\n", I[0]);
    printf("Толщина покрытия  %f мм\n",T[0]);
    printf("Толщина изоляции между магнитопроводом и обмоткой %f мм\n",In[0]);
    printf("Количество витков %f\n", min.N[0]);
    printf("Масса  %f кг\n", min.M[2]);
    printf("Потери %f Вт\n", min.P[2]);
    printf("Сопротивление %f Ом\n", Rz[0]);
    printf("Длина %f м\n", min.L[0]);
    printf("Сечение %f мм**2\n", min.s[0]);
    printf("Площадь поверхности провода %f м**2\n", S[0]);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", min.R[1], min.R[2], min.h);
    printf("Число витков по высоте %f по радиусу %f \n", min.h/sqrt(min.s[0]), (min.R[2] - min.R[1])/sqrt(min.s[0]));

    printf("\nВторая катушка\n");
    printf("Напряжение %f В\n", U[1]);
    printf("Максимальный ток %f А\n", I[1]);
    printf("Толщина покрытия %f мм\n",T[1]);
    printf("Толщина изоляции между обмотками                  %f мм\n",In[1]);
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

void trans1(void)
{
    struct TRANS t;

    //Input parameters
    t.PW = 100;     //Power in kW

    //Magnetic core
    t.m.D = 7.65;                   //Density of transformer steel kg/gm**3
    t.m.Mu = 20000.;                //Magnetic permeability of transformer steel
    t.m.B = 1.85;                   //The max magnetic field
    t.m.Lc = 1;                     //Loss in magnetic core W/kg
    t.m.Sfc = 0.955;                //Core Stacking Factor

    //First insulator
    t.i[0].D = 1.5;                   //Density of insulation kg/gm**3

    //First coil
    t.c[0].U = 315.;                   //10000kV
    t.c[0].I = t.PW*1000./t.c[0].U/3.;   //the max current of 10Kv coil
    t.c[0].T  = 0.06;                  //The thickness of the insulation of first coil
    t.c[0].In = 0.;                    //The thickness of the insulation
    t.c[0].C  = 0.0282;                //Conductivity of aluminum OM*m/mm**2
    t.c[0].D  = 2.6989;                //Density aluminum kg/gm**3

    //Secons insulator
    t.i[1].D = 1.5;                   //Density of insulation kg/gm**3

    //Second coil
    t.c[1].U = 400./sqrt(3.);          //400V
    t.c[1].I = t.PW*1000./t.c[1].U/3.;   //the max current of 10Kv coil
    t.c[1].T  = 0.06;                  //The thickness of the insulation of first coil
    t.c[1].In = 2.;                    //The thickness of the insulation
    t.c[1].C  = 0.0282;                //Conductivity of aluminum OM*m/mm**2
    t.c[1].D  = 2.6989;                //Density aluminum kg/gm**3

    t.NN = t.c[0].U/t.c[1].U/sqrt(3.);        // Transformation coefficient

    double R[5], U[2], M[5], N[2], Nr[2], Nh[2], I[4], s[2], P[2], L[4], S[2], Rz[4], B, T[2], In[2], H[2], W[2], P1[2];
    double LM, mu, h, Rm, Sfc, Lc, PW, NN;
    double MT, PT, Mmin=100000, Pmin = 100000;
    struct save min;
    int i=1;

    PW = 100; //Power in kW
    U[i]  = 315.;      //10000kV
    U[!i] = 400./sqrt(3.); //400V
    I[i]  = PW*1000./U[i]/3.; //the max current of 10Kv coil
    I[!i] = PW*1000/U[!i]/3.; //the max current of 400v coil
    T[i]  = 0.06; //The thickness of the insulation of first coil
    T[!i] = 0.06; //The thickness of the insulation of second coil
    mu = 20000.;     //Magnetic permeability of transformer steel
    B = 1.85;    //The max magnetic field
    Lc = 1; //Loss in magnetic core W/kg
    In[0] = 0.;    //The thickness of the insulation
    In[1] = 2.;    //The thickness of the insulation
    Sfc = 0.955; //Core Stacking Factor
    NN = U[i]/U[!i]/sqrt(3.); // Transformation coefficient

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

    for(N[!i]=10; N[!i] <= 100; N[!i]++){
        for(h=200.; h < 300.; h+=1.){
            //{ s[i] = 90.;
            for(s[i]=50; s[i] < 200.; s[i]+=1){
            //for(s[i]=0.1; s[i] < 5.; s[i]+=0.1){
                //{ s[!i] = 90.;
                for(s[!i]=50.; s[!i] < 200.; s[!i]+=1.){

                    N[i] = N[!i]*NN;

                    R[0] = 1000.*Rc(B, U[i], N[i], 50)*sqrt(2.-Sfc);
                    R[1] = R[0] + In[0];
                    R[2] = R[1] + Rr(h, s[0], N[0]);
                    R[3] = R[2] + In[1];
                    R[4] = R[3] + Rr(h, s[1], N[1]);

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
                    //printf("h = %f s1 = %f s2 = %f N1 = %f N2 = %f PT = %f MT = %f\n",h, s[0], s[1], N[0], N[1], PT, MT);
                    //printf("N0 = %f N1 = %f R0 = %f R1 = %f R2 = %f R3 = %f R4 = %f M0 = %f M1 = %f M2 = %f M3 = %f M4 = %f P0 = %f P2 = %f P4 = %f\n",
                    //       N[0], N[1], R[0], R[1], R[2], R[3], R[4], M[0], M[1], M[2], M[3], M[4], P[0], P[2], P[4]);

                }
            }
        }
    }
    Nr[0] = ceil((min.R[2] - min.R[1])/sqrt(min.s[0]));
    Nh[0] = round(min.N[0]/Nr[0]);
    N[0] = Nr[0]*Nh[0];
    H[0] = (min.R[2] - min.R[1])/Nr[0];
    W[0] = min.h/Nh[0];
    L[2] = Ln1(Nh[0], Nr[0], H[0], min.R[1]);
    Rz[2] = L[2]*Cal/((H[0] - 2.*T[0])*(W[0] - 2.*T[0]));
    P1[0] = Rz[2]*I[0]*I[0];
    printf("Nr[0] = %f Nh[0] = %f N0 = %f H0 = %f W = %f L = %f R = %f P = %f\n", Nr[0], Nh[0], N[0], H[0], W[0], L[2], Rz[2], P1[0]);

    Nr[1] = ceil((min.R[4] - min.R[3])/sqrt(min.s[1]));
    Nh[1] = round(min.N[1]/Nr[1]);
    N[1] = Nr[1]*Nh[1];
    H[1] = (min.R[4] - min.R[3])/Nr[1];
    W[1] = min.h/Nh[1];
    L[3] = Ln1(Nh[1], Nr[1], H[1], min.R[3]);
    Rz[3] = L[3]*Cal/((H[1] - 2.*T[1])*(W[1] - 2.*T[1]));
    P1[1] = Rz[3]*I[1]*I[1];
    printf("Nr[1] = %f Nh[1] = %f N1 = %f H1 = %f W1 = %f L = %f R = %f P = %f\n", Nr[1], Nh[1], N[1], H[1], W[1], L[3], Rz[3], P1[1]);

    printf("Transformation = %f U = %f\n", N[i]/N[!i], U[i]*N[i]/N[!i]);

    S[0] = 4.*min.L[0]*sqrt(min.s[0])/1000.;
    S[1] = 4.*min.L[1]*sqrt(min.s[1])/1000.;
    LM = Lm(min.h, min.R[0], min.R[4])/1000.;
    I[2] = Iid(min.h, min.R[0], min.R[4], U[i], mu, N[i], 50);
    I[3] = Iid1(min.h, min.R[0], min.R[4], B, mu, N[i]);
    Rz[0] = min.L[0]*Cal/min.s[0];
    Rz[1] = min.L[1]*Cal/min.s[1];
    //B = I[2]*mu*MU*N[0]/LM;
    Rm = 1000.*Rc(1.8, U[i], N[0], 50);

    //printf("Power = %f\n", Pw1(min.h, min.s[0], min.N[0], min.R[1], min.R[2], Cal, I[0], T[0]));
    //printf("Power = %f\n", Pw1(min.h, min.s[1], min.N[1], min.R[3], min.R[4], Cal, I[1], T[1]));

    //printf("h = %f s1 = %f s2 = %f PT = %f MT = %f\n",min.h, min.s[0], min.s[1], min.PT, min.MT);
    //printf("I0 = %f I1 = %f L1 = %f L2 = %f S1 = %f S2 = %f R1 = %f R2 = %f Lm = %f Iid = %f \n",
    //       I[0], I[1], min.L[0], min.L[1], S[0], S[1], Rz[0], Rz[1], LM, I[2]);
    //printf("N0 = %f N1 = %f R0 = %f R1 = %f R2 = %f R3 = %f R4 = %f M0 = %f M1 = %f M2 = %f M3 = %f M4 = %f P0 = %f P2 = %f P4 = %f\n",
    //       min.N[0], min.N[1], min.R[0], min.R[1], min.R[2], min.R[3], min.R[4], min.M[0], min.M[1], min.M[2], min.M[3], min.M[4], min.P[0], min.P[2], min.P[4]);

    printf("Входные параметры\n");
    printf("Схема подключения треугольник - звезда\n");
    printf("Мощность %f кВ\n", PW);
    printf("Максимальное поле в магнитном сердечнике %f Тл\n", B);
    printf("Потери в магнитном сердечнике %f Вт/кг\n", Lc);
    printf("Толщина изоляции между катушками %f мм\n", In);
    printf("Коэффициент заполнения магнитопровода %f %%\n",Sfc*100.);
    printf("Коэффициент трансформации обмоток %f \n", NN);

    printf("\nВыходные параметры\n");
    printf("Суммарные потери %f Вт\n", min.PT);
    printf("Масса %f кг\n", min.MT);
    printf("Ток холостого хода %f А %f A\n", I[2], I[3]);

    printf("\nМагнитный сердечник\n");
    printf("Масса  %f кг\n", min.M[0]);
    printf("Потери %f Вт\n", min.P[0]);
    printf("Радиус %f мм\n", min.R[0]);

    printf("\nПервая катушка\n");
    printf("Напряжение %f В\n", U[0]);
    printf("Максимальный ток %f А\n", I[0]);
    printf("Толщина покрытия  %f мм\n",T[0]);
    printf("Толщина изоляции между магнитопроводом и обмоткой %f мм\n",In[0]);
    printf("Количество витков %f\n", min.N[0]);
    printf("Масса  %f кг\n", min.M[2]);
    printf("Потери %f Вт\n", min.P[2]);
    printf("Сопротивление %f Ом\n", Rz[0]);
    printf("Длина %f м\n", min.L[0]);
    printf("Сечение %f мм**2\n", min.s[0]);
    printf("Площадь поверхности провода %f м**2\n", S[0]);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", min.R[1], min.R[2], min.h);
    printf("Число витков по высоте %f по радиусу %f \n", min.h/sqrt(min.s[0]), (min.R[2] - min.R[1])/sqrt(min.s[0]));

    printf("\nВторая катушка\n");
    printf("Напряжение %f В\n", U[1]);
    printf("Максимальный ток %f А\n", I[1]);
    printf("Толщина покрытия %f мм\n",T[1]);
    printf("Толщина изоляции между обмотками                  %f мм\n",In[1]);
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

