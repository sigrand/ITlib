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
    double R;   //The external radius mm
    double s;   //The cross-sectional area of the wire mm**2
    double L;   //The lenght of coil wire m
    double N;   //The number of turns
    double Nr;   //The number of turns in radial derection
    double Nh;   //The number of turns vertical direction
    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double Rz;   //Resistance
    double h;   //The wire height
    double w;   //The wire width
    double C;   //Conductivity
    double D;   //Density aluminum kg/gm**3
    double S;   //The surface area of the wire m**2
    double sp[3]; //The cross-sectional area of the wire range
};

struct MCORE {
    double R;   //The external radius mm
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
    double R;   //The external radius mm
    double M;   //Mass kg
    double V;   //Volume dm**3
    double P;   //The loss power W
    double D;   //Density of insulation kg/gm**3
    double T;   //The thickness of the insulation mm
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
    double H;   //The height of the coil, insulation and magnetic core mm
    double I;   //No-load current
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

void trans1(void)
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

void trans_3_phase(struct TRANS *t, double PT, int p, double N1, double N2, double H1, double H2)
{
    struct TRANS tm;
    double Mmin=1000000, Pmin = 1000000;
    int i, j=0;


    if(t->c[0].U > t->c[1].U)   { t->NN = t->c[0].U/t->c[1].U; i = 1; }
    else                        { t->NN = t->c[1].U/t->c[0].U; i = 0; }

    for(t->c[i].N = N1; t->c[i].N <= N2; t->c[i].N++){
        for(t->H = H1; t->H <= H2; t->H++){
            for(t->c[0].s = t->c[0].sp[0]; t->c[0].s <= t->c[0].sp[1]; t->c[0].s+=t->c[0].sp[2]){
                for(t->c[1].s = t->c[1].sp[0]; t->c[1].s <= t->c[1].sp[1]; t->c[1].s+=t->c[1].sp[2]){
                    //printf("U0 = %f U1 = %f i = %d N = %f H = %f s1 = %f s2 = %f\n", t->c[0].U, t->c[1].U, i, t->c[i].N, t->H, t->c[0].s, t->c[1].s);

                    //t->c[1].s = t->c[0].s;
                    t->c[!i].N = t->c[i].N*t->NN;

                    t->m.R = 1000.*Rc(t->m.B, t->c[p].U, t->c[p].N, 50)*sqrt(2.-t->m.Sfc);
                    t->i[0].R = t->m.R + t->i[0].T;
                    t->c[0].R = t->i[0].R + Rr(t->H, t->c[0].s, t->c[0].N);
                    t->i[1].R = t->c[0].R + t->i[1].T;
                    t->c[1].R = t->i[1].R + Rr(t->H, t->c[1].s, t->c[1].N);

                    t->c[0].L = Ln(t->H, t->c[0].s, t->c[0].N, t->i[0].R);
                    t->c[1].L = Ln(t->H, t->c[1].s, t->c[1].N, t->i[1].R);

                    t->m.V    = Vol(t->H, 0, t->m.R) + Vou(t->m.R, t->c[1].R);
                    t->i[0].V = Vol(t->H, t->m.R, t->i[0].R);
                    t->c[0].V = Vol(t->H, t->i[0].R, t->c[0].R);
                    t->i[1].V = Vol(t->H, t->c[0].R, t->i[1].R);
                    t->c[1].V = Vol(t->H, t->i[1].R, t->c[1].R);

                    t->m.M    = t->m.V*t->m.D/1000;        //Mass steel
                    t->i[0].M = t->i[0].V*t->i[0].D/1000.; //Mass of insulation
                    t->c[0].M = t->c[0].V*t->c[0].D/1000.; //Mass coil
                    t->i[1].M = t->i[1].V*t->i[1].D/1000.; //Mass of insulation
                    t->c[1].M = t->c[1].V*t->c[1].D/1000.; //Mass 1coil

                    t->M = 3.*(t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);

                    t->m.P = t->m.M*t->m.Lc;                         //Loss in steel
                    t->c[0].P = Pw(t->H, t->c[0].s, t->c[0].N, t->i[0].R, t->c[0].C, t->c[0].I, t->c[0].T);
                    t->c[1].P = Pw(t->H, t->c[1].s, t->c[1].N, t->i[1].R, t->c[1].C, t->c[1].I, t->c[1].T);

                    t->P = 3.*(t->m.P + t->c[0].P + t->c[1].P);

                    if(t->P <= PT) {
                         if(t->M < Mmin) {
                         //if(t->P < Pmin){
                            Mmin = t->M;
                            Pmin = t->P;
                            memcpy (&tm, t, sizeof(tm));
                            j++;
                        }
                    }


                }
            }
        }
    }

    if(j) memcpy (t, &tm, sizeof(tm));

}

void trans(void)
{
    struct TRANS t;
    int p = 1;                              //The primary coil

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
    t.i[0].T = 0.;                    //The thickness of the insulation mm

    //First coil
    t.c[p].U = 10000.;                   //10000kV
    t.c[p].I = t.PW*1000./t.c[p].U/3.;   //the max current of 10Kv coil
    t.c[p].T  = 0.06;                  //The thickness of the insulation of first coil
    t.c[p].C  = 0.0282;                //Conductivity of aluminum OM*m/mm**2
    t.c[p].D  = 2.6989;                //Density aluminum kg/gm**3
    t.c[p].sp[0] = 0.1;                //From
    t.c[p].sp[1] = 6;                  //To
    t.c[p].sp[2] = 0.1;                //Spep

    //Secons insulator
    t.i[1].D = 1.5;                   //Density of insulation kg/gm**3
    t.i[1].T = 2.;                    //The thickness of the insulation mm

    //Second coil
    t.c[!p].U = 400.;                    //400V
    t.c[!p].I = t.PW*1000.*sqrt(3.)/t.c[!p].U/3.;   //the max current of 10Kv coil
    t.c[!p].T  = 0.06;                  //The thickness of the insulation of first coil
    t.c[!p].C  = 0.0282;                //Conductivity of aluminum OM*m/mm**2
    t.c[!p].D  = 2.6989;                //Density aluminum kg/gm**3
    t.c[!p].sp[0] = 50;                  //From
    t.c[!p].sp[1] = 200;                 //To
    t.c[!p].sp[2] = 1;                   //Spep

    //t.NN = t.c[0].U/t.c[1].U;        // Transformation coefficient

    double N[2], H[2], s[2][2];

    N[0] = 10; N[1] = 100;
    H[0] = 200; H[1] = 300;
    //s[0][0] = 0.1; s[0][1] = 6;
    s[0][0] = 50; s[0][1] = 200;
    //s[1][0] = 50; s[1][1] = 200;
    s[1][0] = 0.1; s[1][1] = 6;

    trans_3_phase(&t, 1000, p, N[0], N[1], H[0], H[1]);

    t.c[0].S = 4.*t.c[0].L*sqrt(t.c[0].s)/1000.;
    t.c[1].S = 4.*t.c[1].L*sqrt(t.c[1].s)/1000.;
    //LM = Lm(min.h, min.R[0], min.R[4])/1000.;
    //t.I = Iid(t.H, min.R[0], min.R[4], U[i], mu, N[i], 50);
    t.I = Iid1(t.H, t.i[p].R, t.c[p].R, t.m.B, t.m.Mu, t.c[p].N);
    t.c[0].Rz = t.c[0].L*t.c[0].C/t.c[0].s;
    t.c[1].Rz = t.c[1].L*t.c[1].C/t.c[1].s;
    //B = I[2]*mu*MU*N[0]/LM;
    //Rm = 1000.*Rc(1.8, U[i], N[0], 50);

    printf("Входные параметры\n");
    printf("Схема подключения треугольник - звезда\n");
    printf("Мощность %f кВ\n", t.PW);
    printf("Максимальное поле в магнитном сердечнике %f Тл\n", t.m.B);
    printf("Потери в магнитном сердечнике %f Вт/кг\n", t.m.Lc);
    printf("Коэффициент заполнения магнитопровода %f %%\n",t.m.Sfc*100.);
    printf("Коэффициент трансформации обмоток %f \n", t.NN);

    printf("\nВыходные параметры\n");
    printf("Суммарные потери %f Вт\n", t.P);
    printf("Масса %f кг\n", t.M);
    printf("Ток холостого хода %f А\n", t.I);

    printf("\nМагнитный сердечник\n");
    printf("Масса  %f кг\n", t.m.M);
    printf("Потери %f Вт\n", t.m.P);
    printf("Радиус %f мм\n", t.m.R);

    printf("\nПервая катушка\n");
    printf("Напряжение %f В\n", t.c[0].U);
    printf("Максимальный ток %f А\n", t.c[0].I);
    printf("Толщина покрытия  %f мм\n",t.c[0].T);
    printf("Толщина изоляции между магнитопроводом и обмоткой %f мм\n",t.i[0].T);
    printf("Количество витков %f\n", t.c[0].N);
    printf("Масса  %f кг\n", t.c[0].M);
    printf("Потери %f Вт\n", t.c[0].P);
    printf("Сопротивление %f Ом\n", t.c[0].Rz);
    printf("Длина %f м\n", t.c[0].L);
    printf("Сечение %f мм**2\n", t.c[0].s);
    printf("Площадь поверхности провода %f м**2\n", t.c[0].S);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", t.i[0].R, t.c[0].R, t.H);
    printf("Число витков по высоте %f по радиусу %f \n", t.H/sqrt(t.c[0].s), (t.c[0].R - t.i[0].R)/sqrt(t.c[0].s));

    printf("\nВторая катушка\n");
    printf("Напряжение %f В\n", t.c[1].U);
    printf("Максимальный ток %f А\n", t.c[1].I);
    printf("Толщина покрытия  %f мм\n",t.c[1].T);
    printf("Толщина изоляции между обмотками                  %f мм\n",t.i[1].T);
    printf("Количество витков %f\n", t.c[1].N);
    printf("Масса  %f кг\n", t.c[1].M);
    printf("Потери %f Вт\n", t.c[1].P);
    printf("Сопротивление %f Ом\n", t.c[1].Rz);
    printf("Длина %f м\n", t.c[1].L);
    printf("Сечение %f мм**2\n", t.c[1].s);
    printf("Площадь поверхности провода %f м**2\n", t.c[1].S);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", t.i[1].R, t.c[1].R, t.H);
    printf("Число витков по высоте %f по радиусу %f \n", t.H/sqrt(t.c[1].s), (t.c[1].R - t.i[1].R)/sqrt(t.c[1].s));

    t.c[0].Nr = ceil((t.c[0].R - t.i[0].R)/sqrt(t.c[0].s));
    t.c[0].Nh = round(t.c[0].N/t.c[0].Nr);
    t.c[0].N = t.c[0].Nr*t.c[0].Nh;
    t.c[0].h = (t.c[0].R - t.i[0].R)/t.c[0].Nr;
    t.c[0].w = t.H/t.c[0].Nh;
    t.c[0].L = Ln1(t.c[0].Nh, t.c[0].Nr, t.c[0].h, t.i[0].R);
    t.c[0].Rz = t.c[0].L*t.c[0].C/((t.c[0].h - 2.*t.c[0].T)*(t.c[0].w - 2.*t.c[0].T));
    t.c[0].P = t.c[0].Rz*t.c[0].I*t.c[0].I;
    printf("Nr = %f Nh = %f N = %f h = %f w = %f L = %f Rz = %f P = %f\n",
           t.c[0].Nr, t.c[0].Nh, t.c[0].N, t.c[0].h, t.c[0].w, t.c[0].L, t.c[0].Rz, t.c[0].P);

    t.c[1].Nr = ceil((t.c[1].R - t.i[1].R)/sqrt(t.c[1].s));
    t.c[1].Nh = round(t.c[1].N/t.c[1].Nr);
    t.c[1].N = t.c[1].Nr*t.c[1].Nh;
    t.c[1].h = (t.c[1].R - t.i[1].R)/t.c[1].Nr;
    t.c[1].w = t.H/t.c[1].Nh;
    t.c[1].L = Ln1(t.c[1].Nh, t.c[1].Nr, t.c[1].h, t.i[1].R);
    t.c[1].Rz = t.c[1].L*t.c[1].C/((t.c[1].h - 2.*t.c[1].T)*(t.c[1].w - 2.*t.c[1].T));
    t.c[1].P = t.c[1].Rz*t.c[1].I*t.c[1].I;
    printf("Nr = %f Nh = %f N = %f h = %f w = %f L = %f Rz = %f P = %f\n",
           t.c[1].Nr, t.c[1].Nh, t.c[1].N, t.c[1].h, t.c[1].w, t.c[1].L, t.c[1].Rz, t.c[1].P);

    t.M = 3.*(t.m.M + t.i[0].M + t.c[0].M + t.i[1].M + t.c[1].M);
    t.P = 3.*(t.m.P + t.c[0].P + t.c[1].P);

    if(t.c[0].U > t.c[1].U)   t.NN = t.c[0].N/t.c[1].N;
    else                      t.NN = t.c[1].N/t.c[0].N;

    if(t.c[p].U < t.c[!p].U)  t.c[!p].U = t.NN*t.c[p].U;
    else                      t.c[!p].U = t.c[p].U/t.NN;


    //printf("Transformation = %f U = %f\n", t.NN, t.c[!p].U);

    printf("Входные параметры\n");
    printf("Схема подключения треугольник - звезда\n");
    printf("Мощность %f кВ\n", t.PW);
    printf("Максимальное поле в магнитном сердечнике %f Тл\n", t.m.B);
    printf("Потери в магнитном сердечнике %f Вт/кг\n", t.m.Lc);
    printf("Коэффициент заполнения магнитопровода %f %%\n",t.m.Sfc*100.);
    printf("Коэффициент трансформации обмоток %f \n", t.NN);

    printf("\nВыходные параметры\n");
    printf("Суммарные потери %f Вт\n", t.P);
    printf("Масса %f кг\n", t.M);
    printf("Ток холостого хода %f А\n", t.I);

    printf("\nМагнитный сердечник\n");
    printf("Масса  %f кг\n", t.m.M);
    printf("Потери %f Вт\n", t.m.P);
    printf("Радиус %f мм\n", t.m.R);

    printf("\nПервая катушка\n");
    printf("Напряжение %f В\n", t.c[0].U);
    printf("Максимальный ток %f А\n", t.c[0].I);
    printf("Толщина покрытия  %f мм\n",t.c[0].T);
    printf("Толщина изоляции между магнитопроводом и обмоткой %f мм\n",t.i[0].T);
    printf("Количество витков %f\n", t.c[0].N);
    printf("Масса  %f кг\n", t.c[0].M);
    printf("Потери %f Вт\n", t.c[0].P);
    printf("Сопротивление %f Ом\n", t.c[0].Rz);
    printf("Длина %f м\n", t.c[0].L);
    printf("Сечение %f мм**2\n", t.c[0].s);
    printf("Площадь поверхности провода %f м**2\n", t.c[0].S);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", t.i[0].R, t.c[0].R, t.H);
    printf("Число витков по высоте %f по радиусу %f \n", t.c[0].Nh, t.c[0].Nr);

    printf("\nВторая катушка\n");
    printf("Напряжение %f В\n", t.c[1].U);
    printf("Максимальный ток %f А\n", t.c[1].I);
    printf("Толщина покрытия  %f мм\n",t.c[1].T);
    printf("Толщина изоляции между обмотками                  %f мм\n",t.i[1].T);
    printf("Количество витков %f\n", t.c[1].N);
    printf("Масса  %f кг\n", t.c[1].M);
    printf("Потери %f Вт\n", t.c[1].P);
    printf("Сопротивление %f Ом\n", t.c[1].Rz);
    printf("Длина %f м\n", t.c[1].L);
    printf("Сечение %f мм**2\n", t.c[1].s);
    printf("Площадь поверхности провода %f м**2\n", t.c[1].S);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", t.i[1].R, t.c[1].R, t.H);
    printf("Число витков по высоте %f по радиусу %f \n", t.c[1].Nh, t.c[1].Nr);



}

