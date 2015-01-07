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
    double sp[2]; //The cross-sectional area of the wire range
    double Np[2]; //The number of turns range and step
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
    double Hp[2];   //The height of the coilrange mm
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
    return sqrt(U/(sqrt(2.)*PI*PI*F*B*N));
}

//The first coil turns
//B - Max Magnetic field
//U  - Input Voltage
//N - number of coils
//F - frequency
inline double Nc(double B, double U, double R, double F)
{
    return U*1000000./(sqrt(2.)*PI*PI*F*B*R*R);
}

int trans_3_phase(struct TRANS *t, double PT, int p)
{
    struct TRANS tm;
    double Mmin=1000000, Pmin = 1000000, s, s1, H, N[2];
    int i, j=0, k;

    k = !t->m.R ? 1 : 0;
    if(t->c[0].U > t->c[1].U)   { t->NN = t->c[0].U/t->c[1].U; i = 1; }
    else                        { t->NN = t->c[1].U/t->c[0].U; i = 0; }

    if(t->m.R) {
        k = 0;
        N[p] = Nc(t->m.B, t->c[p].U, t->m.R*sqrt(t->m.Sfc), 50);
        if(p == i) N[!i] = N[i]*t->NN;
        else       N[i] = N[!i]/t->NN;
        t->c[i].Np[0] = N[i];
        t->c[i].N = N[i];
    } else k = 1;

    printf("i = %d N0 = %f N1 = %f H = %f Hp1 = %f Hp2 = %f  Np1 = %f  Np2 = %f \n", i, N[0], N[1], t->H, t->Hp[0], t->Hp[1], t->c[i].Np[0], t->c[i].Np[1]);
    for(N[i] = t->c[i].N; N[i] <= t->c[i].Np[0]; N[i] += t->c[i].Np[1]){
        for(H = t->H; H <= t->Hp[0]; H += t->Hp[1]){
            for(s = t->c[0].s; s <= t->c[0].sp[0]; s += t->c[0].sp[1]){
                {//for(s1 = t->c[1].s ; s1 <= t->c[1].sp[0]; s1 += t->c[1].sp[1]){
                    //printf("N = %f H = %f s = %f s1 = %f \n", N[i], H, s, s1);
                    //printf("U0 = %f U1 = %f i = %d N = %f H = %f s1 = %f s2 = %f\n", t->c[0].U, t->c[1].U, i, t->c[i].N, t->H, t->c[0].s, t->c[1].s);
                    s1 = s;
                    if(k){
                        N[!i] = N[i]*t->NN;
                        t->m.R = 1000.*Rc(t->m.B, t->c[p].U, N[p], 50)/sqrt(t->m.Sfc);
                    }

                    t->i[0].R = t->m.R + t->i[0].T;
                    t->c[0].R = t->i[0].R + Rr(H, s , N[0]);
                    t->i[1].R = t->c[0].R + t->i[1].T;
                    t->c[1].R = t->i[1].R + Rr(H, s1, N[1]);

                    t->c[0].L = Ln(H, s , N[0], t->i[0].R);
                    t->c[1].L = Ln(H, s1, N[1], t->i[1].R);

                    t->m.V    = Vol(H, 0, t->m.R) + Vou(t->m.R, t->c[1].R);
                    t->i[0].V = Vol(H, t->m.R, t->i[0].R);
                    t->c[0].V = Vol(H, t->i[0].R, t->c[0].R);
                    t->i[1].V = Vol(H, t->c[0].R, t->i[1].R);
                    t->c[1].V = Vol(H, t->i[1].R, t->c[1].R);

                    t->m.M    = t->m.V*t->m.D/1000;        //Mass steel
                    t->i[0].M = t->i[0].V*t->i[0].D/1000.; //Mass of insulation
                    t->c[0].M = t->c[0].V*t->c[0].D/1000.; //Mass coil
                    t->i[1].M = t->i[1].V*t->i[1].D/1000.; //Mass of insulation
                    t->c[1].M = t->c[1].V*t->c[1].D/1000.; //Mass 1coil

                    t->M = 3.*(t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);

                    t->m.P = t->m.M*t->m.Lc;                         //Loss in steel
                    t->c[0].P = Pw(H, s , N[0], t->i[0].R, t->c[0].C, t->c[0].I, t->c[0].T);
                    t->c[1].P = Pw(H, s1, N[1], t->i[1].R, t->c[1].C, t->c[1].I, t->c[1].T);

                    t->P = 3.*(t->m.P + t->c[0].P + t->c[1].P);

                    if(t->P <= PT) {
                         if(t->M < Mmin) {
                         //if(t->P < Pmin){
                            Mmin = t->M;
                            Pmin = t->P;
                            memcpy (&tm, t, sizeof(tm));
                            tm.c[0].s  = s ;
                            tm.c[1].s  = s1;
                            tm.c[0].N = N[0];
                            tm.c[1].N = N[1];
                            tm.H = H;
                            j++;
                        }
                    }
                }
            }
        }
    }

    if(j) {
        memcpy (t, &tm, sizeof(tm));
        return 0;
    } else return 1;
}

void trans(void)
{
    struct TRANS t;
    int p = 1;                              //The primary coil

    //Input parameters
    t.PW = 100;     //Power in kW
    t.H = 200;
    t.Hp[0] = 400;                  //Max value
    t.Hp[1] = 1;

    //Magnetic core
    t.m.D = 7.65;                    //Density of transformer steel kg/gm**3
    t.m.Mu = 20000.;                 //Magnetic permeability of transformer steel
    t.m.B = 1.85;                    //The max magnetic field
    t.m.Lc = 1.;                     //Loss in magnetic core W/kg
    t.m.Sfc = 0.955;                 //Core Stacking Factor
    t.m.R = 70;                       //Core Radius

    //First insulator
    t.i[0].D = 1.5;                   //Density of insulation kg/gm**3
    t.i[0].T = 0.;                    //The thickness of the insulation mm

    //First coil
    t.c[p].T  = 0.06;                  //The thickness of the insulation of first coil
    t.c[p].C  = 0.0282;                //Conductivity of aluminum OM*m/mm**2
    t.c[p].D  = 2.6989;                //Density aluminum kg/gm**3

    t.c[p].U = 315.;                   //10000kV
    t.c[p].I = t.PW*1000./t.c[p].U/3.;   //the max current of 10Kv coil
    t.c[p].s    = 50;//0.1;                //From
    t.c[p].sp[0] = 200; //6;                  //To
    t.c[p].sp[1] = 1; //0.1;                //Spep
    t.c[p].N     = 30;                  //From
    t.c[p].Np[0] = 100;                 //To
    t.c[p].Np[1] = 1;                   //Spep
    /*
    t.c[p].U = 10000.;                   //10000kV
    t.c[p].I = t.PW*1000./t.c[p].U/3.;   //the max current of 10Kv coil
    t.c[p].s    = 0.1;                //From
    t.c[p].sp[0] = 6;                  //To
    t.c[p].sp[1] = 0.1;                //Spep
    */
    //Secons insulator
    t.i[1].D = 1.5;                   //Density of insulation kg/gm**3
    t.i[1].T = 2.;                    //The thickness of the insulation mm

    //Second coil
    t.c[!p].U = 400;                    //400V
    t.c[!p].I = t.PW*1000.*sqrt(3.)/t.c[!p].U/3.;   //the max current of 10Kv coil
    t.c[!p].T  = 0.06;                  //The thickness of the insulation of first coil
    t.c[!p].C  = 0.0282;                //Conductivity of aluminum OM*m/mm**2
    t.c[!p].D  = 2.6989;                //Density aluminum kg/gm**3
    t.c[!p].s     = 50;                  //From
    t.c[!p].sp[0] = 200;                 //To
    t.c[!p].sp[1] = 1;                   //Spep
    t.c[!p].N     = 30;                  //From
    t.c[!p].Np[0] = 100;                 //To
    t.c[!p].Np[1] = 1;                   //Spep

    //t.NN = t.c[0].U/t.c[1].U;        // Transformation coefficient

    //double N[2], H[2], s[2][2];

    //N[0] = 10; N[1] = 100;
    //H[0] = 200; H[1] = 300;
    //s[0][0] = 0.1; s[0][1] = 6;
    //s[0][0] = 50; s[0][1] = 200;
    //s[1][0] = 50; s[1][1] = 200;
    //s[1][0] = 0.1; s[1][1] = 6;

    if(trans_3_phase(&t, 1000, p)){
        printf("No any result!!!\n");
        return;
    }

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

