#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libit/types.h"
#include "./trans.h"

#define PI 3.1415926535
#define MU 1.256637E-6
#define POWER 100


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

//Calculate square the top and bottom of the cylinder
//R1 - internal radius in mm
//R2 - External radius
inline double Sup(double R1, double R2)
{
    return 2*PI*(R2*R2 - R1*R1)/1000000.;
}

//Calculate square the sides of the cylinder
//h - vertical size in mm
//R2 - External radius
inline double Ssd(double h, double R2)
{
    return 2*PI*h*R2/1000000.;
}

//Calculate square of 3 phase magnetic core
//R1 - internal radius in mm
//R2 - External radius
//W  - The distance between the coils
inline double Sou(double R1, double R2, double W)
{
    return 4.*((PI+2)*PI*R1*R1*5/12. + (R2-R1)*(PI+2)*R1 + (PI+2)*R1*W/2.)/1000000.;
}

//Calculate square of 1 phase magnetic core
//h - vertical size in mm
//R1 - internal radius in mm
//R2 - External radius
inline double Sou1(double h, double R1, double R2)
{
    return (8*(PI+2)*PI*R1*R1/3. + 4*(R2-R1)*(PI+2)*R1 + 2*h*(PI+2)*R1)/1000000.;
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
//W  - The distance between the coils
inline double Vou(double R1, double R2, double W)
{
    return 4.*(PI*R1*R1*PI*R1*5./(4*6.) + (R2-R1)*PI*R1*R1/2. + PI*R1*R1*W/4.)/1000.;
}

//Calculate vоlume outside
//R1 - internal radius in mm
//R2 - External radius
inline double Vou1(double R1, double R2)
{
    return (8.*PI*R1*R1*PI*R1*2./(4*3.) + 4*(R2-R1)*PI*R1*R1/2.)/1000.;
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

int trans_optim(TRANS *t, double PT, int p, int ph)
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

    //printf("i = %d N0 = %f N1 = %f H = %f Hp1 = %f Hp2 = %f  Np1 = %f  Np2 = %f \n", i, N[0], N[1], t->H, t->Hp[0], t->Hp[1], t->c[i].Np[0], t->c[i].Np[1]);
    for(N[i] = t->c[i].N; N[i] <= t->c[i].Np[0]; N[i] += t->c[i].Np[1]){
        for(H = t->H; H <= t->Hp[0]; H += t->Hp[1]){
            for(s = t->c[0].s; s <= t->c[0].sp[0]; s += t->c[0].sp[1]){
                for(s1 = t->c[1].s ; s1 <= t->c[1].sp[0]; s1 += t->c[1].sp[1]){
                    //printf("N = %f H = %f s = %f s1 = %f \n", N[i], H, s, s1);
                    //printf("U0 = %f U1 = %f i = %d N = %f H = %f s1 = %f s2 = %f\n", t->c[0].U, t->c[1].U, i, t->c[i].N, t->H, t->c[0].s, t->c[1].s);
                    //s1 = s;
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


                    if(ph == 3) t->m.S = Sou(t->m.R, t->c[1].R, t->W);
                    else        t->m.S = Sou1(H, t->m.R, t->c[1].R);

                    t->i[0].S = Sup(t->m.R, t->i[0].R);
                    t->c[0].S = Sup(t->i[0].R, t->c[0].R);
                    t->i[1].S = Sup(t->c[0].R, t->i[1].R);
                    t->c[1].S = Sup(t->i[1].R, t->c[1].R) + Ssd(H, t->c[1].R);

                    t->S = 3*(t->m.S + t->i[0].S + t->c[0].S + t->i[1].S + t->c[1].S);
                    t->T = PT/(10*t->S);

                    if(ph == 3) t->m.V = Vol(H, 0, t->m.R) + Vou(t->m.R, t->c[1].R, t->W);
                    else        t->m.V = Vol(H, 0, t->m.R)*2. + Vou1(t->m.R, t->c[1].R);
                    t->i[0].V = Vol(H, t->m.R, t->i[0].R);
                    t->c[0].V = Vol(H, t->i[0].R, t->c[0].R);
                    t->i[1].V = Vol(H, t->c[0].R, t->i[1].R);
                    t->c[1].V = Vol(H, t->i[1].R, t->c[1].R);

                    t->m.M    = t->m.V*t->m.D/1000;        //Mass steel
                    t->i[0].M = t->i[0].V*t->i[0].D/1000.; //Mass of insulation
                    t->c[0].M = t->c[0].V*t->c[0].D/1000.; //Mass coil
                    t->i[1].M = t->i[1].V*t->i[1].D/1000.; //Mass of insulation
                    t->c[1].M = t->c[1].V*t->c[1].D/1000.; //Mass 1coil

                    if(ph == 3) t->M = 3.*(t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);
                    else        t->M =    (t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);

                    t->m.P = 3.*t->m.M*t->m.Lc;                         //Loss in steel
                    t->c[0].P = Pw(H, s , N[0], t->i[0].R, t->c[0].C, t->c[0].I, t->c[0].T);
                    t->c[1].P = Pw(H, s1, N[1], t->i[1].R, t->c[1].C, t->c[1].I, t->c[1].T);

                    if(ph == 3) t->P = t->m.P + 3.*(t->c[0].P + t->c[1].P);
                    else        t->P =    (t->m.P + t->c[0].P + t->c[1].P);

                    if((t->P <= PT) && (t->m.P <= PT*0.3)) {
                    //if(t->P <= PT) {
                         if(t->M < Mmin) {
                         //if(t->m.P < Pmin){
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

        t->c[0].Sw = 4.*t->c[0].L*sqrt(t->c[0].s)/1000.;
        t->c[1].Sw = 4.*t->c[1].L*sqrt(t->c[1].s)/1000.;
        //t->I = Iid1(t->H, t->i[p].R, t->c[p].R, t->m.B, t->m.Mu, t->c[p].N);
        t->I = Iid(t->H, t->m.R, t->c[1].R, t->c[p].U, t->m.Mu, t->c[p].N, 50);
        t->c[0].Rz = t->c[0].L*t->c[0].C/t->c[0].s;
        t->c[1].Rz = t->c[1].L*t->c[1].C/t->c[1].s;

        return 0;
    } else return 1;

}

void trans_calc(TRANS *t, int p, int ph)
{
    struct TRANS tm;
    int i;

    if(t->c[0].U > t->c[1].U)   { t->NN = t->c[0].U/t->c[1].U; i = 1; }
    else                        { t->NN = t->c[1].U/t->c[0].U; i = 0; }

    //t->c[p].N = Nc(t->m.B, t->c[p].U, t->m.R*sqrt(t->m.Sfc), 50);
    //if(p == i) t->c[!i].N = t->c[i].N*t->NN;
    //else       t->c[i].N = t->c[!i].N/t->NN;

    t->i[0].R = t->m.R + t->i[0].T;
    t->c[0].R = t->i[0].R + Rr(t->H, t->c[0].s, t->c[0].N);
    t->i[1].R = t->c[0].R + t->i[1].T;
    t->c[1].R = t->i[1].R + Rr(t->H, t->c[1].s, t->c[1].N);

    t->c[0].L = Ln1(t->c[0].Nh, t->c[0].Nr, t->c[0].h, t->i[0].R);
    t->c[1].L = Ln1(t->c[1].Nh, t->c[1].Nr, t->c[1].h, t->i[1].R);

    if(ph == 3) t->m.V = Vol(t->H, 0, t->m.R) + Vou(t->m.R, t->c[1].R, t->W);
    else        t->m.V = Vol(t->H, 0, t->m.R)*2. + Vou1(t->m.R, t->c[1].R);
    t->i[0].V = Vol(t->H, t->m.R, t->i[0].R);
    t->c[0].V = Vol(t->H, t->i[0].R, t->c[0].R);
    t->i[1].V = Vol(t->H, t->c[0].R, t->i[1].R);
    t->c[1].V = Vol(t->H, t->i[1].R, t->c[1].R);

    t->m.M    = t->m.V*t->m.D/1000;        //Mass steel
    t->i[0].M = t->i[0].V*t->i[0].D/1000.; //Mass of insulation
    t->c[0].M = t->c[0].V*t->c[0].D/1000.; //Mass coil
    t->i[1].M = t->i[1].V*t->i[1].D/1000.; //Mass of insulation
    t->c[1].M = t->c[1].V*t->c[1].D/1000.; //Mass 1coil

    if(ph == 3) t->M = 3.*(t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);
    else        t->M =    (t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);

    t->m.P = 3*t->m.M*t->m.Lc;                         //Loss in steel
    t->c[0].P = Pw(t->H, t->c[0].s, t->c[0].N, t->i[0].R, t->c[0].C, t->c[0].I, t->c[0].T);
    t->c[1].P = Pw(t->H, t->c[1].s, t->c[1].N, t->i[1].R, t->c[1].C, t->c[1].I, t->c[1].T);

    if(ph == 3) t->P = t->m.P + 3*(t->c[0].P + t->c[1].P);
    else        t->P =    (t->m.P + t->c[0].P + t->c[1].P);


    t->c[0].Sw = 2.*t->c[0].L*(t->c[0].h + t->c[0].w)/1000.;
    t->c[1].Sw = 2.*t->c[1].L*(t->c[1].h + t->c[1].w)/1000.;
    t->I = Iid1(t->H, t->i[p].R, t->c[p].R, t->m.B, t->m.Mu, t->c[p].N);
    t->c[0].Rz = t->c[0].L*t->c[0].C/t->c[0].s;
    t->c[1].Rz = t->c[1].L*t->c[1].C/t->c[1].s;

}

void real_coil(TRANS *t, int p)
{
    t->c[0].Nr = ceil((t->c[0].R - t->i[0].R)/sqrt(t->c[0].s));
    t->c[0].Nh = round(t->c[0].N/t->c[0].Nr);
    if(t->c[0].Nr*t->c[0].Nh != t->c[0].N) {
        t->c[0].N = t->c[0].Nr*t->c[0].Nh;
        t->c[0].s = (t->c[0].R - t->i[0].R)*t->H/t->c[0].N;
    }
    t->c[0].h = (t->c[0].R - t->i[0].R)/t->c[0].Nr;
    t->c[0].w = t->H/t->c[0].Nh;
    //t->c[0].L = Ln1(t->c[0].Nh, t->c[0].Nr, t->c[0].h, t->i[0].R);
    //t->c[0].Rz = t->c[0].L*t->c[0].C/((t->c[0].h - 2.*t->c[0].T)*(t->c[0].w - 2.*t->c[0].T));
    //t->c[0].P = t->c[0].Rz*t->c[0].I*t->c[0].I;
    //printf("Nr = %f Nh = %f N = %f h = %f w = %f L = %f Rz = %f P = %f\n",
    //       t->c[0].Nr, t->c[0].Nh, t->c[0].N, t->c[0].h, t->c[0].w, t->c[0].L, t->c[0].Rz, t->c[0].P);

    t->c[1].Nr = ceil((t->c[1].R - t->i[1].R)/sqrt(t->c[1].s));
    t->c[1].Nh = round(t->c[1].N/t->c[1].Nr);
    if(t->c[1].Nr*t->c[1].Nh != t->c[1].N) {
        t->c[1].N = t->c[1].Nr*t->c[1].Nh;
        t->c[1].s = (t->c[1].R - t->i[1].R)*t->H/t->c[1].N;
    }
    t->c[1].h = (t->c[1].R - t->i[1].R)/t->c[1].Nr;
    t->c[1].w = t->H/t->c[1].Nh;
    //t->c[1].L = Ln1(t->c[1].Nh, t->c[1].Nr, t->c[1].h, t->i[1].R);
    //t->c[1].Rz = t->c[1].L*t->c[1].C/((t->c[1].h - 2.*t->c[1].T)*(t->c[1].w - 2.*t->c[1].T));
    //t->c[1].P = t->c[1].Rz*t->c[1].I*t->c[1].I;
    //printf("Nr = %f Nh = %f N = %f h = %f w = %f L = %f Rz = %f P = %f\n",
    //       t->c[1].Nr, t->c[1].Nh, t->c[1].N, t->c[1].h, t->c[1].w, t->c[1].L, t->c[1].Rz, t->c[1].P);

    t->M = 3.*(t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);
    t->P = 3.*(t->m.P + t->c[0].P + t->c[1].P);

    if(t->c[0].U > t->c[1].U)   t->NN = t->c[0].N/t->c[1].N;
    else                        t->NN = t->c[1].N/t->c[0].N;

    if(t->c[p].U < t->c[!p].U)  t->c[!p].U = t->NN*t->c[p].U;
    else                        t->c[!p].U = t->c[p].U/t->NN;

}

void print_result(TRANS *t, int p)
{
    printf("Входные параметры\n");
    printf("Схема подключения треугольник - звезда\n");
    printf("Мощность %f кВ\n", t->PW);
    printf("Максимальное поле в магнитном сердечнике %f Тл\n", t->m.B);
    printf("Потери в магнитном сердечнике %f Вт/кг\n", t->m.Lc);
    printf("Коэффициент заполнения магнитопровода %f %%\n",t->m.Sfc*100.);
    printf("Коэффициент трансформации обмоток %f \n", t->NN);

    printf("\nВыходные параметры\n");
    printf("Суммарные потери %f Вт\n", t->P);
    printf("Потери холостого хода %f Вт\n", t->m.P);
    printf("Потери короткого замыкания %f Вт\n", t->P - t->m.P);
    printf("Масса %f кг\n", t->M);
    printf("Площадь поверхности %f м**2\n", t->S);
    printf("Темпрература над окружением %f K\n", t->T);
    printf("Ток холостого хода %f А\n", t->I);
    printf("Ток холостого хода %f %%\n", t->I*100./t->c[p].I);

    printf("\nМагнитный сердечник\n");
    printf("Масса  %f кг\n", t->m.M);
    printf("Потери %f Вт\n", t->m.P);
    printf("Радиус %f мм\n", t->m.R);

    printf("\nПервая катушка\n");
    printf("Напряжение %f В\n", t->c[0].U);
    printf("Максимальный ток %f А\n", t->c[0].I);
    printf("Толщина покрытия  %f мм\n",t->c[0].T);
    printf("Толщина изоляции между магнитопроводом и обмоткой %f мм\n",t->i[0].T);
    printf("Количество витков %f\n", t->c[0].N);
    printf("Масса  %f кг\n", t->c[0].M);
    printf("Потери %f Вт\n", t->c[0].P);
    printf("Сопротивление %f Ом\n", t->c[0].Rz);
    printf("Длина %f м\n", t->c[0].L);
    printf("Сечение %f мм**2\n", t->c[0].s);
    printf("Площадь поверхности провода %f м**2\n", t->c[0].Sw);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", t->i[0].R, t->c[0].R, t->H);
    printf("Число витков по высоте %f по радиусу %f \n", t->c[0].Nh, t->c[0].Nr);
    printf("Высота сечения провода %f ширина %f \n", t->c[0].h, t->c[0].w);

    printf("\nВторая катушка\n");
    printf("Напряжение %f В\n", t->c[1].U);
    printf("Максимальный ток %f А\n", t->c[1].I);
    printf("Толщина покрытия  %f мм\n",t->c[1].T);
    printf("Толщина изоляции между обмотками                  %f мм\n",t->i[1].T);
    printf("Количество витков %f\n", t->c[1].N);
    printf("Масса  %f кг\n", t->c[1].M);
    printf("Потери %f Вт\n", t->c[1].P);
    printf("Сопротивление %f Ом\n", t->c[1].Rz);
    printf("Длина %f м\n", t->c[1].L);
    printf("Сечение %f мм**2\n", t->c[1].s);
    printf("Площадь поверхности провода %f м**2\n", t->c[1].Sw);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", t->i[1].R, t->c[1].R, t->H);
    printf("Число витков по высоте %f по радиусу %f \n", t->c[1].Nh, t->c[1].Nr);
    printf("Высота сечения провода %f ширина %f \n", t->c[1].h, t->c[1].w);
}

void trans(void)
{
    int i, n = 7;
    TRANS tr[n], t;
    COIL coil[9];
    int p = 1;                              //The primary coil

    for(i=0; i < n; i++) {
        tr[i] = (TRANS) {.PW = POWER, .W = 2, .H = 200, .Hp[0] = 400, .Hp[1] = 1 };
        tr[i].m = (MCORE) {.D = 7.65, .Mu = 10000., .B = 1.85, .Lc = 1., .Sfc = 0.955, .R = 0 };
        tr[i].i[0] = (INS) {.D = 1.5, .T = 1.};
        tr[i].i[1] = (INS) {.D = 1.5, .T = 2.};
    }

    //Tree phase
    coil[0] = (COIL) {.U = 400,   .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 50,  .sp[0] = 200, .sp[1] = 1,  .N = 30, .Np[0] = 100, .Np[1] = 1};
    coil[0].I = POWER*1000*sqrt(3.)/coil[0].U/3.;
    coil[1] = (COIL) {.U = 315,   .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 50,  .sp[0] = 200, .sp[1] = 1,  .N = 30, .Np[0] = 100, .Np[1] = 1};
    coil[1].I = POWER*1000/coil[1].U/3.;
    coil[2] = (COIL) {.U = 10000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 0.1, .sp[0] = 6,  .sp[1] = 0.1, .N = 30, .Np[0] = 100, .Np[1] = 1};
    coil[2].I = POWER*1000/coil[2].U/3.;
    coil[3] = (COIL) {.U = 12000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 0.1, .sp[0] = 6,  .sp[1] = 0.1, .N = 30, .Np[0] = 100, .Np[1] = 1};
    coil[3].I = POWER*1000/coil[3].U/3.;
    coil[4] = (COIL) {.U = 24000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 0.1, .sp[0] = 6,  .sp[1] = 0.1, .N = 30, .Np[0] = 100, .Np[1] = 1};
    coil[4].I = POWER*1000/coil[4].U/3.;
    coil[5] = (COIL) {.U = 35000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 0.1, .sp[0] = 6,  .sp[1] = 0.1, .N = 30, .Np[0] = 100, .Np[1] = 1};
    coil[5].I = POWER*1000/coil[5].U/3.;
    coil[6] = (COIL) {.U = 36000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 0.1, .sp[0] = 6,  .sp[1] = 0.1, .N = 30, .Np[0] = 100, .Np[1] = 1};
    coil[6].I = POWER*1000/coil[6].U/3.;

    for(i=0; i < n; i++) {
        memcpy (&tr[i].c[0], &coil[0],   sizeof(coil[0]));
        memcpy (&tr[i].c[1], &coil[i+1], sizeof(coil[1]));
    }
    memcpy (&t, &tr[2], sizeof(t));


    if(trans_optim(&t, 1000, p, 3)){
        printf("No any result!!!\n");
        return;
    }
    /*
    //One Phase
    coil[7] = (COIL) {.U = 230, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 0.1, .sp[0] = 10,  .sp[1] = 0.1};
    coil[7].I = 1.5*1000/coil[7].U;
    coil[8] = (COIL) {.U = 40, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 1, .sp[0] = 100,  .sp[1] = 1, .N = 20, .Np[0] = 100, .Np[1] = 1};
    coil[8].I = 1.5*1000/coil[8].U;

    t = (TRANS) {.PW = 1500, .H = 20, .Hp[0] = 200, .Hp[1] = 1 };
    t.m = (MCORE) {.D = 7.65, .Mu = 20000., .B = 1.9, .Lc = 1., .Sfc = 0.955, .R = 24 };
    t.i[0] = (INS) {.D = 1.5, .T = 1.};
    t.i[1] = (INS) {.D = 1.5, .T = 1.};
    memcpy (&t.c[0], &coil[7], sizeof(coil[0]));
    memcpy (&t.c[1], &coil[8], sizeof(coil[1]));
    p = 0;

    if(trans_optim(&t, 1500*0.04, p, 1)){
        printf("No any result!!!\n");
        return;
    }
    */
    print_result(&t, p);

    real_coil(&t, p);
    trans_calc(&t, p, 3);

    print_result(&t, p);
    printf("H = %f R1 = %f R2 = %f R3 = %f R4 = %f R5 = %f \n", t.H, t.m.R, t.i[0].R, t.c[0].R, t.i[1].R, t.c[1].R);

}

