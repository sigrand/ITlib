#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libit/types.h"
#include "./trans.h"

#define PI 3.1415926535
//#define E 2.71828182
#define MU 1.256637E-6

//#define POWER 100


//Calculate square of the coil,
//s - square of one wire in mm**2
//n - number of coils
inline double sq(double s, double n)
{
    return n*s;
}

//Calculate square of the magnetic core,
//R - Core radius
//k - coefficient
inline double SQ(double R2, double R1)
{
    return PI*R2*R2;
}

//Calculate square of the magnetic core,
//R - Core radius
//k - coefficient
inline double SQ1(double R2, double R1)
{
    return 2*R2*R2*PI/3 + 1.73*R1*R1/2 + 2*(R2-R1)*R1 + PI*(R2-R1)*(R2-R1)/3;
}

//Calculate square elips

inline double SQL(double a, double b)
{
    return PI*a*b;
}

//Calculate  length of the magnetic core,
//R - Core radius
//k - coefficient
//h - the length from magnetic core
inline double LN(double R2, double R1)
{
    return 2*PI*R2;
}

//Calculate  length of the magnetic core,
//R - Core radius
//k - coefficient
//h - the length from magnetic core
inline double LN1(double R2, double R1)
{
    return 4*PI*R2/3. + 2*R1 + 2*PI*(R2-R1)/3;
}

//Calculate  length elips

inline double LNL(double a, double b)
{
    return PI*(3*(a + b) - sqrt((3*a + b)*(a + 3*b)));
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

//Calculate spiral h correction
//h - vertical one wire size
//Nh - number of coils
//R - radius of spiral
inline double Hspi(double h, double Nh, double R)
{
    return (Nh+1)*h/sqrt(1 - h*h/(4.*PI*PI*R*R));
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
inline double Ln2(double h, double s, double n, double r)
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

//Calculate lenght of spiral
//st - step
//fi - angel
inline double LnSp(double st, double fi)
{
    return st*(fi*sqrt(1 + fi*fi) + ln(fi + sqrt(1 + fi*fi)))/2.;
}

//Calculate lenght of coil in m
//Nh - number of turns in vertical
//Nr - number of turns in radial
//H - the height of the wire
//R1 - inside radius of coil
//w  - spiral width
inline double LnS(double Nh, double Nr, double h, double R1, double w, fp ln)
{
    int i;
    double L = 0, L1, L2;
    for(i = 0; i < Nr; i++){
        //L1 =  2.*PI*(R1 + h/2. + h*i);
        L1 = (ln)(R1 + h/2. + h*i, 0);
        L2 = L1*sqrt(1. + 1./(L1*L1/(w*w) - 1));
        //printf("%d   L1 = %f  L2 = %f\n", i, L1, L2);
        L+= L2*Nh/1000.;
    }
    return L;
}

//Calculate radial size of the coil
//h - vertical size in mm
//s - square of one wire in mm**2
//n - number of coils
inline double Rr(double h, double s, double n)
{
    return sq(s,n)/h;
}

//Calculate radial size of the coil with spiral wire
//h - vertical size in mm
//s - square of one wire in mm**2
//n - number of coils
inline double Rrsp(double h, double s, double n)
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
//Ln - lenght ow wire
//h - the vertical size of wire
//w - the horisontal size of wire
//C - Conductivity
//I - Current
//Ti - The thickness of the insulation
//Ks - Skin effect factor
inline double Pw2(double Ln, double h, double w, double C, double I, double Ti, double Ks)
{
    return C*Ln*I*I*Ks/((h - 2.*Ti)*(w - 2.*Ti));
}

inline double Pw3(double Ln, double h, double w, double C, double I, double Ti, double Ks)
{
    return C*Ln*I*I*Ks/((h)*(w));
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
//inline double Vol(double h, double R1, double R2)
//{
//    return PI*h*(R2*R2 - R1*R1)/1000.;
//}

//Calculate volume of tor
//R - Tor radius
//S - the square
inline double Vol_tor(double R, double S)
{
    return 2*S*PI*R/1000.;
}

//Calculate volume
//h - vertical size in mm
//R1 - internal radius in mm
//R2 - External radius
inline double Vol(double h, double S1, double S2)
{
    return (S2 - S1)*h/1000.;
}

//Calculate vоlume outside
//R1 - internal radius in mm
//R2 - External radius
//W  - The distance between the coils
//inline double Vou(double R1, double R2, double W)
//{
//    return 4.*(PI*R1*R1*PI*R1*5./(4*6.) + (R2-R1)*PI*R1*R1/2. + PI*R1*R1*W/4.)/1000.;
//}

//Calculate vоlume outside
//R1 - internal radius in mm
//R2 - External radius
//W  - The distance between the coils
inline double Vou(double S, double R1,  double R2, double W)
{
    return 4.*(S*PI*R1*5./(4*6.) + (R2-R1)*S/2. + S*W/4.)/1000.;
}

//Calculate vоlume outside
//R1 - internal radius in mm
//R2 - External radius
//W  - The distance between the coils
inline double Vou2(double S, double R1,  double R2, double W)
{
    return 4.*(S*PI*R1*5./(4*6.) + (R2- 2*R1 + 0.866*R1)*S/2. + S*W/4.)/1000.;
}

//Calculate vоlume outside
//R1 - internal radius in mm
//R2 - External radius
inline double Vou1(double R1, double R2)
{
    return (8.*PI*R1*R1*PI*R1*2./(4*3.) + 4*(R2-R1)*PI*R1*R1/2.)/1000.;
}

//The cross section of tor wire
//R - Tor radius
//S1 - Square of cross section core
//S2 - Square of cross section core + first coil
//N - Number of windings
//Ln - length of wire
inline double CS_tor(double R, double S1, double S2, double N, double Ln)
{
    return 2*(S2-S1)*PI*R/((N-1)*Ln);
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
    return (B/1.414)*Lm(h, R1, R2)/(mu*MU*N*1000.);
}

//The magnetic core radius
//B - Max Magnetic field
//U  - Input Voltage
//N - number of coils
//F - frequency
inline double Rc(double B, double U, double N, double F)
{
    return sqrt(U/(1.414*PI*PI*F*B*N));
}

//The first coil turns
//B - Max Magnetic field
//U  - Input Voltage
//S - The square of magnetic core
//F - frequency
inline double Nc(double B, double U, double Sq, double F)
{
    return U*1000000./(1.414*PI*F*B*Sq);
}

inline double Bc(double N, double U, double Sq, double F)
{
    return U*1000000./(1.414*PI*F*N*Sq);
}

int trans_optim_tor(TRANS *t, double PW, double M, double FR)
{
    struct TRANS tm;
    double Mmin=10000000, Pmin = 10000000, Imax = 0, s, s1, N, Nr[2], l, SP; // Ln[2];
    double tmp, Ks, R1, R, r, sq, a, b, a1, b1, mR, xl, N1, L;
    int j=0;
    double U, UP = 1.8;
    double L0, L1, L2, fi, V, Ms, C, h, E, cosf, st = 0.35/(2*PI);
    //fp sq, ln;
    //fp1 vou;

    //sq = SQ; ln = LN;
    //R1 = t->m.R;

    t->NN = t->c[1].U/t->c[0].U;

    for(mR = t->m.R; mR <= t->m.Rp[0]; mR += t->m.Rp[1]) {
        for(R1 = t->m.R1; R1 <= t->m.R1p[0]; R1 += t->m.R1p[1]){
            //for(N1 = t->c[0].N; N1 <= t->c[0].Np[0]; N1 += t->c[0].Np[1]){
                for(N = t->c[1].N; N  <= t->c[1].Np[0]; N += t->c[1].Np[1]){
                    for(a = R1; a <=  R1; a += 0.2) {
                        b = R1;
                        t->c[0].Rl = sqrt(t->c[0].U*t->c[0].U/(t->c[0].I*t->c[0].I) - t->c[0].Rz*t->c[0].Rz) - 0.3;
                        L = t->c[0].Rl/(2*PI*FR);
                        N1 = sqrt(2*PI*L*mR/(1.26e-6*7.4*SQL(a,b)));
                        sq = t->c[0].s*(N1+1)/0.9;
                        if(PI*(mR - b)*(mR - b) <= sq) continue;

                        //R = t->m.R - sqrt((PI*(t->m.R - R1)*(t->m.R - R1) - sq)/PI);
                        //t->c[0].L = LN((R1 + R)/2., 0)*(N1+1)/1000.;

                        R = mR - sqrt((PI*(mR - b)*(mR - b) - sq)/PI);
                        b1 = (R + b)/2; a1 = a + (R - b)/2;
                        t->c[0].L = LNL(a1, b1)*(N1+1)/1000.;

                        //printf(" L = %f L = %f R = %f a1 = %f b1 = %f a = %f b = %f R1 = %f R = %f m.R = %f\n", tmp, t->c[0].L, R, a1, b1, a , b, R1, R, t->m.R);
                        t->c[0].V = t->c[0].s*t->c[0].L/1000.;

                        //printf("N = %f s = %f Sq = %f R = %f L = %f V = %f\n", N1, t->c[0].s, sq, R, t->c[0].L, t->c[0].V);
                        if (R + b >= mR) continue;

                        t->c[0].h =  R - b;
                        t->c[0].w =  t->c[0].s/t->c[0].h;

                        //Remove thickness of the wire insulation
                        t->c[1].s = SQ(mR-R,0)*0.9/N;
                        r = sqrt(t->c[1].s/PI);
                        //t->c[1].L = LN(R + r/2., 0)*N/1000. + 0.3;
                        t->c[1].L = LNL(R + r/2., a + R - b + r/2.)*N/1000. + 0.3;
                        

                        //t->m.V = Vol_tor(t->m.R, SQ(R1,0))/1000;
                        t->m.V = Vol_tor(mR, SQL(a,b))/1000;
                        //t->c[0].V = (Vol_tor(t->m.R, SQ(R,0)) - Vol_tor(t->m.R, SQ(R1,0)))/1000;
                        
                        //The volume of spiral
                        E = sqrt(1-(b*b)/(a*a));
                        l = 0;
                        cosf =  l/sqrt(4*b*b - E*E*l*l);
                        h = b - b*sqrt(1 - cosf*cosf)/sqrt(1 - E*E*cosf*cosf);
                        h = 0;

                        fi = (mR - b + h)/st;
                        L0 = LnSp(st, fi);
                        fi = (mR)/st;
                        L1 = LnSp(st, fi) - L0;
                        fi = (mR + b - h)/st;
                        L2 = LnSp(st, fi) - L0 - L1;
                        V = ((a+l/2)*(L1 + L2))*0.35/1000000;
                        //V = (L1*a + L2*a)*0.35/1000000;
                        Ms = V*t->m.D*(t->m.Sfc + 1)/(2);

                        //b/sqrt(1 - E*E*cos(f)*cos(f))

                        t->c[1].V = t->c[1].s*t->c[1].L/1000.;
                        
                        t->m.M    = t->m.V*t->m.D*(t->m.Sfc + 1)/(2);       //Mass steel
                        t->c[0].M = t->c[0].V*t->c[0].D;                    //Mass coil
                        t->c[1].M = t->c[1].V*t->c[1].D;                    //Mass coil

                        t->M =   t->m.M + t->c[0].M + t->c[1].M;

                        t->m.P = t->m.M*t->m.Lc*t->m.Sfc;                         //Loss in steel
                        Ks = t->Ks[0]*t->c[0].s + t->Ks[1];
                        t->c[0].Rz = t->c[0].L*t->c[0].C*Ks/t->c[0].s;

                        //L = 7.4*1.26e-6*N1*N1*SQL(a,b)*0.001/(2*PI*mR);
                        //t->c[0].Rl = L*2*PI*FR;
                        //t->c[0].I = t->c[0].U/sqrt(t->c[0].Rz*t->c[0].Rz + (t->c[0].Rl + 0.3)*(t->c[0].Rl + 0.3));
                        U = t->c[0].Rl*t->c[0].I;
                        //if(t->c[0].I >100) continue;

                        t->c[0].P = t->c[0].Rz*t->c[0].I*t->c[0].I;
                        //t->c[0].P = Pw3(t->c[0].L, t->c[0].h, t->c[0].w, t->c[0].C, t->c[0].I, t->c[0].T, Ks);
                        Ks = t->Ks[0]*t->c[1].s + t->Ks[1];
                        t->c[1].Rz = t->c[1].L*t->c[1].C*Ks/t->c[1].s;
                        //t->c[1].U = (t->c[0].U - t->c[0].I*t->c[0].Rz)*N/t->c[0].N;
                        //t->c[1].U = (t->c[0].U - t->c[0].I*sqrt(t->c[0].Rz*t->c[0].Rz + 0.09))*N/t->c[0].N;
                        //xl = 1.414*PI*FR*SQL(a, b)*B*t->c[0].N/(t->c[0].I*1000000);
                        //U = (sqrt((t->c[0].U*t->c[0].U)/(t->c[0].I*t->c[0].I) - t->c[0].Rz*t->c[0].Rz) - 0.03)*t->c[0].I;
                        //U = t->c[0].I*t->c[0].Rl;
                        //B = 1.414*L*t->c[0].I*1000000/(N1*SQL(a, b));
                        t->c[1].U = U*N/N1;
                        t->m.B = Bc(N, t->c[1].U, SQL(a,b), FR);
                        t->m.H = t->c[0].I*t->c[0].N*1000/(1.4*PI*mR);
                        t->m.B1 = 1.26e-6*7.4*t->m.H;
                        t->m.Mu = t->m.B/(1.26e-6*t->m.H);
                        t->c[0].N1 = Nc(t->m.B, t->c[1].U, SQL(a,b), FR);
                        t->c[0].N2 = Nc(t->m.B1, t->c[1].U, SQL(a,b), FR);

                        //printf("xl = %f R = %f U = %f U1 = %f\n", xl, t->c[0].Rz,t->c[1].U, tmp);

                        if(t->c[1].U <= 2) continue;

                        //t->c[1].U = (t->c[0].U - t->c[0].I*t->c[0].Rz)*t->c[1].N/t->c[0].N;

                        t->c[1].I = (t->c[1].U - UP)/t->c[1].Rz;
                        t->c[1].P = t->c[1].Rz*t->c[1].I*t->c[1].I;
                        //if(t->c[1].P > PW*1000 - (t->m.P + t->c[0].P)){
                        //    t->c[1].P = PW*1000 - (t->m.P + t->c[0].P);
                        //    t->c[1].I = t->c[1].P/t->c[1].U;
                        //}
                        //t->c[1].P = t->c[1].C*t->c[1].L*t->c[1].I*t->c[1].I*Ks/t->c[1].s;

                        t->P =    (t->m.P + t->c[0].P + t->c[1].P + UP*t->c[1].I);
                        //printf("P = %f I = %f M = %f\n", t->P, t->c[1].I, t->M );
                       //printf("m.R = %f m.R1 = %f c.R = %f N = %f s = %f L = %f  m.V = %f  c[0].V = %f m.M = %f c.M = %f P = %f I = %f I = %f U = %f L = %f RL = %f\n",
                       //        t->m.R, R1, R, N1, t->c[0].s, t->c[0].L, t->m.V, t->c[0].V, t->m.M, t->c[0].M, t->P, t->c[0].I,t->c[1].I, t->c[1].U, L, RL);
                        if(t->P > t->c[0].I*t->c[0].U) continue;

                        if(t->M <= M ) {
                            if(t->c[1].I > Imax) {
                            //if(t->M < Mmin){
                                Mmin = t->M;
                                Pmin = t->P;
                                Imax = t->c[1].I;
                                memcpy (&tm, t, sizeof(tm));
                                tm.m.R1 = R1;
                                tm.c[0].R = R;
                                tm.c[1].N = N;
                                tm.c[0].N = N1;
                                tm.m.R = mR;
                                tm.m.a = a; tm.m.b = b;
                                j++;
                            }
                        }
                    }
                //}
                //}
            }
        }
    }

    if(j) {
        memcpy (t, &tm, sizeof(tm));

        t->c[0].Sw = 2.*t->c[0].L*(t->c[0].h + t->c[0].w)/1000.;
        //t->c[1].Sw = 2.*t->c[1].L*(t->c[1].h + t->c[1].w)/1000.;
        //t->I = Iid1(t->H, t->i[p].R, t->c[p].R, t->m.B, t->m.Mu, t->c[p].N);
        t->I = Iid(t->H, t->m.R, t->c[1].R, t->c[0].U, t->m.Mu, t->c[0].N, FR);
        //Ks = t->Ks[0]*t->c[0].s + t->Ks[1];
        //t->c[0].Rz = t->c[0].L*t->c[0].C*Ks/t->c[0].s;
        //Ks = t->Ks[0]*t->c[1].s + t->Ks[1];
        //t->c[1].Rz = t->c[1].L*t->c[1].C*Ks/t->c[1].s;

        printf("m.R = %f m.R1 = %f c.R = %f m.V = %f m.M = %f m.P = %f P = %f a = %f b = %f\n",
               t->m.R, t->m.R1, t->c[0].R, t->m.V,  t->m.M,  t->m.P, t->P, t->m.a, t->m.b);
        printf("N = %f s = %f L = %f  V = %f M = %f P = %f Rz = %f I = %f w = %f h = %f RL = %f\n",
               t->c[0].N, t->c[0].s, t->c[0].L, t->c[0].V, t->c[0].M, t->c[0].P, t->c[0].Rz, t->c[0].I, t->c[0].w,  t->c[0].h, t->c[0].Rl);
        printf("N = %f s = %f L = %f  V = %f M = %f P = %f Rz = %f I = %f U = %f PER = %f Pv = %f M = %f\n",
               t->c[1].N, t->c[1].s, t->c[1].L, t->c[1].V, t->c[1].M, t->c[1].P, t->c[1].Rz, t->c[1].I, t->c[1].U, U*t->c[1].I/t->P, UP*t->c[1].I, t->M);
        printf("B = %f B1 = %f H = %f mu = %f N1 = %f N2 = %f\n", tm.m.B, tm.m.B1, tm.m.H, tm.m.Mu, t->c[0].N1, t->c[0].N2);
        printf("h = %f L1 = %f L2 = %f V = %f M = %f\n", h, L1, L2, V, Ms);

        return 0;
    } else return 1;

}

int trans_optim(TRANS *t, double PT, int p, double FR, int type)
{
    struct TRANS tm;
    double Mmin=10000000, Pmin = 10000000, s, s1, H, N[2], Nr[2], l, SP; // Ln[2];
    double Hb, sb, s1b, Nb, Ks, R1;
    int i, j=0, k;
    fp sq, ln;
    fp1 vou;

    if(type == Round_3phase)
    {
        sq = SQ; ln = LN; vou = Vou;
    }
    else if(type == Round_120_3phase)
    {
        sq = SQ1; ln = LN1; vou = Vou2;
    }
    else if(type == Tore_1phase)
    {
        sq = SQ; ln = LN;
    }
    R1 = t->m.R;

    k = !t->m.R ? 1 : 0;
    if(t->c[0].U > t->c[1].U)   { t->NN = t->c[0].U/t->c[1].U; i = 1; }
    else                        { t->NN = t->c[1].U/t->c[0].U; i = 0; }

    Nb = t->c[i].N; Hb = t->H; sb = t->c[0].s; s1b = t->c[1].s, Nr[0] = t->c[0].Nr, Nr[1] = t->c[1].Nr;

    //printf("i = %d N0 = %f N1 = %f H = %f Hp1 = %f Hp2 = %f  Np1 = %f  Np2 = %f \n", i, N[0], N[1], t->H, t->Hp[0], t->Hp[1], t->c[i].Np[0], t->c[i].Np[1]);

    for( ; t->m.R <= t->m.Rp[0]; t->m.R += t->m.Rp[1]){
        for(H = Hb; H <= t->Hp[0]; H += t->Hp[1]){
            for(s = sb; s <= t->c[0].sp[0]; s += t->c[0].sp[1]){
                for(s1 = s1b; s1 <= t->c[1].sp[0]; s1 += t->c[1].sp[1]){
                    //printf("R = %f  H = %f s = %f s1 = %f\n", t->m.R, H, s, s1);


                    N[i] = Nc(t->m.B, t->c[i].U, (sq)(t->m.R*sqrt(t->m.Sfc), R1), FR);
                    //N[i] = Nc(t->m.B, t->c[i].U, t->m.R*sqrt(t->m.Sfc), FR);
                    N[!i] = N[i]*t->NN;

                    t->i[0].R = t->m.R + t->i[0].T;
                    t->c[0].R = t->i[0].R + Rr(H, s , N[0]);
                    t->i[1].R = t->c[0].R + t->i[1].T;
                    t->c[1].R = t->i[1].R + Rr(H, s1, N[1]);

                    t->c[0].Nr = 1;
                    //t->c[0].Nr = Nr[0] ? Nr[0] : ceil((t->c[0].R - t->i[0].R)/sqrt(s));
                    t->c[0].Nh = round(N[0]/t->c[0].Nr);
                    t->c[0].N = t->c[0].Nr*t->c[0].Nh;
                    //t->c[0].s = (t->c[0].R - t->i[0].R)*H/t->c[0].N;
                    t->c[0].h = (t->c[0].R - t->i[0].R)/t->c[0].Nr;
                    t->c[0].w1 = H/t->c[0].Nh;

                    //Spiral calc
                     SP = (ln)(t->i[0].R, R1)*(t->c[0].Nh+1)/H;
                    t->c[0].w = (ln)(t->i[0].R, R1)*sqrt(1/(1+SP*SP));
                    t->c[0].s = t->c[0].h*t->c[0].w;
                    //printf("w0 = %f w1 = %f\n", t->c[0].w1, t->c[0].w);

                    t->c[1].Nr = 1;
                    //t->c[1].Nr = Nr[1] ? Nr[1] : ceil((t->c[1].R - t->i[1].R)/sqrt(s1));
                    t->c[1].Nh = round(N[1]/t->c[1].Nr);
                    t->c[1].N = t->c[1].Nr*t->c[1].Nh;
                    //t->c[1].s = (t->c[1].R - t->i[1].R)*H/t->c[1].N;
                    t->c[1].h = (t->c[1].R - t->i[1].R)/t->c[1].Nr;
                    t->c[1].w1 = H/t->c[1].Nh;

                    //Spiral calc
                    SP = (ln)(t->i[1].R, 0)*(t->c[1].Nh+1)/H;
                    t->c[1].w = (ln)(t->i[1].R, 0)*sqrt(1/(1+SP*SP));
                    t->c[1].s = t->c[1].h*t->c[1].w;

                    t->c[0].L = LnS(t->c[0].Nh, t->c[0].Nr, t->c[0].h, t->i[0].R, t->c[0].w, ln);
                    t->c[1].L = LnS(t->c[1].Nh, t->c[1].Nr, t->c[1].h, t->i[1].R, t->c[1].w, ln);


                    if(type == Round_3phase || type == Round_120_3phase) t->m.S = Sou(t->m.R, t->c[1].R, t->W);
                    else        t->m.S = Sou1(H, t->m.R, t->c[1].R);

                    t->i[0].S = Sup(t->m.R, t->i[0].R);
                    t->c[0].S = Sup(t->i[0].R, t->c[0].R);
                    t->i[1].S = Sup(t->c[0].R, t->i[1].R);
                    t->c[1].S = Sup(t->i[1].R, t->c[1].R) + Ssd(H, t->c[1].R);

                    t->S = 3*(t->m.S + t->i[0].S + t->c[0].S + t->i[1].S + t->c[1].S);
                    t->T = PT/(10*t->S);

                    if(type == Round_3phase || type == Round_120_3phase)
                        t->m.V = Vol(H, 0, (sq)(t->m.R, R1)) + (vou)((sq)(t->m.R, R1), t->m.R, t->c[1].R, t->W);
                    else        t->m.V = Vol(H, 0, (sq)(t->m.R, R1))*2. + Vou1(t->m.R, t->c[1].R);

                    t->i[0].V = Vol(H, (sq)(t->m.R, R1)   , (sq)(t->i[0].R, R1));
                    t->c[0].V = Vol(H, (sq)(t->i[0].R, R1), (sq)(t->c[0].R, R1));
                    t->i[1].V = Vol(H, (sq)(t->c[0].R, R1), (sq)(t->i[1].R, R1));
                    t->c[1].V = Vol(H, (sq)(t->i[1].R, R1), (sq)(t->c[1].R, R1));

                    t->m.M    = t->m.V*t->m.D*(t->m.Sfc + 1)/(2*1000);        //Mass steel
                    t->i[0].M = t->i[0].V*t->i[0].D/1000.; //Mass of insulation
                    t->c[0].M = t->c[0].V*t->c[0].D/1000.; //Mass coil
                    t->i[1].M = t->i[1].V*t->i[1].D/1000.; //Mass of insulation
                    t->c[1].M = t->c[1].V*t->c[1].D/1000.; //Mass 1coil

                    if(type == Round_3phase || type == Round_120_3phase) t->M = 3.*(t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);
                    else        t->M =    (t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);


                    t->m.P = 3.*t->m.M*t->m.Lc*t->m.Sfc;                         //Loss in steel
                    Ks = t->Ks[0]*t->c[0].s + t->Ks[1];
                    t->c[0].P = Pw2(t->c[0].L, t->c[0].h, t->c[0].w, t->c[0].C, t->c[0].I, t->c[0].T, Ks);
                    Ks = t->Ks[0]*t->c[1].s + t->Ks[1];
                    t->c[1].P = Pw2(t->c[1].L, t->c[1].h, t->c[1].w, t->c[1].C, t->c[1].I, t->c[1].T, Ks);

                    if(type == Round_3phase || type == Round_120_3phase) t->P = t->m.P + 3.*(t->c[0].P + t->c[1].P);
                    else        t->P =    (t->m.P + t->c[0].P + t->c[1].P);

                    if((t->P <= PT) && (t->m.P <= PT*1)) {
                    //if(t->P <= PT) {
                         if(t->M < Mmin) {
                         //if(t->m.P < Pmin){
                            Mmin = t->M;
                            Pmin = t->P;
                            memcpy (&tm, t, sizeof(tm));
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

        t->c[0].Sw = 2.*t->c[0].L*(t->c[0].h + t->c[0].w)/1000.;
        t->c[1].Sw = 2.*t->c[1].L*(t->c[1].h + t->c[1].w)/1000.;
        //t->I = Iid1(t->H, t->i[p].R, t->c[p].R, t->m.B, t->m.Mu, t->c[p].N);
        t->I = Iid(t->H, t->m.R, t->c[1].R, t->c[p].U, t->m.Mu, t->c[p].N, FR);
        Ks = t->Ks[0]*t->c[0].s + t->Ks[1];
        t->c[0].Rz = t->c[0].L*t->c[0].C*Ks/t->c[0].s;
        Ks = t->Ks[0]*t->c[1].s + t->Ks[1];
        t->c[1].Rz = t->c[1].L*t->c[1].C*Ks/t->c[1].s;

        if(t->c[0].U > t->c[1].U)   t->NN = t->c[0].N/t->c[1].N;
        else                        t->NN = t->c[1].N/t->c[0].N;

        if(t->c[p].U < t->c[!p].U)  t->c[!p].U = t->NN*t->c[p].U;
        else                        t->c[!p].U = t->c[p].U/t->NN;


        return 0;
    } else return 1;

}


void trans_param(TRANS *t, double PT, int p, double FR, int type)
{
    struct TRANS tm;
    int i;
    double SP, Ks, R1;
    fp sq, ln;
    fp1 vou;

    if(type == Round_3phase)
    {
        sq = SQ; ln = LN; vou = Vou;
    }
    else if(type == Round_120_3phase)
    {
        sq = SQ1; ln = LN1; vou = Vou2;
    }
    else if(type == Tore_1phase)
    {

    }
    R1 = t->m.R;

    if(t->c[0].U > t->c[1].U)   t->NN = t->c[0].N/t->c[1].N;
    else                        t->NN = t->c[1].N/t->c[0].N;

    if(t->c[p].U < t->c[!p].U)  t->c[!p].U = t->NN*t->c[p].U;
    else                        t->c[!p].U = t->c[p].U/t->NN;

    t->i[0].R = t->m.R + t->i[0].T;
    t->i[1].R = t->c[0].R + t->i[1].T;

    t->c[0].h = (t->c[0].R - t->i[0].R)/t->c[0].Nr;
    t->c[0].w1 = t->H/t->c[0].Nh;

    //Spiral calc
    SP = (ln)(t->i[0].R, R1)*(t->c[0].Nh+1)/t->H;
    t->c[0].w = (ln)(t->i[0].R, R1)*sqrt(1/(1+SP*SP));
    t->c[0].s = t->c[0].h*t->c[0].w;

    //SP = 2*PI*t->i[0].R*(t->c[0].Nh+1)/t->H;
    //t->c[0].w = 2*PI*t->i[0].R*sqrt(1/(1+SP*SP));
    //t->c[0].s = t->c[0].h*t->c[0].w;
    t->c[0].s1 = (t->c[0].h - 2*t->c[0].T)*(t->c[0].w - 2*t->c[0].T);

    t->c[1].h = (t->c[1].R - t->i[1].R)/t->c[1].Nr;
    t->c[1].w1 = t->H/t->c[1].Nh;

    //Spiral calc
    SP = (ln)(t->i[1].R, 0)*(t->c[1].Nh+1)/t->H;
    t->c[1].w = (ln)(t->i[1].R, 0)*sqrt(1/(1+SP*SP));
    t->c[1].s = t->c[1].h*t->c[1].w;

    //SP = 2*PI*t->i[1].R*(t->c[1].Nh+1)/t->H;
    //t->c[1].w = 2*PI*t->i[1].R*sqrt(1/(1+SP*SP));
    //t->c[1].s = t->c[1].h*t->c[1].w;
    t->c[1].s1 = (t->c[1].h - 2*t->c[1].T)*(t->c[1].w - 2*t->c[1].T);

    t->c[0].L = LnS(t->c[0].Nh, t->c[0].Nr, t->c[0].h, t->i[0].R, t->c[0].w, ln);
    t->c[1].L = LnS(t->c[1].Nh, t->c[1].Nr, t->c[1].h, t->i[1].R, t->c[1].w, ln);

    if(type == Round_3phase || type == Round_120_3phase)  t->m.S = Sou(t->m.R, t->c[1].R, t->W);
    else        t->m.S = Sou1(t->H, t->m.R, t->c[1].R);

    t->i[0].S = Sup(t->m.R, t->i[0].R);
    t->c[0].S = Sup(t->i[0].R, t->c[0].R);
    t->i[1].S = Sup(t->c[0].R, t->i[1].R);
    t->c[1].S = Sup(t->i[1].R, t->c[1].R) + Ssd(t->H, t->c[1].R);

    t->S = 3*(t->m.S + t->i[0].S + t->c[0].S + t->i[1].S + t->c[1].S);
    t->T = PT/(10*t->S);

    if(type == Round_3phase || type == Round_120_3phase)  t->m.V = Vol(t->H, 0, (sq)(t->m.R, R1)) + (vou)((sq)(t->m.R, R1), t->m.R, t->c[1].R, t->W);
    else        t->m.V = Vol(t->H, 0, (sq)(t->m.R, R1))*2. + Vou1(t->m.R, t->c[1].R);
    t->i[0].V = Vol(t->H, (sq)(t->m.R, R1)   , (sq)(t->i[0].R, R1));
    t->c[0].V = Vol(t->H, (sq)(t->i[0].R, R1), (sq)(t->c[0].R, R1));
    t->i[1].V = Vol(t->H, (sq)(t->c[0].R, R1), (sq)(t->i[1].R, R1));
    t->c[1].V = Vol(t->H, (sq)(t->i[1].R, R1), (sq)(t->c[1].R, R1));

    t->m.M    = t->m.V*t->m.D*(t->m.Sfc + 1)/(2*1000);        //Mass steel
    t->i[0].M = t->i[0].V*t->i[0].D/1000.; //Mass of insulation
    t->c[0].M = t->c[0].V*t->c[0].D/1000.; //Mass coil
    t->i[1].M = t->i[1].V*t->i[1].D/1000.; //Mass of insulation
    t->c[1].M = t->c[1].V*t->c[1].D/1000.; //Mass 1coil

    if(type == Round_3phase || type == Round_120_3phase)  t->M = 3.*(t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);
    else        t->M =    (t->m.M + t->i[0].M + t->c[0].M + t->i[1].M + t->c[1].M);

    t->m.P = 3.*t->m.M*t->m.Lc*t->m.Sfc;  //Loss in steel
    Ks = t->Ks[0]*t->c[0].s + t->Ks[1];
    t->c[0].P = Pw2(t->c[0].L, t->c[0].h, t->c[0].w, t->c[0].C, t->c[0].I, t->c[0].T, Ks);
    Ks = t->Ks[0]*t->c[1].s + t->Ks[1];
    t->c[1].P = Pw2(t->c[1].L, t->c[1].h, t->c[1].w, t->c[1].C, t->c[1].I, t->c[1].T, Ks);

    if(type == Round_3phase || type == Round_120_3phase)  t->P = t->m.P + 3*(t->c[0].P + t->c[1].P);
    else        t->P =   (t->m.P + t->c[0].P + t->c[1].P);


    t->c[0].Sw = 2.*t->c[0].L*(t->c[0].h + t->c[0].w)/1000.;
    t->c[1].Sw = 2.*t->c[1].L*(t->c[1].h + t->c[1].w)/1000.;
    //t->I = Iid1(t->H, t->i[p].R, t->c[p].R, t->m.B, t->m.Mu, t->c[p].N);
    t->I = Iid(t->H, t->m.R, t->c[1].R, t->c[p].U, t->m.Mu, t->c[p].N, FR);
    Ks = t->Ks[0]*t->c[0].s + t->Ks[1];
    t->c[0].Rz = t->c[0].L*t->c[0].C*Ks/t->c[0].s;
    Ks = t->Ks[0]*t->c[1].s + t->Ks[1];
    t->c[1].Rz = t->c[1].L*t->c[1].C*Ks/t->c[1].s;

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
    printf("КПД %f %%\n", (t->PW*1000-t->P)*100/(t->PW*1000));
    printf("Масса %f кг\n", t->M);
    printf("Площадь поверхности %f м**2\n", t->S);
    printf("Температура над окружением %f K\n", t->T);
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
    printf("Сечение провода %f мм**2  без изоляции %f мм**2\n", t->c[0].s, t->c[0].s1);
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
    printf("Сечение провода %f мм**2  без изоляции %f мм**2\n", t->c[1].s, t->c[1].s1);
    printf("Площадь поверхности провода %f м**2\n", t->c[1].Sw);
    printf("Радиус внутренний %f мм внешний %f мм высота %f мм\n", t->i[1].R, t->c[1].R, t->H);
    printf("Число витков по высоте %f по радиусу %f \n", t->c[1].Nh, t->c[1].Nr);
    printf("Высота сечения провода %f ширина %f \n", t->c[1].h, t->c[1].w);

    printf("H = %f R1 = %f R2 = %f R3 = %f R4 = %f R5 = %f W = %f\n", t->H, t->m.R, t->i[0].R, t->c[0].R, t->i[1].R, t->c[1].R, t->W);
}

void trans(void)
{
    int i, n = 7;
    TRANS tr[n], t;
    //COIL coil[9];
    double LOST, POWER = 23;
    double FR = 50;
    int p = 1;

    //3 Phase power treasformer
    if(POWER == 25){
        i=0;
        LOST = 400;
        tr[i] = (TRANS) {.PW = POWER, .W = 2, .H = 100, .Hp[0] = 300, .Hp[1] = 1, .Ks[0] = 0.000155, .Ks[1] = 0.999845 };
        tr[i].m = (MCORE) {.D = 7.65, .Mu = 10000., .B = 1.85, .Lc = 1., .Sfc = 0.95, .R = 40, .Rp[0] = 60, .Rp[1] = 1 };
        tr[i].i[0] = (INS) {.D = 3.26, .T = 1.};
        tr[i].c[0] = (COIL) {.U = 400,   .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 10,  .sp[0] = 40, .sp[1] = 1, .Nr = 0}; //25KW
        tr[i].c[0].I = POWER*1000*sqrt(3.)/tr[i].c[0].U/3.;
        tr[i].i[1] = (INS) {.D = 3.26, .T = 2.};
        tr[i].c[1] = (COIL) {.U = 10000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 0.5, .sp[0] = 4,  .sp[1] = 0.1, .Nr = 20}; //25KW
        tr[i].c[1].I = POWER*1000/tr[i].c[1].U/3.;
    }
    else if(POWER == 50){
        i=1;
        LOST = 600;
        tr[i] = (TRANS) {.PW = POWER, .W = 3, .H = 100, .Hp[0] = 300, .Hp[1] = 1, .Ks[0] = 0.000155, .Ks[1] = 0.999845 };
        tr[i].m = (MCORE) {.D = 7.65, .Mu = 10000., .B = 1.85, .Lc = 1., .Sfc = 0.95, .R = 50, .Rp[0] = 70, .Rp[1] = 1 };
        tr[i].i[0] = (INS) {.D = 3.26, .T = 1.};
        tr[i].c[0] = (COIL) {.U = 400,   .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 50,  .sp[0] = 100, .sp[1] = 1, .Nr = 0}; //50KW
        tr[i].c[0].I = POWER*1000*sqrt(3.)/tr[i].c[0].U/3.;
        tr[i].i[1] = (INS) {.D = 3.26, .T = 2.};
        tr[i].c[1] = (COIL) {.U = 10000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 1, .sp[0] = 6,  .sp[1] = 0.2, .Nr = 20}; //50KW
        tr[i].c[1].I = POWER*1000/tr[i].c[1].U/3.;
    }
    else if(POWER == 100){
        i=2;
        LOST = 980;
        tr[i] = (TRANS) {.PW = POWER, .W = 3, .H = 200, .Hp[0] = 400, .Hp[1] = 1, .Ks[0] = 0.000155, .Ks[1] = 0.999845 };
        tr[i].m = (MCORE) {.D = 7.65, .Mu = 10000., .B = 1.9, .Lc = 1.55, .Sfc = 0.95, .R = 65, .Rp[0] = 80, .Rp[1] = 1 }; //1.9 1.55
        tr[i].i[0] = (INS) {.D = 3.26, .T = 2.5};
        tr[i].c[0] = (COIL) {.U = 400,   .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 100,  .sp[0] = 200, .sp[1] = 1, .Nr = 0}; //100KW
        tr[i].c[0].I = POWER*1000*sqrt(3.)/tr[i].c[0].U/3.;
        tr[i].i[1] = (INS) {.D = 3.26, .T = 2.5};
        //tr[i].c[1] = (COIL) {.U = 10000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 2, .sp[0] = 8,  .sp[1] = 0.2, .Nr = 20}; //100KW
        tr[i].c[1] = (COIL) {.U = 315,   .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 100,  .sp[0] = 200, .sp[1] = 1, .Nr = 0};//100KW
        tr[i].c[1].I = POWER*1000/tr[i].c[1].U/3.;
    }
    else if(POWER == 250){
         i=3;

    }
    else if(POWER == 500){
        i=4;
        LOST = 3100;
        tr[i] = (TRANS) {.PW = POWER, .W = 5, .H = 300, .Hp[0] = 500, .Hp[1] = 1, .Ks[0] = 0.000155, .Ks[1] = 0.999845 };
        tr[i].m = (MCORE) {.D = 7.65, .Mu = 10000., .B = 1.85, .Lc = 1., .Sfc = 0.95, .R = 100, .Rp[0] = 120, .Rp[1] = 1 };
        tr[i].i[0] = (INS) {.D = 3.26, .T = 1.};
        tr[i].c[0] = (COIL) {.U = 400,   .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 600,  .sp[0] = 800, .sp[1] = 1, .Nr = 0}; //500KW
        tr[i].c[0].I = POWER*1000*sqrt(3.)/tr[i].c[0].U/3.;
        tr[i].i[1] = (INS) {.D = 3.26, .T = 2.};
        tr[i].c[1] = (COIL) {.U = 10000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 20, .sp[0] = 50,  .sp[1] = 1, .Nr = 20}; //500KW
        //tr[i].c[1] = (COIL) {.U = 315,   .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 600,  .sp[0] = 800, .sp[1] = 1, .Nr = 0};//500KW
        tr[i].c[1].I = POWER*1000/tr[i].c[1].U/3.;
    }
    else if(POWER == 1000){
        i=5;
        LOST = 6000;
        tr[i] = (TRANS) {.PW = POWER, .W = 5, .H = 400, .Hp[0] = 600, .Hp[1] = 1, .Ks[0] = 0.000155, .Ks[1] = 0.999845 };
        tr[i].m = (MCORE) {.D = 7.65, .Mu = 10000., .B = 1.85, .Lc = 1., .Sfc = 0.95, .R = 120, .Rp[0] = 140, .Rp[1] = 1 };
        tr[i].i[0] = (INS) {.D = 3.26, .T = 1.};
        tr[i].c[0] = (COIL) {.U = 400,   .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 1200,  .sp[0] = 1300, .sp[1] = 2, .Nr = 0}; //1000KW
        tr[i].c[0].I = POWER*1000*sqrt(3.)/tr[i].c[0].U/3.;
        tr[i].i[1] = (INS) {.D = 3.26, .T = 2.};
        tr[i].c[1] = (COIL) {.U = 10000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 35, .sp[0] = 50,  .sp[1] = 1, .Nr = 20}; //1000KW
        tr[i].c[1].I = POWER*1000/tr[i].c[1].U/3.;
    }
        else if(POWER == 2000){
            i=6;
            LOST = 12000;
            tr[i] = (TRANS) {.PW = POWER, .W = 5, .H = 400, .Hp[0] = 600, .Hp[1] = 1, .Ks[0] = 0.000155, .Ks[1] = 0.999845 };
            tr[i].m = (MCORE) {.D = 7.65, .Mu = 10000., .B = 1.85, .Lc = 1., .Sfc = 0.95, .R = 140, .Rp[0] = 160, .Rp[1] = 1 };
            tr[i].i[0] = (INS) {.D = 3.26, .T = 1.};
            tr[i].c[0] = (COIL) {.U = 400,   .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 1500,  .sp[0] = 2300, .sp[1] = 2, .Nr = 0}; //2000KW
            tr[i].c[0].I = POWER*1000*sqrt(3.)/tr[i].c[0].U/3.;
            tr[i].i[1] = (INS) {.D = 3.26, .T = 2.};
            tr[i].c[1] = (COIL) {.U = 10000, .T = 0.06,.C = 0.0282, .D = 2.6989, .s = 60, .sp[0] = 85,  .sp[1] = 1, .Nr = 20}; //2000KW
            tr[i].c[1].I = POWER*1000/tr[i].c[1].U/3.;
        }


    //if(trans_optim(&(tr[i]), LOST, p, FR, Round_120_3phase)) { printf("No any result!!!\n"); return; }

    //One phase tor transformer
            if(POWER == 23){
                i=10;
                LOST = POWER*1000.*(1 - 0.8);
                tr[i] = (TRANS) {.PW = POWER, .W = 5, .H = 400, .Hp[0] = 600, .Hp[1] = 1, .Ks[0] = 0.000155, .Ks[1] = 0.999845 };
                tr[i].m = (MCORE) {.D = 7.65, .Mu = 10000., .B = 2.7, .Lc = 3, .Sfc = 0.95, .R = 100, .Rp[0] = 200, .Rp[1] = 0.1, .R1 = 50, .R1p[0] = 100, .R1p[1] = 0.1 };
                //tr[i].i[0] = (INS) {.D = 3.26, .T = 1.};
                tr[i].c[0] = (COIL) {.U = 230,   .T = 0.06,.C = 0.0282, .D = 2.6989, .R = 20.370513, .Rp[0] = 20.370513, .Rp[1] = 0.1, .s = 2.3, .N = 100, .Np[0] = 600, .Np[1] = 1};
                tr[i].c[0].I = 100.;
                //tr[i].i[1] = (INS) {.D = 3.26, .T = 2.};
                tr[i].c[1] = (COIL) {.U = 3, .T = 0.06,.C = 0.0282, .D = 2.6989, .N = 1, .Np[0] = 15, .Np[1] = 1}; //8.92 //0.0175
                //tr[i].c[1].I = POWER*1000/tr[i].c[1].U/3.;
            }

                if(trans_optim_tor(&(tr[i]), POWER, 10 , FR)) { printf("No any result!!!\n"); return; }



    //10000kV
    /*
    tr[2].c[1].R = 133;
    tr[2].c[1].Nh = 83;
    tr[2].c[1].Nr = 20;
    tr[2].c[1].N = tr[2].c[1].Nh*tr[2].c[1].Nr;
    tr[2].c[1].I = POWER*1000/tr[2].c[1].U/3.;
    */
    //12000kV
    /*
    tr[2].c[1].R = 129;
    tr[2].c[1].Nh = 81;
    tr[2].c[1].Nr = 24;
    tr[2].c[1].N = tr[2].c[1].Nh*tr[2].c[1].Nr;
    tr[2].c[1].I = POWER*1000/tr[2].c[1].U/3.;
    */
    /*
    tr[2].H = 256;
    tr[2].m.R = 77;

    //400V
    tr[2].c[0].R = 114;
    tr[2].c[0].Nh = 56;
    tr[2].c[0].Nr = 1;
    tr[2].c[0].N = tr[2].c[0].Nh*tr[2].c[0].Nr;
    tr[2].c[0].I = POWER*1000*sqrt(3.)/tr[2].c[0].U/3.;
    //315V
    tr[2].c[1].R = 143;
    tr[2].c[1].Nh = 44;
    tr[2].c[1].Nr = 1;
    tr[2].c[1].N = tr[2].c[1].Nh*tr[2].c[1].Nr;
    tr[2].c[1].I = POWER*1000/tr[2].c[1].U/3.;


    trans_param(&(tr[i]), LOST, p, FR, Round_120_3phase);
    */
    print_result(&(tr[i]), p);

}

