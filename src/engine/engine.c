#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libit/types.h"
#include "./engine.h"

inline double BH(double H)
{
    int i;
    double a, b;
    for(i=0; i < 9; i++){
        //printf("H = %f HB0 = %f HB1 = %f\n",H, HB[i][0], HB[i+1][0]);
        if(H > HB[i][0] && H <= HB[i+1][0]){
            a = (HB[i+1][1] - HB[i][1])/(HB[i+1][0] - HB[i][0]);
            b = HB[i+1][1] - a*HB[i+1][0];
            return a*H + b;
        }
    }
    //a = (HB[0][1] - HB[1][1])/(HB[0][0] - HB[1][0]);
    //b = HB[1][1] - a*HB[1][0];
    //return a*H + b;
    //return 0;
}


double resistance_COIL(COIL *cl)
{
    double s, l, h;
    s = PI*cl->R1*cl->R1*cl->FF/(2*cl->N);
    h = s*cl->N/(PI*cl->R1);
    printf("h = %f s = %f d = %f\n", h, s, 2*sqrt(s/PI));
    l = (2*cl->L + 2*(cl->R2 - cl->R1))*0.001;
    return cl->C*l*cl->N/s;
}

//Caculation geometry paraments
double get_sizes(ENGINE *en)
{
    double al;
    //Stator
    al = 2*PI/(2*en->ST.N);
    en->ST.CR.al = al*180/PI;
    en->ST.CR.R2 = en->ST.CR.R1*cos(al)/(1 + (2*PI*(2 + cos(al)))/(en->ST.N*(en->Ph-1)));
    en->ST.CR.a = 2*PI*en->ST.CR.R2/(en->ST.N*(en->Ph-1));
    en->ST.CR.w = 2*en->ST.CR.a;
    en->ST.CR.p = 2*PI*en->ST.CR.R2/en->ST.N;
    en->ST.CR.h = 2*en->ST.CR.a;

    //Rotor
    en->RT.CR.al = 360/(2*en->RT.N);
    en->RT.CR.R1 = en->ST.CR.R2 - en->h;
    en->RT.CR.a = 2*PI*en->RT.CR.R1/(en->RT.N*en->Ph);
    en->RT.CR.R2 = en->RT.CR.R1 - 3*en->RT.CR.a;
    en->RT.CR.w = 3*en->RT.CR.a;
    en->RT.CR.p = 2*PI*en->RT.CR.R1/en->RT.N;
    en->RT.CR.h = 2*en->RT.CR.a;
}

double volume_COIL(ENGINE *en)
{
    double R = (en->ST.CR.p - en->ST.CR.w)/2;
    double S = en->ST.CR.w*R;
    double L = 2*(en->ST.CR.w + en->L);
    en->ST.CL.S = en->ST.CR.w*R;
    return (S*L + PI*R*R*en->ST.CR.w)*0.000001;
}

double volume_CORE_STATOR(ENGINE *en)
{
    CORE *cr = &en->ST.CR;
    double sg, R1, R2;
    R1 = cr->R2 + cr->w;
    R2 = cr->R2;
    sg = PI*(R1*R1 - R2*R2)*(cr->p - cr->w)/cr->p;
    R1 = cr->R1;
    R2 = cr->R2;
    return (PI*(R1*R1 - R2*R2) - sg)*en->L*0.000001;
}

double volume_CORE_ROTOR(ENGINE *en)
{
    CORE *cr = &en->RT.CR;
    double sg, R1, R2;
    R1 = cr->R1;
    R2 = cr->R1 - en->ST.CR.w;
    sg = PI*(R1*R1 - R2*R2)*(cr->p - cr->w)/cr->p;
    R1 = cr->R1;
    R2 = cr->R2;
    return (PI*(R1*R1 - R2*R2) - sg)*en->L*0.000001;
}

//Calculate volume of engine
double volume_ENGINE(ENGINE *en)
{

    en->ST.CR.V = volume_CORE_STATOR(en);
    printf("Stator core volume = %f\n", en->ST.CR.V);

    en->RT.CR.V = volume_CORE_ROTOR(en);
    printf("Rotor core volume = %f\n", en->RT.CR.V);

    //Coils
    en->ST.CL.V = volume_COIL(en);
    printf("Coil volume = %f total = %f S = %f\n", en->ST.CL.V, en->ST.N*en->ST.CL.V, en->ST.CL.S);
    /*
    vcl = en->ST.N*volume_COIL(&en->ST.CL);

    printf("Coil volume = %f\n", volume_COIL(&en->ST.CL));

    vs = volume_STATOR(&en->ST);

    printf("Stator volume = %f\n", vs);

    vr = volume_ROTOR(&en->RT);

    printf("Rotor volume = %f\n", vr);
    */

    return en->ST.CR.V + en->RT.CR.V + en->ST.CL.V*en->ST.N;
}

//Calculate mass of engine
double mass_ENGINE(ENGINE *en)
{
    //double ms, mr, mcr, mcl, R1, R2;

    en->ST.CR.M = en->ST.CR.V*en->ST.CR.D;
    printf("Stator core mass = %f\n", en->ST.CR.M);

    en->RT.CR.M = en->RT.CR.V*en->RT.CR.D;
    printf("Rotor core mass = %f\n", en->RT.CR.M);

    en->ST.CL.M = en->ST.CL.V*en->ST.CL.D;
    printf("Coil mass = %f total = %f\n", en->ST.CL.M, en->ST.CL.M*en->ST.N);

    /*
    printf("Core mass = %f  %f\n", volume_CORE_STATOR(&en->ST.CR)*en->ST.CR.D, mcr);

    //Coils
    mcl = en->ST.N*volume_COIL(&en->ST.CL)*en->ST.CL.D;

    printf("Coil mass = %f  %f\n", volume_COIL(&en->ST.CL)*en->ST.CL.D, mcl);

    ms = volume_STATOR(&en->ST)*en->ST.D;

    printf("Stator mass = %f\n", ms);

    mr = volume_ROTOR(&en->RT)*en->ST.D;

    printf("Rotor mass = %f\n", mr);
    */

    return en->ST.CR.M + en->RT.CR.M + en->ST.CL.M*en->ST.N;

}

void B_calc(ENGINE *en)
{
    int i = 0;
    double Hst[3], Hrt[3], Hsy[2], Hry[2], Hh[3];
    double Bst[3], Brt[3], Bsy[2], Bry[2], Bh[3];
    double Lst, Lrt, Lsy, Lry, Lh;
    double Sst, Srt, Ssy, Sry, Sh;
    double S[3];
    double I = 20; //Current density
    double d0 = 100, d1 = 100, l0, l1;

    //Set up lentgh
    Lst = (en->ST.CR.h + en->ST.CR.a/2)*0.001;
    Lrt = (en->RT.CR.h + en->RT.CR.a/2)*0.001;
    Lh = en->h*0.001;
    Lsy = en->ST.CR.p*0.001;
    Lry = en->ST.CR.p*0.001;

    //Set up square
    Sst = en->ST.CR.w*en->L;
    Sst = en->ST.CR.w*en->L;
    Sh = en->ST.CR.w*en->L;
    Ssy = en->ST.CR.a*en->L;
    Sry = en->RT.CR.a*en->L;
    S[0] = Sst; S[1] = Sst;  S[2] = Sst;

    l0 = Lst + Lrt + Lh + Lsy + Lry;
    printf("Lst = %f Lrt = %f Lh = %f Lsy = %f Lry = %f l = %f Ssy = %f\n",Lst, Lrt, Lh, Lsy, Lry, l0, Ssy);

    for(i=0, Hsy[0] = 100, Hsy[1] = 100; abs(d0) > 2 && abs(d1) > 2; i++){

        Hry[0] = Hsy[0];
        Hry[1] = Hsy[1];
        Hst[1] = (Hsy[0] + Hsy[1])*Ssy/S[1];
        Hrt[1] = (Hry[0] + Hry[1])*Ssy/S[1];
        Hst[2] = Hsy[1]*Ssy/S[2];
        Hrt[2] = Hry[1]*Ssy/S[2];
        Hst[0] = Hsy[0]*Ssy/S[0];
        Hrt[0] = Hry[0]*Ssy/S[0];
        Hh[1] =  BH(Hst[1])/MU;
        Hh[0] =  BH(Hst[0])/MU;
        Hh[2] =  BH(Hst[2])/MU;

        d0 = (Lst*(Hst[0] + Hst[1]) + Lrt*(Hrt[0] + Hrt[1]) + Lh*(Hh[0] + Hh[1]) + Lsy*Hsy[0] + Lry*Hry[0]) - en->ST.CL.S*I;
        d1 = (Lst*(Hst[1] + Hst[2]) + Lrt*(Hrt[1] + Hrt[2]) + Lh*(Hh[1] + Hh[2]) + Lsy*Hsy[1] + Lry*Hry[1]) - en->ST.CL.S*I;

        if(d0 > 0) Hsy[0]--;
        else Hsy[0]++;

        if(d1 > 0) Hsy[1]--;
        else Hsy[1]++;

        //printf("I = %f d0 = %f d1 = %f H0 = %f H1 = %f B = %f\n", en->ST.CL.S*I, d0, d1, Hsy[0], Hsy[1], BH(Hst[1]));
        //if(i > 10) break;
    }
    printf("I = %f d0 = %f d1 = %f H0 = %f H1 = %f H2 = %f B = %f\n", en->ST.CL.S*I, d0, d1, Hsy[0], Hsy[1], Hh[1], BH(Hst[1]));
}

double get_square(ENGINE *en, double Ts, double Tr)
{
    double S, d, di;
    if(Tr < Ts) {
        d = Ts - Tr;
        di = en->RT.CR.w - d;
        S = di > en->ST.CR.w ? en->ST.CR.w : di;
        //printf("d = %f di = %f S = %f\n", d, di, S);
    } else {
        d = Tr - Ts;
        S = en->ST.CR.w - d;
    }
    return S;
}

void position_calc(ENGINE *en)
{
    double Ts[6], Tr[6];
    double S[3];
    double d, di;
    //Set up
    Ts[0] = -en->ST.CR.w/2 - en->ST.CR.p;
    Ts[1] =  en->ST.CR.w/2 - en->ST.CR.p;
    Ts[2] = -en->ST.CR.w/2;
    Ts[3] =  en->ST.CR.w/2;
    Ts[4] = -en->ST.CR.w/2 + en->ST.CR.p;
    Ts[5] =  en->ST.CR.w/2 + en->ST.CR.p;

    Tr[0] = -en->RT.CR.w/2 - en->RT.CR.p;
    Tr[1] =  en->RT.CR.w/2 - en->RT.CR.p;
    Tr[2] = -en->RT.CR.w/2;
    Tr[3] =  en->RT.CR.w/2;
    Tr[4] = -en->RT.CR.w/2 + en->RT.CR.p;
    Tr[5] =  en->RT.CR.w/2 + en->RT.CR.p;

    printf("Ts[0] = %f Ts[1] = %f Ts[2] = %f Ts[3] = %f Ts[4] = %f Ts[5] = %f \n", Ts[0],Ts[1],Ts[2],Ts[3],Ts[4],Ts[5]);
    printf("Tr[0] = %f Tr[1] = %f Tr[2] = %f Tr[3] = %f Tr[4] = %f Tr[5] = %f \n", Tr[0],Tr[1],Tr[2],Tr[3],Tr[4],Tr[5]);
    /*
    if(Tr[0] < Ts[0]) {
        d = Ts[0] - Tr[0];
        di = en->RT.CR.p - d;
        S[0] = di > en->ST.CR.w ? en->ST.CR.w : di;
    } else {
        d = Tr[0] - Ts[0];
        S[0] = en->ST.CR.p - d;
    }
    */
    S[0] = get_square(en, Ts[0], Tr[0]);
    S[1] = get_square(en, Ts[2], Tr[2]);
    S[2] = get_square(en, Ts[4], Tr[4]);
    printf("S[0] = %f S[1] = %f S[2] = %f \n", S[0], S[1], S[2]);
}

void engine(void)
{
    double vl, ms, rz;
    ENGINE en;
    en.L = 50;
    en.R = 150;
    en.Ph = 6;
    en.Pn = 4;
    en.ST.N = en.Pn*en.Ph;
    en.RT.N = en.Pn*(en.Ph-1);
    en.h = 0.5;
    en.EX = 1.5;

    en.ST.D = 1.830;        //Magnesium 1.830
    en.RT.D = 1.830;        //Magnesium
    en.ST.CL.D = 2.6989;    //Aluminium
    //en.ST.CL.D = 8.92;    //Copper
    en.ST.CL.C = 0.0282;    //Conductivity
    en.ST.CL.N = 60;        //The numner of turns
    en.ST.CL.FF = 0.9;      //60% filling
    en.ST.CR.D = 7.65;      //Magnetic steel
    en.RT.CR.D = 7.65;      //Magnetic steel

    //Calculation core size
    en.ST.CR.R1 = en.R;

    get_sizes(&en);

    printf("Phaese = %f Groups = %f  ST.N = %f RT.N = %f \n", en.Ph, en.Pn, en.ST.N, en.RT.N);
    printf("ST.R1 = %f ST.R2 = %f ST.w = %f ST.p = %f al = %f ST.l = %f\n", en.ST.CR.R1, en.ST.CR.R2, en.ST.CR.w, en.ST.CR.p, en.ST.CR.al, 2*PI*en.ST.CR.R2);
    printf("RT.R1 = %f RT.R2 = %f RT.w = %f RT.p = %f al = %f RT.l = %f\n", en.RT.CR.R1, en.RT.CR.R2, en.RT.CR.w, en.RT.CR.p, en.RT.CR.al, 2*PI*en.RT.CR.R1);

    vl = volume_ENGINE(&en);
    printf("Engine voulue = %f\n", vl);

    ms = mass_ENGINE(&en);
    printf("Engine mass = %f\n", ms);

    B_calc(&en);

    position_calc(&en);

    /*
    rz = resistance_COIL(&en.ST.CL);
    printf("Resistance = %f   %f\n", rz, rz*en.ST.N);
    */

    //Conductivity

}
