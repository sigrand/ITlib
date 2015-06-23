#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libit/types.h"
#include "./engine.h"

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
    en->ST.CR.R2 = en->ST.CR.R1*cos(al)/(1 + (2*PI*(2 + cos(al)))/(en->ST.N*(en->Ph-1)));
    en->ST.CR.a = 2*PI*en->ST.CR.R2/(en->ST.N*(en->Ph-1));
    en->ST.CR.w = 2*en->ST.CR.a;
    en->ST.CR.p = 2*PI*en->ST.CR.R2/en->ST.N;

    //Rotor
    en->RT.CR.R1 = en->ST.CR.R2 - en->h;
    en->RT.CR.a = 2*PI*en->RT.CR.R1/(en->RT.N*en->Ph);
    en->RT.CR.R2 = en->RT.CR.R1 - 3*en->RT.CR.a;
    en->RT.CR.w = 3*en->RT.CR.a;
    en->RT.CR.p = 2*PI*en->RT.CR.R1/en->RT.N;
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

void engine(void)
{
    double vl, ms, rz;
    ENGINE en;
    en.L = 50;
    en.R = 150;
    en.Ph = 6;
    en.Pn = 3;
    en.ST.N = en.Pn*en.Ph;
    en.RT.N = en.Pn*(en.Ph-1);
    en.h = 0.4;
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
    printf("ST.R1 = %f ST.R2 = %f ST.w = %f ST.p = %f ST.l = %f\n", en.ST.CR.R1, en.ST.CR.R2, en.ST.CR.w, en.ST.CR.p, 2*PI*en.ST.CR.R2);
    printf("RT.R1 = %f RT.R2 = %f RT.w = %f RT.p = %f RT.l = %f\n", en.RT.CR.R1, en.RT.CR.R2, en.RT.CR.w, en.RT.CR.p, 2*PI*en.RT.CR.R1);

    vl = volume_ENGINE(&en);
    printf("Engine voulue = %f\n", vl);

    ms = mass_ENGINE(&en);
    printf("Engine mass = %f\n", ms);

    /*
    rz = resistance_COIL(&en.ST.CL);
    printf("Resistance = %f   %f\n", rz, rz*en.ST.N);
    */

    //Conductivity

}
