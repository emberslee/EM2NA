/* IdealRNA.hpp - Idealized RNA nucleotide conformations */
#include <string>
#include <map>
#include <vector>

using namespace std;

map<string, map<string,vector<double> > >parse_ideal_rna()
{
    vector<double> tmp(3,0);
    map<string, map<string,vector<double> > >ideal_rna;

    map<string,vector<double> >P;
    P[" O3'"]=tmp; P[" O3'"][0]=   4.210; P[" O3'"][1]=   6.832; P[" O3'"][2]=  -0.518;
    P[" P  "]=tmp; P[" P  "][0]=   5.130; P[" P  "][1]=   7.667; P[" P  "][2]=  -1.527;
    P[" OP1"]=tmp; P[" OP1"][0]=   5.914; P[" OP1"][1]=   8.669; P[" OP1"][2]=  -0.770;
    P[" OP2"]=tmp; P[" OP2"][0]=   4.303; P[" OP2"][1]=   8.192; P[" OP2"][2]=  -2.635;
    P[" O5'"]=tmp; P[" O5'"][0]=   6.107; P[" O5'"][1]=   6.526; P[" O5'"][2]=  -2.080;
    ideal_rna["P"]=P;
    map<string,vector<double> >().swap(P);

    map<string,vector<double> >A;
    A[" P  "]=tmp; A[" P  "][0]=  -0.356; A[" P  "][1]=   9.218; A[" P  "][2]=   1.848;
    A[" OP1"]=tmp; A[" OP1"][0]=  -0.311; A[" OP1"][1]=  10.489; A[" OP1"][2]=   2.605;
    A[" OP2"]=tmp; A[" OP2"][0]=  -1.334; A[" OP2"][1]=   9.156; A[" OP2"][2]=   0.740;
    A[" O5'"]=tmp; A[" O5'"][0]=   1.105; A[" O5'"][1]=   8.869; A[" O5'"][2]=   1.295;
    A[" C5'"]=tmp; A[" C5'"][0]=   2.021; A[" C5'"][1]=   8.156; A[" C5'"][2]=   2.146;
    A[" C4'"]=tmp; A[" C4'"][0]=   2.726; A[" C4'"][1]=   7.072; A[" C4'"][2]=   1.355;
    A[" O4'"]=tmp; A[" O4'"][0]=   1.986; A[" O4'"][1]=   5.817; A[" O4'"][2]=   1.352;
    A[" C3'"]=tmp; A[" C3'"][0]=   2.952; A[" C3'"][1]=   7.370; A[" C3'"][2]=  -0.127;
    A[" O3'"]=tmp; A[" O3'"][0]=   4.210; A[" O3'"][1]=   6.832; A[" O3'"][2]=  -0.518;
    A[" C2'"]=tmp; A[" C2'"][0]=   1.848; A[" C2'"][1]=   6.598; A[" C2'"][2]=  -0.850;
    A[" C1'"]=tmp; A[" C1'"][0]=   1.913; A[" C1'"][1]=   5.344; A[" C1'"][2]=   0.016;
    A[" N9 "]=tmp; A[" N9 "][0]=   0.711; A[" N9 "][1]=   4.472; A[" N9 "][2]=  -0.101;
    A[" C8 "]=tmp; A[" C8 "][0]=  -0.589; A[" C8 "][1]=   4.841; A[" C8 "][2]=  -0.292;
    A[" N7 "]=tmp; A[" N7 "][0]=  -1.415; A[" N7 "][1]=   3.843; A[" N7 "][2]=  -0.354;
    A[" C5 "]=tmp; A[" C5 "][0]=  -0.604; A[" C5 "][1]=   2.728; A[" C5 "][2]=  -0.192;
    A[" C6 "]=tmp; A[" C6 "][0]=  -0.881; A[" C6 "][1]=   1.351; A[" C6 "][2]=  -0.162;
    A[" N6 "]=tmp; A[" N6 "][0]=  -2.113; A[" N6 "][1]=   0.841; A[" N6 "][2]=  -0.301;
    A[" N1 "]=tmp; A[" N1 "][0]=   0.158; A[" N1 "][1]=   0.514; A[" N1 "][2]=   0.016;
    A[" C2 "]=tmp; A[" C2 "][0]=   1.380; A[" C2 "][1]=   1.027; A[" C2 "][2]=   0.154;
    A[" N3 "]=tmp; A[" N3 "][0]=   1.758; A[" N3 "][1]=   2.286; A[" N3 "][2]=   0.143;
    A[" C4 "]=tmp; A[" C4 "][0]=   0.700; A[" C4 "][1]=   3.103; A[" C4 "][2]=  -0.037;
    ideal_rna["  A"]=A;
    map<string,vector<double> >().swap(A);

    map<string,vector<double> >C;
    C[" P  "]=tmp; C[" P  "][0]=   5.130; C[" P  "][1]=   7.667; C[" P  "][2]=  -1.527;
    C[" OP1"]=tmp; C[" OP1"][0]=   5.914; C[" OP1"][1]=   8.669; C[" OP1"][2]=  -0.770;
    C[" OP2"]=tmp; C[" OP2"][0]=   4.303; C[" OP2"][1]=   8.192; C[" OP2"][2]=  -2.635;
    C[" O5'"]=tmp; C[" O5'"][0]=   6.107; C[" O5'"][1]=   6.526; C[" O5'"][2]=  -2.080;
    C[" C5'"]=tmp; C[" C5'"][0]=   6.430; C[" C5'"][1]=   5.410; C[" C5'"][2]=  -1.229;
    C[" C4'"]=tmp; C[" C4'"][0]=   6.362; C[" C4'"][1]=   4.119; C[" C4'"][2]=  -2.020;
    C[" O4'"]=tmp; C[" O4'"][0]=   5.026; C[" O4'"][1]=   3.539; C[" O4'"][2]=  -2.023;
    C[" C3'"]=tmp; C[" C3'"][0]=   6.720; C[" C3'"][1]=   4.227; C[" C3'"][2]=  -3.502;
    C[" O3'"]=tmp; C[" O3'"][0]=   7.422; C[" O3'"][1]=   3.053; C[" O3'"][2]=  -3.893;
    C[" C2'"]=tmp; C[" C2'"][0]=   5.374; C[" C2'"][1]=   4.251; C[" C2'"][2]=  -4.225;
    C[" C1'"]=tmp; C[" C1'"][0]=   4.689; C[" C1'"][1]=   3.199; C[" C1'"][2]=  -3.359;
    C[" N1 "]=tmp; C[" N1 "][0]=   3.204; C[" N1 "][1]=   3.200; C[" N1 "][2]=  -3.476;
    C[" C2 "]=tmp; C[" C2 "][0]=   2.545; C[" C2 "][1]=   1.977; C[" C2 "][2]=  -3.385;
    C[" O2 "]=tmp; C[" O2 "][0]=   3.212; C[" O2 "][1]=   0.950; C[" O2 "][2]=  -3.212;
    C[" N3 "]=tmp; C[" N3 "][0]=   1.191; C[" N3 "][1]=   1.956; C[" N3 "][2]=  -3.490;
    C[" C4 "]=tmp; C[" C4 "][0]=   0.504; C[" C4 "][1]=   3.090; C[" C4 "][2]=  -3.677;
    C[" N4 "]=tmp; C[" N4 "][0]=  -0.815; C[" N4 "][1]=   3.018; C[" N4 "][2]=  -3.773;
    C[" C5 "]=tmp; C[" C5 "][0]=   1.163; C[" C5 "][1]=   4.359; C[" C5 "][2]=  -3.774;
    C[" C6 "]=tmp; C[" C6 "][0]=   2.516; C[" C6 "][1]=   4.358; C[" C6 "][2]=  -3.667;
    ideal_rna["  C"]=C;
    map<string,vector<double> >().swap(C);

    map<string,vector<double> >G;
    G[" P  "]=tmp; G[" P  "][0]=  -5.706; G[" P  "][1]=  -7.249; G[" P  "][2]=  -5.223;
    G[" OP1"]=tmp; G[" OP1"][0]=  -6.417; G[" OP1"][1]=  -8.303; G[" OP1"][2]=  -5.980;
    G[" OP2"]=tmp; G[" OP2"][0]=  -6.461; G[" OP2"][1]=  -6.623; G[" OP2"][2]=  -4.115;
    G[" O5'"]=tmp; G[" O5'"][0]=  -4.320; G[" O5'"][1]=  -7.825; G[" O5'"][2]=  -4.670;
    G[" C5'"]=tmp; G[" C5'"][0]=  -3.159; G[" C5'"][1]=  -7.787; G[" C5'"][2]=  -5.521;
    G[" C4'"]=tmp; G[" C4'"][0]=  -1.951; G[" C4'"][1]=  -7.323; G[" C4'"][2]=  -4.730;
    G[" O4'"]=tmp; G[" O4'"][0]=  -1.813; G[" O4'"][1]=  -5.874; G[" O4'"][2]=  -4.727;
    G[" C3'"]=tmp; G[" C3'"][0]=  -1.943; G[" C3'"][1]=  -7.697; G[" C3'"][2]=  -3.248;
    G[" O3'"]=tmp; G[" O3'"][0]=  -0.610; G[" O3'"][1]=  -8.002; G[" O3'"][2]=  -2.857;
    G[" C2'"]=tmp; G[" C2'"][0]=  -2.383; G[" C2'"][1]=  -6.424; G[" C2'"][2]=  -2.525;
    G[" C1'"]=tmp; G[" C1'"][0]=  -1.593; G[" C1'"][1]=  -5.448; G[" C1'"][2]=  -3.391;
    G[" N9 "]=tmp; G[" N9 "][0]=  -2.053; G[" N9 "][1]=  -4.036; G[" N9 "][2]=  -3.274;
    G[" C8 "]=tmp; G[" C8 "][0]=  -3.327; G[" C8 "][1]=  -3.548; G[" C8 "][2]=  -3.081;
    G[" N7 "]=tmp; G[" N7 "][0]=  -3.395; G[" N7 "][1]=  -2.240; G[" N7 "][2]=  -3.021;
    G[" C5 "]=tmp; G[" C5 "][0]=  -2.072; G[" C5 "][1]=  -1.833; G[" C5 "][2]=  -3.185;
    G[" C6 "]=tmp; G[" C6 "][0]=  -1.514; G[" C6 "][1]=  -0.528; G[" C6 "][2]=  -3.210;
    G[" O6 "]=tmp; G[" O6 "][0]=  -2.081; G[" O6 "][1]=   0.554; G[" O6 "][2]=  -3.091;
    G[" N1 "]=tmp; G[" N1 "][0]=  -0.125; G[" N1 "][1]=  -0.568; G[" N1 "][2]=  -3.400;
    G[" C2 "]=tmp; G[" C2 "][0]=   0.625; G[" C2 "][1]=  -1.719; G[" C2 "][2]=  -3.547;
    G[" N2 "]=tmp; G[" N2 "][0]=   1.938; G[" N2 "][1]=  -1.545; G[" N2 "][2]=  -3.719;
    G[" N3 "]=tmp; G[" N3 "][0]=   0.101; G[" N3 "][1]=  -2.942; G[" N3 "][2]=  -3.524;
    G[" C4 "]=tmp; G[" C4 "][0]=  -1.245; G[" C4 "][1]=  -2.921; G[" C4 "][2]=  -3.340;
    ideal_rna["  G"]=G;
    map<string,vector<double> >().swap(G);

    map<string,vector<double> >U;
    U[" P  "]=tmp; U[" P  "][0]=  -0.356; U[" P  "][1]=  -9.218; U[" P  "][2]=  -1.848;
    U[" OP1"]=tmp; U[" OP1"][0]=  -0.311; U[" OP1"][1]= -10.489; U[" OP1"][2]=  -2.605;
    U[" OP2"]=tmp; U[" OP2"][0]=  -1.334; U[" OP2"][1]=  -9.156; U[" OP2"][2]=  -0.740;
    U[" O5'"]=tmp; U[" O5'"][0]=   1.105; U[" O5'"][1]=  -8.869; U[" O5'"][2]=  -1.295;
    U[" C5'"]=tmp; U[" C5'"][0]=   2.021; U[" C5'"][1]=  -8.156; U[" C5'"][2]=  -2.146;
    U[" C4'"]=tmp; U[" C4'"][0]=   2.726; U[" C4'"][1]=  -7.072; U[" C4'"][2]=  -1.355;
    U[" O4'"]=tmp; U[" O4'"][0]=   1.986; U[" O4'"][1]=  -5.817; U[" O4'"][2]=  -1.352;
    U[" C3'"]=tmp; U[" C3'"][0]=   2.952; U[" C3'"][1]=  -7.370; U[" C3'"][2]=   0.127;
    U[" O3'"]=tmp; U[" O3'"][0]=   4.210; U[" O3'"][1]=  -6.832; U[" O3'"][2]=   0.518;
    U[" C2'"]=tmp; U[" C2'"][0]=   1.848; U[" C2'"][1]=  -6.598; U[" C2'"][2]=   0.850;
    U[" C1'"]=tmp; U[" C1'"][0]=   1.913; U[" C1'"][1]=  -5.344; U[" C1'"][2]=  -0.016;
    U[" N1 "]=tmp; U[" N1 "][0]=   0.711; U[" N1 "][1]=  -4.472; U[" N1 "][2]=   0.101;
    U[" C2 "]=tmp; U[" C2 "][0]=   0.911; U[" C2 "][1]=  -3.116; U[" C2 "][2]=   0.009;
    U[" O2 "]=tmp; U[" O2 "][0]=   2.010; U[" O2 "][1]=  -2.617; U[" O2 "][2]=  -0.161;
    U[" N3 "]=tmp; U[" N3 "][0]=  -0.226; U[" N3 "][1]=  -2.339; U[" N3 "][2]=   0.123;
    U[" C4 "]=tmp; U[" C4 "][0]=  -1.513; U[" C4 "][1]=  -2.797; U[" C4 "][2]=   0.316;
    U[" O4 "]=tmp; U[" O4 "][0]=  -2.454; U[" O4 "][1]=  -2.004; U[" O4 "][2]=   0.403;
    U[" C5 "]=tmp; U[" C5 "][0]=  -1.622; U[" C5 "][1]=  -4.234; U[" C5 "][2]=   0.400;
    U[" C7 "]=tmp; U[" C7 "][0]=  -2.989; U[" C7 "][1]=  -4.816; U[" C7 "][2]=   0.610;
    U[" C6 "]=tmp; U[" C6 "][0]=  -0.533; U[" C6 "][1]=  -5.013; U[" C6 "][2]=   0.293;
    ideal_rna["  U"]=U;
    map<string,vector<double> >().swap(U);

    map<string,vector<double> >A0;
    A0[" P  "]=tmp; A0[" P  "][0]=  -0.356; A0[" P  "][1]=   9.218; A0[" P  "][2]=   1.848;
    A0[" C4'"]=tmp; A0[" C4'"][0]=   2.726; A0[" C4'"][1]=   7.072; A0[" C4'"][2]=   1.355;
    A0[" C1'"]=tmp; A0[" C1'"][0]=   1.913; A0[" C1'"][1]=   5.344; A0[" C1'"][2]=   0.016;
    map<string,vector<double> >A1;
    A1[" P  "]=tmp; A1[" P  "][0]=   5.130; A1[" P  "][1]=   7.667; A1[" P  "][2]=  -1.527;
    A1[" C4'"]=tmp; A1[" C4'"][0]=   6.362; A1[" C4'"][1]=   4.119; A1[" C4'"][2]=  -2.020;
    A1[" C1'"]=tmp; A1[" C1'"][0]=   4.689; A1[" C1'"][1]=   3.199; A1[" C1'"][2]=  -3.359;
    map<string,vector<double> >A2;
    A2[" P  "]=tmp; A2[" P  "][0]=   8.657; A2[" P  "][1]=   3.187; A2[" P  "][2]=  -4.902;
    A2[" C4'"]=tmp; A2[" C4'"][0]=   7.568; A2[" C4'"][1]=  -0.407; A2[" C4'"][2]=  -5.395;
    A2[" C1'"]=tmp; A2[" C1'"][0]=   5.674; A2[" C1'"][1]=  -0.168; A2[" C1'"][2]=  -6.734;
    map<string,vector<double> >B0;
    B0[" P  "]=tmp; B0[" P  "][0]=  -8.877; B0[" P  "][1]=  -2.510; B0[" P  "][2]=  -8.598;
    B0[" C4'"]=tmp; B0[" C4'"][0]=  -5.883; B0[" C4'"][1]=  -4.778; B0[" C4'"][2]=  -8.105;
    B0[" C1'"]=tmp; B0[" C1'"][0]=  -4.491; B0[" C1'"][1]=  -3.471; B0[" C1'"][2]=  -6.766;
    map<string,vector<double> >B1;
    B1[" P  "]=tmp; B1[" P  "][0]=  -5.706; B1[" P  "][1]=  -7.249; B1[" P  "][2]=  -5.223;
    B1[" C4'"]=tmp; B1[" C4'"][0]=  -1.951; B1[" C4'"][1]=  -7.323; B1[" C4'"][2]=  -4.730;
    B1[" C1'"]=tmp; B1[" C1'"][0]=  -1.593; B1[" C1'"][1]=  -5.448; B1[" C1'"][2]=  -3.391;
    map<string,vector<double> >B2;
    B2[" P  "]=tmp; B2[" P  "][0]=  -0.356; B2[" P  "][1]=  -9.218; B2[" P  "][2]=  -1.848;
    B2[" C4'"]=tmp; B2[" C4'"][0]=   2.726; B2[" C4'"][1]=  -7.072; B2[" C4'"][2]=  -1.355;
    B2[" C1'"]=tmp; B2[" C1'"][0]=   1.913; B2[" C1'"][1]=  -5.344; B2[" C1'"][2]=  -0.016;

    ideal_rna["  A  A0"]=A0; ideal_rna["  A  A1"]=A1; ideal_rna["  A  A2"]=B1; ideal_rna["  A  A3"]=B2;
    ideal_rna["  A  C0"]=A0; ideal_rna["  A  C1"]=A1; ideal_rna["  A  C2"]=B1; ideal_rna["  A  C3"]=B2;
    ideal_rna["  A  G0"]=A0; ideal_rna["  A  G1"]=A1; ideal_rna["  A  G2"]=B1; ideal_rna["  A  G3"]=B2;
    ideal_rna["  A  U0"]=A0; ideal_rna["  A  U1"]=A1; ideal_rna["  A  U2"]=B1; ideal_rna["  A  U3"]=B2;
    ideal_rna["  C  A0"]=A0; ideal_rna["  C  A1"]=A1; ideal_rna["  C  A2"]=B1; ideal_rna["  C  A3"]=B2;
    ideal_rna["  C  C0"]=A0; ideal_rna["  C  C1"]=A1; ideal_rna["  C  C2"]=B1; ideal_rna["  C  C3"]=B2;
    ideal_rna["  C  G0"]=A0; ideal_rna["  C  G1"]=A1; ideal_rna["  C  G2"]=B1; ideal_rna["  C  G3"]=B2;
    ideal_rna["  C  U0"]=A0; ideal_rna["  C  U1"]=A1; ideal_rna["  C  U2"]=B1; ideal_rna["  C  U3"]=B2;
    ideal_rna["  G  A0"]=A0; ideal_rna["  G  A1"]=A1; ideal_rna["  G  A2"]=B1; ideal_rna["  G  A3"]=B2;
    ideal_rna["  G  C0"]=A0; ideal_rna["  G  C1"]=A1; ideal_rna["  G  C2"]=B1; ideal_rna["  G  C3"]=B2;
    ideal_rna["  G  G0"]=A0; ideal_rna["  G  G1"]=A1; ideal_rna["  G  G2"]=B1; ideal_rna["  G  G3"]=B2;
    ideal_rna["  G  U0"]=A0; ideal_rna["  G  U1"]=A1; ideal_rna["  G  U2"]=B1; ideal_rna["  G  U3"]=B2;
    ideal_rna["  U  A0"]=A0; ideal_rna["  U  A1"]=A1; ideal_rna["  U  A2"]=B1; ideal_rna["  U  A3"]=B2;
    ideal_rna["  U  C0"]=A0; ideal_rna["  U  C1"]=A1; ideal_rna["  U  C2"]=B1; ideal_rna["  U  C3"]=B2;
    ideal_rna["  U  G0"]=A0; ideal_rna["  U  G1"]=A1; ideal_rna["  U  G2"]=B1; ideal_rna["  U  G3"]=B2;
    ideal_rna["  U  U0"]=A0; ideal_rna["  U  U1"]=A1; ideal_rna["  U  U2"]=B1; ideal_rna["  U  U3"]=B2;
    ideal_rna["  A  A  A0"]=A0; ideal_rna["  A  A  A1"]=A1; ideal_rna["  A  A  A2"]=A2;  ideal_rna["  A  A  A3"]=B0; ideal_rna["  A  A  A4"]=B1; ideal_rna["  A  A  A5"]=B2;
    ideal_rna["  A  A  C0"]=A0; ideal_rna["  A  A  C1"]=A1; ideal_rna["  A  A  C2"]=A2;  ideal_rna["  A  A  C3"]=B0; ideal_rna["  A  A  C4"]=B1; ideal_rna["  A  A  C5"]=B2;
    ideal_rna["  A  A  G0"]=A0; ideal_rna["  A  A  G1"]=A1; ideal_rna["  A  A  G2"]=A2;  ideal_rna["  A  A  G3"]=B0; ideal_rna["  A  A  G4"]=B1; ideal_rna["  A  A  G5"]=B2;
    ideal_rna["  A  A  U0"]=A0; ideal_rna["  A  A  U1"]=A1; ideal_rna["  A  A  U2"]=A2;  ideal_rna["  A  A  U3"]=B0; ideal_rna["  A  A  U4"]=B1; ideal_rna["  A  A  U5"]=B2;
    ideal_rna["  A  C  A0"]=A0; ideal_rna["  A  C  A1"]=A1; ideal_rna["  A  C  A2"]=A2;  ideal_rna["  A  C  A3"]=B0; ideal_rna["  A  C  A4"]=B1; ideal_rna["  A  C  A5"]=B2;
    ideal_rna["  A  C  C0"]=A0; ideal_rna["  A  C  C1"]=A1; ideal_rna["  A  C  C2"]=A2;  ideal_rna["  A  C  C3"]=B0; ideal_rna["  A  C  C4"]=B1; ideal_rna["  A  C  C5"]=B2;
    ideal_rna["  A  C  G0"]=A0; ideal_rna["  A  C  G1"]=A1; ideal_rna["  A  C  G2"]=A2;  ideal_rna["  A  C  G3"]=B0; ideal_rna["  A  C  G4"]=B1; ideal_rna["  A  C  G5"]=B2;
    ideal_rna["  A  C  U0"]=A0; ideal_rna["  A  C  U1"]=A1; ideal_rna["  A  C  U2"]=A2;  ideal_rna["  A  C  U3"]=B0; ideal_rna["  A  C  U4"]=B1; ideal_rna["  A  C  U5"]=B2;
    ideal_rna["  A  G  A0"]=A0; ideal_rna["  A  G  A1"]=A1; ideal_rna["  A  G  A2"]=A2;  ideal_rna["  A  G  A3"]=B0; ideal_rna["  A  G  A4"]=B1; ideal_rna["  A  G  A5"]=B2;
    ideal_rna["  A  G  C0"]=A0; ideal_rna["  A  G  C1"]=A1; ideal_rna["  A  G  C2"]=A2;  ideal_rna["  A  G  C3"]=B0; ideal_rna["  A  G  C4"]=B1; ideal_rna["  A  G  C5"]=B2;
    ideal_rna["  A  G  G0"]=A0; ideal_rna["  A  G  G1"]=A1; ideal_rna["  A  G  G2"]=A2;  ideal_rna["  A  G  G3"]=B0; ideal_rna["  A  G  G4"]=B1; ideal_rna["  A  G  G5"]=B2;
    ideal_rna["  A  G  U0"]=A0; ideal_rna["  A  G  U1"]=A1; ideal_rna["  A  G  U2"]=A2;  ideal_rna["  A  G  U3"]=B0; ideal_rna["  A  G  U4"]=B1; ideal_rna["  A  G  U5"]=B2;
    ideal_rna["  A  U  A0"]=A0; ideal_rna["  A  U  A1"]=A1; ideal_rna["  A  U  A2"]=A2;  ideal_rna["  A  U  A3"]=B0; ideal_rna["  A  U  A4"]=B1; ideal_rna["  A  U  A5"]=B2;
    ideal_rna["  A  U  C0"]=A0; ideal_rna["  A  U  C1"]=A1; ideal_rna["  A  U  C2"]=A2;  ideal_rna["  A  U  C3"]=B0; ideal_rna["  A  U  C4"]=B1; ideal_rna["  A  U  C5"]=B2;
    ideal_rna["  A  U  G0"]=A0; ideal_rna["  A  U  G1"]=A1; ideal_rna["  A  U  G2"]=A2;  ideal_rna["  A  U  G3"]=B0; ideal_rna["  A  U  G4"]=B1; ideal_rna["  A  U  G5"]=B2;
    ideal_rna["  A  U  U0"]=A0; ideal_rna["  A  U  U1"]=A1; ideal_rna["  A  U  U2"]=A2;  ideal_rna["  A  U  U3"]=B0; ideal_rna["  A  U  U4"]=B1; ideal_rna["  A  U  U5"]=B2;
    ideal_rna["  C  A  A0"]=A0; ideal_rna["  C  A  A1"]=A1; ideal_rna["  C  A  A2"]=A2;  ideal_rna["  C  A  A3"]=B0; ideal_rna["  C  A  A4"]=B1; ideal_rna["  C  A  A5"]=B2;
    ideal_rna["  C  A  C0"]=A0; ideal_rna["  C  A  C1"]=A1; ideal_rna["  C  A  C2"]=A2;  ideal_rna["  C  A  C3"]=B0; ideal_rna["  C  A  C4"]=B1; ideal_rna["  C  A  C5"]=B2;
    ideal_rna["  C  A  G0"]=A0; ideal_rna["  C  A  G1"]=A1; ideal_rna["  C  A  G2"]=A2;  ideal_rna["  C  A  G3"]=B0; ideal_rna["  C  A  G4"]=B1; ideal_rna["  C  A  G5"]=B2;
    ideal_rna["  C  A  U0"]=A0; ideal_rna["  C  A  U1"]=A1; ideal_rna["  C  A  U2"]=A2;  ideal_rna["  C  A  U3"]=B0; ideal_rna["  C  A  U4"]=B1; ideal_rna["  C  A  U5"]=B2;
    ideal_rna["  C  C  A0"]=A0; ideal_rna["  C  C  A1"]=A1; ideal_rna["  C  C  A2"]=A2;  ideal_rna["  C  C  A3"]=B0; ideal_rna["  C  C  A4"]=B1; ideal_rna["  C  C  A5"]=B2;
    ideal_rna["  C  C  C0"]=A0; ideal_rna["  C  C  C1"]=A1; ideal_rna["  C  C  C2"]=A2;  ideal_rna["  C  C  C3"]=B0; ideal_rna["  C  C  C4"]=B1; ideal_rna["  C  C  C5"]=B2;
    ideal_rna["  C  C  G0"]=A0; ideal_rna["  C  C  G1"]=A1; ideal_rna["  C  C  G2"]=A2;  ideal_rna["  C  C  G3"]=B0; ideal_rna["  C  C  G4"]=B1; ideal_rna["  C  C  G5"]=B2;
    ideal_rna["  C  C  U0"]=A0; ideal_rna["  C  C  U1"]=A1; ideal_rna["  C  C  U2"]=A2;  ideal_rna["  C  C  U3"]=B0; ideal_rna["  C  C  U4"]=B1; ideal_rna["  C  C  U5"]=B2;
    ideal_rna["  C  G  A0"]=A0; ideal_rna["  C  G  A1"]=A1; ideal_rna["  C  G  A2"]=A2;  ideal_rna["  C  G  A3"]=B0; ideal_rna["  C  G  A4"]=B1; ideal_rna["  C  G  A5"]=B2;
    ideal_rna["  C  G  C0"]=A0; ideal_rna["  C  G  C1"]=A1; ideal_rna["  C  G  C2"]=A2;  ideal_rna["  C  G  C3"]=B0; ideal_rna["  C  G  C4"]=B1; ideal_rna["  C  G  C5"]=B2;
    ideal_rna["  C  G  G0"]=A0; ideal_rna["  C  G  G1"]=A1; ideal_rna["  C  G  G2"]=A2;  ideal_rna["  C  G  G3"]=B0; ideal_rna["  C  G  G4"]=B1; ideal_rna["  C  G  G5"]=B2;
    ideal_rna["  C  G  U0"]=A0; ideal_rna["  C  G  U1"]=A1; ideal_rna["  C  G  U2"]=A2;  ideal_rna["  C  G  U3"]=B0; ideal_rna["  C  G  U4"]=B1; ideal_rna["  C  G  U5"]=B2;
    ideal_rna["  C  U  A0"]=A0; ideal_rna["  C  U  A1"]=A1; ideal_rna["  C  U  A2"]=A2;  ideal_rna["  C  U  A3"]=B0; ideal_rna["  C  U  A4"]=B1; ideal_rna["  C  U  A5"]=B2;
    ideal_rna["  C  U  C0"]=A0; ideal_rna["  C  U  C1"]=A1; ideal_rna["  C  U  C2"]=A2;  ideal_rna["  C  U  C3"]=B0; ideal_rna["  C  U  C4"]=B1; ideal_rna["  C  U  C5"]=B2;
    ideal_rna["  C  U  G0"]=A0; ideal_rna["  C  U  G1"]=A1; ideal_rna["  C  U  G2"]=A2;  ideal_rna["  C  U  G3"]=B0; ideal_rna["  C  U  G4"]=B1; ideal_rna["  C  U  G5"]=B2;
    ideal_rna["  C  U  U0"]=A0; ideal_rna["  C  U  U1"]=A1; ideal_rna["  C  U  U2"]=A2;  ideal_rna["  C  U  U3"]=B0; ideal_rna["  C  U  U4"]=B1; ideal_rna["  C  U  U5"]=B2;
    ideal_rna["  G  A  A0"]=A0; ideal_rna["  G  A  A1"]=A1; ideal_rna["  G  A  A2"]=A2;  ideal_rna["  G  A  A3"]=B0; ideal_rna["  G  A  A4"]=B1; ideal_rna["  G  A  A5"]=B2;
    ideal_rna["  G  A  C0"]=A0; ideal_rna["  G  A  C1"]=A1; ideal_rna["  G  A  C2"]=A2;  ideal_rna["  G  A  C3"]=B0; ideal_rna["  G  A  C4"]=B1; ideal_rna["  G  A  C5"]=B2;
    ideal_rna["  G  A  G0"]=A0; ideal_rna["  G  A  G1"]=A1; ideal_rna["  G  A  G2"]=A2;  ideal_rna["  G  A  G3"]=B0; ideal_rna["  G  A  G4"]=B1; ideal_rna["  G  A  G5"]=B2;
    ideal_rna["  G  A  U0"]=A0; ideal_rna["  G  A  U1"]=A1; ideal_rna["  G  A  U2"]=A2;  ideal_rna["  G  A  U3"]=B0; ideal_rna["  G  A  U4"]=B1; ideal_rna["  G  A  U5"]=B2;
    ideal_rna["  G  C  A0"]=A0; ideal_rna["  G  C  A1"]=A1; ideal_rna["  G  C  A2"]=A2;  ideal_rna["  G  C  A3"]=B0; ideal_rna["  G  C  A4"]=B1; ideal_rna["  G  C  A5"]=B2;
    ideal_rna["  G  C  C0"]=A0; ideal_rna["  G  C  C1"]=A1; ideal_rna["  G  C  C2"]=A2;  ideal_rna["  G  C  C3"]=B0; ideal_rna["  G  C  C4"]=B1; ideal_rna["  G  C  C5"]=B2;
    ideal_rna["  G  C  G0"]=A0; ideal_rna["  G  C  G1"]=A1; ideal_rna["  G  C  G2"]=A2;  ideal_rna["  G  C  G3"]=B0; ideal_rna["  G  C  G4"]=B1; ideal_rna["  G  C  G5"]=B2;
    ideal_rna["  G  C  U0"]=A0; ideal_rna["  G  C  U1"]=A1; ideal_rna["  G  C  U2"]=A2;  ideal_rna["  G  C  U3"]=B0; ideal_rna["  G  C  U4"]=B1; ideal_rna["  G  C  U5"]=B2;
    ideal_rna["  G  G  A0"]=A0; ideal_rna["  G  G  A1"]=A1; ideal_rna["  G  G  A2"]=A2;  ideal_rna["  G  G  A3"]=B0; ideal_rna["  G  G  A4"]=B1; ideal_rna["  G  G  A5"]=B2;
    ideal_rna["  G  G  C0"]=A0; ideal_rna["  G  G  C1"]=A1; ideal_rna["  G  G  C2"]=A2;  ideal_rna["  G  G  C3"]=B0; ideal_rna["  G  G  C4"]=B1; ideal_rna["  G  G  C5"]=B2;
    ideal_rna["  G  G  G0"]=A0; ideal_rna["  G  G  G1"]=A1; ideal_rna["  G  G  G2"]=A2;  ideal_rna["  G  G  G3"]=B0; ideal_rna["  G  G  G4"]=B1; ideal_rna["  G  G  G5"]=B2;
    ideal_rna["  G  G  U0"]=A0; ideal_rna["  G  G  U1"]=A1; ideal_rna["  G  G  U2"]=A2;  ideal_rna["  G  G  U3"]=B0; ideal_rna["  G  G  U4"]=B1; ideal_rna["  G  G  U5"]=B2;
    ideal_rna["  G  U  A0"]=A0; ideal_rna["  G  U  A1"]=A1; ideal_rna["  G  U  A2"]=A2;  ideal_rna["  G  U  A3"]=B0; ideal_rna["  G  U  A4"]=B1; ideal_rna["  G  U  A5"]=B2;
    ideal_rna["  G  U  C0"]=A0; ideal_rna["  G  U  C1"]=A1; ideal_rna["  G  U  C2"]=A2;  ideal_rna["  G  U  C3"]=B0; ideal_rna["  G  U  C4"]=B1; ideal_rna["  G  U  C5"]=B2;
    ideal_rna["  G  U  G0"]=A0; ideal_rna["  G  U  G1"]=A1; ideal_rna["  G  U  G2"]=A2;  ideal_rna["  G  U  G3"]=B0; ideal_rna["  G  U  G4"]=B1; ideal_rna["  G  U  G5"]=B2;
    ideal_rna["  G  U  U0"]=A0; ideal_rna["  G  U  U1"]=A1; ideal_rna["  G  U  U2"]=A2;  ideal_rna["  G  U  U3"]=B0; ideal_rna["  G  U  U4"]=B1; ideal_rna["  G  U  U5"]=B2;
    ideal_rna["  U  A  A0"]=A0; ideal_rna["  U  A  A1"]=A1; ideal_rna["  U  A  A2"]=A2;  ideal_rna["  U  A  A3"]=B0; ideal_rna["  U  A  A4"]=B1; ideal_rna["  U  A  A5"]=B2;
    ideal_rna["  U  A  C0"]=A0; ideal_rna["  U  A  C1"]=A1; ideal_rna["  U  A  C2"]=A2;  ideal_rna["  U  A  C3"]=B0; ideal_rna["  U  A  C4"]=B1; ideal_rna["  U  A  C5"]=B2;
    ideal_rna["  U  A  G0"]=A0; ideal_rna["  U  A  G1"]=A1; ideal_rna["  U  A  G2"]=A2;  ideal_rna["  U  A  G3"]=B0; ideal_rna["  U  A  G4"]=B1; ideal_rna["  U  A  G5"]=B2;
    ideal_rna["  U  A  U0"]=A0; ideal_rna["  U  A  U1"]=A1; ideal_rna["  U  A  U2"]=A2;  ideal_rna["  U  A  U3"]=B0; ideal_rna["  U  A  U4"]=B1; ideal_rna["  U  A  U5"]=B2;
    ideal_rna["  U  C  A0"]=A0; ideal_rna["  U  C  A1"]=A1; ideal_rna["  U  C  A2"]=A2;  ideal_rna["  U  C  A3"]=B0; ideal_rna["  U  C  A4"]=B1; ideal_rna["  U  C  A5"]=B2;
    ideal_rna["  U  C  C0"]=A0; ideal_rna["  U  C  C1"]=A1; ideal_rna["  U  C  C2"]=A2;  ideal_rna["  U  C  C3"]=B0; ideal_rna["  U  C  C4"]=B1; ideal_rna["  U  C  C5"]=B2;
    ideal_rna["  U  C  G0"]=A0; ideal_rna["  U  C  G1"]=A1; ideal_rna["  U  C  G2"]=A2;  ideal_rna["  U  C  G3"]=B0; ideal_rna["  U  C  G4"]=B1; ideal_rna["  U  C  G5"]=B2;
    ideal_rna["  U  C  U0"]=A0; ideal_rna["  U  C  U1"]=A1; ideal_rna["  U  C  U2"]=A2;  ideal_rna["  U  C  U3"]=B0; ideal_rna["  U  C  U4"]=B1; ideal_rna["  U  C  U5"]=B2;
    ideal_rna["  U  G  A0"]=A0; ideal_rna["  U  G  A1"]=A1; ideal_rna["  U  G  A2"]=A2;  ideal_rna["  U  G  A3"]=B0; ideal_rna["  U  G  A4"]=B1; ideal_rna["  U  G  A5"]=B2;
    ideal_rna["  U  G  C0"]=A0; ideal_rna["  U  G  C1"]=A1; ideal_rna["  U  G  C2"]=A2;  ideal_rna["  U  G  C3"]=B0; ideal_rna["  U  G  C4"]=B1; ideal_rna["  U  G  C5"]=B2;
    ideal_rna["  U  G  G0"]=A0; ideal_rna["  U  G  G1"]=A1; ideal_rna["  U  G  G2"]=A2;  ideal_rna["  U  G  G3"]=B0; ideal_rna["  U  G  G4"]=B1; ideal_rna["  U  G  G5"]=B2;
    ideal_rna["  U  G  U0"]=A0; ideal_rna["  U  G  U1"]=A1; ideal_rna["  U  G  U2"]=A2;  ideal_rna["  U  G  U3"]=B0; ideal_rna["  U  G  U4"]=B1; ideal_rna["  U  G  U5"]=B2;
    ideal_rna["  U  U  A0"]=A0; ideal_rna["  U  U  A1"]=A1; ideal_rna["  U  U  A2"]=A2;  ideal_rna["  U  U  A3"]=B0; ideal_rna["  U  U  A4"]=B1; ideal_rna["  U  U  A5"]=B2;
    ideal_rna["  U  U  C0"]=A0; ideal_rna["  U  U  C1"]=A1; ideal_rna["  U  U  C2"]=A2;  ideal_rna["  U  U  C3"]=B0; ideal_rna["  U  U  C4"]=B1; ideal_rna["  U  U  C5"]=B2;
    ideal_rna["  U  U  G0"]=A0; ideal_rna["  U  U  G1"]=A1; ideal_rna["  U  U  G2"]=A2;  ideal_rna["  U  U  G3"]=B0; ideal_rna["  U  U  G4"]=B1; ideal_rna["  U  U  G5"]=B2;
    ideal_rna["  U  U  U0"]=A0; ideal_rna["  U  U  U1"]=A1; ideal_rna["  U  U  U2"]=A2;  ideal_rna["  U  U  U3"]=B0; ideal_rna["  U  U  U4"]=B1; ideal_rna["  U  U  U5"]=B2;

    map<string,vector<double> >().swap(A0);
    map<string,vector<double> >().swap(B0);
    map<string,vector<double> >().swap(A1);
    map<string,vector<double> >().swap(B1);
    map<string,vector<double> >().swap(A2);
    map<string,vector<double> >().swap(B2);

    vector<vector<double> > xyz_list1(3,tmp);
    vector<vector<double> > xyz_list2(3,tmp);
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t
    for (map<string, map<string,vector<double> > >::iterator iter = ideal_rna.begin();
        iter != ideal_rna.end(); iter++)
    {
        string key =  iter->first;
        if (key.size()<=3) continue;
        int idx=(char)(key[key.size()-1])-'0';
        bool reverse=false;
        if (key.size()==7 && idx>=2)
        {
            reverse=true;
            idx=(3-idx);
        }
        else if (key.size()==10 && idx>=3)
        {
            reverse=true;
            idx=(5-idx);
        }
        string nt=key.substr(idx*3,3);
        if (reverse)
        {
            if      (nt=="  A") nt="  U";
            else if (nt=="  C") nt="  G";
            else if (nt=="  G") nt="  C";
            else if (nt=="  U") nt="  A";
        }
    
        xyz_list1[0][0]=ideal_rna[nt][" P  "][0];
        xyz_list1[0][1]=ideal_rna[nt][" P  "][1];
        xyz_list1[0][2]=ideal_rna[nt][" P  "][2];
        xyz_list1[1][0]=ideal_rna[nt][" C4'"][0];
        xyz_list1[1][1]=ideal_rna[nt][" C4'"][1];
        xyz_list1[1][2]=ideal_rna[nt][" C4'"][2];
        xyz_list1[2][0]=ideal_rna[nt][" C1'"][0];
        xyz_list1[2][1]=ideal_rna[nt][" C1'"][1];
        xyz_list1[2][2]=ideal_rna[nt][" C1'"][2];
        
        xyz_list2[0][0]=ideal_rna[key][" P  "][0];
        xyz_list2[0][1]=ideal_rna[key][" P  "][1];
        xyz_list2[0][2]=ideal_rna[key][" P  "][2];
        xyz_list2[1][0]=ideal_rna[key][" C4'"][0];
        xyz_list2[1][1]=ideal_rna[key][" C4'"][1];
        xyz_list2[1][2]=ideal_rna[key][" C4'"][2];
        xyz_list2[2][0]=ideal_rna[key][" C1'"][0];
        xyz_list2[2][1]=ideal_rna[key][" C1'"][1];
        xyz_list2[2][2]=ideal_rna[key][" C1'"][2];
        
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
        
        for (map<string, vector<double> >::iterator i = ideal_rna[nt].begin();
            i != ideal_rna[nt].end(); i++)
        {
            string k=i-> first;
            if (ideal_rna[key].count(k)) continue;
            ideal_rna[key][k]=tmp;
            ChangeCoor(ideal_rna[nt][k], RotMatix, TranVect, ideal_rna[key][k]);
        }
    }

    vector<vector<double> >().swap(xyz_list1);
    vector<vector<double> >().swap(xyz_list2);
    vector<vector<double> >().swap(RotMatix);  // U
    vector<double> ().swap(TranVect);  // t
    vector<double> ().swap(tmp);
    return ideal_rna;
}
