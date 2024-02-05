#include <math.h>
#include <stdio.h>
extern double fabs(double x);
extern double acos(double x);
extern double acosh(double x);
extern double atan(double x);
extern double atanh(double x);
extern double asin(double x);
extern double asinh(double x);
extern double acos(double x);
extern double acosh(double x);
extern double asin(double x);
extern double asinh(double x);
extern double atan(double x);
extern double atanh(double x);
extern double ceil(double x);
extern double cos(double x);
extern double cosh(double x);
extern double tan(double x);
extern double tanh(double x);
extern double sin(double x);
extern double sinh(double x);
extern double exp(double x);
extern double floor(double x);
extern double pow(double x, double y);
extern double factorial(double x);
extern double log(double x);
extern double arbitrary_log(double x, double base);
extern double gcd_pair(double a, double b);
extern double lcm_pair(double a, double b);
extern double gcd_multi(unsigned int size, ...);
extern double lcm_multi(unsigned int size, ...);
extern double multi_min(unsigned int size, ...);
extern double multi_max(unsigned int size, ...);
extern void NR_MINIMISE(double (*func)(double VOI, double *C, double *R,
                                       double *S, double *A),
                        double VOI, double *C, double *R, double *S, double *A,
                        double *V);

void OC_CellML_RHS_routine(double VOI, double *OC_STATE, double *OC_RATE,
                           double *OC_WANTED, double *OC_KNOWN) {

  double DUMMY_ASSIGNMENT;
  double CONSTANTS[110], ALGEBRAIC[70];

  /* Constant C_m */
  CONSTANTS[0] = 0.58;
  /* Constant gam */
  CONSTANTS[1] = 2.79;
  /* Constant R_a */
  CONSTANTS[2] = 150;
  /* Constant tsi */
  CONSTANTS[3] = 0.000001;
  /* Constant tsi2 */
  CONSTANTS[4] = 0.0025;
  /* Constant tsi3 */
  CONSTANTS[5] = 0.0005;
  /* Constant FF */
  CONSTANTS[6] = 96485;
  /* Constant tau_K */
  CONSTANTS[7] = 559;
  /* Constant tau_Na */
  CONSTANTS[8] = 559;
  /* Constant f_T */
  CONSTANTS[9] = 0.00174;
  /* Constant tau_K2 */
  CONSTANTS[10] = 40229.885;
  /* Constant tau_Na2 */
  CONSTANTS[11] = 40229.885;
  /* Constant I_K_rest */
  CONSTANTS[12] = 0.34;
  /* Constant I_Na_rest */
  CONSTANTS[13] = -0.43;
  /* Constant alpha_h_bar */
  CONSTANTS[14] = 0.0081;
  /* Constant alpha_m_bar */
  CONSTANTS[15] = 0.288;
  /* Constant alpha_n_bar */
  CONSTANTS[16] = 0.0131;
  /* Constant beta_h_bar */
  CONSTANTS[17] = 4.38;
  /* Constant beta_m_bar */
  CONSTANTS[18] = 1.38;
  /* Constant beta_n_bar */
  CONSTANTS[19] = 0.067;
  /* Constant V_m */
  CONSTANTS[20] = -46;
  /* Constant V_n */
  CONSTANTS[21] = -40;
  /* Constant V_h */
  CONSTANTS[22] = -45;
  /* Constant V_a */
  CONSTANTS[23] = 70;
  /* Constant V_S_inf */
  CONSTANTS[24] = -68;
  /* Constant V_h_K_inf */
  CONSTANTS[25] = -40;
  /* Constant A_a */
  CONSTANTS[26] = 150;
  /* Constant A_S_inf */
  CONSTANTS[27] = 7.1;
  /* Constant A_h_K_inf */
  CONSTANTS[28] = 7.5;
  /* Constant K_alpha_h */
  CONSTANTS[29] = 14.7;
  /* Constant K_beta_h */
  CONSTANTS[30] = 9;
  /* Constant K_alpha_m */
  CONSTANTS[31] = 10;
  /* Constant K_alpha_n */
  CONSTANTS[32] = 7;
  /* Constant K_beta_m */
  CONSTANTS[33] = 18;
  /* Constant K_beta_n */
  CONSTANTS[34] = 40;
  /* Constant RR */
  CONSTANTS[35] = 8314.41;
  /* Constant TT */
  CONSTANTS[36] = 293;
  /* Constant g_Cl_bar */
  CONSTANTS[37] = 3.275;
  /* Constant g_K_bar */
  CONSTANTS[38] = 10.8;
  /* Constant g_Na_bar */
  CONSTANTS[39] = 134;
  /* Constant G_K */
  CONSTANTS[40] = 1.85;
  /* Constant del */
  CONSTANTS[41] = 0.4;
  /* Constant K_K */
  CONSTANTS[42] = 950;
  /* Constant K_S */
  CONSTANTS[43] = 1;
  /* Constant K_m_K */
  CONSTANTS[44] = 1;
  /* Constant K_m_Na */
  CONSTANTS[45] = 13;
  /* Constant S_i */
  CONSTANTS[46] = 10;
  /* Constant J_NaK_bar */
  CONSTANTS[47] = 0.0001656;
  /* Constant V_tau */
  CONSTANTS[48] = 70;
  /* Constant vS */
  DUMMY_ASSIGNMENT /*OC_STATE[0]*/ = -79.974;
  /* Constant vT */
  DUMMY_ASSIGNMENT /*OC_STATE[1]*/ = -80.2;
  /* Constant K_t */
  DUMMY_ASSIGNMENT /*OC_STATE[2]*/ = 5.9;
  /* Constant K_i */
  DUMMY_ASSIGNMENT /*OC_STATE[3]*/ = 150.9;
  /* Constant K_e */
  DUMMY_ASSIGNMENT /*OC_STATE[4]*/ = 5.9;
  /* Constant Na_i */
  DUMMY_ASSIGNMENT /*OC_STATE[5]*/ = 12.7;
  /* Constant Na_t */
  DUMMY_ASSIGNMENT /*OC_STATE[6]*/ = 132.0;
  /* Constant Na_e */
  DUMMY_ASSIGNMENT /*OC_STATE[7]*/ = 133.0;
  /* Constant eta_Cl */
  CONSTANTS[49] = 0.1;
  /* Constant eta_IR */
  CONSTANTS[50] = 1.0;
  /* Constant eta_DR */
  CONSTANTS[51] = 0.45;
  /* Constant eta_Na */
  CONSTANTS[52] = 0.1;
  /* Constant eta_NaK */
  CONSTANTS[53] = 0.1;
  /* Constant I_HH */
  DUMMY_ASSIGNMENT /*OC_KNOWN[0]*/ = 0.0;
  /* Constant n */
  DUMMY_ASSIGNMENT /*OC_STATE[8]*/ = 0.009466;
  /* Constant h_K */
  DUMMY_ASSIGNMENT /*OC_STATE[9]*/ = 0.9952;
  /* Constant m */
  DUMMY_ASSIGNMENT /*OC_STATE[10]*/ = 0.0358;
  /* Constant h */
  DUMMY_ASSIGNMENT /*OC_STATE[11]*/ = 0.4981;
  /* Constant S */
  DUMMY_ASSIGNMENT /*OC_STATE[12]*/ = 0.581;
  /* Constant n_t */
  DUMMY_ASSIGNMENT /*OC_STATE[13]*/ = 0.009466;
  /* Constant h_K_t */
  DUMMY_ASSIGNMENT /*OC_STATE[14]*/ = 0.9952;
  /* Constant m_t */
  DUMMY_ASSIGNMENT /*OC_STATE[15]*/ = 0.0358;
  /* Constant h_t */
  DUMMY_ASSIGNMENT /*OC_STATE[16]*/ = 0.4981;
  /* Constant S_t */
  DUMMY_ASSIGNMENT /*OC_STATE[17]*/ = 0.581;
  /* Constant O_0 */
  DUMMY_ASSIGNMENT /*OC_STATE[18]*/ = 0.0;
  /* Constant O_1 */
  DUMMY_ASSIGNMENT /*OC_STATE[19]*/ = 0.0;
  /* Constant O_2 */
  DUMMY_ASSIGNMENT /*OC_STATE[20]*/ = 0.0;
  /* Constant O_3 */
  DUMMY_ASSIGNMENT /*OC_STATE[21]*/ = 0.0;
  /* Constant O_4 */
  DUMMY_ASSIGNMENT /*OC_STATE[22]*/ = 0.0;
  /* Constant k_L */
  CONSTANTS[54] = 0.002;
  /* Constant k_Lm */
  CONSTANTS[55] = 1000;
  /* Constant f */
  CONSTANTS[56] = 0.2;
  /* Constant alpha1 */
  CONSTANTS[57] = 0.2;
  /* Constant K */
  CONSTANTS[58] = 4.5;
  /* Constant Vbar */
  CONSTANTS[59] = -20;
  /* Constant C_0 */
  DUMMY_ASSIGNMENT /*OC_STATE[23]*/ = 1.0;
  /* Constant C_1 */
  DUMMY_ASSIGNMENT /*OC_STATE[24]*/ = 0.0;
  /* Constant C_2 */
  DUMMY_ASSIGNMENT /*OC_STATE[25]*/ = 0.0;
  /* Constant C_3 */
  DUMMY_ASSIGNMENT /*OC_STATE[26]*/ = 0.0;
  /* Constant C_4 */
  DUMMY_ASSIGNMENT /*OC_STATE[27]*/ = 0.0;
  /* Constant nu_SR */
  CONSTANTS[60] = 2.4375;
  /* Constant K_SR */
  CONSTANTS[61] = 1;
  /* Constant L_e */
  CONSTANTS[62] = 0.00004;
  /* Constant tau_R */
  CONSTANTS[63] = 0.75;
  /* Constant tau_SR_R */
  CONSTANTS[64] = 0.75;
  /* Constant L_S_0 */
  CONSTANTS[65] = 1.0;
  /* Constant L_S */
  DUMMY_ASSIGNMENT /*OC_KNOWN[1]*/ = 1.0;
  /* Constant R_R */
  CONSTANTS[66] = 0.5;
  /* Constant k_T_on */
  CONSTANTS[67] = 0.0885;
  /* Constant k_T_off */
  CONSTANTS[68] = 0.115;
  /* Constant T_tot_0 */
  CONSTANTS[69] = 140;
  /* Constant k_P_on */
  CONSTANTS[70] = 0;
  /* Constant k_P_off */
  CONSTANTS[71] = 0;
  /* Constant P_tot */
  CONSTANTS[72] = 1500;
  /* Constant k_Mg_on */
  CONSTANTS[73] = 0;
  /* Constant k_Mg_off */
  CONSTANTS[74] = 0;
  /* Constant k_Cs_on */
  CONSTANTS[75] = 0.000004;
  /* Constant k_Cs_off */
  CONSTANTS[76] = 0.005;
  /* Constant Cs_tot */
  CONSTANTS[77] = 31000;
  /* Constant k_CATP_on */
  CONSTANTS[78] = 0.15;
  /* Constant k_CATP_off */
  CONSTANTS[79] = 30;
  /* Constant k_MATP_on */
  CONSTANTS[80] = 0.0015;
  /* Constant k_MATP_off */
  CONSTANTS[81] = 0.15;
  /* Constant tau_ATP */
  CONSTANTS[82] = 0.375;
  /* Constant tau_Mg */
  CONSTANTS[83] = 1.5;
  /* Constant k_0_on */
  CONSTANTS[84] = 0;
  /* Constant k_0_off */
  CONSTANTS[85] = 0.15;
  /* Constant k_Ca_on */
  CONSTANTS[86] = 0.15;
  /* Constant k_Ca_off */
  CONSTANTS[87] = 0.05;
  /* Constant f_o */
  CONSTANTS[88] = 0.5;
  /* Constant f_p */
  CONSTANTS[89] = 5;
  /* Constant h_o */
  CONSTANTS[90] = 0.08;
  /* Constant h_p */
  CONSTANTS[91] = 0.06;
  /* Constant g_o */
  CONSTANTS[92] = 0.04;
  /* Constant b_p */
  CONSTANTS[93] = 0.00000394;
  /* Constant k_p */
  CONSTANTS[94] = 0.00000362;
  /* Constant A_p */
  CONSTANTS[95] = 1;
  /* Constant B_p */
  CONSTANTS[96] = 0.0001;
  /* Constant PP */
  CONSTANTS[97] = 6;
  /* Constant x_0 */
  CONSTANTS[98] = 0.05;
  /* Constant x_1 */
  CONSTANTS[99] = 0.0;
  /* Constant x_2 */
  CONSTANTS[100] = 0.05;
  /* Constant eta */
  CONSTANTS[101] = 0.000107;
  /* Constant dummy */
  DUMMY_ASSIGNMENT /*OC_STATE[28]*/ = 0.0;
  /* Constant zeta */
  CONSTANTS[102] = 0.0021;
  /* Constant Ca_1 */
  DUMMY_ASSIGNMENT /*OC_STATE[29]*/ = 0.1;
  /* Constant Ca_SR1 */
  DUMMY_ASSIGNMENT /*OC_STATE[30]*/ = 1500.0;
  /* Constant Ca_2 */
  DUMMY_ASSIGNMENT /*OC_STATE[31]*/ = 0.1;
  /* Constant Ca_SR2 */
  DUMMY_ASSIGNMENT /*OC_STATE[32]*/ = 1500.0;
  /* Constant Ca_T_2 */
  DUMMY_ASSIGNMENT /*OC_STATE[33]*/ = 25;
  /* Constant Ca_P1 */
  DUMMY_ASSIGNMENT /*OC_STATE[34]*/ = 615.000000;
  /* Constant Ca_P2 */
  DUMMY_ASSIGNMENT /*OC_STATE[35]*/ = 615.000000;
  /* Constant Mg_P1 */
  DUMMY_ASSIGNMENT /*OC_STATE[36]*/ = 811.000000;
  /* Constant Mg_P2 */
  DUMMY_ASSIGNMENT /*OC_STATE[37]*/ = 811.000000;
  /* Constant Ca_Cs1 */
  DUMMY_ASSIGNMENT /*OC_STATE[38]*/ = 16900.0;
  /* Constant Ca_Cs2 */
  DUMMY_ASSIGNMENT /*OC_STATE[39]*/ = 16900.0;
  /* Constant Ca_ATP1 */
  DUMMY_ASSIGNMENT /*OC_STATE[40]*/ = 0.4;
  /* Constant Ca_ATP2 */
  DUMMY_ASSIGNMENT /*OC_STATE[41]*/ = 0.4;
  /* Constant Mg_ATP1 */
  DUMMY_ASSIGNMENT /*OC_STATE[42]*/ = 7200.0;
  /* Constant Mg_ATP2 */
  DUMMY_ASSIGNMENT /*OC_STATE[43]*/ = 7200.0;
  /* Constant ATP1 */
  DUMMY_ASSIGNMENT /*OC_STATE[44]*/ = 799.6;
  /* Constant ATP2 */
  DUMMY_ASSIGNMENT /*OC_STATE[45]*/ = 799.6;
  /* Constant Mg1 */
  DUMMY_ASSIGNMENT /*OC_STATE[46]*/ = 1000.0;
  /* Constant Mg2 */
  DUMMY_ASSIGNMENT /*OC_STATE[47]*/ = 1000.0;
  /* Constant Ca_CaT2 */
  DUMMY_ASSIGNMENT /*OC_STATE[48]*/ = 3.0;
  /* Constant D_0 */
  DUMMY_ASSIGNMENT /*OC_STATE[49]*/ = 0.8;
  /* Constant D_1 */
  DUMMY_ASSIGNMENT /*OC_STATE[50]*/ = 1.2;
  /* Constant D_2 */
  DUMMY_ASSIGNMENT /*OC_STATE[51]*/ = 3.0;
  /* Constant A_1 */
  DUMMY_ASSIGNMENT /*OC_STATE[52]*/ = 0.3;
  /* Constant A_2 */
  DUMMY_ASSIGNMENT /*OC_STATE[53]*/ = 0.23;
  /* Constant P */
  DUMMY_ASSIGNMENT /*OC_STATE[54]*/ = 0.23;
  /* Constant P_SR */
  DUMMY_ASSIGNMENT /*OC_STATE[55]*/ = 0.23;
  /* Constant P_C_SR */
  DUMMY_ASSIGNMENT /*OC_STATE[56]*/ = 0.23;
  /* Constant i2 */
  CONSTANTS[103] = 60;
  /* Constant V_o_Eqn */
  CONSTANTS[104] =
      0.950000 * CONSTANTS[65] * 3.14159265358979 * pow(CONSTANTS[66], 2.00000);
  /* Constant V_SR_Eqn */
  CONSTANTS[105] = 0.0500000 * CONSTANTS[65] * 3.14159265358979 *
                   pow(CONSTANTS[66], 2.00000);
  /* Constant V_1_Eqn */
  CONSTANTS[106] = 0.0100000 * CONSTANTS[104];
  /* Constant V_2_Eqn */
  CONSTANTS[107] = 0.990000 * CONSTANTS[104];
  /* Constant V_SR1_Eqn */
  CONSTANTS[108] = 0.0100000 * CONSTANTS[105];
  /* Constant V_SR2_Eqn */
  CONSTANTS[109] = 0.990000 * CONSTANTS[105];
  /* dCa1_Eqn */
  OC_RATE[29] =
      (((((CONSTANTS[103] * (OC_STATE[18] + OC_STATE[19] + OC_STATE[20] +
                             OC_STATE[21] + OC_STATE[22])) *
              ((OC_STATE[30] - OC_STATE[29]) / CONSTANTS[106]) -
          CONSTANTS[60] * ((OC_STATE[29] / (OC_STATE[29] + CONSTANTS[61])) /
                           CONSTANTS[106])) +
         CONSTANTS[62] * ((OC_STATE[30] - OC_STATE[29]) / CONSTANTS[106])) +
        -CONSTANTS[63] * ((OC_STATE[29] - OC_STATE[31]) / CONSTANTS[106])) +
       -((CONSTANTS[70] * OC_STATE[29]) *
             ((CONSTANTS[72] + -OC_STATE[34]) + -OC_STATE[36]) +
         -CONSTANTS[71] * OC_STATE[34])) +
      -((CONSTANTS[78] * OC_STATE[29]) * OC_STATE[44] +
        -CONSTANTS[79] * OC_STATE[40]);
  /* dCaSR1_Eqn */
  OC_RATE[30] =
      (((-(CONSTANTS[103] * (OC_STATE[18] + OC_STATE[19] + OC_STATE[20] +
                             OC_STATE[21] + OC_STATE[22])) *
             ((OC_STATE[30] - OC_STATE[29]) / CONSTANTS[108]) +
         CONSTANTS[60] * ((OC_STATE[29] / (OC_STATE[29] + CONSTANTS[61])) /
                          CONSTANTS[108])) +
        -CONSTANTS[62] * ((OC_STATE[30] - OC_STATE[29]) / CONSTANTS[108])) +
       -CONSTANTS[64] * ((OC_STATE[30] - OC_STATE[32]) / CONSTANTS[108])) +
      -((CONSTANTS[75] * OC_STATE[30]) * (CONSTANTS[77] - OC_STATE[38]) +
        -CONSTANTS[76] * OC_STATE[38]);
  /* dCa_SR2_Eqn */
  OC_RATE[32] =
      (((CONSTANTS[60] * ((OC_STATE[31] / (OC_STATE[31] + CONSTANTS[61])) /
                          CONSTANTS[109]) +
         -CONSTANTS[62] * ((OC_STATE[32] + -OC_STATE[31]) / CONSTANTS[109])) +
        CONSTANTS[64] * ((OC_STATE[30] + -OC_STATE[32]) / CONSTANTS[109])) +
       -((CONSTANTS[75] * OC_STATE[32]) * (CONSTANTS[77] + -OC_STATE[39]) +
         -CONSTANTS[76] * OC_STATE[39])) -
      (1000.00 / 1.00000) *
          (CONSTANTS[95] *
               (OC_STATE[55] * (0.00100000 / 1.00000) * OC_STATE[32] -
                CONSTANTS[97]) *
               (OC_STATE[55] * (0.00100000 / 1.00000) * OC_STATE[32] -
                            CONSTANTS[97] >
                        0.00000
                    ? 1.00000
                    : 0.00000) *
               (0.00100000 / 1.00000) * OC_STATE[55] * OC_STATE[32] -
           CONSTANTS[96] * OC_STATE[56] *
               (CONSTANTS[97] -
                OC_STATE[55] * (0.00100000 / 1.00000) * OC_STATE[32]) *
               (CONSTANTS[97] - OC_STATE[55] * (0.00100000 / 1.00000) *
                                    OC_STATE[32] >
                        0.00000
                    ? 1.00000
                    : 0.00000));
  /* dCa_P1_Eqn */
  OC_RATE[34] = (CONSTANTS[70] * OC_STATE[29]) *
                    ((CONSTANTS[72] + -OC_STATE[34]) + -OC_STATE[36]) +
                -CONSTANTS[71] * OC_STATE[34];
  /* dCa_P2_Eqn */
  OC_RATE[35] = (CONSTANTS[70] * OC_STATE[31]) *
                    ((CONSTANTS[72] + -OC_STATE[35]) + -OC_STATE[37]) +
                -CONSTANTS[71] * OC_STATE[35];
  /* dMg_P1_Eqn */
  OC_RATE[36] =
      (CONSTANTS[73] * (CONSTANTS[72] + -OC_STATE[34] + -OC_STATE[36])) *
          OC_STATE[46] +
      -CONSTANTS[74] * OC_STATE[36];
  /* dMP2_Eqn */
  OC_RATE[37] =
      (CONSTANTS[73] * (CONSTANTS[72] + -OC_STATE[35] + -OC_STATE[37])) *
          OC_STATE[47] +
      -CONSTANTS[74] * OC_STATE[37];
  /* dCa_Cs1_Eqn */
  OC_RATE[38] =
      (CONSTANTS[75] * OC_STATE[30]) * (CONSTANTS[77] + -OC_STATE[38]) +
      -CONSTANTS[76] * OC_STATE[38];
  /* dCs_Cs2_Eqn */
  OC_RATE[39] =
      (CONSTANTS[75] * OC_STATE[32]) * (CONSTANTS[77] + -OC_STATE[39]) +
      -CONSTANTS[76] * OC_STATE[39];
  /* dCa_ATP1_Eqn */
  OC_RATE[40] =
      ((CONSTANTS[78] * OC_STATE[29]) * OC_STATE[44] +
       -CONSTANTS[79] * OC_STATE[40]) +
      -CONSTANTS[82] * ((OC_STATE[40] + -OC_STATE[41]) / CONSTANTS[106]);
  /* dCa_ATP2_Eqn */
  OC_RATE[41] =
      ((CONSTANTS[78] * OC_STATE[31]) * OC_STATE[45] +
       -CONSTANTS[79] * OC_STATE[41]) +
      CONSTANTS[82] * ((OC_STATE[40] + -OC_STATE[41]) / CONSTANTS[107]);
  /* dMg_ATP1_Eqn */
  OC_RATE[42] =
      ((CONSTANTS[80] * OC_STATE[46]) * OC_STATE[44] +
       -CONSTANTS[81] * OC_STATE[42]) +
      -CONSTANTS[82] * ((OC_STATE[42] + -OC_STATE[43]) / CONSTANTS[106]);
  /* dMg_ATP2_Eqn */
  OC_RATE[43] =
      ((CONSTANTS[80] * OC_STATE[47]) * OC_STATE[45] +
       -CONSTANTS[81] * OC_STATE[43]) +
      CONSTANTS[82] * ((OC_STATE[42] + -OC_STATE[43]) / CONSTANTS[107]);
  /* dATP1_Eqn */
  OC_RATE[44] =
      (-((CONSTANTS[78] * OC_STATE[29]) * OC_STATE[44] +
         -CONSTANTS[79] * OC_STATE[40]) +
       -((CONSTANTS[80] * OC_STATE[46]) * OC_STATE[44] +
         -CONSTANTS[81] * OC_STATE[42])) +
      -CONSTANTS[82] * ((OC_STATE[44] + -OC_STATE[45]) / CONSTANTS[106]);
  /* dATP2_Eqn */
  OC_RATE[45] =
      (-((CONSTANTS[78] * OC_STATE[31]) * OC_STATE[45] +
         -CONSTANTS[79] * OC_STATE[41]) +
       -((CONSTANTS[80] * OC_STATE[47]) * OC_STATE[45] +
         -CONSTANTS[81] * OC_STATE[43])) +
      CONSTANTS[82] * ((OC_STATE[44] + -OC_STATE[45]) / CONSTANTS[107]);
  /* dMg1_Eqn */
  OC_RATE[46] =
      (-((CONSTANTS[73] * (CONSTANTS[72] + -OC_STATE[34] + -OC_STATE[36])) *
             OC_STATE[46] +
         -CONSTANTS[74] * OC_STATE[36]) +
       -((CONSTANTS[80] * OC_STATE[46]) * OC_STATE[44] +
         -CONSTANTS[81] * OC_STATE[42])) +
      -CONSTANTS[83] * ((OC_STATE[46] + -OC_STATE[47]) / CONSTANTS[106]);
  /* dMg2_Eqn */
  OC_RATE[47] =
      (-((CONSTANTS[73] * (CONSTANTS[72] + -OC_STATE[35] + -OC_STATE[37])) *
             OC_STATE[47] +
         -CONSTANTS[74] * OC_STATE[37]) +
       -((CONSTANTS[80] * OC_STATE[47]) * OC_STATE[45] +
         -CONSTANTS[81] * OC_STATE[43])) +
      CONSTANTS[83] * ((OC_STATE[46] + -OC_STATE[47]) / CONSTANTS[107]);
  /* dCa_CaT2_Eqn */
  OC_RATE[48] = (((CONSTANTS[67] * OC_STATE[31]) * OC_STATE[33] +
                  -CONSTANTS[68] * OC_STATE[48]) +
                 -CONSTANTS[86] * OC_STATE[48]) +
                CONSTANTS[87] * OC_STATE[51];
  /* dD_1_Eqn */
  OC_RATE[50] = ((((CONSTANTS[67] * OC_STATE[31] * OC_STATE[49] +
                    -CONSTANTS[68] * OC_STATE[50]) +
                   CONSTANTS[84] * OC_STATE[33]) +
                  -CONSTANTS[85] * OC_STATE[50]) +
                 (-CONSTANTS[67] * OC_STATE[31]) * OC_STATE[50]) +
                CONSTANTS[68] * OC_STATE[51];
  /* dD_2_Eqn */
  OC_RATE[51] = (((((CONSTANTS[67] * OC_STATE[31] * OC_STATE[50] +
                     -CONSTANTS[68] * OC_STATE[51]) +
                    CONSTANTS[86] * OC_STATE[48]) +
                   -CONSTANTS[87] * OC_STATE[51]) +
                  -CONSTANTS[88] * OC_STATE[51]) +
                 CONSTANTS[89] * OC_STATE[52]) +
                CONSTANTS[92] * OC_STATE[53];
  /* dA_1_Eqn */
  OC_RATE[52] =
      ((CONSTANTS[88] * OC_STATE[51] + -CONSTANTS[89] * OC_STATE[52]) +
       CONSTANTS[91] * OC_STATE[53]) +
      -CONSTANTS[90] * OC_STATE[52];
  /* dA_2_Eqn */
  OC_RATE[53] = (-CONSTANTS[91] * OC_STATE[53] + CONSTANTS[90] * OC_STATE[52]) +
                -CONSTANTS[92] * OC_STATE[53];
  /* dP_Eqn */
  OC_RATE[54] = (0.00100000 / 1.00000) * (CONSTANTS[90] * OC_STATE[52] -
                                          CONSTANTS[91] * OC_STATE[53]) +
                -1.00000 * CONSTANTS[93] * OC_STATE[54] +
                -1.00000 * CONSTANTS[94] *
                    ((OC_STATE[54] - OC_STATE[55]) / CONSTANTS[107]);
  /* dP_SR_Eqn */
  OC_RATE[55] =
      CONSTANTS[94] * ((OC_STATE[54] - OC_STATE[55]) / CONSTANTS[109]) -
      1.00000 * (CONSTANTS[95] *
                     (OC_STATE[55] * (0.00100000 / 1.00000) * OC_STATE[32] -
                      CONSTANTS[97]) *
                     (OC_STATE[55] * (0.00100000 / 1.00000) * OC_STATE[32] -
                                  CONSTANTS[97] >
                              0.00000
                          ? 1.00000
                          : 0.00000) *
                     (0.00100000 / 1.00000) * OC_STATE[55] * OC_STATE[32] -
                 CONSTANTS[96] * OC_STATE[56] *
                     (CONSTANTS[97] -
                      OC_STATE[55] * (0.00100000 / 1.00000) * OC_STATE[32]) *
                     (CONSTANTS[97] - OC_STATE[55] * (0.00100000 / 1.00000) *
                                          OC_STATE[32] >
                              0.00000
                          ? 1.00000
                          : 0.00000));
  /* dP_C_SR_Eqn */
  OC_RATE[56] =
      1.00000 * (CONSTANTS[95] *
                     (OC_STATE[55] * (0.00100000 / 1.00000) * OC_STATE[32] -
                      CONSTANTS[97]) *
                     (OC_STATE[55] * (0.00100000 / 1.00000) * OC_STATE[32] -
                                  CONSTANTS[97] >
                              0.00000
                          ? 1.00000
                          : 0.00000) *
                     (0.00100000 / 1.00000) * OC_STATE[55] * OC_STATE[32] -
                 CONSTANTS[96] * OC_STATE[56] *
                     (CONSTANTS[97] -
                      OC_STATE[55] * (0.00100000 / 1.00000) * OC_STATE[32]) *
                     (CONSTANTS[97] - OC_STATE[55] * (0.00100000 / 1.00000) *
                                          OC_STATE[32] >
                              0.00000
                          ? 1.00000
                          : 0.00000));
  /* T_0_Eqn */
  ALGEBRAIC[12] =
      (CONSTANTS[69] + -OC_STATE[33] + -OC_STATE[48] + -OC_STATE[49] +
                   -OC_STATE[50] + -OC_STATE[51] + -OC_STATE[52] +
                   -OC_STATE[53] >
               0.00000
           ? CONSTANTS[69] + -OC_STATE[33] + -OC_STATE[48] + -OC_STATE[49] +
                 -OC_STATE[50] + -OC_STATE[51] + -OC_STATE[52] + -OC_STATE[53]
           : 0.00000);
  /* dCa2_Eqn */
  OC_RATE[31] =
      ((((-CONSTANTS[60] * ((OC_STATE[31] / (OC_STATE[31] + CONSTANTS[61])) /
                            CONSTANTS[107]) +
          CONSTANTS[62] * ((OC_STATE[32] + -OC_STATE[31]) / CONSTANTS[107])) +
         CONSTANTS[63] * ((OC_STATE[29] - OC_STATE[31]) / CONSTANTS[107])) +
        -(((((((CONSTANTS[67] * OC_STATE[31] * ALGEBRAIC[12] +
                -CONSTANTS[68] * OC_STATE[33]) +
               CONSTANTS[67] * OC_STATE[31] * OC_STATE[33]) +
              -CONSTANTS[68] * OC_STATE[48]) +
             CONSTANTS[67] * OC_STATE[31] * OC_STATE[49]) +
            -CONSTANTS[68] * OC_STATE[50]) +
           CONSTANTS[67] * OC_STATE[31] * OC_STATE[50]) +
          -CONSTANTS[68] * OC_STATE[51])) +
       -((CONSTANTS[70] * OC_STATE[31]) *
             (CONSTANTS[72] + -OC_STATE[35] + -OC_STATE[37]) +
         -CONSTANTS[71] * OC_STATE[35])) +
      -((CONSTANTS[78] * OC_STATE[31]) * OC_STATE[45] +
        -CONSTANTS[79] * OC_STATE[41]);
  /* dCa_T_2_Eqn */
  OC_RATE[33] = (((((CONSTANTS[67] * OC_STATE[31]) * ALGEBRAIC[12] +
                    -CONSTANTS[68] * OC_STATE[33]) +
                   (-CONSTANTS[67] * OC_STATE[31]) * OC_STATE[33]) +
                  CONSTANTS[68] * OC_STATE[48]) +
                 -CONSTANTS[84] * OC_STATE[33]) +
                CONSTANTS[85] * OC_STATE[50];
  /* dD_0_Eqn */
  OC_RATE[49] = (((-CONSTANTS[67] * OC_STATE[31]) * OC_STATE[49] +
                  CONSTANTS[68] * OC_STATE[50]) +
                 CONSTANTS[84] * ALGEBRAIC[12]) +
                -CONSTANTS[85] * OC_STATE[49];
  /* stress_Eqn */
  OC_WANTED[0] =
      ((((OC_STATE[52] / CONSTANTS[69]) * CONSTANTS[99] +
         (OC_STATE[53] / CONSTANTS[69]) * CONSTANTS[100]) -
        CONSTANTS[101]) /
       CONSTANTS[102]) *
      (OC_KNOWN[1] >= 0.635000 && OC_KNOWN[1] <= 0.850000
           ? (0.700000 / (0.850000 - 0.635000)) * (OC_KNOWN[1] - 0.635000)
       : OC_KNOWN[1] > 0.850000 && OC_KNOWN[1] <= 1.17000
           ? 0.700000 +
                 (0.300000 / (1.17000 - 0.850000)) * (OC_KNOWN[1] - 0.850000)
       : OC_KNOWN[1] > 1.17000 && OC_KNOWN[1] <= 1.25500 ? 1.00000
       : OC_KNOWN[1] > 1.25500 && OC_KNOWN[1] <= 1.97000
           ? 1.00000 - (1.00000 / (1.97000 - 1.25500)) * (OC_KNOWN[1] - 1.25500)
           : 0.00000);
  /* dummy_Eqn */
  OC_RATE[28] = OC_WANTED[0];
  /* alpha_n_Eqn */
  ALGEBRAIC[1] =
      CONSTANTS[16] *
      ((OC_STATE[0] - CONSTANTS[21]) /
       (1.00000 - exp(-((OC_STATE[0] - CONSTANTS[21]) / CONSTANTS[32]))));
  /* beta_n_Eqn */
  ALGEBRAIC[14] =
      CONSTANTS[19] * exp(-((OC_STATE[0] - CONSTANTS[21]) / CONSTANTS[34]));
  /* dn_Eqn */
  OC_RATE[8] =
      ALGEBRAIC[1] * (1.00000 - OC_STATE[8]) - ALGEBRAIC[14] * OC_STATE[8];
  /* h_K_inf_Eqn */
  ALGEBRAIC[2] =
      1.00000 / (1.00000 + exp((OC_STATE[0] - CONSTANTS[25]) / CONSTANTS[28]));
  /* tau_h_K_Eqn */
  ALGEBRAIC[15] = 1000.00 * exp(-((OC_STATE[0] + 40.0000) / 25.7500));
  /* dh_K_Eqn */
  OC_RATE[9] = (ALGEBRAIC[2] - OC_STATE[9]) / ALGEBRAIC[15];
  /* alpha_m_Eqn */
  ALGEBRAIC[4] =
      CONSTANTS[15] *
      ((OC_STATE[0] - CONSTANTS[20]) /
       (1.00000 - exp(-((OC_STATE[0] - CONSTANTS[20]) / CONSTANTS[31]))));
  /* beta_m_Eqn */
  ALGEBRAIC[17] =
      CONSTANTS[18] * exp(-((OC_STATE[0] - CONSTANTS[20]) / CONSTANTS[33]));
  /* dm_Eqn */
  OC_RATE[10] =
      ALGEBRAIC[4] * (1.00000 - OC_STATE[10]) - ALGEBRAIC[17] * OC_STATE[10];
  /* alpha_h_Eqn */
  ALGEBRAIC[3] =
      CONSTANTS[14] * exp(-((OC_STATE[0] - CONSTANTS[22]) / CONSTANTS[29]));
  /* beta_h_Eqn */
  ALGEBRAIC[16] =
      CONSTANTS[17] /
      (1.00000 + exp(-((OC_STATE[0] - CONSTANTS[22]) / CONSTANTS[30])));
  /* dh_Eqn */
  OC_RATE[11] =
      ALGEBRAIC[3] * (1.00000 - OC_STATE[11]) - ALGEBRAIC[16] * OC_STATE[11];
  /* S_inf_Eqn */
  ALGEBRAIC[5] =
      1.00000 / (1.00000 + exp((OC_STATE[0] - CONSTANTS[24]) / CONSTANTS[27]));
  /* tau_S_Eqn */
  ALGEBRAIC[18] =
      8571.00 /
      (0.200000 +
       5.65000 * pow((OC_STATE[0] + CONSTANTS[48]) / 100.000, 2.00000));
  /* dS_Eqn */
  OC_RATE[12] = (ALGEBRAIC[5] - OC_STATE[12]) / ALGEBRAIC[18];
  /* alpha_n_t_Eqn */
  ALGEBRAIC[6] =
      CONSTANTS[16] *
      ((OC_STATE[1] - CONSTANTS[21]) /
       (1.00000 - exp(-((OC_STATE[1] - CONSTANTS[21]) / CONSTANTS[32]))));
  /* beta_n_t_Eqn */
  ALGEBRAIC[19] =
      CONSTANTS[19] * exp(-((OC_STATE[1] - CONSTANTS[21]) / CONSTANTS[34]));
  /* dn_t_Eqn */
  OC_RATE[13] =
      ALGEBRAIC[6] * (1.00000 - OC_STATE[13]) - ALGEBRAIC[19] * OC_STATE[13];
  /* h_K_inf_t_Eqn */
  ALGEBRAIC[7] =
      1.00000 / (1.00000 + exp((OC_STATE[1] - CONSTANTS[25]) / CONSTANTS[28]));
  /* tau_h_K_t_Eqn */
  ALGEBRAIC[20] = 1.00000 * exp(-((OC_STATE[1] + 40.0000) / 25.7500));
  /* dh_K_t_Eqn */
  OC_RATE[14] = (ALGEBRAIC[7] - OC_STATE[14]) / ALGEBRAIC[20];
  /* alpha_m_t_Eqn */
  ALGEBRAIC[9] =
      CONSTANTS[15] *
      ((OC_STATE[1] - CONSTANTS[20]) /
       (1.00000 - exp(-((OC_STATE[1] - CONSTANTS[20]) / CONSTANTS[31]))));
  /* beta_m_t_Eqn */
  ALGEBRAIC[22] =
      CONSTANTS[18] * exp(-((OC_STATE[1] - CONSTANTS[20]) / CONSTANTS[33]));
  /* dm_t_Eqn */
  OC_RATE[15] =
      ALGEBRAIC[9] * (1.00000 - OC_STATE[15]) - ALGEBRAIC[22] * OC_STATE[15];
  /* alpha_h_t_Eqn */
  ALGEBRAIC[8] =
      CONSTANTS[14] * exp(-((OC_STATE[1] - CONSTANTS[22]) / CONSTANTS[29]));
  /* beta_h_t_Eqn */
  ALGEBRAIC[21] =
      CONSTANTS[17] /
      (1.00000 + exp(-((OC_STATE[1] - CONSTANTS[22]) / CONSTANTS[30])));
  /* dh_t_Eqn */
  OC_RATE[16] =
      ALGEBRAIC[8] * (1.00000 - OC_STATE[16]) - ALGEBRAIC[21] * OC_STATE[16];
  /* S_inf_t_Eqn */
  ALGEBRAIC[10] =
      1.00000 / (1.00000 + exp((OC_STATE[1] - CONSTANTS[24]) / CONSTANTS[27]));
  /* tau_S_t_Eqn */
  ALGEBRAIC[23] =
      8571.00 /
      (0.200000 +
       5.65000 * pow((OC_STATE[1] + CONSTANTS[48]) / 100.000, 2.00000));
  /* dS_t_Eqn */
  OC_RATE[17] = (ALGEBRAIC[10] - OC_STATE[17]) / ALGEBRAIC[23];
  /* k_C_Eqn */
  ALGEBRAIC[11] =
      0.500000 * CONSTANTS[57] *
      exp((OC_STATE[1] - CONSTANTS[59]) / (8.00000 * CONSTANTS[58]));
  /* k_Cm_Eqn */
  ALGEBRAIC[24] =
      0.500000 * CONSTANTS[57] *
      exp((CONSTANTS[59] - OC_STATE[1]) / (8.00000 * CONSTANTS[58]));
  /* dC_0_Eqn */
  OC_RATE[23] = -CONSTANTS[54] * OC_STATE[23] + CONSTANTS[55] * OC_STATE[18] +
                -4.00000 * ALGEBRAIC[11] * OC_STATE[23] +
                ALGEBRAIC[24] * OC_STATE[24];
  /* dO_0_Eqn */
  OC_RATE[18] = CONSTANTS[54] * OC_STATE[23] + -CONSTANTS[55] * OC_STATE[18] +
                (-4.00000 * ALGEBRAIC[11] * OC_STATE[18]) / CONSTANTS[56] +
                CONSTANTS[56] * ALGEBRAIC[24] * OC_STATE[19];
  /* dC_1_Eqn */
  OC_RATE[24] = 4.00000 * ALGEBRAIC[11] * OC_STATE[23] +
                -ALGEBRAIC[24] * OC_STATE[24] +
                (-CONSTANTS[54] * OC_STATE[24]) / CONSTANTS[56] +
                CONSTANTS[56] * CONSTANTS[55] * OC_STATE[19] +
                -3.00000 * ALGEBRAIC[11] * OC_STATE[24] +
                2.00000 * ALGEBRAIC[24] * OC_STATE[25];
  /* dO_1_Eqn */
  OC_RATE[19] = (CONSTANTS[54] * OC_STATE[24]) / CONSTANTS[56] +
                -CONSTANTS[55] * CONSTANTS[56] * OC_STATE[19] +
                (4.00000 * ALGEBRAIC[11] * OC_STATE[18]) / CONSTANTS[56] +
                -CONSTANTS[56] * ALGEBRAIC[24] * OC_STATE[19] +
                (-3.00000 * ALGEBRAIC[11] * OC_STATE[19]) / CONSTANTS[56] +
                2.00000 * CONSTANTS[56] * ALGEBRAIC[24] * OC_STATE[20];
  /* dC_2_Eqn */
  OC_RATE[25] = 3.00000 * ALGEBRAIC[11] * OC_STATE[24] +
                -2.00000 * ALGEBRAIC[24] * OC_STATE[25] +
                (-CONSTANTS[54] * OC_STATE[25]) / pow(CONSTANTS[56], 2.00000) +
                pow(CONSTANTS[56], 2.00000) * CONSTANTS[55] * OC_STATE[20] +
                -2.00000 * ALGEBRAIC[11] * OC_STATE[25] +
                3.00000 * ALGEBRAIC[24] * OC_STATE[26];
  /* dO_2_Eqn */
  OC_RATE[20] = (3.00000 * ALGEBRAIC[11] * OC_STATE[19]) / CONSTANTS[56] +
                -2.00000 * CONSTANTS[56] * ALGEBRAIC[24] * OC_STATE[20] +
                (CONSTANTS[54] * OC_STATE[25]) / pow(CONSTANTS[56], 2.00000) +
                -CONSTANTS[55] * pow(CONSTANTS[56], 2.00000) * OC_STATE[20] +
                (-2.00000 * ALGEBRAIC[11] * OC_STATE[20]) / CONSTANTS[56] +
                3.00000 * CONSTANTS[56] * ALGEBRAIC[24] * OC_STATE[21];
  /* dC_3_Eqn */
  OC_RATE[26] = 2.00000 * ALGEBRAIC[11] * OC_STATE[25] +
                -3.00000 * ALGEBRAIC[24] * OC_STATE[26] +
                (-CONSTANTS[54] * OC_STATE[26]) / pow(CONSTANTS[56], 3.00000) +
                CONSTANTS[55] * pow(CONSTANTS[56], 3.00000) * OC_STATE[21] +
                -ALGEBRAIC[11] * OC_STATE[26] +
                4.00000 * ALGEBRAIC[24] * OC_STATE[27];
  /* dO_3_Eqn */
  OC_RATE[21] = (CONSTANTS[54] * OC_STATE[26]) / pow(CONSTANTS[56], 3.00000) +
                -CONSTANTS[55] * pow(CONSTANTS[56], 3.00000) * OC_STATE[21] +
                (2.00000 * ALGEBRAIC[11] * OC_STATE[20]) / CONSTANTS[56] +
                -3.00000 * ALGEBRAIC[24] * CONSTANTS[56] * OC_STATE[21] +
                (-ALGEBRAIC[11] * OC_STATE[21]) / CONSTANTS[56] +
                4.00000 * CONSTANTS[56] * ALGEBRAIC[24] * OC_STATE[22];
  /* dC_4_Eqn */
  OC_RATE[27] = ALGEBRAIC[11] * OC_STATE[26] +
                -4.00000 * ALGEBRAIC[24] * OC_STATE[27] +
                (-CONSTANTS[54] * OC_STATE[27]) / pow(CONSTANTS[56], 4.00000) +
                CONSTANTS[55] * pow(CONSTANTS[56], 4.00000) * OC_STATE[22];
  /* dO_4_Eqn */
  OC_RATE[22] = (ALGEBRAIC[11] * OC_STATE[21]) / CONSTANTS[56] +
                -4.00000 * CONSTANTS[56] * ALGEBRAIC[24] * OC_STATE[22] +
                (CONSTANTS[54] * OC_STATE[27]) / pow(CONSTANTS[56], 4.00000) +
                -CONSTANTS[55] * pow(CONSTANTS[56], 4.00000) * OC_STATE[22];
  /* J_K_eqn */
  ALGEBRAIC[30] = OC_STATE[0] *
                  ((OC_STATE[3] -
                    OC_STATE[4] * exp((-1.00000 * CONSTANTS[6] * OC_STATE[0]) /
                                      (CONSTANTS[35] * CONSTANTS[36]))) /
                   (1.00000 - exp((-1.00000 * CONSTANTS[6] * OC_STATE[0]) /
                                  (CONSTANTS[35] * CONSTANTS[36]))));
  /* E_K_eqn */
  ALGEBRAIC[13] = ((CONSTANTS[35] * CONSTANTS[36]) / CONSTANTS[6]) *
                  log(OC_STATE[4] / OC_STATE[3]);
  /* K_R_Eqn */
  ALGEBRAIC[36] =
      OC_STATE[4] * exp((-CONSTANTS[41] * ALGEBRAIC[13]) *
                        (CONSTANTS[6] / (CONSTANTS[35] * CONSTANTS[36])));
  /* g_IR_bar_Eqn */
  ALGEBRAIC[37] =
      CONSTANTS[40] * (pow(ALGEBRAIC[36], 2.00000) /
                       (CONSTANTS[42] + pow(ALGEBRAIC[36], 2.00000)));
  /* y_Eqn */
  ALGEBRAIC[38] =
      1.00000 -
      pow(1.00000 + (CONSTANTS[43] *
                     (1.00000 + pow(ALGEBRAIC[36], 2.00000) / CONSTANTS[42])) /
                        (pow(CONSTANTS[46], 2.00000) *
                         exp((2.00000 * (1.00000 - CONSTANTS[41]) *
                              OC_STATE[0] * CONSTANTS[6]) /
                             (CONSTANTS[35] * CONSTANTS[36]))),
          -1.00000);
  /* g_IR_Eqn */
  ALGEBRAIC[39] = ALGEBRAIC[37] * ALGEBRAIC[38];
  /* I_IR_Eqn */
  ALGEBRAIC[40] = ALGEBRAIC[39] *
                  (ALGEBRAIC[30] > 0.00000 ? 1.00000 : 0.00000) *
                  (ALGEBRAIC[30] / 50.0000);
  /* g_DR_Eqn */
  ALGEBRAIC[41] = (CONSTANTS[38] * pow(OC_STATE[8], 4.00000)) * OC_STATE[9];
  /* I_DR_Eqn */
  ALGEBRAIC[42] = ALGEBRAIC[41] * (ALGEBRAIC[30] / 50.0000);
  /* sig_Eqn */
  ALGEBRAIC[46] = (1.00000 / 7.00000) * (exp(OC_STATE[7] / 67.3000) - 1.00000);
  /* f1_Eqn */
  ALGEBRAIC[47] =
      pow(1.00000 +
              0.120000 * exp(-0.100000 * OC_STATE[0] *
                             (CONSTANTS[6] / (CONSTANTS[35] * CONSTANTS[36]))) +
              0.0400000 * ALGEBRAIC[46] *
                  exp(-(OC_STATE[0] *
                        (CONSTANTS[6] / (CONSTANTS[35] * CONSTANTS[36])))),
          -1.00000);
  /* I_NaK_bar_Eqn */
  ALGEBRAIC[48] =
      CONSTANTS[6] *
      (CONSTANTS[47] / (pow(1.00000 + CONSTANTS[44] / OC_STATE[4], 2.00000) *
                        pow(1.00000 + CONSTANTS[45] / OC_STATE[5], 3.00000)));
  /* I_NaK_Eqn */
  ALGEBRAIC[49] = ALGEBRAIC[48] * ALGEBRAIC[47];
  /* dK_e_Eqn */
  OC_RATE[4] = (ALGEBRAIC[40] + ALGEBRAIC[42] + CONSTANTS[12] +
                -2.00000 * ALGEBRAIC[49]) /
                   ((1000.00 / 1.00000) * CONSTANTS[6] * CONSTANTS[5]) +
               (OC_STATE[2] - OC_STATE[4]) / CONSTANTS[10];
  /* g_Na_Eqn */
  ALGEBRAIC[44] =
      ((CONSTANTS[39] * pow(OC_STATE[10], 3.00000)) * OC_STATE[11]) *
      OC_STATE[12];
  /* J_Na_eqn */
  ALGEBRAIC[43] = OC_STATE[0] *
                  ((OC_STATE[5] -
                    OC_STATE[7] * exp((-1.00000 * CONSTANTS[6] * OC_STATE[0]) /
                                      (CONSTANTS[35] * CONSTANTS[36]))) /
                   (1.00000 - exp((-1.00000 * CONSTANTS[6] * OC_STATE[0]) /
                                  (CONSTANTS[35] * CONSTANTS[36]))));
  /* I_Na_Eqn */
  ALGEBRAIC[45] = ALGEBRAIC[44] * (ALGEBRAIC[43] / 75.0000);
  /* dNa_e_Eqn */
  OC_RATE[7] = (ALGEBRAIC[45] + CONSTANTS[13] + 3.00000 * ALGEBRAIC[49]) /
                   ((1000.00 / 1.00000) * CONSTANTS[6] * CONSTANTS[5]) +
               (OC_STATE[6] - OC_STATE[7]) / CONSTANTS[11];
  /* I_T_Eqn */
  ALGEBRAIC[0] =
      (1000.00 / 1.00000) * ((OC_STATE[0] - OC_STATE[1]) / CONSTANTS[2]);
  /* Cl_i_eqn */
  ALGEBRAIC[26] = 156.500 / (5.00000 + exp((-CONSTANTS[6] * ALGEBRAIC[13]) /
                                           (CONSTANTS[35] * CONSTANTS[36])));
  /* Cl_o_eqn */
  ALGEBRAIC[27] = 156.500 - 5.00000 * ALGEBRAIC[26];
  /* J_Cl_eqn */
  ALGEBRAIC[33] =
      OC_STATE[0] *
      ((ALGEBRAIC[26] - ALGEBRAIC[27] * exp((CONSTANTS[6] * OC_STATE[0]) /
                                            (CONSTANTS[35] * CONSTANTS[36]))) /
       (1.00000 -
        exp((CONSTANTS[6] * OC_STATE[0]) / (CONSTANTS[35] * CONSTANTS[36]))));
  /* a_Eqn */
  ALGEBRAIC[32] =
      1.00000 / (1.00000 + exp((OC_STATE[0] - CONSTANTS[23]) / CONSTANTS[26]));
  /* g_Cl_Eqn */
  ALGEBRAIC[34] = CONSTANTS[37] * pow(ALGEBRAIC[32], 4.00000);
  /* I_Cl_Eqn */
  ALGEBRAIC[35] = ALGEBRAIC[34] * (ALGEBRAIC[33] / 45.0000);
  /* I_ionic_s_Eqn */
  ALGEBRAIC[50] = ALGEBRAIC[35] + ALGEBRAIC[40] + ALGEBRAIC[42] +
                  ALGEBRAIC[45] + ALGEBRAIC[49] + -OC_KNOWN[0];
  /* vS_diff_calculation */
  OC_RATE[0] = -((ALGEBRAIC[50] + ALGEBRAIC[0]) / CONSTANTS[0]);
  /* J_K_t_eqn */
  ALGEBRAIC[31] = OC_STATE[1] *
                  ((OC_STATE[3] -
                    OC_STATE[2] * exp((-1.00000 * CONSTANTS[6] * OC_STATE[1]) /
                                      (CONSTANTS[35] * CONSTANTS[36]))) /
                   (1.00000 - exp((-1.00000 * CONSTANTS[6] * OC_STATE[1]) /
                                  (CONSTANTS[35] * CONSTANTS[36]))));
  /* E_K_t_eqn */
  ALGEBRAIC[25] = ((CONSTANTS[35] * CONSTANTS[36]) / CONSTANTS[6]) *
                  log(OC_STATE[2] / OC_STATE[3]);
  /* K_R_t_Eqn */
  ALGEBRAIC[55] =
      OC_STATE[2] * exp((-CONSTANTS[41] * ALGEBRAIC[25]) *
                        (CONSTANTS[6] / (CONSTANTS[35] * CONSTANTS[36])));
  /* g_IR_bar_t_Eqn */
  ALGEBRAIC[56] =
      CONSTANTS[40] * (pow(ALGEBRAIC[55], 2.00000) /
                       (CONSTANTS[42] + pow(ALGEBRAIC[55], 2.00000)));
  /* y_t_Eqn */
  ALGEBRAIC[57] =
      1.00000 -
      pow(1.00000 + (CONSTANTS[43] *
                     (1.00000 + pow(ALGEBRAIC[55], 2.00000) / CONSTANTS[42])) /
                        (pow(CONSTANTS[46], 2.00000) *
                         exp((2.00000 * (1.00000 - CONSTANTS[41]) *
                              OC_STATE[1] * CONSTANTS[6]) /
                             (CONSTANTS[35] * CONSTANTS[36]))),
          -1.00000);
  /* g_IR_t_Eqn */
  ALGEBRAIC[58] = ALGEBRAIC[56] * ALGEBRAIC[57];
  /* I_IR_t_Eqn */
  ALGEBRAIC[59] = CONSTANTS[50] * ALGEBRAIC[58] * (ALGEBRAIC[31] / 50.0000);
  /* g_DR_t_Eqn */
  ALGEBRAIC[60] = (CONSTANTS[38] * pow(OC_STATE[13], 4.00000)) * OC_STATE[14];
  /* I_DR_t_Eqn */
  ALGEBRAIC[61] = CONSTANTS[51] * ALGEBRAIC[60] * (ALGEBRAIC[31] / 50.0000);
  /* sig_t_Eqn */
  ALGEBRAIC[65] = (1.00000 / 7.00000) * (exp(OC_STATE[6] / 67.3000) - 1.00000);
  /* f1_t_Eqn */
  ALGEBRAIC[66] =
      pow(1.00000 +
              0.120000 * exp(-0.100000 * OC_STATE[1] *
                             (CONSTANTS[6] / (CONSTANTS[35] * CONSTANTS[36]))) +
              0.0400000 * ALGEBRAIC[65] *
                  exp(-(OC_STATE[1] *
                        (CONSTANTS[6] / (CONSTANTS[35] * CONSTANTS[36])))),
          -1.00000);
  /* I_NaK_bar_t_Eqn */
  ALGEBRAIC[67] =
      CONSTANTS[6] *
      (CONSTANTS[47] / (pow(1.00000 + CONSTANTS[44] / OC_STATE[2], 2.00000) *
                        pow(1.00000 + CONSTANTS[45] / OC_STATE[5], 3.00000)));
  /* I_NaK_t_Eqn */
  ALGEBRAIC[68] = CONSTANTS[53] * ALGEBRAIC[67] * ALGEBRAIC[66];
  /* dK_i_Eqn */
  OC_RATE[3] =
      -CONSTANTS[9] * ((ALGEBRAIC[59] + ALGEBRAIC[61] + CONSTANTS[12] +
                        -2.00000 * ALGEBRAIC[68]) /
                       ((1000.00 / 1.00000) * CONSTANTS[6] * CONSTANTS[3])) -
      (ALGEBRAIC[40] + ALGEBRAIC[42] + CONSTANTS[12] +
       -2.00000 * ALGEBRAIC[49]) /
          ((1000.00 / 1.00000) * CONSTANTS[6] * CONSTANTS[4]);
  /* dK_t_Eqn */
  OC_RATE[2] = (ALGEBRAIC[59] + ALGEBRAIC[61] + CONSTANTS[12] +
                -2.00000 * ALGEBRAIC[68]) /
                   ((1000.00 / 1.00000) * CONSTANTS[6] * CONSTANTS[3]) -
               (OC_STATE[2] - OC_STATE[4]) / CONSTANTS[7];
  /* g_Na_t_Eqn */
  ALGEBRAIC[63] =
      ((CONSTANTS[39] * pow(OC_STATE[15], 3.00000)) * OC_STATE[16]) *
      OC_STATE[17];
  /* J_Na_t_eqn */
  ALGEBRAIC[62] = OC_STATE[1] *
                  ((OC_STATE[5] -
                    OC_STATE[6] * exp((-1.00000 * CONSTANTS[6] * OC_STATE[1]) /
                                      (CONSTANTS[35] * CONSTANTS[36]))) /
                   (1.00000 - exp((-1.00000 * CONSTANTS[6] * OC_STATE[1]) /
                                  (CONSTANTS[35] * CONSTANTS[36]))));
  /* I_Na_t_Eqn */
  ALGEBRAIC[64] = CONSTANTS[52] * ALGEBRAIC[63] * (ALGEBRAIC[62] / 75.0000);
  /* dNa_i_Eqn */
  OC_RATE[5] = -CONSTANTS[9] *
                   ((ALGEBRAIC[64] + CONSTANTS[13] + 3.00000 * ALGEBRAIC[68]) /
                    ((1000.00 / 1.00000) * CONSTANTS[6] * CONSTANTS[3])) -
               (ALGEBRAIC[45] + CONSTANTS[13] + 3.00000 * ALGEBRAIC[49]) /
                   ((1000.00 / 1.00000) * CONSTANTS[6] * CONSTANTS[4]);
  /* dNa_t_Eqn */
  OC_RATE[6] = (ALGEBRAIC[64] + CONSTANTS[13] + 3.00000 * ALGEBRAIC[68]) /
                   ((1000.00 / 1.00000) * CONSTANTS[6] * CONSTANTS[3]) -
               (OC_STATE[6] - OC_STATE[7]) / CONSTANTS[8];
  /* Cl_i_t_eqn */
  ALGEBRAIC[28] = 156.500 / (5.00000 + exp((-CONSTANTS[6] * ALGEBRAIC[25]) /
                                           (CONSTANTS[35] * CONSTANTS[36])));
  /* Cl_o_t_eqn */
  ALGEBRAIC[29] = 156.500 - 5.00000 * ALGEBRAIC[28];
  /* J_Cl_t_eqn */
  ALGEBRAIC[52] =
      OC_STATE[1] *
      ((ALGEBRAIC[28] - ALGEBRAIC[29] * exp((CONSTANTS[6] * OC_STATE[1]) /
                                            (CONSTANTS[35] * CONSTANTS[36]))) /
       (1.00000 -
        exp((CONSTANTS[6] * OC_STATE[1]) / (CONSTANTS[35] * CONSTANTS[36]))));
  /* a_t_Eqn */
  ALGEBRAIC[51] =
      1.00000 / (1.00000 + exp((OC_STATE[1] - CONSTANTS[23]) / CONSTANTS[26]));
  /* g_Cl_t_Eqn */
  ALGEBRAIC[53] = CONSTANTS[37] * pow(ALGEBRAIC[51], 4.00000);
  /* I_Cl_t_Eqn */
  ALGEBRAIC[54] = CONSTANTS[49] * ALGEBRAIC[53] * (ALGEBRAIC[52] / 45.0000);
  /* I_ionic_t_Eqn */
  ALGEBRAIC[69] = ALGEBRAIC[54] + ALGEBRAIC[59] + ALGEBRAIC[61] +
                  ALGEBRAIC[64] + ALGEBRAIC[68];
  /* dvT_Eqn */
  OC_RATE[1] = -((ALGEBRAIC[69] - ALGEBRAIC[0] / CONSTANTS[1]) / CONSTANTS[0]);

} // OC_CellML_RHS_routine()

;