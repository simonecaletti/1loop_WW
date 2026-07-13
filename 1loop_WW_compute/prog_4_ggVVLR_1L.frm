#-

#:maxtermsize 1000000
#:workspace 1G


* 1) Compute the helicity coefficients
**************************************************************************************************************** Declare symbols for names
autodeclare s d, x, y, PL, NP, PLA, PLB, PLC, NPA, NPB, Sec, n, box0, [d-4], [d-3], [d-2], [d-1]                        ;
**************************************************************************************************************** Declare kinematic variables
autodeclare s s, u                                                                                                        ;
s t, mt, mw2 ,m32, m42, pi, eps, S, U                                                                               ;
**************************************************************************************************************** Declare couplings colour factors and other parameters for the amplitude
s vq, aq, g, Nc, NA, CF, d, Nf, Nh, NfV, NhV, NfA, NfV2, NhV2, NfA2[], [Nc^2-1], [Nc-1], [gh], I, XXX, NffcV, NffcA     ;
**************************************************************************************************************** Declare some extra symbols for numbers
s n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15                                                      ;
**************************************************************************************************************** Declare indices
autodeclare i mu=d, nu=d, ro=d, sigma=d, a=NA, b=NA, i=Nc, j=Nc , J=Nc                                                  ;
**************************************************************************************************************** Declare vectors for momenta and polarisations
autodeclare v k, mk, p, mp, q, va                                                                                  ;
*************************************************************************************************************** Declare extrasymbol for gamma_5
s a, x                                                                                                                ;

**************************************************************************************************************** Declare functions for particle names
cf Quark1, Quark2, Quark3, Qbar1, Qbar2, Qbar3, gluon, ghost, Ghost, Vboson, Vboson1, Vboson2                                          ;
**************************************************************************************************************** Declare functions for propagators and vertices
cf prop, propmt, vrtx, Vqqbv, Vqqba, Vqqbg, Vttbv, Vttba, Vttbg, Vggg, Vgggg, Vgggg1, Vgggg2, Vgggg3, VghGgh, VHpb, Vvvv;
**************************************************************************************************************** Declare functions for gamma matrices and spinor bilinears
cf gaQ, Bi, BiVU, BiUV                                                                                                  ;
**************************************************************************************************************** Declare functions for polarisations:
*cf pol, U, V, Ub, Vb                                                                                                    ;
**************************************************************************************************************** Declare polarisation vectors as functions for individual particles
*cf Eps1,Eps2,Eps3,Eps4,Eps5                                                                                             ;
**************************************************************************************************************** Declare colour objects
cf T,f,Dd,Tr,delC, TR                                                                                                       ;
**************************************************************************************************************** Declare loops
cf LoopF,LoopGh, PL00M, PL00Mx12, PL0M0, PL0M0x12, NPL0M, NPL0Mx12                                                       ;
**************************************************************************************************************** Declare functions for factors
cf Num, Den, coeff, Ratio                                                                                               ;
**************************************************************************************************************** Declare functions for write helicity form factors
cf FtoHe                                                                                                                ;
**************************************************************************************************************** Declare label type objects
cf Label, INT, DiaMatch, Shift, Sector, Save, ACC,                                                                      ;

**************************************************************************************************************** Declare functions for spinor helicity calculus
cf spa, sps, spaInv, spsInv, spas, spsa, spaa, spss                                                                     ;

**************************************************************************************************************** Declare objects for projectors and form factors
autodeclare s c, F, He, PH                                                                                              ;

**************************************************************************************************************** Declare propagators
s D1,...,D15                                                                                                            ;
**************************************************************************************************************** Declare symbols for one loop reduction to masters
cf INT1L, M                                                                                                             ;
s P12, P13, P13x12, root, heps                                                                                          ;

s AmpT1d, AmpT2d, AmpT3d, AmpT4d, AmpT5d, AmpT6d, AmpT7d, AmpT8d, AmpT9d, AmpT10d, AmpT11d, AmpT12d, AmpT13d, AmpT14d, AmpT15d, AmpT16d, AmpT17d, AmpT18d, AmpT19d, AmpT20d, AmpT21d, AmpT22d, AmpT23d, AmpT24d, AmpT25d, AmpT26d, AmpT27d, AmpT28d, AmpT29d, AmpT30d, AmpT31d, AmpT32d, AmpT33d, AmpT34d, AmpT35d, AmpT36d ;
**************************************************************************************************************** Declare propagator symbols
#include ../Procedures/Procedures.proc
#include ../Procedures/Procedures_2.proc



off statistics ;
.global
#message Importing the projectors coefficients and the change of basis coefficient



*#do i1=1,36
*g [F'i1'] =
*+ c'i2'p'i1' * [AmpT'i2'd]
*#enddo
*;
*#enddo
*.sort

#do i1=1,18
  g [He'i1'] =
  #do i2=1,36
  #do i3=1,36
   + FtoHe('i1','i2')*c'i3'p'i2' * AmpT'i3'd
  #enddo
  #enddo
  ;
#enddo
.sort

#include ../Tensor_decomposition/projectorscoeffs.inc
    .sort
#include ../Tensor_decomposition/helicity_coeff_LR.inc
.sort


#message calling polyratfun to simplify the expressions


#call FactArgDen

id Den (-1+d) = 1/[d-1];
id Den(-3+d)  = 1/[d-3];
id Den(-2+d)  = 1/[d-2];
id Den(d-4)   = 1/[d-4];


b  Den, Num, Ratio, d, [d-1],[d-2], [d-3], [d-4], AmpT1d, AmpT2d, AmpT3d, AmpT4d, AmpT5d, AmpT6d, AmpT7d, AmpT8d, AmpT9d, AmpT10d, AmpT11d, AmpT12d, AmpT13d, AmpT14d, AmpT15d, AmpT16d, AmpT17d, AmpT18d, AmpT19d, AmpT20d, AmpT21d, AmpT22d, AmpT23d, AmpT24d, AmpT25d, AmpT26d, AmpT27d, AmpT28d, AmpT29d, AmpT30d, AmpT31d, AmpT32d, AmpT33d, AmpT34d, AmpT35d, AmpT36d;
.sort

Collect Num;
.sort
#call PFD

#message polyratfun in sij and mj
id Num(n1?) = Ratio(n1,1)            ;
id Den(n1?) = Ratio(1,n1)            ;
.sort
polyratfun Ratio                     ;
.sort
polyratfun                           ;
.sort

id Ratio(n1?,n2?) =  Num(n1) * Den(n2) ;
.sort
#call FactArgNum
#call FactArgDen
.sort

#do i=1,36
#include Tensors_x_amp/Simplified_AmpT'i'd.inc
#enddo;
.sort

*b i_, d_, CF, Nf, Nh, NfV, NhV, NfV2, NhV2, g, vq, aq, Nc, Den, Num, Ratio, NffcV, NffcA, Tr, d, [d-1], [d-2], [d-3], [d-4], PL00M, PL00Mx12, PL0M0, PL0M0x12, NPL0M, NPL0Mx12;
*.sort
*Collect Num;
*.sort
#call PFD

#message polyratfun in sij and mj
id Num(n1?) = Ratio(n1,1)            ;
id Den(n1?) = Ratio(1,n1)            ;
.sort
polyratfun Ratio                     ;
.sort
polyratfun                           ;
.sort

id Ratio(n1?,n2?) =  Num(n1) * Den(n2) ;
.sort
#call FactArgNum
#call FactArgDen
.sort

#message Exporting into file.inc

#do i=1,18
#write <Helicity_form_factors/He'i'LR.inc> "g He'i' = %E ; " [He'i']
#enddo

#message Export in Mathematica

*Perform epsilon expansion for the d dimensional factors

id [d-1]^n1?    =  3^(n1)-2*eps*n1*3^(n1-1) + 2*eps^2*3^(n1-2)*n1*(n1-1);
id [d-2]^n1?    =  eps^2*2^(n1-1)*(n1-1)*n1-eps*2^n1*n1+2^n1;
id [d-3]^n1?    =  2*eps^2*(n1^2-n1)-2*eps*n1+1;
id [d-4]^n1?    =  2^n1* (-eps)^(n1);
    .sort

*AG: change notation for S@M

id s^n1? = S^n1          ;
id u^n1? = U^n1          ;
    .sort
argument Num, Den ;
id s^n1? = S^n1          ;
id u^n1? = U^n1          ;
endargument       ;
    .sort

id Tr(a99,a97) = TR;
    .sort
#do i=1,18
format mathematica;
#write <../Numerical_Evaluation_1loop/gg1l_results_AG/He'i'_LR.m> "HeLR['i'] = ( %E ) ; " [He'i'] ;
#write <../Numerical_Evaluation_1loop/gg1l_results_AG/He'i'_LR.m> "" ;
format normal;
#enddo
print+s;
    .end
