#-

#:maxtermsize 1000000
#:workspace 1G


* 1) Substitute the master from reduze and simplify the tensors applied to the amplitutes

**************************************************************************************************************** Declare symbols for names
autodeclare s d, x, y, PL, NP, PLA, PLB, PLC, NPA, NPB, Sec, n, box0, [d-4], [d-3], [d-2], [d-1], PL00M, PL0M0, NPL0M          ; *********AG: added extra topology for 1 loop qq -> VV
**************************************************************************************************************** Declare kinematic variables
autodeclare s s                                                                                                         ;
s t, u, mt, mw2 ,m32, m42, pi                                                                                           ;
**************************************************************************************************************** Declare couplings colour factors and other parameters for the amplitude
s vq, aq, g, Nc, NA, CF, d, Nf, Nh, NfV, NhV, NfA, NfV2, NhV2, NfA2[], [Nc^2-1], [Nc-1], [gh], I, XXX, NffcV, NffcA     ;
**************************************************************************************************************** Declare some extra symbols for numbers
s n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15                                                      ;
**************************************************************************************************************** Declare indices
autodeclare i mu=d, nu=d, ro=d, sigma=d, a=NA, b=NA, i=Nc, j=Nc , J=Nc                                                  ;
**************************************************************************************************************** Declare vectors for momenta and polarisations
autodeclare v k, mk, p, mp, q, eps, va                                                                                  ;
*************************************************************************************************************** Declare extrasymbol for gamma_5
s a, x                                                                                                                ;

**************************************************************************************************************** Declare functions for particle names
cf Quark1, Quark2, Quark3, Qbar1, Qbar2, Qbar3, gluon, ghost, Ghost, Vboson, Vboson1, Vboson2                                          ;
**************************************************************************************************************** Declare functions for propagators and vertices
cf prop, propmt, vrtx, Vqqbv, Vqqba, Vqqbg, Vttbv, Vttba, Vttbg, Vggg, Vgggg, Vgggg1, Vgggg2, Vgggg3, VghGgh, VHpb, Vvvv;
**************************************************************************************************************** Declare functions for gamma matrices and spinor bilinears
cf gaQ, Bi, BiVU, BiUV                                                                                                  ;
**************************************************************************************************************** Declare functions for polarisations:
cf pol, U, V, Ub, Vb                                                                                                    ;
**************************************************************************************************************** Declare polarisation vectors as functions for individual particles
cf Eps1,Eps2,Eps3,Eps4,Eps5                                                                                             ;
**************************************************************************************************************** Declare colour objects
cf T,f,Dd,Tr,delC                                                                                                       ;
**************************************************************************************************************** Declare loops
cf LoopF,LoopGh, PL00M, PL00Mx12, PL0M0, PL0M0x12, NPL0M, NPL0Mx12                                                       ;
**************************************************************************************************************** Declare functions for factors
cf Num, Den, coeff, Ratio                                                                                               ;
**************************************************************************************************************** Declare functions for write helicity form factors
cf FtoHe                                                                                                                ;
**************************************************************************************************************** Declare label type objects
cf Label, INT, DiaMatch, Shift, Sector, Save, ACC                                                                       ;

**************************************************************************************************************** Declare functions for spinor helicity calculus
cf spa, sps, spaInv, spsInv, spas, spsa, spaa, spss                                                                     ;

**************************************************************************************************************** Declare objects for projectors and form factors
autodeclare s c, F, He                                                                                                  ;

**************************************************************************************************************** Declare propagators
s D1,...,D15                                                                                                            ;
**************************************************************************************************************** Declare symbols for one loop reduction to masters
cf INT1L, M                                                                                                             ;
s P12, P13, P13x12, root, heps                                                                                          ;
**************************************************************************************************************** Declare propagator symbols

#include ../Procedures/Procedures.proc
#include ../Procedures/Procedures_2.proc



off statistics ;
.global
#message including Tensor x amp

#do i=1,36
#include Tensors_x_amp/AmpT'i'd.inc
#enddo;

*#message computing form factors

*#do i1=1,36
*g [F'i1'] =
*#do i2=1,36
*+ c'i2'p'i1' * [AmpT'i2'd]
*#enddo
*;
*#enddo
*.sort

*#include ../Tensor_decomposition/projectorscoeffs.inc
    .sort
#message calling polyratfun to simplify the expressions



#message Substituting master integrals from reduze

b INT ;
.sort
keep brackets;
#include ../reduze_results/1loop/List_1loop.inc
.sort

#call FactArgDen

id Den (-1+d) = 1/[d-1];
id Den(-3+d)  = 1/[d-3];
id Den(-2+d)  = 1/[d-2];
id Den(d-4)   = 1/[d-4];

b i_, d_, CF, Nf, Nh, NfV, NhV, NfV2, NhV2, g, vq, aq, Nc, Den, Num, Ratio, NffcV, NffcA, INT, Tr, d, [d-1],[d-2], [d-3], [d-4];
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
#message change in notation to match AMFlow and printing on external file (now commented)

*AG: for Leonardo code
*id INT(x?,n1?,n2?,n3?,n4?,[],n5?,n6?,n7?,n8?)    = INT(x,n5,n6,n7,n8);
*.sort



*AG: using AMFLow notation for numerical evaluation
id INT(PL00M,n1?,n2?,n3?,n4?,[],n5?,n6?,n7?,n8?)    = PL00M(n5,n6,n7,n8);
id INT(PL00Mx12,n1?,n2?,n3?,n4?,[],n5?,n6?,n7?,n8?) = PL00Mx12(n5,n6,n7,n8);
id INT(PL0M0,n1?,n2?,n3?,n4?,[],n5?,n6?,n7?,n8?)    = PL0M0(n5,n6,n7,n8);
id INT(PL0M0x12,n1?,n2?,n3?,n4?,[],n5?,n6?,n7?,n8?) = PL0M0x12(n5,n6,n7,n8);
id INT(NPL0M,n1?,n2?,n3?,n4?,[],n5?,n6?,n7?,n8?)    = NPL0M(n5,n6,n7,n8);
id INT(NPL0Mx12,n1?,n2?,n3?,n4?,[],n5?,n6?,n7?,n8?) = NPL0Mx12(n5,n6,n7,n8);
.sort

    #do i1=1,36
    #write <Tensors_x_amp/Simplified_AmpT'i1'd.inc> "id AmpT'i1'd = %E ; " [AmpT'i1'd]
    #enddo
    .sort
    print +s;
    .end
