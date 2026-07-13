#-
#:maxtermsize 1000000
#:workspace 1G

* 1) analyze the diagrams
* 2) do color algebra
* 3) substitute Feynman Rules and prepare contractions

**************************************************************************************************************** Declare symbols for names
autodeclare s d, x, y, PL, NP, PLA, PLB, PLC, NPA, NPB, Sec, n, box0, [d-4], [d-3], [d-2]                                                          ; *********AG: added extra topology for 1 loop qq -> VV
**************************************************************************************************************** Declare kinematic variables
autodeclare s s                                                                                                                                    ;
s t, u, mt, mw2 ,m32, m42, pi                                                                                                                      ;
**************************************************************************************************************** Declare couplings colour factors and other parameters for the amplitude
s vq, aq, g, Nc, NA, CF, d, Nf, Nh, NfV, NhV, NfA, NfV2, NhV2, NfA2, NffcV, NffcA, NffcV2, NffcA2, NffcVA,[Nc^2-1], [Nc-1], [gh], I, XXX                  ;
**************************************************************************************************************** Declare some extra symbols for numbers
s n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15                                                                                 ;
**************************************************************************************************************** Declare indices
autodeclare i mu=d, nu=d, ro=d, sigma=d, a=NA, b=NA, i=Nc, j=Nc , J=Nc                                                                             ;
**************************************************************************************************************** Declare vectors for momenta and polarisations
autodeclare v k, mk, p, mp, q, eps, va                                                                                                             ;
*************************************************************************************************************** Declare extrasymbol for gamma_5
s a, x                                                                                                                                             ;

**************************************************************************************************************** Declare functions for particle names
cf Quark1, Quark2, Quark3, Qbar1, Qbar2, Qbar3, gluon, ghost, Ghost, Vboson, Vboson1, Vboson2                                                      ;
**************************************************************************************************************** Declare functions for propagators and vertices
cf prop, propmt, vrtx, Vqqbv, Vqqba, Vqqbg, Vttbv, Vttba, Vttbg, Vtqbv, Vtqba, Vqtba ,Vqtbv, Vggg, Vgggg, Vgggg1, Vgggg2, Vgggg3, VghGgh, VHpb, Vvvv;
**************************************************************************************************************** Declare functions for gamma matrices and spinor bilinears
cf gaQ, Bi, BiVU, BiUV                                                                                                                             ;
**************************************************************************************************************** Declare functions for polarisations:
cf pol, U, V, Ub, Vb                                                                                                                               ;
**************************************************************************************************************** Declare polarisation vectors as functions for individual particles
cf Eps1,Eps2,Eps3,Eps4,Eps5                                                                                                                        ;
**************************************************************************************************************** Declare colour objects
cf T,f,Dd,Tr,delC                                                                                                                                  ;
**************************************************************************************************************** Declare loops
cf LoopF,LoopGh                                                                                                                                    ;
**************************************************************************************************************** Declare functions for factors
cf Num, Den, coeff, Ratio                                                                                                                          ;
**************************************************************************************************************** Declare functions for write helicity form factors
cf FtoHe                                                                                                                                           ;
**************************************************************************************************************** Declare label type objects
cf Label, INT, DiaMatch, Shift, Sector, Save, ACC                                                                                                  ;

**************************************************************************************************************** Declare functions for spinor helicity calculus
cf spa, sps, spaInv, spsInv, spas, spsa, spaa, spss                                                                                                ;
**************************************************************************************************************** Declare extra symbol for flavour changing vertex
cf fc                                                                                                                                              ;

**************************************************************************************************************** Declare objects for projectors and form factors
autodeclare s c, F, He                                                                                                                             ;

**************************************************************************************************************** Declare propagators
s D1,...,D15                                                                                                                                       ;
**************************************************************************************************************** Declare symbols for one loop reduction to masters
cf INT1L, M                                                                                                                                        ;
s P12, P13, P13x12, root, heps                                                                                                                     ;
**************************************************************************************************************** Declare propagator symbols

****************************************************************************************************************
#include ../Procedures/Procedures.proc
#include ../Procedures/Procedures_2.proc

off statistics ;
.global
#message Importing diagrams from QGraf

#do i = 1,8

*#message ******* ----- ~~~~~ process diagram 'i' ~~~~~ ----- *******

g diag'i'  = d'i'*DiaMatch('i')                        ; ****AG: for one loop: map the integral that are defined up to a shift of the loop momenta

#enddo


#include ggVV.1loop

**
#message Plugg in momentum conservation

multiply replace_(p4,-p1-p2-p3)                       ;
.sort

** Feynman Rules
#message Start Feynman rules

*Set negative external particle indices to their values plus one hundred
argument pol, vrtx                 ;
id gluon(-3,p1?)  = gluon(97,p1) ;
id gluon(-1,p1?)   = gluon(99,p1)  ;
id Vboson1(-5,p1?)  = Vboson1(95,p1) ;
id Vboson2(-7,p1?)  = Vboson2(93,p1) ;
endargument                        ;
.sort
** Propagators
#message Doing propagators

#do i=1,20
#do j=1,20
* we assume Feynman gauge
id prop(gluon('i','j',p1?))  = - i_*d_(a'i',a'j')*d_(mu'i',mu'j')*prop(p1)                  ;
id prop(ghost('i','j',p1?))  = + i_*ghost('i','j')*d_(a'i',a'j')*prop(p1)                   ;
#enddo
#enddo
id prop(Quark1(n1?,n2?,p1?)) = i_*gaQ(n1,n2,p1)*delC(n1,n2)*prop(p1)                        ;
id prop(Quark2(n1?,n2?,p1?)) = i_*(gaQ(n1,n2,p1) + mt*gaQ(n1,n2,1) )*delC(n1,n2)*propmt(p1) ;
.sort

** Relabelling the indices of all particles in vertices and polarisations
argument vrtx,pol                            ;
#do j1=99,90,-1
id Vboson1('j1',p1?) = Vboson1(mu'j1',p1)    ;
id Vboson2('j1',p1?) = Vboson2(mu'j1',p1)    ;
id gluon('j1',p1?)  = gluon(mu'j1',a'j1',p1) ;
id ghost('j1',p1?)  = ghost('j1',a'j1',p1)   ;
id Ghost('j1',p1?)  = Ghost('j1',a'j1',p1)   ;
#enddo
#do j1=1,20
id Vboson1('j1',p1?) = Vboson1(mu'j1',p1)    ;
id Vboson2('j1',p1?) = Vboson2(mu'j1',p1)    ;
id gluon('j1',p1?)  = gluon(mu'j1',a'j1',p1) ;
id ghost('j1',p1?)  = ghost('j1',a'j1',p1)   ;
id Ghost('j1',p1?)  = Ghost('j1',a'j1',p1)   ;
#enddo
endargument                                  ;
.sort

** Vertices
#message Doing vertices (including gamma5)

id vrtx(Qbar1(n1?,p1?),Vboson2(mu1?,p3?),Quark2(n2?,p2?))
    = (vq*Vtqbv(n1,n2,mu1) + aq*Vtqba(n1,n2,mu1))*( -i_*delC(n1,n2) ) ;

id vrtx(Qbar2(n1?,p1?),Vboson1(mu1?,p3?),Quark1(n2?,p2?))
    = (vq*Vqtbv(n1,n2,mu1) + aq*Vqtba(n1,n2,mu1))*(-i_*delC(n1,n2) ) ;

id vrtx(Qbar1(n1?,p1?),gluon(mu1?,a1?,p3?),Quark1(n2?,p2?))
    = Vqqbg(n1,n2,mu1)*( -i_*g*T(a1,n1,n2) )                          ;

id vrtx(Qbar2(n1?,p1?),gluon(mu1?,a1?,p3?),Quark2(n2?,p2?))
    = Vttbg(n1,n2,mu1)*( -i_*g*T(a1,n1,n2) )                          ;

id vrtx(Ghost(n1?,a1?,p1?),gluon(mu1?,a2?,p2?),ghost(n3?,a3?,p3?))
    = VghGgh(n1,n3,p1,mu1)*( + g*f(a2,a1,a3) )                        ;

id vrtx(gluon(mu1?,a1?,p1?),gluon(mu2?,a2?,p2?),gluon(mu3?,a3?,p3?))
    = Vggg(a1,a2,a3,mu1,mu2,mu3,p1,p2,p3)                             ;

id vrtx(gluon(mu1?,a1?,p1?),gluon(mu2?,a2?,p2?),gluon(mu3?,a3?,p3?),gluon(mu4?,a4?,p4?))
    = Vgggg(a1,a2,a3,a4,mu1,mu2,mu3,mu4)                              ;

** External gluons Vboson

id pol(gluon(mu1?,a1?,p1?))      = gluon(mu1,a1,p1);
id pol(Vboson1(mu1?,p1?))        = Eps3(mu1) ;
id pol(Vboson2(mu1?,p1?))        = Eps4(mu1) ;
id gluon(mu?,a1?,p1)             = Eps1(mu)  ;
id gluon(mu?,a1?,p2)             = Eps2(mu)  ;
.sort

** AG: check that all polarisation vectors are substituted
id pol(?x) = 1/(1-1) ;
.sort

********* THIS PART MAKES THE FERMION LOOPS, CAN BE IMPROVED BY USING DOLLAR VARIABLES
********* TAKING INSPIRATION FROM THE GHOST LOOP PART
********* RECALL THAT YOU MUST CHANGE TO ACCOMODATE THE AXIAL PART


** AG: My implementation
#message Build closed fermion loops

*************************************************************************** AG: 1 loop.

*********AG: For flavour changing currents. Nffc is the number of flavour in the flavour changing vertex

*id once  Vtqbv(n1?,n2?,mu?)*gaQ(n2?,n3?,p1?)    = NffcV*LoopF(1,n1,mu,p1,n3)                                  ;
*id once  Vtqbv(n1?,n2?,mu?)*gaQ(n2?,n3?,1)      = NffcV*LoopF(1,n1,mu,n3)                                     ;
*id once  Vqtbv(n1?,n2?,mu?)*gaQ(n2?,n3?,p1?)    = NffcV*LoopF(1,n1,mu,p1,n3)                                  ;
*id once  Vqtbv(n1?,n2?,mu?)*gaQ(n2?,n3?,1)      = NffcV*LoopF(1,n1,mu,n3)                                     ;

*********AG:  Light quark loop with axial coupling

*id once  Vtqba(n1?,n2?,mu?)*gaQ(n2?,n3?,p1?)    = -1/6*NffcA*LoopF(1,n1,nu1,ro1,sigma1,p1,n3)*e_(mu,nu1,ro1,sigma1) ;
*id once  Vtqba(n1?,n2?,mu?)*gaQ(n2?,n3?,1)      = -1/6*NffcA*LoopF(1,n1,nu2,ro2,sigma2,n3)*e_(mu,nu2,ro2,sigma2)  ;

*********AG:  Top quark loop with axial coupling

*id once  Vqtba(n1?,n2?,mu?)*gaQ(n2?,n3?,p1?)    = -1/6*NffcA*LoopF(1,n1,nu3,ro3,sigma3,p1,n3)*e_(mu,nu3,ro3,sigma3) ;
*id once  Vqtba(n1?,n2?,mu?)*gaQ(n2?,n3?,1)      = -1/6*NffcA*LoopF(1,n1,nu4,ro4,sigma4,n3)*e_(mu,nu4,ro4,sigma4)  ;

*repeat                                                                  ;

*id LoopF(1,n1?,?p1,n2?)*gaQ(n2?,n3?,p2?)   = LoopF(1,n1,?p1,p2,n3)      ;
*id LoopF(1,n1?,?p1,n2?)*gaQ(n2?,n3?,1)     = LoopF(1,n1,?p1,n3)         ;

*id LoopF(1,n1?,?p1,n2?)*Vqqbg(n2?,n3?,mu?) = LoopF(1,n1,?p1,mu,n3)      ;
*id LoopF(1,n1?,?p1,n2?)*Vttbg(n2?,n3?,mu?) = LoopF(1,n1,?p1,mu,n3)      ;

*id LoopF(1,n1?,?p1,n2?)*Vqtbv(n2?,n3?,mu?) = NffcV*LoopF(1,n1,?p1,mu,n3);
*id LoopF(1,n1?,?p1,n2?)*Vtqbv(n2?,n3?,mu?) = NffcV*LoopF(1,n1,?p1,mu,n3);

*id LoopF(1,n1?,?p1,n2?)*Vtqba(n2?,n3?,mu?) = -1/6*NffcA*LoopF(1,n1,?p1,nu5,ro5, sigma5,n3)*e_(mu,nu5,ro5,sigma5);
*id LoopF(1,n1?,?p1,n2?)*Vqtba(n2?,n3?,mu?) = -1/6*NffcA*LoopF(1,n1,?p1,nu6,ro6, sigma6,n3)*e_(mu,nu6,ro6,sigma6);

*id LoopF(1,n1?,?p1,n2?)*Vtqba(n2?,n3?,mu?) = -1/6*NffcA*LoopF(1,n1,?p1,nu7,ro7, sigma7,n3)*e_(mu,nu7,ro7,sigma7);
*id LoopF(1,n1?,?p1,n2?)*Vqtba(n2?,n3?,mu?) = -1/6*NffcA*LoopF(1,n1,?p1,nu8,ro8, sigma8,n3)*e_(mu,nu8,ro8,sigma8);

*id LoopF(1,n1?,?p1,n2?)*LoopF(1,n2?,?p3,n3?) = LoopF(1,n1,?p1,?p3,n3)   ;
*endrepeat                                                               ;

*.sort



*AG: Leonardo implementation:


id once Vqtbv(n1?,n2?,mu?)*gaQ(n2?,n3?,p1?) = LoopF(1,n1,mu,p1,n3) ;
id once Vqtbv(n1?,n2?,mu?)*gaQ(n2?,n3?,1)   = LoopF(1,n1,mu,n3) ;
repeat ;
id LoopF(1,n1?,?p1,n2?)*gaQ(n2?,n3?,p2?)   = LoopF(1,n1,?p1,p2,n3) ;
id LoopF(1,n1?,?p1,n2?)*gaQ(n2?,n3?,1)     = LoopF(1,n1,?p1,n3) ;
id LoopF(1,n1?,?p1,n2?)*Vqqbg(n2?,n3?,mu?) = LoopF(1,n1,?p1,mu,n3) ;
id LoopF(1,n1?,?p1,n2?)*Vttbg(n2?,n3?,mu?) = LoopF(1,n1,?p1,mu,n3) ;
* if both vectors on the same loop
id LoopF(1,n1?,?p1,n2?)*Vtqbv(n2?,n3?,mu?) = LoopF(1,n1,?p1,mu,n3) ;
id LoopF(1,n1?,?p1,n2?)*Vtqba(n2?,n3?,mu?) = -1/6*LoopF(1,n1,?p1,nu1,ro1,sigma1,n3)*e_(mu,nu1,ro1,sigma1) ;
endrepeat ;
*for axial-axial contributions
id once Vqtba(n1?,n2?,mu?)*gaQ(n2?,n3?,p1?) = -1/6*LoopF(1,n1,nu2,ro2,sigma2,p1,n3)* e_(mu,nu2,ro2,sigma2);
repeat ;
id LoopF(1,n1?,?p1,n2?)*gaQ(n2?,n3?,p2?)    = LoopF(1,n1,?p1,p2,n3) ;
id LoopF(1,n1?,?p1,n2?)*gaQ(n2?,n3?,1)      = LoopF(1,n1,?p1,n3) ;
id LoopF(1,n1?,?p1,n2?)*Vqqbg(n2?,n3?,mu?)  = LoopF(1,n1,?p1,mu,n3) ;
id LoopF(1,n1?,?p1,n2?)*Vttbg(n2?,n3?,mu?)  = LoopF(1,n1,?p1,mu,n3) ;
* if both photons on the same loop
id LoopF(1,n1?,?p1,n2?)*Vtqbv(n2?,n3?,mu?) = LoopF(1,n1,?p1,mu,n3) ;
id LoopF(1,n1?,?p1,n2?)*Vtqba(n2?,n3?,mu?) = -1/6*LoopF(1,n1,?p1,nu3,ro3,sigma3,n3)*e_(mu,nu3,ro3,sigma3) ;
endrepeat ;
id LoopF(n1?,n2?,?x,n2?) = LoopF(n1,?x) ;
.sort



* double check that everything closes to itself

id LoopF(n1?,n2?,?x,n2?)                       = LoopF(n1,?x)           ;
.sort

*check that we caught all fermion loops
id gaQ(?x) = 1/(1-1) ;
.sort

*************************************************************************** Color Algebra
#message Color algebra 1/2: contract deltas and Ts

b delC,T                            ;
.sort
keep brackets                       ;
#do i=1,100
#do j=1,100
id delC('i','j') = d_(j'i',j'j')    ;
id T(a1?,'i','j') = T(a1,j'i',j'j') ;
#enddo
#enddo
.sort

#message Color algebra 2/2: Build Color Traces

antisymmetrize f;
.sort

#call coloralgebra


**       Contract with polarisation vectors
#message Contract with polarisation vectors
#do i=1,4
id Eps'i'(mu?) = eps'i'(mu) ;
#enddo
.sort

#call kinematics
#call transversality
#call gaugechoice
.sort


* Print intermediate results

b  vq, aq, d_, i_, DiaMatch, NffcV, NffcA ;
#do i = 1,8
  #write <Diagrams/diagram'i'.inc> "id diag'i' = %E;", diag'i'
#enddo
.sort
print+s;

#message Done!
.end
