#-

#:maxtermsize 1000000
#:workspace 1G

* analyze the diagrams
* 1) take diagrams from prog_1
* 2) apply conjugate tensors to the amplitude
* 3) substitute reduze notations needed for the master integrals reduction

**************************************************************************************************************** Declare symbols for names
autodeclare s d, x, y, PL, NP, PLA, PLB, PLC, NPA, NPB, Sec, n, box0, [d-4], [d-3], [d-2]                               ; *********AG: added extra topology for 1 loop qq -> VV
**************************************************************************************************************** Declare kinematic variables
autodeclare s s                                                                                                         ;
s t, u, mt, mw2 ,m32, m42, pi                                                                                               ;
**************************************************************************************************************** Declare couplings colour factors and other parameters for the amplitude
s vq, aq, g, Nc, NA, CF, d, Nf, Nh, NfV, NhV, NfA, NfV2, NhV2, NffcV, NffcA, NfA2[], [Nc^2-1], [Nc-1], [gh], I, XXX                   ;
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
cf LoopF,LoopGh                                                                                                         ;
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
***** AG: need to declare also denominators for the three loop case, otherwise the procedures cannot be called
*one loop

*massless propagators
#do j=1,3
s
Dk'j'                                ,
#do i1=1,3
Dk'j'p'i1',Dk'j'm'i1'                ,
#do i2='i1',3
Dk'j'p'i1''i2',Dk'j'm'i1''i2'        ,
#do i3='i2',3
Dk'j'p'i1''i2''i3',Dk'j'm'i1''i2''i3',
#enddo
#enddo
#enddo
#enddo

*massive propagators
#do j=1,3
s
Dmk'j',
#do i1=1,3
Dmk'j'p'i1',Dmk'j'm'i1',
#do i2='i1',3
Dmk'j'p'i1''i2',Dmk'j'm'i1''i2',
#do i3='i2',3
Dmk'j'p'i1''i2''i3',Dmk'j'm'i1''i2''i3',
#enddo
#enddo
#enddo
;
#enddo


*two loop

*massless propagators
#do j1=1,2
#do j2='j1',3
s
Dk`j1'`j2',

#do i1=1,3
Dk`j1'`j2'p'i1',Dk`j1'`j2'm'i1',

#do i2='i1',3
Dk`j1'`j2'p'i1''i2',Dk`j1'`j2'm'i1''i2'

#do i3='i2',3
Dk`j1'`j2'p'i1''i2''i3',Dk`j1'`j2'm'i1''i2''i3',

#enddo
#enddo
#enddo
;
#enddo
#enddo

*massive propagators
#do j1=1,2
#do j2='j1',3
s
Dmk`j1'`j2',

#do i1=1,3
Dmk`j1'`j2'p'i1',Dmk`j1'`j2'm'i1',

#do i2='i1',3
Dmk`j1'`j2'p'i1''i2',Dmk`j1'`j2'm'i1''i2'

#do i3='i2',3
Dmk`j1'`j2'p'i1''i2''i3',Dmk`j1'`j2'm'i1''i2''i3',

#enddo
#enddo
#enddo
;
#enddo
#enddo



*three loop

*massless propagators
s
Dk123,

#do i1=1,3
Dk123p'i1',Dk123m'i1',

#do i2='i1',3
Dk123p'i1''i2',Dk123m'i1''i2',

#do i3='i2',3
Dk123p'i1''i2''i3',Dk123m'i1''i2''i3',

#enddo
#enddo
#enddo
;


*massive propagators
s
Dmk123,

#do i1=1,3
Dmk123p'i1',Dmk123m'i1',

#do i2='i1',3
Dmk123p'i1''i2',Dmk123m'i1''i2',

#do i3='i2',3
Dmk123p'i1''i2''i3',Dmk123m'i1''i2''i3',

#enddo
#enddo
#enddo
;


set color: Nc,Nf,NfV,NfV2,CF,NA ;

set props:

#do j=1,3

Dk'j',
#do i1=1,3
Dk'j'p'i1',Dk'j'm'i1',
#do i2='i1',3
Dk'j'p'i1''i2',Dk'j'm'i1''i2'
#do i3='i2',3
Dk'j'p'i1''i2''i3',Dk'j'm'i1''i2''i3',
#enddo
#enddo
#enddo

#enddo

#do j=1,3

Dmk'j',
#do i1=1,3
Dmk'j'p'i1',Dmk'j'm'i1',
#do i2='i1',3
Dmk'j'p'i1''i2',Dmk'j'm'i1''i2'
#do i3='i2',3
Dmk'j'p'i1''i2''i3',Dmk'j'm'i1''i2''i3',
#enddo
#enddo
#enddo

#enddo

Dk12,
#do i1=1,3
Dk12p'i1',Dk12m'i1',
#do i2='i1',3
Dk12p'i1''i2',Dk12m'i1''i2'
#do i3='i2',3
Dk12p'i1''i2''i3',Dk12m'i1''i2''i3',
#enddo
#enddo
#enddo

Dmk12,
#do i1=1,3
Dmk12p'i1',Dmk12m'i1',
#do i2='i1',3
Dmk12p'i1''i2',Dmk12m'i1''i2'
#do i3='i2',3
Dmk12p'i1''i2''i3',Dmk12m'i1''i2''i3',
#enddo
#enddo
#enddo

Dk13,
#do i1=1,3
Dk13p'i1',Dk13m'i1',
#do i2='i1',3
Dk13p'i1''i2',Dk13m'i1''i2'
#do i3='i2',3
Dk13p'i1''i2''i3',Dk13m'i1''i2''i3',
#enddo
#enddo
#enddo

Dmk13,
#do i1=1,3
Dmk13p'i1',Dmk13m'i1',
#do i2='i1',3
Dmk13p'i1''i2',Dmk13m'i1''i2'
#do i3='i2',3
Dmk13p'i1''i2''i3',Dmk13m'i1''i2''i3',
#enddo
#enddo
#enddo

Dk23,
#do i1=1,3
Dk23p'i1',Dk23m'i1',
#do i2='i1',3
Dk23p'i1''i2',Dk23m'i1''i2'
#do i3='i2',3
Dk23p'i1''i2''i3',Dk23m'i1''i2''i3',
#enddo
#enddo
#enddo

Dmk23,
#do i1=1,3
Dmk23p'i1',Dmk23m'i1',
#do i2='i1',3
Dmk23p'i1''i2',Dmk23m'i1''i2'
#do i3='i2',3
Dmk23p'i1''i2''i3',Dmk23m'i1''i2''i3',
#enddo
#enddo
#enddo

Dk123,
#do i1=1,3
Dk123p'i1',Dk123m'i1',
#do i2='i1',3
Dk123p'i1''i2',Dk123m'i1''i2'
#do i3='i2',3
Dk123p'i1''i2''i3',Dk123m'i1''i2''i3',
#enddo
#enddo
#enddo

Dmk123,
#do i1=1,3
Dmk123p'i1',Dmk123m'i1',
#do i2='i1',3
Dmk123p'i1''i2',Dmk123m'i1''i2'
#do i3='i2',3
Dmk123p'i1''i2''i3',Dmk123m'i1''i2''i3',
#enddo
#enddo
#enddo
;

set propsnum: D1,...,D15 ;
****************************************************************************************************************

#include ../Procedures/Procedures.proc
#include ../Procedures/Procedures_2.proc

off statistics ;
.global

********AG: Start analysis of diagrams

#message importing diagrams from prog_1

g amp =
#do i1=1,8
+ diag'i1'
#enddo
;
#do i=1,8
#include Diagrams/diagram'i'.inc
#enddo;

#message Importing shifts and sector mappings from Reduze
multiply replace_(k1,q1)                        ;
.sort

*************************************************************************** AG: for 1loop

#message include shifts from jobs_2 in reduze
#include ../reduze_results/mytmp/1loop_reduction/ggWW.1loop_shift.inc

id Shift(k1,p1?,[])            = replace_(q1,p1)       ;
*id Shift(k1,p1?,[],k2,p2?,[]) = replace_(q1,p1,q2,p2) ;
*id Shift(k2,p2?,[],k1,p1?,[]) = replace_(q2,p2,q1,p1) ;
.sort

#message Cleaning up propagators
factarg prop, propmt                                                                                                                    ;
.sort

***************************************************************************
*************************************************************************** AG: I leave prop just the propagators with loop momenta
repeat id prop(k1?,-1,1)       = prop(k1)                                                                                               ;
repeat id prop(-1,p?)          = prop(p)                                                                                                ;
repeat id propmt(k1?,-1,1)     = propmt(k1)                                                                                             ; *For mt propagators
repeat id propmt(-1,p?)        = propmt(p)                                                                                              ;
.sort

#do i1=1,3
id prop(p'i1')                  = Den( ( p'i1'(mu) ) * ( p'i1'(mu) ) )                                                                   ;
id propmt(p'i1')                = Den( ( p'i1'(mu) ) * ( p'i1'(mu) ) - mt^2)                                                             ; *For mt propagators
#do i2='i1'+1,3
id prop(p'i1'+p'i2')            = Den( ( p'i1'(mu)+p'i2'(mu) ) * ( p'i1'(mu)+p'i2'(mu) ) )                                               ;
id propmt(p'i1'+p'i2')          = Den( ( p'i1'(mu)+p'i2'(mu) ) * ( p'i1'(mu)+p'i2'(mu) ) - mt^2)                                         ; *For mt propagators
#do i3='i2'+1,3
id prop(p'i1'+p'i2'+p'i3')      = Den( ( p'i1'(mu)+p'i2'(mu)+p'i3'(mu) ) * ( p'i1'(mu)+p'i2'(mu)+p'i3'(mu) ) )                           ;
id propmt(p'i1'+p'i2'+p'i3')    = Den( ( p'i1'(mu)+p'i2'(mu)+p'i3'(mu) ) * ( p'i1'(mu)+p'i2'(mu)+p'i3'(mu) ) -mt^2)                      ;
#enddo
#enddo
#enddo
.sort
*************************************************************************** AG: commented scalar product with m3=m4
argument Den             ;
id p4 = -p1-p2-p3        ;
id p1.p1 = 0             ;
id p2.p2 = 0             ;
id p3.p3 = m32           ;
id p1.p2 = s/2           ;
id p1.p3 = (-s-u+m42)/2  ;
id p2.p3 = (u-m32)/2     ;
endargument              ;
.sort
#call FactArgDen
.sort
********* Substitute propagators
#message Identifying propagators and renaming momenta for traces

id k.va = 0                 ;

argument prop, propmt, LoopF;
id k    = k1                ;
id p1?.k = p1.k1            ;
endargument                 ;
id p6?.va = e_(p1,p2,p3,p6) ;

#call identifyprop
#call namemomenta

*check everything got subsistuted
id prop(?x)   = 1/(1-1) ;
id propmt(?x) = 1/(1-1) ;
.sort
id p6?.va = e_(p1,p2,p3,p6) ;

#message declaration denominators


********* Store compact intermediate result

b vq, aq, g, d_, i_, Sector ;

.store

************ AG: compute here Td_j * Amp. Then we need to combine all the contributions weighted with the correspoinding coefficient


#include ../Tensor_decomposition/tensors.inc
.sort
.store
#message apply conj. tensors to the amplitude.
#do i1=1,36
g [AmpT'i1'd]       = [T'i1'd]  * amp ;
#enddo
.sort
#message rename momenta inside fermion loop to perform traces
#call namemomentaLoopF
.sort

id LoopF(1,?x) = g_(1,?x);
.sort
*Check that all Fermion loops got converted into trace of gamma matrices
id LoopF(?x) = 1/(1-1) ;
.sort
tracen 1;
.sort

#message First step of opening numerator momenta

*repeat ;
*id mp12 = -p12 ;
*id mp13 = -p13 ;
*id mp23 = -p23 ;
*id mp123 = -p123 ;
*id mk1 = -k1;

*id p12 = p1+p2 ;
*id p13 = p1+p3 ;
*id p23 = p2+p3 ;
*id p123 = p1+p2+p3 ;
*endrepeat ;
*.sort

#do j=1,1
repeat id q?.mk'j'= -q.k'j'  ;
#do i1=1,3
repeat id q?.k'j'p'i1'= q.k'j'+q.p'i1'  ;
repeat id q?.k'j'm'i1'= q.k'j'-q.p'i1'  ;
repeat id q?.mk'j'p'i1'= -q.k'j'-q.p'i1';
repeat id q?.mk'j'm'i1'= -q.k'j'+q.p'i1';
#do i2='i1',3
repeat id q?.k'j'p'i1''i2' = q.k'j'+q.p'i1'+q.p'i2' ;
repeat id q?.k'j'm'i1''i2' = q.k'j'-q.p'i1'-q.p'i2' ;
repeat id q?.mk'j'p'i1''i2'= -q.k'j'-q.p'i1'-q.p'i2';
repeat id q?.mk'j'm'i1''i2'= -q.k'j'+q.p'i1'+q.p'i2';
#do i3='i2',3
repeat id q?.k'j'p'i1''i2''i3'= q.k'j'+q.p'i1'+q.p'i2'+q.p'i3'   ;
repeat id q?.k'j'm'i1''i2''i3'= q.k'j'-q.p'i1'-q.p'i2'-q.p'i3'   ;
repeat id q?.mk'j'p'i1''i2''i3'=-q.k'j'-q.p'i1'-q.p'i2'-q.p'i3'  ;
repeat id q?.mk'j'm'i1''i2''i3'= -q.k'j'+q.p'i1'+q.p'i2'+q.p'i3' ;
#enddo
#enddo
#enddo
#enddo

.sort

#call kinematics
*#call transversality
#call gaugechoice
.sort

#call namemomenta    ;
#call allsimplnosort ;
.sort
#message Deal with LEVI-CIVITA
*#call Levi_Civita

id e_(va,p4?,p5?,p6?) = e_(p1,p2,p3,ro20) * e_(ro20, p4,p5,p6);
Contract 0;
    .sort

#call kinematics
#call transversality
#call gaugechoice
.sort

#message Computing polarisation sums

#message First polarisation sum

id eps1.p6?*eps1c.p7? = (-p6.p7 +2*(p1.p6*p2.p7+p1.p7*p2.p6)/(s) )  ;
id eps1.p6?*eps1c.p6? = (-p6.p6 +4*(p1.p6*p2.p6)/(s) )              ;
id eps1.eps1c = 2-d                                                 ;

id e_(eps1,p6?,p7?,p8?)*eps1c.p9? = -e_(p9,p6,p7,p8) + 2*(e_(p1,p6,p7,p8)*p2.p9+e_(p2,p6,p7,p8)*p1.p9)/(s) ;
id e_(eps1c,p6?,p7?,p8?)*eps1.p9? = -e_(p9,p6,p7,p8) + 2*(e_(p1,p6,p7,p8)*p2.p9+e_(p2,p6,p7,p8)*p1.p9)/(s) ;
id e_(eps1,eps1c,p6?,p7?)         = 0                                                                      ;
id eps1.va                        = e_(p1,p2,p3,eps1)                                                      ;
id eps1c.va                       = e_(p1,p2,p3,eps1c)                                                     ;

#call allsimplnosort
.sort
#call kinematics
*#call transversality
#call gaugechoice

#message Second polarization sum
id eps2.p6?*eps2c.p7? = (-p6.p7 +2*(p1.p6*p2.p7+p1.p7*p2.p6)/(s) )  ;
id eps2.p6?*eps2c.p6? = (-p6.p6 +4*(p1.p6*p2.p6)/(s) )              ;
id eps2.eps2c = 2-d                                                 ;

id e_(eps2,p6?,p7?,p8?)*eps2c.p9? = -e_(p9,p6,p7,p8) + 2*(e_(p1,p6,p7,p8)*p2.p9+e_(p2,p6,p7,p8)*p1.p9)/(s) ;
id e_(eps2c,p6?,p7?,p8?)*eps2.p9? = -e_(p9,p6,p7,p8) + 2*(e_(p1,p6,p7,p8)*p2.p9+e_(p2,p6,p7,p8)*p1.p9)/(s) ;
id e_(eps2,eps2c,p6?,p7?)         = 0                                                                      ;
id eps2.va                        = e_(p1,p2,p3,eps2)                                                      ;
id eps2c.va                       = e_(p1,p2,p3,eps2c)                                                     ;

#call allsimplnosort
.sort
#call kinematics
*#call transversality
#call gaugechoice

#message Third polarization sum
id eps3.p6?*eps3c.p7?             = (-p6.p7 + (p3.p6*p3.p7)/(m32) )                ;
id eps3.p6?*eps3c.p6?             = (-p6.p6 + (p3.p6*p3.p6)/(m32) )                ;
id eps3.eps3c                     = 1-d                                            ;

id e_(eps3,p6?,p7?,p8?)*eps3c.p9? = -e_(p9,p6,p7,p8) + e_(p3,p6,p7,p8)*p3.p9/(m32) ;
id e_(eps3c,p6?,p7?,p8?)*eps3.p9? = -e_(p9,p6,p7,p8) + e_(p3,p6,p7,p8)*p3.p9/(m32) ;
id e_(eps3,eps3c,p6?,p7?)         = 0                                              ;
id e_(eps3,eps3c,p6?,p7?)         = 0                                              ;
id eps3.va                        = e_(p1,p2,p3,eps3)                              ;
id eps3c.va                       = e_(p1,p2,p3,eps3c)                             ;

#call allsimplnosort
.sort
#call kinematics
*#call transversality
#call gaugechoice

#message Fourth polarization sum

id eps4.p6?*eps4c.p7?             = (-p6.p7 + (p4.p6*p4.p7)/(m42) )                ;
id eps4.p6?*eps4c.p6?             = (-p6.p6 + (p4.p6*p4.p6)/(m42) )                ;
id eps4.eps4c                     = 1-d                                            ;

id e_(eps4,p6?,p7?,p8?)*eps4c.p9? = -e_(p9,p6,p7,p8) + e_(p4,p6,p7,p8)*p4.p9/(m42) ;
id e_(eps4c,p6?,p7?,p8?)*eps4.p9? = -e_(p9,p6,p7,p8) + e_(p4,p6,p7,p8)*p4.p9/(m42) ;
id e_(eps4,eps4c,p6?,p7?)         = 0                                              ;
id eps4.va                        = e_(p1,p2,p3,eps4)                              ;
id eps4c.va                       = e_(p1,p2,p3,eps4c)                             ;


#message Deal with final levi civita tensors
*#call Levi_Civita

id e_(va, p6?,p7?,p8?) = e_(p1,p2,p3,ro21) * e_(ro21,p6,p7,p8) ;
Contract 0 ;
.sort

#call kinematics
#call transversality
#call gaugechoice
.sort

.sort
#do j=1,1
repeat id q?.mk'j'= -q.k'j'  ;
#do i1=1,3
repeat id q?.k'j'p'i1'= q.k'j'+q.p'i1'  ;
repeat id q?.k'j'm'i1'= q.k'j'-q.p'i1'  ;
repeat id q?.mk'j'p'i1'= -q.k'j'-q.p'i1';
repeat id q?.mk'j'm'i1'= -q.k'j'+q.p'i1';
#do i2='i1',3
repeat id q?.k'j'p'i1''i2' = q.k'j'+q.p'i1'+q.p'i2' ;
repeat id q?.k'j'm'i1''i2' = q.k'j'-q.p'i1'-q.p'i2' ;
repeat id q?.mk'j'p'i1''i2'= -q.k'j'-q.p'i1'-q.p'i2';
repeat id q?.mk'j'm'i1''i2'= -q.k'j'+q.p'i1'+q.p'i2';
#do i3='i2',3
repeat id q?.k'j'p'i1''i2''i3'= q.k'j'+q.p'i1'+q.p'i2'+q.p'i3'   ;
repeat id q?.k'j'm'i1''i2''i3'= q.k'j'-q.p'i1'-q.p'i2'-q.p'i3'   ;
repeat id q?.mk'j'p'i1''i2''i3'=-q.k'j'-q.p'i1'-q.p'i2'-q.p'i3'  ;
repeat id q?.mk'j'm'i1''i2''i3'= -q.k'j'+q.p'i1'+q.p'i2'+q.p'i3' ;
#enddo
#enddo
#enddo
#enddo

.sort
#call allsimplnosort
.sort
#call kinematics
*#call transversality
#call gaugechoice
    .sort

*AG: at this point only k, p_i are left inside the levi civita. Thus we can set them to 0.
id k1.va =0;
    .sort

id e_(p1?,p2?,p3?,p4?) =0;
.sort


#message Mapping propagators onto integral families
*
id Sector(PL?,?x) = Label(PL) ;

****************************************************** maps propagators onto the propagators of the integral families and their crossings
******** AG:  map each integral into his topologies
#include Topos1L.rec
**check that all propagators got mapped
id 1/Dk1?props = 1/(1-1) ;
id 1/Dmk1?props = 1/(1-1) ;
.sort
*

*#message Third step of opening numerator momenta

*invert shorthand notation and go back to sums of momenta
*#call mktok
*#call openmomenta
*#call kinematics
*#call allsimplnosort

*id k1.va           = 0;
*id e_(p1,p2,p3,k1) = 0;

#message Rewriting numerators in terms of inverse propagators
****************************************************** Write scalar products in terms of inverse propagators

#include Scalar_products/ScalProd_1L.inc

.sort

*check that all scalar products got subsituted
id k1 = 1/(1-1) ;

.sort


*AG: Arrived here
id Label(PL?)/<D1^n1?>/.../<D4^n4?> = INT(PL,n1,...,n4) ;

******* AG: simplify rational coefficients
b i_, d_, CF, Nf, Nh, NfV, NhV, NfV2, NhV2, g, vq, aq, Nc, INT, Den, Num, Ratio, d, [d-2], [d-3], [d-4], cd, NffcV, NffcA, Tr ;
.sort

Collect Num                                             ;
.sort
id Num(n1?) = Ratio(n1,1)                               ;
id Den(n1?) = Ratio(1,n1)                               ;
.sort
polyratfun Ratio                                        ;
.sort
polyratfun                                              ;
.sort
id Ratio(n1?,n2?) =  Num(n1)*Den(n2)  ;
.sort
#call FactArgNum
#call FactArgDen
.sort

#message Converting integrals to standard Reduze format
id INT(box0?,n1?,n2?,n3?,n4?) =
  INT(box0, coeff(n1) + coeff(n2) + coeff(n3) + coeff(n4),
          coeff(n1)*2^0 + coeff(n2)*2^1 + coeff(n3)*2^2 + coeff(n4)*2^3,
          (abs_(n1) + n1 )/2 + (abs_(n2) + n2 )/2 + (abs_(n3) + n3 )/2 + (abs_(n4) + n4 )/2 ,
          (abs_(n1) - n1 )/2 + (abs_(n2) - n2 )/2 + (abs_(n3) - n3 )/2 + (abs_(n4) - n4 )/2 ,
          [],
          n1,...,n4)                  ;
.sort


argument INT                          ;
id coeff(0) = 0                       ;
endargument                           ;
.sort
argument INT                          ;
id coeff(n1?) = (abs_(n1) + n1 )/2/n1 ;
endargument                           ;
.sort
b i_, d_, g, CF, vq, aq, Num, Den ;

#do i1=1,36
#write <Tensors_x_amp/AmpT'i1'd.inc> "g [AmpT'i1'd] = %E ; " [AmpT'i1'd]
#enddo
.sort
**** AG: printing integrals for jobs_4 in reduze
b INT;

print[];
.end
