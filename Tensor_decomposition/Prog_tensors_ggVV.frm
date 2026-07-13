#-
#:maxtermsize 1000000
#:workspace 1G
*************************************************************************** All tensor structures for q qb --> two vector bosons


*************************************************************************** AG: declaration of the symbols
s [],d, [d-1],[d-2],[d-3],[d-4],[d-5],[d-6] ;
autodeclare s F,c,n,m                       ;
* AG: Mandelstam invariants
s s,t,u,m32,m42,mw2                         ;
* AG: spinor indices
s n1,n2                                     ;
autodeclare i mu = d, nu = d                ;
autodeclare v eps,epsc,p,q,va               ;
cf Eps,Num,Den,Ratio,Label                  ;
cf BiUV, BiVU                               ; * AG: BiUV = ubar(p2)....v(p1), BiVU = vbar(p1)....u(p2)

*************************************************************************** AG: Factorization of numrators and denominators
#Procedure FactArgDen
factarg Den ;
.sort
repeat id Den(n1?,n2?,?x) = Den(n1)*Den(n2,?x) ;
.sort
repeat id Den(n1?number_)   = 1/n1             ; * AG: if n1 is a number: Den(n1)->1/n1
repeat id Den(n1?symbol_)   = 1/n1             ; * AG: if n1 is a symbol: Den(n1)->1/n1
repeat id Den(1/n1?symbol_) = n1               ; * AG: if n1 is a symbol: Den(1/n1)->n1
.sort
#EndProcedure


#Procedure FactArgNum
factarg Num ;
.sort
repeat id Num(n1?,n2?,?x) = Num(n1)*Num(n2,?x) ;
.sort
repeat id Num(n1?number_) = n1                 ;
repeat id Num(n1?symbol_) = n1                 ;
repeat id Num(1/n1?symbol_) = 1/n1             ;
.sort
#EndProcedure


*************************************************************************** AG: Kinematics
#Procedure allsimplnosort
id p1.p1  = 0                                       ;
id p2.p2  = 0                                       ;
id p3.p3  = m32                                     ;
id p1.p2. = s/2                                     ;
id p1.p3  = (-s-u+m42)/2                            ;
id p2.p3  = (u-m32)/2                               ;
id p1.va  = 0                                       ;
id p2.va  = 0                                       ;
id p3.va  = 0                                       ;
id va.va  = (d-3)/4*(s*(-s-u+m42)*(u-m32)-s*s*m32)  ;
#EndProcedure


#Procedure kinematics
id p4     = -p1-p2-p3                               ;
.sort
id p1.p1  = 0                                       ;
id p2.p2  = 0                                       ;
id p3.p3  = m32                                     ;
id p1.p2 = s/2                                      ;
id p1.p3  = (-s-u+m42)/2                            ;
id p2.p3  = (u-m32)/2                               ;
id p1.va  = 0                                       ;
id p2.va  = 0                                       ;
id p3.va  = 0                                       ;
*id va.va = (d-3)/4*(s*(-s-u+mw2)*(u-mw2)-s*s*mw2)  ;
id va.va  = (d-3)/4*(s*(-s-u+m42)*(u-m32)-s*s*m32)  ;
.sort
#EndProcedure

*************************************************************************** Gauge condition

#Procedure gaugechoice
id eps1.p1  = 0                   ;
id eps1c.p1 = 0                   ;
id eps1.p2  = 0                   ;
id eps1c.p2 = 0                   ;
id eps2.p1  = 0                   ;
id eps2c.p1 = 0                   ;
id eps2.p2  = 0                   ;
id eps2c.p2 = 0                   ;
id eps3.p3  = 0                   ;
id eps4.p3  = -eps4.p1 -eps4.p2   ;
id eps3c.p3 = 0                   ;
id eps4c.p3 = -eps4c.p1 -eps4c.p2 ;
.sort
#EndProcedure


off statistics                    ;
.global

*************************************************************************** Define the tensors and their hermitian conjugates

#message The tensors and their hermitian conjugates are

*parity even
g [T1] = eps1.p3*p3.eps2*eps4.p1*p1.eps3 ;
g [T2] = eps1.p3*p3.eps2*eps4.p1*eps3.p2 ;
g [T3] = eps1.p3*p3.eps2*eps4.p2*p1.eps3 ;
g [T4] = eps1.p3*p3.eps2*eps4.p2*eps3.p2 ;
g [T5] = eps1.p3*p3.eps2*eps4.eps3 ;
g [T6] = eps1.eps2*eps4.p1*p1.eps3 ;
g [T7] = eps1.eps2*eps4.p1*eps3.p2 ;
g [T8] = eps1.eps2*eps4.p2*p1.eps3 ;
g [T9] = eps1.eps2*eps4.p2*eps3.p2 ;
g [T10] = eps1.eps4*p3.eps2*p1.eps3 ;
g [T11] = eps1.eps4*p3.eps2*eps3.p2 ;
g [T12] = eps1.eps3*p3.eps2*eps4.p1 ;
g [T13] = eps1.eps3*p3.eps2*eps4.p2 ;
g [T14] = eps1.p3*eps2.eps4*p1.eps3 ;
g [T15] = eps1.p3*eps2.eps4*eps3.p2 ;
g [T16] = eps1.p3*eps2.eps3*eps4.p1 ;
g [T17] = eps1.p3*eps2.eps3*eps4.p2 ;
g [T18] = eps1.eps2*eps4.eps3 + eps1.eps4*eps2.eps3 + eps1.eps3*eps2.eps4 ;

g [T1d] = eps1c.p3*p3.eps2c*eps4c.p1*p1.eps3c ;
g [T2d] = eps1c.p3*p3.eps2c*eps4c.p1*eps3c.p2 ;
g [T3d] = eps1c.p3*p3.eps2c*eps4c.p2*p1.eps3c ;
g [T4d] = eps1c.p3*p3.eps2c*eps4c.p2*eps3c.p2 ;
g [T5d] = eps1c.p3*p3.eps2c*eps4c.eps3c ;
g [T6d] = eps1c.eps2c*eps4c.p1*p1.eps3c ;
g [T7d] = eps1c.eps2c*eps4c.p1*eps3c.p2 ;
g [T8d] = eps1c.eps2c*eps4c.p2*p1.eps3c ;
g [T9d] = eps1c.eps2c*eps4c.p2*eps3c.p2 ;
g [T10d] = eps1c.eps4c*p3.eps2c*p1.eps3c ;
g [T11d] = eps1c.eps4c*p3.eps2c*eps3c.p2 ;
g [T12d] = eps1c.eps3c*p3.eps2c*eps4c.p1 ;
g [T13d] = eps1c.eps3c*p3.eps2c*eps4c.p2 ;
g [T14d] = eps1c.p3*eps2c.eps4c*p1.eps3c ;
g [T15d] = eps1c.p3*eps2c.eps4c*eps3c.p2 ;
g [T16d] = eps1c.p3*eps2c.eps3c*eps4c.p1 ;
g [T17d] = eps1c.p3*eps2c.eps3c*eps4c.p2 ;
g [T18d] = eps1c.eps2c*eps4c.eps3c + eps1c.eps4c*eps2c.eps3c + eps1c.eps3c*eps2c.eps4c ;

* parity odd
g [T19] = eps1.va*eps2.p3*eps4.p1*p1.eps3 ;
g [T20] = eps1.va*eps2.p3*eps4.p1*eps3.p2 ;
g [T21] = eps1.va*eps2.p3*eps4.p2*p1.eps3 ;
g [T22] = eps1.va*eps2.p3*eps4.p2*eps3.p2 ;
g [T23] = eps1.p3*va.eps2*eps4.p1*p1.eps3 ;
g [T24] = eps1.p3*va.eps2*eps4.p1*eps3.p2 ;
g [T25] = eps1.p3*va.eps2*eps4.p2*p1.eps3 ;
g [T26] = eps1.p3*va.eps2*eps4.p2*eps3.p2 ;
g [T27] = eps1.p3*va.eps4*eps2.p3*p1.eps3 ;
g [T28] = eps1.p3*va.eps4*eps2.p3*eps3.p2 ;
g [T29] = eps1.p3*va.eps3*eps2.p3*eps4.p1 ;
g [T30] = eps1.p3*va.eps3*eps2.p3*eps4.p2 ;
g [T31] = (eps1.eps2*va.eps4*p1.eps3 + eps1.va*eps2.eps4*p1.eps3 + eps1.eps4*eps2.va*p1.eps3) ;
g [T32] = (eps1.eps2*va.eps4*p2.eps3 + eps1.va*eps2.eps4*p2.eps3 + eps1.eps4*eps2.va*p2.eps3) ;
g [T33] = (eps1.eps2*va.eps3*p1.eps4 + eps1.va*eps2.eps3*p1.eps4 + eps1.eps3*eps2.va*p1.eps4) ;
g [T34] = (eps1.eps2*va.eps3*p2.eps4 + eps1.va*eps2.eps3*p2.eps4 + eps1.eps3*eps2.va*p2.eps4) ;
g [T35] = (eps1.eps4*va.eps3*p3.eps2 + eps3.eps4*va.eps1*p3.eps2 + eps1.eps3*eps4.va*p3.eps2) ;
g [T36] = (eps2.eps4*va.eps3*p3.eps1 + eps3.eps4*va.eps2*p3.eps1 + eps2.eps3*eps4.va*p3.eps1) ;

g [T19d] =  - va.eps1c*p3.eps2c*p1.eps4c*p1.eps3c ;
g [T20d] =  - va.eps1c*p3.eps2c*p1.eps4c*p2.eps3c ;
g [T21d] =  - va.eps1c*p3.eps2c*p1.eps3c*p2.eps4c ;
g [T22d] =  - va.eps1c*p3.eps2c*p2.eps4c*p2.eps3c ;
g [T23d] =  - va.eps2c*p3.eps1c*p1.eps4c*p1.eps3c ;
g [T24d] =  - va.eps2c*p3.eps1c*p1.eps4c*p2.eps3c ;
g [T25d] =  - va.eps2c*p3.eps1c*p1.eps3c*p2.eps4c ;
g [T26d] =  - va.eps2c*p3.eps1c*p2.eps4c*p2.eps3c ;
g [T27d] =  - va.eps4c*p3.eps1c*p3.eps2c*p1.eps3c ;
g [T28d] =  - va.eps4c*p3.eps1c*p3.eps2c*p2.eps3c ;
g [T29d] =  - va.eps3c*p3.eps1c*p3.eps2c*p1.eps4c ;
g [T30d] =  - va.eps3c*p3.eps1c*p3.eps2c*p2.eps4c ;
g [T31d] = -(eps1c.eps2c*va.eps4c*p1.eps3c + eps1c.va*eps2c.eps4c*p1.eps3c + eps1c.eps4c*eps2c.va*p1.eps3c) ;
g [T32d] = -(eps1c.eps2c*va.eps4c*p2.eps3c + eps1c.va*eps2c.eps4c*p2.eps3c + eps1c.eps4c*eps2c.va*p2.eps3c) ;
g [T33d] = -(eps1c.eps2c*va.eps3c*p1.eps4c + eps1c.va*eps2c.eps3c*p1.eps4c + eps1c.eps3c*eps2c.va*p1.eps4c) ;
g [T34d] = -(eps1c.eps2c*va.eps3c*p2.eps4c + eps1c.va*eps2c.eps3c*p2.eps4c + eps1c.eps3c*eps2c.va*p2.eps4c) ;
g [T35d] = -(eps1c.eps4c*va.eps3c*p3.eps2c + eps3c.eps4c*va.eps1c*p3.eps2c + eps1c.eps3c*eps4c.va*p3.eps2c) ;
g [T36d] = -(eps2c.eps4c*va.eps3c*p3.eps1c + eps3c.eps4c*va.eps2c*p3.eps1c + eps2c.eps3c*eps4c.va*p3.eps1c) ;

#define NUM "36"

.sort

print +s ;
.sort
*************************************************************************** Write them to a file
#do i=1,'NUM'
#write <tensors.inc> "g [T'i'] = %E ;" [T'i'] *AG: %E print the tensor ith as it appear in the code
#enddo

#write <tensors.inc> ""

#do i=1,'NUM'
#write <tensors.inc> "g [T'i'd] = %E ;" [T'i'd]
#enddo

.store

*************************************************************************** Define the decomposition of the amplitude and the projectors
*AG: definition of the amplitude: tensor times form factor
g [Amp] =
#do i = 1,'NUM'
+ F'i'*[T'i']
#enddo
;

*AG: definition of the projectors: see 2.15 of https://arxiv.org/pdf/1503.08835. Notice that the expression of the projection is completely general
g [Proj] =
#do i = 1,'NUM'
+ c'i'*[T'i'd]
#enddo
;

b c1,...,c'NUM' ;
.store

*************************************************************************** Define tensor metric, i.e. the product T_i^\dagger\cdot T_j. Notice that (M_ij)^-1=c_{ij}

#message The tensor metric is given by

#do i=1,'NUM'
#do j=1,'NUM'
g [M'i'p'j'] = [T'i'd]*[T'j'] ;
#enddo
#enddo
;

*************************************************************************** Apply the general projector on the tensor structure and obtain the system of equations

g expr = [Proj]*[Amp] ;

*************************************************************************** Compute fermion polarisation sums.

#call kinematics
#call gaugechoice
id eps1.p6?*eps1c.p7? = (-p6.p7 +2*(p1.p6*p2.p7+p1.p7*p2.p6)/(s) )  ;
id eps1.p6?*eps1c.p6? = (-p6.p6 +4*(p1.p6*p2.p6)/(s) )              ;
id eps1.eps1c = 2-d                                                 ;

#call kinematics
#call gaugechoice
id eps2.p6?*eps2c.p7? = (-p6.p7 +2*(p1.p6*p2.p7+p1.p7*p2.p6)/(s) )  ;
id eps2.p6?*eps2c.p6? = (-p6.p6 +4*(p1.p6*p2.p6)/(s) )              ;
id eps2.eps2c = 2-d                                                 ;

#call kinematics
#call gaugechoice
id eps3.p6?*eps3c.p7? = (-p6.p7 + (p3.p6*p3.p7)/(m32) )  ;
id eps3.p6?*eps3c.p6? = (-p6.p6 + (p3.p6*p3.p6)/(m32) )  ;
id eps3.eps3c = 1-d ;
#call kinematics
#call gaugechoice
id eps4.p6?*eps4c.p7? = (-p6.p7 + (p4.p6*p4.p7)/(m42) )  ;
id eps4.p6?*eps4c.p6? = (-p6.p6 + (p4.p6*p4.p6)/(m42) )  ;
id eps4.eps4c = 1-d ;
#call kinematics
#call gaugechoice

b F1,...,F'NUM' ;
print +s ;
.store

*************************************************************************** Export the tensor metric for use in Mathematica
#message Exporting tensor metric for use in Mathematica

#do i=1,'NUM'
#do j=1,'NUM'
l [TM'i''j'] = [M'i'p'j'] ;
#enddo
#enddo
;
.sort

format mathematica ;

*#write <FormFactors/TensorMetric.m> "TM[1][1] = ( %E ) ; " [TM11] ;
*#close <FormFactors/TensorMetric.m>

#do i=1,36
#do j=1,36
#write <TensorMetric_ggVV.m> "TM['i']['j'] = ( %E ) ; " [TM'i''j'] ;
#write <TensorMetric_ggVV.m> "" ;
#enddo
#enddo
#close <TensorMetric_ggVV.m>

format normal ;

.store

*************************************************************************** Construct the equations for the coefficients

#do i=1,'NUM'
b F1,...,F'NUM' ;
.sort

g eq'i' =Label('i')*( expr[F'i'] ) ;
.sort
#enddo

b Label,c1,...,c'NUM' ;
.sort
collect Num ;
.sort
id Num(n1?) = Ratio(n1,1) ;
.sort

b Label ;
.store

*************************************************************************** Solve now the equations and construct the projectors

#message We construct the projectors that single out the form factors

*system number 'N'
#do N = 1,'NUM'

#message Start Projector 'N'

#do i=1,'NUM'
#if('i'== 'N')
l zero'i' = eq'i' - Label('i')*1 ;
#else
l zero'i' = eq'i' ;
#endif
#enddo
.sort

polyratfun Ratio ;
.sort

#do i=1,'NUM'
l [c'i'] = c'i' ;
#enddo
.sort

*solve the equations

#do i=1,'NUM'

#$NN = 1 ;
#$DD = 1 ;
#$Eq = -7 ;
.sort
id Label(n1?$Eq)*c'i'*Ratio(n?$NN,m?$DD) = Label(n1)*c'i'*Ratio(n,m) ;
.sort

*#message use equation '$Eq' for coefficient 'i' of projector 'N'
id c'i' = 0 - zero'$Eq'*[]*Ratio($DD,$NN) ;

id Label('$Eq') = 1 ;
id []*c'i' = 0 ;
id [] = 1 ;
.sort

#enddo

polyratfun ;
.sort
id Ratio(n1?,n2?) = Num(n1)*Den(n2) ;
.sort
#call FactArgDen
#call FactArgNum

b Num;
.sort
id Num(n1?) = n1;
.sort

*Write the projector to a file
#write <projectorscoeffs.inc> "**** Projector 'N' ****"
#write <projectorscoeffs.inc> ""
#do i=1,'NUM'
#write <projectorscoeffs.inc> "id c'i'p'N' = %E ; " [c'i']
#write <projectorscoeffs.inc> ""
#enddo

print +s ;
.store

#enddo
.store

******** Now define the projectors with the coefficients calculated
#message We found the projectors

#do i1=1,'NUM'
g [P'i1'] =
#do i2=1,'NUM'
+ c'i2'p'i1'*[T'i2'd]
#enddo
;
#enddo

#include projectorscoeffs.inc

b Den ;
print +s ;
.store

*************************************************************************** We perform a series of consistency checks
#message We check that the coefficients defining the projectors define indeed the inverse of the tensor metric
#do i1=1,'NUM'
#do i2=1,'NUM'
g [cM'i1'p'i2'] =
#do i3=1,'NUM'
+ c'i1'p'i3'*[M'i3'p'i2']
#enddo
;
#enddo
#enddo

#include projectorscoeffs.inc

b Den,Num ;
.sort
collect Num ;
.sort
id Num(n1?) = Ratio(n1,1) ;
id Den(n1?) = Ratio(1,n1) ;
.sort
polyratfun Ratio ;
.sort
polyratfun ;
.sort

id Ratio(n1?,n1?) = 1 ;
id Ratio(n1?,n2?) = Num(n1)*Den(n2) ;
.sort
#call FactArgNum
#call FactArgDen

print +s;
.store

*Check that the projectors single out the correct form factors
#message We check that the projectors single out the correct form factors

#do i=1,'NUM'

l [FormFactor'i'] = [P'i']*[Amp] ;
.sort
#call kinematics
#call gaugechoice
id eps1.p6?*eps1c.p7? = (-p6.p7 +2*(p1.p6*p2.p7+p1.p7*p2.p6)/(s) )  ;
id eps1.p6?*eps1c.p6? = (-p6.p6 +4*(p1.p6*p2.p6)/(s) )              ;
id eps1.eps1c = 2-d                                                 ;
#call kinematics
#call gaugechoice
id eps2.p6?*eps2c.p7? = (-p6.p7 +2*(p1.p6*p2.p7+p1.p7*p2.p6)/(s) )  ;
id eps2.p6?*eps2c.p6? = (-p6.p6 +4*(p1.p6*p2.p6)/(s) )              ;
id eps2.eps2c = 2-d                                                 ;
#call kinematics
#call gaugechoice
id eps3.p6?*eps3c.p7? = (-p6.p7 + (p3.p6*p3.p7)/(m32) )  ;
id eps3.p6?*eps3c.p6? = (-p6.p6 + (p3.p6*p3.p6)/(m32) )  ;
id eps3.eps3c = 1-d ;
#call kinematics
#call gaugechoice
id eps4.p6?*eps4c.p7? = (-p6.p7 + (p4.p6*p4.p7)/(m42) )  ;
id eps4.p6?*eps4c.p6? = (-p6.p6 + (p4.p6*p4.p6)/(m42) )  ;
id eps4.eps4c = 1-d ;
#call kinematics
#call gaugechoice
#enddo

b Den,Num ;
.sort
collect Num ;
.sort
id Num(n1?) = Ratio(n1,1) ;
id Den(n1?) = Ratio(1,n1) ;
.sort
polyratfun Ratio ;
.sort
polyratfun ;
.sort

id Ratio(n1?,n1?) = 1 ;
id Ratio(n1?,n2?) = Num(n1)*Den(n2) ;
.sort
#call FactArgNum
#call FactArgDen

b F1,...,F'NUM' ;
print +s ;

.end
