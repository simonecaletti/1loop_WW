repeat id Label(PL00M)*k1.k1 = Label(PL00M)*(D1)                                             ;
repeat id Label(PL00M)*k1.p1 = Label(PL00M)*(  1/2*D1 - 1/2*D2)                              ;
repeat id Label(PL00M)*k1.p2 = Label(PL00M)*(  1/2*D2 - 1/2*D3 + 1/2*s)                      ;
repeat id Label(PL00M)*k1.p3 = Label(PL00M)*(  1/2*D3 - 1/2*D4 - 1/2*s + 1/2*m42 - 1/2*mt^2) ;

repeat id Label(PL00Mx12)*k1.k1 = Label(PL00Mx12)*(D1)                                                ;
repeat id Label(PL00Mx12)*k1.p1 = Label(PL00Mx12)*(  1/2*D2 - 1/2*D3 + 1/2*s)                         ;
repeat id Label(PL00Mx12)*k1.p2 = Label(PL00Mx12)*(  1/2*D1 - 1/2*D2)                                 ;
repeat id Label(PL00Mx12)*k1.p3 = Label(PL00Mx12)*(  1/2*D3 - 1/2*D4 - 1/2*s+1/2*m42 - 1/2*mt^2)      ;

repeat id Label(PL0M0)*k1.k1 = Label(PL0M0)*(mt^2+D1)                                               ;
repeat id Label(PL0M0)*k1.p1 = Label(PL0M0)*( -1/2*D2 + 1/2*D1)                                     ;
repeat id Label(PL0M0)*k1.p2 = Label(PL0M0)*(  1/2*D2 - 1/2*D3 + 1/2*s)                             ;
repeat id Label(PL0M0)*k1.p3 = Label(PL0M0)*(  1/2*D3 - 1/2*D4 - 1/2*s + 1/2*m42 + 1/2*mt^2)        ;

repeat id Label(PL0M0x12)*k1.k1 = Label(PL0M0x12)*(mt^2+D1)                                         ;
repeat id Label(PL0M0x12)*k1.p2 = Label(PL0M0x12)*( -1/2*D2 + 1/2*D1)                               ;
repeat id Label(PL0M0x12)*k1.p1 = Label(PL0M0x12)*(  1/2*D2 - 1/2*D3 + 1/2*s)                       ;
repeat id Label(PL0M0x12)*k1.p3 = Label(PL0M0x12)*(  1/2*D3 - 1/2*D4 - 1/2*s + 1/2*m42 + 1/2*mt^2)  ;

repeat id Label(NPL0M)*k1.k1 = Label(NPL0M)*(D1)                                                    ;
repeat id Label(NPL0M)*k1.p1 = Label(NPL0M)*( -1/2*D2 + 1/2*D1)                                     ;
repeat id Label(NPL0M)*k1.p2 = Label(NPL0M)*(  1/2*D3 - 1/2*D4 - 1/2*(m32-s-u) )      ;
repeat id Label(NPL0M)*k1.p3 = Label(NPL0M)*(  1/2*D2 - 1/2*D3 + 1/2*(m42+m32-s-u)-1/2*mt^2)        ;


repeat id Label(NPL0Mx12)*k1.k1 = Label(NPL0Mx12)*(D1)                                        ;
repeat id Label(NPL0Mx12)*k1.p1 = Label(NPL0Mx12)*(  1/2*D3 - 1/2*D4 - 1/2*u + 1/2*m42)      ;
repeat id Label(NPL0Mx12)*k1.p2 = Label(NPL0Mx12)*( -1/2*D2 + 1/2*D1)                         ;
repeat id Label(NPL0Mx12)*k1.p3 = Label(NPL0Mx12)*(  1/2*D2 - 1/2*D3 + 1/2*u-1/2*mt^2)        ;
