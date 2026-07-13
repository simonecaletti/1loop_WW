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
g [T19] = eps1.va*p3.eps2*eps4.p1*p1.eps3 ;
g [T20] = eps1.va*p3.eps2*eps4.p1*eps3.p2 ;
g [T21] = eps1.va*p3.eps2*eps4.p2*p1.eps3 ;
g [T22] = eps1.va*p3.eps2*eps4.p2*eps3.p2 ;
g [T23] = eps1.p3*eps2.va*eps4.p1*p1.eps3 ;
g [T24] = eps1.p3*eps2.va*eps4.p1*eps3.p2 ;
g [T25] = eps1.p3*eps2.va*eps4.p2*p1.eps3 ;
g [T26] = eps1.p3*eps2.va*eps4.p2*eps3.p2 ;
g [T27] = eps1.p3*p3.eps2*eps4.va*p1.eps3 ;
g [T28] = eps1.p3*p3.eps2*eps4.va*eps3.p2 ;
g [T29] = eps1.p3*p3.eps2*eps4.p1*eps3.va ;
g [T30] = eps1.p3*p3.eps2*eps4.p2*eps3.va ;
g [T31] = eps1.eps2*eps4.va*p1.eps3 + eps1.eps4*eps2.va*p1.eps3 + eps1.va*
      eps2.eps4*p1.eps3 ;
g [T32] = eps1.eps2*eps4.va*eps3.p2 + eps1.eps4*eps2.va*eps3.p2 + eps1.va*
      eps2.eps4*eps3.p2 ;
g [T33] = eps1.eps2*eps4.p1*eps3.va + eps1.eps3*eps2.va*eps4.p1 + eps1.va*
      eps2.eps3*eps4.p1 ;
g [T34] = eps1.eps2*eps4.p2*eps3.va + eps1.eps3*eps2.va*eps4.p2 + eps1.va*
      eps2.eps3*eps4.p2 ;
g [T35] = eps1.eps4*p3.eps2*eps3.va + eps1.eps3*p3.eps2*eps4.va + eps1.va*
      p3.eps2*eps4.eps3 ;
g [T36] = eps1.p3*eps2.eps4*eps3.va + eps1.p3*eps2.eps3*eps4.va + eps1.p3*
      eps2.va*eps4.eps3 ;

g [T1d] = p3.eps1c*p3.eps2c*p1.eps4c*p1.eps3c ;
g [T2d] = p3.eps1c*p3.eps2c*p1.eps4c*p2.eps3c ;
g [T3d] = p3.eps1c*p3.eps2c*p1.eps3c*p2.eps4c ;
g [T4d] = p3.eps1c*p3.eps2c*p2.eps4c*p2.eps3c ;
g [T5d] = p3.eps1c*p3.eps2c*eps4c.eps3c ;
g [T6d] = p1.eps4c*p1.eps3c*eps1c.eps2c ;
g [T7d] = p1.eps4c*p2.eps3c*eps1c.eps2c ;
g [T8d] = p1.eps3c*p2.eps4c*eps1c.eps2c ;
g [T9d] = p2.eps4c*p2.eps3c*eps1c.eps2c ;
g [T10d] = p3.eps2c*p1.eps3c*eps1c.eps4c ;
g [T11d] = p3.eps2c*p2.eps3c*eps1c.eps4c ;
g [T12d] = p3.eps2c*p1.eps4c*eps1c.eps3c ;
g [T13d] = p3.eps2c*p2.eps4c*eps1c.eps3c ;
g [T14d] = p3.eps1c*p1.eps3c*eps2c.eps4c ;
g [T15d] = p3.eps1c*p2.eps3c*eps2c.eps4c ;
g [T16d] = p3.eps1c*p1.eps4c*eps2c.eps3c ;
g [T17d] = p3.eps1c*p2.eps4c*eps2c.eps3c ;
g [T18d] = eps1c.eps2c*eps4c.eps3c + eps1c.eps4c*eps2c.eps3c + eps1c.eps3c*
      eps2c.eps4c ;
g [T19d] =  - p3.eps2c*p1.eps4c*p1.eps3c*eps1c.va ;
g [T20d] =  - p3.eps2c*p1.eps4c*p2.eps3c*eps1c.va ;
g [T21d] =  - p3.eps2c*p1.eps3c*p2.eps4c*eps1c.va ;
g [T22d] =  - p3.eps2c*p2.eps4c*p2.eps3c*eps1c.va ;
g [T23d] =  - p3.eps1c*p1.eps4c*p1.eps3c*eps2c.va ;
g [T24d] =  - p3.eps1c*p1.eps4c*p2.eps3c*eps2c.va ;
g [T25d] =  - p3.eps1c*p1.eps3c*p2.eps4c*eps2c.va ;
g [T26d] =  - p3.eps1c*p2.eps4c*p2.eps3c*eps2c.va ;
g [T27d] =  - p3.eps1c*p3.eps2c*p1.eps3c*eps4c.va ;
g [T28d] =  - p3.eps1c*p3.eps2c*p2.eps3c*eps4c.va ;
g [T29d] =  - p3.eps1c*p3.eps2c*p1.eps4c*eps3c.va ;
g [T30d] =  - p3.eps1c*p3.eps2c*p2.eps4c*eps3c.va ;
g [T31d] =  - p1.eps3c*eps1c.eps2c*eps4c.va - p1.eps3c*eps1c.eps4c*eps2c.va - 
      p1.eps3c*eps1c.va*eps2c.eps4c ;
g [T32d] =  - p2.eps3c*eps1c.eps2c*eps4c.va - p2.eps3c*eps1c.eps4c*eps2c.va - 
      p2.eps3c*eps1c.va*eps2c.eps4c ;
g [T33d] =  - p1.eps4c*eps1c.eps2c*eps3c.va - p1.eps4c*eps1c.eps3c*eps2c.va - 
      p1.eps4c*eps1c.va*eps2c.eps3c ;
g [T34d] =  - p2.eps4c*eps1c.eps2c*eps3c.va - p2.eps4c*eps1c.eps3c*eps2c.va - 
      p2.eps4c*eps1c.va*eps2c.eps3c ;
g [T35d] =  - p3.eps2c*eps1c.eps4c*eps3c.va - p3.eps2c*eps1c.eps3c*eps4c.va - 
      p3.eps2c*eps1c.va*eps4c.eps3c ;
g [T36d] =  - p3.eps1c*eps2c.eps4c*eps3c.va - p3.eps1c*eps2c.eps3c*eps4c.va - 
      p3.eps1c*eps2c.va*eps4c.eps3c ;
