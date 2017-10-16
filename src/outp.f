!
! Copyright 1996-2017 the Authors
!
! Licensed under the EUPL, Version 1.1 only (the "Licence");
!
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!   https://joinup.ec.europa.eu/software/page/eupl
!
! Unless required by applicable law or agreed to in writing,
! software distributed under the Licence is distributed on an
! "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
! either express or implied.
!
! See the Licence for the specific language governing permissions
! and limitations under the Licence.



! outp.f : output
! box: the initial (const*) and overview (prof*) output SRs are not adjusted, 
!      i.e. they produce a heap of output that's irrelevant. To avoid huge 
!     output files, the size of the plou* SRs has been adjusted by using
!     n_bl and n_bln. Output starts at k=1 to save the deposited/surface 
!     values as well.

! This file contains the following subroutines:
!
!     - outm: restart files, meteorological data
!     - outc: restart files, chemical data
!
!     - ploutj: output of photolysis rates up to level n_bln
!
!     - ploutc: utilitary subroutine to call the following subroutines:
!              - ploutcg: output of gas phase chemical species
!              - ploutci:
!              - ploutcl:
!              - ploutcr:
!              - ploutcgr: output of instantaneous rates
!              - ploutcgs:
!     -
! (to be continued)

      subroutine outm
!
! Description:
!    profiles of meteorological data needed to restart the program
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      08/2016   Use module for parameters                <Josue Bock>
!                    Comments / header
!                    Removed one unused argument
!
! 1.1       ?        Output every 12h with increasing names   <Roland von Glasow>
!                    Changed most common blocks
!
! 1.0       ?        Original code.                           <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:
      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     n1,
     &     nb,
     &     nka,
     &     nkt,
     &     mb

      implicit double precision (a-h,o-z)

! Local scalars:
      character (len=10) fname
      character (len=1) sub

! Common blocks:
      common /cb11/ totrad (mb,n1)
      double precision totrad

      common /cb18/ alat,declin                ! for the SZA calculation
      double precision alat,declin

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
      common /cb44/ r0,r1,g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhow,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      common /cb45/ u(n),v(n),w(n)
      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb48/ sk,sl,dtrad(n),dtcon(n)
      double precision sk, sl, dtrad, dtcon

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      common /cb55/ dtrad0(n),dtrad1(n),sk0,sl0,sk1,sl1,time0,time2
      common /cb63/ fcs(nka),xmol3(nka)

!- End of header ---------------------------------------------------------------

      fname='rstm .dat'
      sub='_'

      if (lday*24+lst.eq.12) sub='A'
      if (lday*24+lst.eq.24) sub='B'
      if (lday*24+lst.eq.36) sub='C'
      if (lday*24+lst.eq.48) sub='D'
      if (lday*24+lst.eq.60) sub='E'
      if (lday*24+lst.eq.72) sub='F'
      if (lday*24+lst.eq.84) sub='G'
      if (lday*24+lst.eq.96) sub='H'
      if (lday*24+lst.eq.108) sub='I'
      if (lday*24+lst.eq.120) sub='J'
      if (lday*24+lst.eq.132) sub='K'
      if (lday*24+lst.eq.144) sub='L'
      if (lday*24+lst.eq.156) sub='M'
      if (lday*24+lst.eq.168) sub='N'
      if (lday*24+lst.eq.180) sub='O'
      if (lday*24+lst.eq.192) sub='P'
      fname(5:5)=sub
 3000 continue
      open (15,file=fname,status='unknown',form='unformatted',err=3000)
! double precision arrays
      write (15)
     &     atkm,atkh,b0m,dfddt,dtrad,dtrad0,dtrad1,eb,ff,fcs,feu,fsum,
     &     gh,p,rho,t,talt,tb,tke,tkep,theta,totrad,u,v,w,xl,xm1,xm1a,
     &     xm2,xmol3,
! double precision single vars
     &     a0m,alat,declin,ds1,ds2,reif,sk,sk0,sk1,sl,sl0,sl1,tau,time0,
     &     time2,trdep,
! integer arrays
     &     nar,
! integer single vars
     &     it,lcl,lct,lday,lmin,lst
      close (15)
 
      end subroutine outm
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine outc
!
! Description:
!    output of chemical data needed to restart the program
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      08/2016   Use module for parameters                <Josue Bock>
!                    Comments / header
!                    Removed one unused argument
!                    Removed unnecessary and outdated output
!
! 1.1       ?        Output every 12h with increasing names   <Roland von Glasow>
!                    Changed most common blocks
!
! 1.0       ?        Original code.                           <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:
      USE gas_common, ONLY:
! Imported Array Variables with intent (in):
     &     s1,
     &     es1,
     &     s3

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nka,
     &     nkt,
     &     nkc,
     &     nphrxn,
     &     nlev,
     &     nrxn

      implicit double precision (a-h,o-z)

! Local scalars:
      character (len=10) fname
      character (len=1) sub

! Common blocks:
      common /band_rat/ photol_j(nphrxn,n)
      double precision photol_j

      common /blck01/ am3(n),cm3(n)
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /blck78/ sa1(nka,j2),sac1(nka,j2)

      common /budg/ bg(2,nrxn,nlev),il(nlev)
      double precision bg
      integer il

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /kinv_i/ kinv
      common /kpp_crys/ xcryssulf,xcrysss,xdelisulf,xdeliss
      common /kpp_l1/ cloudt(nkc,n)
      logical cloudt

      common /kpp_mol/ xgamma(nf,j6,nkc)
      common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)

!- End of header ---------------------------------------------------------------

      fname='rstc .dat'
      sub='_'
      if (lday*24+lst.eq.12) sub='A'
      if (lday*24+lst.eq.24) sub='B'
      if (lday*24+lst.eq.36) sub='C'
      if (lday*24+lst.eq.48) sub='D'
      if (lday*24+lst.eq.60) sub='E'
      if (lday*24+lst.eq.72) sub='F'
      if (lday*24+lst.eq.84) sub='G'
      if (lday*24+lst.eq.96) sub='H'
      if (lday*24+lst.eq.108) sub='I'
      if (lday*24+lst.eq.120) sub='J'
      if (lday*24+lst.eq.132) sub='K'
      if (lday*24+lst.eq.144) sub='L'
      if (lday*24+lst.eq.156) sub='M'
      if (lday*24+lst.eq.168) sub='N'
      if (lday*24+lst.eq.180) sub='O'
      if (lday*24+lst.eq.192) sub='P'
      fname(5:5)=sub
 3000 continue
      open (16,file=fname,status='unknown',form='unformatted',err=3000)
! double precision arrays
      write (16) am3,cm,cm3,conv2,cw,es1,photol_j,rc,s1,s3,sa1,
     &     sac1,sl1,sion1,vd,vdm,vt,xgamma,
! double precision, single values
     &     xcryssulf,xcrysss,xdelisulf,xdeliss,
! logicals
     &     cloudt,
! integers
     &     il,kinv,lday,lmin,lst
      close (16)

      end subroutine outc
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploutj (fogtype,n_bln)
!
! Description:
!    output of photolysis rates
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      08/2016   Use module for parameters                <Josue Bock>
!                    Comments / header
!                    Removed one unused argument
!                    Cleaning
!
! 1.0       ?        Original code.                          <Roland von Glasow>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:
      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nphrxn

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      character (len=1) fogtype
      integer n_bln

! Local scalars:
      character (len=10) fname
      integer j,k
      double precision xday,xst,xmin

! Local arrays:
      double precision i0                    ! the local array where are written the photolysis rates,
      dimension i0(nphrxn,n_bln) !    up to the selected level (1:n_bln)

! Common blocks:
      common /band_rat/ photol_j(nphrxn,n)
      double precision photol_j

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

!- End of header ---------------------------------------------------------------

      fname='jra .out'
      fname(4:4)=fogtype

      do k=1,n_bln
         do j=1,nphrxn
            i0(j,k)=photol_j(j,k)
         enddo
      enddo

      xday=float(lday)
      xst=float(lst)
      xmin=float(lmin)

 3000 continue
      open (69, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
!      write (69) lday,lst,lmin,i0
      write (69) xday,xst,xmin,i0
      close (69)

      end subroutine ploutj
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploutc (fogtype,mic,n_bl,n_bl8)
      implicit none
      character (len=1) fogtype
      logical mic
      integer n_bl,n_bl8
      call ploutcg (fogtype,n_bl)
      if (mic) call ploutci (fogtype,n_bl)
      if (mic) call ploutcl (fogtype,n_bl)
      call ploutcr (fogtype,n_bl)
      call ploucgr (fogtype,n_bl8)
      call ploucgs (fogtype,n_bl)

      end subroutine ploutc
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploutcg (fogtype,n_bl)
!
! Description:
!    output of gas phase chemical species s1:
!
!    1. NO      2. NO2     3. HNO3     4. NH3     5. SO2      6. H2SO4
!    7. O3      8. CH4     9. C2H6    10. C3H8   11. ALKA    12. ETHE
!   13. ALKE   14. AROM   15. HCOOH   16. ACTA   17. HCHO    18. ALD2
!   19. H2O2   20. CH3OOH 21. HONO    22. PAN    23. TPAN    24. KET
!   25. CRES   26. DIAL   27. GLYX    28. MGLY   29. NH4NO3  30. HCl
!   31. R3N2   32. RAN2   33. RAN1    34. N2O5   35. HNO4    36. NO3
!   37. DMS    38. HOCl   39. ClNO2   40. ClNO3  41. Cl2     42. HBr
!   43. HOBr   44. BrNO2  45. BrNO3   46. Br2    47. BrCl    48. HI
!   49. HOI    50. I2O2   51. INO2    52. INO3   53. I2      54. ICl
!   55. IBr    56. CH3I   57. CH2I2   58. CH2ClI 59. C3H7I   60. DMSO
!   61. CH3SO2 62. CH3SO3 63. CH3SO3H 64. CO     65. Cl2O2   66. DMOO    
!   67. CH3S   68. CH3SO  69. MSIA    70. DMSO2  71. CH2BrI  72. CHBr2I
!   73. C2H5I  74. HIO3   75. NUCV    76. SO3    77. HOSO2   78. CO2
!   79. I2O    80. I2O3   81. I2O4    82. I2O5   83. INO     84. Br2O
!   85. ClONO  86. ClO3   87. Cl2O3   88. CH3OH  89. C2H5OH  90. H2
!   91. NHS    92. RCl    93. RBr     94. XOR    95. SOR     96. SPAN
!   97. Hg     98. HgO    99. HgCl   100. HgCl2 101. HgBr   102. HgBr2
!   102 -- 130 undefined
!
! DMOO = CH3SCH2OO, MSIA = CH3S(O)OH
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      08/2016   Use module          <Josue Bock>
!                    Comments / header
!
! 1.0       ?        Original code.      <Roland von Glasow>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:
      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
! Imported Array Variables with intent (in):
     &     s1

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      character (len=1) fogtype
      integer n_bl

! Local scalars:
      character (len=10) fname
      integer j,k, nmax

! Local arrays:
      double precision i0
      dimension i0(n_bl,j1)

! Common blocks:
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

!- End of header ---------------------------------------------------------------

      nmax = n_bl
      if (n_bl.gt.3) nmax = n_bl-2
      do k=1,nmax
         do j=1,j1
            i0(k,j)=s1(j,k)
         enddo
      enddo
! only for 1D:
      if (n_bl.gt.3) then
         do j=1,j1
            i0(nf-1,j)=etw(lcl)
            i0(nf,j)=etw(lct)
         enddo
      endif
      fname='sg1 .out'
      fname(4:4)=fogtype
 3000 continue
      open (61, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
      write (61) lday,lst,lmin,i0
      close (61)

      end subroutine ploutcg
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploutci (fogtype,n_bl)
!
! Description:
!    output of ions sion1:
!
!    1. H+       2. NH4+     3. OH-      4. CH2OHSO3-  5. HSO3-
!    6. SO3=     7. SO4-     8. SO4=     9. HCO3-     10. CO3-
!   11. O2-     12. NO2-    13. NO3-    14. Cl-       15. Cl2-
!   16. HCOO-   17. Fe3+    18. Mn2+    19. HSO4-     20. Na+ (check electroneg)
!   21. NO4-    22. ClO-    23. ClOH-   24. Br-       25. Br2-
!   26. BrO-    27. BrOH-   28. BrCl2-  29. Br2Cl-    30. CH3SO3-     
!   31. HSO5-   32. SO3-    33. SO5-    34. I-        35. IO2-        
!   36. IO3-    37. ICl2-   38. IBr2-   39. MS-       40. Hg+
!   41. Hg2+    42. HgOH+   43. HgCl+   44. HgCl3-    45. HgCl42-
!   46. HgBr+   47. HgBr3-  48. HgBr42- 49. Hg(SO3)22- 50. --
!
! MS- = CH3S(O)O-, CH3SO3- = CH3S(OO)O-
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      08/2016   Use module          <Josue Bock>
!                    Comments / header
!
! 1.0       ?        Original code.      <Roland von Glasow>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:
      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nkc,
     &     nf,
     &     n

      implicit double precision (a-h,o-z)

! Subroutine arguments
! Scalar arguments with intent(in):
      character (len=1) fogtype
      integer n_bl

! Local scalars:
      character (len=10) fname

! Local arrays:
      double precision i0
      dimension i0(n_bl,j6,nkc_l)

! Common blocks:
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /liq_pl/ nkc_l

!- End of header ---------------------------------------------------------------

      nmax = n_bl
      if (n_bl.gt.3) nmax = n_bl-2
      do k=1,nmax
         do j=1,j6
            do i=1,nkc_l
               i0(k,j,i)=sion1(j,i,k)
            enddo
         enddo
      enddo
! only for 1D:
      if (n_bl.gt.3) then
         do j=1,j6
            do i=1,nkc_l
               i0(nf-1,j,i)=etw(lcl)
               i0(nf,j,i)=etw(lct)
            enddo
         enddo
      endif
      fname='ion .out'
      fname(4:4)=fogtype
 3000 continue
      open (62, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
      write (62) lday,lst,lmin,i0
      close (62)

      end subroutine ploutci
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploutcl (fogtype,n_bl)
!
! Description:
!    output of liquid phase species sl1:
!
!    1. NO      2. NO2     3. HNO3     4. NH3     5. SO2      6. H2SO4
!    7. O3      8. CH4     9. C2H6    10. C3H8   11. ALKA    12. ETHE
!   13. ALKE   14. AROM   15. HCOOH   16. ACTA   17. HCHO    18. ALD2
!   19. H2O2   20. CH3OOH 21. HONO    22. PAN    23. TPAN    24. KET
!   25. CRES   26. DIAL   27. GLYX    28. MGLY   29. NH4NO3  30. HCl
!   31. R3N2   32. RAN2   33. RAN1    34. N2O5   35. HNO4    36. NO3
!   37. DMS    38. HOCl   39. ClNO2   40. ClNO3  41. Cl2     42. HBr
!   43. HOBr   44. BrNO2  45. BrNO3   46. Br2    47. BrCl    48. HI
!   49. HOI    50. I2O2   51. INO2    52. INO3   53. I2      54. ICl
!   55. IBr    56. CH3I   57. CH2I2   58. CH2ClI 59. C3H7I   60. DMSO
!   61. CH3SO2 62. CH3SO3 63. CH3SO3H 64. CO     65. Cl2O2   66. DMOO    
!   67. CH3S   68. CH3SO  69. MSIA    70. DMSO2  71. CH2BrI  72. CHBr2I
!   73. C2H5I  74. HIO3   75. NUCV    76. SO3    77. HOSO2   78. CO2
!   79. I2O    80. I2O3   81. I2O4    82. I2O5   83. INO     84. Br2O
!   85. ClONO  86. ClO3   87. Cl2O3   88. CH3OH  89. C2H5OH  90. H2
!   91. NHS    92. RCl    93. RBr     94. XOR    95. SOR     96. SPAN
!   97. Hg     98. HgO    99. HgCl   100. HgCl2 101. HgBr   102. HgBr2
!
! DMOO = CH3SCH2OO, MSIA = CH3S(O)OH.
!
!   j2-j3+1. --     j2-j3+2. OH    j2-j3+3. HO2    j2-j3+4. DOM     j2-j3+5.  HIO2
!   j2-j3+6. CH3OO  j2-j3+7. IO    j2-j3+8. Cl     j2-j3+9. Br      j2-j3+10. --
!   j2-j3+11.O2     j2-j3+12.OIO   j2-j3+13. HgOH2 j2-j3+14. HgOHCl j2-j3+15. HgSO3
!   j2-j3+16. HgOHBr j2-j3+17. --  j2-j3+18. --    j2-j3+19. --     j2-j3+20. --
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      08/2016   Use module          <Josue Bock>
!                    Comments / header
!
! 1.0       ?        Original code.      <Roland von Glasow>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:
      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nkc,
     &     nf,
     &     n,
     &     nka,
     &     nkt

 
      implicit double precision (a-h,o-z)

! Subroutine arguments
! Scalar arguments with intent(in):
      character (len=1) fogtype
      integer n_bl

! Local scalars:
      character (len=10) fname
      integer i,j,k, nmax

! Local arrays:
      double precision i0,irc,icw
      dimension i0(n_bl,j2,nkc_l),irc(n_bl,nkc_l),icw(n_bl,nkc_l)

! Common blocks:
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      double precision ff, fsum
      integer nar

      common /liq_pl/ nkc_l

!- End of header ---------------------------------------------------------------

      nmax = n_bl
      if (n_bl.gt.3) nmax = n_bl-2
      do k=1,nmax
         do j=1,j2
            do i=1,nkc_l
               i0(k,j,i)=sl1(j,i,k)
            enddo
         enddo
         x0=0.
         do ia=1,nka
            do jt=1,nkt
               x0=x0+ff(jt,ia,k)*en(ia)
            enddo
         enddo
         i0(k,j2,1)=x0
         i0(k,j2,2)=fsum(k)
      enddo
! only for 1D:
      if (n_bl.gt.3) then
         do j=1,j2
            do i=1,nkc_l
               i0(nf-1,j,i)=etw(lcl)
               i0(nf,j,i)=etw(lct)
            enddo
         enddo
      endif
      fname='sl1 .out'
      fname(4:4)=fogtype

      do kc=1,nkc_l
         do k=1,n_bl
            irc(k,kc)=rc(kc,k) ! jjb this transposition is maybe useless, just to stick to old write
            icw(k,kc)=cw(kc,k)
         enddo      
      enddo
 3000 continue
      open (63, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
      write (63) lday,lst,lmin,i0,irc,icw
      close (63)

      end subroutine ploutcl
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploutcr (fogtype,n_bl)


! output of radical species s3:
! species 1-3 are treated as long lived 
!    1. ----    2. ----    3. ---     4. OH      5. HO2
!    6. AHO2    7. MCO3    8. CH3OO   9. ETO2   10. KO2
!   11. R3O2   12. RAO2   13. TO2    14. TCO3   15. ZO2
!   16. EO2    17. PO2    18. CHO2   19. CRO2   20. PRN1
!   21. O(1D)  22. Cl     23. ClO    24. OClO   25. Br
!   26. BrO    27. I      28. IO     29. OIO    30. O3P
!   31. ClRO2  32. BrRO2  33. IRO2

! Modifications:
!     30-11-2016   J. Bock   implicit none, header
!     30-11-2016   J. Bock   correction in 1D case: indexes n_bl / n_bl-1
!                            instead of nf / nf-1
!     30-11-2016   J. Bock   loop order: priority for write (efficiency)
!     14-02-2017   J. Bock   renamed MO2 -> CH3OO
!                            MO2 was only used in the gas phase, while CH3OOlz was used in aq. phase,
!                            thus it could not be recognised by the newly developed routine searching
!                            for exchanged species.

      USE gas_common, ONLY:
! Imported Parameters:
     &     j5,
! Imported Array Variables with intent (in):
     &     s3

      USE global_params, ONLY :
! Imported Parameters:
     &     n

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      character (len=1) fogtype
      integer n_bl

! Local scalars:
      character (len=10) fname
      integer nmax, j, k

! Local arrays:
      double precision i0(n_bl,j5)

! Common blocks:
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

!- End of header ---------------------------------------------------------------

      nmax = n_bl
      if (n_bl.gt.3) nmax = n_bl-2

      do j=1,j5
         do k=1,nmax
            i0(k,j)=s3(j,k)
         enddo
      enddo
! only for 1D:
      if (n_bl.gt.3) then
         do j=1,j5
            i0(n_bl-1,j)=etw(lcl)
            i0(n_bl,j)=etw(lct)
         enddo
      endif

      fname='sr1 .out'
      fname(4:4)=fogtype
 3000 continue
      open (64, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
      write (64) lday,lst,lmin,i0
      close (64)

      end subroutine ploutcr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploucgr (fogtype,n_bl8)
!
! Description:
!    output of instantaneous rates [mol/(m3 s)]

! Modifications:
!     30-11-2016   J. Bock   USE module instead of local parameters
!     30-11-2016   J. Bock   implicit none

! Declarations:
! Modules used:
      USE global_params, ONLY :
! Imported Parameters:
     &     nlev,
     &     nrxn

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      character (len=1) fogtype
      integer n_bl8

! Local scalars:
      character (len=10) fname
      integer j, k

! Local arrays:
      double precision i0(n_bl8,nrxn)  ! the local array where are written the photolysis
                                       !   rates, up to the selected level (1:n_bl8)

! Common blocks:
      common /budg/ bg(2,nrxn,nlev),il(nlev)
      double precision bg
      integer il

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

!- End of header ---------------------------------------------------------------

      do j=1,nrxn
         do k=1,n_bl8
            i0(k,j)=bg(1,j,k)
         enddo
      enddo
      fname='gr .out'
      fname(3:3)=fogtype
 3000 continue
      open (66, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
      write (66) lday,lst,lmin,il,i0
      close (66)

      end subroutine ploucgr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploucgs (fogtype,n_bl)

!
! Description:
!     output of gas phase rates

! Modifications:
!     30-11-2016   J. Bock   USE module instead of local parameters
!     30-11-2016   J. Bock   implicit none, header


! Declarations:
! Modules used:
      USE global_params, ONLY :
! Imported Parameters:
     &     n

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      character (len=1) fogtype
      integer n_bl

! Local scalars:
      character (len=10) fname
      integer j, k

! Local arrays:
      double precision i0(n_bl,122),i1(n_bl,122)

! Common blocks:
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /budgs/ bgs(2,122,n)
      double precision bgs

!- End of header ---------------------------------------------------------------
      do j=1,122
         do k=1,n_bl
            i0(k,j)=bgs(1,j,k)
            i1(k,j)=bgs(2,j,k)
         enddo
      enddo
      fname='gs .out'
      fname(3:3)=fogtype
 3000 continue
      open (67, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
      write (67) lday,lst,lmin,i0,i1
      close (67)

      end subroutine ploucgs
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploutm (fogtype,n_bln)
! output of meteorological data

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     n1,
     &     nb,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)

      common /cb15/ fnseb,flgeg,hr(n1)
      double precision fnseb, flgeg, hr

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb48/ sk,sl,dtrad(n),dtcon(n)
      double precision sk, sl, dtrad, dtcon

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      double precision i0
      dimension i0(12,n_bln)
      character *10 fname
      character *1 fogtype
      do k=1,n_bln
         i0(1,k)=rho(k)
         i0(2,k)=atkh(k)
         i0(3,k)=theta(k)
         i0(4,k)=t(k)
         i0(5,k)=thetl(k)
         i0(6,k)=dtcon(k)
         i0(7,k)=dtrad(k)
         i0(8,k)=fsum(k)
         i0(9,k)=xm1(k)
         i0(10,k)=xm2(k)
         i0(11,k)=feu(k)
      enddo
! only for 1D:
      if (n_bln.gt.3) then
         do j=1,12
            i0(j,n-1)=etw(lcl)
            i0(j,n)=etw(lct)
         enddo
      endif
      fname='pm .out'
      fname(3:3)=fogtype
 3000 continue
      open (17, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
      write (17) lday,lst,lmin,i0
      close (17)
      fname(2:2)='b'
 3010 continue
      write (77,*) 'in pba'
      open (19, file=fname,status='old',form='unformatted',
     & position='append',err=3010)
      write (19) lday,lst,lmin,tb,eb,ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,
     & ajm,reif,tau,trdep,sl,sk,fnseb,flgeg
      close (19)

      end subroutine ploutm
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploutp (fogtype)
! output of particle distributions
      double precision ff,fsum,xm1,xm2,feu,dfddt,xm1a,xm2a

      parameter (nf=100,n=nf+50,nka=70,nkt=70)
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      dimension ff2(nka,nkt)
      character *10 fname
      character *1 fogtype
      fname='f1 .out'
      fname(3:3)=fogtype
 3000 continue
      open (41, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
      do kk=1,3
!         if (kk.eq.1) k=lcl-5
!         if (kk.eq.2) k=lcl
!         if (kk.eq.3) k=(lct+lcl)/2
!         k=max0(k,5)
         if (kk.eq.1) k=2
         if (kk.eq.2) k=12
         if (kk.eq.3) k=22
         do ia=1,nka
         do jt=1,nkt
            ff2(ia,jt)=ff(jt,ia,k)
         enddo
         enddo
         x0=xm2(k)
         x1=feu(k)
         x2=eta(k)
         write (41) lday,lst,lmin,k,x0,x1,x2,ff2
      enddo
      close (41)
      fname(2:2)='2'
 3010 continue
      open (42, file=fname,status='old',form='unformatted',
     & position='append',err=3010)
      do kk=1,3
!         if (kk.eq.1) k=lct-4
!         if (kk.eq.2) k=lct-2
!         if (kk.eq.3) k=lct
!         k=max0(k,5)
         if (kk.eq.1) k=32
         if (kk.eq.2) k=42
         if (kk.eq.3) k=52
         do ia=1,nka
         do jt=1,nkt
            ff2(ia,jt)=ff(jt,ia,k)
         enddo
         enddo
         x0=xm2(k)
         x1=feu(k)
         x2=eta(k)
         write (42) lday,lst,lmin,k,x0,x1,x2,ff2
      enddo
      close (42)
 3020 continue
      fname(2:2)='3'
      open (43, file=fname,status='old',form='unformatted',
     & position='append',err=3020)
      do kk=1,3
!         if (kk.eq.1) k=lct+2
!         if (kk.eq.2) k=lct+4
!         if (kk.eq.3) k=lct+6
!         k=max0(k,5)
         if (kk.eq.1) k=62
         if (kk.eq.2) k=72
         if (kk.eq.3) k=82
         do ia=1,nka
         do jt=1,nkt
            ff2(ia,jt)=ff(jt,ia,k)
         enddo
         enddo
         x0=xm2(k)
         x1=feu(k)
         x2=eta(k)
         write (43) lday,lst,lmin,k,x0,x1,x2,ff2
      enddo
      close (43)

      end subroutine ploutp
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploutr (fogtype,n_bln)
! output of radiation data

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     n1,
     &     n4,
     &     nka,
     &     nkt,
     &     mb,
     &     mbs

      implicit double precision (a-h,o-z)

      common /cb11/ totrad (mb,n1)
      double precision totrad

      common /cb15/ fnseb,flgeg,hr(n1)
      double precision fnseb, flgeg, hr

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /kurz/ fs1(n4),fs2(n4),totds(n4),ss(n4),fsn(n4),dtdts(n1)
      double precision fs1, fs2, totds, ss, fsn, dtdts

      common /lang/ fl1(n4),fl2(n4),fln(n4),dtdtl(n1)
      double precision fl1, fl2, fln, dtdtl

      double precision i0
      dimension i0(12,n_bln)
      character *10 fname
      character *1 fogtype
      do k=1,n_bln
         i0(1,k)=fs1(k)
         i0(2,k)=fs2(k)
         i0(3,k)=totds(k)
         i0(4,k)=dtdts(k)
         i0(5,k)=fl1(k)
         i0(6,k)=fl2(k)
         i0(7,k)=dtdtl(k)
         do ib=1,mbs
            i0(8,k)=i0(8,k)+totrad(ib,k)
         enddo
         do ib=mbs+1,mb
            i0(9,k)=i0(9,k)+totrad(ib,k)
         enddo
         r1=0.
!         r2=0.
         r3=0.
         do ia=1,nka
         do jt=2,nkt
            x0=rq(jt,ia)
            x1=(rw(jt,ia)-rw(jt-1,ia))
            r1=r1+ff(jt,ia,k)*x0*x1
!            r2=r2+ff(jt,ia,k)*x0**2*x1
            r3=r3+ff(jt,ia,k)*x0**3*x1
         enddo
         enddo
         i0(10,k)=r1
!         i0(11,k)=r2
         i0(12,k)=r3
         i0(11,k)=totrad(1,k)
      enddo
! only for 1D:
      if (n_bln.gt.3) then
         do j=1,12
            i0(j,n-1)=etw(lcl)
            i0(j,n)=etw(lct)
         enddo
         i0(1,n-2)=fnseb
         i0(1,n-3)=flgeg
      endif
      fname='pr .out'
      fname(3:3)=fogtype
 3000 continue
      open (14, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
      write (14) lday,lst,lmin,i0
      close (14)

      end subroutine ploutr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ploutt (fogtype,n_bln)
! output of turbulence data
      implicit double precision (a-h,o-z)
      parameter (nf=100,n=nf+50)
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb42a/ tkeps(n),tkepb(n),tkepd(n)
      common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
      common /cb45/ u(n),v(n),w(n)
      double precision i0
      dimension i0(12,n_bln)
      character *10 fname
      character *1 fogtype
      do k=1,n_bln
         i0(1,k)=atkm(k)
         i0(2,k)=u(k)
         i0(3,k)=v(k)
         i0(4,k)=w(k)
         i0(5,k)=tke(k)
         i0(6,k)=tkeps(k)
         i0(7,k)=tkepb(k)
         i0(8,k)=tkepd(k)
         i0(9,k)=-atkh(k)*buoy(k)
         i0(10,k)=buoy(k)
         i0(11,k)=xl(k)
      enddo
! only for 1D:
      if (n_bln.gt.3) then
         do j=1,12
            i0(j,n-1)=etw(lcl)
            i0(j,n)=etw(lct)
         enddo
      endif
      fname='pt .out'
      fname(3:3)=fogtype
 3000 continue
      open (18, file=fname,status='old',form='unformatted',
     & position='append',err=3000)
      write (18) lday,lst,lmin,i0
      close (18)

      end subroutine ploutt
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine constm (chem,mic,rst)

      USE constants, ONLY :
! Imported Parameters:
     &     cp              ! Specific heat of dry air, in J/(kg.K)


      implicit double precision (a-h,o-z)
! output of constants and parameters used in the current run
      parameter (nf=100,n=nf+50,nb=20,nkt=70,nka=70,n1=n+12,n2=n1-1)
      logical chem,mic,rst
      common /band_o3/ scaleo3_m
      common /blck06/ kw(nka),ka
      common /cb18/ alat,declin                ! for the SZA calculation
      double precision alat,declin

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb44/ r0,r1,g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhow,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb51/ dlgew,dlgenw,dlne
      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)

      dimension xsum(n)
      write (26,6000)
 6000 format (16x,'constants and parameters of the current run'
     & ,///,6x,'numerical grid',/,6x,'eta:')
      write (26,6010) (eta(k),k=1,n)
 6010 format (1x,13e12.5)
      write (26,6020)
 6020 format (6x,'etw:')
      write (26,6010) (etw(k),k=1,n)
      write (26,6030)
 6030 format (6x,'deta:')
      write (26,6010) (deta(k),k=1,n)
      write (26,6040)
 6040 format (6x,'detw:')
      write (26,6010) (detw(k),k=1,n)
      write (26,6050)
 6050 format (6x,'enw:')
      write (26,6010) (enw(k),k=1,nka)
      write (26,6060)
 6060 format (6x,'en:')
      write (26,6010) (en(k),k=1,nka)
      write (26,6070)
 6070 format (6x,'rn:')
      write (26,6010) (rn(k),k=1,nka)
      write (26,6090)
 6090 format (6x,'ew:')
      write (26,6010) (ew(k),k=1,nkt)
      write (26,6100)
 6100 format (6x,'e:')
      write (26,6010) (e(k),k=1,nkt)
      write (26,6110)
 6110 format (6x,'dew:')
      write (26,6010) (dew(k),k=1,nkt)
      write (26,6120)
 6120 format (6x,'rq(k,1):')
      write (26,6010) (rq(k,1),k=1,nkt)
      write (26,6130)
 6130 format (6x,'rq(k,nka):')
      write (26,6010) (rq(k,nka),k=1,nkt)
      write (26,6140)
 6140 format (6x,'zb:')
      write (26,6010) (zb(k),k=1,nb)
      write (26,6150)
 6150 format (6x,'dzb:')
      write (26,6010) (dzb(k),k=1,nb)
      write (26,6160)
 6160 format (6x,'dzbw:')
      write (26,6010) (dzbw(k),k=1,nb)
      write (26,6170)
 6170 format (//,16x,'constants and parameters'/)
      write (26,6173)
 6173 format (6x,'aerosol type: 1=urban 2=rural 3=ocean 4=background')
      write (26,6176) (nar(k),k=1,n)
 6176 format (6x,60i2)
      write (26,6177)
 6177 format (6x,'solution term b0m of droplet growth equation')
      write (26,6010) (b0m(k),k=1,nka)
      write (26,6180)
 6180 format (6x,'r0,r1,g,cp,a0m,ug,vg,dlgew,dlgenw,dlne')
      write (26,6010) r0,r1,g,cp,a0m,ug,vg,dlgew,dlgenw,dlne
      write (26,6190)
 6190 format (6x,'z0')
      write (26,6010) z0
      write (26,6195) ka,kw
 6195 format (//,16x,'ka and kw for aqueous phase reactions',21i6)
      write (26,6200)
 6200 format (//,16x,'dimensions of arrays: n,n1,n2,nb,nka,nkt')
      write (26,6210) n,n1,n2,nb,nka,nkt
 6210 format (/,1x,6i10)
      write (26,6220) alat,declin,wmin,wmax,tw,scaleo3_m
 6220 format (//,6x,'geogr. latitude ',f9.1,' declination ',f9.1,
     &' large scale subsidence in m/s',2f9.5,
     &' water temperature',f8.2,' O3 column (hv only) ',f4.0,//)
      write (26,6230) chem,mic,rst
 6230 format (6x,'current program evaluation: ','   chem: ',l1,
     &' mic: ',l1,'   rst: ',l1,//)
      xxsum=0.
      do k=2,nf
         xsum(k)=0.
         do ia=1,nka
            do jt=1,nkt
               xsum(k)=xsum(k)+ff(jt,ia,k)*en(ia)
            enddo
         enddo
         xsum(k)=xsum(k)*1.e+09
         xxsum=xxsum+xsum(k)*detw(k)
      enddo
      write (26,6240) 
 6240 format (/,6x,'aerosol mass in ug m**-3 in layers 2 - nf')
      write (26,6250) xsum
 6250 format (1x,15f8.3)
      write (26,6260) xxsum
 6260 format(6x,'total aerosol mass in ug m**-2 of layers 2 - nf',f12.3)

      end subroutine constm
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine constc
!
! Description:
!    output of chemical constants and parameters used in current run
!    This subroutine is called only once during initialisation (no restart case)
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      11/2016   Cleaned + header + implicit none         <Josue Bock>
!                    Removed outdated features
!
! 1.1       ?        Added advected species in the output     <Roland von Glasow>
!
! 1.0       ?        Original code.                           <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!

! Declarations:
! Modules used:
      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
     &     j5

      USE global_params, ONLY :
! Imported Parameters:
     &     j2

      implicit none

! Local scalars:
      integer j

! Common blocks:
      common /kpp_eul/ xadv(10),nspec(10),neula
      double precision xadv
      integer nspec, neula

!- End of header ---------------------------------------------------------------

      write (60,5900) neula
      write (60,5910) (nspec(j),xadv(j),j=1,10)
 5900 format ('euler (=0) or lagrangean view (=1): ',i3,' the following'
     &     ,' species are advected only if neula=0')
 5910 format (i3,d12.3)

      write (60,6090)
 6090 format (//,16x,'dimension of arrays: j1, j2, j5')
      write (60,6100) j1,j2,j5
 6100 format (/,1x,5i10,//)

      end subroutine constc
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine profm (dt)

! output of meteorological profiles

! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     n1,
     &     nb,
     &     nka,
     &     nkt,
     &     mbs

      implicit double precision (a-h,o-z)

      common /cb16/ u0,albedo(mbs),thk(n1)
      double precision u0, albedo, thk

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
      common /cb44/ r0,r1,g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhow,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      common /cb45/ u(n),v(n),w(n)
      common /cb46/ ustern,gclu,gclt
      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb48/ sk,sl,dtrad(n),dtcon(n)
      double precision sk, sl, dtrad, dtcon

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)

      character *10 srname
      dimension xsum(n)
      srname='          '
      write (26,6000) it,dt,lday,lst,lmin
 6000 format (//,6x,i8,'-th. timestep dt = ',f4.1,' sec ',i2,' day ',
     & i2,' hour ',i2,' min '/)
      write (26,6010)
 6010 format (1x,'k  height  d/dz(thetl)*100  atkm  atkh  xl',
     &'  tke  tkep  thetl   d/dt(dtrad)*3600   d/dt(dtcon)*3600')
      write (26,6020) (k,etw(k),(thetl(k+1)-thetl(k))/deta(k)*100.,
     & atkm(k),atkh(k),xl(k),
     & tke(k),tkep(k),thetl(k),dtrad(k)*3600.,dtcon(k)*3600.,k=n-1,1,-1)
 6020 format (1x,i3,f7.1,3f9.3,f7.1,f9.4,e12.4,f9.3,2f12.4)
      write (26,6030)
 6030 format (/,1x,'k  height       u         v         t          p',
     & '      theta      feu          q        m2        fsum')
      write (26,6040) (k,eta(k),u(k),v(k),
     & t(k)-273.15,p(k)/100.,theta(k),
     & feu(k)*100.,1000.*xm1(k),1000.*xm2(k),
     & fsum(k),k=n,1,-1)
      xxm1=0.
      xxm2=0.
      do k=1,n
!         xxm1=xm1(k)*detw(k)+xxm1    
         xxm1=xm1(k)*detw(k)*rho(k)*1000+xxm1
         xxm2=xm2(k)*detw(k)*1000+xxm2
! xxm1 in g/m**2 vapour content of atm., xxm2 in g/m**2 liquid water content of atm.
      enddo
 6040 format (1x,i3,f10.1,9f10.3)
      write (26,6050) xxm1,xxm2,tau*1000.,reif*1000.,trdep*1000.
     & ,ds1*1000.,ds2*1000.
 6050 format (/,1x,'sum(xm1*detw*rho)',f10.3,1x,'sum(xm2*detw)',f10.3,
     & 1x,'dew',f10.3,5x,'rime',f10.3,5x,'particles',f10.3,5x,'aerosol',
     & f10.3,5x,'droplets',f10.3)
      write (26,6060) u0,sk,sl,-5.6697d-8*t(1)**4
 6060 format (1x,'surface radiative fluxes'/,10x,'u0: ',f10.4,
     & 10x,'solar: ',e11.4,3x,'infrared: ',e11.4,3x,'emission: ',e11.4)
      write (26,6070) z0,ustern,ajq,ajs,ajm,ajd
 6070 format (1x,'roughness length ',e10.3,' friction velocity u*',
     & f8.3,/,1x,'surface moisture fluxes'/,
     & 10x,'water vapor: ',e11.4,3x,'droplet sedimentation: ',e11.4,3x,
     & 'ground moisture: ',e11.4,3x,'dew storage: ',e11.4)
      write (26,6080) ajb,ajl,ajt,sk+sl-5.669d-8*t(1)**4
 6080 format (1x,'surface heat fluxes'/,
     & 10x,'ground heat: ',e11.4,3x,'latent heat: ',e11.4,3x,
     & 'sensible heat: ',e11.4,3x,'net radiation: ',e11.4)
      write (26,6090)
 6090 format (/,1x,'temperature and volumetric moisture content in',
     &' ground:')
      write (26,6100) (zb(k),k=1,nb)
      write (26,6110) (tb(k)-273.15,k=1,nb)
      write (26,6120) (eb(k),k=1,nb)
 6100 format (4x,'zb:',/,10f10.3,/,10f10.3)
 6110 format (4x,'tb:',/,10f10.3,/,10f10.3)
 6120 format (4x,'eb:',/,10f10.3,/,10f10.3)
      xxsum=0.
      do k=2,nf
         xsum(k)=0.
         do ia=1,nka
            do jt=1,nkt
               xsum(k)=xsum(k)+ff(jt,ia,k)*en(ia)
            enddo
         enddo
         xsum(k)=xsum(k)*1.e+09
         xxsum=xxsum+xsum(k)*detw(k)
      enddo
      write (26,6240) 
 6240 format (/,6x,'aerosol mass in ug m**-3 in layers 2 - nf')
      write (26,6250) xsum
 6250 format (1x,15f8.3)
      write (26,6260) xxsum
 6260 format(6x,'total aerosol mass in ug m**-2 of layers 2 - nf',f12.3)

! 141  format (2i3,9d12.4)
! 142  format (6x,9d12.4)

      call ion_mass (srname)
      
!      do k=2,nf
!         do kc=1,nkc_l
!            write (*,141) (k,kc,(dss(k,l,kc),l=1,lsp))
!            write (*,142) (svc(k,kc,kkc),kkc=1,nkc_l)
!c        write (*,142) (fss(k,kc,1),fss(k,kc,2),svc(k,kc,1),svc(k,kc,2))
!         enddo
!      enddo

      end subroutine profm
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine profc (dt,mic)

! vertical profiles of chemical species

! technical comment:
!     when writing the blocks headers ('height' + 'species_names') the format
!     can include 10a12 even if less than 10 species names are written.
!     However, when writting the arrays (one line per model level), the exact
!     number of species has to be used in the format descriptor (fmt), because
!     it is repeated for each line.

      USE gas_common, ONLY :
! Imported parameter
     &     j1,
     &     j5,
! Imported Array Variables with intent (in):
     &     s1,
     &     s3,
     &     gas_name,
     &     rad_name

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j3,
     &     j6,
     &     nkc,
     &     nlev,
     &     nrxn,
     &     nf,
     &     n

      implicit none

      double precision dt
      logical mic

      double precision si(10,n),xfac(nf)

      character (len=17) :: fmt
      integer :: k, kc, l
      integer :: jblock, jlay, jspec, jspec_max
      integer :: nrate
      double precision :: zx0(n) ! conversion factor

      common /blck01/ am3(n),cm3(n)
      double precision am3, cm3

      common /blck11/ rc(nkc,n)
      double precision rc

      common /blck12/ cw(nkc,n),cm(nkc,n)
      double precision cw, cm

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1

      common /budg/ bg(2,nrxn,nlev),il(nlev) ! jjb corrected, was written il(nrxn)
      double precision bg
      integer il

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /liq_pl/ nkc_l
      integer nkc_l


! ==============================================================================
!  -- 1.1 -- Gas phase species: header and common conversion factor
! ==============================================================================

!     Write header
!     ------------
      write (60,6000) it,dt,lday,lst,lmin
 6000 format (//,6x,i8,'-th. timestep dt = ',f4.1,' sec ',i2,' day ',
     & i2,' hour ',i2,' min '/)
      write (60,6010)
 6010 format (10x,'gas phase species in ppb; at the ground total',
     &' deposition in molecules/cm**2')


!     output in ppb
!     define conversion factor (zx0):
      do jlay=2,n
         zx0(jlay) = 1.e+09/am3(jlay)
      end do



! ==============================================================================
!  -- 1.2 -- Non radical, gas phase species
! ==============================================================================

!     Write data blocks of 10 species
!     -------------------------------
      jblock = 0
      do while (10*jblock < j1)
         jspec_max = MIN(10,j1-10*jblock) ! Number of output species in the current block (10 or less)

         ! No conversion for layer 1
         do jspec=1,jspec_max
            si(jspec,1)=s1(10*jblock+jspec,1)
         enddo

         do jlay=2,n
            do jspec=1,jspec_max
              si(jspec,jlay)=s1(10*jblock+jspec,jlay)*zx0(jlay)
            end do
         end do

         ! write gas names
         write (60,6015)'height',
     &      (TRIM(gas_name(10*jblock+jspec)),jspec=1,jspec_max)
 6015    format (3x,a7,10a12)
         ! write gas concentrations (and deposition at the ground)
         write (fmt,'( "(f10.1,",i2,"e12.5)" )')jspec_max
         write (60,fmt) (eta(jlay),
     &                   (si(jspec,jlay),jspec=1,jspec_max),jlay=n,1,-1)

         jblock = jblock + 1
      end do



! ==============================================================================
!  -- 1.3 -- Radical, gas phase species
! ==============================================================================

!     Write data blocks of 10 species
!     -------------------------------
      jblock = 0
      do while (10*jblock < j5)
         jspec_max = MIN(10,j5-10*jblock) ! Number of output species in the current block (10 or less)

         ! no deposition for radicals, start at layer 2
         do jlay=2,n
            do jspec=1,jspec_max
              si(jspec,jlay)=s3(10*jblock+jspec,jlay)*zx0(jlay)
            end do
         end do

         ! write radicals names
         write (60,6015)'height',
     &      (TRIM(rad_name(10*jblock+jspec)),jspec=1,jspec_max)
         ! write radicals concentrations
         write (fmt,'( "(f10.1,",i2,"e12.5)" )')jspec_max
         write (60,fmt) (eta(jlay),
     &                   (si(jspec,jlay),jspec=1,jspec_max),jlay=n,2,-1)

         jblock = jblock + 1
      end do



! ==============================================================================
!  -- 2345 -- to be done
! ==============================================================================

!     output of reaction rates integrated (later: over 1 hour), converted to mol/(mol*h)
      write (60,6170)
         write (60,6180) eta(il(1)),eta(il(2)),eta(il(3))
         nrate=600
         if (nkc_l.eq.4) nrate=1120
         do l=1,nrate
           write (60,6190) l,bg(2,l,1)/am3(il(1)),bg(1,l,1)/am3(il(1)),
     &           bg(2,l,2)/am3(il(2)),bg(1,l,2)/am3(il(2))
     &           ,bg(2,l,3)/am3(il(3)),bg(1,l,3)/am3(il(3))
         enddo
! 6170 format (//,'reaction rates integrated over 1 hour, converted to' 
 6170 format (//,'accumulated reaction rates [mol/(m^3(air) s)]')
 6180 format (/,'height = ',19x,f10.2,' m',4x,f10.2,' m',4x,f10.2,' m')

! 6190 format ('reaction no. ',i4,' rate : ',4d16.4)
 6190 format ('no. ',i4,' : ',6d16.8)



!      if (.not.mic.or.lct.le.1) return
      if (.not.mic) return
! liquid phase chemistry
      do 1033 kc=1,nkc_l
      write (60,6040) kc
 6040 format (/,10x,'aqueous phase species in mole/m**3;',
     &' at the ground total deposition in mole/m**2 for bin:',i3,//,
     & 4x,'height',5x,'cw',7x,'rc',7x,'hno3',6x,'nh3',7x,'so2',5x,
     & 'h2so4',7x,'o3',9x,'h2o2',7x,'fe(3)',6x,'mn(2)')
      write (60,6020) (eta(k),cw(kc,k),rc(kc,k),sl1(3,kc,k),
     & sl1(4,kc,k),sl1(5,kc,k),sl1(6,kc,k),sl1(7,kc,k),
     & sl1(19,kc,k),sl1(44,kc,k),sl1(45,kc,k),k=nf,1,-1)
!     & sl1(19,kc,k),sl1(44,kc,k),sl1(45,kc,k),k=lct,lcl,-1)
      write (60,6041)
 6041 format (4x,'height',5x,'OH',7x,'HO2',7x,'NO3',6x,'NO',7x,'NO2',5x,
     & 'HONO',7x,'HCHO')
      write (60,6042) (eta(k),sl1(j2-j3+2,kc,k),sl1(j2-j3+3,kc,k),
     & sl1(36,kc,k),sl1(1,kc,k),sl1(2,kc,k),sl1(21,kc,k),sl1(17,kc,k)
     & ,k=nf,1,-1)
!     6 ,k=lct,lcl,-1)
 6042 format (f10.1,7e10.3)
      write (60,6050) kc
! 6050 format (/,20x,'ion concentrations in mole/liter',//,
 6050 format (/,20x,'ion concentrations in mol/m^3 for bin:',i3,//,
!     & 4x,'height',7x,'h+',7x,'nh4+',6x,'cl-',2x,'ch2ohso3-',2x,
     & 4x,'height',7x,'h+',7x,'nh4+',6x,'cl-',2x,'Br-',2x,
!     & 3x,'hso3-',6x,'so3=',6x,'so4-',5x,'so4=',5x,'no3-',6x,'fe')
     & 3x,'hso3-',6x,'so3=',6x,'so4-',5x,'so4=',5x,'no3-',6x,'tracer')
      do k=1,nf
         xfac(k)=1.             !mol/m^3_air       --> mol/m^3_air
      enddo
!      write (60,6020) (eta(k),sion1(1,kc,k)*xfac(k),sion1(2,kc,k)*
      write (60,6020) (eta(k),sion1(1,kc,k)*xfac(k),sion1(2,kc,k)*
     & xfac(k),sion1(14,kc,k)*xfac(k),sion1(24,kc,k)*xfac(k),
     & (sion1(l,kc,k)*xfac(k),l=5,8),
!     & sion1(17,kc,k)*xfac(k),k=lct,lcl,-1)
!     & sion1(13,kc,k)*xfac(k),sion1(3,kc,k)*xfac(k),k=lct,lcl,-1)
     & sion1(13,kc,k)*xfac(k),sion1(21,kc,k),k=nf,1,-1)
 6020 format (f10.1,10e10.3)
 1033 continue
      
      write (60,*) 'done with profc'

      end subroutine profc
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine profr

      USE global_params, ONLY :
! Imported Parameters:
     &     n1,
     &     n4,
     &     mb,
     &     mbs,
     &     nrfl

      implicit double precision (a-h,o-z)

      common /cb02/ t(n4),p(n4),rho(n4),xm1(n4),rho2(n1),frac(n1),
     & ts,ntypa(n1),ntypd(n1)
      double precision t,p,rho,xm1,rho2,frac,ts
      integer ntypa,ntypd

      common /cb11/ totrad (mb,n1)
      double precision totrad

      common /cb15/ fnseb,flgeg,hr(n1)
      double precision fnseb, flgeg, hr

      common /cb16/ u0,albedo(mbs),thk(n1)
      double precision u0, albedo, thk

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /kurz/ fs1(n4),fs2(n4),totds(n4),ss(n4),fsn(n4),dtdts(n1)
      double precision fs1, fs2, totds, ss, fsn, dtdts

      common /lang/ fl1(n4),fl2(n4),fln(n4),dtdtl(n1)
      double precision fl1, fl2, fln, dtdtl

      write (40,6000) lday,lst,lmin,u0
 6000 format (10x,'day: ',i5,10x,'hour: ',i5,10x,'minute: ',i5,
     & 10x,'cosine of zenith distance: ',f8.2,/)
      write (40,6010)
! 6020 format (7f14.3,e14.5,f14.1)
! 6030 format (/,10x,' fs1, fs2, ss, f1l, fl2, f1f, f2f, p') ! jjb
 6030 format (/,10x,' fs1, fs2, ss, f1l, fl2')
! 6040 format (7f14.3,f14.1)
 6010 format (10x,' totrad(l,i) l=1,6, fn, hr, p')
      do i=1,nrfl
         write (40,6021) (totrad(ib,i),ib=1,mbs),hr(i),p(i)
      enddo
      write (40,6030)
      write (40,6030)
      do i=1,nrfl
         write (40,6021) (totrad(ib,i),ib=mbs+1,mb)
      enddo
      write (40,6030)
 6021 format(12f10.3)
 6041 format(6f14.3)
      do i=1,nrfl
         write (40,6041) fs1(i),fs2(i),ss(i),fl1(i),fl2(i)
      enddo

      end subroutine profr
