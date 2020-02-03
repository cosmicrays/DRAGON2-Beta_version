 

cc#include "config.h"

      subroutine dmspec(mxin,modein,sigmavin,egev,yieldk,
     & spectrum, nE)
      implicit none
c     
c      real*8 mxin,sigmavin,dshaloyield,egev
      integer modein,yieldk,istat,i, nE
      double precision mxin,sigmavin,dshaloyield,egev(nE),spectrum(nE)
      real*8 jpsidelta,egath,fluxgac,dshrgacont,dshrgacsusy
      integer, save :: counter = 0

c      STATIC counter      
#ifdef HAVE_DS
c      include 'dssusy.h'
      include 'dsge.h'
      include 'dsdirver.h'
      include 'dsmssm.h'
      include 'dsidtag.h'
      include 'dshacom.h'
      include 'dshmcom.h'
      include 'dsprep.h'
     
      if (counter.eq.0) call dsinit
       counter=1

c      prtlevel=0
c      dmspec = 0.
   
ccc   tipo higgsino puro
ccc   set the WIMP mass:
ccc   mxin=1000.d0 ! GeV
ccc   select one final state in WIMP pair annihilations  
ccc   modein=13 ! w+w-
ccc   initialize them into the code (see subroutine below):
ccc   sigmavin=3.d-26 ! in cm^3 / s
      call setoneannmod(mxin,sigmavin,modein)
c      write(*,*) mxin
c      write(*,*) sigmavin
c      write(*,*) modein
      
ccc      yieldk =  151             ! positrons
      do i=1,nE
         spectrum(i) = dshaloyield(egev(i),yieldk,istat) ! yield per annihil.
c        dmspec=dshaloyield(egev,yieldk,istat) ! yield per annihil.
c        write(*,*) egev(i),spectrum(i)
      enddo
      return 
      end
      
      subroutine setoneannmod(mxin,sigmavin,modein)
c     
c mxin is the value (in GeV) you assign to the WIMP mass
c
c sigmavin is the value (in cm^3 / s) you assign to the WIMP pair 
c annihilation rate
c
c modein selects one (only one in this simple setup, assigning a 100% 
c branching ratio to it) among the mode list (the list includes the
c final states which can be non-zero in the pair annihilation of
c neutralino, at zero temperature and within the MSSM; the higgs boson
c channels 5 to 11 can be used only if you define branching ratios for 
c their decay chains; this is not implemented at this stage): 
c
c Ch No  Particles       Old Ch No 
c  5     h10 h30         7
c  6     h20 h30         11
c  8     z0 h10          8
c  9     z0 h20          9
c  11    w+ h- / w- h+   10
c  12    z0 z0           6
c  13    w+ w-           5
c  17    mu+ mu-         13
c  19    tau+ tau-       4
c  22    cc-bar          1
c  24    tt-bar          3
c  25    bb-bar          2
c  26    gluon gluon     12
c  29    z gamma         14
c
      implicit none
      include 'dshacom.h'
      include 'dsprep.h'
      integer i
      real*8 mxin,sigmavin
      integer modein
c
      do i=1,29
        habr(i) = 0.d0
      enddo
      habr(modein)=1.d0
 
c no additional yield from internal bremsstrahlung, but keep the 
c model-independent final state radiation contribution that is 
c automatically included in the tabulated PYTHIA results)
      haib='none'

c WIMP mass for halo yield routines
      hamwimp=mxin ! GeV

      dshasetupcalled=.true.

c WIMP parameters for halo rates routines: 
      mx=mxin ! GeV
      sigmav = sigmavin ! in cm^3 / s
      do i=1,29
        sigv(i) = 0.d0
      enddo
      newmodelsigmav=.false.
      dsprepcalled=.true.
c
c#else
c      dmspec = 0
#endif
      return
      end

