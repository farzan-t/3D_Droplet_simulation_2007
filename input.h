       implicit none
c*******************************************************************
c     This section defines the pointers and the memory space needed
c     by the multigrid pressure subroutines
c***
      INTEGER NGM,NN3,NN2,NN1,MEMLEN,mem
  
c***
c     memory necessary for a cube of side: 2**NGM+2
c***
      PARAMETER (NGM=6)
      PARAMETER (NN3=(2**(3*NGM)-1)/7,NN2=2**(2*NGM),NN1=2**NGM)
      PARAMETER (MEMLEN=64*NN3 + 64*NN2 + 192*NN1 + 72*NGM)
c***
c     memory for a box : 130x66x66: 5246832
c     memory for a box : 258x130x130: 40128440
c***
c******************************************************************
      INTEGER icc(NGM),iu(NGM),ires(NGM),irhs(NGM),ia(NGM),ib(NGM),
     %         ic(NGM),iei(NGM),ng
      DOUBLE PRECISION z(MEMLEN)
      COMMON /pointers/ icc,iu,ires,irhs,ia,ib,ic,iei,ng
      COMMON /memory/ z, mem
c******************************************************************
      INTEGER NGMX,NMX,NGMY,NMY,NGMZ,NMZ
      PARAMETER (NGMX=6,NMX=2**NGMX+2)
      PARAMETER (NGMY=5,NMY=2**NGMY+2)
      PARAMETER (NGMZ=5,NMZ=2**NGMZ+2)
      DOUBLE PRECISION u(NMX,NMY,NMZ),v(NMX,NMY,NMZ),w(NMX,NMY,NMZ) 
      DOUBLE PRECISION S11(NMX,NMY,NMZ),S22(NMX,NMY,NMZ),
     %                 S33(NMX,NMY,NMZ),S12(NMX,NMY,NMZ),
     %                 S13(NMX,NMY,NMZ),S23(NMX,NMY,NMZ)
      common /surfcom/ u,v,w,S11,S22,S33,S12,S13,S23
c**********************************************************************
      DOUBLE PRECISION mut(NMX,NMY,NMZ),lmix(NMX),cmu(NMX,NMY,NMZ),
     & turbst(NMX,NMY,NMZ),ktur(NMX,NMY,NMZ),eps(NMX,NMY,NMZ),
     & u2p(NMX,NMY,NMZ),v2p(NMX,NMY,NMZ),w2p(NMX,NMY,NMZ),
     & uvp(NMX,NMY,NMZ),uwp(NMX,NMY,NMZ),vwp(NMX,NMY,NMZ),
     & hx(NMX,NMY,NMZ),hy(NMX,NMY,NMZ),hz(NMX,NMY,NMZ),
     & ax(NMX,NMY,NMZ),ay(NMX,NMY,NMZ),az(NMX,NMY,NMZ),
     & bx(NMX,NMY,NMZ),by(NMX,NMY,NMZ),bz(NMX,NMY,NMZ)
      common /turbulence/mut,lmix,turbst,ktur,eps,cmu,
     & u2p,v2p,w2p,uvp,uwp,vwp,hx,hy,hz,
     & ax,ay,az,bx,by,bz

c***********************************************************************
      INTEGER nx,ny,nz
      common /dimension/ nx,ny,nz
c***********************************************************************
c     common block containing the constants controlling different
c     algorithms
c***
      DOUBLE PRECISION ANGLEDIF,MAXCFL,OMEGA,EPSC,EPSURR,tauh
      INTEGER CYCLES,NPRE,NPOST,NTOP,restart,icase,ijet,idrop,
     % tauh2,ibubble,isheet,iturb,nfilter,ibc,tswap,convectype,
     % hswitch,surfacet
      COMMON /DIFF/ ANGLEDIF, EPSC, EPSURR,tauh,tauh2
      COMMON /CYCLING/MAXCFL,OMEGA,CYCLES,NPRE,NPOST,NTOP,
     % restart,icase,ijet,idrop,ibubble,isheet,iturb
     & ,nfilter,ibc,tswap,convectype,hswitch,surfacet

c***********************************************************************
c     common block containing the initial and final time step, the
c     storing,checking and dumping rates and the names of the I/O files
c***
      INTEGER t0,tmax,tstore,ntecplot,tsave,tprint
      COMMON /TIMES/ t0, tmax,tstore,ntecplot,tsave,tprint
      
c***********************************************************************
c     common blocks containing the dimensional parameters
c***
      DOUBLE PRECISION tau_SI,timescale,lengthscale,length_SI, 
     %     rho1_SI,rho2_SI,visc1_SI,visc2_SI,sigma_SI,grav_SI, 
     %     lengthsize,numer,grav1_SI,xvel1_SI,yvel1_SI,zvel1_SI
     %     ,xvel2_SI,yvel2_SI,zvel2_SI,radius_si,z0_SI,y0_SI
      COMMON /dinem/tau_SI,timescale,lengthscale,length_SI, 
     %     rho1_SI,rho2_SI,visc1_SI,visc2_SI,sigma_SI,grav_SI, 
     %     lengthsize,numer,grav1_SI,xvel1_SI,yvel1_SI,zvel1_SI
     %     ,xvel2_SI,yvel2_SI,zvel2_SI,radius_si,z0_SI,y0_SI
c***********************************************************************
c     common blocks containing the non-dimensional physical constants
c***
      DOUBLE PRECISION tau,time,radius,h,rho1,rho2,mu1,mu2,sigma,g,g1,
     % y0,z0,xvel1,yvel1,zvel1,xvel2,yvel2,zvel2,rhon(NMX,NMY,NMZ),
     % mun(NMX,NMY,NMZ)
      COMMON /PHYSICS/ tau,time,radius,h,rho1,rho2,mu1,mu2,sigma,g,g1, 
     % y0,z0,xvel1,yvel1,zvel1,xvel2,yvel2,zvel2,rhon,mun

c***********************NEW*******************************************

      DOUBLE PRECISION hxmax(NMX,NMY,NMZ),hxmin(NMX,NMY,NMZ),
     & hymax(NMX,NMY,NMZ),hymin(NMX,NMY,NMZ),
     & hzmax(NMX,NMY,NMZ),hzmin(NMX,NMY,NMZ)

       COMMON /NEWW/ hxmax,hxmin,hymax,hymin,hzmax,hzmin

c***********************************************************************
