c***********************************************************************
      PROGRAM THREE DIM
c***
      INCLUDE 'input.h'
      INTEGER t
c***
c     SUBROUTINES: initialize,printout,VOFmethod,surfacetens,solvelpres
c***
c     Open files foe input-output data
c***
       call openfiles()
c***
c    Read input data and creat module
c***
      call initialize()
c***
c    MAIN DO LOOP
c***
      do t = t0,tmax
c***
c    Output some fields
c***
        call printout(t,z(icc(ng)),z(iu(ng)))
c***
c    VOF calculation
c***
      call VOFmethod(t)


      if(hswitch.eq.1) call CalcH3d(z(icc(ng)))

c***
c    eddy viscosity calculation
c***
	if (iturb.eq.1.or.iturb.eq.2) call eddyvisc (z(icc(ng)),z(iu(ng)))
c***
c    k-epsilon turbulence model
c***
	if (iturb.eq.3.or.iturb.eq.4) 
     &	call kepsmodel (t,z(iu(ng)),z(ires(ng)))
c***
      call checktime(t)
c***
c     get surface tension contribution to the stress tensor
c***
c      call surfacetens(z(icc(ng)),z(iu(ng)),z(ires(ng)),z(irhs(ng))
c     %        ,z(ia(ng)),z(ib(ng)),z(ic(ng)),z(iei(ng)))	
	if (surfacet.eq.1) then
      call surfacetens(z(icc(ng)),z(iu(ng)),z(ires(ng)),z(irhs(ng))
     %        ,z(ia(ng)),z(ib(ng)),z(ic(ng)),z(iei(ng)))
	else if (surfacet.eq.2) then
	call csf(z(icc(ng)),z(iu(ng)),z(ires(ng)),z(irhs(ng)))
	else if (surfacet.eq.3) then
	call csfbkz(z(icc(ng)),z(iu(ng)),z(ires(ng)),z(irhs(ng)))
      else if (surfacet.eq.4) then
	call meier(z(icc(ng)),z(iu(ng)),z(ires(ng)),z(irhs(ng)))
      else if (surfacet.eq.5) then
	call SGIP(z(icc(ng)),z(iu(ng)),z(ires(ng)),z(irhs(ng)))
      else if (surfacet.eq.6) then
	call SGIP2(z(icc(ng)),z(iu(ng)),z(ires(ng)),z(irhs(ng)))

      else
      STOP 'capalg'
	endif
c***
c     N-S + multigrid
c***
      call solvelpres(z(icc(ng)),z(iu(ng)),z(ires(ng)),z(irhs(ng))
     %        ,z(ia(ng)),z(ib(ng)),z(ic(ng)),z(iei(ng)),t)

c		if (mod(t,10).eq.0 .and. t.gt.50) call newcc(z(icc(ng)))

      enddo
c***
      stop
      end
c **************	test new idea for new c	*********
	subroutine newcc(cc)
	include 'input.h'

	double precision cc(nx,ny,nz),ccf(nx,ny,nz),temp(nx,ny,nz)
	integer i,j,k,nf
	nf=3
	write(23,*) 'new c '
	call copy(ccf,cc,nx,ny,nz)
	do i=1,nf
	call filterc27(ccf,temp,nx,ny,nz)
	call copy (ccf,temp,nx,ny,nz)
	call bc_out(ccf,nx,ny,nz)
	enddo
c	call filterc27(cc,ccf,nx,ny,nz)
	do i=1,nx
	do j=1,ny
	do k=1,nz
		if (cc(i,j,k).gt.0.001 .and. cc(i,j,k).lt.0.999) then
		cc(i,j,k)=ccf(i,j,k)
		
		endif
	enddo
	enddo
	enddo
			
	end
c***********************************************************************
      SUBROUTINE openfiles () 
	OPEN(1,FILE='inputdata.in',STATUS='OLD')
	OPEN(23,FILE='outputdata.out')
	open(24,file='spurious.plt')
	open(26,file='normv.plt')
	open(27,file='hysteresis.plt')
	open(28,file='velocity.plt')
	
			write(24,*) 'VARIABLES="time","Vmax"'
			write(26,*) 'VARIABLES="time","norm V"'
	write(27,*) 'VARIABLES="time","thetaH","theta A","theta R"'
			write(28,*) 'VARIABLES="time (ms)","u-center"'
	return
	end

c***********************************************************************
      SUBROUTINE initialize ()
c***
      INCLUDE 'input.h'
      DOUBLE PRECISION vscale
      DOUBLE PRECISION FINDNINF
      EXTERNAL FINDNINF

c***
c     SUBROUTINES: readinput,  writeminmax, nondim,
c          nondimnumbers, inijet,inidrop,inisheet,iniget,bc_v,
c             bc_c,memoryaloc, printdata
c   read input data
c***
      call readinput ()
c***
c     check max dimension of MEMLEN and initialize indices
c***
      call memoryaloc()
c***
c   non-dimensionalze the input data and calculate the non-dim numbers
c***
      call nondim(vscale)
c***
      call nondimnumbers(vscale)
c***  
      if (t0 .eq. 0) then
        write (6,*) 'initial velocity max  = ',
     %        FINDNINF(u,v,w,nx,ny,nz)
        write (23,*) 'initial velocity max  = ',
     %        FINDNINF(u,v,w,nx,ny,nz)

c***
c     read IC and make a few checks
c***
         write (*,*) '--- RESTART = ', restart
          if(restart.eq.0) then
	        t0=0
			write (*,*) 't0 =', t0
			write (23,*) 't0 =', t0
			  if(ijet.eq.1) then
	             call inijet(z(icc(ng)),z(iu(ng)))
	          elseif(idrop.eq.1) then
				call smooth2(z(icc(ng)),z(iu(ng)))
		!		call smoothdrop(z(icc(ng)),z(iu(ng)))
	    !         call inidrop(z(icc(ng)),z(iu(ng)))
                elseif(isheet.eq.1) then
	             call inisheet(z(icc(ng)),z(iu(ng)))
                else
                   stop
	          endif
	     else
			write (*,*) '--- RESTARTING WITH FILE dump.dat ---'
			write (23,*) '--- RESTARTING WITH FILE dump.dat ---'
			call iniget(z(icc(ng)))
           endif
         call bc_c(z(icc(ng)))
         call bc_v()
         call calcrho(z(icc(ng)))
      endif 
      call printdata()
c***
      return 
      end
c************************** READ INITIAL DATA ***********************       
      SUBROUTINE readinput ()
      INCLUDE 'input.h'
c      ***** numerical properties *********  
      read (1,*) nx        ! max points in x-direction
      read (1,*) ny        ! max points in y-direction
      read (1,*) nz        ! max points in z-direction
      read (1,*) tmax      ! maximum number of time steps
      read (1,*) nfilter   ! filter volume fraction (def.is zero) 
      read (1,*) cycles    ! max no. of cycles in multigrid method(def.is 60 )
      read (1,*) OMEGA     ! relaxation coefficient(def.is 1.0)
      read (1,*) NPRE      ! pre-relaxation coefficient(def.is 1.0)
      read (1,*) NTOP      ! post-relaxation coefficient(def.is 1.0)
      read (1,*) NPOST     ! top-relaxation coefficient(def.is 1.0)
      read (1,*) restart   ! =0 for new data, =1 for continueing the caculation (def.is 0)
      read (1,*) maxcfl    ! maximum value for CFL(def.is 0.8)
	read (1,*) iturb     ! =0 laminar,=1 mixing length with sigmtau model, =2
     % 					 !  mixing length with H'K' model, =3 k-e model, =4 anisotropic k-e model
      read (1,*) convectype !convection term, 1=upwind, 2=central, otherwise=no conv.(def.is 1)
      read (1,*) hswitch    ! =0 off, =1, on

c      ***** physical and flow properties **********
      read (1,*) xvel1_SI   !initial x-component of velocity of fluid 1
      read (1,*) yvel1_SI   !initial y-component of velocity of fluid 1
      read (1,*) zvel1_SI   !initial z-component of velocity of fluid 1
      read (1,*) xvel2_SI   !initial x-component of velocity of fluid 2
      read (1,*) yvel2_SI   !initial y-component of velocity of fluid 2
      read (1,*) zvel2_SI   !initial z-component of velocity of fluid 2
      read (1,*) rho2_SI    !initial density of fluid 2
      read (1,*) rho1_SI    !initial density of fluid 1
      read (1,*) sigma_SI   !surface tension coefficient
      read (1,*) visc2_SI   !initial viscosity of fluid 2
      read (1,*) visc1_SI   !initial viscosity of fluid 1
      read (1,*) grav_SI    !garvity in z-direction
      read (1,*) grav1_SI   !garvity in x-direction

c      ****** geometry and flow type *******
      read (1,*) icase      ! module (1-jet, 2=drop, 3=bubble)
      read (1,*) ibc        ! bourdary condition type (=1, original, =2, standard) 
	read (1,*) length_SI  ! length of the domain in x-direction
      read (1,*) radius_si  ! module radius (NOTE: MUST BE CORRECTED!!) 
c      read (1,*) x0_SI     ! x location of center of jet, drop of bubble 
	read (1,*) y0_SI      ! y location of center of jet, drop of bubble
 	read (1,*) z0_SI      ! z location of center of jet, drop of bubble
c	(if z0_SI=1000 z0 will be in center of channel automatically)

c     ****** print data ********
      read (1,*) tprint    ! print data on screen
      read (1,*) ntecplot  ! print data for ntecplot
      read (1,*) tsave     ! print unformatted data
	read (1,*) surfacet	!switch to method of surface tension
	close(1)
c***
c     check nx, ny, nz against the maximum values in input.h
c***
      write (*,*)
      if (nx .ne. NMX) write (*,*) 'STOP: nx .ne. NMX:', nx, NMX
      if (ny .ne. NMY) write (*,*) 'STOP: ny .ne. NMY:', ny, NMY
      if (nz .ne. NMZ) write (*,*) 'STOP: nz .ne. NMZ:', nz, NMZ
c***
      if (ntecplot.eq.0) then
         write (*,*) 'STOP: missing all storing rates'
         stop
      endif
     	if (icase.eq. 1) then
	ijet=1
      idrop=0
	ibubble=0
	isheet=0
	elseif(icase.eq.2) then
	ijet=0
	idrop=1
	ibubble=0
	isheet=0
      elseif(icase.eq.3) then	
	ijet=0
	idrop=0
	ibubble=1
	isheet=0
      elseif(icase.eq.4) then	
	ijet=0
	idrop=0
	ibubble=0
      isheet=1
	else
	stop "incorrct module number"
	endif
      return
      end
c***********************************************************************
      SUBROUTINE memoryaloc()
c***
      INCLUDE 'input.h'
      INTEGER MYMALLOC,nnx,nny,nnz
      INTEGER j, nn
      INTRINSIC MIN0
      EXTERNAL MYMALLOC
c***
c     initialize memory allocation scheme with MYMALLOC, compute ng
c***
      mem = 0
      nn = MIN0(nx,ny,nz)
      j =  nn - 2
c***
c     ng is the number of multigrid levels 
c***
      ng = 1
      do while (j .gt. 2)
         j = j/2
         ng = ng + 1
      enddo
      if ((2**ng + 2) .ne. nn) then 
         write (*,*) "STOP: nn .ne. 2**ng + 2", nn,2**ng + 2
         stop
      endif
c***
c     dynamically set indices, only up to ng levels
c***
      nnx = nx
      nny = ny
      nnz = nz
      do j=ng,1,-1
         irhs(j) = MYMALLOC(nnx*nny*nnz)
         ires(j) = MYMALLOC(nnx*nny*nnz)
         iu(j)   = MYMALLOC(nnx*nny*nnz)
         icc(j)  = MYMALLOC(nnx*nny*nnz)
         ia(j)   = MYMALLOC(nnx*nny*nnz)
         ib(j)   = MYMALLOC(nnx*nny*nnz)
         ic(j)   = MYMALLOC(nnx*nny*nnz)
         iei(j)  = MYMALLOC(nnx*nny*nnz)
         nn  = nn/2 + 1
         nnx = nnx/2 + 1
         nny = nny/2 + 1
         nnz = nnz/2 + 1
      enddo
      print*,'MEMLEN,mem:',MEMLEN,mem
c***
      return
      end
c***********************************************************************
      SUBROUTINE writeminmax(string,a,nx,ny,nz)
c***
      IMPLICIT none
      INTEGER nx,ny,nz,i,j,k,i1,j1,k1
      DOUBLE PRECISION a(nx,ny,nz), FINDMAX,FINDMIN,xmax,xmin
      CHARACTER *(*) string
      EXTERNAL FINDMAX,FINDMIN
c***
      xmax = FINDMAX(a,nx,ny,nz,i,j,k)
      xmin = FINDMIN(a,nx,ny,nz,i1,j1,k1)
      write(*,10) string,xmax,i,j,k,xmin,i1,j1,k1
   10 format(1x,a10,1p,1e12.5,'(',3i3,')',1x,'***',1x,
     %     1p,1e12.5,'(',3i3,')')
c***
      return
      end
c***********************************************************************
      SUBROUTINE nondim(vscale)
c***
      INCLUDE 'input.h'
	INTEGER n0,i,j,k
	DOUBLE PRECISION denom,vscale,x0,C_drag,gtotal_si,vel2_si,Re_2
c***
      do i=1,nx
           do j=1,ny
              do k=1,nz
              cmu(i,j,k)=0.09d0    
              enddo
           enddo
      enddo


	lengthscale=length_si
      vscale=dsqrt((xvel1_si-xvel2_si)*(xvel1_si-xvel2_si)+
     &       (yvel1_si-yvel2_si)*(yvel1_si-yvel2_si)+
     &       (zvel1_si-zvel2_si)*(zvel1_si-zvel2_si))
      C_drag = 1.
	vel2_si=dsqrt(xvel2_si*xvel2_si+yvel2_si*yvel2_si
     &	+zvel2_si*zvel2_si)
	gtotal_si=dsqrt(grav_si*grav_si+grav1_si*grav1_si)
      if((vscale.eq.0.).and.(gtotal_SI.ne.0.)) then
         vscale=sqrt((abs(gtotal_SI)*lengthscale)*ABS(rho1_SI-rho2_SI)/
     &    rho2_SI*(1./C_drag))
         Re_2=vel2_si*lengthscale*rho2_SI/visc2_SI
	   if(Re_2.le.10.and.vscale.eq.0.) C_drag = 10./Re_2
	   vscale=sqrt((abs(gtotal_SI)*lengthscale)*dabs(rho1_SI-rho2_SI)/
     &       rho2_SI*(1./C_drag))
      write (23,*) 'C_drag=',C_drag
       endif
c	vscale=1.
      if(vscale.eq.0.) vscale=1.
	timescale=lengthscale/vscale
      numer = length_SI
      denom = lengthscale
      lengthsize = numer/denom
      h = lengthsize/(nx-2)
      time=0.
c***
c     now derive adimensional physical variables to be stored in
c     the common physics, (be sure time in sec and length in meters)
c***
	tau_si=MAXCFL*0.8*(length_si/(nx-2))/vscale
      numer = tau_SI
      denom = timescale
      tau = numer/denom
	
c***
	tauh=tau/h
      tauh2 = tau/(h*h)
      tswap = 0
c***
      radius=radius_si/lengthscale

c    make velocities non-dimensional and rescale them according
c     to : v --> v1 = v*timescale/lengthsize --> v2 = v1*tau/h;

      xvel1=xvel1_si/vscale
      yvel1=yvel1_si/vscale
      zvel1=zvel1_si/vscale
      xvel2=xvel2_si/vscale
      yvel2=yvel2_si/vscale
      zvel2=zvel2_si/vscale

      n0 = INT(radius/h) + 1
      x0 = (2+n0)*h  
      y0 = (ny-2)*h/2
      z0 = (nz-2)*h/2
      write (*,*) 'n0,h,x0,y0,z0=',n0,h,x0,y0,z0
c***
      numer = rho2_SI
      denom = rho1_SI
      rho1  = 1.d0
      rho2  = numer/denom
c***
      denom = rho1_SI * lengthscale**2 / timescale 
      numer = visc1_SI
      mu1 = numer/denom
      numer = visc2_SI
      mu2   =  numer/denom
c***
      numer = sigma_SI
      denom = rho1_SI * lengthscale**3 / timescale**2
      sigma = numer/denom
c***
      numer = grav_SI
      denom = lengthscale / timescale**2
      g = numer/denom

      numer = grav1_SI
      denom = lengthscale / timescale**2
      g1 = numer/denom

	call fillzero(u2p,nx,ny,nz)   
	call fillzero(v2p,nx,ny,nz)
	call fillzero(w2p,nx,ny,nz)
	call fillzero(uvp,nx,ny,nz)
	call fillzero(uwp,nx,ny,nz)
	call fillzero(vwp,nx,ny,nz)   
      call fillzero(mut,nx,ny,nz)
	call fillzero(turbst,nx,ny,nz)
c      call fillzero(sigmturb,nx,ny,nz)
      return
      end
c     *********** NON-DIMENSIONAL numbers *************
      SUBROUTINE nondimnumbers(vscale)
      include 'input.h'
      REAL length_covered,La_1,oh_1,oh_2
      REAL tau_sigm1,tau_sigm2,tau_tvisc1,tau_tvisc2,CFL     
      DOUBLE PRECISION Re_1,Re_2, We_1,We_2,We_c,vel1_si,vel2_si,
     & vscale,h_si

	h_si=lengthscale/(nx-2)
      tau_sigm1=tau_SI*sqrt(sigma_SI/(rho1_SI*h_si*h_si*h_si))
      tau_sigm2=tau_SI*sqrt(sigma_SI/(rho2_SI*h_si*h_si*h_si))
      tau_tvisc1=tau_SI*visc1_SI/(rho1_SI*h_si*h_si)
      tau_tvisc2=tau_SI*visc2_SI/(rho2_SI*h_si*h_si)

	vel2_si=sqrt(xvel2_SI*xvel2_SI+yvel2_SI*yvel2_SI
     &	+zvel2_SI*zvel2_SI)
      vel1_si=sqrt(xvel1_SI*xvel1_SI+yvel1_SI*yvel1_SI
     &	+zvel1_SI*zvel1_SI)

      write (23,*)  '-------------------------'
      write (23,*) 'Non-Dim No. controlling the stability proprerties'
      write (23,*)  '-------------------------'
      write (23,*) 'tau/sim 1=',tau_sigm1
      write (23,*) 'tau/sim 2=',tau_sigm2
      write (23,*) 'tauSI/t_visc1=',tau_tvisc1
      write (23,*) 'tauSI/t_visc2=',tau_tvisc2
      if(vscale.gt.0.) then
      CFL=vscale*tau_SI/(h*lengthscale)
      write (23,*)   'CFL =',CFL
      length_covered = tmax*tau_SI*vscale/lengthscale
      write (23,*) 'distance covered',length_covered
      endif

      write (23,*)  '-------------------------'
      write (23,*)  'Other Non-dimensional Numbers'
      write (23,*)  '-------------------------'

      Re_1=vel1_si*lengthscale*rho1_SI/visc1_SI
      Re_2=vel2_si*lengthscale*rho2_SI/visc2_SI
      We_2=rho2_SI*vel2_si*vel2_si*(radius_si-y0_si)/sigma_SI
      We_1=rho1_SI*vel1_si*vel1_si*lengthscale/sigma_SI
      We_c=We_2*rho1_SI/((rho2_SI+rho1_SI)*2*3.14159265)
      La_1=lengthscale*sigma_SI*rho1_SI/(visc1_SI*visc1_SI)
      Oh_1=visc1_SI/sqrt(lengthscale*sigma_SI*rho1_SI)
      Oh_2=visc2_SI/sqrt(lengthscale*sigma_SI*rho2_SI)

      write (23,*) 'vscale, lengthscale', vscale , lengthscale 
	write (23,*) 'vel_2SI, vel_1SI',vel2_Si,vel1_SI  
      write (23,*)  'Re_1=',Re_1
      write (23,*)  'Re_2=',Re_2
      write (23,*)  'We_1=',We_1
      write (23,*)  'We_2=',We_2
      write (23,*)  'La_1=1/Oh_1**2',La_1
      write (23,*)  'Oh_1=',Oh_1
      write (23,*)  'Oh_2=',Oh_2

      write (6,*)  'vscale, lengthscale', vscale , lengthscale  
      write (6,*)  'Re_1=',Re_1
      write (6,*)  'Re_2=',Re_2
      write (6,*)  'We_1=',We_1
      write (6,*)  'We_2=',We_2
      write (6,*)  'La_1=',La_1
      write (6,*)  'Oh_1=',Oh_1
      write (6,*)  'Oh_2=',Oh_2

      return
      end

c    ********  creat initial conditions for jet *********

      SUBROUTINE inijet2(cc,dist)
      include 'input.h'
      INTEGER i,j,k
      DOUBLE PRECISION cc(nx,ny,nz),dist(nx,ny,nz)
      DOUBLE PRECISION yy,zz,rr

      do k=1,nz
         zz = k*h
         do j=1,ny
            yy = j*h
            rr= DSQRT((yy-y0)*(yy-y0) + (zz-z0)*(zz-z0))
              do i=1,nx
c     put pressure field in dist
                dist(i,j,k)=0.d0
               if(rr.le.radius) then
               cc(i,j,k)=1.d0
               u(i,j,k)=xvel1*(1.0-(dist(1,j,k)/radius)**4)**0.25
               v(i,j,k)=yvel1
		   	 w(i,j,k)=zvel1
			 if(iturb.eq.3.or.iturb.eq.4) then
			 ktur(i,j,k)=0.00015*xvel1*xvel1
			 eps(i,j,k)=0.09*rho1*ktur(i,j,k)*ktur(i,j,k)/mu1
	         mut(i,j,k)=0.09*rho1*ktur(i,j,k)*ktur(i,j,k)/eps(i,j,k)	
			 endif
           else
               u(i,j,k)=xvel2
               v(i,j,k)=yvel2
               w(i,j,k)=zvel2
               cc(i,j,k)=0.d0
			 if(iturb.eq.3.or.iturb.eq.4) then
			 ktur(i,j,k)=0.00015*(xvel2+dabs(xvel1-xvel2)*0.01)**2
			 eps(i,j,k)=0.09*rho2*ktur(i,j,k)*ktur(i,j,k)/mu2
	         mut(i,j,k)=0.09*rho2*ktur(i,j,k)*ktur(i,j,k)/eps(i,j,k)	
			 endif
            endif
         enddo
	  enddo
      enddo
      return
      end
c***********************************************************************
      SUBROUTINE inijet(cc,dist)
c***
      INCLUDE 'input.h'
      INTEGER i,j,k,n0
      DOUBLE PRECISION cc(nx,ny,nz),dist(nx,ny,nz)
      DOUBLE PRECISION zz,yy,outlet
      INTRINSIC DSQRT,DFLOAT,INT,DABS,DEXP
c***
      outlet=0.
      do k=1,nz
         zz = k*h
         do j=1,ny
            yy = j*h
            do i=1,nx
               dist(i,j,k) = DSQRT((yy-y0)*(yy-y0) + (zz-z0)*(zz-z0))
               if (dist(1,j,k).le.radius) then
                 u(1,j,k) = xvel1*(1.0-(dist(1,j,k)/radius)**4)**0.25
                 v(1,j,k) = yvel1
                 cc(1,j,k) = 1.d0
c              outlet=outlet+xvel1*(1.0-(dist(1,j,k)/radius)**4)**0.25
			 if(iturb.eq.3.or.iturb.eq.4) then
			 ktur(i,j,k)=0.00015*xvel1*xvel1
			 eps(i,j,k)=0.09*rho1*ktur(i,j,k)*ktur(i,j,k)/mu1
	         mut(i,j,k)=0.09*rho1*ktur(i,j,k)*ktur(i,j,k)/eps(i,j,k)	
			 endif
               else
                 u(1,j,k) = xvel2
                 v(1,j,k) = yvel2
                 cc(1,j,k) = 0.d0
c                 outlet = outlet + xvel2
			 if(iturb.eq.3.or.iturb.eq.4) then
			 ktur(i,j,k)=0.00015*(xvel2+dabs(xvel1-xvel2)*0.01)**2
			 eps(i,j,k)=0.09*rho2*ktur(i,j,k)*ktur(i,j,k)/mu2
	         mut(i,j,k)=0.09*rho2*ktur(i,j,k)*ktur(i,j,k)/eps(i,j,k)	
			 endif
			 endif
               if (dist(2,j,k).le.radius) then
                 u(2,j,k) = xvel1*(1.0-(dist(2,j,k)/radius)**4)**0.25
                 v(2,j,k) = yvel1
                 cc(2,j,k) = 1.d0
               else
                 u(2,j,k) = xvel2
                 v(2,j,k) = yvel2
                 cc(2,j,k) = 0.d0
               endif
               w(i,j,k) = 0.d0
            enddo 
         enddo 
      enddo 
      do k=1,nz
        do j=1,ny
c            u(nx,j,k) = outlet / (nx*ny*nz)
            u(nx,j,k) = u(nx-1,j,k)
        enddo
      enddo      
c     put pressure field in dist
c***
      do k=1,nz
         do j=1,ny
            do i=1,nx
               dist(i,j,k) = 0.d0
            enddo 
         enddo 
      enddo 
c***
      return
      end
c    ********  creat initial conditions for drop *********
      SUBROUTINE inidrop(cc,dist)
      include 'input.h'
      INTEGER i,j,k
      DOUBLE PRECISION cc(nx,ny,nz),dist(nx,ny,nz),ccf(nx,ny,nz)
      DOUBLE PRECISION x0,xx,yy,zz,rr
      x0=(nx-2)*h/2.0

c***	if z0_SI equal 1000 z0 will be in center of channel
	if (z0_si.ne.1000) z0= z0_SI/length_SI

	do k=1,nz
         zz = k*h
         do j=1,ny
            yy = j*h
            do i=1,nx
            xx = i*h
            rr= DSQRT((xx-x0)*(xx-x0)+(yy-y0)*(yy-y0)+(zz-z0)*(zz-z0))
c     put pressure field in dist
                dist(i,j,k) = 0.d0
               if(rr.le.radius) then
               cc(i,j,k)=1.d0
               u(i,j,k)=xvel1
               v(i,j,k)=yvel1
		   	 w(i,j,k)=zvel1
           else
c              u(i,j,k)=xvel2
c	*************new**********
			u(i,j,k)=xvel2*(-((k*1.d0-1.d0)/(nz*1.d0)-1.d0)**2+1.d0)

               v(i,j,k)=yvel2
               w(i,j,k)=zvel2
               cc(i,j,k)=0.0
            endif
         enddo
	  enddo
      enddo
      return
      end

c    ********  creat initial conditions for sheet *********

      SUBROUTINE inisheet(cc,dist)
      include 'input.h'
      INTEGER i,j,k
      DOUBLE PRECISION cc(nx,ny,nz),dist(nx,ny,nz)
      DOUBLE PRECISION yy,zz,rr

      do k=1,nz
         zz = k*h
         do j=1,ny
            yy = j*h
            
              do i=1,nx
c     put pressure field in dist
                dist(i,j,k) = 0.d0
               if(zz.le.radius) then
               cc(i,j,k)=1.d0
               u(i,j,k)=xvel1
               v(i,j,k)=yvel1
		   	 w(i,j,k)=zvel1
           else
               u(i,j,k)=xvel2
               v(i,j,k)=yvel2
               w(i,j,k)=zvel2
               cc(i,j,k)=0.d0
            endif
         enddo
	  enddo
      enddo
      return
      end
c    ********  READ initial data from file *********
      SUBROUTINE iniget(cc)
      include 'input.h'
      INTEGER i,j,k
      DOUBLE PRECISION cc(nx,ny,nz)
	open(2,form='unformatted',FILE='save.dat',STATUS='OLD') 
      read(2) t0
      read (2) nx        ! max points in x-direction
      read (2) ny        ! max points in y-direction
      read (2) nz        ! max points in z-direction
      read (2) tmax      ! maximum number of time steps
      read (2) nfilter   ! filter volume fraction (def.is zero) 
      read (2) cycles    ! max no. of cycles in multigrid method(def.is 60 )
      read (2) OMEGA     ! relaxation coefficient(def.is 1.0)
      read (2) NPRE      ! pre-relaxation coefficient(def.is 1.0)
      read (2) NTOP      ! post-relaxation coefficient(def.is 1.0)
      read (2) NPOST     ! top-relaxation coefficient(def.is 1.0)
      read (2) restart   ! =0 for new data, =1 for continueing the caculation(def.is 0)
      read (2) maxcfl    ! maximum value for CFL(def.is 0.8)
	read (2) iturb     ! =0 laminar,=1 mixing length with sigmtau model, =2
      				   !  mixing length with H'K' model, =3 k-epsilon model
      read (2) convectype !convection term, 1=upwind, 2=central, otherwise=no conv.(def.is 1)
      read (2) hswitch    !=0 off, =1 on

c      ***** physical and flow properties **********
      read (2) xvel1_SI   !initial x-component of velocity of fluid 1
      read (2) yvel1_SI   !initial y-component of velocity of fluid 1
      read (2) zvel1_SI   !initial z-component of velocity of fluid 1
      read (2) xvel2_SI   !initial x-component of velocity of fluid 2
      read (2) yvel2_SI   !initial y-component of velocity of fluid 2
      read (2) zvel2_SI   !initial z-component of velocity of fluid 2
      read (2) rho2_SI    !initial density of fluid 2
      read (2) rho1_SI    !initial density of fluid 1
      read (2) sigma_SI   !surface tension coefficient
      read (2) visc2_SI   !initial viscosity of fluid 2
      read (2) visc1_SI   !initial viscosity of fluid 1
      read (2) grav_SI    !garvity in z-direction
      read (2) grav1_SI   !garvity in x-direction

c      ****** geometry and flow type *******
      read (2) icase      ! module (1-jet, 2=drop, 3=bubble)
      read (2) ibc        ! bourdary condition type (=1, original, =2, standard) 
	read (2) length_SI  ! length of the domain in x-direction
      read (2) radius_si  ! module radius (NOTE: MUST BE CORRECTED!!) 
c      read (2) x0_SI     ! x location of center of jet, drop of bubble 
c	read (2) y0_SI      ! y location of center of jet, drop of bubble
c	read (2) z0_SI      ! z location of center of jet, drop of bubble

c     ****** print data ********
      read (2) tprint    ! print data on screen
      read (2) ntecplot  ! print data for ntecplot
      read (2) tsave     ! print unformatted data
c***
c     check nx, ny, nz against the maximum values in input.h
c***
      write (*,*)
      if (nx .ne. NMX) write (*,*) 'STOP: nx .ne. NMX:', nx, NMX
      if (ny .ne. NMY) write (*,*) 'STOP: ny .ne. NMY:', ny, NMY
      if (nz .ne. NMZ) write (*,*) 'STOP: nz .ne. NMZ:', nz, NMZ
c***
      if (ntecplot.eq.0) then
         write (*,*) 'STOP: missing all storing rates'
         stop
      endif
     	if (icase.eq. 1) then
	ijet=1
      idrop=0
	ibubble=0
	isheet=0
	elseif(icase.eq.2) then
	ijet=0
	idrop=1
	ibubble=0
	isheet=0
      elseif(icase.eq.3) then	
	ijet=0
	idrop=0
	ibubble=1
	isheet=0
      elseif(icase.eq.4) then	
	ijet=0
	idrop=0
	ibubble=0
      isheet=1
	else
	stop "incorrct module number"
	endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
	      read(2) u(i,j,k),v(i,j,k),w(i,j,k),cc(i,j,k)
            enddo
         enddo
      enddo
	tmax=t0+tmax
	return
	end
c***********************************************************************
      SUBROUTINE VOFmethod(t)
c***  
      INCLUDE 'input.h'
      INTEGER t, d
	INTRINSIC MOD
c***  
c     alternate first direction of interface motion
c***  
      tswap = tswap + 1
      if (tswap .gt. 3) tswap = 1
      if (MOD(tswap,3) .eq. 0) then
         call swp(w,z(icc(ng)),z(iei(ng)),z(ires(ng)),z(irhs(ng)),3)
         call swp(u,z(icc(ng)),z(iei(ng)),z(ires(ng)),z(irhs(ng)),1)
         call swp(v,z(icc(ng)),z(iei(ng)),z(ires(ng)),z(irhs(ng)),2)
      elseif (MOD(tswap,2) .eq. 0) then
         call swp(v,z(icc(ng)),z(iei(ng)),z(ires(ng)),z(irhs(ng)),2)
         call swp(w,z(icc(ng)),z(iei(ng)),z(ires(ng)),z(irhs(ng)),3)
         call swp(u,z(icc(ng)),z(iei(ng)),z(ires(ng)),z(irhs(ng)),1)
      else 
         call swp(u,z(icc(ng)),z(iei(ng)),z(ires(ng)),z(irhs(ng)),1)
         call swp(v,z(icc(ng)),z(iei(ng)),z(ires(ng)),z(irhs(ng)),2)
         call swp(w,z(icc(ng)),z(iei(ng)),z(ires(ng)),z(irhs(ng)),3)
      endif 
      call bc_c(z(icc(ng)))
c  ********* adding turbulent diffusion in F equation
      if (iturb.ne.0) call fturbdif(z(icc(ng)))
c***
      call calcrho (z(icc(ng)))
      return
	end
c***********************************************************************
      SUBROUTINE minmaxvcp (p,cc,t)
c***
c     debugging routine to write min and max and their position
c     of main fields
c***
      INCLUDE 'input.h'
      INTEGER t
      DOUBLE PRECISION p(nx,ny,nz), cc(nx,ny,nz)
c***
      call writeminmax ('u', u, nx, ny, nz)
      call writeminmax ('v', v, nx, ny, nz)
      call writeminmax ('w', w, nx, ny, nz)
      call writeminmax ('p', p, nx, ny, nz)
      call writeminmax ('cc', cc, nx, ny, nz)
c***
      return 
      end

c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
c Split advection of the interface along the x (d=1), y (d=2) and z (d=3)
c direction
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
      SUBROUTINE swp(us,cc,vof1,vof2,vof3,d)
c***
      INCLUDE 'input.h'
      INTEGER i,j,k,invx,invy,invz,d
      DOUBLE PRECISION us(nx,ny,nz),cc(nx,ny,nz),mx,my,mz,mm1,mm2
      DOUBLE PRECISION a1,a2,alpha,AL3D,FL3D
      DOUBLE PRECISION vof1(nx,ny,nz),vof2(nx,ny,nz),vof3(nx,ny,nz)
      INTRINSIC DMAX1,DMIN1
      EXTERNAL AL3D,FL3D
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
              a1 = us(i,j,k)*tauh
              if (d.eq.1) then
                a2 = us(i+1,j,k)*tauh
              elseif (d.eq.2) then
                a2 = us(i,j+1,k)*tauh
              elseif (d.eq.3) then
                a2 = us(i,j,k+1)*tauh
              endif
               
c***
c     3 cases: 1: DEFAULT (c=0. and fluxes=0.); 2: c=1.; 3:c>0.
c***
               vof1(i,j,k) = 0.0d0
               vof2(i,j,k) = 0.0d0
               vof3(i,j,k) = 0.0d0

               if (cc(i,j,k) .EQ. 1.0d0) then
                  vof1(i,j,k) = DMAX1(-a1,0.d0)
                  vof2(i,j,k) = 1.d0 - DMAX1(a1,0.d0) + DMIN1(a2,0.d0)
                  vof3(i,j,k) = DMAX1(a2,0.d0)

               else if (cc(i,j,k) .GT. 0.d0) then
c***
c     (1) normal vector: mx,my,mz; (2) mx,my,mz>0. and mx+my+mz = 1.;
c     (3) get alpha;               (4) back to original plane;
c     (5) lagrangian advection;    (6) get fluxes
c*(1)*
                mm1 = cc(i-1,j-1,k-1)+cc(i-1,j-1,k+1)+cc(i-1,j+1,k-1)
     %               +cc(i-1,j+1,k+1)+2.0d0*(cc(i-1,j-1,k)+cc(i-1,j+1,k)
     %               +cc(i-1,j,k-1)+cc(i-1,j,k+1))+4.0d0*cc(i-1,j,k)
                mm2 = cc(i+1,j-1,k-1)+cc(i+1,j-1,k+1)+cc(i+1,j+1,k-1)
     %               +cc(i+1,j+1,k+1)+2.0d0*(cc(i+1,j-1,k)+cc(i+1,j+1,k)
     %               +cc(i+1,j,k-1)+cc(i+1,j,k+1))+4.0d0*cc(i+1,j,k)
                mx = mm1 - mm2
                  
                mm1 = cc(i-1,j-1,k-1)+cc(i-1,j-1,k+1)+cc(i+1,j-1,k-1)
     %               +cc(i+1,j-1,k+1)+2.0d0*(cc(i-1,j-1,k)+cc(i+1,j-1,k)
     %               +cc(i,j-1,k-1)+cc(i,j-1,k+1))+4.0d0*cc(i,j-1,k)
                mm2 = cc(i-1,j+1,k-1)+cc(i-1,j+1,k+1)+cc(i+1,j+1,k-1)
     %               +cc(i+1,j+1,k+1)+2.0d0*(cc(i-1,j+1,k)+cc(i+1,j+1,k)
     %               +cc(i,j+1,k-1)+cc(i,j+1,k+1))+4.0d0*cc(i,j+1,k)
                my = mm1 - mm2
                  
                mm1 = cc(i-1,j-1,k-1)+cc(i-1,j+1,k-1)+cc(i+1,j-1,k-1)
     %               +cc(i+1,j+1,k-1)+2.0d0*(cc(i-1,j,k-1)+cc(i+1,j,k-1)
     %               +cc(i,j-1,k-1)+cc(i,j+1,k-1))+4.0d0*cc(i,j,k-1)
                mm2 = cc(i-1,j-1,k+1)+cc(i-1,j+1,k+1)+cc(i+1,j-1,k+1)
     %               +cc(i+1,j+1,k+1)+2.0d0*(cc(i-1,j,k+1)+cc(i+1,j,k+1)
     %               +cc(i,j-1,k+1)+cc(i,j+1,k+1))+4.0d0*cc(i,j,k+1)
                mz = mm1 - mm2
c*(2)*  
                invx = 1
                invy = 1
                invz = 1
                if (mx .LT. 0.0d0) then
                   mx = -mx
                   invx = -1
                endif
                if (my .LT. 0.0d0) then
                   my = -my
                   invy = -1
                endif
                if (mz .LT. 0.0d0) then
                   mz = -mz
                   invz = -1
                endif
                mm2 = mx+my+mz+1.e-30
                mx = mx/mm2
                my = my/mm2
                mz = mz/mm2
c*(3)*  
                  alpha = AL3D(mx,my,mz,cc(i,j,k))
c*(4)*  
                  mx = invx*mx
                  my = invy*my
                  mz = invz*mz
                  alpha = alpha + DMIN1(0.d0,mx) + DMIN1(0.d0,my) +
     %                 DMIN1(0.d0,mz)
c*(5)*  

                  mm1 = DMAX1(a1,0.0d0)
                  mm2 = 1.d0 - mm1 + DMIN1(0.d0,a2)
                  if (d.eq.1) then
                    mx = mx/(1.0d0 - a1 + a2)
                    alpha = alpha + mx*a1
                    if (a1 .LT. 0.d0) 
     %                 vof1(i,j,k) = FL3D(mx,my,mz,alpha,a1  ,-a1)
                    if (a2 .GT. 0.d0) 
     %                 vof3(i,j,k) = FL3D(mx,my,mz,alpha,1.d0,a2)
                    vof2(i,j,k) = FL3D(mx,my,mz,alpha,mm1,mm2)
                  elseif (d.eq.2) then
                    my = my/(1.0d0 - a1 + a2)
                    alpha = alpha + my*a1
                    if (a1 .LT. 0.d0) 
     %                 vof1(i,j,k) = FL3D(my,mz,mx,alpha,a1  ,-a1)
                    if (a2 .GT. 0.d0) 
     %                 vof3(i,j,k) = FL3D(my,mz,mx,alpha,1.d0,a2)
                    vof2(i,j,k) = FL3D(my,mz,mx,alpha,mm1,mm2)
                  elseif (d.eq.3) then
                    mz = mz/(1.0d0 - a1 + a2)
                    alpha = alpha + mz*a1
                    if (a1 .LT. 0.d0) 
     %                 vof1(i,j,k) = FL3D(mz,mx,my,alpha,a1  ,-a1)
                    if (a2 .GT. 0.d0) 
     %                 vof3(i,j,k) = FL3D(mz,mx,my,alpha,1.d0,a2)
                    vof2(i,j,k) = FL3D(mz,mx,my,alpha,mm1,mm2)
                  endif
               endif
            enddo
         enddo
      enddo
c***  
c     (1) apply proper boundary conditions to fluxes
c     (2) new values of cc and  clip it: 0. <= cc <= 1.
c     (3) apply proper boundary conditions to cc
c*(1)*  
      call bc_flux(vof1,vof3,us,d)
c*(2)*  
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
              if (d.eq.1) then
                cc(i,j,k) = vof3(i-1,j,k) + vof2(i,j,k) + vof1(i+1,j,k)
              elseif (d.eq.2) then
                cc(i,j,k) = vof3(i,j-1,k) + vof2(i,j,k) + vof1(i,j+1,k)
              elseif (d.eq.3) then
                cc(i,j,k) = vof3(i,j,k-1) + vof2(i,j,k) + vof1(i,j,k+1)
              endif
              cc(i,j,k) = DMAX1(0.0d0,DMIN1(1.0d0,cc(i,j,k)))
           enddo
         enddo
      enddo
c*(3)*
      call bc_c(cc)
c***



      return
      end

*********************************
      SUBROUTINE printtheta(cc,t)
      INCLUDE 'input.h'
c	IMPLICIT none
      INTEGER i,j,k,nf,ii,iback,ifront,t
      DOUBLE PRECISION cc(nx,ny,nz),ccf(nx,ny,nz),theta1,theta2
	DOUBLE PRECISION theta,xx,zz,phi,xc,zc,mx,my,mz,mm1,mm2
	CHARACTER (128) FileName 
	nf=3
	call copy(ccf,cc,nx,ny,nz)
c	do i=1,nf
c	call filterc27(ccf,temp,nx,ny,nz)
c	call copy (ccf,temp,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	enddo
	do i=1,nf
	call filterc27(ccf,ccf,nx,ny,nz)
	call bc_out(ccf,nx,ny,nz)
	enddo
	write (FileName, '(A,I4.4,A)') 'theta', t/10, '.plt' 
      OPEN(5,FILE=FileName) 
c		open(5,file='theta.plt')
		write(5,*) 'VARIABLES="phi","theta"'

	theta1=datand(y0/radius_si)
	theta2=datand((y0+h*lengthscale)/radius_si)
		
		j=2

		k=(nz-2)/2
		ii=0
		do i=1,nx
			if (cc(i,j,k).gt.0.d0 .and. cc(i,j,k).lt.0.99) then
			ii=ii+1
			if (ii.eq.1) iback=i
			if (ii.gt.1) ifront=i
			endif
		enddo
	xc=(ifront+iback)/2

c		xc=(nx-2)/2.0
		zc=(nz-2)/2.0

		do k=2,nz-1
		do i=2,nx-1
		
		if (cc(i,j,k).gt.0.d0 .and. cc(i,j,k).lt.1.d0) then

                mm1 = ccf(i-1,j-1,k-1)+ccf(i-1,j-1,k+1)+ccf(i-1,j+1,k-1)
     %            +ccf(i-1,j+1,k+1)+2.0d0*(ccf(i-1,j-1,k)+ccf(i-1,j+1,k)
     %               +ccf(i-1,j,k-1)+ccf(i-1,j,k+1))+4.0d0*ccf(i-1,j,k)
                mm2 = ccf(i+1,j-1,k-1)+ccf(i+1,j-1,k+1)+ccf(i+1,j+1,k-1)
     %            +ccf(i+1,j+1,k+1)+2.0d0*(ccf(i+1,j-1,k)+ccf(i+1,j+1,k)
     %                +ccf(i+1,j,k-1)+ccf(i+1,j,k+1))+4.0d0*ccf(i+1,j,k)
                mx = mm1 - mm2
                  
                mm1 = ccf(i-1,j-1,k-1)+ccf(i-1,j-1,k+1)+ccf(i+1,j-1,k-1)
     %            +ccf(i+1,j-1,k+1)+2.0d0*(ccf(i-1,j-1,k)+ccf(i+1,j-1,k)
     %               +ccf(i,j-1,k-1)+ccf(i,j-1,k+1))+4.0d0*ccf(i,j-1,k)
                mm2 = ccf(i-1,j+1,k-1)+ccf(i-1,j+1,k+1)+ccf(i+1,j+1,k-1)
     %            +ccf(i+1,j+1,k+1)+2.0d0*(ccf(i-1,j+1,k)+ccf(i+1,j+1,k)
     %               +ccf(i,j+1,k-1)+ccf(i,j+1,k+1))+4.0d0*ccf(i,j+1,k)
                my = mm1 - mm2
                  
                mm1 = ccf(i-1,j-1,k-1)+ccf(i-1,j+1,k-1)+ccf(i+1,j-1,k-1)
     %            +ccf(i+1,j+1,k-1)+2.0d0*(ccf(i-1,j,k-1)+ccf(i+1,j,k-1)
     %               +ccf(i,j-1,k-1)+ccf(i,j+1,k-1))+4.0d0*ccf(i,j,k-1)
                mm2 = ccf(i-1,j-1,k+1)+ccf(i-1,j+1,k+1)+ccf(i+1,j-1,k+1)
     %            +ccf(i+1,j+1,k+1)+2.0d0*(ccf(i-1,j,k+1)+ccf(i+1,j,k+1)
     %               +ccf(i,j-1,k+1)+ccf(i,j+1,k+1))+4.0d0*ccf(i,j,k+1)
                mz = mm1 - mm2




			theta=datand(dsqrt(mx*mx+mz*mz)/my)
			if (theta .lt.0) theta=180+theta

			xx=i-xc
			zz=k-zc
			phi=datand(zz/(xx+0.000000001))
			if (xx.lt.0) phi=phi+180
			if (xx.ge.0 .and. zz.lt.0) phi=phi+360
			theta=theta+theta2-theta1
			write(5,122) phi,theta
			

		endif






	enddo
	enddo
	close(5)
 122	format (2(1x,e15.7))
	
	return 
	end

*********************************
      SUBROUTINE thetahysteresis(cc,nx,ny,nz,tsi)
c      INCLUDE 'input.h'
	IMPLICIT none
      INTEGER i,j,k,nf,nx,ny,nz,ii
      DOUBLE PRECISION cc(nx,ny,nz),ccf(nx,ny,nz) !,temp(nx,ny,nz)
	DOUBLE PRECISION theta,xx,zz,phi,xc,yc,mx,my,mz,mm1,mm2
	DOUBLE PRECISION thetah,thetaa,thetar,tsi
	nf=3
c	call copy(ccf,cc,nx,ny,nz)
c	do i=1,nf
c	call filterc27(ccf,temp,nx,ny,nz)
c	call copy (ccf,temp,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	enddo
	call copy(ccf,cc,nx,ny,nz)
	do i=1,nf
	call filterc27(ccf,ccf,nx,ny,nz)
	call bc_out(ccf,nx,ny,nz)
	enddo		
		
		j=2		
		k=(nz-2)/2
		ii=0
		do i=2,nx-1
		
		if (cc(i,j,k).gt.0.d0 .and. cc(i,j,k).lt.1.d0) then
	ii=ii+1

                mm1 = ccf(i-1,j-1,k-1)+ccf(i-1,j-1,k+1)+ccf(i-1,j+1,k-1)
     %            +ccf(i-1,j+1,k+1)+2.0d0*(ccf(i-1,j-1,k)+ccf(i-1,j+1,k)
     %               +ccf(i-1,j,k-1)+ccf(i-1,j,k+1))+4.0d0*ccf(i-1,j,k)
                mm2 = ccf(i+1,j-1,k-1)+ccf(i+1,j-1,k+1)+ccf(i+1,j+1,k-1)
     %            +ccf(i+1,j+1,k+1)+2.0d0*(ccf(i+1,j-1,k)+ccf(i+1,j+1,k)
     %                +ccf(i+1,j,k-1)+ccf(i+1,j,k+1))+4.0d0*ccf(i+1,j,k)
                mx = mm1 - mm2
                  
                mm1 = ccf(i-1,j-1,k-1)+ccf(i-1,j-1,k+1)+ccf(i+1,j-1,k-1)
     %            +ccf(i+1,j-1,k+1)+2.0d0*(ccf(i-1,j-1,k)+ccf(i+1,j-1,k)
     %               +ccf(i,j-1,k-1)+ccf(i,j-1,k+1))+4.0d0*ccf(i,j-1,k)
                mm2 = ccf(i-1,j+1,k-1)+ccf(i-1,j+1,k+1)+ccf(i+1,j+1,k-1)
     %            +ccf(i+1,j+1,k+1)+2.0d0*(ccf(i-1,j+1,k)+ccf(i+1,j+1,k)
     %               +ccf(i,j+1,k-1)+ccf(i,j+1,k+1))+4.0d0*ccf(i,j+1,k)
                my = mm1 - mm2
                  
                mm1 = ccf(i-1,j-1,k-1)+ccf(i-1,j+1,k-1)+ccf(i+1,j-1,k-1)
     %            +ccf(i+1,j+1,k-1)+2.0d0*(ccf(i-1,j,k-1)+ccf(i+1,j,k-1)
     %               +ccf(i,j-1,k-1)+ccf(i,j+1,k-1))+4.0d0*ccf(i,j,k-1)
                mm2 = ccf(i-1,j-1,k+1)+ccf(i-1,j+1,k+1)+ccf(i+1,j-1,k+1)
     %            +ccf(i+1,j+1,k+1)+2.0d0*(ccf(i-1,j,k+1)+ccf(i+1,j,k+1)
     %               +ccf(i,j-1,k+1)+ccf(i,j+1,k+1))+4.0d0*ccf(i,j,k+1)
                mz = mm1 - mm2




			theta=datand(dsqrt(mx*mx+mz*mz)/my)
			if (theta .lt.0) theta=180+theta
			if (ii.eq.1) thetaR=theta
			if (ii.gt.1) thetaA=theta
			thetaH=thetaA-thetaR
			
			
		
		endif
	enddo
	write(27,122) tsi,thetaH,thetaA,thetaR
 122	format (4(1x,e15.7))
	return 
	end
************************ print theta 2 dimentional*******************
      SUBROUTINE printtheta2d(cc)
      INCLUDE 'input.h'
      INTEGER i,j,k
      DOUBLE PRECISION cc(nx,ny,nz),mx,my,mz,mm1,mm2,theta
	  		
		open(5,file='theta.plt')
		write(5,*) 'VARIABLES="x","y","theta"'
		write(5,*) 'zone i=',nx-2,'j=',ny-2
		
		k=3
		do i=2,nx-1
		do j=2,ny-1
	
		if (cc(i,j,k).gt.0.d0 .and. cc(i,j,k).lt.1.d0) then

                mm1 = cc(i-1,j-1,k-1)+cc(i-1,j-1,k+1)+cc(i-1,j+1,k-1)
     %               +cc(i-1,j+1,k+1)+2.0d0*(cc(i-1,j-1,k)+cc(i-1,j+1,k)
     %               +cc(i-1,j,k-1)+cc(i-1,j,k+1))+4.0d0*cc(i-1,j,k)
                mm2 = cc(i+1,j-1,k-1)+cc(i+1,j-1,k+1)+cc(i+1,j+1,k-1)
     %               +cc(i+1,j+1,k+1)+2.0d0*(cc(i+1,j-1,k)+cc(i+1,j+1,k)
     %               +cc(i+1,j,k-1)+cc(i+1,j,k+1))+4.0d0*cc(i+1,j,k)
                mx = mm1 - mm2
                  
                mm1 = cc(i-1,j-1,k-1)+cc(i-1,j-1,k+1)+cc(i+1,j-1,k-1)
     %               +cc(i+1,j-1,k+1)+2.0d0*(cc(i-1,j-1,k)+cc(i+1,j-1,k)
     %               +cc(i,j-1,k-1)+cc(i,j-1,k+1))+4.0d0*cc(i,j-1,k)
                mm2 = cc(i-1,j+1,k-1)+cc(i-1,j+1,k+1)+cc(i+1,j+1,k-1)
     %               +cc(i+1,j+1,k+1)+2.0d0*(cc(i-1,j+1,k)+cc(i+1,j+1,k)
     %               +cc(i,j+1,k-1)+cc(i,j+1,k+1))+4.0d0*cc(i,j+1,k)
                my = mm1 - mm2
                  
                mm1 = cc(i-1,j-1,k-1)+cc(i-1,j+1,k-1)+cc(i+1,j-1,k-1)
     %               +cc(i+1,j+1,k-1)+2.0d0*(cc(i-1,j,k-1)+cc(i+1,j,k-1)
     %               +cc(i,j-1,k-1)+cc(i,j+1,k-1))+4.0d0*cc(i,j,k-1)
                mm2 = cc(i-1,j-1,k+1)+cc(i-1,j+1,k+1)+cc(i+1,j-1,k+1)
     %               +cc(i+1,j+1,k+1)+2.0d0*(cc(i-1,j,k+1)+cc(i+1,j,k+1)
     %               +cc(i,j-1,k+1)+cc(i,j+1,k+1))+4.0d0*cc(i,j,k+1)
                mz = mm1 - mm2




			theta=datand(sqrt(mx*mx+my*my)/mz)
		endif
		if(cc(i,j,k).eq.0.d0 .or. cc(i,j,k).eq.1.d0) then
			theta=0.d0

		endif
	if (theta .lt.0) then
	theta=180+theta
	endif
	write(5,122) i*h,j*h,theta


 122	format (3(1x,e15.7))

	enddo
	enddo
	close(5)
	return
	end


c  ********* adding turbulent diffusion in F equation
      subroutine fturbdif(cc)
      include 'input.h'
      integer i,j,k
      DOUBLE PRECISION ftr,const,cmin,cmax,FMODMAX,FINDMIN
	DOUBLE PRECISION cc(nx,ny,nz)
      
        const=1.0d0
        ftr=1. ! (default is 1)if it is zero, no turbulence in F equation 
	  ftr=ftr*tauh2*const/2.d0
           do i=2,nx-1
	     do j=2,ny-1
		 do k=2,nz-1
	     cc(i,j,k)=cc(i,j,k)+ftr*(
     &	 (mut(i,j,k)/rhon(i,j,k)+mut(i,j+1,k)/rhon(i,j+1,k))*
     &	 (cc(i,j+1,k)-cc(i,j,k))-(mut(i,j,k)/rhon(i,j,k)+
     &	 mut(i,j-1,k)/rhon(i,j-1,k))*(cc(i,j,k)-cc(i,j-1,k))+
     &	 (mut(i,j,k)/rhon(i,j,k)+mut(i,j,k+1)/rhon(i,j,k+1))*
     &     (cc(i,j,k+1)-cc(i,j,k))-(mut(i,j,k)/rhon(i,j,k)+mut(i,j,k-1)/
     &     rhon(i,j,k-1))*(cc(i,j,k)-cc(i,j,k-1)))
	     enddo
		 enddo
           enddo
	      call bc_c(cc)
	cmin=FINDMIN(cc,nx,ny,nz,i,j,k)
	cmax=FMODMAX(cc,nx,ny,nz,i,j,k)
	if(cmax.gt.1.0.or.cmin.lt.0.0)
     &		write(*,*)' WARNING, cmin, cmax= ',cmin,cmax
	return
      end
c  **********************************************************************
c  ***** filter velocity field *******      
      SUBROUTINE filter(ufilt)
      INCLUDE 'input.h'
	INTEGER i,j,k
      DOUBLE PRECISION ufilt(nx,ny,nz)
      do j=2,ny-1
         do i=2,nx-1
	      do k=2,nz-1
            ufilt(i,j,k) = 0.25*u(i,j,k)
     *  + 0.125*(u(i,j,k-1)+u(i,j,k+1)+u(i,j-1,k)+u(i,j+1,k))
     *  + 0.0625*(u(i,j-1,k-1)+u(i,j-1,k+1)+u(i,j+1,k-1)+u(i,j+1,k+1))
            enddo
	   enddo
      enddo
      return
      end
c***********************************************************************
      SUBROUTINE fillzero(u,nx,ny,nz)
c***
c     fill of zeroes a 3D array
c***
      IMPLICIT none
      INTEGER nx,ny,nz,i,j,k
      DOUBLE PRECISION u(nx,ny,nz)
c***
      do k=1,nz
         do j=1,ny
            do i=1,nx
               u(i,j,k)=0.d0
            enddo
         enddo
      enddo
c***
      return
      end







c----------	
	subroutine CalcH3d(cc)
c
c     ******************************************************************
c

	include 'input.h'
	  integer i,j,k
	  DOUBLE PRECISION nix,niy,niz,cc(nx,ny,nz)
c        DOUBLE PRECISION hxmin(nx,ny,nz),hxmax(nx,ny,nz),hymin(nx,ny,nz)
c        DOUBLE PRECISION hymax(nx,ny,nz),hzmin(nx,ny,nz),hzmax(nx,ny,nz)
        DOUBLE PRECISION absxn,absyn,abszn,s1,s2,s3
        DOUBLE PRECISION h1a,h1b,h2a,h2b,h3a,h3b,hmintot,hmaxtot  
	  DOUBLE PRECISION n1,n2,n3,clim1,diag,rsixn3,loop,delta
        DOUBLE PRECISION mm1,mm2,normx,normy,normz,ccnew
	  INTRINSIC DABS,DSQRT      
        
	do k=2,nz-1
	   do j=2,ny-1
		do i=2,nx-1
		
!" nix, niy, niz are unit normal vector components"


                mm1 = cc(i-1,j-1,k-1)+cc(i-1,j-1,k+1)+cc(i-1,j+1,k-1)
     %               +cc(i-1,j+1,k+1)+2.0d0*(cc(i-1,j-1,k)+cc(i-1,j+1,k)
     %               +cc(i-1,j,k-1)+cc(i-1,j,k+1))+4.0d0*cc(i-1,j,k)
                mm2 = cc(i+1,j-1,k-1)+cc(i+1,j-1,k+1)+cc(i+1,j+1,k-1)
     %               +cc(i+1,j+1,k+1)+2.0d0*(cc(i+1,j-1,k)+cc(i+1,j+1,k)
     %               +cc(i+1,j,k-1)+cc(i+1,j,k+1))+4.0d0*cc(i+1,j,k)
                normx = mm1 - mm2

                mm1 = cc(i-1,j-1,k-1)+cc(i-1,j-1,k+1)+cc(i+1,j-1,k-1)
     %               +cc(i+1,j-1,k+1)+2.0d0*(cc(i-1,j-1,k)+cc(i+1,j-1,k)
     %               +cc(i,j-1,k-1)+cc(i,j-1,k+1))+4.0d0*cc(i,j-1,k)
               mm2 = cc(i-1,j+1,k-1)+cc(i-1,j+1,k+1)+cc(i+1,j+1,k-1)
     %               +cc(i+1,j+1,k+1)+2.0d0*(cc(i-1,j+1,k)+cc(i+1,j+1,k)
     %               +cc(i,j+1,k-1)+cc(i,j+1,k+1))+4.0d0*cc(i,j+1,k)
                normy = mm1 - mm2

                mm1 = cc(i-1,j-1,k-1)+cc(i-1,j+1,k-1)+cc(i+1,j-1,k-1)
     %               +cc(i+1,j+1,k-1)+2.0d0*(cc(i-1,j,k-1)+cc(i+1,j,k-1)
     %               +cc(i,j-1,k-1)+cc(i,j+1,k-1))+4.0d0*cc(i,j,k-1)
                mm2 = cc(i-1,j-1,k+1)+cc(i-1,j+1,k+1)+cc(i+1,j-1,k+1)
     %               +cc(i+1,j+1,k+1)+2.0d0*(cc(i-1,j,k+1)+cc(i+1,j,k+1)
     %               +cc(i,j-1,k+1)+cc(i,j+1,k+1))+4.0d0*cc(i,j,k+1)
                normz = mm1 - mm2

        nix=normx/DSQRT(normx**2+normy**2+normz**2+1.d-50)
        niy=normy/DSQRT(normx**2+normy**2+normz**2+1.d-50)
        niz=normz/DSQRT(normx**2+normy**2+normz**2+1.d-50)


		   if (cc(i,j,k).gt.0.5d0) then
			ccnew=1.d0-cc(i,j,k)
		   else
              ccnew=cc(i,j,k)
		   endif

c
c................ order the normals
c
		   absxn=dabs(nix)
		   absyn=dabs(niy)
		   abszn=dabs(niz)
		   n1=dmin1(absxn,absyn,abszn)
		   n3=dmax1(absxn,absyn,abszn)
		   n2=absxn+absyn+abszn-n1-n3
		   

		   

      clim1=dsqrt(2.d0*ccnew*n3*n2)

      if ((n1.eq.0d0).and.(n2.gt.0d0).and.(ccnew.le.clim1))  then ! or if (n1.le.1d-50)

      s2=dsqrt(2*ccnew*n3/n2)
      s3=dsqrt(2*ccnew*n2/n3)

      h1a=s2*s3/2.d0
      h1b=s2*s3/2.d0
      h2a=s3
      h2b=0.d0
      h3a=s2
      h3b=0.d0

		go to 199

      endif



      if ((n1.eq.0d0).and.(n2.gt.0d0).and.(ccnew.gt.clim1)) then  !"????????"

      h1a=cc(i,j,k)   !"????????"
      h1b=cc(i,j,k)   !"????????"
      h2a=cc(i,j,k)+n2/(2*n3)  !"????????"
      h2b=cc(i,j,k)-n2/(2*n3)  !"????????"
      h3a=1.d0
      h3b=0.d0

		go to 199

      endif


      if ((n1.eq.0d0).and.(n2.eq.0d0))  then
      h1a=ccnew
      h1b=ccnew
      h2a=ccnew
      h2b=ccnew
      h3a=1.d0
      h3b=0.d0

		go to 199

      endif

c
c................ triangular section, figures 5,11
c


		   if (6.d0*n2*n3*ccnew.lt.n1**2) then
			diag=(6.d0*n1*n2*n3*ccnew)**(1.d0/3.d0)

         s1=diag/n1		
         s2=diag/n2		
         s3=diag/n3		

! Case (a)  Fig.(5)
      h1a=(s2*s3)/2.d0
      h1b=0.d0
      h2a=(s1*s3)/2.d0
      h2b=0.d0
      h3a=(s1*s2)/2.d0
      h3b=0.d0

c                 write (9,*) nsubcyc,ijk,"triangle"
			go to 199
		   endif
	
c
c................ quadrilateral section A, figures 6,12,13
c

		   if (6.d0*n2*n3*ccnew.lt.3.d0*n2**2-3.d0*n1*n2+n1**2) then
			diag=(3.d0*n1+dsqrt(72.d0*n2*n3*ccnew-3.d0*n1**2))
     &                 /6.d0

         s1=diag/n1		
         s2=diag/n2		
         s3=diag/n3		

! Case (b)  Fig.(6)
      h1a=(s2*s3)/2.d0
      h1b=(s2*s3)/2.d0*((s1-1.d0)/s1)**2
      h2a=s3*(2.d0*s1-1.d0)/(2.d0*s1)
      h2b=0.d0
      h3a=s2*(2.d0*s1-1.d0)/(2.d0*s1)
      h3b=0.d0

c                 write (9,*) nsubcyc,ijk,"quadrilateral A"
			go to 199
		   endif
		
c
c................ pentagonal section, figures 7,14
c
		   if ((n1+n2.ge.n3.and.6.d0*n1*n2*n3*ccnew.lt.n3**3
     1            -(n3-n1)**3-(n3-n2)**3).or.(n1+n2.lt.n3.and.
     2            2.d0*n3*ccnew.lt.n1+n2)) then
			rsixn3=1.d0/(6.d0*n1*n2*n3)
			diag=n2
  130             do 131 loop=1,1000
			   delta=-(rsixn3*(diag**3-(diag-n1)**3-(diag-n2)**3)       !"ALI: Newton-Raphson iteration"
     &                     -ccnew)/(3.d0*rsixn3*(diag**2-(diag-n1)**2
     &                     -(diag-n2)**2))
			   diag=diag+delta
			   if (dabs(delta).lt.1.d-6) goto 132
  131             continue
  132             continue

         s1=diag/n1		
         s2=diag/n2		
         s3=diag/n3		

! Case (c)  Fig.(7)
      h1a=s3*(2*s2-1.d0)/(2.d0*s2)
      h1b=(s2*s3)/2.d0*((s1-1.d0)/s1)**2
      h2a=s3*(2*s1-1.d0)/(2.d0*s1)
      h2b=(s1*s3)/2.d0*((s2-1.d0)/s2)**2
      h3a=2.d0*(s1+s2)-s1/s2-s2/s1-s1*s2-1.d0
      h3b=0.d0

			go to 199
		   endif
		
c
c................ hexagonal section, figures 8,15
c
		   if ((n1+n2.gt.n3).and.6.d0*n1*n2*n3*ccnew.ge.n3**3
     &             -(n3-n1)**3-(n3-n2)**3) then
			rsixn3=1.d0/(6.d0*n1*n2*n3)
			diag=n3
  140             do 141 loop=1,1000
			   delta=-(rsixn3*(diag**3-(diag-n1)**3-(diag-n2)**3        !"ALI: Newton-Raphson iteration"
     1                     -(diag-n3)**3)-ccnew)/(3.d0*rsixn3*(diag**2
     2                     -(diag-n1)**2-(diag-n2)**2-(diag-n3)**2))
			   diag=diag+delta
			   if (dabs(delta).lt.1.d-6) go to 142
  141             continue
  142             continue

         s1=diag/n1		
         s2=diag/n2		
         s3=diag/n3		

! Case (d)    Fig.(8)
      h1a=2.d0*(s2+s3)-s2/s3-s3/s2-s2*s3-1.d0
      h1b=(s2*s3)/2.d0*((s1-1.d0)/s1)**2
      h2a=2.d0*(s1+s3)-s1/s3-s3/s1-s1*s3-1.d0
      h2b=(s1*s3)/2.d0*((s2-1.d0)/s2)**2
      h3a=2.d0*(s1+s2)-s1/s2-s2/s1-s1*s2-1.d0
      h3b=(s1*s2)/2.d0*((s3-1.d0)/s3)**2

  			go to 199
		   endif
		
c
c................ quadrilateral section B, figure 9
c
		   if ((n1+n2.lt.n3).and.6.d0*n1*n2*n3*ccnew.ge.
     &             (n1+n2)**3-n1**3-n2**3) then

  150          diag=0.5d0*(n1+n2)+ccnew*n3

         s1=diag/n1		
         s2=diag/n2		
         s3=diag/n3		

! case (e)          Fig.(9)
      h1a=s3*(2.d0*s2-1.d0)/(2.d0*s2)
      h1b=s3*(2.d0*s2-1.d0)*(s1-1.d0)/(2.d0*s1*s2)
      h2a=s3*(2.d0*s1-1.d0)/(2.d0*s1)
      h2b=s3*(2.d0*s1-1.d0)*(s2-1.d0)/(2.d0*s1*s2)
      h3a=1.d0
      h3b=0.d0		

c              write (9,*) nsubcyc,ijk,"quadrilateral B"
		   go to 199
		   
		   endif
		
!      "HERE if none of the above cases are not found, the program STOPS"		

		   Write (*,*) "NO MATCH CASES FOUND"
		   stop
		

  199          continue

          if (cc(i,j,k).gt.0.5) then      

          h1a=1.d0-h1a
          h1b=1.d0-h1b
          h2a=1.d0-h2a
          h2b=1.d0-h2b
          h3a=1.d0-h3a
          h3b=1.d0-h3b
          
          endif

                      if(absxn.eq.n1) then

                         if (nix.lt.0) then
                         hxmin(i,j,k)=h1a
                         hxmax(i,j,k)=h1b
                         else
                         hxmin(i,j,k)=h1b
                         hxmax(i,j,k)=h1a
                         endif

                                if(dabs(absyn-n2).lt.1.d-6) then

                                        if (niy.lt.0) then
                                        hymin(i,j,k)=h2a
                                        hymax(i,j,k)=h2b
                                        else
                                        hymin(i,j,k)=h2b
                                        hymax(i,j,k)=h2a
                                        endif

                                        if (niz.lt.0) then
                                        hzmin(i,j,k)=h3a
                                        hzmax(i,j,k)=h3b
                                        else
                                        hzmin(i,j,k)=h3b
                                        hzmax(i,j,k)=h3a
                                        endif

                                elseif(dabs(abszn-n2).lt.1.d-6) then

                                        if (niy.lt.0) then
                                        hymin(i,j,k)=h3a
                                        hymax(i,j,k)=h3b
                                        else
                                        hymin(i,j,k)=h3b
                                        hymax(i,j,k)=h3a
                                        endif

                                        if (niz.lt.0) then
                                        hzmin(i,j,k)=h2a
                                        hzmax(i,j,k)=h2b
                                        else
                                        hzmin(i,j,k)=h2b
                                        hzmax(i,j,k)=h2a
                                        endif

                                else

                                write (*,*) "PROBLEM finding MIN and MAX
     &                           unit vector components"
                                stop

                                endif

                      elseif(absyn.eq.n1) then

                         if (niy.lt.0) then
                         hymin(i,j,k)=h1a
                         hymax(i,j,k)=h1b
                         else
                         hymin(i,j,k)=h1b
                         hymax(i,j,k)=h1a
                         endif

                                if(dabs(absxn-n2).lt.1.d-6) then

                                        if (nix.lt.0) then
                                        hxmin(i,j,k)=h2a
                                        hxmax(i,j,k)=h2b
                                        else
                                        hxmin(i,j,k)=h2b
                                        hxmax(i,j,k)=h2a
                                        endif

                                        if (niz.lt.0) then
                                        hzmin(i,j,k)=h3a
                                        hzmax(i,j,k)=h3b
                                        else
                                        hzmin(i,j,k)=h3b
                                        hzmax(i,j,k)=h3a
                                        endif

                                elseif(dabs(abszn-n2).lt.1.d-6) then

                                        if (nix.lt.0) then
                                        hxmin(i,j,k)=h3a
                                        hxmax(i,j,k)=h3b
                                        else
                                        hxmin(i,j,k)=h3b
                                        hxmax(i,j,k)=h3a
                                        endif

                                        if (niz.lt.0) then
                                        hzmin(i,j,k)=h2a
                                        hzmax(i,j,k)=h2b
                                        else
                                        hzmin(i,j,k)=h2b
                                        hzmax(i,j,k)=h2a
                                        endif

                                else

                                write (*,*) "PROBLEM finding MIN and MAX
     &                           unit vector components"
                                stop

                                endif

                      elseif(abszn.eq.n1) then

                         if (niz.lt.0) then
                         hzmin(i,j,k)=h1a
                         hzmax(i,j,k)=h1b
                         else
                         hzmin(i,j,k)=h1b
                         hzmax(i,j,k)=h1a
                         endif

                                if(dabs(absxn-n2).lt.1.d-6) then

                                        if (nix.lt.0) then
                                        hxmin(i,j,k)=h2a
                                        hxmax(i,j,k)=h2b
                                        else
                                        hxmin(i,j,k)=h2b
                                        hxmax(i,j,k)=h2a
                                        endif

                                        if (niz.lt.0) then
                                        hymin(i,j,k)=h3a
                                        hymax(i,j,k)=h3b
                                        else
                                        hymin(i,j,k)=h3b
                                        hymax(i,j,k)=h3a
                                        endif

                                elseif(dabs(absyn-n2).lt.1.d-6) then

                                        if (nix.lt.0) then
                                        hxmin(i,j,k)=h3a
                                        hxmax(i,j,k)=h3b
                                        else
                                        hxmin(i,j,k)=h3b
                                        hxmax(i,j,k)=h3a
                                        endif

                                        if (niz.lt.0) then
                                        hymin(i,j,k)=h2a
                                        hymax(i,j,k)=h2b
                                        else
                                        hymin(i,j,k)=h2b
                                        hymax(i,j,k)=h2a
                                        endif

                                else

                                write (*,*) "PROBLEM finding MIN and MAX
     &                           unit vector components"
                                stop

                                endif

                      endif 



      if (cc(i,j,k).le.0.5) then     ! "EXCEPTIONS are the CASES I added to YOUNGS method "

      call exchange(hxmin(i,j,k),hxmax(i,j,k))
      call exchange(hymin(i,j,k),hymax(i,j,k))
      call exchange(hzmin(i,j,k),hzmax(i,j,k))

      endif

	if (hxmin(i,j,k).lt.0.d0) hxmin(i,j,k)=0.d0
	if (hymin(i,j,k).lt.0.d0) hymin(i,j,k)=0.d0
	if (hzmin(i,j,k).lt.0.d0) hzmin(i,j,k)=0.d0
	if (hxmax(i,j,k).lt.0.d0) hxmax(i,j,k)=0.d0
	if (hymax(i,j,k).lt.0.d0) hymax(i,j,k)=0.d0
	if (hzmax(i,j,k).lt.0.d0) hzmax(i,j,k)=0.d0

	if (hxmin(i,j,k).gt.1.d0) hxmin(i,j,k)=1.d0
	if (hymin(i,j,k).gt.1.d0) hymin(i,j,k)=1.d0
	if (hzmin(i,j,k).gt.1.d0) hzmin(i,j,k)=1.d0
	if (hxmax(i,j,k).gt.1.d0) hxmax(i,j,k)=1.d0
	if (hymax(i,j,k).gt.1.d0) hymax(i,j,k)=1.d0
	if (hzmax(i,j,k).gt.1.d0) hzmax(i,j,k)=1.d0




      hmintot=dmin1(hxmin(i,j,k),hxmax(i,j,k),hymin(i,j,k),hymax(i,j,k)
     %	,hzmin(i,j,k),hzmax(i,j,k))  

      hmaxtot=dmax1(hxmin(i,j,k),hxmax(i,j,k),hymin(i,j,k),hymax(i,j,k)
     %	,hzmin(i,j,k),hzmax(i,j,k))  



      if(hmintot.lt.0.d0.or.hmaxtot.gt.1.d0) then 
	write (*,*) "H problem" 
     & ,"i,j,k,hmintot,hmaxtot",i,j,k,hmintot,hmaxtot
	
	write (*,*) "nix,niy,niz",nix,niy,niz
	write (*,*)	"hx, hy,hz",hxmin(i,j,k),hxmax(i,j,k),hymin(i,j,k)
     %	,hymax(i,j,k),hzmin(i,j,k),hzmax(i,j,k)


      endif

cc	write (*,*) 'i,j,k,hxmin,hxmax', i,j,k,hxmin(i,j,k),hxmax(i,j,k)	
cc 	write (*,*) 'i,j,k,hymin,hymax', i,j,k,hymin(i,j,k),hymax(i,j,k)	
cc	write (*,*) 'i,j,k,hzmin,hzmax', i,j,k,hzmin(i,j,k),hzmax(i,j,k)	

   10 format(i3,1x,i3,1x,i3,1x,f7.3,2x,f7.3)

c	open (12,file='test.txt')
c	if (hxmax(i,j,k).lt.hxmin(i,j,k)) write (12,*) hxmax(i,j,k)
c     *,hxmin(i,j,k)


           enddo
         enddo
      enddo
c	close(12)
	return
	end

c	***********		exchange	***************


      subroutine exchange (a,b)		
      DOUBLE PRECISION a,b,temp

      temp=a
      a=b
      b=temp

      return
      end



c  ******** Calculation of Turbulent Eddy Viscosity  **********************
      SUBROUTINE eddyvisc (cc,ufilt)
      INCLUDE 'input.h'
      INTEGER i,j,k,j1,k1,j2
	REAL alfa,cmu1
	DOUBLE PRECISION FINDMAX,umin,umax,dudy,dudz
	DOUBLE PRECISION cc(nx,ny,nz),ufilt(nx,ny,nz)
      EXTERNAL FINDMAX
c  note that calculation of DELTA and value of alfa is for a round jet in x-direction.
      alfa=0.08
	cmu1=1.0
c  note the following calculation of mixing length is solly for jet flow with jet of 
c  liquid 1 being at some distance away from boundaries y and z. 
      call filter(ufilt)
	j1=y0/h
	k1=z0/h
      do i=1,nx
	if (xvel1.gt.xvel2) then
	umax=ufilt(i,j1,k1)
	umin=xvel2
         j2 = 0
	     do j=2,ny-1
	      if (ufilt(i,j,k1)-umin.gt.0.1*(umax-umin))then
             j2=j
             go to 10
		  endif   
     	     enddo
  10        continue
	else
	umax=xvel2
	umin=ufilt(i,j1,k1)
	     do j=2,ny-1
	      if (ufilt(i,j,k1)-umin.le.0.1*(umax-umin))then
             j2=j
             go to 12
		  endif   
           enddo
  12        continue
	endif
	if (j2.eq.2) j2=j1
      lmix(i)=alfa*h*abs(j2-j1)
      enddo
c  note: mut is calculated at the cell center (i,j,k) where cc(i,j,k) is calculated

	do i=2,nx-1
	do k=2,nz-1
	do j=2,ny-1
      dudy=0.25d0*(u(i,j,k)-u(i,j-1,k)+u(i+1,j,k)-u(i+1,j-1,k)
     &	        +u(i,j+1,k)-u(i,j,k)+u(i+1,j+1,k)-u(i+1,j,k))/h
	dudz=0.25d0*(u(i,j,k)-u(i,j,k-1)+u(i+1,j,k)-u(i+1,j,k-1)
     &	        +u(i,j,k+1)-u(i,j,k)+u(i+1,j,k+1)-u(i+1,j,k))/h
      mut(i,j,k)=cmu1*rhon(i,j,k)*0.25*(lmix(i)+lmix(i+1))*
     &	(lmix(i)+lmix(i+1))*dsqrt(dudy*dudy+dudz*dudz)
	enddo
	enddo
	enddo
	call bc_scalar(mut,ibc,nx,ny,nz)
      return
      end
c***********************************************************************
c ********* calculation of k-epsilon *******
	SUBROUTINE kepsmodel (t,kturnew,epsnew)
      INCLUDE 'input.h'
      INTEGER i,j,k,t
	DOUBLE PRECISION kturnew(nx,ny,nz),epsnew(nx,ny,nz)
	DOUBLE PRECISION conv,diff,prod,sigk,sige,c1,c2,h2,h22
	DOUBLE PRECISION dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
	DOUBLE PRECISION cip,cim,cjp,cjm,ckp,ckm,muim,muip,mujm,mujp
      sigk=1.0d0
	sige=1.3d0
	c1=1.44d0
	c2=1.92d0

	h2=h*h
	h22=h2*2.d0
      do j=2,ny-1
        do i=2,nx-1
          do k=2,nz-1
c  ****** k-equation
c  ****** convective term  ***********************
	  conv=(
     *   +dmax1(u(i+1,j,k),0.0d0)*ktur(i,j,k)
     *   +dmin1(u(i+1,j,k),0.0d0)*ktur(i+1,j,k)
     *   -dmax1(u(i,j,k),0.0d0)*ktur(i-1,j,k)
     *   -dmin1(u(i,j,k),0.0d0)*ktur(i,j,k)
     *   +dmax1(v(i,j+1,k),0.0d0)*ktur(i,j,k)
     *   +dmin1(v(i,j+1,k),0.0d0)*ktur(i,j+1,k)
     *   -dmax1(v(i,j,k),0.0d0)*ktur(i,j-1,k)
     *   -dmin1(v(i,j,k),0.0d0)*ktur(i,j,k)
     *   +dmax1(w(i,j,k+1),0.0d0)*ktur(i,j,k)
     *   +dmin1(w(i,j,k+1),0.0d0)*ktur(i,j,k+1)
     *   -dmax1(w(i,j,k),0.d0)*ktur(i,j,k-1)
     *   -dmin1(w(i,j,k),0.d0)*ktur(i,j,k))/h

c  **** diffusion term ****************
	    cip=mun(i+1,j,k)+mun(i,j,k)+(mut(i,j,k)+mut(i+1,j,k))/sigk
	    cim=mun(i-1,j,k)+mun(i,j,k)+(mut(i,j,k)+mut(i-1,j,k))/sigk
	    cjp=mun(i,j+1,k)+mun(i,j,k)+(mut(i,j,k)+mut(i,j+1,k))/sigk
	    cjm=mun(i,j-1,k)+mun(i,j,k)+(mut(i,j,k)+mut(i,j-1,k))/sigk
	    ckp=mun(i,j,k+1)+mun(i,j,k)+(mut(i,j,k)+mut(i,j,k+1))/sigk
	    ckm=mun(i,j,k-1)+mun(i,j,k)+(mut(i,j,k)+mut(i,j,k-1))/sigk
	   
	   diff=1.d0/(rhon(i,j,k)*h22)*(
     *        cip*(ktur(i+1,j,k)-ktur(i,j,k))
     *	   -cim*(ktur(i,j,k)-ktur(i-1,j,k))
     *	   +cjp*(ktur(i,j+1,k)-ktur(i,j,k))
     *	   -cjm*(ktur(i,j,k)-ktur(i,j-1,k))
     *	   +ckp*(ktur(i,j,k+1)-ktur(i,j,k))
     *	   -ckm*(ktur(i,j,k)-ktur(i,j,k-1)))
c  **** production term  **************
		dudx=u(i+1,j,k)-u(i,j,k)
          dudy=0.25d0*(u(i,j,k)-u(i,j-1,k)+u(i+1,j,k)-u(i+1,j-1,k)
     &	        +u(i,j+1,k)-u(i,j,k)+u(i+1,j+1,k)-u(i+1,j,k))
		dudz=0.25d0*(u(i,j,k)-u(i,j,k-1)+u(i+1,j,k)-u(i+1,j,k-1)
     &	        +u(i,j,k+1)-u(i,j,k)+u(i+1,j,k+1)-u(i+1,j,k))
		dvdx=0.25d0*(v(i,j,k)-v(i-1,j,k)+v(i+1,j,k)-v(i,j,k)
     &        	+v(i,j+1,k)-v(i-1,j+1,k)+v(i+1,j+1,k)-v(i,j+1,k))
		dvdy=v(i,j+1,k)-v(i,j,k)
		dvdz=0.25d0*(v(i,j,k)-v(i,j,k-1)+v(i,j,k+1)-v(i,j,k)
     &        	+v(i,j+1,k)-v(i,j+1,k-1)+v(i,j+1,k+1)-v(i,j+1,k))
		dwdx=0.25d0*(w(i,j,k)-w(i-1,j,k)+w(i+1,j,k)-w(i,j,k)
     &        	+w(i,j,k+1)-w(i-1,j,k+1)+w(i+1,j,k+1)-w(i,j,k+1))
		dwdy=0.25d0*(w(i,j,k)-w(i,j-1,k)+w(i,j,k+1)-w(i,j-1,k+1)
     &	        +w(i,j+1,k)-w(i,j,k)+w(i,j+1,k+1)-w(i,j,k+1))
		dwdz=w(i,j,k+1)-w(i,j,k)



		     if(t.eq.0.or.iturb.eq.3) then

		prod=mut(i,j,k)/rhon(i,j,k)/h2*(2.d0*
     &	(dudx*dudx+dvdy*dvdy+dwdz*dwdz+dvdx*dudy+dwdx*dudz+dvdz*dwdy)+
     &	dudy*dudy+dudz*dudz+dvdx*dvdx+dvdz*dvdz+dwdx*dwdx+dwdy*dwdy)

               elseif (iturb.eq.4) then
               prod=-(u2p(i,j,k)*dudx+v2p(i,j,k)*dvdy+w2p(i,j,k)*dwdz
     &        +uvp(i,j,k)*(dudy+dvdx)+uwp(i,j,k)*(dudz+dwdx)
     &        +vwp(i,j,k)*(dvdz+dwdy))/h

               else 
			 write (*,*) 'iturb problem'
			 stop  

	         endif


     		kturnew(i,j,k)=ktur(i,j,k)+tau*(-conv+diff+prod-eps(i,j,k))

c  ***** epsilon Equation  *******************
c  ****** convective term  ***********************
	   conv=(
     *    +dmax1(u(i+1,j,k),0.0d0)*eps(i,j,k)
     *    +dmin1(u(i+1,j,k),0.0d0)*eps(i+1,j,k)
     *	-dmax1(u(i,j,k),0.0d0)*eps(i-1,j,k)
     *    -dmin1(u(i,j,k),0.0d0)*eps(i,j,k)
     *    +dmax1(v(i,j+1,k),0.0d0)*eps(i,j,k)
     *    +dmin1(v(i,j+1,k),0.0d0)*eps(i,j+1,k)
     *    -dmax1(v(i,j,k),0.0d0)*eps(i,j-1,k)
     *    -dmin1(v(i,j,k),0.d0)*eps(i,j,k)
     *    +dmax1(w(i,j,k+1),0.0d0)*eps(i,j,k)
     *    +dmin1(w(i,j,k+1),0.0d0)*eps(i,j,k+1)
     *    -dmax1(w(i,j,k),0.d0)*eps(i,j,k-1)
     *    -dmin1(w(i,j,k),0.d0)*eps(i,j,k))/h

c  **** diffusion term ****************
	    cip=mun(i+1,j,k)+mun(i,j,k)+(mut(i,j,k)+mut(i+1,j,k))/sige
	    cim=mun(i-1,j,k)+mun(i,j,k)+(mut(i,j,k)+mut(i-1,j,k))/sige
	    cjp=mun(i,j+1,k)+mun(i,j,k)+(mut(i,j,k)+mut(i,j+1,k))/sige
	    cjm=mun(i,j-1,k)+mun(i,j,k)+(mut(i,j,k)+mut(i,j-1,k))/sige
	    ckp=mun(i,j,k+1)+mun(i,j,k)+(mut(i,j,k)+mut(i,j,k+1))/sige
	    ckm=mun(i,j,k-1)+mun(i,j,k)+(mut(i,j,k)+mut(i,j,k-1))/sige
	   diff=1.d0/(rhon(i,j,k)*h22)*(
     *	cip*(eps(i+1,j,k)-eps(i,j,k))-cim*(eps(i,j,k)-eps(i-1,j,k))+
     *	cjp*(eps(i,j+1,k)-eps(i,j,k))-cjm*(eps(i,j,k)-eps(i,j-1,k))+
     *	ckp*(eps(i,j,k+1)-eps(i,j,k))-ckm*(eps(i,j,k)-eps(i,j,k-1)))
c  **** production term  **************
		prod=prod*c1*eps(i,j,k)/ktur(i,j,k)
		epsnew(i,j,k) =eps(i,j,k)+tau*
     &		(-conv+diff+prod-c2*eps(i,j,k)*eps(i,j,k)/ktur(i,j,k))
	   enddo
	 enddo
      enddo
      do j=2,ny-1
       do i=2,nx-1
         do k=2,nz-1
		ktur(i,j,k)=kturnew(i,j,k)
		eps(i,j,k)=epsnew(i,j,k)
          mut(i,j,k)=cmu(i,j,k)*rhon(i,j,k)*ktur(i,j,k)*ktur(i,j,k)
     &		/eps(i,j,k)
	   enddo
	 enddo
      enddo
      call bc_scalar(ktur,ibc,nx,ny,nz)
	call bc_scalar(eps,ibc,nx,ny,nz)
	call bc_scalar(mut,ibc,nx,ny,nz)

      	if(iturb.eq.4) call uv2p ()	

	return
      end



c ********* calculation of Reynolds stress terms, (SHIH ZHU and LUMLEY Model) *******

      SUBROUTINE uv2p ()

      include 'input.h'
      INTEGER i,j,k
	DOUBLE PRECISION ss11,ss12,ss13,ss22,ss23,ss33,sskk
	DOUBLE PRECISION dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz	
      DOUBLE PRECISION sstar11,sstar22,sstar33,omeg12,omeg13,omeg23
      DOUBLE PRECISION wstar,ustar,omegastar,sstar,sixsqws
	DOUBLE PRECISION sixsq,Asstar,c2,keps,k2eps,k2eps2,k3eps2

      do k=2,nz-1
        do j=2,ny-1
          do i=2,nx-1
       	keps=ktur(i,j,k)/eps(i,j,k)
		k2eps=keps*ktur(i,j,k)
		k2eps2=keps*keps
		k3eps2=k2eps2*ktur(i,j,k)
		
		dudx=(u(i+1,j,k)-u(i,j,k))/h
                dudy=0.25d0*(u(i,j,k)-u(i,j-1,k)+u(i+1,j,k)-u(i+1,j-1,k)
     &	        +u(i,j+1,k)-u(i,j,k)+u(i+1,j+1,k)-u(i+1,j,k))/h
		dudz=0.25d0*(u(i,j,k)-u(i,j,k-1)+u(i+1,j,k)-u(i+1,j,k-1)
     &	        +u(i,j,k+1)-u(i,j,k)+u(i+1,j,k+1)-u(i+1,j,k))/h
		dvdx=0.25d0*(v(i,j,k)-v(i-1,j,k)+v(i+1,j,k)-v(i,j,k)
     &        	+v(i,j+1,k)-v(i-1,j+1,k)+v(i+1,j+1,k)-v(i,j+1,k))/h
		dvdy=(v(i,j+1,k)-v(i,j,k))/h
		dvdz=0.25d0*(v(i,j,k)-v(i,j,k-1)+v(i,j,k+1)-v(i,j,k)
     &        	+v(i,j+1,k)-v(i,j+1,k-1)+v(i,j+1,k+1)-v(i,j+1,k))/h
		dwdx=0.25d0*(w(i,j,k)-w(i-1,j,k)+w(i+1,j,k)-w(i,j,k)
     &        	+w(i,j,k+1)-w(i-1,j,k+1)+w(i+1,j,k+1)-w(i,j,k+1))/h
		dwdy=0.25d0*(w(i,j,k)-w(i,j-1,k)+w(i,j,k+1)-w(i,j-1,k+1)
     &	        +w(i,j+1,k)-w(i,j,k)+w(i,j+1,k+1)-w(i,j,k+1))/h
		dwdz=(w(i,j,k+1)-w(i,j,k))/h
		
		ss11=dudx
		ss12=0.5*(dudy+dvdx)
		ss13=0.5*(dudz+dwdx)
		ss22=dvdy
		ss23=0.5*(dvdz+dwdy)		
		ss33=dwdz
		
		sskk=ss11+ss22+ss33
		sstar11=ss11-1./3.*sskk
		sstar22=ss22-1./3.*sskk
                sstar33=ss33-1./3.*sskk
		sstar=dsqrt(sstar11*sstar11+sstar22*sstar22+sstar33*sstar33
     *		        +2.*ss12*ss12+2.*ss13*ss13+2.*ss23*ss23+1.E-30)
	
            omeg12=0.5*(dudy-dvdx)
            omeg13=0.5*(dudz-dwdx)
            omeg23=0.5*(dvdz-dwdy)
            omegastar=dsqrt(2.*omeg12*omeg12+
     *      2.*omeg13*omeg13+2.*omeg23*omeg23)

		wstar=(sstar11*sstar11*sstar11+sstar22*sstar22*sstar22+
     *		sstar33*sstar33*sstar33+3.*sstar11*ss12*ss12+
     *		3.*ss12*ss12*sstar22+3.*sstar11*ss13*ss13+
     *		6.*ss12*ss13*ss23+3.*ss13*ss13*sstar33+
     *		3.*sstar22*ss23*ss23+3.*ss23*ss23*sstar33)/sstar/sstar/sstar

		ustar=dsqrt(sstar*sstar+omegastar*omegastar)

		sixsq=sqrt(6.)
		sixsqws=sixsq*wstar
		
		if((1-dabs(sixsqws)).lt.1.d-4) then
		
		    if (sixsqws.lt.0) then
		    sixsqws=-1.
		    else
		    sixsqws=1.
		    endif
				
		endif
		
		Asstar=sixsq*dcos(dacos(sixsqws)/3.)
          cmu(i,j,k)=1.0/(6.5+asstar*ustar*keps)

          mut(i,j,k)=cmu(i,j,k)*rhon(i,j,k)*ktur(i,j,k)*ktur(i,j,k)
     &		/eps(i,j,k)

		c2=1.-9.*cmu(i,j,k)*cmu(i,j,k)*sstar*sstar*k2eps2
		if (c2.lt.0.)c2=0.
		c2=dsqrt(c2)/(1.+6.*sstar*k2eps2*omegastar)

	    u2p(i,j,k)=2./3.*ktur(i,j,k)-2.*cmu(i,j,k)*k2eps*sstar11+
     &            4.*c2*k3eps2*(ss12*omeg12+ss13*omeg13)

	    v2p(i,j,k)=2./3.*ktur(i,j,k)-2.*cmu(i,j,k)*k2eps*sstar22-
     &            4.*c2*k3eps2*(ss12*omeg12-ss23*omeg23)

            w2p(i,j,k)=2./3.*ktur(i,j,k)-2.*cmu(i,j,k)*k2eps*sstar33-
     &            4.*c2*k3eps2*(ss13*omeg13+ss23*omeg23)

            uvp(i,j,k)=-2.*cmu(i,j,k)*k2eps*ss12+2.*c2*k3eps2*
     &       (-sstar11*omeg12+sstar22*omeg12+ss13*omeg23+ss23*omeg13)
     			
            uwp(i,j,k)=-2.*cmu(i,j,k)*k2eps*ss13+2.*c2*k3eps2*
     &       (-sstar11*omeg13+sstar33*omeg13-ss12*omeg23+ss23*omeg12)

            vwp(i,j,k)=-2.*cmu(i,j,k)*k2eps*ss23+2.*c2*k3eps2*
     &       (-sstar22*omeg23+sstar33*omeg23-ss12*omeg13-ss13*omeg12)

		  enddo
	    enddo
	enddo	

      call bc_scalar(mut,ibc,nx,ny,nz)

      call bc_scalar(u2p,ibc,nx,ny,nz)
	call bc_scalar(v2p,ibc,nx,ny,nz)
      call bc_scalar(w2p,ibc,nx,ny,nz)
	call bc_scalar(uvp,ibc,nx,ny,nz)
	call bc_scalar(uwp,ibc,nx,ny,nz)
	call bc_scalar(vwp,ibc,nx,ny,nz)
	
	return
      end




c***********************************************************************

	SUBROUTINE checktime(t)
	INCLUDE 'input.h'
	INTEGER t
	DOUBLE PRECISION CFL,cfl1
	cfl1=CFL()
	if(cfl1.gt.MAXCFL)then
	tau=tau*(MAXCFL*0.8)/cfl1
	tauh=tau/h
      tauh2 = tau/(h*h)
	endif
	if(cfl1.le.MAXCFL*0.25)then
	cfl1=maxcfl*0.5	!new because when CFL is very small tau become very big
	write (*,*) "*****NEW CFL*****"
	tau=tau*(MAXCFL*0.5)/cfl1
	tauh=tau/h
      tauh2 = tau/(h*h)
	endif
	time=time+tau
	return
      end
c***********************************************************************
      SUBROUTINE calcrho (cc)
      INCLUDE 'input.h'
      INTEGER i,j,k
	DOUBLE PRECISION cc(nx,ny,nz)
	do i=1,nx
		do k=1,nz
			do j=1,ny
			rhon(i,j,k)=rho2+cc(i,j,k)*(rho1-rho2)
			mun(i,j,k)=mu2+cc(i,j,k)*(mu1-mu2)
			enddo
		enddo
	enddo
      return
      end

c***********************************************************************
c     This subroutine calculates the part of the stress tensor T due to
c     surface tension. Temporary values are calculated in the center
c     of the cell shown below, then the appropriate averages are
c     performed. The equation for T is : T = - sigh*|dc|*(I - n@n) 
c        
c                            ON OUTPUT:
c
c     S11, S22, S33, S12, S13, S23:        Sij 
c             
c***   
c
c                                       u               cc,p,S11
c             p(i-1) *------------------*------------------*
c                   /.                                    /|
c                  / .                S12                / |
c          v(i-1) *  .               *                v *  |
c                /   .                                 /   |
c               /    .       u(j-1)            p(j-1) /    |
c   p(i-1,j-1) *------------------*------------------*     |
c              |     .                 S13           |     |
c              |     + w(i-1)           +            |     * w
c              |     .                               |     |
c              |     .                               | S23 |
c     S23(i-1) |  +  .                               |  *  |
c              |     .                               |     |
c              |     .         S13(j-1)              |     |
c   w(i-1,j-1) *     .            *           w(j-1) *     |
c              |     .                               |     |
c              |     + . . . . . . . . .+ . . . . . .| . . * p(k-1) 
c              |    . p(i-1,k-1)      u(k-1)         |    /         
c              |   .                                 |   /         
c              |  + v(i-1,k-1)       +               |  * v(k-1)   
c              | .                  S12(k-1)         | /
c              |.                                    |/
c              *------------------*------------------* p(j-1,k-1)
c     p(i-1,j-1,k-1)          u(j-1,k-1)
c
c
c***********************************************************************
c***
      SUBROUTINE surfacetens( cc,S11t,S22t,S33t,S12t,S13t,S23t,mag)
c***
      INCLUDE 'input.h'
      INTEGER i, j, k
      DOUBLE PRECISION cc(nx,ny,nz), mag(nx,ny,nz),deltasctr(nx,ny,nz)
      DOUBLE PRECISION S12t(nx,ny,nz),S13t(nx,ny,nz),S23t(nx,ny,nz)
      DOUBLE PRECISION S11t(nx,ny,nz),S22t(nx,ny,nz),S33t(nx,ny,nz)
      DOUBLE PRECISION sigh,MINGRAD,diff,deltasij,ftr,mu
      PARAMETER (MINGRAD=1.0d-50)
      INTRINSIC DSQRT,DABS
c***
      ftr=0. ! (default is 1)if it is zero, no turbulence in surface tension 
      if (sigma .eq. 0.d0.and.iturb.eq.0) then
         call fillzero(S11,nx,ny,nz)
         call fillzero(S22,nx,ny,nz)
         call fillzero(S33,nx,ny,nz)
         call fillzero(S12,nx,ny,nz)
         call fillzero(S13,nx,ny,nz)
         call fillzero(S23,nx,ny,nz)
      else
      sigh = .25d0*tau**2*sigma/h**3
c***
c     calculate the components of the gradient of cc in the center
c     of the above cell. (Stored in S11,S22,S33 to save memory space)
c***
      do k=2,nz
         do j=2,ny
            do i=2,nx
               S11(i,j,k) = (cc(i,j,k) + cc(i,j-1,k) + cc(i,j,k-1)
     %                + cc(i,j-1,k-1)) - (cc(i-1,j,k) + cc(i-1,j-1,k)
     %                + cc(i-1,j,k-1) + cc (i-1,j-1,k-1))
               S22(i,j,k) = (cc(i,j,k) + cc(i-1,j,k) + cc(i,j,k-1)
     %                + cc(i-1,j,k-1)) - (cc(i,j-1,k) + cc(i-1,j-1,k)
     %                + cc(i,j-1,k-1) + cc(i-1,j-1,k-1))
               S33(i,j,k) = (cc(i,j,k) + cc(i,j-1,k) + cc(i-1,j,k)
     %                + cc(i-1,j-1,k)) - (cc(i,j,k-1) + cc(i,j-1,k-1)
     %                + cc(i-1,j,k-1) + cc(i-1,j-1,k-1))
            enddo 
         enddo 
      enddo 
c***
c     calculate the unit normal vector n and the modulus of grad cc as:
c     mag  = |dc| = MAX(|dc|,MINGRAD), to avoid division by zero.
c***
      do k=2,nz
         do j=2,ny
            do i=2,nx
               mag(i,j,k) = DSQRT(S11(i,j,k)**2 + S22(i,j,k)**2 +
     %              S33(i,j,k)**2)
               diff = MINGRAD - mag(i,j,k)
               mag(i,j,k) = mag(i,j,k) + 0.5d0*(DABS(diff) + diff)
               deltasctr(i,j,k)=0.25d0/h*mag(i,j,k)
               S11(i,j,k) = S11(i,j,k)/mag(i,j,k)
               S22(i,j,k) = S22(i,j,k)/mag(i,j,k)
               S33(i,j,k) = S33(i,j,k)/mag(i,j,k)
               mag(i,j,k) = sigh*mag(i,j,k)
            enddo 
         enddo 
      enddo 
c***
c     get temporary stress tensor coefficients: S11t, S22t, S33t,
c     S12t, S13t, S23t, from: T = - sigh*|dc|*(I - n@n) 
c***
      do k=2,nz
         do j=2,ny
            do i=2,nx
               S11t(i,j,k)  = - mag(i,j,k)*
     %              (1.d0 - S11(i,j,k)*S11(i,j,k))
               S22t(i,j,k)  = - mag(i,j,k)*
     %              (1.d0 - S22(i,j,k)*S22(i,j,k))
               S33t(i,j,k)  = - mag(i,j,k)*
     %              (1.d0 - S33(i,j,k)*S33(i,j,k))
               S12t(i,j,k) =  mag(i,j,k)*
     %                    S11(i,j,k)*S22(i,j,k)
               S13t(i,j,k) =  mag(i,j,k)*
     %                    S11(i,j,k)*S33(i,j,k)
               S23t(i,j,k) =  mag(i,j,k)*
     %                    S22(i,j,k)*S33(i,j,k)
            enddo
         enddo
      enddo
c***
c     stress tensor coefficients in the vertices of the cell
c     are given as averages of the surrounding 8 points
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               S11(i,j,k) = 0.125d0*(S11t(i,j,k) + S11t(i+1,j,k)
     %              +S11t(i,j+1,k)+S11t(i+1,j+1,k)+S11t(i,j,k+1)
     %              +S11t(i+1,j,k+1)+S11t(i,j+1,k+1)+S11t(i+1,j+1,k+1))
               S22(i,j,k) = 0.125d0*(S22t(i,j,k) + S22t(i+1,j,k)
     %              +S22t(i,j+1,k)+S22t(i+1,j+1,k)+S22t(i,j,k+1)
     %              +S22t(i+1,j,k+1)+S22t(i,j+1,k+1)+S22t(i+1,j+1,k+1))
               S33(i,j,k) = 0.125d0*(S33t(i,j,k) + S33t(i+1,j,k)
     %              +S33t(i,j+1,k)+S33t(i+1,j+1,k)+S33t(i,j,k+1)
     %              +S33t(i+1,j,k+1)+S33t(i,j+1,k+1)+S33t(i+1,j+1,k+1))
c
        if(iturb.eq.2.or.iturb.eq.3.or.iturb.eq.4) then
        	deltasij=0.125d0*(deltasctr(i,j,k)+deltasctr(i+1,j,k)
     &          +deltasctr(i,j+1,k)+deltasctr(i+1,j+1,k)
     &		  +deltasctr(i,j,k+1)+deltasctr(i+1,j,k+1)
     &          +deltasctr(i,j+1,k+1)+deltasctr(i+1,j+1,k+1))
	     mu=mu2+cc(i,j,k)*( mu1 - mu2)
        turbst(i,j,k)=ftr*sigma*(dsqrt(0.1d0*cmu(i,j,k)*mut(i,j,k)/mu))
     &		 *tauh*deltasij
        endif
            enddo
         enddo
      enddo
      if(iturb.eq.2.or.iturb.eq.3.or.iturb.eq.4) 
     &	call bc_scalar(turbst,ibc,nx,ny,nz)
c***
c     stress tensor coefficients in the faces of the cell are
c     given as averages along the normal to the surface
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               S12(i,j,k) =  0.5d0*(S12t(i,j,k) + S12t(i,j,k+1))
               S13(i,j,k) =  0.5d0*(S13t(i,j,k) + S13t(i,j+1,k))
               S23(i,j,k) =  0.5d0*(S23t(i,j,k) + S23t(i+1,j,k))
            enddo
         enddo
      enddo
      endif
c***
      return
      end
c     ************csf*************	
      SUBROUTINE csf(ccc,w11,w22,w33)
      include 'input.h'
      DOUBLE PRECISION sigh,d1x,d2x,d1y,d2y,d1z,d2z
      DOUBLE PRECISION nxu,nxd,nyu,nyd,nzu,nzd
      DOUBLE PRECISION ccf(nx,ny,nz),ccc(nx,ny,nz)
c      DOUBLE PRECISION S12t(nx,ny,nz),S13t(nx,ny,nz),S23t(nx,ny,nz)
      DOUBLE PRECISION w11(nx,ny,nz),w22(nx,ny,nz),w33(nx,ny,nz)

	INTEGER i, j, k,nf
	nf=1	

	call copy(ccf,ccc,nx,ny,nz)
	do i=1,nf
	call filterc27(ccf,w22,nx,ny,nz)
	call copy (ccf,w22,nx,ny,nz)
	call bc_out(ccf,nx,ny,nz)
	enddo
c	if (nf.ge.1)then	 
c	call filterc27(ccc,ccf,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	do i=1,nf-1
c	call copy (w22,ccf,nx,ny,nz)
c	call filterc27(w22,ccf,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	enddo
c	else
c	call copy(ccf,ccc,nx,ny,nz)
c	endif

	if(sigma.eq.0.d0) then
         call fillzero(S11,nx,ny,nz)
         call fillzero(S22,nx,ny,nz)
	   call fillzero(S33,nx,ny,nz)
	   call fillzero(S12,nx,ny,nz)
	   call fillzero(S13,nx,ny,nz)			
         call fillzero(S23,nx,ny,nz)
      else
      sigh =0.5*tau**2*sigma/h**3

      do j=2,ny-1
         do i=2,nx-1
			do k=2,nz-1



            nxu=dsqrt((ccf(i+1,j,k)-ccf(i,j,k))**2
     *          +0.0625d0*(ccf(i,j+1,k)-ccf(i,j-1,k)+ccf(i+1,j+1,k)
     *          -ccf(i+1,j-1,k))**2+0.0625d0*(ccf(i,j,k+1)
     *          -ccf(i,j,k-1)+ccf(i+1,j,k+1)-ccf(i+1,j,k-1))**2+1.d-50)

            nxd=dsqrt((ccf(i,j,k)-ccf(i-1,j,k))**2
     *          +0.0625d0*(ccf(i,j+1,k)-ccf(i,j-1,k)+ccf(i-1,j+1,k)
     *          -ccf(i-1,j-1,k))**2+0.0625d0*(ccf(i,j,k+1)
     *          -ccf(i,j,k-1)+ccf(i-1,j,k+1)-ccf(i-1,j,k-1))**2+1.d-50)

            nyu=dsqrt((ccf(i,j+1,k)-ccf(i,j,k))**2
     *          +0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j+1,k)
     *          -ccf(i-1,j+1,k))**2+0.0625d0*(ccf(i,j,k+1)
     *          -ccf(i,j,k-1)+ccf(i,j+1,k+1)-ccf(i,j+1,k-1))**2+1.d-50)

            nyd=dsqrt((ccf(i,j,k)-ccf(i,j-1,k))**2
     *          +0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j-1,k)
     *          -ccf(i-1,j-1,k))**2+0.0625d0*(ccf(i,j,k+1)
     *          -ccf(i,j,k-1)+ccf(i,j-1,k+1)-ccf(i,j-1,k-1))**2+1.d-50)

            nzu=dsqrt((ccf(i,j,k+1)-ccf(i,j,k))**2
     *           +0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j,k+1)
     *           -ccf(i-1,j,k+1))**2+0.0625d0*(ccf(i,j+1,k)
     *           -ccf(i,j-1,k)+ccf(i,j+1,k+1)-ccf(i,j-1,k+1))**2+1.d-50)

            nzd=dsqrt((ccf(i,j,k)-ccf(i,j,k-1))**2
     *           +0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j,k-1)
     *           -ccf(i-1,j,k-1))**2+0.0625d0*(ccf(i,j+1,k)
     *           -ccf(i,j-1,k)+ccf(i,j+1,k-1)-ccf(i,j-1,k-1))**2+1.d-50)

            d1x=(ccf(i+1,j,k)-ccf(i,j,k))/nxu
            d2x=(ccf(i,j,k)-ccf(i-1,j,k))/nxd
            d1y=(ccf(i,j+1,k)-ccf(i,j,k))/nyu
            d2y=(ccf(i,j,k)-ccf(i,j-1,k))/nyd
            d1z=(ccf(i,j,k+1)-ccf(i,j,k))/nzu
            d2z=(ccf(i,j,k)-ccf(i,j,k-1))/nzd
c	here we temporary use w11 instead of divergence
            w11(i,j,k) = - (d1x-d2x+d1y-d2y+d1z-d2z)

			enddo
		enddo
	enddo
	!!!call bcp(ccf,nx,ny)
	call bc_out(w11,nx,ny,nz)

	do	k=2,nz-1
        do j=2,ny-1
	     do i=2,nx-1

	s11(i,j,k)=sigh*(ccf(i,j,k)-ccf(i-1,j,k))*(w11(i,j,k)
     *)
      s22(i,j,k)=sigh*(ccf(i,j,k)-ccf(i,j-1,k))*(w11(i,j,k)
     *)
	s33(i,j,k)=sigh*(ccf(i,j,k)-ccf(i,j,k-1))*(w11(i,j,k)
     *)
     	

		 enddo
	    enddo
	enddo

	
	call copy (w11,s11,nx,ny,nz)
	call copy (w22,s22,nx,ny,nz)
	call copy (w33,s33,nx,ny,nz)	   

c      !???  calls boundary conditions for a homogenous field

       !???  call bcvechom(div,cc,nx,ny) !=bcv for homogenous condition
c	call bc_w(w11,w22,w33)
      call bc_tensor()  
	endif
      return
      end

c     ************csfbkz*************
c   Our implementation of the original Brackbill-Kothe-Zemach
c   algorithm
c
c$$$@ARTICLE{brackbill92,
c$$$	author="J.U. Brackbill  and  D. B. Kothe and C. Zemach ",
c$$$        title="A Continuum Method for Modeling Surface Tension",
c$$$        journal="J. Comput. Phys.",
c$$$        year="1992",
c$$$        volume="100",
c$$$        pages="335-354"}
c$$$
c
c
c   On input:
c
c        cc          (filtered) volume fractions
c        S1         arbitrary nx*ny array;
c        S2              "
c        w1             "
c        w2
c        sigma     surface tension coefficient.
c
c
c   On output
c
c        w2 first component of force multiplied by tau/h
c        cc  second "   "  "  "  "  "
c        in momentum.f  at locations described on figure.
c	
      SUBROUTINE csfbkz(ccc,w11,w22,w33)
      include 'input.h'
      DOUBLE PRECISION sigh,d1x,d2x,d1y,d2y,d1z,d2z
      DOUBLE PRECISION nxu,nxd,nyu,nyd,nzu,nzd
      DOUBLE PRECISION ccf(nx,ny,nz),ccc(nx,ny,nz),rho
c      DOUBLE PRECISION S12t(nx,ny,nz),S13t(nx,ny,nz),S23t(nx,ny,nz)
      DOUBLE PRECISION w11(nx,ny,nz),w22(nx,ny,nz),w33(nx,ny,nz)
      DOUBLE PRECISION theta,mx,my,mz,mmx,mmy,mmz,normv
	INTEGER i, j, k,nf
	nf=1
	
	call copy(ccf,ccc,nx,ny,nz)
	do i=1,nf
	call filterc27(ccf,w22,nx,ny,nz)
	call copy (ccf,w22,nx,ny,nz)
	call bc_out(ccf,nx,ny,nz)
	enddo
		
c	if (nf.ge.1)then	 
c	call filterc27(ccc,ccf,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	do i=1,nf-1
c	call copy (w22,ccf,nx,ny,nz)
c	call filterc27(w22,ccf,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	enddo
c	else
c	call copy(ccf,ccc,nx,ny,nz)
c	endif

	if(sigma.eq.0.d0) then
         call fillzero(S11,nx,ny,nz)
         call fillzero(S22,nx,ny,nz)
	   call fillzero(S33,nx,ny,nz)
	   call fillzero(S12,nx,ny,nz)
	   call fillzero(S13,nx,ny,nz)			
         call fillzero(S23,nx,ny,nz)
      else
      sigh =0.5*tau**2*sigma/h**3

      do j=2,ny-1
         do i=2,nx-1
			do k=2,nz-1



      nxu=dsqrt((ccf(i+1,j,k)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i,j+1,k)-ccf(i,j-1,k)+ccf(i+1,j+1,k)-ccf(i+1,j-1,k
     *))**2
     *+0.0625d0*(ccf(i,j,k+1)-ccf(i,j,k-1)+ccf(i+1,j,k+1)-ccf(i+1,j,k-1
     &))**2
     *+1.d-50)
      nxd=dsqrt((ccf(i,j,k)-ccf(i-1,j,k))**2
     *+0.0625d0*(ccf(i-1,j+1,k)-ccf(i-1,j-1,k)+ccf(i,j+1,k)-ccf(i,j-1,k
     &))**2
     *+0.0625d0*(ccf(i-1,j,k+1)-ccf(i-1,j,k-1)+ccf(i,j,k+1)-ccf(i,j,k-1
     &))**2
     *+1.d-50)
      nyu=dsqrt((ccf(i,j+1,k)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i+1,j+1,k)-ccf(i-1,j+1,k)+ccf(i+1,j,k)-ccf(i-1,j,k
     &))**2
     *+0.0625d0*(ccf(i,j+1,k+1)-ccf(i,j+1,k-1)+ccf(i,j,k+1)-ccf(i,j,k-1
     &))**2
     *+1.d-50)
      nyd=dsqrt((ccf(i,j,k)-ccf(i,j-1,k))**2
     *+0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j-1,k)-ccf(i-1,j-1,k
     &))**2
     *+0.0625d0*(ccf(i,j,k+1)-ccf(i,j,k-1)+ccf(i,j-1,k+1)-ccf(i,j-1,k-1
     &))**2
     *+1.d-50)
      nzu=dsqrt((ccf(i,j,k+1)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j,k+1)-ccf(i-1,j,k+1
     &))**2
     *+0.0625d0*(ccf(i,j+1,k)-ccf(i,j-1,k)+ccf(i,j+1,k+1)-ccf(i,j-1,k+1
     &))**2
     *+1.d-50)
      nzd=dsqrt((ccf(i,j,k)-ccf(i,j,k-1))**2
     *+0.0625d0*(ccf(i+1,j,k-1)-ccf(i-1,j,k-1)+ccf(i+1,j,k)-ccf(i-1,j,k
     &))**2
     *+0.0625d0*(ccf(i,j+1,k-1)-ccf(i,j-1,k-1)+ccf(i,j+1,k)-ccf(i,j-1,k
     &))**2
     *+1.d-50)

            d1x=(ccf(i+1,j,k)-ccf(i,j,k))/nxu
            d2x=(ccf(i,j,k)-ccf(i-1,j,k))/nxd
            d1y=(ccf(i,j+1,k)-ccf(i,j,k))/nyu
            d2y=(ccf(i,j,k)-ccf(i,j-1,k))/nyd
            d1z=(ccf(i,j,k+1)-ccf(i,j,k))/nzu
            d2z=(ccf(i,j,k)-ccf(i,j,k-1))/nzd
c	here we temporary use w11 instead of divergence
            w11(i,j,k) = - (d1x-d2x+d1y-d2y+d1z-d2z)

			enddo
		enddo
	enddo
	!!!call bcp(ccf,nx,ny)
	call bc_out(w11,nx,ny,nz)

	do	k=2,nz-1
        do j=2,ny-1
	     do i=2,nx-1

	s11(i,j,k)=sigh*(ccf(i,j,k)-ccf(i-1,j,k))*w11(i,j,k)
      s22(i,j,k)=sigh*(ccf(i,j,k)-ccf(i,j-1,k))*w11(i,j,k)
	s33(i,j,k)=sigh*(ccf(i,j,k)-ccf(i,j,k-1))*w11(i,j,k)

c	**************************	new  **********************
     	if (j.eq.2) then

	mx=(ccf(i+1,j,k)-ccf(i-1,j,k))
	mz=(ccf(i,j,k+1)-ccf(i,j,k-1))
	mmx=mx/dsqrt(mz**2+mx**2+1.d-50)
	mmz=mz/dsqrt(mz**2+mx**2+1.d-50)
	normv=mmx*u(i,j,k)+mmz*w(i,j,k)
	normv=-normv
	if (dabs(normv).gt.0.1) then
		if (normv.gt.0)	theta=115
		if	(normv.lt.0) theta=45
	else
	theta=350.0*normv+80


c	**** determination of theta as a function of Vcl *****
c	sign=(ccf(i+1,j,k)-ccf(i-1,j,k))*u(i,j,k)+
c     &(ccf(i,j,k+1)-ccf(i,j,k-1))*w(i,j,k)
c	normv=dsqrt(u(i,j,k)**2+w(i,j,k)**2)
c	sign=-sign
c	if (sign.gt.0) sign=1
c	if (sign.lt.0) sign=-1
c	if (normv.gt.0.1) then
c		if (sign.eq.1)	theta=115
c		if	(sign.eq.-1) theta=45
c	else
c	theta=350.0*normv*sign+80.0
	
	endif
c	theta=135.0
c	******************************************************
	mx=(ccf(i+1,j,k)-ccf(i-1,j,k))
	mz=(ccf(i,j,k+1)-ccf(i,j,k-1))
	my=(ccf(i,j+1,k)-ccf(i,j-1,k))
	mmx=mx/dsqrt(mz**2+mx**2+1.d-50)
	mmz=mz/dsqrt(mz**2+mx**2+1.d-50)
	mmx=mmx*dsind(theta)
	mmz=mmz*dsind(theta)

	mmy=dcosd(theta)

	s11(i,j,k)=sigh*(w11(i,j,k)+w11(i-1,j,k))*mx*dsind(theta)/2.d0
	s22(i,j,k)=sigh*(w11(i,j,k)+w11(i,j-1,k))*mmy
     & *dsqrt(mz**2+mx**2+my**2)*-1.d0
	s33(i,j,k)=sigh*(w11(i,j,k)+w11(i,j,k-1))*mz*dsind(theta)/2.d0	
	endif
c ***********************************************     	
	rho=rho2+ccc(i,j,k)*(1.d0-rho2)
	s11(i,j,k)=s11(i,j,j)*rho/(1-rho2)
	s22(i,j,k)=s22(i,j,j)*rho/(1-rho2)
	s33(i,j,k)=s33(i,j,j)*rho/(1-rho2)

	w11(i,j,k)=0.5d0*(s11(i,j,k)+s11(i-1,j,k))
	w22(i,j,k)=0.5d0*(s22(i,j,k)+s22(i,j-1,k))
	w33(i,j,k)=0.5d0*(s33(i,j,k)+s33(i,j,k-1))


		 enddo
	    enddo
	enddo
	


c	call copy (w11,s11,nx,ny,nz)
c	call copy (w22,s22,nx,ny,nz)
c	call copy (w33,s33,nx,ny,nz)	   

c      !???  calls boundary conditions for a homogenous field

       !???  call bcvechom(div,cc,nx,ny) !=bcv for homogenous condition
c	call bc_w(w11,w22,w33)
	endif
      return
      end

c     ************meier*************	
      SUBROUTINE meier(ccc,w11,w22,w33)
      include 'input.h'
      DOUBLE PRECISION sigh,d1x,d2x,d1y,d2y,d1z,d2z
      DOUBLE PRECISION nxu,nxd,nyu,nyd,nzu,nzd
      DOUBLE PRECISION ccf(nx,ny,nz),ccc(nx,ny,nz)
c      DOUBLE PRECISION S12t(nx,ny,nz),S13t(nx,ny,nz),S23t(nx,ny,nz)
      DOUBLE PRECISION w11(nx,ny,nz),w22(nx,ny,nz),w33(nx,ny,nz)
	INTEGER i, j, k,nf
	nf=1	
	call copy(ccf,ccc,nx,ny,nz)
	do i=1,nf
	call filterc27(ccf,w22,nx,ny,nz)
	call copy (ccf,w22,nx,ny,nz)
	call bc_out(ccf,nx,ny,nz)
	enddo

c	if (nf.ge.1)then	 
c	call filterc27(ccc,ccf,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	do i=1,nf-1
c	call copy (w22,ccf,nx,ny,nz)
c	call filterc27(w22,ccf,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	enddo
c	else
c	call copy(ccf,ccc,nx,ny,nz)
c	endif


	if(sigma.eq.0.d0) then
         call fillzero(S11,nx,ny,nz)
         call fillzero(S22,nx,ny,nz)
	   call fillzero(S33,nx,ny,nz)
	   call fillzero(S12,nx,ny,nz)
	   call fillzero(S13,nx,ny,nz)			
         call fillzero(S23,nx,ny,nz)
      else
      sigh =0.5*tau*sigma/h**2

      do j=2,ny-1
         do i=2,nx-1
			do k=2,nz-1



      nxu=dsqrt((ccf(i+1,j,k)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i,j+1,k)-ccf(i,j-1,k)+ccf(i+1,j+1,k)-ccf(i+1,j-1,k
     *))**2
     *+0.0625d0*(ccf(i,j,k+1)-ccf(i,j,k-1)+ccf(i+1,j,k+1)-ccf(i+1,j,k-1
     &))**2
     *+1.d-50)
      nxd=dsqrt((ccf(i,j,k)-ccf(i-1,j,k))**2
     *+0.0625d0*(ccf(i-1,j+1,k)-ccf(i-1,j-1,k)+ccf(i,j+1,k)-ccf(i,j-1,k
     &))**2
     *+0.0625d0*(ccf(i-1,j,k+1)-ccf(i-1,j,k-1)+ccf(i,j,k+1)-ccf(i,j,k-1
     &))**2
     *+1.d-50)
      nyu=dsqrt((ccf(i,j+1,k)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i+1,j+1,k)-ccf(i-1,j+1,k)+ccf(i+1,j,k)-ccf(i-1,j,k
     &))**2
     *+0.0625d0*(ccf(i,j+1,k+1)-ccf(i,j+1,k-1)+ccf(i,j,k+1)-ccf(i,j,k-1
     &))**2
     *+1.d-50)
      nyd=dsqrt((ccf(i,j,k)-ccf(i,j-1,k))**2
     *+0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j-1,k)-ccf(i-1,j-1,k
     &))**2
     *+0.0625d0*(ccf(i,j,k+1)-ccf(i,j,k-1)+ccf(i,j-1,k+1)-ccf(i,j-1,k-1
     &))**2
     *+1.d-50)
      nzu=dsqrt((ccf(i,j,k+1)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j,k+1)-ccf(i-1,j,k+1
     &))**2
     *+0.0625d0*(ccf(i,j+1,k)-ccf(i,j-1,k)+ccf(i,j+1,k+1)-ccf(i,j-1,k+1
     &))**2
     *+1.d-50)
      nzd=dsqrt((ccf(i,j,k)-ccf(i,j,k-1))**2
     *+0.0625d0*(ccf(i+1,j,k-1)-ccf(i-1,j,k-1)+ccf(i+1,j,k)-ccf(i-1,j,k
     &))**2
     *+0.0625d0*(ccf(i,j+1,k-1)-ccf(i,j-1,k-1)+ccf(i,j+1,k)-ccf(i,j-1,k
     &))**2
     *+1.d-50)

            d1x=(ccf(i+1,j,k)-ccf(i,j,k))/nxu
            d2x=(ccf(i,j,k)-ccf(i-1,j,k))/nxd
            d1y=(ccf(i,j+1,k)-ccf(i,j,k))/nyu
            d2y=(ccf(i,j,k)-ccf(i,j-1,k))/nyd
            d1z=(ccf(i,j,k+1)-ccf(i,j,k))/nzu
            d2z=(ccf(i,j,k)-ccf(i,j,k-1))/nzd
c	here we temporary use w11 instead of divergence
            w11(i,j,k) = - (d1x-d2x+d1y-d2y+d1z-d2z)

			enddo
		enddo
	enddo
	!!!call bcp(ccf,nx,ny)
	call bc_out(w11,nx,ny,nz)

	do	k=2,nz-1
        do j=2,ny-1
	     do i=2,nx-1
c		    S11t(i,j,k)=(0.5+Ox(i-1,j,k))*(sigh*div(i-1,j,k)*Ax(i-1,j,k))+
c     *                (0.5-Ox(i,j,k))*(sigh*div(i,j,k)*Ax(i,j,k))
c             S22t(i,j,k)=(0.5+Oy(i,j-1))*(sigh*div(i,j-1)*Ay(i,j-1))+
c     *                (0.5-Oy(i,j,k))*(sigh*div(i,j,k)*Ay(i,j,k))
c		    S33t(i,j,k)=(0.5+Oz(i-1,j,k))*(sigh*div(i-1,j,k)*Az(i-1,j,k))+
c     *                (0.5-Oz(i,j,k))*(sigh*div(i,j,k)*Az(i,j,k))
	
c		    S11t(i,j,k)=(0.5)*(sigh*div(i-1,j,k)*Ax(i-1,j,k))+
c     *                (0.5)*(sigh*div(i,j,k)*Ax(i,j,k))
c			S22t(i,j,k)=(0.5)*(sigh*div(i,j-1,k)*Ay(i,j-1,k))+
c     *                (0.5)*(sigh*div(i,j,k)*Ay(i,j,k))
c		    S33t(i,j,k)=(0.5)*(sigh*div(i-1,j,k)*Az(i-1,j,k))+
c     *                (0.5)*(sigh*div(i,j,k)*Az(i,j,k))
	s11(i,j,k)=sigh*(ccf(i,j,k)-ccf(i-1,j,k))*(w11(i,j,k)+w11(i-1,j,k
     *))
      s22(i,j,k)=sigh*(ccf(i,j,k)-ccf(i,j-1,k))*(w11(i,j,k)+w11(i,j-1,k
     *))
	s33(i,j,k)=sigh*(ccf(i,j,k)-ccf(i,j,k-1))*(w11(i,j,k)+w11(i,j,k-1
     *))
     	

		 enddo
	    enddo
	enddo

	
	call copy (w11,s11,nx,ny,nz)
	call copy (w22,s22,nx,ny,nz)
	call copy (w33,s33,nx,ny,nz)	   

c      !???  calls boundary conditions for a homogenous field

       !???  call bcvechom(div,cc,nx,ny) !=bcv for homogenous condition
c	call bc_w(w11,w22,w33)
	endif
      return
      end

c     ************SGIP*************	
      SUBROUTINE SGIP(ccc,w11,w22,w33)
      include 'input.h'
      DOUBLE PRECISION sigh,d1x,d2x,d1y,d2y,d1z,d2z
      DOUBLE PRECISION nxu,nxd,nyu,nyd,nzu,nzd
      DOUBLE PRECISION ccf(nx,ny,nz),ccc(nx,ny,nz)
c      DOUBLE PRECISION S12t(nx,ny,nz),S13t(nx,ny,nz),S23t(nx,ny,nz)
      DOUBLE PRECISION w11(nx,ny,nz),w22(nx,ny,nz),w33(nx,ny,nz)
      DOUBLE PRECISION Lx,Ly,Lz,Lxd,Lyd,Lzd,theta,sign,normv
	DOUBLE PRECISION mx,my,mz,mmx,mmy,mmz,tempa
	INTEGER i, j, k,nf,ii,ifront,iback
	nf=1
	
	call copy(ccf,ccc,nx,ny,nz)
	do i=1,nf
	call filterc27(ccf,w22,nx,ny,nz)
	call copy (ccf,w22,nx,ny,nz)
	call bc_out(ccf,nx,ny,nz)
	enddo	
c	if (nf.ge.1)then	 
c	call filterc27(ccc,ccf,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	do i=1,nf-1
c	call copy (w22,ccf,nx,ny,nz)
c	call filterc27(w22,ccf,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	enddo
c	else
c	call copy(ccf,ccc,nx,ny,nz)
c	endif


	if(sigma.eq.0.d0) then
         call fillzero(S11,nx,ny,nz)
         call fillzero(S22,nx,ny,nz)
	   call fillzero(S33,nx,ny,nz)
	   call fillzero(S12,nx,ny,nz)
	   call fillzero(S13,nx,ny,nz)			
         call fillzero(S23,nx,ny,nz)
      else
      sigh =0.5*tau**2*sigma/h**3

      do j=2,ny-1
         do i=2,nx-1
			do k=2,nz-1



      nxu=dsqrt((ccf(i+1,j,k)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i,j+1,k)-ccf(i,j-1,k)+ccf(i+1,j+1,k)-ccf(i+1,j-1,k
     *))**2
     *+0.0625d0*(ccf(i,j,k+1)-ccf(i,j,k-1)+ccf(i+1,j,k+1)-ccf(i+1,j,k-1
     &))**2
     *+1.d-50)
      nxd=dsqrt((ccf(i,j,k)-ccf(i-1,j,k))**2
     *+0.0625d0*(ccf(i-1,j+1,k)-ccf(i-1,j-1,k)+ccf(i,j+1,k)-ccf(i,j-1,k
     &))**2
     *+0.0625d0*(ccf(i-1,j,k+1)-ccf(i-1,j,k-1)+ccf(i,j,k+1)-ccf(i,j,k-1
     &))**2
     *+1.d-50)
      nyu=dsqrt((ccf(i,j+1,k)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i+1,j+1,k)-ccf(i-1,j+1,k)+ccf(i+1,j,k)-ccf(i-1,j,k
     &))**2
     *+0.0625d0*(ccf(i,j+1,k+1)-ccf(i,j+1,k-1)+ccf(i,j,k+1)-ccf(i,j,k-1
     &))**2
     *+1.d-50)
      nyd=dsqrt((ccf(i,j,k)-ccf(i,j-1,k))**2
     *+0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j-1,k)-ccf(i-1,j-1,k
     &))**2
     *+0.0625d0*(ccf(i,j,k+1)-ccf(i,j,k-1)+ccf(i,j-1,k+1)-ccf(i,j-1,k-1
     &))**2
     *+1.d-50)
      nzu=dsqrt((ccf(i,j,k+1)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j,k+1)-ccf(i-1,j,k+1
     &))**2
     *+0.0625d0*(ccf(i,j+1,k)-ccf(i,j-1,k)+ccf(i,j+1,k+1)-ccf(i,j-1,k+1
     &))**2
     *+1.d-50)
      nzd=dsqrt((ccf(i,j,k)-ccf(i,j,k-1))**2
     *+0.0625d0*(ccf(i+1,j,k-1)-ccf(i-1,j,k-1)+ccf(i+1,j,k)-ccf(i-1,j,k
     &))**2
     *+0.0625d0*(ccf(i,j+1,k-1)-ccf(i,j-1,k-1)+ccf(i,j+1,k)-ccf(i,j-1,k
     &))**2
     *+1.d-50)

            d1x=(ccf(i+1,j,k)-ccf(i,j,k))/nxu
            d2x=(ccf(i,j,k)-ccf(i-1,j,k))/nxd
            d1y=(ccf(i,j+1,k)-ccf(i,j,k))/nyu
            d2y=(ccf(i,j,k)-ccf(i,j-1,k))/nyd
            d1z=(ccf(i,j,k+1)-ccf(i,j,k))/nzu
            d2z=(ccf(i,j,k)-ccf(i,j,k-1))/nzd
c	here we temporary use w11 instead of divergence
            w11(i,j,k) = - (d1x-d2x+d1y-d2y+d1z-d2z)

			enddo
		enddo
	enddo
	call bc_out(w11,nx,ny,nz)

	do	k=2,nz-1
        do j=2,ny-1
	     do i=2,nx-1

	Lx=(hxmax(i,j,k)+hxmin(i,j,k))/2.d0
	Lxd=(hxmax(i-1,j,k)+hxmin(i-1,j,k))/2.d0
	Ly=(hymax(i,j,k)+hymin(i,j,k))/2.d0
	Lyd=(hymax(i,j-1,k)+hymin(i,j-1,k))/2.d0
	Lz=(hzmax(i,j,k)+hzmin(i,j,k))/2.d0
	Lzd=(hzmax(i,j,k-1)+hzmin(i,j,k-1))/2.d0

	
	s11(i,j,k)=sigh*(w11(i,j,k)+w11(i-1,j,k))*	
     &(ccf(i+1,j,k)-ccf(i-1,j,k))*dabs(Lx-Lxd)/2.d0
      s22(i,j,k)=sigh*(w11(i,j,k)+w11(i,j-1,k))*
     &(ccf(i,j+1,k)-ccf(i,j-1,k))*dabs(Ly-Lyd)/2.d0
	s33(i,j,k)=sigh*(w11(i,j,k)+w11(i,j,k-1))*
     &(ccf(i,j,k+1)-ccf(i,j,k-1))*dabs(Lz-Lzd)/2.d0
c	if (s11(i,j,k).ne.0.d0)  write(25,*) s11(i,j,k),i,j,k


     	if (j.eq.2) then
c	********************************************************
	mx=(ccf(i+1,j,k)-ccf(i-1,j,k))
	mz=(ccf(i,j,k+1)-ccf(i,j,k-1))
	mmx=mx/dsqrt(mz**2+mx**2+1.d-50)
	mmz=mz/dsqrt(mz**2+mx**2+1.d-50)
	normv=mmx*u(i,j,k)+mmz*w(i,j,k)
	normv=-normv
	if (dabs(normv).gt.0.1) then
		if (normv.gt.0)	theta=115
		if	(normv.lt.0) theta=45
	else
	theta=350.0*normv+80


c	**** determination of theta as a function of Vcl *****
c	sign=(ccf(i+1,j,k)-ccf(i-1,j,k))*u(i,j,k)+
c     &(ccf(i,j,k+1)-ccf(i,j,k-1))*w(i,j,k)
c	normv=dsqrt(u(i,j,k)**2+w(i,j,k)**2)
c	sign=-sign
c	if (sign.gt.0) sign=1
c	if (sign.lt.0) sign=-1
c	if (normv.gt.0.1) then
c		if (sign.eq.1)	theta=115
c		if	(sign.eq.-1) theta=45
c	else
c	theta=350.0*normv*sign+80.0
	
	endif
c	theta=135.0
c	******************************************************
	mx=(ccf(i+1,j,k)-ccf(i-1,j,k))
	mz=(ccf(i,j,k+1)-ccf(i,j,k-1))
	my=(ccf(i,j+1,k)-ccf(i,j-1,k))
	mmx=mx/dsqrt(mz**2+mx**2+1.d-50)
	mmz=mz/dsqrt(mz**2+mx**2+1.d-50)
	mmx=mmx*dsind(theta)
	mmz=mmz*dsind(theta)

	mmy=dcosd(theta)

	s11(i,j,k)=sigh*(w11(i,j,k)+w11(i-1,j,k))*mx*dsind(theta)/2.d0
	s22(i,j,k)=sigh*(w11(i,j,k)+w11(i,j-1,k))*mmy
     & *dsqrt(mz**2+mx**2+my**2)*-1.d0
	s33(i,j,k)=sigh*(w11(i,j,k)+w11(i,j,k-1))*mz*dsind(theta)/2.d0	
	endif
		 enddo
	    enddo
	enddo
		
	call bc_out(s11,nx,ny,nz)
	call bc_out(s22,nx,ny,nz)
	call bc_out(s33,nx,ny,nz)
	tempa=-10000
	j=2
	do i=2,nx-1
	do k=2,ny-1
	if (dabs(s11(i,j,k)).gt.tempa)	tempa=s11(i,j,k)
	enddo
	enddo
c	write(*,*) tempa
	
	call copy (w11,s11,nx,ny,nz)
	call copy (w22,s22,nx,ny,nz)
	call copy (w33,s33,nx,ny,nz)	   

c      !???  calls boundary conditions for a homogenous field

       !???  call bcvechom(div,cc,nx,ny) !=bcv for homogenous condition
c	call bc_w(w11,w22,w33)
	endif
      return
      end

c     ************SGIP2 secound order*************	
      SUBROUTINE SGIP2(ccc,w11,w22,w33)
      include 'input.h'
      DOUBLE PRECISION sigh,d1x,d2x,d1y,d2y,d1z,d2z
      DOUBLE PRECISION nxu,nxd,nyu,nyd,nzu,nzd
      DOUBLE PRECISION ccf(nx,ny,nz),ccc(nx,ny,nz)
c      DOUBLE PRECISION S12t(nx,ny,nz),S13t(nx,ny,nz),S23t(nx,ny,nz)
      DOUBLE PRECISION w11(nx,ny,nz),w22(nx,ny,nz),w33(nx,ny,nz)
      DOUBLE PRECISION Lxu,Lyu,Lzu,Lxd,Lyd,Lzd


	INTEGER i, j, k,nf
	nf=3	

	call copy(ccf,ccc,nx,ny,nz)
	do i=1,nf
	call filterc27(ccf,w22,nx,ny,nz)
	call copy (ccf,w22,nx,ny,nz)
	call bc_out(ccf,nx,ny,nz)
	enddo

c	if (nf.ge.1)then	 
c	call filterc27(ccc,ccf,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	do i=1,nf-1
c	call copy (w22,ccf,nx,ny,nz)
c	call filterc27(w22,ccf,nx,ny,nz)
c	call bc_out(ccf,nx,ny,nz)
c	enddo
c	else
c	call copy(ccf,ccc,nx,ny,nz)
c	endif


	if(sigma.eq.0.d0) then
         call fillzero(S11,nx,ny,nz)
         call fillzero(S22,nx,ny,nz)
	   call fillzero(S33,nx,ny,nz)
	   call fillzero(S12,nx,ny,nz)
	   call fillzero(S13,nx,ny,nz)			
         call fillzero(S23,nx,ny,nz)
      else
      sigh =0.5*tau*sigma/h**2

      do j=2,ny-1
         do i=2,nx-1
			do k=2,nz-1



      nxu=dsqrt((ccf(i+1,j,k)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i,j+1,k)-ccf(i,j-1,k)+ccf(i+1,j+1,k)-ccf(i+1,j-1,k
     *))**2
     *+0.0625d0*(ccf(i,j,k+1)-ccf(i,j,k-1)+ccf(i+1,j,k+1)-ccf(i+1,j,k-1
     &))**2
     *+1.d-50)
      nxd=dsqrt((ccf(i,j,k)-ccf(i-1,j,k))**2
     *+0.0625d0*(ccf(i-1,j+1,k)-ccf(i-1,j-1,k)+ccf(i,j+1,k)-ccf(i,j-1,k
     &))**2
     *+0.0625d0*(ccf(i-1,j,k+1)-ccf(i-1,j,k-1)+ccf(i,j,k+1)-ccf(i,j,k-1
     &))**2
     *+1.d-50)
      nyu=dsqrt((ccf(i,j+1,k)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i+1,j+1,k)-ccf(i-1,j+1,k)+ccf(i+1,j,k)-ccf(i-1,j,k
     &))**2
     *+0.0625d0*(ccf(i,j+1,k+1)-ccf(i,j+1,k-1)+ccf(i,j,k+1)-ccf(i,j,k-1
     &))**2
     *+1.d-50)
      nyd=dsqrt((ccf(i,j,k)-ccf(i,j-1,k))**2
     *+0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j-1,k)-ccf(i-1,j-1,k
     &))**2
     *+0.0625d0*(ccf(i,j,k+1)-ccf(i,j,k-1)+ccf(i,j-1,k+1)-ccf(i,j-1,k-1
     &))**2
     *+1.d-50)
      nzu=dsqrt((ccf(i,j,k+1)-ccf(i,j,k))**2
     *+0.0625d0*(ccf(i+1,j,k)-ccf(i-1,j,k)+ccf(i+1,j,k+1)-ccf(i-1,j,k+1
     &))**2
     *+0.0625d0*(ccf(i,j+1,k)-ccf(i,j-1,k)+ccf(i,j+1,k+1)-ccf(i,j-1,k+1
     &))**2
     *+1.d-50)
      nzd=dsqrt((ccf(i,j,k)-ccf(i,j,k-1))**2
     *+0.0625d0*(ccf(i+1,j,k-1)-ccf(i-1,j,k-1)+ccf(i+1,j,k)-ccf(i-1,j,k
     &))**2
     *+0.0625d0*(ccf(i,j+1,k-1)-ccf(i,j-1,k-1)+ccf(i,j+1,k)-ccf(i,j-1,k
     &))**2
     *+1.d-50)

            d1x=(ccf(i+1,j,k)-ccf(i,j,k))/nxu
            d2x=(ccf(i,j,k)-ccf(i-1,j,k))/nxd
            d1y=(ccf(i,j+1,k)-ccf(i,j,k))/nyu
            d2y=(ccf(i,j,k)-ccf(i,j-1,k))/nyd
            d1z=(ccf(i,j,k+1)-ccf(i,j,k))/nzu
            d2z=(ccf(i,j,k)-ccf(i,j,k-1))/nzd
c	here we temporary use w11 instead of divergence
            w11(i,j,k) = - (d1x-d2x+d1y-d2y+d1z-d2z)

			enddo
		enddo
	enddo
	!!!call bcp(ccf,nx,ny)
	call bc_out(w11,nx,ny,nz)

	do	k=2,nz-1
        do j=2,ny-1
	     do i=2,nx-1

	Lxu=(hxmax(i+1,j,k)+hxmin(i+1,j,k))/2
	Lxd=(hxmax(i-1,j,k)+hxmin(i-1,j,k))/2
	Lyu=(hymax(i,j+1,k)+hymin(i,j+1,k))/2
	Lyd=(hymax(i,j-1,k)+hymin(i,j-1,k))/2
	Lzu=(hzmax(i,j,k+1)+hzmin(i,j,k+1))/2
	Lzd=(hzmax(i,j,k-1)+hzmin(i,j,k-1))/2

	
	s11(i,j,k)=sigh*(w11(i+1,j,k)+w11(i-1,j,k))*	
     &(ccf(i+1,j,k)-ccf(i-1,j,k))*(Lxu-Lxd)*0.25
      s22(i,j,k)=sigh*(w11(i,j+1,k)+w11(i,j-1,k))*
     &(ccf(i,j+1,k)-ccf(i,j-1,k))*(Lyu-Lyd)*0.25
	s33(i,j,k)=sigh*(w11(i,j,k+1)+w11(i,j,k-1))*
     &(ccf(i,j,k+1)-ccf(i,j,k-1))*(Lzu-Lzd)*0.25
c	if (s11(i,j,k).ne.0.d0)  write(25,*) s11(i,j,k),i,j,k
     	

		 enddo
	    enddo
	enddo

	
	call copy (w11,s11,nx,ny,nz)
	call copy (w22,s22,nx,ny,nz)
	call copy (w33,s33,nx,ny,nz)	   

c      !???  calls boundary conditions for a homogenous field

       !???  call bcvechom(div,cc,nx,ny) !=bcv for homogenous condition
c	call bc_w(w11,w22,w33)
	endif
      return
      end

c***********************************************************************
c                      SKETCH OF THE C-GRID in 3-D
c     
c     S11,S22,S33,divergence,u**2,v**2,cc are defined in the same 
c     position as p, while (i,j,k) are the default indices 
c
c                                       u                  p
c             p(i-1) *------------------*------------------*
c                   /.                                    /|
c                  / .                S12                / |
c          v(i-1) *  .               *                v *  |
c                /   .                                 /   |
c               /    .       u(j-1)            p(j-1) /    |
c   p(i-1,j-1) *------------------*------------------*	   |
c              |     .                 S13           |     |
c              |     + w(i-1)           +            |     * w
c              |     .                               |     |
c              |     .                               | S23 |
c     S23(i-1) |  +  .                               |  *  |
c              |     .                               |     |
c              |     .         S13(j-1)              |     |
c   w(i-1,j-1) *     .            *           w(j-1) *     |
c              |     .                               |     |
c              |     + . . . . . . . . .+ . . . . . .| . . * p(k-1) 
c              |    . p(i-1,k-1)      u(k-1)         |    /         
c              |   .                                 |   /         
c              |  + v(i-1,k-1)       +               |  * v(k-1)   
c              | .                  S12(k-1)         | /
c              |.                                    |/
c              *------------------*------------------* p(j-1,k-1)
c     p(i-1,j-1,k-1)          u(j-1,k-1)
c
c
c***
c
c
c        SKETCH OF THE TOP-BOTTOM BOUNDARY CONDITIONS
c
c     
c      1       2       3         N-2     N-1      N
c      *-------*-------*--      --*-------*-------*  k=N;  p(i,j,N)
c      |       |       |          |       |       |
c  ....^.......^.......^..      ..^.......^.......^        w(i,j,N)=0
c      |       |       |          | u(N-1,j,N-1)  |
c      *-------*-------*--      --*--->---*-------*  k=N-1
c      |       |       |          |       |       |
c      |       |       |          |       |       |      nx=ny=nz=N
c      |       |       |          |       |       |
c      *-------*-------*--      --*-------*-------*
c      |       |       |          |       |       |     
c
c
c      | S12(2,j,3)    |          |       |       |
c      *---x---*-------*--      --*-------*-------*  k=3
c      |       |       |          |       |       |     
c      |       |       |          |       |       |     
c      |       |       |          |       |       |     
c      *-------*-------*--      --*-------*-------* p(N,j,2)
c      |   S23(2,j,2)  |          |       |       |     
c  ....^.......x.......^..      ..^.......^.......^  z=0;  w(i,j,2)=0
c      |       |       |          |       |       |     
c      *-------*-------*--      --*-------*-------*  k=1;  p(i,j,1)
c      |       |       |          |       |       |     
c      ^   x   ^       ^          ^       ^       ^
c      S13(2,j,1)                     w(N-1,j,1)
c
c
c***********************************************************************
      SUBROUTINE solvelpres(cc,w1,w2,w3,a,b,c,ei,t)
c***
      INCLUDE 'input.h'
      INTEGER i,j,k,t,imax,jmax,kmax
      DOUBLE PRECISION w1(nx,ny,nz), w2(nx,ny,nz), w3(nx,ny,nz)
      DOUBLE PRECISION a(nx,ny,nz),  b(nx,ny,nz),  c(nx,ny,nz)
      DOUBLE PRECISION cc(nx,ny,nz), ei(nx,ny,nz)
      DOUBLE PRECISION r1,r2,r10,r20,hi2,stdgrav,temp
      DOUBLE PRECISION horgrav, horgr1, volliq
      DOUBLE PRECISION epsdiv,MINDIV, MAXDIV
      PARAMETER (MINDIV=2.5d0, MAXDIV=25.d0)
      DOUBLE PRECISION SUMFIELD,FMODMAX,turbstx,turbsty,turbstz,ftr
	DOUBLE PRECISION ftrctr,coef
      EXTERNAL SUMFIELD,FMODMAX
      INTRINSIC MOD
c     SUBROUTINES: abcei, bc_tensor, bc_v, bc_w, mglin
c***
c     define constants, renormalize gravity
c***
      ftr=0. ! (default is 1)if it is zero, no turbulence in momentum eq.
      ftrctr=0.25d0*tauh2*ftr
      ftr=2.d0*ftr*tauh2

      r20 = tauh2*mu2
      r10 = 0.25d0*tauh2*(mu1-mu2)
      r1 = 2.d0*tauh2*(mu1-mu2)
      r2 = 2.d0*tauh2*mu2
      stdgrav= g * tau**2/h
      horgrav = g1 * tau**2/h
c***
c     total liquid volume fraction and lift force to compensate
c     gravity and keep droplets/bubble steady
c***
      volliq = SUMFIELD(cc,nx,ny,nz)/((nx-2)*(ny-2)*(nz-2))
      horgr1 = horgrav*(rho1*volliq + rho2*(1.d0 - volliq))

	coef=1.
	if (iturb.eq.4) coef=.0


c***
c     add the viscous part to the stress tensor
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               temp = r2 + r1*cc(i,j,k)+ftr*coef*mut(i,j,k)
               S11(i,j,k) = - S11(i,j,k) +
     %              temp * (u(i+1,j,k) - u(i,j,k))
               S22(i,j,k) = - S22(i,j,k) +
     %              temp * (v(i,j+1,k) - v(i,j,k))
               S33(i,j,k) = - S33(i,j,k) +
     %              temp * (w(i,j,k+1) - w(i,j,k))
	    temp=r20+r10*(cc(i,j,k)+cc(i-1,j,k)+cc(i,j-1,k)+cc(i-1,j-1,k))
     &+ftrctr*coef*(mut(i,j,k)+mut(i-1,j,k)+mut(i,j-1,k)+mut(i-1,j-1,k))
               S12(i,j,k) = - S12(i,j,k) + (temp) *
     %              (u(i,j,k) - u(i,j-1,k) + v(i,j,k) - v(i-1,j,k))
		temp=r20+r10*(cc(i,j,k)+cc(i-1,j,k)+cc(i,j,k-1)+cc(i-1,j,k-1))
     &+ftrctr*coef*(mut(i,j,k)+mut(i-1,j,k)+mut(i,j,k-1)+mut(i-1,j,k-1))
			 S13(i,j,k) = - S13(i,j,k) + (temp) *
     %              (u(i,j,k) - u(i,j,k-1) + w(i,j,k) - w(i-1,j,k))
		temp=r20+r10*(cc(i,j,k)+cc(i,j-1,k)+cc(i,j,k-1)+cc(i,j-1,k-1))
     &+ftrctr*coef*(mut(i,j,k)+mut(i,j-1,k)+mut(i,j,k-1)+mut(i,j-1,k-1))
	         S23(i,j,k) = - S23(i,j,k) + (temp) *
     %              (v(i,j,k) - v(i,j,k-1) + w(i,j,k) - w(i,j-1,k))
            enddo 
         enddo 
      enddo              
c***
      call bc_tensor()
c***
      r2 = rho2
      r1 = .5d0*(rho1 - rho2)
c***
      call abcei(a,b,c,ei,cc,r1,r2,ibc,nx,ny,nz)
c***
c     advance velocity field in time: viscous and gravity terms
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
	      turbstx=(turbst(i,j,k)+turbst(i-1,j,k))/2.0d0
		  turbsty=(turbst(i,j,k)+turbst(i,j-1,k))/2.0d0
		  turbstz=(turbst(i,j,k)+turbst(i,j,k-1))/2.0d0
               w1(i,j,k) = u(i,j,k) + a(i-1,j,k) * (S11(i,j,k) -
     %              S11(i-1,j,k) + S12(i,j+1,k) - S12(i,j,k) +
     %              S13(i,j,k+1) - S13(i,j,k)-turbstx) + horgrav
c     %              - horgr1/(r1*(cc(i,j,k)+cc(i,j-1,k)) + r2)
               w2(i,j,k) = v(i,j,k) + b(i,j-1,k) * (S12(i+1,j,k) -
     %              S12(i,j,k) + S22(i,j,k) - S22(i,j-1,k) +
     %              S23(i,j,k+1) - S23(i,j,k)-turbsty)+ stdgrav
               w3(i,j,k) = w(i,j,k) + c(i,j,k-1) * (S13(i+1,j,k) -
     %              S13(i,j,k) + S23(i,j+1,k) - S23(i,j,k) +
     %              S33(i,j,k) - S33(i,j,k-1)-turbstz) 
            enddo
         enddo
      enddo
c***
      call bc_w(w1,w2,w3)

c *********** compute the convection terms in momentum equation

      call convective (w1,w2,w3)

c***
      call bc_v()
c***
c     now let w3 = (div v) / h
c***
      hi2=1.d0/(h*h)
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               w3(i,j,k) = (u(i+1,j,k) - u(i,j,k) + v(i,j+1,k) -
     %                      v(i,j,k) + w(i,j,k+1) - w(i,j,k)) * hi2
            enddo
         enddo
      enddo
c***
c     check calls
c***
      if (t .le. 1) then
         write (6,*) "divergence sum = ", SUMFIELD(w3,nx,ny,nz)
         write (23,*) "divergence sum = ", SUMFIELD(w3,nx,ny,nz)
      endif
      if (MOD(t,tprint) .eq. 0) then
         write (*,*) "max div. before = ",
     %        h*FMODMAX(w3,nx,ny,nz,imax,jmax,kmax)
         write (23,*) "max div. before = ",
     %        h*FMODMAX(w3,nx,ny,nz,imax,jmax,kmax),
     %        imax,jmax,kmax
      endif 
c***
c     call the pressure solver
c***
      call mglin(r1,r2)
c***  
c     make the velocity field about divergence free
c***  
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               u(i,j,k) = u(i,j,k) - (w1(i,j,k)-w1(i-1,j,k))*a(i-1,j,k)
               v(i,j,k) = v(i,j,k) - (w1(i,j,k)-w1(i,j-1,k))*b(i,j-1,k)
               w(i,j,k) = w(i,j,k) - (w1(i,j,k)-w1(i,j,k-1))*c(i,j,k-1)
            enddo
         enddo
      enddo
c***
      call bc_v()
c***
c     check if more or fewer relaxations are needed
c***
      if (t .le. 1 .or. MOD(t,10) .eq. 0) then
         do k=2,nz-1
            do j=2,ny-1
               do i=2,nx-1
                  w3(i,j,k) = (u(i+1,j,k) - u(i,j,k) + v(i,j+1,k) -
     %                         v(i,j,k) + w(i,j,k+1) - w(i,j,k))*hi2
               enddo
            enddo
         enddo
         epsdiv = h*FMODMAX(w3,nx,ny,nz,imax,jmax,kmax)
         if (epsdiv .gt. MAXDIV.and. ntop.lt.30) then
            NTOP = NTOP + 2
            NPRE = NPRE + 2
            NPOST = NPOST + 2
            write (*,*) 'more relaxations: t,NPRE,NTOP,NPOST:',
     %           t,NPRE,NTOP,NPOST
            write (23,*) 'more relaxations: t,NPRE,NTOP,NPOST:',
     %           t,NPRE,NTOP,NPOST
         elseif (epsdiv .lt. MINDIV.and. ntop.gt.4) then
            NPRE = NPRE - 1
            NTOP = NTOP - 1
            NPOST = NPOST - 1
            write (*,*) 'fewer relaxations: t,NPRE,NTOP,NPOST:',
     %           t,NPRE,NTOP,NPOST
            write (23,*) 'fewer relaxations: t,NPRE,NTOP,NPOST:',
     %           t,NPRE,NTOP,NPOST
         endif
      endif 
c***
c     check calls
c***
      if (MOD(t,tprint) .eq. 0) then
         do k=2,nz-1
            do j=2,ny-1
               do i=2,nx-1
                  w3(i,j,k) = (u(i+1,j,k) - u(i,j,k) + v(i,j+1,k) -
     %                         v(i,j,k) + w(i,j,k+1) - w(i,j,k))*hi2
               enddo
            enddo
         enddo
         write (6,*) "divergence sum = ", SUMFIELD(w3,nx,ny,nz)
         write (23,*) "divergence sum = ", SUMFIELD(w3,nx,ny,nz)
         write (*,*) "max div. after = ",
     %                h*FMODMAX(w3,nx,ny,nz,imax,jmax,kmax)
         write (23,*) "max div. after = ",
     %                h*FMODMAX(w3,nx,ny,nz,imax,jmax,kmax),
     %                imax,jmax,kmax
      endif
c***
      return
      end
c *********** compute the convection terms in momentum equation
      subroutine convective (w1,w2,w3)
      include 'input.h'
      integer i,j,k
      DOUBLE PRECISION w1(nx,ny,nz),w2(nx,ny,nz),w3(nx,ny,nz)
      DOUBLE PRECISION s1,s2,t1,t2,q1,q2

      if(convectype.eq.1) then
c
c  UPWIND SCHEME
c
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1

               t1=.5d0*(w2(i,j,k) + w2(i-1,j,k))
               t2=.5d0*(w2(i,j+1,k) + w2(i-1,j+1,k))
               s1=.5d0*(w1(i-1,j,k) + w1(i,j,k))
               s2=.5d0*(w1(i,j,k) + w1(i+1,j,k))
               q1=.5d0*(w3(i,j,k) + w3(i-1,j,k))
               q2=.5d0*(w3(i,j,k+1) + w3(i-1,j,k+1))

                u(i,j,k) = w1(i,j,k) + tauh* (
     %         + dmax1(t1,0.d0)*w1(i,j-1,k) + dmin1(t1,0.d0)*w1(i,j,k)
     %         - dmax1(t2,0.d0)*w1(i,j,k) - dmin1(t2,0.d0)*w1(i,j+1,k)
     %         + dmax1(s1,0.d0)*w1(i-1,j,k) + dmin1(s1,0.d0)*w1(i,j,k)
     %         - dmax1(s2,0.d0)*w1(i,j,k) - dmin1(s2,0.d0)*w1(i+1,j,k)
     %         + dmax1(q1,0.d0)*w1(i,j,k-1) + dmin1(q1,0.d0)*w1(i,j,k)
     %         - dmax1(q2,0.d0)*w1(i,j,k) - dmin1(q2,0.d0)*w1(i,j,k+1)
     %         +u2p(i-1,j,k)-u2p(i,j,k)+0.25d0*( 
     %         +uvp(i,j,k)-uvp(i,j+1,k)+uvp(i,j-1,k)-uvp(i,j,k)
     %        +uvp(i-1,j,k)-uvp(i-1,j+1,k)+uvp(i-1,j-1,k)-uvp(i-1,j,k) 
     %         +uwp(i,j,k)-uwp(i,j,k+1)+uwp(i,j,k-1)-uwp(i,j,k)
     %        +uwp(i-1,j,k)-uwp(i-1,j,k+1)+uwp(i-1,j,k-1)-uwp(i-1,j,k)))



               t1=.5d0*(w1(i,j,k) + w1(i,j-1,k))
               t2=.5d0*(w1(i+1,j,k) + w1(i+1,j-1,k))
               s1=.5d0*(w2(i,j-1,k) + w2(i,j,k))
               s2=.5d0*(w2(i,j,k) + w2(i,j+1,k))
               q1=.5d0*(w3(i,j,k) + w3(i,j-1,k))
               q2=.5d0*(w3(i,j,k+1) + w3(i,j-1,k+1))

               v(i,j,k) = w2(i,j,k) + tauh* (
     %         + dmax1(t1,0.d0)*w2(i-1,j,k) + dmin1(t1,0.d0)*w2(i,j,k)
     %         - dmax1(t2,0.d0)*w2(i,j,k) - dmin1(t2,0.d0)*w2(i+1,j,k)
     %         + dmax1(s1,0.d0)*w2(i,j-1,k) + dmin1(s1,0.d0)*w2(i,j,k)
     %         - dmax1(s2,0.d0)*w2(i,j,k) - dmin1(s2,0.d0)*w2(i,j+1,k)
     %         + dmax1(q1,0.d0)*w2(i,j,k-1) + dmin1(q1,0.d0)*w2(i,j,k)
     %         - dmax1(q2,0.d0)*w2(i,j,k) - dmin1(q2,0.d0)*w2(i,j,k+1)
     %         +v2p(i,j-1,k)-v2p(i,j,k)+0.25d0*(
     %        +uvp(i-1,j,k)-uvp(i,j,k)+uvp(i,j,k)-uvp(i+1,j,k)
     %        +uvp(i-1,j-1,k)-uvp(i,j-1,k)+uvp(i,j-1,k)-uvp(i+1,j-1,k)
     %        +vwp(i,j,k-1)-vwp(i,j,k)+vwp(i,j,k)-vwp(i,j,k+1)
     %        +vwp(i,j-1,k-1)-vwp(i,j-1,k)+vwp(i,j-1,k)-vwp(i,j-1,k+1)))


               t1=.5d0*(w1(i,j,k) + w1(i,j,k-1))
               t2=.5d0*(w1(i+1,j,k) + w1(i+1,j,k-1))
               s1=.5d0*(w3(i,j,k-1) + w3(i,j,k))
               s2=.5d0*(w3(i,j,k) + w3(i,j,k+1))
               q1=.5d0*(w2(i,j,k) + w2(i,j,k-1))
               q2=.5d0*(w2(i,j+1,k) + w2(i,j+1,k-1))

               w(i,j,k) = w3(i,j,k) + tauh* (
     %         + dmax1(t1,0.d0)*w3(i-1,j,k) + dmin1(t1,0.d0)*w3(i,j,k)
     %         - dmax1(t2,0.d0)*w3(i,j,k) - dmin1(t2,0.d0)*w3(i+1,j,k)
     %         + dmax1(s1,0.d0)*w3(i,j,k-1) + dmin1(s1,0.d0)*w3(i,j,k)
     %         - dmax1(s2,0.d0)*w3(i,j,k) - dmin1(s2,0.d0)*w3(i,j,k+1)
     %         + dmax1(q1,0.d0)*w3(i,j-1,k) + dmin1(q1,0.d0)*w3(i,j,k)
     %         - dmax1(q2,0.d0)*w3(i,j,k) - dmin1(q2,0.d0)*w3(i,j+1,k)
     %         +w2p(i,j,k-1)-w2p(i,j,k)+0.25d0*( 
     %        +uwp(i-1,j,k)-uwp(i,j,k)+uwp(i,j,k)-uwp(i+1,j,k)
     %        +uwp(i-1,j,k-1)-uwp(i,j,k-1)+uwp(i,j,k-1)-uwp(i+1,j,k-1)     
     %        +vwp(i,j-1,k)-vwp(i,j,k)+vwp(i,j,k)-vwp(i,j+1,k)
     %        +vwp(i,j-1,k-1)-vwp(i,j,k-1)+vwp(i,j,k-1)-vwp(i,j+1,k-1)))    

            enddo
         enddo
      enddo
      else if(convectype.eq.2) then
c
c  Centeral difference
c
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1

               u(i,j,k) = w1(i,j,k) + tauh*(0.25d0 * (
     %         (w1(i,j,k)+w1(i-1,j,k))**2 - (w1(i+1,j,k)+w1(i,j,k))**2
     %         + (w1(i,j-1,k)+w1(i,j,k))*(w2(i,j,k)  +w2(i-1,j,k))
     %         - (w1(i,j+1,k)+w1(i,j,k))*(w2(i,j+1,k)+w2(i-1,j+1,k))
     %         + (w1(i,j,k-1)+w1(i,j,k))*(w3(i,j,k)  +w3(i-1,j,k))
     %         - (w1(i,j,k+1)+w1(i,j,k))*(w3(i,j,k+1)+w3(i-1,j,k+1)) 
     %         +uvp(i,j,k)-uvp(i,j+1,k)+uvp(i,j-1,k)-uvp(i,j,k)
     %        +uvp(i-1,j,k)-uvp(i-1,j+1,k)+uvp(i-1,j-1,k)-uvp(i-1,j,k) 
     %         +uwp(i,j,k)-uwp(i,j,k+1)+uwp(i,j,k-1)-uwp(i,j,k)
     %        +uwp(i-1,j,k)-uwp(i-1,j,k+1)+uwp(i-1,j,k-1)-uwp(i-1,j,k) )
     %         +u2p(i-1,j,k)-u2p(i,j,k))                                 

               v(i,j,k) = w2(i,j,k) + tauh*(0.25d0 * (
     %         (w2(i,j,k)+w2(i,j-1,k))**2 - (w2(i,j+1,k)+w2(i,j,k))**2
     %         + (w2(i-1,j,k)+w2(i,j,k))*(w1(i,j,k)  +w1(i,j-1,k))
     %         - (w2(i+1,j,k)+w2(i,j,k))*(w1(i+1,j,k)+w1(i+1,j-1,k))
     %         + (w2(i,j,k-1)+w2(i,j,k))*(w3(i,j,k)  +w3(i,j-1,k))
     %         - (w2(i,j,k+1)+w2(i,j,k))*(w3(i,j,k+1)+w3(i,j-1,k+1)) 
     %        +uvp(i-1,j,k)-uvp(i,j,k)+uvp(i,j,k)-uvp(i+1,j,k)
     %        +uvp(i-1,j-1,k)-uvp(i,j-1,k)+uvp(i,j-1,k)-uvp(i+1,j-1,k)
     %        +vwp(i,j,k-1)-vwp(i,j,k)+vwp(i,j,k)-vwp(i,j,k+1)
     %        +vwp(i,j-1,k-1)-vwp(i,j-1,k)+vwp(i,j-1,k)-vwp(i,j-1,k+1) )
     %        +v2p(i,j-1,k)-v2p(i,j,k))                         
                
               w(i,j,k) = w3(i,j,k) + tauh*(0.25d0 * (
     %         (w3(i,j,k)+w3(i,j,k-1))**2 - (w3(i,j,k+1)+w3(i,j,k))**2
     %         + (w3(i-1,j,k)+w3(i,j,k))*(w1(i,j,k)  +w1(i,j,k-1))
     %         - (w3(i+1,j,k)+w3(i,j,k))*(w1(i+1,j,k)+w1(i+1,j,k-1))
     %         + (w3(i,j-1,k)+w3(i,j,k))*(w2(i,j,k)  +w2(i,j,k-1))
     %         - (w3(i,j+1,k)+w3(i,j,k))*(w2(i,j+1,k)+w2(i,j+1,k-1)) 
     %        +uwp(i-1,j,k)-uwp(i,j,k)+uwp(i,j,k)-uwp(i+1,j,k)
     %        +uwp(i-1,j,k-1)-uwp(i,j,k-1)+uwp(i,j,k-1)-uwp(i+1,j,k-1)     
     %        +vwp(i,j-1,k)-vwp(i,j,k)+vwp(i,j,k)-vwp(i,j+1,k)
     %        +vwp(i,j-1,k-1)-vwp(i,j,k-1)+vwp(i,j,k-1)-vwp(i,j+1,k-1) )    
     %         +w2p(i,j,k-1)-w2p(i,j,k))                    
            enddo
         enddo
      enddo
      else
c zero convection term case
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               u(i,j,k) = w1(i,j,k)
               v(i,j,k) = w2(i,j,k)
               w(i,j,k) = w3(i,j,k)
            enddo
         enddo
      enddo
      endif
      return
      end
c***********************************************************************
      SUBROUTINE abcei (a,b,c,ei,cc,r1,r2,ibc,nx,ny,nz)
c***
      IMPLICIT none
      INTEGER nx,ny,nz,i,j,k,ibc
      DOUBLE PRECISION a(nx,ny,nz),b(nx,ny,nz),c(nx,ny,nz),
     %                 ei(nx,ny,nz),cc(nx,ny,nz)
      DOUBLE PRECISION r1,r2
c     SUBROUTINES: bc_coef
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               a(i,j,k) = 1.d0 / (r1*(cc(i,j,k)+cc(i+1,j,k)) + r2)
               b(i,j,k) = 1.d0 / (r1*(cc(i,j,k)+cc(i,j+1,k)) + r2)
               c(i,j,k) = 1.d0 / (r1*(cc(i,j,k)+cc(i,j,k+1)) + r2)
            enddo
         enddo
      enddo
c***
c     need indices i=1 for a, j=1 for b and k=1 for c to calculate ei
c***
      call bc_coef(a,b,c,cc,r1,r2,ibc,nx,ny,nz)
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               ei(i,j,k) = 1.d0/(a(i,j,k)   + a(i-1,j,k) + b(i,j,k) +
     %                           b(i,j-1,k) + c(i,j,k)   + c(i,j,k-1))
            enddo
         enddo
      enddo
c***
      return
      end
c***********************************************************************
c     This subroutine invert the pressure equation:
c
c     div (1/rho grad p) = rhs
c***
c     Li's ALGEBRAIC VERSION
c***  
c     ON INPUT:
c     
c     3D arrays (nx,ny,nz):
c
c       u: initial right hand side (divergence field)
c       p: anything
c     res: memory reserved for residual on the top grid
c       a: coefficients for equation at each level (computed in acfromac)
c       b: coefficients "   "
c       c: coefficients "   "
c      ei: inverse coefficients  "  "
c   
c     scalar parameters:
c 
c       r1, r2: parameters used in the physical version only
c     nx,ny,nz: dimensions of arrays
c       cycles: number of cycles at each V level
c
c     ON OUTPUT:
c     
c     everything the same, except p: pressure
c***********************************************************************
      SUBROUTINE mglin(r1,r2)
c***
      INCLUDE 'input.h'
      INTEGER j,jcycle,jj,jpost,jpre,nf,ngrid,nn,nxmin,nymin,nzmin
      INTEGER nnx,nny,nnz,nfx,nfy,nfz
      DOUBLE PRECISION r1, r2
      INTRINSIC MIN0
c     SUBROUTINES: copy, rstrct, abcei, fillzero, slvsml, interp,
c     relax, resid, addint
c***
c     get minimum dimension
c***
      nn  = MIN0(nx,ny,nz)
      nnx = nx
      nny = ny
      nnz = nz
      ngrid = NG
c***
c     restrict rhs and get coefficients
c***
      do while (nn .gt. 4)
         nn  =  nn/2 + 1
         nnx = nnx/2 + 1
         nny = nny/2 + 1
         nnz = nnz/2 + 1
         ngrid = ngrid - 1
         call rstrct(z(irhs(ngrid)),z(irhs(ngrid+1)),ibc,
     %        nnx,nny,nnz)
c***
         call acfromac(z(ia(ngrid)),z(ib(ngrid)),z(ic(ngrid)),
     %        z(iei(ngrid)),z(ia(ngrid+1)),z(ib(ngrid+1)),
     %        z(ic(ngrid+1)),nnx,nny,nnz)
      enddo
c***
c     define minimum values, just for clarity
c***
      nn = 4
      nxmin = nnx
      nymin = nny
      nzmin = nnz
      
      call fillzero(z(iu(1)),nnx,nny,nnz)
      call slvsml(z(iu(1)),z(irhs(1)),z(ia(1)),
     %     z(ib(1)),z(ic(1)),z(iei(1)),ibc,nnx,nny,nnz,
     %     NTOP,OMEGA,lengthscale)
c***
c     up with the V cycles
c***
      ngrid = NG
      do j=2,ngrid
         nn  = 2*nn  - 2
         nnx = 2*nnx - 2
         nny = 2*nny - 2
         nnz = 2*nnz - 2
         call interp(z(iu(j)),z(iu(j-1)),ibc,nnx,nny,nnz)
c***
c     cycles within a fixed V level
c***
         do jcycle=1,cycles
            nf  = nn
            nfx = nnx
            nfy = nny
            nfz = nnz
c***
c     down in the fixed V cycle
c***
            do jj=j,2,-1
               do jpre=1,NPRE
                  call relax(z(iu(jj)),lengthsize,z(irhs(jj)),z(ia(jj)),
     %                 z(ib(jj)),z(ic(jj)),z(iei(jj)),ibc,
     &                 nfx,nfy,nfz,OMEGA)
               enddo

               call resid(lengthsize,z(ires(jj)),z(iu(jj)),z(irhs(jj)),
     %         z(ia(jj)),z(ib(jj)),z(ic(jj)),ibc,nfx,nfy,nfz)
   
               nf  = nf/2  + 1
               nfx = nfx/2 + 1
               nfy = nfy/2 + 1
               nfz = nfz/2 + 1

               call rstrct(z(irhs(jj-1)),z(ires(jj)),ibc,nfx,nfy,nfz)
               call fillzero(z(iu(jj-1)),nfx,nfy,nfz)
            enddo
c***
c     lowest level in the fixed V cycle, indices to lowest value
c     (jusr for clarity, not really necessary)
c***            
            call slvsml(z(iu(1)),z(irhs(1)),z(ia(1)),z(ib(1)),z(ic(1)),
     %           z(iei(1)),ibc,nxmin,nymin,nzmin,NTOP,OMEGA,lengthsize)
            nf = 4
            nfx = nxmin
            nfy = nymin
            nfz = nzmin
c***
c     up in the fixed V cycle
c***  
            do jj=2,j
               nf  = 2*nf  - 2
               nfx = nfx*2 - 2
               nfy = nfy*2 - 2
               nfz = nfz*2 - 2
c***  
c     Use res for temporary storage inside addint
c***  
               call addint(z(iu(jj)),z(iu(jj-1)),z(ires(jj)),
     %              ibc,nfx,nfy,nfz)

               do jpost=1,NPOST
                  call relax(z(iu(jj)),lengthsize,z(irhs(jj)),z(ia(jj)),
     %                 z(ib(jj)),z(ic(jj)),z(iei(jj)),ibc,
     %                 nfx,nfy,nfz,OMEGA)
               enddo
            enddo   
         enddo
      enddo
c*** 
      return
      end
c***********************************************************************
      SUBROUTINE rstrct(uc,uf,ibc,ncx,ncy,ncz)
c***
      IMPLICIT none
      INTEGER ncx,ncy,ncz,ic,iff,jc,jf,kc,kf,ibc
      DOUBLE PRECISION uc(ncx,ncy,ncz),uf(2*ncx-2,2*ncy-2,2*ncz-2)
c     SUBROUTINES: bc_scalar
c***
c     each point inside the cube equal 1/8 of the sum of the vertices
c***
      do kc=2,ncz-1
         kf=2*kc-2
         do jc=2,ncy-1
            jf=2*jc-2
            do ic=2,ncx-1
               iff = 2*ic - 2
               uc(ic,jc,kc)= .125d0*( uf(iff,jf,kf) + uf(iff+1,jf,kf) 
     %            + uf(iff,jf+1,kf) + uf(iff+1,jf+1,kf)
     %            + uf(iff,jf,kf+1) + uf(iff+1,jf,kf+1)
     %            + uf(iff,jf+1,kf+1) + uf(iff+1,jf+1,kf+1) )
            enddo
         enddo
      enddo
c***
      call bc_scalar(uc,ibc,ncx,ncy,ncz)
c***
      return
      end
c***********************************************************************
      SUBROUTINE acfromac(a,b,c,ei,a2,b2,c2,nx,ny,nz)
c***
      IMPLICIT none
      INTEGER nx,ny,nz,i,j,k
      DOUBLE PRECISION a(nx,ny,nz),b(nx,ny,nz)
      DOUBLE PRECISION c(nx,ny,nz),ei(nx,ny,nz)
      DOUBLE PRECISION a2(2*nx-2,2*ny-2,2*nz-2)
      DOUBLE PRECISION b2(2*nx-2,2*ny-2,2*nz-2)
      DOUBLE PRECISION c2(2*nx-2,2*ny-2,2*nz-2)
c***
c     with this numerics, a is on the E face, b on the N face, 
c     c on the top, ei in the center of the cell
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=1,nx-1
               a(i,j,k) = 0.25d0*(a2(2*i-1,2*j-1,2*k-1) + 
     %              a2(2*i-1,2*j-2,2*k-1) + a2(2*i-1,2*j-1,2*k-2) +
     %              a2(2*i-1,2*j-2,2*k-2))
            enddo
         enddo
      enddo
      do k=2,nz-1
         do j=1,ny-1
            do i=2,nx-1
               b(i,j,k) = 0.25d0*(b2(2*i-1,2*j-1,2*k-1) +
     %              b2(2*i-2,2*j-1,2*k-1) + b2(2*i-1,2*j-1,2*k-2) +
     %              b2(2*i-2,2*j-1,2*k-2))
            enddo
         enddo
      enddo
      do k=1,nz-1
         do j=2,ny-1
            do i=2,nx-1
               c(i,j,k) = 0.25d0*(c2(2*i-1,2*j-1,2*k-1) +
     %              c2(2*i-2,2*j-1,2*k-1) + c2(2*i-1,2*j-2,2*k-1) +
     %              c2(2*i-2,2*j-2,2*k-1))
            enddo
         enddo
      enddo
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               ei(i,j,k) = 1.d0/(a(i,j,k) + a(i-1,j,k) + b(i,j,k) +
     %                           b(i,j-1,k) + c(i,j,k) + c(i,j,k-1))
            enddo
         enddo
      enddo
c***
      return
      end
c***********************************************************************
      SUBROUTINE slvsml(u,rhs,a,b,c,ei,ibc,nx,ny,nz,NTOP,OMEGA,
     &     	lengthsize)
c***
      IMPLICIT none
      INTEGER NTOP,j,nx,ny,nz,ibc
      DOUBLE PRECISION a(nx,ny,nz),b(nx,ny,nz),c(nx,ny,nz)
      DOUBLE PRECISION rhs(nx,ny,nz), u(nx,ny,nz), 
     %     OMEGA,lengthsize,ei(nx,ny,nz)
c     SUBROUTINES: relax
c***
c     call NTOP times relax at the lowest grid
c***  
      do j=1,NTOP
         call relax(u,lengthsize,rhs,a,b,c,ei,ibc,nx,ny,nz,OMEGA)
      enddo
      
c***
      return
      end
c***********************************************************************
      SUBROUTINE interp(uf,uc,ibc,nfx,nfy,nfz)
c***
      IMPLICIT none
      INTEGER nfx,nfy,nfz,ibc
      DOUBLE PRECISION uc(nfx/2+1,nfy/2+1,nfz/2+1),uf(nfx,nfy,nfz)
      INTEGER ic,if,jc,jf,kc,kf,ncx,ncy,ncz
c     SUBROUTINES: bc_scalar
c***
c     each point of the cube equal to the point inside
c***
      ncx=nfx/2 + 1
      ncy=nfy/2 + 1
      ncz=nfz/2 + 1
      do kc=2,ncz-1
         kf=2*kc-2
         do jc=2,ncy-1
            jf=2*jc-2
            do ic=2,ncx-1
               if=2*ic-2
               uf(if,jf,kf) = uc(ic,jc,kc)
               uf(if,jf+1,kf) = uc(ic,jc,kc)
               uf(if+1,jf,kf) = uc(ic,jc,kc)
               uf(if+1,jf+1,kf) = uc(ic,jc,kc)
               uf(if,jf,kf+1) = uc(ic,jc,kc)
               uf(if,jf+1,kf+1) = uc(ic,jc,kc)
               uf(if+1,jf,kf+1) = uc(ic,jc,kc)
               uf(if+1,jf+1,kf+1) = uc(ic,jc,kc)
            enddo
         enddo
      enddo
c***
      call bc_scalar(uf,ibc,nfx,nfy,nfz)
c***
      return
      end
c***********************************************************************
      SUBROUTINE resid(lengthsize,res,u,rhs,a,b,c,ibc,nx,ny,nz)
c***
      IMPLICIT none
      INTEGER nx,ny,nz,i,j,k,ibc
      DOUBLE PRECISION res(nx,ny,nz),rhs(nx,ny,nz),u(nx,ny,nz)
      DOUBLE PRECISION a(nx,ny,nz),b(nx,ny,nz),c(nx,ny,nz)
      DOUBLE PRECISION h1,h2i,lengthsize
c***
c     calculate the residue of the pressure equation
c***
      h1=lengthsize/(nx-2)
      h2i=1.d0/(h1*h1)
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               res(i,j,k)= - h2i*( a(i,j,k)*(u(i+1,j,k)-u(i,j,k)) 
     %                     + a(i-1,j,k)*(u(i-1,j,k)-u(i,j,k))    
     %                     + b(i,j,k)*(u(i,j+1,k)-u(i,j,k))     
     %                     + b(i,j-1,k)*(u(i,j-1,k) - u(i,j,k))
     %                     + c(i,j,k)*(u(i,j,k+1)-u(i,j,k))   
     %                     + c(i,j,k-1)*(u(i,j,k-1) - u(i,j,k))
     %                     ) + rhs(i,j,k)                    
            enddo
         enddo
      enddo
c***
      call bc_resid(lengthsize,res,ibc,nx,ny,nz)
c***
      return
      end

c***********************************************************************
      SUBROUTINE addint(uf,uc,res,ibc,nfx,nfy,nfz)
c***
      IMPLICIT none
      INTEGER nfx,nfy,nfz,i,j,k,ibc
      DOUBLE PRECISION res(nfx,nfy,nfz),uc(nfx/2+1,nfy/2+1,nfz/2+1),
     %                 uf(nfx,nfy,nfz)
c     SUBROUTINES: interp
c***
c     add the interpolated residue to the field uf
c***
      call interp(res,uc,ibc,nfx,nfy,nfz)
      do k=1,nfz
         do j=1,nfy
            do i=1,nfx
               uf(i,j,k)=uf(i,j,k)+res(i,j,k)
            enddo
         enddo
      enddo
c***
      return
      end
c***********************************************************************
      SUBROUTINE relax(u,lengthsize,rhs,a,b,c,ei,ibc,nx,ny,nz,OMEGA)
c***
      IMPLICIT none
      INTEGER nx,ny,nz,i,ipass,isw,j,jsw,k,ksw,ibc
      DOUBLE PRECISION rhs(nx,ny,nz),u(nx,ny,nz)
      DOUBLE PRECISION a(nx,ny,nz),b(nx,ny,nz),c(nx,ny,nz),
     %     ei(nx,ny,nz)
      DOUBLE PRECISION h1, h2, OMEGA, temp,lengthsize
c     SUBROUTINES: bc_scalar
c***
      h1=lengthsize/(nx-2)
      h2=h1*h1
      ksw=1
c***
      do ipass=1,2
c***
      call bc_scalar(u,ibc,nx,ny,nz)
c***
c     solve pressure equation with an alternate sweep
c***
         jsw=ksw
         do k=2,nz-1
           isw=jsw
           do j=2,ny-1
              do i=isw+1,nx-1,2
                 temp = (a(i,j,k)*u(i+1,j,k)+ a(i-1,j,k)*u(i-1,j,k)
     %                  + b(i,j,k)*u(i,j+1,k) + b(i,j-1,k)*u(i,j-1,k)
     %                  + c(i,j,k)*u(i,j,k+1) + c(i,j,k-1)*u(i,j,k-1)
     %                  - h2*rhs(i,j,k))*ei(i,j,k)
                 u(i,j,k) = OMEGA*temp + (1.d0 - OMEGA)*u(i,j,k)
              enddo
              isw=3-isw
           enddo
           jsw=3-jsw
         enddo
         ksw=3-ksw
      enddo
c***
      call bc_scalar(u,ibc,nx,ny,nz)
c***
      return
      end

c***********************************************************************
      SUBROUTINE finddiff(u,v,w,y,nx,ny,nz,imax,jmax,kmax)
c***
c     get the max difference, and its indices, between the components 
c     of a vector field
c***
      IMPLICIT none
      INTEGER imax, jmax, kmax, nx, ny, nz, i, j, k 
      DOUBLE PRECISION u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz),
     %                 y, res, t1, t2, t3
      INTRINSIC DABS,DMAX1 
c***
      y=0.d0
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               t1=DABS(u(i,j,k)-v(i,j,k))
               t2=DABS(u(i,j,k)-w(i,j,k))
               t3=DABS(v(i,j,k)-w(i,j,k))
               res = DMAX1(t1,t2,t3)
               if(res .gt. y) then
                  y = res
                  imax = i
                  jmax = j
                  kmax = k
               endif
            enddo
         enddo
      enddo
c***
      return
      end
c***********************************************************************
      SUBROUTINE copy(aout,ain,nx,ny,nz)
c***
c     copy ain in aout
c***
      IMPLICIT none
      INTEGER nx,ny,nz,i,j,k
      DOUBLE PRECISION ain(nx,ny,nz),aout(nx,ny,nz)
c***
      do k=1,nz
         do j=1,ny
            do i=1,nx
               aout(i,j,k)=ain(i,j,k)
            enddo
         enddo
      enddo
c***
      return
      end
c***********************************************************************


c***********************************************************************
c     The following sunroutines are all routines for the following B.C.'s:
c     periodicity in the x direction;
c     in the y and z directions: free-slip conditions on the wall (at
c     j=1,k=1), symmetry at the center (at j=ny, k=nz).
c     To include directly the four vertices of each horizontal plane,
c     the B.C.'s are imposed in the following way:
c     y-z limiting planes: k: 2->nz-1; j: 2->ny-1;
c     x-z limiting planes: k: 2->nz-1; i: 1->nx;
c     x-y limiting planes: j: 1->ny;   i: 1->nx;
c***********************************************************************
      SUBROUTINE bc_v()
c***
      INCLUDE 'input.h'
      INTEGER i, j, k
      DOUBLE PRECISION r,zz,yy
      DOUBLE PRECISION outlet, vout, outlet2, vout2

c***
c     periodicity of the 3 velocity components in the x direction
ctemp note that periodicity is not used here!!
c***

      do k=1,nz
         do j=ny,1,-1
	u(1,j,k)=xvel2*(-((j*1.d0-1.5d0)/(ny*1.d0-1.5d0)-1.d0)**2+1.d0)
     
	if (j.eq.1) u(1,1,k)=-u(1,2,k)
c		u(nx,j,k)=xvel2*(-((j-1.d0)/(ny*1.d0)-1.d0)**2+1.d0)
		v(nx,j,k)= v(nx-1,j,k) 
          w(nx,j,k)= w(nx-1,j,k) 
         enddo
      enddo

c***
c     y direction: free-slip on the wall for the parallel components,
c:    to zero the perpendicular one: 
c     grad u = grad w = 0;      v(.,1,.) = v(.,2,.) = 0;
c     symmetry at the center (w becomes the radial component!): 
c     grad u = grad w = 0;      v(.,ny,.) = 0 
c***
      do k=2,nz-1
         do i=2,nx-1
			u(i,1,k ) = -u(i,2,k)
			u(i,ny,k) = u(i,ny-1,k)
			v(i,2 ,k) = 0.d0
			v(i,1 ,k) = 0.d0
			v(i,ny,k) = 0.d0
			w(i,1,k ) = -w(i,2,k )
			w(i,ny,k) = w(i,ny-1,k)	
         enddo 
      enddo 
c***
c     z direction: same as before with v <-> w 
c     grad u = grad v = 0;      w(.,.,1) = w(.,.,2) = 0;
c     grad u = grad v = 0;      w(.,.,nz) = 0 
c***
      do j=1,ny
         do i=2,nx-1
			u(i,j,1 ) = u(i,j,2)
			u(i,j,nz) = u(i,j,nz-1)
			v(i,j,1 ) = v(i,j,2)
		    v(i,j,nz) = v(i,j,nz-1)
			w(i,j,2 ) = 0.d0
			w(i,j,nz)=0.d0
         enddo 
      enddo 
c***
      outlet=0.
      outlet2=0.
      do k=2,nz-1
         do j=2,ny-1
               u(2,j,k)= u(1,j,k)
               outlet = outlet + u(1,j,k)
c compute total flux
               outlet2 = outlet2 + u(nx-1,j,k)
         enddo 
      enddo
c average entrance flux
      vout = outlet/ (ny-2)/(nz-2)
c average exit flux
      vout2 = outlet2/ (ny-2)/(nz-2)
      do k=2,nz-1
        do j=2,ny-1
            u(nx,j,k) = u(nx-1,j,k) + vout - vout2
        enddo
      enddo
      return
      end
c***********************************************************************
      SUBROUTINE bc_w(w1,w2,w3)
c***
      INCLUDE 'input.h'
      INTEGER i, j, k
      DOUBLE PRECISION w1(nx,ny,nz), w2(nx,ny,nz), w3(nx,ny,nz)
      DOUBLE PRECISION r,zz,yy
      DOUBLE PRECISION outlet, vout, outlet2, vout2

c***
c     periodicity of the 3 velocity components in the x direction
c***
      do k=2,nz-1
         do j=2,ny-1
          w2(1,j,k)=v(1,j,k)
          w2(2,j,k)=w2(1,j,k) 
          w2(nx,j,k)=v(nx,j,k)
          w3(1,j,k)=w(1,j,k)
          w3(2,j,k)=w3(1,j,k) 
          w3(nx,j,k)=w(nx,j,k)
         enddo
      enddo
c***
c     y direction: free-slip on the wall for the parallel components,
c:    to zero the perpendicular one: 
c     grad u = grad w = 0;      v(.,1,.) = v(.,2,.) = 0;
c     symmetry at the center (w becomes the radial component!): 
c     grad u = grad w = 0;      v(.,ny,.) = 0 
c***
      do k=2,nz-1
         do i=2,nx-1
            w1(i,1,k) = -w1(i,2,k)
            w1(i,ny,k) = w1(i,ny-1,k)
            w2(i,1 ,k) =0.d0
            w2(i,2 ,k) =0.d0
            w2(i,ny,k) =v(i,ny,k)
            w3(i,1 ,k) = -w3(i,2,k)
            w3(i,ny,k) = w3(i,ny-1,k)
         enddo 
      enddo 
c***
c     z direction: same as before with v <-> w 
c     grad u = grad v = 0;      w(.,.,1) = w(.,.,2) = 0;
c     grad u = grad v = 0;      w(.,.,nz) = 0 
c***
      do j=1,ny
         do i=2,nx-1
            w1(i,j,1 ) = w1(i,j,2)
            w1(i,j,nz) = w1(i,j,nz-1)
            w2(i,j,1 ) = w2(i,j,2)
            w2(i,j,nz) = w2(i,j,nz-1)
            w3(i,j,1 ) = w(i,j,1)
            w3(i,j,2 ) = w3(i,j,1)
            w3(i,j,nz) = w(i,j,nz)
         enddo 
      enddo 
c***
      outlet=0.
      outlet2=0.
      do k=2,nz-1
         do j=2,ny-1
               w1(1,j,k)=u(1,j,k)
               w1(2,j,k)=w1(1,j,k)
               outlet = outlet + u(1,j,k)
c compute total flux
               outlet2 = outlet2 + w1(nx-1,j,k)
         enddo 
      enddo
c average entrance flux
      vout = outlet/ (ny-2)/(nz-2)
c average exit flux
      vout2 = outlet2/ (ny-2)/(nz-2)
      do k=2,nz-1
        do j=2,ny-1
            w1(nx,j,k) = w1(nx-1,j,k) + vout - vout2
        enddo
      enddo
      return
      end
c***********************************************************************
      SUBROUTINE bc_tensor()
c***
      INCLUDE 'input.h'
      INTEGER i, j, k
c***
c     periodicity in the x direction of all components of the stress 
c     tensor [ Note: there are a few commented out, this is because 
c                    they are not used in the code ]
c***
      do k=2,nz-1
         do j=2,ny-1
            S11(1 ,j,k) = 0.d0
c---            S11(nx,j,k) = S11(2   ,j,k)
c---            S22(1 ,j,k) = S22(nx-1,j,k)
c---            S22(nx,j,k) = S22(2   ,j,k)
c---            S33(1 ,j,k) = S33(nx-1,j,k)
c---            S33(nx,j,k) = S33(2   ,j,k)
c---            S12(1 ,j,k) = S12(nx-1,j,k)
            S12(nx,j,k) = 0.d0
c---            S13(1 ,j,k) = S13(nx-1,j,k)
            S13(nx,j,k) = 0.d0
c---            S23(1 ,j,k) = S23(nx-1,j,k)
c---            S23(nx,j,k) = S23(2   ,j,k)
         enddo 
      enddo
c***
c     y direction: free slip at the wall and simmetry at the center 
c     for indices 11, 33, 13 (this means zero gradient for both 
c     conditions); to zero all the rest
c***
      do k=2,nz-1
         do i=1,nx
c---            S11(i,1 ,k) = S11(i,2   ,k)
c---            S11(i,ny,k) = S11(i,ny-1,k)
            S22(i,1 ,k) =  S22(i,2 ,k)
c---            S22(i,ny,k) = 0.d0
c---            S33(i,1 ,k) = S33(i,2,k)
c---            S33(i,ny,k) = S33(i,ny-1,k)
c---            S12(i,1 ,k) = 0.d0
c            S12(i,2 ,k) = S12(i,2 ,k)
            S12(i,ny,k) = 0.d0
c---            S13(i,1 ,k) = S13(i,2,k)
c---            S13(i,ny,k) = S13(i,ny-1,k)
c---            S23(i,1 ,k) = 0.d0
c            S23(i,2 ,k) = 0.d0
            S23(i,ny,k) = 0.d0
         enddo 
      enddo
c***
c     z direction: free slip at the wall and simmetry at the center 
c     for indices 11, 22, 12 (this means zero gradient for both 
c     conditions); to zero all the rest
c***
      do j=1,ny
         do i=1,nx
c---            S11(i,j,1 ) = S11(i,j,2)
c---            S11(i,j,nz) = S11(i,j,nz-1)
c---            S22(i,j,1 ) = S22(i,j,2)
c---            S22(i,j,nz) = S22(i,j,nz-1)
            S33(i,j,1 ) = 0.d0
			
c---            S33(i,j,nz) = 0.d0
c---            S12(i,j,1 ) = S12(i,j,2)
c---            S12(i,j,nz) = S12(i,j,nz-1)
c---            S13(i,j,1 ) = 0.d0
            S13(i,j,2 ) = 0.d0
            S13(i,j,nz) = 0.d0
c---            S23(i,j,1 ) = 0.d0
            S23(i,j,2 ) = 0.d0
            S23(i,j,nz) = 0.d0
         enddo 
      enddo
c***
      return
      end
c***********************************************************************
      SUBROUTINE bc_scalar(scal,ibc,nx,ny,nz)
c***
      IMPLICIT none
      INTEGER nx,ny,nz, i, j, k,ibc
      DOUBLE PRECISION scal(nx,ny,nz)
c***  
c     periodicity in the x direction
c***  
      do k=2,nz-1
         do j=2,ny-1
            scal(1 ,j,k) = scal(2,j,k)
            scal(nx,j,k) = 0.d0 !scal(nx-1,j,k)
         enddo
      enddo
c***  
c     gradient equal to zero in the y and z directions
c***  
      do k=2,nz-1
         do i=1,nx
            scal(i,1,k ) = scal(i,2   ,k)
            scal(i,ny,k) = scal(i,ny-1,k)
         enddo
      enddo 
      do j=1,ny
         do i=1,nx
            scal(i,j,1 ) = scal(i,j,2)
            scal(i,j,nz) = scal(i,j,nz-1)
         enddo
      enddo 
c***
      return
      end
c***********************************************************************
      SUBROUTINE bc_c(cc)
c***
      INCLUDE 'input.h'
      INTEGER i, j, k
      DOUBLE PRECISION cc(nx,ny,nz),ccf(nx,ny,nz),w22(nx,ny,nz)
	DOUBLE PRECISION zz,yy,r,normv,theta,nf,mmy,mmx,mmz,mx,mz,sign
c***  
c     periodicity in the x direction
c***  
      do k=1,nz
         do j=1,ny
            cc(2,j,k) = cc(1,j,k)
            cc(nx,j,k) = cc(nx-1,j,k)
         enddo 
      enddo

	call copy(ccf,cc,nx,ny,nz)
	nf=1
	do i=1,nf
	call filterc27(ccf,w22,nx,ny,nz)
	call copy (ccf,w22,nx,ny,nz)
	call bc_out(ccf,nx,ny,nz)
	enddo	 

c***  
c     gradient equal to zero in the y and z directions
c*** 	
	do k=2,nz-1
         do i=2,nx-1
c	**** determination of theta as a function of Vcl *****
	j=2
	sign=(ccf(i+1,j,k)-ccf(i-1,j,k))*u(i,j,k)+
     &(ccf(i,j,k+1)-ccf(i,j,k-1))*w(i,j,k)
	normv=dsqrt(u(i,j,k)**2+w(i,j,k)**2)
	sign=-sign
	if (sign.gt.0) sign=1
	if (sign.lt.0) sign=-1
	if (normv.gt.0.1) then
		if (sign.eq.1)	theta=115
		if	(sign.eq.-1) theta=45
	else
	theta=350.0*normv*sign+80	
	endif
c	******************************************************
	mmy=dcosd(theta)
c            cc(i,1,k ) = cc(i,2   ,k)
            
		if (cc(i,2,k).lt.0.1 .OR.cc(i,2,k).gt.0.9 ) then
			cc(i,1,k)=cc(i,2,k)
		else	
		cc(i,1,k)=cc(i,2,k) !+h*mmy
		endif
		if (cc(i,1,k).gt.1 ) cc(i,1,k)=1
		if (cc(i,1,k).lt.0 ) cc(i,1,k)=0
		cc(i,ny,k) = cc(i,ny-1,k)
         enddo
      enddo 
      do j=1,ny
         do i=2,nx-1
            cc(i,j,1 ) = cc(i,j,2)
            cc(i,j,nz) = cc(i,j,nz-1)
         enddo
      enddo 
c***
      return
      end

c	*****************************************************
      SUBROUTINE bc_out(u1,nx,ny,nz)
c***

      INTEGER i, j, k,nx,ny,nz
      DOUBLE PRECISION u1(nx,ny,nz)
c***  
c     periodicity in the x direction
c***  
      do k=1,nz
         do j=1,ny
            u1(1,j,k) = u1(2,j,k)
            u1(nx,j,k) = u1(nx-1,j,k)
         enddo 
      enddo
c***  
c     gradient equal to zero in the y and z directions
c***  
      do k=2,nz-1
         do i=1,nx
            u1(i,1,k ) = u1(i,2   ,k)
            u1(i,ny,k) = u1(i,ny-1,k)
         enddo
      enddo 
      do j=1,ny
         do i=2,nx-1
            u1(i,j,1 ) = u1(i,j,2)
            u1(i,j,nz) = u1(i,j,nz-1)
         enddo
      enddo 
c***
      return
      end
c***********************************************************************
      SUBROUTINE bc_resid(lengthsize,scal,ibc,nx,ny,nz)
c***
      IMPLICIT none
      INTEGER nx,ny,nz, i, j, k,ibc
      DOUBLE PRECISION scal(nx,ny,nz),lengthsize
c***  
c     periodicity in the x direction
c***  
      do k=1,nz
         do j=1,ny
            scal(1 ,j,k) = scal(2,j,k)
            scal(nx,j,k) = scal(nx-1,j,k)
         enddo
      enddo
c***  
c     residue equal to zero in the y and z directions
c***  
      do k=2,nz-1
         do i=2,nx-1
            scal(i,1 ,k) = 0.d0
            scal(i,ny,k) = 0.d0
         enddo
      enddo
      do j=1,ny
         do i=2,nx-1
            scal(i,j,1 ) = 0.d0
            scal(i,j,nz) = 0.d0
         enddo
      enddo
c***  
      return
      end
c***********************************************************************
      SUBROUTINE bc_coef(a,b,c,cc,r1,r2,ibc,nx,ny,nz)
c***
      IMPLICIT none
      INTEGER nx, ny, nz, i, j, k,ibc
      DOUBLE PRECISION a(nx,ny,nz), b(nx,ny,nz), c(nx,ny,nz)
      DOUBLE PRECISION cc(nx,ny,nz), r1, r2
c***
c     need to define indices i=1 for a, j=1 for b and k=1 for c: use
c     periodic boundary conditions in x direction, simple average
c     in y and z directions
c***
      do k=1,nz
         do j=1,ny
            a(1 ,j,k) = a(2,j,k)
            a(nx,j,k) = a(nx-1,j,k)
         enddo
      enddo
      do k=2,nz-1
         do i=2,nx-1
            b(i,1,k) = 1.d0 / (r1*(cc(i,1,k) + cc(i,2,k)) + r2)
         enddo
      enddo
      do j=1,ny
         do i=2,nx-1
            c(i,j,1) = 1.d0 / (r1*(cc(i,j,1) + cc(i,j,2)) + r2)
         enddo
      enddo
c***
      return 
      end
c***********************************************************************
      SUBROUTINE bc_flux(vof1,vof3,vel,indx)
c***
      INCLUDE 'input.h'
      INTEGER indx, i, j, k
      DOUBLE PRECISION vof1(nx,ny,nz), vof3(nx,ny,nz), vel(nx,ny,nz)
c***
c     along x direction
c***
      if (indx .eq. 1) then
         do k=2,nz-1
            do j=2,ny-1
               vof3(1 ,j,k) = 0.d0
               vof1(nx,j,k) = 0.d0
            enddo 
         enddo 
c***
c     along y direction
c***
      elseif (indx .eq. 2) then
         do k=2,nz-1
            do i=2,nx-1
               vof3(i,1 ,k) = 0.d0 !DMAX1(0.0d0,vel(i,2,k)*tauh)
               vof1(i,ny,k) = 0.d0
            enddo
         enddo
c***
c     along z direction
c***
      elseif (indx .eq. 3) then
         do j=2,ny-1
            do i=2,nx-1
               vof3(i,j,1) = 0.d0 !DMAX1(0.0d0,vel(i,j,2)*tauh)
               vof1(i,j,nz) = 0.0d0
            enddo
         enddo
c***  
c     wrong value
c***  
      else 
         stop 'wrong value for indx in one of the swaps'
      endif 
c***
      return
      end
c***********************************************************************
      INTEGER FUNCTION MYMALLOC(len)
c***
      INCLUDE 'input.h'
      INTEGER len
c***
      if ((mem+len+1) .gt. MEMLEN)  then 
         write(*,*) 'MEMLEN=', MEMLEN,', mem+len+1:', mem+len+1
         write(*,*) 'insufficient memory in MYMALLOC'
         stop
      endif
      z(mem+1) = len
      MYMALLOC = mem + 2
      mem = mem + len + 1
c***
      return
      end

c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
c PROGRAM TO FIND alpha IN: m1 x1 + m2 x2 + m3 x3 = alpha,
c GIVEN m1+m2+m3=1 (all > 0) AND THE VOLUMETRIC FRACTION cc
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
      DOUBLE PRECISION FUNCTION AL3D(b1,b2,b3,cc)
c***
      IMPLICIT none
      DOUBLE PRECISION m1,m2,m3,cc,b1,b2,b3,tmp,pr,ch,mm,m12
      DOUBLE PRECISION p,p12,q,teta,cs
      DOUBLE PRECISION UNTIER,V1,V2,V3,const
      PARAMETER (UNTIER=1.d0/3.d0)
      INTRINSIC DMAX1,DMIN1,DSQRT,DACOS,DCOS
c***  
c     (1) order coefficients: m1<m2<m3; (2) get ranges: V1<V2<v3;
c     (3) limit ch (0.d0 < ch < 0.5d0); (4) calculate alpha
c*(1)* 
      m1 = DMIN1(b1,b2)
      m3 = DMAX1(b1,b2)
      m2 = b3
      if (m2 .LT. m1) then
         tmp = m1
         m1 = m2
         m2 = tmp
      else if (m2 .GT. m3) then
         tmp = m3
         m3 = m2
         m2 = tmp
      endif
c*(2)*
      m12 = m1 + m2 
      pr  = DMAX1(6.d0*m1*m2*m3,1.d-50)
      V1  = m1*m1*m1/pr
      V2  = V1 + 0.5d0*(m2-m1)/m3
      if (m3 .LT. m12) then
         mm = m3
         V3 = (m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) +
     %        m2*m2*(m2-3.d0*m3))/pr
      else
         mm = m12
         V3 = 0.5d0*mm/m3
      endif
c*(3)*
      ch = DMIN1(cc,1.d0-cc)
c*(4)*      
      if (ch .LT. V1) then
c***         AL3D = cbrt(pr*ch)
         AL3D = (pr*ch)**UNTIER
      else if (ch .LT. V2) then
         AL3D = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(ch-V1)))
      else if (ch .LT. V3) then
         p = 2.d0*m1*m2
         q = 1.5d0*m1*m2*(m12 - 2.d0*m3*ch)
         p12 = DSQRT(p)
	   const=q/(p*p12)
	   if(const.gt.1.0) const=1.
	   if(const.lt.-1.0) const=-1.
         teta = DACOS(const)/3.d0
         cs = DCOS(teta)
         AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + m12
      else if (mm .EQ. m12) then
         AL3D = m3*ch + 0.5d0*mm
      else 
         p = m1*(m2+m3) + m2*m3 - 0.25d0
         q = 1.5d0*m1*m2*m3*(0.5d0-ch)
         p12 = DSQRT(p)
	   const=q/(p*p12)
	   if(const.gt.1.0) const=1.
	   if(const.lt.-1.0) const=-1.
         teta = DACOS(const)/3.d0
         cs = DCOS(teta)
         AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + 0.5d0
      endif

      if (cc .GT. 0.5d0)  AL3D = 1.d0 - AL3D
c***
      return
      end
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
c PROGRAM TO FIND THE "CUT VOLUME" V0 GIVEN r0, dr0 AND
c m1 x1 + m2 x2 + m3 x3 = alpha
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
      DOUBLE PRECISION FUNCTION FL3D(m1,m2,m3,alpha,r0,dr0)
c***
      IMPLICIT none
      DOUBLE PRECISION m1,m2,m3,alpha,r0,dr0
      DOUBLE PRECISION al,al0,n1,n2,n3,b1,b2,b3,b12,bm,tmp,pr
      INTRINSIC DMAX1,DMIN1,DABS
c***
c     (1) move origin to r0 along r ;  (2) reflect parallelepiped;
c     (3) limit alpha (0<= al0 <=0.5); (4) order coefficients: b1<b2<b3;
c     (5) calculate volume (NOTE: it is assumed:s0=t0=0; ds0=dt0=1.)
c*(1)*
      al = alpha - m1*r0
c*(2)*
      al = al + DMAX1(0.d0,-m1*dr0)+DMAX1(0.d0,-m2)+DMAX1(0.d0,-m3)
      tmp = DABS(m1)*dr0 + DABS(m2) + DABS(m3)
      n1 = DABS(m1)/tmp
      n2 = DABS(m2)/tmp
      n3 = DABS(m3)/tmp
      al = DMAX1(0.d0,DMIN1(1.d0,al/tmp))
c*(3)*
      al0 = DMIN1(al,1.d0-al)
c*(4)*
      b1 = DMIN1(n1*dr0,n2)
      b3 = DMAX1(n1*dr0,n2)
      b2 = n3
      if (b2 .LT. b1) then
         tmp = b1
         b1 = b2
         b2 = tmp
      else if (b2 .GT. b3) then
         tmp = b3
         b3 = b2
         b2 = tmp
      endif
      b12 = b1 + b2
      bm = DMIN1(b12,b3)
      pr = DMAX1(6.d0*b1*b2*b3,1.0d-50)
c*5*     
      if (al0 .LT. b1) then
         tmp = al0*al0*al0/pr
      else if (al0 .LT. b2) then
         tmp = 0.5d0*al0*(al0-b1)/(b2*b3) +  b1*b1*b1/pr
      else if (al0 .LT. bm) then
         tmp = (al0*al0*(3.d0*b12-al0) + b1*b1*(b1-3.d0*al0) +
     %        b2*b2*(b2-3.d0*al0))/pr
      else if (bm .EQ. b12) then
         tmp = (al0 - 0.5d0*bm)/b3
      else 
         tmp = (al0*al0*(3.d0-2.d0*al0) + b1*b1*(b1-3.d0*al0) +
     %        b2*b2*(b2-3.d0*al0) + b3*b3*(b3-3.d0*al0))/pr
      endif
      
      if (al .LE. 0.5d0) then
         FL3D = tmp*dr0
      else
         FL3D = (1.d0-tmp)*dr0
      endif
c***  
      return
      end
c***********************************************************************
      DOUBLE PRECISION FUNCTION FINDMAX(u,nx,ny,nz,imax,jmax,kmax)
c***
c     get the max value, and its indices, of a scalar field
c***
      IMPLICIT none
      INTEGER imax, jmax, kmax, nx, ny, nz, i, j, k 
      DOUBLE PRECISION u(nx,ny,nz)
c***
      FINDMAX = u(2,2,2)
      imax = 2
      jmax = 2
      kmax = 2
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               if (u(i,j,k) .gt. FINDMAX) then
                  FINDMAX = u(i,j,k)
                  imax = i
                  jmax = j
                  kmax = k
               endif
            enddo
         enddo
      enddo
c***
      return
      end
c***********************************************************************
      DOUBLE PRECISION FUNCTION FMODMAX(u,nx,ny,nz,imax,jmax,kmax)
c***
c     get the absolute max value, and its indices, of a scalar field
c***
      IMPLICIT none
      INTEGER imax, jmax, kmax, nx, ny, nz, i, j, k 
      DOUBLE PRECISION u(nx,ny,nz), res
      INTRINSIC DABS
c***
      FMODMAX=-1.d0
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               res = DABS(u(i,j,k))
               if (res .gt. FMODMAX) then
                  FMODMAX = res
                  imax = i
                  jmax = j
                  kmax = k
               endif
            enddo
         enddo
      enddo
c***
      return
      end
c***********************************************************************
      DOUBLE PRECISION FUNCTION FINDMIN(u,nx,ny,nz,imin,jmin,kmin)
c***
c     get the min value, and its indices, of a scalar field
c***
      IMPLICIT none
      INTEGER imin, jmin, kmin, nx, ny, nz, i, j, k 
      DOUBLE PRECISION u(nx,ny,nz)
c***
      FINDMIN = u(2,2,2)
      imin = 2
      jmin = 2
      kmin = 2
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               if (u(i,j,k) .lt. FINDMIN) then
                  FINDMIN = u(i,j,k)
                  imin = i
                  jmin = j
                  kmin = k
               endif
            enddo
         enddo
      enddo
c***
      return
      end
c***********************************************************************
      DOUBLE PRECISION FUNCTION SUMFIELD(u,nx,ny,nz)
c***
c     get the sum of a scalar field
c***
      IMPLICIT none
      INTEGER nx, ny, nz, i, j, k
      DOUBLE PRECISION u(nx,ny,nz)
c***
      SUMFIELD=0.d0
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               SUMFIELD = SUMFIELD + u(i,j,k) 
            enddo
         enddo
      enddo
c***
      return
      end

c***********************************************************************
      DOUBLE PRECISION FUNCTION FINDN2(u,v,w,nx,ny,nz)
c***
c     get the average module of a vector field
c***
      IMPLICIT none
      INTEGER nx, ny, nz, i, j, k
      DOUBLE PRECISION u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), tnpt
      INTRINSIC DSQRT
c***
      FINDN2=0.d0
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               FINDN2 = FINDN2 + u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2
            enddo
         enddo
      enddo 
      tnpt = DFLOAT((nx-2)*(ny-2)*(nz-2))
      FINDN2 = DSQRT(FINDN2/tnpt)
c***
      return
      end
c***********************************************************************
      DOUBLE PRECISION FUNCTION FINDNINF(u,v,w,nx,ny,nz)
c***
c     get the max module of a vector field
c***
      IMPLICIT none
      INTEGER nx, ny, nz, i, j, k
      DOUBLE PRECISION u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), vv
      INTRINSIC DSQRT,DMAX1
c***
      FINDNINF=0.d0
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               vv = DSQRT(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2)
               FINDNINF = DMAX1(vv,FINDNINF)
            enddo
         enddo
      enddo
c***
      return
      end
      DOUBLE PRECISION FUNCTION CFL()
c***
c    check the CFL condition for the interface motion
c***
      INCLUDE 'input.h'
      INTEGER i, j, k
      INTRINSIC DMAX1, DABS
c***
      CFL = 0.0d0
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               CFL = DMAX1(CFL,DABS(u(i,j,k)))
            enddo
         enddo
      enddo
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               CFL = DMAX1(CFL,DABS(v(i,j,k)))
            enddo
         enddo
      enddo
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
               CFL = DMAX1(CFL,DABS(w(i,j,k)))
            enddo
         enddo
      enddo
      CFL=CFL*tauh
c***
      return
      end
c***********************************************************************
      SUBROUTINE printdata()
c***
      INCLUDE 'input.h'
c***
      write (*,*)
      write (*,*) '*** input/output files ***'
      write (*,*)
      call writeminmax('cc',z(icc(ng)),nx,ny,nz)
      call writeminmax('u',u,nx,ny,nz)
      call writeminmax('v',v,nx,ny,nz)
      call writeminmax('w',w,nx,ny,nz)
      call writeminmax('p',z(iu(ng)),nx,ny,nz)
      write (*,*)
c***
c     write statemnts for input data on the sreen

      write (6,*) "***** numerical properties *********"     
      write(6,*) "nx,ny,nz=",nx,ny,nz
      write(6,*) "maxtimestep=",tmax
      write(6,*) "nfilter=",nfilter
      write(6,*) "MG_cycles_bound=",cycles
      write(6,*) "OMEGA=",OMEGA
      write(6,*) "NPRE=",NPRE
      write(6,*) "NPOST=",NPOST
      write(6,*) "NTOP=",NTOP
	write(6,*) "restart=", restart
	write(6,*) "maxcfl=", maxcfl
	write(6,*) "iturb=",iturb
	write(6,*) "convectype=",convectype
	write(6,*) "hswitch=",hswitch
			
      write(6,*)"***** physical and flow properties **********"
      write(6,*) "xvel1-si=",xvel1_SI
      write(6,*) "yvel1-si=",yvel1_SI
      write(6,*) "zvel1-si=",zvel1_SI
      write(6,*) "xvel2-si=",xvel2_SI
      write(6,*) "yvel2-si=",yvel2_SI
      write(6,*) "zvel2-si=",zvel2_SI
      write(6,*) "rho2-si=",rho2_si
      write(6,*) "rho1si=",rho1_si
      write(6,*) "sigma-si=",sigma_si
      write(6,*) "visc2_si=",visc2_si
      write(6,*) "visc1_si=",visc1_si
      write(6,*) "z-grav_si=",grav_si
      write(6,*) "x-grav_si=",grav1_si

      write(6,*) "****** geometry and flow type *******"
      write(6,*) "ijet=",ijet
      write(6,*) "idrop=",idrop
      write(6,*) "ibubble=",ibubble
      write(6,*) "length_SI",length_SI   
      write(6,*) "radius_si=",radius_si
      write(6,*) "ibc=",ibc
c      write(6,*) "x0_SI=",x0_SI
c	write(6,*) "y0_SI=",y0_SI
c	write(6,*) "z0_SI=",z0_SI
      
	write(6,*) "****** print data ********"
      write(6,*) "tprint=",tprint
      write(6,*) "ntecplot=",ntecplot
	write(6,*) "tsave=", tsave

c     write statemnts for input data on the data file
      
      write(23,*) "***** numerical properties *********"     
      write(23,*) "nx,ny,nz=",nx,ny,nz
      write(23,*) "tmax=",tmax
      write(23,*) "nfilter=",nfilter
      write(23,*) "MG_cycles_bound=",cycles
c      write(23,*) "div_tolerance=",mineps
      write(23,*) "OMEGA=",OMEGA
      write(23,*) "NPRE=",NPRE
      write(23,*) "NPOST=",NPOST
      write(23,*) "NTOP=",NTOP
	write(23,*) "restart=", restart
	write(23,*) "maxcfl=", maxcfl
      write(23,*) "iturb=",iturb
      write(23,*) "convectype=",convectype
	write(23,*) "hswitch=",hswitch
			
      write(23,*)"***** physical and flow properties **********"
      write(23,*) "xvel1-si=",xvel1_SI
      write(23,*) "yvel1-si=",yvel1_SI
      write(23,*) "zvel1-si=",zvel1_SI
      write(23,*) "xvel2-si=",xvel2_SI
      write(23,*) "yvel2-si=",yvel2_SI
      write(23,*) "zvel2-si=",zvel2_SI
      write(23,*) "rho2-si=",rho2_si
      write(23,*) "rho1si=",rho1_si
      write(23,*) "sigma-si=",sigma_si
      write(23,*) "visc2_si=",visc2_si
      write(23,*) "visc1_si=",visc1_si
      write(23,*) "z-grav_si=",grav_si
      write(23,*) "x-grav_si=",grav1_si

      write(23,*) "****** geometry and flow type *******"
      write(23,*) "ijet=",ijet
      write(23,*) "idrop=",idrop
      write(23,*) "ibubble=",ibubble
      write(23,*) "radius_si=",radius_si
      write(23,*) "ibc=",ibc
c      write(23,*) "x0_SI=",x0_SI
c	write(23,*) "y0_SI=",y0_SI
c	write(23,*) "z0_SI=",z0_SI
      
	write(23,*) "****** print data ********"
      write(23,*) "tprint=",tprint
      write(23,*) "ntecplot=",ntecplot
	write(23,*) "tsave=", tsave
c***
      write (*,*)
      write (*,*) '*** non-dimensional physics parameters ***'
      write (*,*)
c***
      write (*,*) 'time step         =', tau
      write (*,*) 'grid size         =', h
      write (*,*) 'density fluid 1   =', rho1
      write (*,*) 'density fluid 2   =', rho2
      write (*,*) 'viscosity fluid 1 =', mu1
      write (*,*) 'viscosity fluid 2 =', mu2
      write (*,*) 'surface tension   =', sigma
      write (*,*) 'gravity           =', g
c***
      write (23,*) 'time step         =', tau
      write (23,*) 'grid size         =', h
      write (23,*) 'density fluid 1   =', rho1
      write (23,*) 'density fluid 2   =', rho2
      write (23,*) 'viscosity fluid 1 =', mu1
      write (23,*) 'viscosity fluid 2 =', mu2
      write (23,*) 'surface tension   =', sigma
      write (23,*) 'gravity           =', g
      write (*,*)
	write (23,*) 'No. Cell in D=',2*radius/h
c***
      return
      end
c***********************************************************************
      SUBROUTINE printout(t,cc,p)
c***
      INCLUDE 'input.h'
      INTEGER t,i,j,k
	DOUBLE PRECISION cc(nx,ny,nz),p(nx,ny,nz)
      DOUBLE PRECISION velmax,SUMFIELD,FINDNINF,FINDMAX,CFL,CFL1,
     &	maxmut,meanmut,totalcc,totalxmom,pn,yy,xx,norm,findn2,tsi
	CHARACTER (128) FileName 

      EXTERNAL SUMFIELD,FINDNINF,CFL,FINDMAX,findn2
c     Checkpoint every tprint steps
c***
c      if (mod(t,10).eq.0) call uplot(t,cc,p)
	if (MOD(t,tprint) .eq. 0) then
      tsi=time*timescale
c	if (mod(t,10).eq.0) call printtheta(cc,t)
c	call thetahysteresis(cc,nx,ny,nz,tsi)
	call velocityc(cc)	
         write(*,*)'***************************************************
     &*************************'
         write(23,*)'**************************************************
     &*************************'
	   write (*,*) ' timestep=',t
	   write(*,*)   ' time=',time,'  time_SI=', time*timescale
	   write (23,*) ' timestep=',t
	   write(23,*)   ' time=',time,'  time_SI=', time*timescale
	   write (*,*) ' tau=',tau,'   tau_SI=',tau*timescale
	   write (23,*) '  tau=',tau,'   tau_SI=',tau*timescale
	   write (*,*) ' CFL = ',CFL()
	   write (23,*) ' CFL = ',CFL()
         velmax = FINDNINF(u,v,w,nx,ny,nz)
	   totalcc=SUMFIELD(z(icc(ng)),nx,ny,nz)
	   totalxmom=SUMFIELD(u,nx,ny,nz)
         write (*,*) 'max |v|= ',velmax,'total cc = ',totalcc
         write (23,*)'max |v|= ',velmax,'total cc = ',totalcc
         write (*,*) 'total x momentum= ',totalxmom
         write (23,*)'total x momentum= ',totalxmom
c         call minmaxvcp(z(iu(ng)),z(icc(ng)),t)
         maxmut=findmax(mut,nx,ny,nz,i,j,k)
         meanmut=sumfield(mut,nx,ny,nz)
	   maxmut=maxmut/mu1
	   meanmut=meanmut/((nx-2)*(ny-2)*(nz-2))/mu1
	write(6,*) 'max mut/mu1=',maxmut,' mean mut/mu1=',meanmut
	write(23,*) 'max mut/mu1=',maxmut,' mean mut/mu1=',meanmut
	
	write(24,111) time*timescale,velmax*lengthscale/timescale
c	if (velmax*lengthscale/timescale.gt.500) stop
	norm=findn2(u,v,w,nx,ny,nz)
	write(26,111) time*timescale,norm*lengthscale/timescale
	endif
	
      if(mod(t,ntecplot).eq.0) then
      call printtecp(t,z(icc(ng)))
      endif
	
c		************plot pressure***********	
c	open(25,file='pressure.plt')

	if (mod(t,100).eq.0) then
	write (FileName, '(A,I4.4,A)') 'pressure', t/100, '.plt' 
      OPEN(25,FILE=FileName) 
	write(25,*) 'VARIABLES="x","y","p"'
	write(25,*) 'zone i=',nx-2,'j=',ny-2
	k=(nz-2)/2
      do j= 2,ny-1
      do i= 2,nx-1

      yy= h*j
      xx= h*i

      pn=p(i,j,k)/(tauh*tauh)*rho1_SI
      write(25,122) xx,yy,pn
      enddo
      enddo
	close(25)
	endif

c *********************************************
	if (t.eq.tsave) call savedata (z(icc(ng)))
 111  format (2(1x,e15.7)) 
 122	format (3(1x,e15.7)) 
 

      return
      end
c***	*************** print tecplot main ********
      SUBROUTINE printtecp(t,cc) 
      INCLUDE 'input.h' 
      INTEGER i,j,k,t,nn 
      DOUBLE PRECISION cc(nx,ny,nz) 
      DOUBLE PRECISION un,vn,wn,cn,zz,yy,xx 
      CHARACTER (128) FileName 
      WRITE(*,*) 'ntecplot print, t=',t 
      nn=1 
      if (mod(t,ntecplot).ne.0) then 
          return 
      endif 
      write (FileName, '(A,I4.4,A)') 'main', t/ntecplot, '.plt' 
      OPEN(40,FILE=FileName) 
      write(40,*) 'VARIABLES="X","Y","Z","U","V","W","CC"' 
      write(40,*) 'zone i=',nx/nn,' j=',ny/nn-1,' k=',nz/nn 
      do k=1,nz,nn 
        zz = k*h*lengthscale*1000
        do j=2,ny,nn 
          yy = j*h*lengthscale*1000
          do i=1,nx,nn 
            xx = i*h*lengthscale*1000
            if (i.eq.nx) then 
                un=u(i,j,k) 
            else 
                un=.5*(u(i,j,k)+u(i+1,j,k)) 
            endif 

            if (j.eq.ny) then 
                vn=v(i,j,k) 
            else 
                vn=.5*(v(i,j,k)+v(i,j+1,k)) 
            endif 

            if (k.eq.nz) then 
               wn=w(i,j,k) 
            else 
               wn=.5*(w(i,j,k)+w(i,j,k+1)) 
            endif 

            cn=cc(i,j,k) 
            write(40,122) xx,yy,zz,un,vn,wn,cn 
          enddo 
        enddo 
      enddo 
 122  format (7(1x,e15.7)) 
      close(40) 
      RETURN 
      END 
c***

c******************old plot	***********************
      SUBROUTINE printtecp2(t,cc)
      INCLUDE 'input.h'
      INTEGER i,j,k,t,nn
      DOUBLE PRECISION cc(nx,ny,nz)
      DOUBLE PRECISION un,vn,wn,cn,zz,yy,xx
      WRITE(*,*) 'ntecplot print, t=',t
      nn=2
      if(t.eq.0) then
	OPEN(40,FILE='main0.plt')
      write(40,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(40,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.ntecplot) then
	OPEN(41,FILE='main1.plt')
      write(41,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(41,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.2*ntecplot) then
	OPEN(42,FILE='main2.plt')
      write(42,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(42,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.3*ntecplot) then
	OPEN(43,FILE='main3.plt')
      write(43,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(43,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.4*ntecplot) then
	OPEN(44,FILE='main4.plt')
      write(44,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(44,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.5*ntecplot) then
	OPEN(45,FILE='main5.plt')
      write(45,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(45,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.6*ntecplot) then
	OPEN(46,FILE='main6.plt')
      write(46,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(46,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.7*ntecplot) then
	OPEN(47,FILE='main7.plt')
      write(47,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(47,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.8*ntecplot) then
	OPEN(48,FILE='main8.plt')
      write(48,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(48,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.9*ntecplot) then
	OPEN(49,FILE='main9.plt')
      write(49,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(49,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.10*ntecplot) then
	OPEN(50,FILE='main10.plt')
      write(50,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(50,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.11*ntecplot) then
	OPEN(51,FILE='main11.plt')
      write(51,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(51,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.12*ntecplot) then
	OPEN(52,FILE='main12.plt')
      write(52,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(52,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.13*ntecplot) then
	OPEN(53,FILE='main13.plt')
      write(53,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(53,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.14*ntecplot) then
	OPEN(54,FILE='main14.plt')
      write(54,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(54,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.15*ntecplot) then
	OPEN(55,FILE='main15.plt')
      write(55,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(55,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.16*ntecplot) then
	OPEN(56,FILE='main16.plt')
      write(56,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(56,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.17*ntecplot) then
	OPEN(57,FILE='main17.plt')
      write(57,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(57,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.18*ntecplot) then
	OPEN(58,FILE='main18.plt')
      write(58,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(58,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.19*ntecplot) then
	OPEN(59,FILE='main19.plt')
      write(59,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(59,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.20*ntecplot) then
	OPEN(60,FILE='main20.plt')
      write(60,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(60,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.21*ntecplot) then
	OPEN(61,FILE='main21.plt')
      write(61,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(61,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.22*ntecplot) then
	OPEN(62,FILE='main22.plt')
      write(62,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(62,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.23*ntecplot) then
	OPEN(63,FILE='main23.plt')
      write(63,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(63,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.24*ntecplot) then
	OPEN(64,FILE='main24.plt')
      write(64,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(64,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      elseif (t.eq.25*ntecplot) then
	OPEN(65,FILE='main25.plt')
      write(65,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(65,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      else
	OPEN(66,FILE='main26.plt')
      write(66,*) 'VARIABLES="X","Y","Z","U","V","W","CC"'
      write(66,*) 'zone i=',nx/nn,' j=',ny/nn,'k=',nz/nn
      endif
      do k=1,nz,nn
       zz = k*h
        do j=1,ny,nn
         yy = j*h
          do i=1,nx,nn
           xx = i*h

      if (i.eq.nx) then
	un=u(i,j,k)
	else
	un=.5*(u(i,j,k)+u(i+1,j,k))
	endif

      if (j.eq.ny) then
	vn=v(i,j,k)
      else
	vn=.5*(v(i,j,k)+v(i,j+1,k))
	endif

      if (k.eq.nz) then
	wn=w(i,j,k)
      else
	wn=.5*(w(i,j,k)+w(i,j,k+1))
	endif
	
      cn=cc(i,j,k)
      if(t.eq.0) then
      write(40,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.ntecplot) then
      write(41,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.2*ntecplot) then
      write(42,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.3*ntecplot) then
      write(43,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.4*ntecplot) then
      write(44,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.5*ntecplot) then
      write(45,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.6*ntecplot) then
      write(46,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.7*ntecplot) then
      write(47,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.8*ntecplot) then
      write(48,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.9*ntecplot) then
      write(49,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.10*ntecplot) then
      write(50,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.11*ntecplot) then
      write(51,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.12*ntecplot) then
      write(52,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.13*ntecplot) then
      write(53,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.14*ntecplot) then
      write(54,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.15*ntecplot) then
      write(55,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.16*ntecplot) then
      write(56,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.17*ntecplot) then
      write(57,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.18*ntecplot) then
      write(58,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.19*ntecplot) then
      write(59,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.20*ntecplot) then
      write(60,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.21*ntecplot) then
      write(61,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.22*ntecplot) then
      write(62,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.23*ntecplot) then
      write(63,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.24*ntecplot) then
      write(64,122) xx,yy,zz,un,vn,wn,cn
      elseif (t.eq.25*ntecplot) then
      write(65,122) xx,yy,zz,un,vn,wn,cn
      else
      write(66,122) xx,yy,zz,un,vn,wn,cn
      endif
            enddo
         enddo
      enddo
 122  format (7(1x,e15.7))
      if(t.eq.0) then
      close(40)
      elseif (t.eq.ntecplot) then
      close(41)
      elseif (t.eq.2*ntecplot) then
      close(42)
      elseif (t.eq.3*ntecplot) then
      close(43)
      elseif (t.eq.4*ntecplot) then
      close(44)
      elseif (t.eq.5*ntecplot) then
      close(45)
      elseif (t.eq.6*ntecplot) then
      close(46)
      elseif (t.eq.7*ntecplot) then
      close(47)
      elseif (t.eq.8*ntecplot) then
      close(48)
      elseif (t.eq.9*ntecplot) then
      close(49)
      elseif (t.eq.10*ntecplot) then
      close(50)
      elseif (t.eq.11*ntecplot) then
      close(51)
      elseif (t.eq.12*ntecplot) then
      close(52)
      elseif (t.eq.13*ntecplot) then
      close(53)
      elseif (t.eq.14*ntecplot) then
      close(54)
      elseif (t.eq.15*ntecplot) then
      close(55)
      elseif (t.eq.16*ntecplot) then
      close(56)
      elseif (t.eq.17*ntecplot) then
      close(57)
      elseif (t.eq.18*ntecplot) then
      close(58)
      elseif (t.eq.19*ntecplot) then
      close(59)
      elseif (t.eq.20*ntecplot) then
      close(60)
      elseif (t.eq.21*ntecplot) then
      close(61)
      elseif (t.eq.22*ntecplot) then
      close(62)
      elseif (t.eq.23*ntecplot) then
      close(63)
      elseif (t.eq.24*ntecplot) then
      close(64)
      elseif (t.eq.25*ntecplot) then
      close(65)
	endif
      RETURN
      END
c    ********  write unformatted data for restart *********
      SUBROUTINE savedata(cc)
      include 'input.h'
      INTEGER i,j,k
      DOUBLE PRECISION cc(nx,ny,nz)
	open(3,form='unformatted',FILE='saved.dat') 
      write(2) t0
      write (2) nx        ! max points in x-direction
      write (2) ny        ! max points in y-direction
      write (2) nz        ! max points in z-direction
      write (2) tmax      ! maximum number of time steps
      write (2) nfilter   ! filter volume fraction (def.is zero) 
      write (2) cycles    ! max no. of cycles in multigrid method(def.is 60 )
      write (2) OMEGA     ! relaxation coefficient(def.is 1.0)
      write (2) NPRE      ! pre-relaxation coefficient(def.is 1.0)
      write (2) NTOP      ! post-relaxation coefficient(def.is 1.0)
      write (2) NPOST     ! top-relaxation coefficient(def.is 1.0)
      write (2) restart   ! =0 for new data, =1 for continueing the caculation(def.is 0)
      write (2) maxcfl    ! maximum value for CFL(def.is 0.8)
	write (2) iturb     ! =0 laminar,=1 turbulent.
	write (2) convectype ! =1 upwind,=2 central.
	write (2) hswitch   ! =0 off, =1 on

      write (2) xvel1_SI   !initial x-component of velocity of fluid 1
      write (2) yvel1_SI   !initial y-component of velocity of fluid 1
      write (2) zvel1_SI   !initial z-component of velocity of fluid 1
      write (2) xvel2_SI   !initial x-component of velocity of fluid 2
      write (2) yvel2_SI   !initial y-component of velocity of fluid 2
      write (2) zvel2_SI   !initial z-component of velocity of fluid 2
      write (2) rho2_SI    !initial density of fluid 2
      write (2) rho1_SI    !initial density of fluid 1
      write (2) sigma_SI   !surface tension coefficient
      write (2) visc2_SI   !initial viscosity of fluid 2
      write (2) visc1_SI   !initial viscosity of fluid 1
      write (2) grav_SI    !garvity in z-direction
      write (2) grav1_SI   !garvity in x-direction
      write (2) icase      ! module (1-jet, 2=drop, 3=bubble)
      write (2) ibc        ! bourdary condition type (=1, original, =2, standard) 
	write (2) length_SI  ! length of the domain in x-direction
      write (2) radius_si  ! module radius (NOTE: MUST BE CORRECTED!!) 
c      write (2) x0_SI     ! x location of center of jet, drop of bubble 
c	write (2) y0_SI      ! y location of center of jet, drop of bubble
c	write (2) z0_SI      ! z location of center of jet, drop of bubble
      write (2) tprint    ! print data on screen
      write (2) ntecplot  ! print data for ntecplot
      write (2) tsave     ! print unformatted data
      do k=1,nz
         do j=1,ny
            do i=1,nx
	      write(2) u(i,j,k),v(i,j,k),w(i,j,k),cc(i,j,k)
            enddo
         enddo
      enddo
	return
	end

cccccccccccc		end		ccccccccccccc

c  ***** filter volume fraction consider 27 near nodes *******
      
      SUBROUTINE filterc27(cc,cfilt,nx,ny,nz)
      INTEGER nx,ny,nz,i,j,k
      DOUBLE PRECISION cc(nx,ny,nz), cfilt(nx,ny,nz)
      do k=2,nz-1
	do j=2,ny-1
         do i=2,nx-1
	      cfilt(i,j,k) = 0.2729*cc(i,j,k)
     &+ 0.08185*(cc(i-1,j,k)+cc(i+1,j,k)+cc(i,j-1,k)
     &+cc(i,j+1,k)+cc(i,j,k+1)+cc(i,j,k-1))
     
     &+0.0178*(cc(i-1,j-1,k)+cc(i-1,j+1,k)+cc(i+1,j-1,k)
     &+cc(i+1,j+1,k)
     &+cc(i,j+1,k+1)+cc(i,j-1,k+1)+cc(i+1,j,k+1)+cc(i-1,j,k+1)
     &+cc(i,j+1,k-1)+cc(i,j-1,k-1)+cc(i+1,j,k-1)+cc(i-1,j,k-1))
     &+.0028*(cc(i+1,j+1,k+1)+cc(i+1,j+1,k-1)
     &+cc(i+1,j-1,k+1)+cc(i+1,j-1,k-1)+cc(i-1,j+1,k+1)
     &+cc(i-1,j+1,k-1)+cc(i-1,j-1,k+1)+cc(i-1,j-1,k-1))

c	      cfilt(i,j,k) = 1.d0/2.75d0*cc(i,j,k)
c     &+ 1.d0/2.75d0*1.d0/6.d0*(cc(i-1,j,k)+cc(i+1,j,k)+cc(i,j-1,k)
c     &+cc(i,j+1,k)+cc(i,j,k+1)+cc(i,j,k-1))
c     
c     &+1.d0/2.75d0*1.d0/24.d0*(cc(i-1,j-1,k)+cc(i-1,j+1,k)+cc(i+1,j-1,k)
c     &+cc(i+1,j+1,k)
c     &+cc(i,j+1,k+1)+cc(i,j-1,k+1)+cc(i+1,j,k+1)+cc(i-1,j,k+1)
c     &+cc(i,j+1,k-1)+cc(i,j-1,k-1)+cc(i+1,j,k-1)+cc(i-1,j,k-1))
c     &+1.d0/2.75d0*1.d0/32.d0*(cc(i+1,j+1,k+1)+cc(i+1,j+1,k-1)
c     &+cc(i+1,j-1,k+1)+cc(i+1,j-1,k-1)+cc(i-1,j+1,k+1)
c     &+cc(i-1,j+1,k-1)+cc(i-1,j-1,k+1)+cc(i-1,j-1,k-1))

!***************************************************************
      if (cfilt(i,j,k).gt.1.d0) cfilt(i,j,k)=1.d0
!*************************************************************** 
         enddo
      enddo
	enddo

      return
      end

c  ***** filter volume fraction consider 18 near nodes *******
      
      SUBROUTINE filterc18(cc,cfilt,nx,ny,nz)
      INTEGER nx,ny,nz,i,j,k
      DOUBLE PRECISION cc(nx,ny,nz), cfilt(nx,ny,nz)
      do k=2,nz-1
	do j=2,ny-1
         do i=2,nx-1
	      cfilt(i,j,k) = 1.d0/2.5d0*cc(i,j,k)
     *+ 1.d0/2.5d0*1.d0/6.d0*(cc(i-1,j,k)+cc(i+1,j,k)+cc(i,j-1,k)
     *+cc(i,j+1,k)+cc(i,j,k+1)+cc(i,j,k-1))
     
     *+1.d0/2.5d0*1.d0/24.d0*(cc(i-1,j-1,k)+cc(i-1,j+1,k)+cc(i+1,j-1,k)
     *+cc(i+1,j+1,k)
     *+cc(i,j+1,k+1)+cc(i,j-1,k+1)+cc(i+1,j,k+1)+cc(i-1,j,k+1)
     *+cc(i,j+1,k-1)+cc(i,j-1,k-1)+cc(i+1,j,k-1)+cc(i-1,j,k-1))
!***************************************************************
      if (cfilt(i,j,k).gt.1.d0) cfilt(i,j,k)=1.d0
!*************************************************************** 
         enddo
      enddo
	enddo

      return
      end
c  ***** filter volume fraction consider 6 near nodes *******
      
      SUBROUTINE filterc6(cc,cfilt,nx,ny,nz)
      INTEGER nx,ny,nz,i,j,k
      DOUBLE PRECISION cc(nx,ny,nz), cfilt(nx,ny,nz)
      do k=2,nz-1
	do j=2,ny-1
         do i=2,nx-1
	      cfilt(i,j,k) = 0.5d0*cc(i,j,k)
     *+0.5*1.d0/6.d0*(cc(i-1,j,k)+cc(i+1,j,k)+cc(i,j-1,k)+cc(i,j+1,k)
     *+cc(i,j,k+1)+cc(i,j,k-1))
!***************************************************************
      if (cfilt(i,j,k).gt.1.d0) cfilt(i,j,k)=1.d0
!*************************************************************** 
         enddo
      enddo
	enddo
      return
      end
c  ***** filter volume fraction consider 8 near nodes in x-y plane ******
      
      SUBROUTINE filterc8(cc,cfilt,nx,ny,nz)
      INTEGER nx,ny,nz,i,j,k
      DOUBLE PRECISION cc(nx,ny,nz), cfilt(nx,ny,nz)
      do k=1,nz-1
	do j=2,ny-1
         do i=2,nx-1
	      cfilt(i,j,k) = 0.25d0*cc(i,j,k)
     *+0.125d0*(cc(i-1,j,k) + cc(i+1,j,k)+cc(i,j-1,k)+cc(i,j+1,k))
     *+0.0625d0*(cc(i-1,j-1,k)+cc(i+1,j-1,k)+cc(i-1,j+1,k)+cc(i+1,j+1,k)
     *)
!***************************************************************
      if (cfilt(i,j,k).gt.1.d0) cfilt(i,j,k)=1.d0
!*************************************************************** 
         enddo
      enddo
	enddo
      return
      end


c    ********  creat initial conditions for smooth drop  *********
      SUBROUTINE smoothdrop(cc,dist)
      include 'input.h'
      INTEGER i,j,k,nh
      DOUBLE PRECISION cc(nx,ny,nz),dist(nx,ny,nz),ccf(nx,ny,nz)
      DOUBLE PRECISION x0,xx,yy,zz,rr,umax
      x0=(nx-2)*h/2.0

c***	if z0_SI equal 1000 z0 will be in center of channel
	if (z0_si.ne.1000) z0= z0_SI/length_SI
c***	if y0_SI equal 1000 y0 will be in center of channel
	if (y0_si.ne.1000) y0= y0_SI/length_SI
c	y0=y0+h !/2.d0

	do k=1,nz
         zz = k*h
         do j=1,ny
            yy = j*h
            do i=1,nx
            xx = i*h
            rr= DSQRT((xx-x0)*(xx-x0)+(yy-y0)*(yy-y0)+(zz-z0)*(zz-z0))
c     put pressure field in dist
                dist(i,j,k) = 0.d0
               if(rr.le.radius-h/2) then
               cc(i,j,k)=1.d0
           else if (rr.ge.radius+h/2) then
               cc(i,j,k)=0.
		else
			cc(i,j,k)=(radius+h/2-rr)/h			
            endif
         enddo
	  enddo
      enddo
	do k=1,nz
		do i=1,nx
			do j=1,ny
            
			if (cc(i,j,k).eq.0) then 
c			u(i,j,k)=xvel2*(-((j*1.d0-1.5d0)/(ny*1.d0)-1.5d0)**2+1.d0)
c			if (j.eq.1) u(i,j,k)=0.0
			u(i,j,k)=xvel2
              v(i,j,k)=yvel2
              w(i,j,k)=zvel2
			elseif (cc(i,j,k).gt.0) then 
               u(i,j,k)=xvel1
               v(i,j,k)=yvel1
		   	 w(i,j,k)=zvel1
			endif
			if (j.gt.1.and.cc(i,j-1,k).gt.0.and.cc(i,j,k).eq.0) then
			nh=j
			endif
c	if(cc(i,1,k).gt.0 .and.cc(i,j,k).eq.0)then
c	umax=2.d0/3.d0*xvel2/(2.d0/3.d0+dfloat(ny-nh)/dfloat(ny)
c     &-1.0-1.d0/3.d0*((ny-nh)/ny)**3)
c	u(i,j,k)=umax*(-(dfloat(j-nh)/dfloat(ny-nh)-1.d0)**2+1.d0)	
c            endif
         enddo
	  enddo
      enddo
      return
      end

c	*********** plot of u,v,p,c in one line of domain	**********
	subroutine uplot(t,cc,p)
	INCLUDE 'input.h' 
      INTEGER i,j,k,t,nn 
      DOUBLE PRECISION cc(nx,ny,nz),p(nx,ny,nz)

      DOUBLE PRECISION un,vn,wn,cn,zz,yy,xx 
      CHARACTER (128) FileName 
      WRITE(*,*) 'ntecplot print, t=',t 
      nn=1 
      if (mod(t,ntecplot).ne.0) then 
          return 
      endif 
      write (FileName, '(A,I4.4,A)') 'u', t/10, '.plt' 
      OPEN(40,FILE=FileName) 
      write(40,*) 'VARIABLES="Z","U","V","p","CC"' 
      write(40,*) 'zone k=',ny-1 
      k=(nz-2)/2
        zz = k*h
	  i=(nx-2)/2 
        xx = i*h
	  do j=1,ny-1 
          yy = j*h                          
                un=0.5*(u(i,j,k)+u(i+1,j,k))              
                vn=0.5*(v(i,j,k)+v(i,j+1,k)) 
                cn=cc(i,j,k) 
            write(40,122) yy,un,vn,p(i,j,k),cn 
         
        enddo 
      
 122  format (7(1x,e15.7)) 
      close(40) 

      RETURN 
      END 

c	*****************************************
      SUBROUTINE smooth2(cc,dist)

	include 'input.h'
      INTEGER i,j,k,l,m,n,index
      DOUBLE PRECISION cc(nx,ny,nz),dist(nx,ny,nz),ccf(nx,ny,nz)
      DOUBLE PRECISION x0,xx(nx),yy(ny),zz(nz),rr,ctotal
      x0=(nx-2)*h/4.0

c***	if z0_SI equal 1000 z0 will be in center of channel
	if (z0_si.ne.1000) z0= z0_SI/length_SI
c***	if y0_SI equal 1000 y0 will be in center of channel
	if (y0_si.ne.1000) y0= y0_SI/length_SI
	y0=y0+h/2.d0
      
	do i=2,nx-1
         xx(i) = (i-1.5)*h
         do j=2,ny-1
            yy(j) = (j-1.5)*h
			do k=2,nz-1
				zz(k) = (k-1.5)*h
		   if (dmin1((xx(i-1)-x0)**2,(xx(i)-x0)**2)
     1            +dmin1((yy(j-1)-y0)**2,(yy(j)-y0)**2)
	2			+dmin1((zz(k-1)-z0)**2,(zz(k)-z0)**2)
     3            .lt.radius**2) go to 40
			cc(i,j,k) =0.d0
	     		go to 60
c
c................ if cell isn't full
c
   40 if (dmax1((xx(i-1)-x0)**2,(xx(i)-x0)**2)
     1            +dmax1((yy(j-1)-y0)**2,(yy(j)-y0)**2)
	2			+dmax1((zz(k-1)-z0)**2,(zz(k)-z0)**2)
     3            .gt.radius**2) go to 50
			cc(i,j,k) = 1.d0
			go to 60
c
c................ cell neither empty nor full
c
   50          cc(i,j,k)=0.d0
	       index=30
		   do l=1,index
		   			do  m=1,index
						do n=1,index !-index/2.0,index/2.0
				if ((xx(i-1)+(dfloat(l)-.5d0)/dfloat(index)*h-x0)**2
     2               +(yy(j-1)+(dfloat(m)-.5d0)/dfloat(index)*h-y0)**2
     3               +(zz(k-1)+(dfloat(n)-.5d0)/dfloat(index)*h-z0)**2 
     4				.lt.radius**2) then
                       cc(i,j,k)=cc(i,j,k)+1.d0/dfloat((index)**3)
c                       u(i,j,k)=u(i,j,k)+xvelliq/dfloat(index**3)
c                       v(i,j,k)=v(i,j,k)+yvelliq/dfloat(index**3)
						
                  endif
             enddo
              enddo
			enddo
   60      CONTINUE
          
            enddo
         enddo
	enddo

      call bc_c(cc)
      do i=1,nx
         do j=1,ny
		 do k=1,nz
			u(i,j,k)= 0.d0 !xvel1*cc(i,j,k)+xvel2*(1-cc(i,j,k))
			v(i,j,k)=yvel1*cc(i,j,k)+yvel2*(1-cc(i,j,k))
			w(i,j,k)=zvel1*cc(i,j,k)+zvel2*(1-cc(i,j,k))

              
        enddo
      enddo
	enddo
      return
      end
c *****************************************
	subroutine plot(cc)
	include 'input.h'
	integer i,j,k
	double precision cc(nx,ny,nz)
	open(123,file='w.plt')
	write(123,*) 'VARIABLES="x","y","cc"'
	write(123,*) 'zone i=',nx-2,'j=',ny-2
	k=(nz-2)/2
      do j= 2,ny-1
      do i= 2,nx-1

      write(123,122) h*i,h*j,cc(i,j,k)
      enddo
      enddo
	close(123)
 122	format (3(1x,e15.7))
	return
	end 
c	*************** velocity center **********
	subroutine velocityc(cc)
	include 'input.h'
	integer i,j,k
	double precision cc(nx,ny,nz),ii,jj,sum,uuu,uud,uu
	double precision vvu,vvd,vv,vvv
	sum=0.0
	ii=0.0
	jj=0.0
	k=(nz-2)/2
	do i=2,nx-1
      do j= 2,ny-1
	ii=i*cc(i,j,k)+ii
	jj=j*cc(i,j,k)+jj
      enddo
	enddo
	
	do i=2,nx-1
      do j= 2,ny-1
	sum=cc(i,j,k)+sum
      enddo
	enddo
	ii=ii/sum
	jj=jj/sum
	i=int(ii)
	j=int(jj)
	uuu=u(i,j+1,k)+(u(i+1,j+1,k)-u(i,j+1,k))*(ii-i*1.0)
	uud=u(i,j,k)+(u(i+1,j,k)-u(i,j,k))*(ii-i*1.0)
	uu=uud+(uuu-uud)*(jj-j*1.0)

	vvu=v(i,j+1,k)+(v(i+1,j+1,k)-v(i,j+1,k))*(ii-i*1.0)
	vvd=v(i,j,k)+(v(i+1,j,k)-v(i,j,k))*(ii-i*1.0)
	vv=vvd+(vvu-vvd)*(jj-j*1.0)
	vvv=dsqrt(vv**2+uu**2)
	write (28,111) time*timescale*1000,vvv*lengthscale/timescale
 111	format (2(1x,e15.7)) 
	return
	end