c     ================================================================
c     CARDIAC ELECTROPHYSIOLOGY SIMULATION - CORRECTED VERSION
c     ================================================================
c     
c     This code simulates 2D cardiac tissue with calcium dynamics
c     and triggered calcium waves. Models cells with sparse t-tubules
c     where calcium influx is confined to cell periphery.
c     
c     Compilation: ifx ead1.f atissue-cax7.f atissue-vx7.f
c     ================================================================

c     double precision code
c     ifx ead1.f atissue-cax7.f atissue-vx7.f

      implicit double precision (a-h,o-z)
      
c     ================================================================
c     SIMULATION PARAMETERS
c     ================================================================
      parameter (Lx=5,Ly=5)        ! Spatial grid dimensions (5x5 cells)
      parameter(nstim=10)          ! Number of stimulus beats

c     ================================================================
c     VARIABLE DECLARATIONS
c     ================================================================
      
      ! Pacing protocol arrays
      double precision rbcl(nstim)
    
      ! Character arrays for file names (legacy)
      character*1 cam(0:9)
      character*15 filenm1(999),filenm2(999), filenm3(999), filenm4(999)
      character*15 filenm5(999),filenm6(999)

      ! Analysis arrays
      double precision tstim(nstim),ampct(1000),apmax(1000)
      double precision ampvt(1000)

      ! ============================================================
      ! VOLTAGE AND MEMBRANE POTENTIAL ARRAYS
      ! ============================================================
      double precision v(0:Lx+1,0:Ly+1)    ! membrane voltage (mV)
      double precision vnew(0:Lx+1,0:Ly+1) ! updated voltage after diffusion
      double precision dv(0:Lx,0:Ly)       ! voltage time derivative
      double precision vold(0:Lx,0:Ly)     ! voltage from previous time step
      
      ! ============================================================
      ! ION CHANNEL GATING VARIABLES
      ! ============================================================
      double precision xm(Lx,Ly)    ! INa activation gate
      double precision xh(Lx,Ly)    ! INa fast inactivation gate
      double precision xj(Lx,Ly)    ! INa slow inactivation gate
      double precision xr(Lx,Ly)    ! IKr activation gate
      double precision xs1(Lx,Ly)   ! IKs activation gate
      double precision qks(Lx,Ly)   ! IKs calcium-dependent gate
      double precision xkur(Lx,Ly)  ! IKur activation gate
      double precision ykur(Lx,Ly)  ! IKur inactivation gate
      
      ! Transient outward current gates
      double precision xtof(Lx,Ly)  ! Ito fast activation
      double precision ytof(Lx,Ly)  ! Ito fast inactivation
      double precision xtos(Lx,Ly)  ! Ito slow activation
      double precision ytos(Lx,Ly)  ! Ito slow inactivation

      ! ============================================================
      ! CALCIUM CONCENTRATION VARIABLES
      ! ============================================================
      double precision cb(Lx,Ly)    ! boundary cytosolic Ca (μM)
      double precision ci(Lx,Ly)    ! interior cytosolic Ca (μM)
      double precision csrb(Lx,Ly)  ! boundary SR Ca (μM)
      double precision csri(Lx,Ly)  ! interior SR Ca (μM)
      double precision cnsr(Lx,Ly)  ! network SR Ca (μM)
      double precision cit(Lx,Ly)   ! total interior Ca (free + buffered)
      double precision cbt(Lx,Ly)   ! total boundary Ca (free + buffered)

      ! ============================================================
      ! L-TYPE CALCIUM CHANNEL MARKOV STATES
      ! ============================================================
      double precision po(Lx,Ly)    ! open probability (boundary-facing)
      double precision c1(Lx,Ly)    ! closed state 1 (boundary)
      double precision c2(Lx,Ly)    ! closed state 2 (boundary)
      double precision xi1(Lx,Ly)   ! inactivated state 1 (boundary)
      double precision xi2(Lx,Ly)   ! inactivated state 2 (boundary)
      
      ! Spark-facing L-type channel states
      double precision c1s(Lx,Ly)   ! closed state 1 (spark-facing)
      double precision c2s(Lx,Ly)   ! closed state 2 (spark-facing)
      double precision xi1s(Lx,Ly)  ! inactivated state 1 (spark-facing)
      double precision xi2s(Lx,Ly)  ! inactivated state 2 (spark-facing)

      ! Additional calcium current variables
      double precision xicaqz(Lx,Ly) ! L-type Ca current density

      ! ============================================================
      ! CALCIUM SPARK VARIABLES
      ! ============================================================
      double precision pi(Lx,Ly)    ! fraction of interior clusters sparking
      double precision pb(Lx,Ly)    ! fraction of boundary clusters sparking
      double precision pox(Lx,Ly)   ! total LCC open probability
      double precision pos(Lx,Ly)   ! open probability (spark-facing)
      double precision ra(Lx,Ly)    ! spark recruitment variable

      ! Spark cluster counters
      integer nsb(Lx,Ly)   ! number of boundary clusters with sparks
      integer nsi(Lx,Ly)   ! number of interior clusters with sparks
	
      ! ============================================================
      ! ANALYSIS AND OUTPUT ARRAYS
      ! ============================================================
      double precision ctmax(1000)                   ! calcium maxima
      double precision camp(0:Lx,0:Ly)               ! calcium amplitude
      double precision cmax(0:Lx,0:Ly,0:1000)        ! maximum Ca each beat
      double precision uxx(0:Lx,0:Ly)                ! analysis array
      double precision xuu(0:Lx,0:Ly)                ! analysis array

      double precision apd(0:Lx,0:Ly,0:1500)  ! action potential duration
      double precision yapd(0:Lx,0:Ly,0:2)    ! voltage history for APD

      ! ============================================================
      ! HETEROGENEITY ARRAYS
      ! ============================================================
      double precision gicaz(0:Lx,0:Ly)    ! L-type Ca conductance distribution
      double precision gnacaz(0:Lx,0:Ly)   ! Na-Ca exchange scaling distribution
      double precision pbxz(0:Lx,0:Ly)     ! spark threshold distribution

      ! Current arrays for analysis
      double precision xicaz(0:Lx,0:Ly)    ! L-type Ca current
      double precision xinacaz(0:Lx,0:Ly)  ! Na-Ca exchange current
      double precision xik1z(0:Lx,0:Ly)    ! IK1 current

      ! Random number generation
      integer, allocatable :: seed_array(:)
      integer :: i, seed_size, n

      ! Temporary array
      real yzz(3)

c     CRITICAL: Define nbt as a parameter to avoid uninitialized use
      integer, parameter :: nbt = 3000  ! total boundary clusters

c     ================================================================
c     MODEL SELECTION AND FILE SETUP
c     ================================================================
   
      mk = 2    ! Model selection flag
c     mk = 1:   Normal case
c     mk = 2:   Enhanced Ca recruitment due to RyR2 leak

      if(mk.eq.1) then
         Dfu = 0.001d0    ! Diffusion coefficient (cm²/ms)
         open(unit=10,file='v1x.dat',status='unknown')
      endif 

      if(mk.eq.2) then 
         Dfu = 0.001d0    ! Diffusion coefficient (cm²/ms)
         open(unit=10,file='v2x.dat',status='unknown')
      endif 

c     ================================================================
c     PHYSIOLOGICAL CONSTANTS
c     ================================================================

      ! Ion concentrations (mM)
      xnao = 136.0d0      ! extracellular sodium
      xki = 140.0d0       ! intracellular potassium
      xko = 5.40d0        ! extracellular potassium
      cao = 1.8d0         ! extracellular calcium
	
      ! Physical constants
      temp = 308.0d0      ! temperature (K) - 35°C
      xxr = 8.314d0       ! gas constant (J/mol/K)
      xf = 96.485d0       ! Faraday's constant (C/mol)
      frt = xf/(xxr*temp) ! F/RT for Nernst calculations

c     ================================================================
c     RANDOM NUMBER INITIALIZATION
c     ================================================================
      call random_seed(size=n)
      allocate(seed_array(n))
      seed_array = 3422   ! Fixed seed for reproducibility
      call random_seed(put=seed_array)

c     ================================================================
c     MODEL PARAMETERS
c     ================================================================
      rbcl1 = 500.0       ! S1 pacing cycle length (ms)
	
      gicai = 1.4         ! L-type Ca current scaling factor
      gtos = 0.05         ! Ito slow conductance
      gtof = 0.15         ! Ito fast conductance
      gnacai = 2.80       ! Na-Ca exchange scaling factor
      zxr = 0.08          ! Calcium dependence parameter
      cxinit = 1000.0     ! Initial SR calcium concentration (μM)
        
c     ================================================================
c     NUMERICAL INTEGRATION PARAMETERS
c     ================================================================
      dt = 0.1d0          ! Basic time step (ms)
      mstp0 = 10          ! Maximum sub-steps for adaptive integration
  
      ! Calculate intracellular Na based on cycle length (empirical relation)
      xmx = -2.0/250.0
      xnai = xmx*rbcl1 + 16.0

c     ================================================================
c     SPATIAL DISCRETIZATION
c     ================================================================
      dx = 0.015d0        ! Spatial step in x-direction (cm)
      dy = 0.015d0        ! Spatial step in y-direction (cm)

      ! Diffusion stability parameters
      slmbdax = Dfu*dt/(4.0d0*dx*dx)  ! x-direction diffusion parameter
      slmbday = Dfu*dt/(4.0d0*dy*dy)  ! y-direction diffusion parameter

c     ================================================================
c     PACING PROTOCOL SETUP
c     ================================================================
      do i = 1, nstim
         rbcl(i) = rbcl1   ! All beats at basic cycle length
      enddo

c     ================================================================
c     INITIAL CONDITIONS
c     ================================================================
      do ix = 1, Lx
         do iy = 1, Ly
            
c           ============================================
c           CALCIUM INITIAL CONDITIONS
c           ============================================
            cb(ix,iy) = 0.1d0      ! boundary cytosolic Ca (μM)
            ci(ix,iy) = 0.1d0      ! interior cytosolic Ca (μM)

            csrb(ix,iy) = cxinit   ! boundary SR Ca (μM)
            csri(ix,iy) = cxinit   ! interior SR Ca (μM)
            cnsr(ix,iy) = cxinit   ! network SR Ca (μM)
	 
c           ============================================
c           L-TYPE CA CHANNEL INITIAL STATES
c           ============================================
            po(ix,iy) = 0.0d0      ! open probability (boundary-facing)
            c1(ix,iy) = 0.0d0      ! closed state 1
            c2(ix,iy) = 1.0d0      ! closed state 2 (all channels start here)
            xi1(ix,iy) = 0.0d0     ! inactivated state 1
            xi2(ix,iy) = 0.0d0     ! inactivated state 2
	
            ra(ix,iy) = 0.0d0      ! spark recruitment variable
	
            ! Spark-facing L-type channel states
            pos(ix,iy) = 0.0d0     ! open probability (spark-facing)
            c1s(ix,iy) = 0.0d0     ! closed state 1
            c2s(ix,iy) = 0.0d0     ! closed state 2
            xi1s(ix,iy) = 0.0d0    ! inactivated state 1
            xi2s(ix,iy) = 0.0d0    ! inactivated state 2

c           ============================================
c           CONVERT FREE TO TOTAL CALCIUM
c           ============================================
            call total(ci(ix,iy),cit(ix,iy))  ! interior
            call total(cb(ix,iy),cbt(ix,iy))  ! boundary

c           ============================================
c           SPARK CLUSTER INITIALIZATION
c           ============================================
            nsb(ix,iy) = 5         ! initial boundary clusters with sparks

c           ============================================
c           VOLTAGE AND ION CHANNEL GATES
c           ============================================
            v(ix,iy) = -90.0d0     ! resting membrane potential (mV)
            
            ! Sodium channel gates
            xm(ix,iy) = 0.001d0    ! activation (mostly closed at rest)
            xh(ix,iy) = 1.0d0      ! fast inactivation (not inactivated)
            xj(ix,iy) = 1.0d0      ! slow inactivation (not inactivated)
            
            ! Potassium channel gates
            xr(ix,iy) = 0.0d0      ! IKr activation (closed)
            xs1(ix,iy) = 0.9d0     ! IKs activation
            qks(ix,iy) = 0.2d0     ! IKs calcium dependence
	  
            ! Transient outward current gates
            xtos(ix,iy) = 0.01d0   ! Ito slow activation
            ytos(ix,iy) = 0.9d0    ! Ito slow inactivation
            xtof(ix,iy) = 0.02d0   ! Ito fast activation
            ytof(ix,iy) = 0.8d0    ! Ito fast inactivation

            ! Initialize voltage history
            vold(ix,iy) = v(ix,iy)

c           ============================================
c           HETEROGENEITY INITIALIZATION
c           ============================================
            gicaz(ix,iy) = gicai   ! L-type Ca conductance
            gnacaz(ix,iy) = gnacai ! Na-Ca exchange strength
	
        enddo
      enddo

c     ================================================================
c     INITIALIZE ALL ARRAYS PROPERLY
c     ================================================================
      ! Initialize analysis variables
      zci1 = 0.0    ! calcium analysis variables
      zci2 = 0.0
      zci3 = 0.0

      ap1 = 0.0     ! action potential analysis variables
      ap2 = 0.0
      ap3 = 0.0

      ! Initialize APD detection arrays properly
      do ix = 0, Lx
         do iy = 0, Ly
            do i = 0, 2
               yapd(ix,iy,i) = 0.0d0
            enddo
         enddo
      enddo

c     ================================================================
c     CELLULAR VOLUME RATIOS
c     ================================================================
      vi = 0.30     ! interior cytoplasm volume (normalized)
      vb = 0.30     ! boundary cytoplasm volume
      vbi = vb/vi   ! boundary to interior volume ratio

      vq = 30.0     ! SR volume scaling factor
      visr = vq     ! interior SR volume factor
      vsrin = (1.0/visr)*vi   ! interior SR to cytoplasm ratio

      vbsr = vq     ! boundary SR volume factor
      vsrbound = (1.0/vbsr)*vb  ! boundary SR to cytoplasm ratio

      vbisr = (vsrbound/vsrin)  ! boundary to interior SR ratio

      vnsr = vq     ! network SR volume factor
      vnsrin = (1.0/vnsr)*vi    ! network SR to cytoplasm ratio

c     ================================================================
c     MAIN INTEGRATION LOOP
c     ================================================================
      kku = 1       ! counter variable
      t = 0.0       ! total elapsed time (ms)

      ! Loop over stimulus beats
      do iz = 1, nstim
         nstep = int(rbcl(iz)/dt)  ! time steps in this beat

         umax = 0.0   ! analysis variables
         vmax = 0.0
	
         ! Initialize calcium maximum tracking for this beat
         do ix = 1, Lx
            do iy = 1, Ly
               cmax(ix,iy,iz) = 0.0
            enddo 
         enddo

         ! Time integration within each beat
         do ncount = 0, nstep
            time = dfloat(ncount)*dt  ! time within current beat
       
            ! Loop over all spatial grid points
            do iy = 1, Ly
               do ix = 1, Lx

c                 ========================================
c                 APPLY TISSUE HETEROGENEITY
c                 ========================================
                  gica = gicaz(ix,iy)   ! local L-type Ca conductance
                  gnaca = gnacaz(ix,iy) ! local Na-Ca exchange strength

c                 ========================================
c                 CALCIUM BUFFERING CALCULATIONS
c                 ========================================
                  ! Convert total calcium to free calcium concentrations
                  call xfree(cit(ix,iy),ci(ix,iy))  ! interior
                  call xfree(cbt(ix,iy),cb(ix,iy))  ! boundary

c                 ========================================
c                 CALCIUM SPARK DYNAMICS
c                 ========================================
                  ! Calculate fraction of boundary clusters with sparks
                  pb(ix,iy) = dfloat(nsb(ix,iy))/dfloat(nbt)

c                 ========================================
c                 SR CALCIUM UPTAKE (SERCA PUMP)
c                 ========================================
                  vupb = 0.1    ! boundary uptake rate
                  vupi = 0.1    ! interior uptake rate
                  call uptake(cb(ix,iy),vupb,xupb)  ! boundary uptake
                  call uptake(ci(ix,iy),vupi,xupi)  ! interior uptake

c                 ========================================
c                 Na-Ca EXCHANGE CURRENT
c                 ========================================
                  ! Calculate NCX at boundary and interior
                  call inaca(v(ix,iy),frt,xnai,xnao,cao,cb(ix,iy),
     +                      xinacaq1)
                  call inaca(v(ix,iy),frt,xnai,xnao,cao,ci(ix,iy),
     +                      xinacaq2)

                  ru = 1.0  ! weighting factor (1.0 = all boundary)
                  xinacaq = ru*xinacaq1 + (1.0-ru)*xinacaq2
                  xinacaq = gnaca*xinacaq  ! apply scaling

c                 ========================================
c                 L-TYPE CALCIUM CURRENT
c                 ========================================
                  pox(ix,iy) = po(ix,iy) + pos(ix,iy)  ! total open probability
                  call ica(v(ix,iy),frt,cao,cb(ix,iy),pox(ix,iy),
     +                   rca,xicaq)
                  xicaq = gica*130.0*xicaq  ! scale and convert units

c                 ========================================
c                 CALCIUM SPARK RECRUITMENT
c                 ========================================
                  ! Model-dependent spark parameters
                  if(mk.eq.1) then
                     ! Normal case
                     qq = 0.5d0        ! spark sensitivity
                     ab = 35.0*qq      ! base spark rate
                     bts = 1.0/30.0    ! spark termination rate
                     csrx = 600.0d0    ! SR Ca threshold
                     phisr = 1.0/(1.0+(csrx/csrb(ix,iy))**10)  ! SR load dependence
                  endif
	
                  if(mk.eq.2) then
                     ! Enhanced recruitment case (RyR2 leak)
                     qq = 0.6d0        ! increased sensitivity
                     ab = 35.0*qq      ! base spark rate
                     bts = 1.0/30.0    ! spark termination rate
                     csrx = 100.0d0    ! lower SR Ca threshold
                     phisr = 1.0/(1.0+(csrx/csrb(ix,iy))**4)   ! different Hill coefficient
                  endif 

                  ! LCC-triggered spark rate
                  alphab = ab*dabs(rca)*po(ix,iy)*phisr
                  bts = 1.0/30.0  ! spark termination rate

c                 ========================================
c                 L-TYPE Ca CHANNEL MARKOV MODEL UPDATE
c                 ========================================
                  call markov(dt,v(ix,iy),cb(ix,iy),c1(ix,iy),
     +                  c2(ix,iy),xi1(ix,iy),xi2(ix,iy),po(ix,iy),
     +                  c1s(ix,iy),c2s(ix,iy),
     +                  xi1s(ix,iy), xi2s(ix,iy),pos(ix,iy),
     +                  alphab,bts,zxr)

c                 ========================================
c                 RyR CALCIUM RELEASE
c                 ========================================
                  ! Boundary RyR release
                  gsrb = (0.01/1.5)*1.0
                  xryrb = gsrb*csrb(ix,iy)*pb(ix,iy)
	
                  ! Interior RyR release (currently disabled)
                  xryri = 0.0

c                 ========================================
c                 CALCIUM BALANCE EQUATIONS
c                 ========================================
                  ! Volume parameters (redefined for this section)
                  vi = 0.50     ! interior volume
                  vb = 1.0      ! boundary volume
                  vbi = vb/vi   ! volume ratio
                  vbisr = vbi   ! SR volume ratio

                  vq = 30.0     ! SR scaling
                  visr = 30.0   ! interior SR scaling
                  vbsr = vq     ! boundary SR scaling
                  visr = vq     ! interior SR scaling

                  ! Diffusion time constants
                  tau1 = 5.0    ! cytoplasm diffusion
                  tau2 = 5.0    ! SR diffusion

                  ! Diffusive currents
                  dfbi = (cb(ix,iy)-ci(ix,iy))/tau1       ! cyto diffusion
                  dfbisr = (csrb(ix,iy)-csri(ix,iy))/tau2 ! SR diffusion

                  ! Net sarcolemmal calcium flux
                  xsarc = -xicaq + xinacaq

                  ! Time derivatives for calcium concentrations
                  dcbt = xryrb - xupb + xsarc - dfbi              ! boundary cyto
                  dcsrb = vbsr*(-xryrb+xupb) - dfbisr            ! boundary SR
                  dcit = xryri - xupi + vbi*dfbi                 ! interior cyto
                  dcsri = visr*(-xryri+xupi) + vbisr*dfbisr      ! interior SR
        
                  ! Update calcium concentrations
                  cbt(ix,iy) = cbt(ix,iy) + dcbt*dt
                  cit(ix,iy) = cit(ix,iy) + dcit*dt
                  csrb(ix,iy) = csrb(ix,iy) + dcsrb*dt	
                  csri(ix,iy) = csri(ix,iy) + dcsri*dt

c                 ========================================
c                 STOCHASTIC SPARK EVOLUTION
c                 ========================================
                  nsbx = nsb(ix,iy)
                  call binevol(nbt,nsbx,alphab,bts,dt,ndeltapx,ndeltamx)

                  ! Ensure spark numbers stay within bounds
                  if(ndeltamx.gt.nsbx.or.ndeltapx.gt.nbt) then
                     nsb(ix,iy) = 0
                  else
                     nsb(ix,iy) = nsb(ix,iy) + ndeltapx - ndeltamx
                  endif 

c                 ========================================
c                 MEMBRANE VOLTAGE DYNAMICS
c                 ========================================
                  ! Convert calcium fluxes to membrane currents
                  wca = 12.0d0               ! conversion factor
                  xinaca = wca*xinacaq       ! Na-Ca exchange current
                  xica = 2.0d0*wca*xicaq     ! L-type Ca current

c                 ========================================
c                 ADAPTIVE TIME STEPPING
c                 ========================================
                  adq = dabs(dv(ix,iy))
                
                  if(adq.gt.25.0d0) then
                     mstp = 10  ! fine time step for fast voltage changes
                  else
                     mstp = 1   ! normal time step
                  endif 
	              
                  hode = dt/dfloat(mstp)  ! actual integration time step

                  ! Sub-stepping loop for voltage integration
                  do iii = 1, mstp  

c                    ================================
c                    ION CHANNEL CURRENT CALCULATIONS
c                    ================================
                     ! Fast sodium current
                     call ina(hode,v(ix,iy),frt,xh(ix,iy),xj(ix,iy),
     &                       xm(ix,iy),xnai,xnao,xina)

                     ! Rapid delayed rectifier K current
                     call ikr(hode,v(ix,iy),frt,xko,xki,
     &                       xr(ix,iy),xikr)

                     ! Slow delayed rectifier K current
                     call iks(hode,v(ix,iy),frt,cb(ix,iy),xnao,xnai,
     &                       xko,xki,xs1(ix,iy),qks(ix,iy),xiks)

                     ! Inward rectifier K current
                     call ik1(hode,v(ix,iy),frt,xki,xko,xik1)

                     ! Transient outward K current
                     call ito(hode,v(ix,iy),frt,xki,xko,xtof(ix,iy),
     &                       ytof(ix,iy),xtos(ix,iy),ytos(ix,iy),
     &                       xito,gtof,gtos)

                     ! Na-K pump current
                     call inak(v(ix,iy),frt,xko,xnao,xnai,xinak)

c                    ================================
c                    STIMULUS CURRENT APPLICATION
c                    ================================
                     if(time.lt.1.0) then
                        stim = 80.0d0  ! stimulus amplitude (μA/μF)
                     else
                        stim = 0.0     ! no stimulus
                     endif       

c                    ================================
c                    TOTAL MEMBRANE CURRENT
c                    ================================
                     dvh = -(xina+xik1+xikr+xiks+xito+xinaca+xica+xinak) 
     +                    + stim

                     ! Update membrane voltage
                     v(ix,iy) = v(ix,iy) + dvh*hode

                  enddo  ! end sub-stepping loop

               enddo  ! end ix loop
            enddo     ! end iy loop

c           ============================================
c           SPATIAL DIFFUSION OF VOLTAGE
c           ============================================
            ! Apply diffusion operator twice for stability
            call Euler_Forward(Lx,Ly,v,vnew,slmbdax,slmbday)
            call Euler_Forward(Lx,Ly,v,vnew,slmbdax,slmbday)

c           ============================================
c           ANALYSIS AND DATA COLLECTION
c           ============================================
            ctoti = 0.0  ! total calcium (for analysis)

            do iy = 1, Ly
               do ix = 1, Lx
                  ! Calculate voltage time derivative
                  dv(ix,iy) = (vnew(ix,iy)-vold(ix,iy))/dt
                  vold(ix,iy) = vnew(ix,iy)
                  
                  ! Accumulate total calcium
                  ctoti = ctoti + ci(ix,iy)
               enddo
            enddo

c           ============================================
c           DATA OUTPUT
c           ============================================
            ! Write voltage data every 10 time steps (1 ms intervals)
            if(mod(ncount,10).eq.0) then
               write(10,*) t, v(3,3)  ! voltage at center cell (3,3 for 5x5 grid)
            endif 

            t = t + dt  ! increment total time

         enddo  ! end time step loop
      enddo     ! end beat loop

500   stop
      end

c     ================================================================
c     SPATIAL DIFFUSION SUBROUTINE
c     ================================================================
c     
c     Implements voltage diffusion between neighboring cells using
c     explicit finite differences with no-flux boundary conditions.
c     ================================================================
      
      subroutine Euler_forward(LLx,LLy,v,vnew,slmbdax,slmbday)
      implicit double precision (a-h,o-z)
        
      double precision  v(0:LLx+1,0:LLy+1),vnew(0:LLx+1,0:LLy+1)

c     ============================================================
c     SET BOUNDARY CONDITIONS (NO-FLUX)
c     ============================================================
      ! Corner points
      v(0,0) = v(2,2)
      v(0,LLy+1) = v(2,LLy-1)
      v(LLx+1,0) = v(LLx-1,2)
      v(LLx+1,LLy+1) = v(LLx-1,LLy-1)

      ! Top and bottom edges
      do ix = 1, LLx
         v(ix,0) = v(ix,2)           ! bottom edge
         v(ix,LLy+1) = v(ix,LLy-1)   ! top edge
      enddo
  
      ! Left and right edges
      do iy = 1, LLy
         v(0,iy) = v(2,iy)           ! left edge
         v(LLx+1,iy) = v(LLx-1,iy)   ! right edge
      enddo

c     ============================================================
c     FIRST DIFFUSION STEP
c     ============================================================
      ! Apply explicit Euler diffusion: ∂V/∂t = D∇²V
      do ix = 1, LLx
         do iy = 1, LLy
          vnew(ix,iy) = v(ix,iy) + slmbdax*(v(ix+1,iy)+v(ix-1,iy)
     #        -2.0d0*v(ix,iy))
     #       +slmbday*(v(ix,iy+1)+v(ix,iy-1)-2.0d0*v(ix,iy))
         enddo
      enddo

c     ============================================================
c     UPDATE BOUNDARY CONDITIONS FOR VNEW
c     ============================================================
      ! Corner points
      vnew(0,0) = vnew(2,2)
      vnew(0,LLy+1) = vnew(2,LLy-1)
      vnew(LLx+1,0) = vnew(LLx-1,2)
      vnew(LLx+1,LLy+1) = vnew(LLx-1,LLy-1)

      ! Edges
      do ix = 1, LLx
         vnew(ix,0) = vnew(ix,2)
         vnew(ix,LLy+1) = vnew(ix,LLy-1)
      enddo

      do iy = 1, LLy
         vnew(0,iy) = vnew(2,iy)
         vnew(LLx+1,iy) = vnew(LLx-1,iy)
      enddo

c     ============================================================
c     SECOND DIFFUSION STEP (FOR IMPROVED ACCURACY)
c     ============================================================
      ! Apply diffusion operator again using updated values
      do ix = 1, LLx
         do iy = 1, LLy
          v(ix,iy) = vnew(ix,iy) + slmbdax*(vnew(ix+1,iy)+vnew(ix-1,iy)
     #        -2.0d0*vnew(ix,iy))
     #      +slmbday*(vnew(ix,iy+1)+vnew(ix,iy-1)-2.0d0*vnew(ix,iy))
         enddo
      enddo

      return
      end

c     ================================================================
c     END OF CARDIAC ELECTROPHYSIOLOGY SIMULATION
c     ================================================================
c
c     SUMMARY:
c     This code simulates 2D cardiac tissue with detailed calcium
c     handling, stochastic calcium sparks, and electrical propagation.
c     Key features:
c     - Two-zone cellular model (boundary vs interior)
c     - Stochastic calcium spark dynamics
c     - Multiple ion channel types with realistic kinetics
c     - Spatial electrical coupling via gap junctions
c     - Model selection for normal vs enhanced RyR2 leak conditions
c     ================================================================
