c     ================================================================
c     ION CHANNEL SUBROUTINES FOR CARDIAC ELECTROPHYSIOLOGY
c     ================================================================
c     
c     This file contains subroutines for calculating ionic currents
c     in cardiac myocytes. Each subroutine implements the kinetics
c     and conductance of a specific ion channel type.
c     
c     Subroutines included:
c     - ina:   Fast sodium current (INa)
c     - ikr:   Rapid delayed rectifier K+ current (IKr) 
c     - iks:   Slow delayed rectifier K+ current (IKs)
c     - ik1:   Inward rectifier K+ current (IK1)
c     - ito:   Transient outward K+ current (Ito fast + slow)
c     - inak:  Na+/K+ pump current (INaK)
c     ================================================================

c     ================================================================
c     FAST SODIUM CURRENT (INa)
c     ================================================================
c     
c     Calculates the fast sodium current responsible for the upstroke
c     of the cardiac action potential. Uses Hodgkin-Huxley formalism
c     with activation (m) and inactivation (h,j) gates.
c     
c     INPUT:  hode - integration time step (ms)
c             v    - membrane voltage (mV)
c             frt  - F/RT for Nernst calculations
c             xh   - fast inactivation gate (0-1)
c             xj   - slow inactivation gate (0-1) 
c             xm   - activation gate (0-1)
c             xnai - intracellular Na+ concentration (mM)
c             xnao - extracellular Na+ concentration (mM)
c     OUTPUT: xina - sodium current (μA/μF)
c             xh,xj,xm - updated gate variables
c     ================================================================

      subroutine ina(hode,v,frt,xh,xj,xm,xnai,xnao,xina)
      implicit double precision (a-h,o-z)

c     ============================================================
c     SODIUM CHANNEL PARAMETERS
c     ============================================================
      gna = 6.0d0         ! Maximum Na+ conductance (mS/μF)
      XKMCAM = 0.3d0      ! CaM-dependent modulation constant (unused)
      deltax = -0.18d0    ! CaM shift parameter (unused)
       
c     ============================================================
c     REVERSAL POTENTIAL CALCULATION
c     ============================================================
      ena = (1.0d0/frt)*dlog(xnao/xnai)  ! Nernst potential for Na+

c     ============================================================
c     ACTIVATION GATE (m) KINETICS
c     ============================================================
      ! Rate constants independent of voltage
      am = 0.32d0*(v+47.13d0)/(1.0d0-dexp(-0.1d0*(v+47.13d0)))
      bm = 0.08d0*dexp(-v/11.0d0)

c     ============================================================
c     CaM-DEPENDENT MODULATION (CURRENTLY DISABLED)
c     ============================================================
c     Original CaM modulation code (commented out):
c     camfact = 1.0d0/(1.0d0+(XKMCAM/caM)**4)
c     vshift = -3.25*camfact
      
      camfact = 0.0d0     ! CaM modulation factor (disabled)
      vshift = 0.0d0      ! Voltage shift (disabled)
      
      vx = v - vshift     ! Effective voltage after shift

c     ============================================================
c     INACTIVATION GATES (h,j) KINETICS
c     ============================================================
      ! Voltage-dependent rate constants with different formulations
      ! for depolarized vs. hyperpolarized potentials
      
      if(vx.lt.(-40.0d0)) then
         ! Hyperpolarized potential (vx < -40 mV)
         ah = 0.135*dexp((80.0+vx)/(-6.8d0))
         bh = 3.56*dexp(0.079*vx) + 310000.0d0*dexp(0.35d0*vx)

         ! j-gate components for hyperpolarized potentials
         aj1a = -127140.0*dexp(0.2444*vx)
         aj1b = 0.00003474d0*dexp(-0.04391d0*vx)
         aj1c = (vx+37.78)/(1.0d0+dexp(0.311*(vx+79.23)))

         aj = (1.0d0+camfact*deltax)*(aj1a-aj1b)*aj1c
         
         rzx=0.1212*dexp(-0.01052*vx)
         
        bj=rzx/(1.0+dexp(-0.1378d0*(vx+40.14d0)))

      else
         ! Depolarized potential (vx >= -40 mV)
         ah = 0.0d0
         bh = 1.0d0/(0.13d0*(1.0d0+dexp((vx+10.66)/(-11.1d0))))
         aj = 0.0d0

         ! j-gate components for depolarized potentials  
         bja = 0.3*dexp(-0.0000002535d0*vx)
         bjb = 1.0+dexp(-0.1d0*(vx+32.0d0))
         bj = bja/bjb
         
      endif          

c     ============================================================
c     TIME CONSTANTS AND CURRENT CALCULATION
c     ============================================================
      tauh = 1.0d0/(ah+bh)    ! h-gate time constant
      tauj = 1.0d0/(aj+bj)    ! j-gate time constant  
      taum = 1.0d0/(am+bm)    ! m-gate time constant

      ! Sodium current: INa = gNa * m³ * h * j * (V - ENa)
      xina = gna*xh*xj*xm*xm*xm*(v-ena)

c     ============================================================
c     GATE VARIABLE UPDATES (EXPONENTIAL INTEGRATION)
c     ============================================================
      xh = ah/(ah+bh) - ((ah/(ah+bh))-xh)*dexp(-hode/tauh)
      xj = aj/(aj+bj) - ((aj/(aj+bj))-xj)*dexp(-hode/tauj)
      xm = am/(am+bm) - ((am/(am+bm))-xm)*dexp(-hode/taum)

      return
      end 

c     ================================================================
c     RAPID DELAYED RECTIFIER POTASSIUM CURRENT (IKr)
c     ================================================================
c     
c     Calculates the rapidly activating delayed rectifier K+ current,
c     which contributes to action potential repolarization.
c     
c     INPUT:  hode - integration time step (ms)
c             v    - membrane voltage (mV) 
c             frt  - F/RT for Nernst calculations
c             xko  - extracellular K+ concentration (mM)
c             xki  - intracellular K+ concentration (mM)
c             xr   - activation gate variable (0-1)
c     OUTPUT: xikr - IKr current (μA/μF)
c             xr   - updated gate variable
c     ================================================================

      subroutine ikr(hode,v,frt,xko,xki,xr,xikr)
      implicit double precision (a-h,o-z)

c     ============================================================
c     REVERSAL POTENTIAL AND CONDUCTANCE SCALING
c     ============================================================
      ek = (1.0d0/frt)*dlog(xko/xki)  ! K+ reversal potential
      
      gss = dsqrt(xko/5.40)           ! [K+]o-dependent conductance scaling
      gkr = 0.007836d0*1.3            ! Maximum IKr conductance (mS/μF)

c     ============================================================
c     ACTIVATION KINETICS
c     ============================================================
      ! Voltage-dependent rate constants for activation
      xkrv1 = 0.00138d0*(v+7.0d0)/(1.0-dexp(-0.123*(v+7.0d0)))
      xkrv2 = 0.00061d0*(v+10.0d0)/(dexp(0.145d0*(v+10.0d0))-1.0d0)
      taukr = 1.0d0/(xkrv1+xkrv2)     ! Time constant
      
      xkrinf = 1.0d0/(1.0d0+dexp(-(v+50.0d0)/7.5d0))  ! Steady-state

c     ============================================================
c     INACTIVATION (RECTIFICATION) 
c     ============================================================
      rg = 1.0d0/(1.0d0+dexp((v+33.0d0)/22.4d0))  ! Inactivation factor

c     ============================================================
c     CURRENT CALCULATION AND GATE UPDATE
c     ============================================================
      xikr = gkr*gss*xr*rg*(v-ek)  ! IKr = gKr * [Ko]^0.5 * r * rect * (V-EK)
      
      xr = xkrinf - (xkrinf-xr)*dexp(-hode/taukr)  ! Gate update

      return
      end 

c     ================================================================
c     SLOW DELAYED RECTIFIER POTASSIUM CURRENT (IKs)
c     ================================================================
c     
c     Calculates the slowly activating delayed rectifier K+ current.
c     Important for action potential duration and repolarization reserve.
c     
c     INPUT:  hode - integration time step (ms)
c             v    - membrane voltage (mV)
c             frt  - F/RT for Nernst calculations  
c             ci   - intracellular Ca2+ concentration (μM)
c             xnao,xnai - extracellular/intracellular Na+ (mM)
c             xko,xki   - extracellular/intracellular K+ (mM)
c             xs1  - activation gate variable (0-1)
c             qks  - Ca2+-dependent modulation variable (0-1)
c     OUTPUT: xiks - IKs current (μA/μF)
c             xs1,qks - updated variables
c     ================================================================

      subroutine iks(hode,v,frt,ci,xnao,xnai,xko,xki,xs1,qks,xiks)
      implicit double precision (a-h,o-z)

c     ============================================================
c     CHANNEL PARAMETERS
c     ============================================================
      prnak = 0.018330d0      ! Na+ permeability relative to K+
      gksx = 0.200d0*1.0      ! Maximum IKs conductance (mS/μF)
      
c     ============================================================
c     Ca2+-DEPENDENT MODULATION PARAMETERS
c     ============================================================
      qks_inf = 0.6d0*0.30d0  ! Steady-state Ca2+ modulation factor
      tauqks = 1000.0d0       ! Time constant for Ca2+ modulation (ms)

c     ============================================================
c     REVERSAL POTENTIAL (K+ WITH Na+ PERMEABILITY)
c     ============================================================
      eks = (1.0d0/frt)*dlog((xko+prnak*xnao)/(xki+prnak*xnai))

c     ============================================================
c     ACTIVATION KINETICS
c     ============================================================
      xs1ss = 1.0/(1.0+dexp(-(v-1.50d0)/16.70d0))  ! Steady-state activation
      xs2ss = xs1ss                                  ! Second gate (same kinetics)

      ! Time constant calculation
      tauxs=1.0d0/(0.0000719*(v+30.0d0)/(1.0d0-dexp(
     +-0.148d0*(v+30.0)))+0.000131d0
     + *(v+30.0d0)/(dexp(0.0687d0*(v+30.0d0))-1.0d0))
     
c     ============================================================
c     CURRENT CALCULATION AND VARIABLE UPDATES
c     ============================================================
      xiks = gksx*qks*xs1**2*(v-eks)  ! IKs = gKs * qKs * xs1² * (V-EKs)

      xs1 = xs1ss - (xs1ss-xs1)*dexp(-hode/tauxs)     ! Activation gate update
      qks = qks + hode*(qks_inf-qks)/tauqks            ! Ca2+ modulation update

      return
      end 

c     ================================================================
c     INWARD RECTIFIER POTASSIUM CURRENT (IK1)
c     ================================================================
c     
c     Calculates the inward rectifier K+ current that maintains
c     resting potential and provides strong repolarizing drive
c     at negative potentials.
c     
c     INPUT:  hode - integration time step (ms) [not used]
c             v    - membrane voltage (mV)
c             frt  - F/RT for Nernst calculations
c             xki  - intracellular K+ concentration (mM)
c             xko  - extracellular K+ concentration (mM)
c     OUTPUT: xik1 - IK1 current (μA/μF)
c     ================================================================

      subroutine ik1(hode,v,frt,xki,xko,xik1)
      implicit double precision (a-h,o-z)

c     ============================================================
c     REVERSAL POTENTIAL AND CONDUCTANCE
c     ============================================================
      ek = (1.0d0/frt)*dlog(xko/xki)     ! K+ reversal potential
      
      gkix = 0.60d0*1.50                 ! Base IK1 conductance (mS/μF)
      gki = gkix*(dsqrt(xko/5.4))        ! [K+]o-dependent scaling

c     ============================================================
c     RECTIFICATION KINETICS
c     ============================================================
      ! Forward and backward rate constants for rectification
      aki = 1.02/(1.0+dexp(0.2385*(v-ek-59.215)))
      
      bki = (0.49124*dexp(0.08032*(v-ek+5.476))
     +      + dexp(0.061750*(v-ek-594.31)))
     +      / (1.0+dexp(-0.5143*(v-ek+4.753)))
      
      xkin = aki/(aki+bki)               ! Rectification factor

c     ============================================================
c     CURRENT CALCULATION
c     ============================================================
      xik1 = gki*xkin*(v-ek)             ! IK1 = gK1 * rect * (V-EK)

      return
      end 

c     ================================================================
c     TRANSIENT OUTWARD POTASSIUM CURRENT (Ito)
c     ================================================================
c     
c     Calculates both fast and slow components of the transient
c     outward K+ current. Ito contributes to early repolarization
c     and action potential notch formation.
c     
c     INPUT:  hode - integration time step (ms)
c             v    - membrane voltage (mV)
c             frt  - F/RT for Nernst calculations
c             xki,xko - intracellular/extracellular K+ (mM)
c             xtof,ytof - fast Ito activation/inactivation gates
c             xtos,ytos - slow Ito activation/inactivation gates  
c             gtof,gtos - fast/slow Ito conductances (mS/μF)
c     OUTPUT: xito - total Ito current (μA/μF)
c             xtof,ytof,xtos,ytos - updated gate variables
c     ================================================================

      subroutine ito(hode,v,frt,xki,xko,xtof,ytof,xtos,ytos,
     &               xito,gtof,gtos)
      implicit double precision (a-h,o-z)

c     ============================================================
c     REVERSAL POTENTIAL
c     ============================================================
      ek = (1.0d0/frt)*dlog(xko/xki)     ! K+ reversal potential

c     ============================================================
c     COMMON VOLTAGE-DEPENDENT FACTORS
c     ============================================================
      rt1 = -(v+3.0)/15.0d0              ! Activation voltage dependence
      rt2 = (v+33.5)/10.0d0              ! Inactivation voltage dependence
      rt3 = (v+60.0d0)/10.0d0            ! Slow inactivation time constant

c     ============================================================
c     SLOW Ito COMPONENT (Itos)
c     ============================================================
      ! Steady-state values
      xtos_inf = 1.0d0/(1.0+dexp(rt1))   ! Activation steady-state
      ytos_inf = 1.0d0/(1.0d0+dexp(rt2)) ! Inactivation steady-state
      rs_inf = 1.0d0/(1.0d0+dexp(rt2))   ! Additional inactivation factor

      ! Time constants
      txs = 9.0d0/(1.0d0+dexp(-rt1)) + 0.5d0     ! Activation time constant
      tys = 3000.0d0/(1.0+dexp(rt3)) + 30.0d0    ! Inactivation time constant

      ! Current calculation: Itos includes additional rectification term
      xitos = gtos*xtos*(ytos+0.5d0*rs_inf)*(v-ek)

      ! Gate updates
      xtos = xtos_inf - (xtos_inf-xtos)*dexp(-hode/txs)
      ytos = ytos_inf - (ytos_inf-ytos)*dexp(-hode/tys)

c     ============================================================
c     FAST Ito COMPONENT (Itof) - Shannon et al. 2005
c     ============================================================
      ! Steady-state values (same activation as slow component)
      xtof_inf = xtos_inf                ! Same activation kinetics
      ytof_inf = ytos_inf                ! Same inactivation kinetics

      ! Time constants (faster than slow component)
      rt4 = -(v/30.0d0)*(v/30.0d0)      ! Gaussian voltage dependence
      rt5 = (v+33.5d0)/10.0d0           ! Inactivation voltage dependence
      
      txf = 3.5d0*dexp(rt4) + 1.5d0     ! Fast activation time constant
      tyf = 20.0/(1.0+dexp(rt5)) + 20.0d0  ! Fast inactivation time constant

      ! Current calculation: standard Hodgkin-Huxley form
      xitof = gtof*xtof*ytof*(v-ek)

      ! Gate updates  
      xtof = xtof_inf - (xtof_inf-xtof)*dexp(-hode/txf)
      ytof = ytof_inf - (ytof_inf-ytof)*dexp(-hode/tyf)

c     ============================================================
c     TOTAL Ito CURRENT
c     ============================================================
      xito = xitos + xitof               ! Total = slow + fast components

      return
      end 

c     ================================================================
c     SODIUM-POTASSIUM PUMP CURRENT (INaK)
c     ================================================================
c     
c     Calculates the Na+/K+ pump (ATPase) current that maintains
c     ionic gradients. Exchanges 3 Na+ out for 2 K+ in, creating
c     a net outward current.
c     
c     INPUT:  v     - membrane voltage (mV)
c             frt   - F/RT for Nernst calculations
c             xko   - extracellular K+ concentration (mM)
c             xnao  - extracellular Na+ concentration (mM) 
c             xnai  - intracellular Na+ concentration (mM)
c     OUTPUT: xinak - INaK current (μA/μF)
c     ================================================================

      subroutine inak(v,frt,xko,xnao,xnai,xinak)
      implicit double precision (a-h,o-z)

c     ============================================================
c     PUMP PARAMETERS
c     ============================================================
      xibarnak = 1.50d0      ! Maximum pump rate (μA/μF)
      xkmko = 1.5d0          ! K+ half-saturation constant (mM)
      xkmnai = 12.0d0        ! Na+ half-saturation constant (mM)
      hh = 1.0d0             ! Na+ cooperativity exponent

c     ============================================================
c     VOLTAGE-DEPENDENT MODULATION
c     ============================================================
      ! Accounts for effects of membrane potential and [Na+]o on pump activity
      sigma = (dexp(xnao/67.3d0)-1.0d0)/7.0d0
      
      fnak = 1.0d0/(1+0.1245*dexp(-0.1*v*frt)
     +              +0.0365*sigma*dexp(-v*frt))

c     ============================================================
c     PUMP CURRENT CALCULATION
c     ============================================================
      ! INaK depends on:
      ! - Voltage (fnak)
      ! - Intracellular [Na+] with cooperativity
      ! - Extracellular [K+] saturation
      xinak = xibarnak*fnak*(1./(1.+(xkmnai/xnai)**hh))*xko/(xko+xkmko)

      return
      end 

c     ================================================================
c     END OF ION CHANNEL SUBROUTINES
c     ================================================================
