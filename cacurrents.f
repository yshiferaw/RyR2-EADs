c     ================================================================
c     CALCIUM HANDLING AND STOCHASTIC PROCESS SUBROUTINES
c     ================================================================
c     
c     Supporting subroutines for cardiac calcium dynamics simulation.
c     Includes stochastic spark generation, calcium buffering,
c     SR uptake, Na-Ca exchange, L-type Ca channels, and Markov
c     state modeling.
c     
c     Subroutines included:
c     - binom:   Binomial random number generation
c     - binevol: Stochastic evolution of spark clusters  
c     - uptake:  SERCA pump calcium uptake
c     - total:   Free to total calcium conversion (buffering)
c     - xfree:   Total to free calcium conversion  
c     - inaca:   Na-Ca exchanger current
c     - ica:     L-type calcium current and driving force
c     - markov:  L-type channel Markov state dynamics
c     ================================================================

c     ================================================================
c     BINOMIAL RANDOM NUMBER GENERATOR
c     ================================================================
c     
c     Generates binomial random numbers using the inverse transform
c     method. Used for stochastic calcium spark generation.
c     
c     INPUT:  n - number of trials
c             p - probability of success per trial
c     OUTPUT: m - number of successes (binomial random variate)
c     ================================================================

      subroutine binom(n,p,m)
      implicit double precision (a-h,o-z)

c     ============================================================
c     INITIALIZE ALGORITHM PARAMETERS
c     ============================================================
      call random_number(r1)

      xlog = -dlog(1.0d0-p)    ! Threshold for acceptance
      xsum = 0.0d0             ! Cumulative sum
      i = 1                    ! Trial counter

c     ============================================================
c     INVERSE TRANSFORM SAMPLING LOOP
c     ============================================================
10    call random_number(r1)
      
      xx1 = r1                        ! Random number [0,1]
      xx2 = dfloat(n-i+1)            ! Remaining trials
      hx = -dlog(xx1)/xx2            ! Exponential increment
      xsum = xsum + hx               ! Update cumulative sum
      
      if(xsum.gt.xlog) goto 100      ! Accept if threshold exceeded
      
      i = i + 1                      ! Next trial
      goto 10	
      
100   m = i - 1                      ! Number of successes

      return
      end 

c     ================================================================
c     STOCHASTIC SPARK CLUSTER EVOLUTION
c     ================================================================
c     
c     Evolves the number of active calcium spark clusters using
c     binomial processes for spark initiation and termination.
c     
c     INPUT:  nt      - total number of clusters
c             nx      - current number of active clusters
c             alpha   - spark initiation rate (1/ms)
c             beta    - spark termination rate (1/ms)  
c             dt      - time step (ms)
c     OUTPUT: ndeltap - number of new sparks initiated
c             ndeltam - number of sparks terminated
c     ================================================================
      
      subroutine binevol(nt,nx,alpha,beta,dt,ndeltap,ndeltam)                  
      implicit double precision (a-h,o-z)

c     ============================================================
c     SPARK INITIATION PROCESS
c     ============================================================
      ! Skip calculation if initiation rate is very low (numerical stability)
      if(alpha.lt.0.00001) then
         ndeltap = 0
         goto 200
      endif 

      na = nt - nx                   ! Available clusters for sparking
      
      if(na.gt.0) then	
         xrate = alpha*dt            ! Initiation probability per cluster
         call binom(na,xrate,ndeltap) ! Sample new sparks
      else
         ndeltap = 0                 ! No clusters available
      endif 
      
c     ============================================================
c     SPARK TERMINATION PROCESS  
c     ============================================================
200   if(nx.gt.0) then
         xrate = beta*dt             ! Termination probability per spark
         call binom(nx,xrate,ndeltam) ! Sample spark terminations
      else
         ndeltam = 0                 ! No active sparks to terminate
      endif 
      
      return
      end

c     ================================================================
c     SERCA PUMP CALCIUM UPTAKE
c     ================================================================
c     
c     Calculates calcium uptake rate by the sarcoplasmic reticulum
c     Ca-ATPase (SERCA) pump with Hill kinetics.
c     
c     INPUT:  ci  - cytosolic calcium concentration (μM)
c             vup - maximum uptake rate (μM/ms)
c     OUTPUT: xup - actual uptake rate (μM/ms)
c     ================================================================

      subroutine uptake(ci,vup,xup)
      implicit double precision (a-h,o-z)

c     ============================================================
c     SERCA PUMP PARAMETERS
c     ============================================================
      double precision, parameter::Ki = 0.30d0    ! Half-saturation (μM)
      double precision, parameter::Knsr = 800.0d0 ! SR load dependence (μM)
      double precision, parameter::HH = 3.00d0    ! Hill coefficient

c     ============================================================
c     HILL KINETICS CALCULATION
c     ============================================================
      ! Michaelis-Menten kinetics with cooperativity:
      ! J_up = V_max * [Ca]^n / (K_m^n + [Ca]^n)
      xup = vup*ci**HH/(Ki**HH+ci**HH)
      
      return
      end

c     ================================================================
c     FREE TO TOTAL CALCIUM CONVERSION (BUFFERING)
c     ================================================================
c     
c     Converts free calcium concentration to total calcium including
c     calcium bound to buffers (calmodulin and SR binding sites).
c     
c     INPUT:  ci  - free calcium concentration (μM)
c     OUTPUT: cit - total calcium concentration (μM)
c     ================================================================

      subroutine total(ci,cit)
      implicit double precision (a-h,o-z)
       		
c     ============================================================
c     CALCIUM BUFFER PARAMETERS
c     ============================================================
      bcal = 24.0d0      ! Calmodulin total concentration (μM)
      xkcal = 7.0d0      ! Calmodulin binding affinity (μM)
      srmax = 47.0d0     ! SR buffer total concentration (μM)
      srkd = 0.6d0       ! SR buffer binding affinity (μM)
       
c     ============================================================
c     BUFFER BINDING CALCULATIONS
c     ============================================================
      bix = bcal*ci/(xkcal+ci)     ! Calmodulin-bound calcium
      six = srmax*ci/(srkd+ci)     ! SR-bound calcium

c     ============================================================
c     TOTAL CALCIUM CONCENTRATION
c     ============================================================
      cit = ci + bix + six         ! Free + bound calcium

      return
      end 

c     ================================================================
c     TOTAL TO FREE CALCIUM CONVERSION
c     ================================================================
c     
c     Converts total calcium concentration back to free calcium
c     using algebraic solution of buffer equilibrium equations.
c     
c     INPUT:  cit - total calcium concentration (μM)
c     OUTPUT: ci  - free calcium concentration (μM)
c     ================================================================

      subroutine xfree(cit,ci)
      implicit double precision (a-h,o-z)
       		
c     ============================================================
c     ALGEBRAIC SOLUTION COEFFICIENTS
c     ============================================================
      ! Coefficients for quadratic solution of buffer equilibrium
      a = 2.23895        ! Coefficient derived from buffer parameters
      b = 52.0344        ! Coefficient derived from buffer parameters  
      c = 0.666509       ! Coefficient derived from buffer parameters

      y = cit            ! Total calcium input

c     ============================================================
c     QUADRATIC FORMULA SOLUTION
c     ============================================================
      ! Solve: a*ci^2 + b*ci + c*ci - y = 0 for ci
      xa = (b+a*c-y)**2 + 4.0*a*c*y
      ci = (-b-a*c+y+dsqrt(xa))/(2.0*a)

      return
      end 

c     ================================================================
c     SODIUM-CALCIUM EXCHANGER CURRENT
c     ================================================================
c     
c     Calculates Na-Ca exchange current using detailed kinetic model
c     with voltage dependence and ion concentration dependencies.
c     
c     INPUT:  v     - membrane voltage (mV)
c             frt   - F/RT for voltage dependence
c             xnai  - intracellular Na+ concentration (mM)
c             xnao  - extracellular Na+ concentration (mM)
c             cao   - extracellular Ca2+ concentration (mM)
c             ci    - intracellular Ca2+ concentration (μM)
c     OUTPUT: xinacaq - NCX current density (dimensionless)
c     ================================================================

      subroutine inaca(v,frt,xnai,xnao,cao,ci,xinacaq)
      implicit double precision (a-h,o-z)
       	
c     ============================================================
c     CONCENTRATION UNIT CONVERSION
c     ============================================================
      cim = ci/1000.0d0              ! Convert μM to mM

c     ============================================================
c     VOLTAGE-DEPENDENT DRIVING FORCE
c     ============================================================
      ! Forward and reverse fluxes with voltage dependence
      zw3a = xnai**3*cao*dexp(v*0.35*frt)        ! Forward flux
      zw3b = xnao**3*cim*dexp(v*(0.35-1.)*frt)   ! Reverse flux
      zw3 = zw3a - zw3b                          ! Net driving force
      
      zw4 = 1.0d0 + 0.2d0*dexp(v*(0.35d0-1.0d0)*frt) ! Voltage factor

c     ============================================================
c     CALCIUM-DEPENDENT REGULATION
c     ============================================================
      xkdna = 0.3d0                              ! Ca regulation constant (μM)
      aloss = 1.0d0/(1.0d0+(xkdna/ci)**3)       ! Ca-dependent activation

c     ============================================================
c     SATURATION CONSTANTS
c     ============================================================
      xmcao = 1.3d0       ! Extracellular Ca2+ saturation (mM)
      xmnao = 87.5d0      ! Extracellular Na+ saturation (mM)
      xmnai = 12.3d0      ! Intracellular Na+ saturation (mM)
      xmcai = 0.0036d0    ! Intracellular Ca2+ saturation (mM)

c     ============================================================
c     DENOMINATOR TERMS (COMPETITIVE BINDING)
c     ============================================================
      yz1 = xmcao*xnai**3 + xmnao**3*cim
      yz2 = xmnai**3*cao*(1.0d0+cim/xmcai)
      yz3 = xmcai*xnao**3*(1.0d0+(xnai/xmnai)**3)
      yz4 = xnai**3*cao + xnao**3*cim
      zw8 = yz1 + yz2 + yz3 + yz4               ! Total denominator

c     ============================================================
c     NCX CURRENT CALCULATION
c     ============================================================
      xinacaq = aloss*zw3/(zw4*zw8)

      return
      end 

c     ================================================================
c     L-TYPE CALCIUM CURRENT AND DRIVING FORCE
c     ================================================================
c     
c     Calculates L-type calcium current using Goldman-Hodgkin-Katz
c     equation for calcium permeation with voltage dependence.
c     
c     INPUT:  v   - membrane voltage (mV)
c             frt - F/RT for voltage dependence
c             cao - extracellular Ca2+ concentration (mM)
c             ci  - intracellular Ca2+ concentration (μM)
c             pox - total channel open probability (0-1)
c     OUTPUT: rca   - calcium driving force
c             xicaq - L-type Ca current density
c     ================================================================

      subroutine ica(v,frt,cao,ci,pox,rca, xicaq)
      implicit double precision (a-h,o-z)

c     ============================================================
c     PHYSICAL CONSTANTS AND PARAMETERS
c     ============================================================
      xf = 96.485d0        ! Faraday's constant (C/mol)
      pca = 0.00054d0      ! Ca2+ permeability constant (Luo-Rudy)
      za = v*2.0d0*frt     ! Voltage factor (z=2 for Ca2+)
      
      factor1 = 4.0d0*pca*xf*frt          ! GHK equation prefactor
      factor = v*factor1                   ! Voltage-dependent factor

c     ============================================================
c     CONCENTRATION UNIT CONVERSION
c     ============================================================
      cim = ci/1000.0d0    ! Convert μM to mM for consistency
             
c     ============================================================
c     GOLDMAN-HODGKIN-KATZ DRIVING FORCE
c     ============================================================
      ! Handle special case when voltage approaches zero
      if(dabs(za).lt.0.001d0) then
         ! L'Hôpital's rule limit as V→0
         rca = factor1*(cim*dexp(za)-0.341d0*(cao))/(2.0d0*frt)
      else   
         ! Standard GHK equation
         rca = factor*(cim*dexp(za)-0.341d0*(cao))/(dexp(za)-1.0d0)         
      endif 

c     ============================================================
c     TOTAL L-TYPE CALCIUM CURRENT
c     ============================================================
      xicaq = rca*pox      ! Current = driving force × open probability
      
      return
      end 

c     ================================================================
c     L-TYPE CALCIUM CHANNEL MARKOV STATE MODEL
c     ================================================================
c     
c     Models L-type calcium channel gating using a 10-state Markov
c     model with calcium-dependent inactivation and spark coupling.
c     States: C2, C1, O, I1, I2 (boundary) + C2s, C1s, Os, I1s, I2s (spark)
c     
c     INPUT:  hode - time step (ms)
c             v    - membrane voltage (mV)  
c             ci   - intracellular Ca2+ (μM)
c             c1,c2,xi1,xi2,po - boundary-facing channel states
c             c1s,c2s,xi1s,xi2s,pos - spark-facing channel states
c             alpha - spark initiation rate (1/ms)
c             bts   - spark termination rate (1/ms)
c             zxr   - Ca2+ sensitivity parameter
c     OUTPUT: Updated channel state probabilities
c     ================================================================
                   
      subroutine markov(hode,v,ci,c1,c2,xi1,xi2,po,
     + c1s,c2s,xi1s,xi2s,pos,alpha,bts,zxr)
      implicit double precision (a-h,o-z)

c     ============================================================
c     CALCIUM-INDEPENDENT TRANSITION RATES  
c     ============================================================
      a23 = 0.3            ! C1 → O rate
      a32 = 3.0            ! O → C1 rate
      a42 = 0.00224d0      ! I1 → C1 rate

c     ============================================================
c     VOLTAGE-DEPENDENT ACTIVATION (C2 ↔ C1)
c     ============================================================
      vth = 1.0            ! Activation threshold (mV)
      s6 = 6.50            ! Activation slope factor
     
      poinf = 1.0d0/(1.0d0+dexp(-(v-vth)/s6))  ! Steady-state activation
      taupo = 1.0d0                            ! Activation time constant

      a12 = poinf/taupo           ! C2 → C1 rate
      a21 = (1.0d0-poinf)/taupo   ! C1 → C2 rate

c     ============================================================
c     VOLTAGE-DEPENDENT INACTIVATION (I2 ↔ C2)
c     ============================================================
      vy = -40.0d0         ! Inactivation voltage
      sy = 4.0d0           ! Inactivation slope
      prv = 1.0d0-1.0d0/(1.0d0+dexp(-(v-vy)/sy))  ! Voltage factor

      vyr = -40.0d0        ! Recovery voltage threshold
      syr = 10.0d0         ! Recovery slope

      recovx = 10.0d0 + 4954.0d0*dexp(v/15.6d0)   ! Recovery rate
      recov = 1.5*recovx                          ! Scaled recovery
      tauba = (recov-450.0d0)*prv + 450.0d0       ! Time constant
      tauba = tauba/2.0                           ! Final scaling

      poix = 1.0d0/(1.0d0+dexp(-(v-vyr)/syr))     ! Recovery steady-state

      a15 = poix/tauba           ! C2 → I2 rate
      a51 = (1.0d0-poix)/tauba   ! I2 → C2 rate

c     ============================================================
c     ADDITIONAL INACTIVATION PATHWAY (O → I1)
c     ============================================================
      vx = -40.0d0         ! Inactivation voltage
      sx = 3.0d0           ! Slope factor
      tau3 = 3.0d0         ! Time constant
      poi = 1.0d0/(1.0d0+dexp(-(v-vx)/sx))  ! Steady-state
      a45 = (1.0d0-poi)/tau3     ! O → I1 rate

c     ============================================================
c     CALCIUM-DEPENDENT TRANSITION RATES
c     ============================================================
      cat = 1.0d0          ! Ca2+ threshold for inactivation

      fca = 1.0d0/(1.0d0+(cat/ci)**2)     ! Ca2+ dependence factor
      a24 = 0.00413d0 + zxr*fca           ! C1 → I1 rate (Ca-dependent)
      a34 = 0.00195 + zxr*fca             ! O → I1 rate (Ca-dependent)

c     ============================================================
c     DETAILED BALANCE RELATIONSHIPS
c     ============================================================
      a43 = a34*(a23/a32)*(a42/a24)       ! I1 → O rate
      a54 = a45*(a51/a15)*(a24/a42)*(a12/a21)  ! I2 → I1 rate

c     ============================================================
c     SPARK-FACING CHANNEL RATES (DIFFERENT Ca SENSITIVITY)
c     ============================================================
      fcax = 1.0d0         ! Enhanced Ca2+ sensitivity for sparks
      a24s = 0.00413d0 + zxr*fcax         ! C1s → I1s rate
      a34s = 0.00195 + zxr*fcax           ! Os → I1s rate

      a43s = a34s*(a23/a32)*(a42/a24s)    ! I1s → Os rate
      a54s = a45*(a51/a15)*(a24s/a42)*(a12/a21)  ! I2s → I1s rate

c     ============================================================
c     STATE DYNAMICS (BOUNDARY-FACING CHANNELS)
c     ============================================================
      dpo = a23*c1 + a43*xi1 - (a34+a32)*po - alpha*po + bts*pos
      dc2 = a21*c1 + a51*xi2 - (a15+a12)*c2 + bts*c2s
      dc1 = a12*c2 + a42*xi1 + a32*po - (a21+a23+a24)*c1 + bts*c1s
      dxi1 = a24*c1 + a54*xi2 + a34*po - (a45+a42+a43)*xi1 + bts*xi1s

c     ============================================================
c     STATE DYNAMICS (SPARK-FACING CHANNELS)
c     ============================================================
      dpos = a23*c1s + a43s*xi1s - (a34s+a32)*pos + alpha*po - bts*pos
      dc2s = a21*c1s + a51*xi2s - (a15+a12)*c2s - bts*c2s
      dc1s = a12*c2s + a42*xi1s + a32*pos - (a21+a23+a24s)*c1s - bts*c1s
      dxi1s = a24s*c1s + a54s*xi2s + a34s*pos - (a45+a42+a43s)*xi1s 
     +        - bts*xi1s
      dxi2s = a45*xi1s + a15*c2s - (a51+a54s)*xi2s - bts*xi2s

c     ============================================================
c     STATE VARIABLE UPDATES
c     ============================================================
      po = po + dpo*hode
      c1 = c1 + dc1*hode
      c2 = c2 + dc2*hode
      xi1 = xi1 + dxi1*hode

      pos = pos + dpos*hode
      c1s = c1s + dc1s*hode
      c2s = c2s + dc2s*hode
      xi1s = xi1s + dxi1s*hode
      xi2s = xi2s + dxi2s*hode	

c     ============================================================
c     CONSERVATION CONSTRAINT (ALL STATES SUM TO 1)
c     ============================================================
      xi2 = 1.0d0 - c1 - c2 - po - xi1 - pos - c1s - c2s - xi1s - xi2s	
         
      return
      end

c     ================================================================
c     END OF CALCIUM AND STOCHASTIC PROCESS SUBROUTINES
c     ================================================================
