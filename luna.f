      !********************************************************************************************************************************************************************** 
      !
      ! !DESCRIPTION:
      ! Calculates Nitrogen fractionation within the leaf, based on optimum calculated fractions in rubisco, cholorophyll, 
      ! Respiration and Storage. Based on Xu et al. 2012 and Ali et al 2015.In Review 
      
      !
      ! !REVISION HISTORY:
      ! version 1.0, by Chonggang Xu, Ashehad Ali and Rosie Fisher. July 14  2015.
      ! version 0.1, by Chonggang Xu, Ashehad Ali and Rosie Fisher. October 30 2014. 
      
      !********************************************************************************************************************************************************************** 
     
     
     
      !********************************************************************************************************************************************************************** 
      ! this subroutine updates the photosynthetic capacity as determined by Vcmax25 and Jmax25
     
      !subroutine Update_Photosynthesis_Capacity(bounds, fn, filterp, 
      !&dayl_factor, atm2lnd_inst, temperature_inst, canopystate_inst, 
      !&photosyns_inst, 
      !&surfalb_inst, solarabs_inst, waterstate_inst, frictionvel_inst)
     
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      !Use the LUNA model to calculate the Nitrogen partioning 
      !subroutine NitrogenAllocation(FNCa,forc_pbot10, relh10, CO2a10,O2a10, PARi10,PARimx10,rb10, hourpd, tair10, tleafd10, tleafn10, &
      !     Jmaxb0, Jmaxb1, Wc2Wjb0, relhExp,&
      !     PNstoreold, PNlcold, PNetold, PNrespold, PNcbold, &
      !     PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt)
      subroutine LUNA(z,vcmx25_z,jmx25_z,PNlc_z,enzs_z,rb10,
     &sla,lnc,par240d,par240x,CO2a10,o2a10,hourpd,relh10,tair10,  
     &gs_func,ftg0,ftg1,max_daily_pchg,ttype,ftToptV,ftHaV,ftHdV,
     &ftToptJ,ftHaJ,ftHdJ) 
     
      implicit none
     
      !REAL*8, intent (in) :: FNCa                       !Area based functional nitrogen content (g N/m2 leaf)
      REAL*8, intent (in) :: sla                        !Specific leaf area m2 g-1 DW
      REAL*8, intent (in) :: lnc                        !Area based leaf nitrogen content (g N/m2 leaf)
      REAL*8, intent (in) :: relh10                     !10-day mean relative humidity (unitless)
      REAL*8, intent (in) :: CO2a10                     !10-day meanCO2 concentration in the air (Pa)
      REAL*8, intent (in) :: O2a10                      !10-day mean O2 concentration in the air (Pa)
      REAL*8, intent (in) :: par240d                  !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
      REAL*8, intent (in) :: par240x                  !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
      REAL*8, intent (in) :: rb10                       !10-day mean boundary layer resistance
      REAL*8, intent (in) :: hourpd                     !hours of light in a the day (hrs)
      REAL*8, intent (in) :: tair10                     !10-day running mean of the 2m temperature (oC)
      REAL*8, intent (in) :: max_daily_pchg             ! maximum daily proportional change in vcmax or Jmax
      ! the below thre variables are commented out as the t pars are assumed the same as tair and pressure is assumed constant in SDGVM
      !REAL*8, intent (in) :: tleafd10                   !10-day running mean of daytime leaf temperature (oC) 
      !REAL*8, intent (in) :: tleafn10                   !10-day running mean of nighttime leaf temperature (oC) 
      !REAL*8, intent (in) :: forc_pbot10                !10-day mean air pressure (Pa)
      ! the below variables are commented out because they are declared as parameters below
      !REAL*8, intent (in) :: Jmaxb0                     !baseline proportion of nitrogen allocated for electron transport rate (unitless)
      !REAL*8, intent (in) :: Jmaxb1                     !coefficient determining the response of electron transport rate to light availability (unitless) 
      !REAL*8, intent (in) :: Wc2Wjb0                    !the baseline ratio of rubisco-limited rate vs light-limited photosynthetic rate (Wc:Wj)
      !REAL*8, intent (in) :: relhExp                    !specifies the impact of relative humidity on electron transport rate (unitless)
      ! the below variables are commented as they declared as local variables below 
      !REAL*8, intent (in) :: PNstoreold                 !old value of the proportion of nitrogen allocated to storage (unitless)
      !REAL*8, intent (in) :: PNlcold                    !old value of the proportion of nitrogen allocated to light capture (unitless)
      !REAL*8, intent (in) :: PNetold                    !old value of the proportion of nitrogen allocated to electron transport (unitless)
      !REAL*8, intent (in) :: PNrespold                  !old value of the proportion of nitrogen allocated to respiration (unitless)
      !REAL*8, intent (in) :: PNcbold                    !old value of the proportion of nitrogen allocated to carboxylation (unitless)  
      !REAL*8, intent (out):: PNstoreopt                 !optimal proportion of nitrogen for storage 
      !REAL*8, intent (out):: PNlcopt                    !optimal proportion of nitrogen for light capture 
      !REAL*8, intent (out):: PNetopt                    !optimal proportion of nitrogen for electron transport 
      !REAL*8, intent (out):: PNrespopt                  !optimal proportion of nitrogen for respiration 
      !REAL*8, intent (out):: PNcbopt                    !optial proportion of nitrogen for carboxyaltion  
      REAL*8, intent (in)   :: ftToptV                  !optimum temp for vcmax (oC)
      REAL*8, intent (in)   :: ftHaV                    !activation energy for vcmax
      REAL*8, intent (in)   :: ftHdV                    !deactivation energy for vcmax
      REAL*8, intent (in)   :: ftToptJ                  !optimum temp for jmax (oC)
      REAL*8, intent (in)   :: ftHaJ                    !activation energy for jmax
      REAL*8, intent (in)   :: ftHdJ                    !deactivation energy for jmax
      integer, intent (in)  :: z                        !canopy layer
      integer, intent (in)  :: ttype                    !temperature scaling method to use
      REAL*8, intent (inout):: vcmx25_z                 !vcmax25
      REAL*8, intent (inout):: jmx25_z                  !jmx25
      REAL*8, intent (inout):: PNlc_z                   !optimal proportion of nitrogen for carboxylation  
      REAL*8, intent (inout):: enzs_z                   !enzyme decay state
      
      !-------------------------------------------------------------------------------------------------------------------------------
      !intermediate variables
      REAL*8 :: FNCa                       !Area based functional nitrogen content (g N/m2 leaf)
      REAL*8 :: SNCa                       !Area based structural nitrogen content (g N/m2 leaf)
      REAL*8 :: par240d_z                  !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
      REAL*8 :: par240x_z                  !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
      REAL*8 :: PARi10                     !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
      REAL*8 :: PARimx10                   !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
      REAL*8 :: PNlcold                    !old value of the proportion of nitrogen allocated to light capture (unitless)
      REAL*8 :: PNetold                    !old value of the proportion of nitrogen allocated to electron transport (unitless)
      REAL*8 :: PNrespold                  !old value of the proportion of nitrogen allocated to respiration (unitless)
      REAL*8 :: PNcbold                    !old value of the proportion of nitrogen allocated to carboxylation (unitless)  
      REAL*8 :: PNstoreopt                 !optimal proportion of nitrogen for storage 
      REAL*8 :: PNlcopt                    !optimal proportion of nitrogen for light capture 
      REAL*8 :: PNetopt                    !optimal proportion of nitrogen for electron transport 
      REAL*8 :: PNrespopt                  !optimal proportion of nitrogen for respiration 
      REAL*8 :: PNcbopt                    !optial proportion of nitrogen for carboxyaltion  
      REAL*8 :: Carboncost1                             !absolute amount of carbon cost associated with maintenance respiration due to deccrease in light capture nitrogen(g dry mass per day) 
      REAL*8 :: Carboncost2                             !absolute amount of carbon cost associated with maintenance respiration due to increase in light capture nitrogen(g dry mass per day) 
      REAL*8 :: Carbongain1                             !absolute amount of carbon gain associated with maintenance respiration due to deccrease in light capture nitrogen(g dry mass per day) 
      REAL*8 :: Carbongain2                             !absolute amount of carbon gain associated with maintenance respiration due to increase in light capture nitrogen(g dry mass per day) 
      REAL*8 :: Fc                                      !the temperature adjustment factor for Vcmax 
      REAL*8 :: Fj                                      !the temperature adjustment factor for Jmax 
      REAL*8 :: PNlc                                    !the current nitrogen allocation proportion for light capture
      REAL*8 :: Jmax                                    !the maximum electron transport rate (umol/m2/s) 
      REAL*8 :: JmaxCoef                                !coefficient determining the response of electron transport rate to light availability (unitless) and humidity
      REAL*8 :: Jmaxb0act                               !base value of Jmax (umol/m2/s) 
      REAL*8 :: JmaxL                                   !the electron transport rate with maximum daily radiation (umol/m2/s)  
      REAL*8 :: JmeanL                                  !the electron transport rate with mean radiation (umol/m2/s) 
      REAL*8 :: Nstore                                  !absolute amount of nitrogen allocated to storage (gN/m2 leaf)
      REAL*8 :: Nresp                                   !absolute amount of nitrogen allocated to respiration (gN/m2 leaf) 
      REAL*8 :: Nlc                                     !absolute amount of nitrogen allocated to light capture (gN/m2 leaf) 
      REAL*8 :: Net                                     !absolute amount of nitrogen allocated to electron transport (gN/m2 leaf) 
      REAL*8 :: Ncb                                     !absolute amount of nitrogen allocated to carboxylation (gN/m2 leaf) 
      REAL*8 :: Nresp1                                  !absolute amount of nitrogen allocated to respiration due to increase in light capture nitrogen(gN/m2 leaf)  
      REAL*8 :: Nlc1                                    !absolute amount of nitrogen allocated to light capture due to increase in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: Net1                                    !absolute amount of nitrogen allocated to electron transport due to increase in light capture nitrogen(gN/m2 leaf)
      REAL*8 :: Ncb1                                    !absolute amount of nitrogen allocated to carboyxlation due to increase in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: Nresp2                                  !absolute amount of nitrogen allocated to respiration due to decrease in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: Nlc2                                    !absolute amount of nitrogen allocated to light capture due to decrease in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: Net2                                    !absolute amount of nitrogen allocated to electron transport due to decrease in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: Ncb2                                    !absolute amount of nitrogen allocated to carboxylation due to increase in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: PSN                                     !g carbon photosynthesized per day per unit(m2) of leaf
      REAL*8 :: RESP                                    !g carbon respired per day per unit(m2) of leaf due to increase in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: PSN1                                    !g carbon photosynthesized per day per unit(m2) of leaf due to increase in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: RESP1                                   !g carbon respired per day per unit(m2) of leaf due to decrease in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: PSN2                                    !g carbon photosynthesized per day per unit(m2) of leaf due to decrease in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: RESP2                                   !g carbon respired per day per unit(m2) of leaf
      REAL*8 :: Npsntarget                              !absolute amount of target nitrogen for photosynthesis(gN/m2 leaf) 
      REAL*8 :: Npsntarget1                             !absolute amount of target nitrogen for photosynthesis due to increase in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: Npsntarget2                             !absolute amount of target nitrogen for photosynthesis due to decrease in light capture nitrogen(gN/m2 leaf) 
      REAL*8 :: NUEj                                    !nitrogen use efficiency for electron transport under current environmental conditions 
      REAL*8 :: NUEc                                    !nitrogen use efficiency for carboxylation under current environmental conditions  
      REAL*8 :: NUEjref                                 !nitrogen use efficiency for electron transport under reference environmental conditions (25oC and 385ppm Co2) 
      REAL*8 :: NUEcref                                 !nitrogen use efficiency for carboxylation under reference environmental conditions (25oC and 385ppm Co2) 
      REAL*8 :: NUEr                                    !nitrogen use efficiency for respiration 
      REAL*8 :: PARi10c                                 !10-day mean constrained photosynthetic active radiation on in a canopy (umol/m2/s)
      REAL*8 :: PARimx10c                               !10-day mean constrained 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
      REAL*8 :: Kj2Kcref                                !the ratio of rubisco-limited photosynthetic rate (Wc) to light limited photosynthetic rate (Wj)
      REAL*8 :: PNlcoldi                                !old value of the proportion of nitrogen allocated to light capture (unitless) 
      REAL*8 :: Kj2Kc                                   !the ratio of Wc to Wj under changed conditions 
      REAL*8 :: Kc                                      !conversion factors for Vc,max to Wc 
      REAL*8 :: Kj                                      !conversion factor for electron transport rate to Wj 
      REAL*8 :: theta                                   !efficiency of light energy conversion (unitless) 
      REAL*8 :: chg                                     !the nitrogen change per interation
      REAL*8 :: chg_constrn                             !the nitrogen change per interation
      REAL*8 :: chg_per_step                            !the nitrogen change per interation
      REAL*8 :: Vcmaxnight                              !Vcmax during night (umol/m2/s)
      REAL*8 :: ci                                      !inter-cellular CO2 concentration (Pa)
      !REAL*8 :: theta_cj                                !interpolation coefficient
      REAL*8 :: tleafd10                   !10-day running mean of daytime leaf temperature (oC) 
      REAL*8 :: tleafn10                   !10-day running mean of nighttime leaf temperature (oC) 
      REAL*8 :: forc_pbot10                !10-day mean air pressure (Pa)
      REAL*8 :: tleafd10c                               !10-day mean daytime leaf temperature, contrained for physiological range (oC)
      REAL*8 :: tleafn10c                               !10-day mean leaf temperature for night, constrained for physiological range (oC)
      REAL*8 :: Vcmax                                   !the maximum carboxyaltion rate (umol/m2/s) 
      REAL*8 :: vcmx25_opt   
      REAL*8 :: jmx25_opt   
      REAL*8 :: rabsorb
      REAL*8 :: radmax2mean 
      integer :: KcKjFlag                                !flag to indicate whether to update the Kc and Kj using the photosynthesis subroutine; 0--Kc and Kj need to be calculated; 1--Kc and Kj is prescribed.
      integer :: jj                                      !index record fo the number of iterations
      integer :: increase_flag                           !whether to increase or decrease
      
      !-------------------------------------------------------------------------------------------------------------------------------
      !pass thru variables
      REAL*8, intent (in) :: ftg0                       !the intercept of the stomatal conductance equation
      REAL*8, intent (in) :: ftg1                       !the slope of the stomatal conductance equation 
      integer, intent (in):: gs_func                    !flag to indicate which stomatal conductance function to use 
      
      !-------------------------------------------------------------------------------------------------------------------------------
      ! functions
      REAL*8 :: RespTBernacchi
      REAL*8 :: t_scalar                                !temperature scaling factor for Vcmax or Jmax 

      !------------------------------------------------------------------------------ 
      !Constants  
      REAL*8, parameter :: rhol = 0.075                     ! leaf reflectance: 1=vis, 2=nir  
      REAL*8, parameter :: taul = 0.075                     ! leaf transmittance: 1=vis, 2=nir 
      REAL*8, parameter :: Cb   = 1.78d0                    ! nitrogen use effiency for choloraphyll for light capture, see Evans 1989  
      REAL*8, parameter :: Cv   = 1.2d-50 * 3600.d0         ! conversion factor from umol CO2 to g carbon
      REAL*8, parameter :: Kc25 = 40.49d0                   ! Mechalis constant of CO2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      REAL*8, parameter :: Ko25 = 27840d0                   ! Mechalis constant of O2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      REAL*8, parameter :: Cp25 = 4.275d0                   ! CO2 compensation point at 25C (Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      REAL*8, parameter :: Fc25 = 294.2d0                   ! Fc25 = 6.22*47.3 #see Rogers (2014) Photosynthesis Research 
      REAL*8, parameter :: Fj25 = 1257.0d0                  ! Fj25 = 8.06*156 # #see COSTE 2005 and Xu et al 2012
      REAL*8, parameter :: NUEr25 = 33.69d0                 ! nitrogen use efficiency for respiration, see Xu et al 2012
      REAL*8, parameter :: O2ref = 209460.0d0               ! ppm of O2 in the air
      REAL*8, parameter :: CO2ref = 380.0d0                 ! reference CO2 concentration for calculation of reference NUE. 
      REAL*8, parameter :: forc_pbot_ref = 101325.0d0       ! reference air pressure for calculation of reference NUE
      REAL*8, parameter :: Jmaxb0 = 0.0311d0                ! the baseline proportion of nitrogen allocated for electron transport (J)     
      REAL*8, parameter :: Jmaxb1 = 0.1745d0                ! the baseline proportion of nitrogen allocated for electron transport (J)    
      REAL*8, parameter :: Wc2Wjb0 = 0.8054d0               ! the baseline ratio of rubisco limited rate vs light limited photosynthetic rate (Wc:Wj) 
      REAL*8, parameter :: relhExp = 6.0999d0               ! electron transport parameters related to relative humidity
      REAL*8, parameter :: NMCp25 = 0.715d0                 ! estimated by assuming 80% maintenance respiration is used for photosynthesis enzyme maintenance
      REAL*8, parameter :: Trange1 = 5.0d0                  ! lower temperature limit (oC) for nitrogen optimization  
      REAL*8, parameter :: Trange2 = 42.0d0                 ! upper temperature limit (oC) for nitrogen optimization
      REAL*8, parameter :: SNC = 0.004d0                    ! structural nitrogen concentration (g N g-1 dry mass carbon)
      REAL*8, parameter :: mp = 9.0d0                       ! slope of stomatal conductance; this is used to estimate model parameter, but may need to be updated from the physiology file, 
      !SDGVM PAR units are in molm-2s-1
      REAL*8, parameter :: PARLowLim = 200.0d-6             ! minimum photosynthetically active radiation for nitrogen optimization
      REAL*8, parameter :: minrelh = 0.25d0                 ! minimum relative humdity for nitrogen optimization
     
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
   
      !print*, z,vcmx25_z,jmx25_z,PNlc_z,enzs_z,lnc,par240d
   
      ! SDGVM asssumes daily average temperatures for the leaf & constant pressure 
      tleafd10    = tair10
      tleafn10    = tair10
      forc_pbot10 = 101325.d0
 
      !SDGVM PAR units are in molm-2s-1
      !par240d_z   = par240d * 1d6
      !par240x_z   = par240x * 1d6
      par240d_z   = par240d * 1.d0
      par240x_z   = par240x * 1.d0
      
      ! SLA is fixed through the canopy in SDGVM and lnc is pre-determined according to Beer's Law scaling or similar 
      ! - therefore relCLNCa, PARTop are not needed
      ! PAR is already passed to this routine in umol m-2 s-1 (check not day)
      ! - therefore the intermediary steps to convert from wm-2 are not necessary

      !------------------------------------------------------------------
      !SNCa     =  1.0d0/slatop(ft) * SNC !(Eq A3)
      SNCa = 1.0d0/sla/0.48d0 * SNC     !(Eq A3)
      FNCa = lnc - SNCa                !(Eq A4)
      
      !------------------------------------------------------------------
      ! for SDGVM the distinction between sunlit and non-sunlit leaves is not made in the LUNA model
      ! - this is because during the day leaves experience both sunlit and non-sunlit states so  
      ! par240d/x_z is 10-day mean/maximum PAR for leaves in a canopy layer
      ! rabsorb is ratio of absorbed light to incident light (rhol and taul are reflectance and transmittance)
      rabsorb     = 1.0d0-rhol-taul
      radmax2mean = par240x_z / par240d_z
      PARi10      = par240d_z / rabsorb
      PARimx10    = PARi10 * radmax2mean
      
      !------------------------------------------------------------------
      !nitrogen allocation model-start          
      PNlcold     = PNlc_z
      PNetold     = 0.0d0
      PNrespold   = 0.0d0
      PNcbold     = 0.0d0                                     
      !end brought in from above subroutine 
      
      ! this call to N allocation has been replaced by the N allocation code 
      !call NitrogenAllocation(FNCa,forc_pbot10(p), relh10, CO2a10, O2a10, PARi10, PARimx10, rb10, hourpd, &
      !     tair10, tleafd10, tleafn10, &
      !     Jmaxb0, Jmaxb1, Wc2Wjb0, relhExp,  PNstoreold, PNlcold, PNetold, PNrespold, &
      !     PNcbold, PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt)
      
      ! set referfence NUE parameters 
      ! - using sinlge common NUE function fed with default values as arguments
      ! call NUEref(NUEjref, NUEcref, Kj2Kcref)
      !print*, 'LUNA start'
      call NUE(O2ref*0.1013d0, 0.7d0*CO2ref*0.1013d0, 25.d0, 25.d0,
     &NUEjref, NUEcref, Kj2Kcref,ttype,ftToptV,ftHaV,ftHdV,
     &ftToptJ,ftHaJ,ftHdJ)
      
      ! appears theta_cj is not necessary for this subroutine
      !theta_cj = 0.95d0
      ! it's not currently clear to me where PNstoreold is passed from in the original code
      !Nstore   = PNstoreold * FNCa                         !FNCa * proportion of storage nitrogen in functional nitrogen
      Nlc      = PNlcold    * FNCa                         !FNCa * proportion of light capturing nitrogen in functional nitrogen
      Net      = PNetold    * FNCa                         !FNCa * proportion of light harvesting (electron transport) nitrogen in functional nitrogen
      Nresp    = PNrespold  * FNCa                         !FNCa * proportion of respiration nitrogen in functional nitrogen
      Ncb      = PNcbold    * FNCa                         !FNCa * proportion of carboxylation nitrogen in functional nitrogen
      if (Nlc > FNCa * 0.5d0) Nlc = 0.5d0 * FNCa
      PNlc     = PNlcold
      PNlcoldi = PNlcold  - 0.001d0
      
      ! constrain the physiological range of PAR and t
      PARi10c   = max(PARLowLim, PARi10)
      PARimx10c = max(PARLowLim, PARimx10)
      tleafd10c = min(max(tleafd10, Trange1), Trange2)  
      tleafn10c = min(max(tleafn10, Trange1), Trange2) 
      
      ! initialse solver parameters 
      ci            = 0.7d0 * CO2a10 
      JmaxCoef      = Jmaxb1 * ((hourpd / 12.0d0)**2.0d0) * 
     &( 1.0d0 - exp( -relhExp * max(relh10 - minrelh, 0.0d0) / 
     & (1.0d0 - minrelh) ) )
      chg_per_step  = 0.02* FNCa
      increase_flag = 0
      
      !------------------------------------------------------------------
      jj = 1
      do while ( (PNlcoldi .NE. PNlc) .and. (jj < 100) )      
     
       ! Fc is the scaling factor to go from leaf N invested in RuBisCO to Vcmax  
       ! Fj is the scaling factor to go from leaf N invested in elec trans to Jmax  
       !Fc   = VcmxTKattge(tair10, tleafd10c) * Fc25
       !Fj   = JmxTKattge(tair10, tleafd10c)  * Fj25
       Fc=T_SCALAR(tleafd10c,'v',ttype,ftToptV,ftHaV,ftHdV,tair10)*Fc25
       Fj=T_SCALAR(tleafd10c,'j',ttype,ftToptJ,ftHaJ,ftHdJ,tair10)*Fj25
       
       NUEr = Cv * NUEr25 * ( RespTBernacchi(tleafd10c) * hourpd +
     &RespTBernacchi(tleafn10c) * (24.0d0 - hourpd) ) !nitrogen use efficiency for respiration (g biomass/m2/day/g N)
      
       !print*, 'LUNA do while loop'
 
       call NUE(O2a10, ci, tair10, tleafd10c, NUEj, NUEc, Kj2Kc,ttype,
     &ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ)
       
       KcKjFlag = 0
       call LUNA_Ninvestments (KcKjFlag,FNCa, Nlc, forc_pbot10,relh10,
     &CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, tair10, 
     &tleafd10c,tleafn10c, 
     &Kj2Kc, Wc2Wjb0, JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref,NUEjref,NUEr,
     &Kc, Kj, ci, Vcmax, Jmax,JmeanL,JmaxL, Net, Ncb, Nresp, PSN, RESP,
     &gs_func,ftg0,ftg1)      

       !target nitrogen allocated to photosynthesis, which may be lower or higher than Npsn_avail
       Npsntarget = Nlc + Ncb + Net       
       PNlcoldi   = Nlc / FNCa
       Nstore     = FNCa - Npsntarget - Nresp
       
       !test the increase of light capture nitrogen
       if ( ((Nstore > 0.0d0) .and.
     &(increase_flag .eq. 1)) .or. (jj .eq. 1) ) then
        Nlc2 = Nlc + chg_per_step
        if (Nlc2 / FNCa > 0.95d0) Nlc2 = 0.95d0 * FNCa
        
        KcKjFlag = 1
        call LUNA_Ninvestments (KcKjFlag,FNCa, Nlc2, forc_pbot10, 
     &relh10, CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, 
     &tair10, tleafd10c,tleafn10c, 
     &Kj2Kc, Wc2Wjb0, JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref,NUEjref,NUEr,
     &Kc,Kj,ci,Vcmax,Jmax,JmeanL,JmaxL, Net2, Ncb2, Nresp2, PSN2, RESP2,
     &gs_func,ftg0,ftg1)
        
        Npsntarget2 = Nlc2 + Ncb2 + Net2
        
        !update the nitrogen change
        Carboncost2 = (Npsntarget2 - Npsntarget) * NMCp25 * Cv * 
     &( RespTBernacchi(tleafd10c) * hourpd + RespTBernacchi(tleafn10c) *
     &(24.0d0  - hourpd) )
        Carbongain2 =  PSN2 - PSN
        
        if( (Carbongain2 > Carboncost2) .and. 
     &(Npsntarget2 + Nresp2 < 0.95d0 * FNCa) ) then
         Nlc = Nlc2
         Net = Net2
         Ncb = Ncb2
         Nstore = FNCa - Npsntarget2 - Nresp2 
         if (jj == 1) increase_flag = 1
        end if
       end if
       
       !------------------------------------------------------------------------------------
       !test the decrease of light capture nitrogen
       if (increase_flag == 0) then  
       
       if (Nstore < 0.0d0) then
        Nlc1 = Nlc * 0.8d0 !bigger step of decrease if it is negative            
       else
        Nlc1 = Nlc - chg_per_step
       end if
       
       if (Nlc1 < 0.05d0) Nlc1 = 0.05d0
        KcKjFlag = 1
        call LUNA_Ninvestments(KcKjFlag,FNCa, Nlc1,forc_pbot10,relh10,
     &CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, 
     &tair10, tleafd10c,tleafn10c, 
     &Kj2Kc, Wc2Wjb0,JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref,NUEjref,NUEr,
     &Kc,Kj,ci,Vcmax,Jmax,JmeanL,JmaxL,Net1, Ncb1, Nresp1, PSN1, RESP1,
     &gs_func,ftg0,ftg1)
        
        Npsntarget1 = Nlc1 + Ncb1 + Net1
        Carboncost1 = (Npsntarget - Npsntarget1) * NMCp25 * Cv *
     &( RespTBernacchi(tleafd10c) * hourpd + RespTBernacchi(tleafn10c) *
     &(24.0d0  - hourpd) )
        Carbongain1 =  PSN - PSN1
        
        if( ((Carbongain1 < Carboncost1) .and. (Nlc1 > 0.05d0)) .or.
     &(Npsntarget + Nresp > 0.95d0 * FNCa) ) then
         Nlc    = Nlc1 
         Net    = Net1   
         Ncb    = Ncb1
         Nstore = FNCa - Npsntarget1 - Nresp1  
        end if
        
       end if
       PNlc = Nlc / FNCa
       
       jj = jj + 1  
      end do                        
      !------------------------------------------------------------------

      ! record optimum proportions      
      PNlcopt    = Nlc    / FNCa
      PNstoreopt = Nstore / FNCa
      PNcbopt    = Ncb    / FNCa
      PNetopt    = Net    / FNCa
      PNrespopt  = Nresp  / FNCa 
      
      !brought in from parent subroutine
      PNlc_z      = PNlcopt

      ! determine change in vcmax and jmax 
      vcmx25_opt  = PNcbopt * FNCa * Fc25
      jmx25_opt   = PNetopt * FNCa * Fj25
      
      chg         = vcmx25_opt-vcmx25_z
      chg_constrn = min(abs(chg),vcmx25_z*max_daily_pchg)
      vcmx25_z    = vcmx25_z+sign(1.0d0,chg)*chg_constrn
      !if(z.eq.1) print*, vcmx25_opt,chg,vcmx25_z      
      !print*, vcmx25_opt,chg,vcmx25_z      
 
      chg         = jmx25_opt-jmx25_z
      chg_constrn = min(abs(chg),jmx25_z*max_daily_pchg)
      jmx25_z     = jmx25_z+sign(1.0d0,chg)*chg_constrn 
      !print*, jmx25_opt,chg,jmx25_z      
      
      if(enzs_z<1.0) enzs_z = enzs_z * (1.0d0 + max_daily_pchg)

      !nitrogen allocation subroutine end  
      !------------------------------------------------------------------

      if( isnan(vcmx25_z) .or. (vcmx25_z>1000.d0) .or. 
     &(vcmx25_z<0.d0) ) then
       !write(iulog, *) 'Error: Vc,mx25 become unrealistic (NaN,>1000,
       write(*, *) 'Error: Vc,mx25 become unrealistic (NaN,>1000,
     & or negative) for z=', z
       write(*, *) 'Error: Vcx25:', vcmx25_z
       write(*, *) 'Error: Jmx25:', jmx25_z
       !write(iulog, *) 'LUNA env:',FNCa,forc_pbot10, relh10, CO2a10, 
       write(*, *) 'LUNA env:',FNCa,forc_pbot10, relh10, CO2a10, 
     &O2a10, PARi10, PARimx10, rb10, 
     &hourpd, tair10, tleafd10, tleafn10
       !call endrun(msg=errmsg(__FILE__, __LINE__))
       stop
      endif
      if( isnan(jmx25_z) .or. (jmx25_z>1000.d0) .or. 
     &(jmx25_z<0.d0) ) then
       !write(iulog, *) 'Error: Jmx25 become unrealistic (NaN,>1000, 
       write(*, *) 'Error: Jmx25 become unrealistic (NaN,>1000, 
     &or negative)for z=', z
       write(*, *) 'Error: Jmx25:', jmx25_z
       write(*, *) 'Error: Vcx25:', vcmx25_z
       !write(iulog, *) 'LUNA env:', FNCa,forc_pbot10, relh10, CO2a10,
       write(*, *) 'LUNA env:', FNCa,forc_pbot10, relh10, CO2a10,
     &O2a10, PARi10, PARimx10, rb10,hourpd, tair10, tleafd10, tleafn10
       !call endrun(msg=errmsg(__FILE__, __LINE__))
       stop
      endif
      !end brought in from above subroutine 
      
      end subroutine LUNA
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      subroutine LUNA_Ninvestments (KcKjFlag, FNCa, Nlc, forc_pbot10, 
     &relh10,CO2a10,O2a10,PARi10,PARimx10,rb10,hourpd,tair10,tleafd10,
     &tleafn10,Kj2Kc, Wc2Wjb0, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, 
     &NUEjref,NUEr,Kc,Kj,ci,Vcmax,Jmax,JmeanL,JmaxL,Net,Ncb,Nresp,PSN,
     &RESP,gs_func,ftg0,ftg1)
      
      !calculate the nitrogen investment for electron transport, carb10oxylation, respiration given a specified value 
      !of nitrogen allocation in light capture [Nlc]. This equation are based on Ali et al 2015b.
      
      implicit none
      
      integer,intent (in) :: KcKjFlag                   !flag to indicate whether to update the Kc and Kj using the photosynthesis subroutine; 0--Kc and Kj need to be calculated; 1--Kc and Kj is prescribed.
      REAL*8, intent (in) :: FNCa                       !Area based functional nitrogen content (g N/m2 leaf)
      REAL*8, intent (in) :: Nlc                        !nitrogen content for light capture(g N/m2 leaf)
      REAL*8, intent (in) :: forc_pbot10                !10-day mean air pressure (Pa)
      REAL*8, intent (in) :: relh10                     !10-day mean relative humidity (unitless)
      REAL*8, intent (in) :: CO2a10                     !10-day mean CO2 concentration in the air (Pa)
      REAL*8, intent (in) :: O2a10                      !10-day mean O2 concentration in the air (Pa)
      REAL*8, intent (in) :: PARi10                     !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
      REAL*8, intent (in) :: PARimx10                   !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
      REAL*8, intent (in) :: rb10                       !10-day mean boundary layer resistance (s/m)
      REAL*8, intent (in) :: hourpd                     !hours of light in a the day (hrs)
      REAL*8, intent (in) :: tair10                     !10-day running mean of the 2m temperature (oC)
      REAL*8, intent (in) :: tleafd10                   !10-day mean daytime leaf temperature (oC) 
      REAL*8, intent (in) :: tleafn10                   !10-day mean nighttime leaf temperature (oC) 
      REAL*8, intent (in) :: Kj2Kc                      !ratio:  Kj / Kc
      REAL*8, intent (in) :: Wc2Wjb0                    !the baseline ratio of rubisco-limited rate vs light-limited photosynthetic rate (Wc:Wj)
      REAL*8, intent (in) :: JmaxCoef                   !coefficient determining the response of electron transport rate to light availability (unitless) and humidity
      REAL*8, intent (in) :: Fc                         !the temperature adjusted ratio of Vcmax to N invested in Vcmax  
      REAL*8, intent (in) :: Fj                         !the temperature adjusted ratio of  Jmax to N invested in electron trans 
      REAL*8, intent (in) :: NUEc                       !nitrogen use efficiency for carboxylation 
      REAL*8, intent (in) :: NUEj                       !nitrogen use efficiency for electron transport
      REAL*8, intent (in) :: NUEcref                    !nitrogen use efficiency for carboxylation under reference climates
      REAL*8, intent (in) :: NUEjref                    !nitrogen use efficiency for electron transport under reference climates
      REAL*8, intent (in) :: NUEr                       !nitrogen use efficiency for respiration
      REAL*8, intent (inout) :: Kc                      !conversion factors from Vc,max to Wc 
      REAL*8, intent (inout) :: Kj                      !conversion factor from electron transport rate to Wj 
      REAL*8, intent (inout) :: ci                      !inter-cellular CO2 concentration (Pa) 
      REAL*8, intent (out)   :: Vcmax                   !the maximum carboxyaltion rate (umol/m2/s) 
      REAL*8, intent (out)   :: Jmax                    !the maximum electron transport rate (umol/m2/s) 
      REAL*8, intent (out)   :: JmaxL                   !the electron transport rate with maximum daily radiation (umol/m2/s)  
      REAL*8, intent (out)   :: JmeanL                  !the electron transport rate with mean radiation (umol/m2/s) 
      REAL*8, intent (out)   :: Net                     !nitrogen content for electron transport(g N/m2 leaf)
      REAL*8, intent (out)   :: Ncb                     !nitrogen content for carboxylation(g N/m2 leaf)
      REAL*8, intent (out)   :: Nresp                   !nitrogen content for respiration(g N/m2 leaf)
      REAL*8, intent (out)   :: PSN                     !daily photosynthetic rate(g C/day/m2 leaf)
      REAL*8, intent (out)   :: RESP                    !daily respiration rate(g C/day/m2 leaf)
      
      !-------------------------------------------------------------------------------------------------------------------------------
      !intermediate variables
      REAL*8 :: A                                       !Gross photosynthetic rate (umol CO2/m2/s)
      REAL*8 :: Wc2Wj                                   !ratio: Wc/Wj  
      REAL*8 :: ELTRNabsorb                             !absorbed electron rate, umol electron/m2 leaf /s
      REAL*8 :: Jmaxb0act                               !base value of Jmax (umol/m2/s) 
      REAL*8 :: theta_cj                                !interpolation coefficient
      REAL*8 :: theta                                   !light absorption rate (0-1)
      REAL*8 :: Vcmaxnight                              !Vcmax during night (umol/m2/s)
      REAL*8 :: Wc                                      !rubisco-limited photosynthetic rate (umol/m2/s)
      REAL*8 :: Wj                                      !light limited photosynthetic rate (umol/m2/s)
      REAL*8 :: k_c      
      REAL*8 :: k_o       
      REAL*8 :: tau        
      REAL*8 :: c_p        
      REAL*8 :: awc        
      REAL*8 :: gs        
      REAL*8 :: NUECHG                                  !the nitrogen use efficiency change under current conidtions compared to reference climate conditions (25oC and 385 ppm )
      
      !-------------------------------------------------------------------------------------------------------------------------------
      !pass thru variables
      REAL*8, intent (in) :: ftg0                       !the intercept of the stomatal conductance equation
      REAL*8, intent (in) :: ftg1                       !the slope of the stomatal conductance equation 
      integer,intent (in) :: gs_func                    !flag to indicate which stomatal conductance function to use 
      
      !-------------------------------------------------------------------------------------------------------------------------------
      !parameters
      REAL*8, parameter :: Cb = 1.78d0                      ! nitrogen use effiency for choloraphyll for light capture, see Evans 1989  
      REAL*8, parameter :: Cv = 1.2d-50 * 3600.d0           ! conversion factor from umol CO2 to g carbon
      REAL*8, parameter :: Jmaxb0 = 0.0311d0                ! the baseline proportion of nitrogen allocated for electron transport (J)     
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      theta_cj    = 0.95d0
      theta       = 0.292d0 / (1.0d0 + 0.076d0 / (Nlc * Cb))
      ELTRNabsorb = theta  * PARi10
      Jmaxb0act   = Jmaxb0 * FNCa * Fj
      ! adjusted for SDGVM units molm-2s-1
      Jmaxb0act   = Jmaxb0act * 1d-6
      Jmax        = Jmaxb0act + JmaxCoef * ELTRNabsorb
      JmaxL       = theta  * PARimx10 / ( sqrt(1.0d0 + 
     &(theta * PARimx10 / Jmax)**2.0d0) )        
      NUEchg      = (NUEc / NUEcref) * (NUEjref / NUEj)
      Wc2Wj       = Wc2Wjb0 * (NUEchg**0.5d0)
      Wc2Wj       = min(1.0d0, Wc2Wj)
      Vcmax       = Wc2Wj * JmaxL * Kj2Kc
     
      !added the two below lines for stability when not calculating LUNA everyday
      ! - in the lowest canopy layers the solver sometimes returns negative values 
      ! - and with larger intervals between LUNA calcultions (~10 days) the max proportional change can equal 1
      ! - allowing returned Vcmax and Jmax to be below 0
      Vcmax       = max(1d-7,Vcmax)
      Jmax        = max(1d-7,Jmax)
     
      JmeanL      = theta * PARi10 / ( sqrt(1.0d0 +
     &(ELTRNabsorb / Jmax)**2.0d0) )
 
      if(KcKjFlag.eq.0) then      !update the Kc,Kj, anc ci information
       
       ! From SDGVM - for consistency
       ! in LUNA ko, kc, and tau are at the 10-day mean leaf temp
       ! - need to check that units are consistent, in SDGVM these are all in Pa
       k_c  = exp(35.8d0   - 80.5d0 / (0.00831d0*(tleafd10+273.15)) )
       k_o  = exp(9.6d0    - 14.51d0/ (0.00831d0*(tleafd10+273.15)) )
       k_o  = k_o * 1.0d3
       tau  = exp(-3.949d0 + 28.99d0/ (0.00831d0*(tleafd10+273.15)) )
       
       ! this function is over-parameterised due to an old attempt to use the Brent solver that is currently not implemented
       ! Vcmax and J are also assumed to be in molm-2s-1 in SDGVM
       ! ga and gs must be in molm-2s-1 - gs is a local variable calculated during the call to assimilation_calc
       ! rd needs to be unpacked, ko, kc, and tau are also needed here 
       ! z is the LAI layer that is currently under calculation - in SDGVM this only triggers print to screen for LAI layer 1        
       ! call Photosynthesis_luna(forc_pbot10, tleafd10, relh10, CO2a10, O2a10,rb10, Vcmax, JmeanL, ci, Kc, Kj, A) 
       ! CALL ASSIMILATION_CALC(fshade,fsunlit,
       !&vmx,xvmax,rd,jshade,jsunlit,upt,qshade,
       !&qsunlit,t,rn,soil2g,wtwp,ga,rh,C3,kg,ko,kc,tau,p,oi,ca,
       !&a,gs,ci,day,mnth,oday,omnth,.FALSE.,gs_func,ftg0,ftg1,i)
 
       ! should rd be zero below? seems like maybe it should according to the code but by not including rd in the a,ci,gs solution this will underestimate ci
       ! the below function calls the SDGVM calculation of assimilation, simultaneously solving A, gs, & ci
       ! for the purposes of the LUNA model A (umolm-2s-1)  and ci (Pa) are returned and used in further calculations
       ! inputs to this model are in Pa units and molm-2s-1 
       CALL ASSIMILATION_CALC(0.0d0,1.0d0,
     &Vcmax,0.0d0,0.0d0,0.0d0,JmeanL,0.0d0,0.0d0,
     &0.0d0,tleafd10,1.0d0,1.0d0,0.0d0,rb10,relh10,1,1.0d0,k_o,k_c,tau,
     &forc_pbot10,O2a10,CO2a10,A,gs,ci,1,2,1,1,.FALSE.,gs_func,ftg0,
     &ftg1,2)

       c_p = 0.5d0*O2a10/tau
       awc = k_c * (1.0d0 + O2a10 / k_o)
       Kj  = max(ci - c_p, 0.0d0) / (4.0d0 * ci + 8.0d0 * c_p)
       Kc  = max(ci - c_p, 0.0d0) / (ci + awc)
       
      else
       
       Wc = Kc * Vcmax
       Wj = Kj * JmeanL
       !A  = (1.0d0 - theta_cj) * max(Wc, Wj) + theta_cj * min(Wc, Wj) 
       ! Vcmax and J are in molm-2s-1 because that's what SDGVM expects, so the above calculation gives Kc and Kj in molm-2s-1
       A  = ( (1.0d0 - theta_cj) * max(Wc, Wj) + theta_cj * min(Wc, Wj) 
     &) * 1e6
       
      endif
      
      ! Cv converts from umolm-2s-1 to gC m-2 hour-1
      PSN        = Cv * A * hourpd
!      Vcmaxnight = VcmxTKattge(tair10, tleafn10) / 
!     &VcmxTKattge(tair10, tleafd10) * Vcmax
      ! SDGVM assumes mean 24 hr air temp 
      Vcmaxnight = Vcmax
      ! again these calculations are adjusted to suit SDGVM molm-2s-1 units for Vcmax and Jmax 
      RESP       = Cv * 0.015d0 * ( Vcmax * hourpd + Vcmaxnight * 
     &(24.d0 - hourpd) ) * 1e6                
      Net        = Jmax  * 1e6 / Fj
      Ncb        = Vcmax * 1e6 / Fc
      Nresp      = RESP  / NUEr
     
      !print*, A,RESP,Vcmax,Jmax
 
      end subroutine LUNA_Ninvestments
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      !Calculate the Nitrogen use effieciency dependence on CO2 and leaf temperature
      subroutine NUE(o2a,ci,tgrow,tleaf,NUEj,NUEc,Kj2Kc,
     &ttype,ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ)
      
      ! this function can be used to calculate NUEref aswell, it just needs to be called with both temps at 25oC and a ci that maintains a ci:ca of 0.7
      
      ! uses module constants of Fc25, Fj25, tfrz, rgas
      
      implicit none
      REAL*8, intent (in) :: o2a                        !air O2 partial presuure (Pa)
      REAL*8, intent (in) :: ci                         !leaf inter-cellular [CO2] (Pa) - originally incorrectly labelled as PPM
      REAL*8, intent (in) :: tgrow                      !10 day growth temperature (oC), 24 hour mean temperature
      REAL*8, intent (in) :: tleaf                      !leaf temperature (oC)
      REAL*8, intent (in) :: ftToptV                    !optimum temp for vcmax (oC)
      REAL*8, intent (in) :: ftHaV                      !activation energy for vcmax
      REAL*8, intent (in) :: ftHdV                      !deactivation energy for vcmax
      REAL*8, intent (in) :: ftToptJ                    !optimum temp for jmax (oC)
      REAL*8, intent (in) :: ftHaJ                      !activation energy for jmax
      REAL*8, intent (in) :: ftHdJ                      !deactivation energy for jmax
      !integer, intent (in):: ttype                      !temperature scaling method to use
      integer:: ttype                      !temperature scaling method to use
      REAL*8, intent (out):: NUEj                       !nitrogen use efficiency for electron transport under refernce environmental conditions (25oC and 385 ppm co2)
      REAL*8, intent (out):: NUEc                       !nitrogen use efficiency for carboxylation under reference environmental conditions  (25oC and 385 ppm co2)
      REAL*8, intent (out):: Kj2Kc                      !the ratio of Kj to Kc 
      !------------------------------------------------
      !intermediate variables
      REAL*8 :: Fj                                      !the temperature adjusted factor for Jmax 
      REAL*8 :: Fc                                      !the temperature adjusted factor for Vcmax 
      REAL*8 :: Kc                                      !conversion factor from Vcmax to Wc 
      REAL*8 :: Kj                                      !conversion factor from J to W 
      REAL*8 :: k_o                                     !Rubsico O2 specifity
      REAL*8 :: k_c                                     !Rubsico CO2 specifity
      REAL*8 :: awc                                     !second deminator term for rubsico limited carboxylation rate based on Farquhar model
      REAL*8 :: c_p                                     !CO2 compenstation point (Pa)
      REAL*8 :: tau
      REAL*8 :: t_scalar                                !temperature scaling factor for Vcmax or Jmax 
      
      REAL*8, parameter :: Fc25 = 294.2d0              ! Fc25 = 6.22*47.3 #see Rogers (2014) Photosynthesis Research 
      REAL*8, parameter :: Fj25 = 1257.0d0             ! Fj25 = 8.06*156  #see COSTE 2005 and Xu et al 2012
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      !print*, 'LUNA NUE, ttype:', ttype
      !Fc  = VcmxTKattge(tgrow, tleaf) * Fc25
      !Fj  = JmxTKattge(tgrow, tleaf)  * Fj25
      Fc  =T_SCALAR(tleaf,'v',ttype,ftToptV,ftHaV,ftHdV,tgrow)*Fc25
      Fj  =T_SCALAR(tleaf,'j',ttype,ftToptJ,ftHaJ,ftHdJ,tgrow)*Fj25
c      k_c = Kc25 * exp((79430.0d0 / (rgas*1.e-3d0 * 
c     &(25.0d0 + tfrz))) * (1.0d0 - (tfrz + 25.0d0) / (tfrz + tleaf)))
c      k_o = Ko25 * exp((36380.0d0 / (rgas*1.e-3d0 *
c     &(25.0d0 + tfrz))) * (1.0d0 - (tfrz + 25.0d0) / (tfrz + tleaf)))
c      c_p = Cp25 * exp((37830.0d0 / (rgas*1.e-3d0 * 
c     &(25.0d0 + tfrz))) * (1.0d0 - (tfrz + 25.0d0) / (tfrz + tleaf)))
      
      ! From SDGVM - for consistency
      ! in LUNA ko, kc, and tau are at the 10-day mean leaf temp
      ! in SDGVM these are all in Pa
      k_c = exp(35.8d0   - 80.5d0 /(0.00831d0*(tleaf+273.15)))
      k_o = exp(9.6d0    - 14.51d0/(0.00831d0*(tleaf+273.15)))
      k_o = k_o*1.d3
      tau = exp(-3.949d0 + 28.99d0/(0.00831d0*(tleaf+273.15)))
      c_p = 0.5d0*o2a/tau
     
      awc = k_c * ( 1.0d0 + o2a/k_o )
      Kj  = max( ci-c_p,0.0d0 ) / ( 4.0d0*ci + 8.0d0*c_p )
      Kc  = max( ci-c_p,0.0d0 ) / ( ci+awc )
     
      NUEj  = Kj * Fj
      NUEc  = Kc * Fc  
      Kj2Kc = Kj / Kc
     
      end subroutine NUE
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      !Calculate the temperature response for respiration, following Bernacchi PCE 2001
      function  RespTBernacchi(tleaf)
      
      implicit none
      
      REAL*8 RespTBernacchi
      REAL*8, intent(in):: tleaf  !leaf temperature (oC)
      
      !RespTBernacchi= exp(18.72d0-46.39d0/(rgas*1.e-6d0 *
      !&(tleaf+tfrz)))
      RespTBernacchi= exp( 18.72d0-46.39d0/(8.31d-3 *
     &(tleaf+273.15)) )
      
      end function RespTBernacchi
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
     
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      function max_daily_pchg(tleaf10) ! - define this as a function to return max_daily_pchg  
      
      ! !LOCAL VARIABLES:
      !
      ! local pointers to implicit in variables
      
      REAL*8               :: max_daily_pchg                  ! maximum daily percentrage change  for nitrogen allocation
      REAL*8, intent (in)  :: tleaf10                         ! 10-day running mean of leaf temperature (oC)
      REAL*8, parameter    :: Q10Enz = 2.0d0                 ! Q10 value for enzyme decay rate
      REAL*8, parameter    :: Enzyme_turnover_daily = 0.1d0  ! the daily turnover rate for photosynthetic enzyme at 25oC in view of ~7 days of half-life time for Rubisco (Suzuki et al. 2001)
      REAL*8               :: EnzTurnoverTFactor              ! temperature adjust factor for enzyme decay
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      !calculate the enzyme ternover rate
      EnzTurnoverTFactor = Q10Enz**(0.1d0*(min(40.0d0, tleaf10)- 25.d0))
      max_daily_pchg     = EnzTurnoverTFactor * Enzyme_turnover_daily
      
      end function 
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
     
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      subroutine LUNA_nogrowth(max_daily_pchg,vcmx25_z,jmx25_z,enzs_z) 
      
      REAL*8, intent (in)    :: max_daily_pchg      ! maximum daily percentrage change for nitrogen allocation
      REAL*8, intent (inout) :: vcmx25_z            ! leaf Vc,max25 (umol/m2 leaf/s) for canopy layer 
      REAL*8, intent (inout) :: jmx25_z             ! leaf Jmax25 (umol electron/m**2/s) for canopy layer
      REAL*8, intent (inout) :: enzs_z              ! enzyme decay status 1.0-fully active; 0-all decayed during stress
      REAL*8                 :: max_daily_decay     ! maximum daily percentrage change for nitrogen allocation
      !INTEGER                :: nrad                ! number of canopy layers
      !INTEGER                :: z                   ! canopy layer index

      ! need to determine if these really need to have each layer or are passed back to each layer in npp calc
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      !nrad            = 12

      !assume enzyme turnover under maintenance is 10 times lower than enzyme change under growth
      max_daily_decay = min(0.5d0, 0.1d0 * max_daily_pchg)

      !print*, 'LUNA nogrowth:', max_daily_decay
      !print*, vcmx25_z(1)

      !do z = 1 , nrad
       !decay is set at only 50% of original enzyme in view that plant will need to maintain their basic functionality
       if(enzs_z>0.5d0) then 
        enzs_z   = enzs_z   * (1.0d0 - max_daily_decay)
        jmx25_z  = jmx25_z  * (1.0d0 - max_daily_decay) 
        vcmx25_z = vcmx25_z * (1.0d0 - max_daily_decay) 
       endif
      !end do              
      
      !print*, vcmx25_z(1)
      
      end subroutine LUNA_nogrowth 
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
