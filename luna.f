      
     !********************************************************************************************************************************************************************** 
     ! this subroutine updates the photosynthetic capacity as determined by Vcmax25 and Jmax25
     
     !subroutine Update_Photosynthesis_Capacity(bounds, fn, filterp, 
     !&dayl_factor, atm2lnd_inst, temperature_inst, canopystate_inst, 
     !&photosyns_inst, 
     !&surfalb_inst, solarabs_inst, waterstate_inst, frictionvel_inst)
     
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      subroutine LUNA_INIT() ! - define this as a function to return max_daily_pchg  
      !
      ! !DESCRIPTION:
      ! Calculates Nitrogen fractionation within the leaf, based on optimum calculated fractions in rubisco, cholorophyll, 
      ! Respiration and Storage. Based on Xu et al. 2012 and Ali et al 2015.In Review 
      
      !
      ! !REVISION HISTORY:
      ! version 1.0, by Chonggang Xu, Ashehad Ali and Rosie Fisher. July 14  2015.
      ! version 0.1, by Chonggang Xu, Ashehad Ali and Rosie Fisher. October 30 2014. 
      
      ! !LOCAL VARIABLES:
      !
      ! local pointers to implicit in variables
      
      real(r8), parameter :: Q10Enz = 2.0_r8                   ! Q10 value for enzyme decay rate
      real(r8), parameter :: Enzyme_turnover_daily = 0.1_r8    ! the daily turnover rate for photosynthetic enzyme at 25oC in view of ~7 days of half-life time for Rubisco (Suzuki et al. 2001)
      real (r8) :: max_daily_pchg                              ! maximum daily percentrage change  for nitrogen allocation
      real (r8) :: EnzTurnoverTFactor                          ! temperature adjust factor for enzyme decay
      real (r8) :: tleaf10                                     ! 10-day running mean of leaf temperature (oC)
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      !calculate the enzyme ternover rate
      EnzTurnoverTFactor = Q10Enz**(0.1_r8*(min(40.0_r8, tleaf10) - 25.0_r8))            
      max_daily_pchg     = EnzTurnoverTFactor * Enzyme_turnover_daily
      
      end subroutine LUNA_INIT 
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
     
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      subroutine LUNA_nogrowth() 
      
      real (r8) :: max_daily_pchg                              ! maximum daily percentrage change  for nitrogen allocation
      real (r8) :: max_daily_decay                             ! maximum daily percentrage change  for nitrogen allocation
      real (r8) :: vcmx25_z      ! Output: [real(r8) (:,:) ] patch leaf Vc,max25 (umol/m2 leaf/s) for canopy layer 
      real (r8) :: jmx25_z       ! Output: [real(r8) (:,:) ] patch leaf Jmax25 (umol electron/m**2/s) for canopy layer
      real (r8) :: pnlc_z        ! Output: [real(r8) (:,:) ] patch proportion of leaf nitrogen allocated for light capture for 
      real (r8) :: enzs_z        ! Output: [real(r8) (:,:) ] enzyme decay status 1.0-fully active; 0-all decayed during stress

      ! need to determine if these really need to have each layer or are passed back to each layer in npp calc
 
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      max_daily_decay = min(0.5_r8, 0.1_r8 * max_daily_pchg)!assume enzyme turnover under maintenance is 10 times lower than enzyme change under growth
      do z = 1 , nrad
       if(enzs_z(z)>0.5_r8) then !decay is set at only 50% of original enzyme in view that plant will need to maintain their basic functionality
        enzs_z(z)   = enzs_z(z)   * (1.0_r8 - max_daily_decay)
        jmx25_z(z)  = jmx25_z(z)  * (1.0_r8 - max_daily_decay) 
        vcmx25_z(z) = vcmx25_z(z) * (1.0_r8 - max_daily_decay) 
       endif
      end do              
      
      end subroutine LUNA_nogrowth 
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
     
     
     
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      !Use the LUNA model to calculate the Nitrogen partioning 
      !subroutine NitrogenAllocation(FNCa,forc_pbot10, relh10, CO2a10,O2a10, PARi10,PARimx10,rb10, hourpd, tair10, tleafd10, tleafn10, &
      !     Jmaxb0, Jmaxb1, Wc2Wjb0, relhExp,&
      !     PNstoreold, PNlcold, PNetold, PNrespold, PNcbold, &
      !     PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt)
      subroutine LUNA() 
      
      implicit none
      real(r8), intent (in) :: FNCa                       !Area based functional nitrogen content (g N/m2 leaf)
      real(r8), intent (in) :: forc_pbot10                !10-day mean air pressure (Pa)
      real(r8), intent (in) :: relh10                     !10-day mean relative humidity (unitless)
      real(r8), intent (in) :: CO2a10                     !10-day meanCO2 concentration in the air (Pa)
      real(r8), intent (in) :: O2a10                      !10-day mean O2 concentration in the air (Pa)
      real(r8), intent (in) :: PARi10                     !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
      real(r8), intent (in) :: PARimx10                   !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
      real(r8), intent (in) :: rb10                       !10-day mean boundary layer resistance
      real(r8), intent (in) :: hourpd                     !hours of light in a the day (hrs)
      real(r8), intent (in) :: tair10                     !10-day running mean of the 2m temperature (oC)
      real(r8), intent (in) :: tleafd10                   !10-day running mean of daytime leaf temperature (oC) 
      real(r8), intent (in) :: tleafn10                   !10-day running mean of nighttime leaf temperature (oC) 
      real(r8), intent (in) :: Jmaxb0                     !baseline proportion of nitrogen allocated for electron transport rate (unitless)
      real(r8), intent (in) :: Jmaxb1                     !coefficient determining the response of electron transport rate to light availability (unitless) 
      real(r8), intent (in) :: Wc2Wjb0                    !the baseline ratio of rubisco-limited rate vs light-limited photosynthetic rate (Wc:Wj)
      real(r8), intent (in) :: relhExp                    !specifies the impact of relative humidity on electron transport rate (unitless)
      real(r8), intent (in) :: PNstoreold                 !old value of the proportion of nitrogen allocated to storage (unitless)
      real(r8), intent (in) :: PNlcold                    !old value of the proportion of nitrogen allocated to light capture (unitless)
      real(r8), intent (in) :: PNetold                    !old value of the proportion of nitrogen allocated to electron transport (unitless)
      real(r8), intent (in) :: PNrespold                  !old value of the proportion of nitrogen allocated to respiration (unitless)
      real(r8), intent (in) :: PNcbold                    !old value of the proportion of nitrogen allocated to carboxylation (unitless)  
      real(r8), intent (out):: PNstoreopt                 !optimal proportion of nitrogen for storage 
      real(r8), intent (out):: PNlcopt                    !optimal proportion of nitrogen for light capture 
      real(r8), intent (out):: PNetopt                    !optimal proportion of nitrogen for electron transport 
      real(r8), intent (out):: PNrespopt                  !optimal proportion of nitrogen for respiration 
      real(r8), intent (out):: PNcbopt                    !optial proportion of nitrogen for carboxyaltion  
      
      !-------------------------------------------------------------------------------------------------------------------------------
      !intermediate variables
      real(r8) :: Carboncost1                             !absolute amount of carbon cost associated with maintenance respiration due to deccrease in light capture nitrogen(g dry mass per day) 
      real(r8) :: Carboncost2                             !absolute amount of carbon cost associated with maintenance respiration due to increase in light capture nitrogen(g dry mass per day) 
      real(r8) :: Carbongain1                             !absolute amount of carbon gain associated with maintenance respiration due to deccrease in light capture nitrogen(g dry mass per day) 
      real(r8) :: Carbongain2                             !absolute amount of carbon gain associated with maintenance respiration due to increase in light capture nitrogen(g dry mass per day) 
      real(r8) :: Fc                                      !the temperature adjustment factor for Vcmax 
      real(r8) :: Fj                                      !the temperature adjustment factor for Jmax 
      real(r8) :: PNlc                                    !the current nitrogen allocation proportion for light capture
      real(r8) :: Jmax                                    !the maximum electron transport rate (umol/m2/s) 
      real(r8) :: JmaxCoef                                !coefficient determining the response of electron transport rate to light availability (unitless) and humidity
      real(r8) :: Jmaxb0act                               !base value of Jmax (umol/m2/s) 
      real(r8) :: JmaxL                                   !the electron transport rate with maximum daily radiation (umol/m2/s)  
      real(r8) :: JmeanL                                  !the electron transport rate with mean radiation (umol/m2/s) 
      real(r8) :: Nstore                                  !absolute amount of nitrogen allocated to storage (gN/m2 leaf)
      real(r8) :: Nresp                                   !absolute amount of nitrogen allocated to respiration (gN/m2 leaf) 
      real(r8) :: Nlc                                     !absolute amount of nitrogen allocated to light capture (gN/m2 leaf) 
      real(r8) :: Net                                     !absolute amount of nitrogen allocated to electron transport (gN/m2 leaf) 
      real(r8) :: Ncb                                     !absolute amount of nitrogen allocated to carboxylation (gN/m2 leaf) 
      real(r8) :: Nresp1                                  !absolute amount of nitrogen allocated to respiration due to increase in light capture nitrogen(gN/m2 leaf)  
      real(r8) :: Nlc1                                    !absolute amount of nitrogen allocated to light capture due to increase in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: Net1                                    !absolute amount of nitrogen allocated to electron transport due to increase in light capture nitrogen(gN/m2 leaf)
      real(r8) :: Ncb1                                    !absolute amount of nitrogen allocated to carboyxlation due to increase in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: Nresp2                                  !absolute amount of nitrogen allocated to respiration due to decrease in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: Nlc2                                    !absolute amount of nitrogen allocated to light capture due to decrease in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: Net2                                    !absolute amount of nitrogen allocated to electron transport due to decrease in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: Ncb2                                    !absolute amount of nitrogen allocated to carboxylation due to increase in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: PSN                                     !g carbon photosynthesized per day per unit(m2) of leaf
      real(r8) :: RESP                                    !g carbon respired per day per unit(m2) of leaf due to increase in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: PSN1                                    !g carbon photosynthesized per day per unit(m2) of leaf due to increase in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: RESP1                                   !g carbon respired per day per unit(m2) of leaf due to decrease in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: PSN2                                    !g carbon photosynthesized per day per unit(m2) of leaf due to decrease in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: RESP2                                   !g carbon respired per day per unit(m2) of leaf
      real(r8) :: Npsntarget                              !absolute amount of target nitrogen for photosynthesis(gN/m2 leaf) 
      real(r8) :: Npsntarget1                             !absolute amount of target nitrogen for photosynthesis due to increase in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: Npsntarget2                             !absolute amount of target nitrogen for photosynthesis due to decrease in light capture nitrogen(gN/m2 leaf) 
      real(r8) :: NUEj                                    !nitrogen use efficiency for electron transport under current environmental conditions 
      real(r8) :: NUEc                                    !nitrogen use efficiency for carboxylation under current environmental conditions  
      real(r8) :: NUEjref                                 !nitrogen use efficiency for electron transport under reference environmental conditions (25oC and 385ppm Co2) 
      real(r8) :: NUEcref                                 !nitrogen use efficiency for carboxylation under reference environmental conditions (25oC and 385ppm Co2) 
      real(r8) :: NUEr                                    !nitrogen use efficiency for respiration 
      real(r8) :: PARi10c                                 !10-day mean constrained photosynthetic active radiation on in a canopy (umol/m2/s)
      real(r8) :: PARimx10c                               !10-day mean constrained 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
      real(r8) :: Kj2Kcref                                !the ratio of rubisco-limited photosynthetic rate (Wc) to light limited photosynthetic rate (Wj)
      real(r8) :: PNlcoldi                                !old value of the proportion of nitrogen allocated to light capture (unitless) 
      real(r8) :: Kj2Kc                                   !the ratio of Wc to Wj under changed conditions 
      real(r8) :: Kc                                      !conversion factors for Vc,max to Wc 
      real(r8) :: Kj                                      !conversion factor for electron transport rate to Wj 
      real(r8) :: theta                                   !efficiency of light energy conversion (unitless) 
      real(r8) :: chg_per_step                            !the nitrogen change per interation
      real(r8) :: Vcmaxnight                              !Vcmax during night (umol/m2/s)
      real(r8) :: ci                                      !inter-cellular CO2 concentration (Pa)
      real(r8) :: theta_cj                                !interpolation coefficient
      real(r8) :: tleafd10c                               !10-day mean daytime leaf temperature, contrained for physiological range (oC)
      real(r8) :: tleafn10c                               !10-day mean leaf temperature for night, constrained for physiological range (oC)
      real(r8) :: Vcmax                                   !the maximum carboxyaltion rate (umol/m2/s) 
      real(r8) :: rhol                                    ! leaf reflectance: 1=vis, 2=nir  
      real(r8) :: taul                                    ! leaf transmittance: 1=vis, 2=nir 
      integer  :: KcKjFlag                                !flag to indicate whether to update the Kc and Kj using the photosynthesis subroutine; 0--Kc and Kj need to be calculated; 1--Kc and Kj is prescribed.
      integer  :: jj                                      !index record fo the number of iterations
      integer  :: increase_flag                           !whether to increase or decrease
      
      !------------------------------------------------------------------------------ 
      !Constants  
      real(r8), parameter :: Cv = 1.2e-5_r8 * 3600.0           ! conversion factor from umol CO2 to g carbon
      real(r8), parameter :: Kc25 = 40.49_r8                   ! Mechanis constant of CO2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      real(r8), parameter :: Ko25 = 27840_r8                   ! Mechanis constant of O2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      real(r8), parameter :: Cp25 = 4.275_r8                   ! CO2 compensation point at 25C (Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      real(r8), parameter :: Fc25 = 294.2_r8                   ! Fc25 = 6.22*47.3 #see Rogers (2014) Photosynthesis Research 
      real(r8), parameter :: Fj25 = 1257.0_r8                  ! Fj25 = 8.06*156 # #see COSTE 2005 and Xu et al 2012
      real(r8), parameter :: NUEr25 = 33.69_r8                 ! nitrogen use efficiency for respiration, see Xu et al 2012
      real(r8), parameter :: Cb = 1.78_r8                      ! nitrogen use effiency for choloraphyll for light capture, see Evans 1989  
      real(r8), parameter :: O2ref = 209460.0_r8               ! ppm of O2 in the air
      real(r8), parameter :: CO2ref = 380.0_r8                 ! reference CO2 concentration for calculation of reference NUE. 
      real(r8), parameter :: forc_pbot_ref = 101325.0_r8       ! reference air pressure for calculation of reference NUE
      real(r8), parameter :: Jmaxb0 = 0.0311_r8                ! the baseline proportion of nitrogen allocated for electron transport (J)     
      real(r8), parameter :: Jmaxb1 = 0.1745_r8                ! the baseline proportion of nitrogen allocated for electron transport (J)    
      real(r8), parameter :: Wc2Wjb0 = 0.8054_r8               ! the baseline ratio of rubisco limited rate vs light limited photosynthetic rate (Wc:Wj) 
      real(r8), parameter :: relhExp = 6.0999_r8               ! electron transport parameters related to relative humidity
      real(r8), parameter :: NMCp25 = 0.715_r8                 ! estimated by assuming 80% maintenance respiration is used for photosynthesis enzyme maintenance
      real(r8), parameter :: Trange1 = 5.0_r8                  ! lower temperature limit (oC) for nitrogen optimization  
      real(r8), parameter :: Trange2 = 42.0_r8                 ! upper temperature limit (oC) for nitrogen optimization
      real(r8), parameter :: SNC = 0.004_r8                    ! structural nitrogen concentration (g N g-1 dry mass carbon)
      real(r8), parameter :: mp = 9.0_r8                       ! slope of stomatal conductance; this is used to estimate model parameter, but may need to be updated from the physiology file, 
      real(r8), parameter :: PARLowLim = 200.0_r8              ! minimum photosynthetically active radiation for nitrogen optimization
      real(r8), parameter :: minrelh = 0.25_r8                 ! minimum relative humdity for nitrogen optimization
     
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
     
      ! SLA is fixed through the canopy in SDGVM and lnc is pre-determined according to Beer's Law scaling or similar 
      ! - therefore relCLNCa, PARTop are not needed
      ! PAR is already passed to this routine in umol m-2 s-1 (check not day)
      ! - therefore the intermediary steps to convert from wm-2 are not necessary

      !------------------------------------------------------------------
      !SNCa     =  1.0_r8/slatop(ft) * SNC !(Eq A3)
      SNCa = 1.0_r8/sla/0.48_r8 * SNC     !(Eq A3)
      FNCa = lnc(z) - SNCa                !(Eq A4)
      
      !------------------------------------------------------------------
      ! for SDGVM the distinction between sunlit and non-sunlit leaves is not made in the LUNA model
      ! - this is because during the day leaves experience both sunlit and non-sunlit states so  
      ! par240x_z is 10-day mean, maximum PAR for leaves 
      ! ratio of absorbed light to incident light (rhol and taul are reflectance and transmittance)
      rabsorb     = 1.0_r8-rhol-taul
      radmax2mean = par240x_z(z) / par240d_z(z)
      PARi10      = par240d_z(z) / rabsorb
      PARimx10    = PARi10*radmax2mean
      
      !------------------------------------------------------------------
      !nitrogen allocation model-start          
      PNlcold     = PNlc_z(z)
      PNetold     = 0.0_r8
      PNrespold   = 0.0_r8
      PNcbold     = 0.0_r8                                     
      !end brought in from above subroutine 
      
      ! this call to N allocation has been replaced by the N allocation code 
      !call NitrogenAllocation(FNCa,forc_pbot10(p), relh10, CO2a10, O2a10, PARi10, PARimx10, rb10v, hourpd, &
      !     tair10, tleafd10, tleafn10, &
      !     Jmaxb0, Jmaxb1, Wc2Wjb0, relhExp,  PNstoreold, PNlcold, PNetold, PNrespold, &
      !     PNcbold, PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt)
      
      ! set referfence NUE parameters 
      ! - using sinlge common NUE function fed with default values as arguments
      ! call NUEref(NUEjref, NUEcref, Kj2Kcref)
      call NUE(O2ref*0.1013_r8, 0.7d0*CO2ref*0.1013_r8, 25.d0, 25.d0,
     &NUEjref, NUEcref, Kj2Kcref)
      
      ! appears theta_cj is not necessary for this subroutine
      !theta_cj = 0.95_r8
      ! it's not currently clear to me where PNstoreold is passed from in the original code
      Nstore   = PNstoreold * FNCa                         !FNCa * proportion of storage nitrogen in functional nitrogen
      Nlc      = PNlcold    * FNCa                         !FNCa * proportion of light capturing nitrogen in functional nitrogen
      Net      = PNetold    * FNCa                         !FNCa * proportion of light harvesting (electron transport) nitrogen in functional nitrogen
      Nresp    = PNrespold  * FNCa                         !FNCa * proportion of respiration nitrogen in functional nitrogen
      Ncb      = PNcbold    * FNCa                         !FNCa * proportion of carboxylation nitrogen in functional nitrogen
      if (Nlc > FNCa * 0.5_r8) Nlc = 0.5_r8 * FNCa
      PNlc     = PNlcold
      PNlcoldi = PNlcold  - 0.001_r8
      
      ! constrain the physiological range of PAR and t
      PARi10c   = max(PARLowLim, PARi10)
      PARimx10c = max(PARLowLim, PARimx10)
      tleafd10c = min(max(tleafd10, Trange1), Trange2)  
      tleafn10c = min(max(tleafn10, Trange1), Trange2) 
      
      ! initialse solver parameters 
      ci            = 0.7_r8 * CO2a10 
      JmaxCoef      = Jmaxb1 * ((hourpd / 12.0_r8)**2.0_r8) * 
     &(1.0_r8 - exp(-relhExp * max(relh10 - minrelh, 0.0_r8) / 
     &(1.0_r8 - minrelh)))
      chg_per_step  = 0.02* FNCa
      increase_flag = 0
      
      !------------------------------------------------------------------
      jj = 1
      do while ( (PNlcoldi .NE. PNlc) .and. (jj < 100) )      
     
       ! Fc is the scaling factor to go from leaf N invested in RuBisCO to Vcmax  
       ! Fj is the scaling factor to go from leaf N invested in elec trans to Jmax  
       Fc   = VcmxTKattge(tair10, tleafd10c) * Fc25
       Fj   = JmxTKattge(tair10, tleafd10c)  * Fj25
       NUEr = Cv * NUEr25 * (RespTBernacchi(tleafd10c) * hourpd +
     &RespTBernacchi(tleafn10c) * (24.0_r8 - hourpd)) !nitrogen use efficiency for respiration (g biomass/m2/day/g N)
       
       call NUE(O2a10, ci, tair10, tleafd10c, NUEj, NUEc, Kj2Kc)
       
       KcKjFlag = 0
       call LUNA_Ninvestments (KcKjFlag,FNCa, Nlc, forc_pbot10,relh10,
     &CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, tair10, 
     &tleafd10c,tleafn10c, 
     &Kj2Kc, Wc2Wjb0, JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref,NUEjref,NUEr,
     &Kc, Kj, ci, Vcmax, Jmax,JmeanL,JmaxL, Net, Ncb, Nresp, PSN, RESP)
       
       !target nitrogen allocated to photosynthesis, which may be lower or higher than Npsn_avail
       Npsntarget = Nlc + Ncb + Net       
       PNlcoldi   = Nlc / FNCa
       Nstore     = FNCa - Npsntarget - Nresp
       
       !test the increase of light capture nitrogen
       if ( ((Nstore > 0.0_r8) .and.
     &(increase_flag .eq. 1)) .or. (jj .eq. 1) ) then
        Nlc2 = Nlc + chg_per_step
        if (Nlc2 / FNCa > 0.95_r8) Nlc2 = 0.95_r8 * FNCa
        
        KcKjFlag = 1
        call LUNA_Ninvestments (KcKjFlag,FNCa, Nlc2, forc_pbot10, 
     &relh10, CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, 
     &tair10, tleafd10c,tleafn10c, 
     &Kj2Kc, Wc2Wjb0, JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref,NUEjref,NUEr,
     &Kc,Kj,ci,Vcmax,Jmax,JmeanL,JmaxL, Net2, Ncb2, Nresp2, PSN2, RESP2)
        
        Npsntarget2 = Nlc2 + Ncb2 + Net2
        
        !update the nitrogen change
        Carboncost2 = (Npsntarget2 - Npsntarget) * NMCp25 * Cv * 
     &(RespTBernacchi(tleafd10c) * hourpd + RespTBernacchi(tleafn10c) *
     &(24.0_r8  - hourpd))
        Carbongain2 =  PSN2 - PSN
        
        if( (Carbongain2 > Carboncost2) .and. 
     &(Npsntarget2 + Nresp2 < 0.95_r8 * FNCa) ) then
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
       
       if (Nstore < 0.0_r8) then
        Nlc1 = Nlc * 0.8_r8 !bigger step of decrease if it is negative            
       else
        Nlc1 = Nlc - chg_per_step
       end if
       
       if (Nlc1 < 0.05_r8) Nlc1 = 0.05_r8
        KcKjFlag = 1
        call LUNA_Ninvestments(KcKjFlag,FNCa, Nlc1,forc_pbot10,relh10,
     &CO2a10,O2a10, PARi10c, PARimx10c,rb10, hourpd, 
     &tair10, tleafd10c,tleafn10c, 
     &Kj2Kc, Wc2Wjb0,JmaxCoef, Fc,Fj, NUEc, NUEj, NUEcref,NUEjref,NUEr,
     &Kc,Kj,ci,Vcmax,Jmax,JmeanL,JmaxL,Net1, Ncb1, Nresp1, PSN1, RESP1)
        
        Npsntarget1 = Nlc1 + Ncb1 + Net1
        Carboncost1 = (Npsntarget - Npsntarget1) * NMCp25 * Cv *
     &( RespTBernacchi(tleafd10c) * hourpd + RespTBernacchi(tleafn10c) *
     &(24.0_r8  - hourpd) )
        Carbongain1 =  PSN - PSN1
        
        if( ((Carbongain1 < Carboncost1) .and. (Nlc1 > 0.05_r8)) .or.
     &(Npsntarget + Nresp > 0.95_r8 * FNCa) ) then
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
      
      !brought in from above subroutine

      ! determine change in vcmax and jmax 
      vcmx25_opt  = PNcbopt * FNCa * Fc25
      jmx25_opt   = PNetopt * FNCa * Fj25
      
      chg         = vcmx25_opt-vcmx25_z(p, z)
      chg_constrn = min(abs(chg),vcmx25_z(p, z)*max_daily_pchg)
      vcmx25_z(z) = vcmx25_z(p, z)+sign(1.0_r8,chg)*chg_constrn
       
      chg         = jmx25_opt-jmx25_z(p, z)
      chg_constrn = min(abs(chg),jmx25_z(p, z)*max_daily_pchg)
      jmx25_z(z)  = jmx25_z(p, z)+sign(1.0_r8,chg)*chg_constrn 
      
      PNlc_z(z)   = PNlcopt
      
      if(enzs_z(z)<1.0) enzs_z(z) = enzs_z(z)* (1.0_r8 + max_daily_pchg)

      !nitrogen allocation model-end  
      !------------------------------------------------------------------

      if(isnan(vcmx25_z(z)).or. vcmx25_z(z)>1000._r8 .or. 
     &vcmx25_z(z)<0._r8) then
       write(iulog, *) 'Error: Vc,mx25 become unrealistic (NaN,>1000,
     & or negative) for z=', z
       write(iulog, *) 'LUNA env:',FNCa,forc_pbot10, relh10, CO2a10, 
     &O2a10, PARi10, PARimx10, rb10v, 
     &hourpd, tair10, tleafd10, tleafn10
       call endrun(msg=errmsg(__FILE__, __LINE__))
      endif
      if(isnan(jmx25_z(u)).or.jmx25_z(z)>1000._r8 .or. 
     &jmx25_z(z)<0._r8)then
       write(iulog, *) 'Error: Jmx25 become unrealistic (NaN,>1000, 
     &or negative)for z=', z
       write(iulog, *) 'LUNA env:', FNCa,forc_pbot10, relh10, CO2a10,
     &O2a10, PARi10, PARimx10, rb10vi,hourpd, tair10, tleafd10, tleafn10
       call endrun(msg=errmsg(__FILE__, __LINE__))
      endif
      !end brought in from above subroutine 
      
      end subroutine LUNA
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      subroutine LUNA_Ninvestments (KcKjFlag, FNCa, Nlc, forc_pbot10, 
     &relh10,CO2a10,O2a10,PARi10,PARimx10,rb10,hourpd,tair10,tleafd10,
     &tleafn10,Kj2Kc, Wc2Wjb0, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, 
     &NUEjref,NUEr,Kc,Kj,ci,Vcmax,Jmax,JmeanL,JmaxL,Net,Ncb,Nresp,PSN,
     &RESP)
      
      !calculate the nitrogen investment for electron transport, carb10oxylation, respiration given a specified value 
      !of nitrogen allocation in light capture [Nlc]. This equation are based on Ali et al 2015b.
      
      implicit none
      integer,  intent (in) :: KcKjFlag                   !flag to indicate whether to update the Kc and Kj using the photosynthesis subroutine; 0--Kc and Kj need to be calculated; 1--Kc and Kj is prescribed.
      real(r8), intent (in) :: FNCa                       !Area based functional nitrogen content (g N/m2 leaf)
      real(r8), intent (in) :: Nlc                        !nitrogen content for light capture(g N/m2 leaf)
      real(r8), intent (in) :: forc_pbot10                !10-day mean air pressure (Pa)
      real(r8), intent (in) :: relh10                     !10-day mean relative humidity (unitless)
      real(r8), intent (in) :: CO2a10                     !10-day mean CO2 concentration in the air (Pa)
      real(r8), intent (in) :: O2a10                      !10-day mean O2 concentration in the air (Pa)
      real(r8), intent (in) :: PARi10                     !10-day mean photosynthetic active radiation on in a canopy (umol/m2/s)
      real(r8), intent (in) :: PARimx10                   !10-day mean 24hr maximum photosynthetic active radiation on in a canopy (umol/m2/s)
      real(r8), intent (in) :: rb10                       !10-day mean boundary layer resistance (s/m)
      real(r8), intent (in) :: hourpd                     !hours of light in a the day (hrs)
      real(r8), intent (in) :: tair10                     !10-day running mean of the 2m temperature (oC)
      real(r8), intent (in) :: tleafd10                   !10-day mean daytime leaf temperature (oC) 
      real(r8), intent (in) :: tleafn10                   !10-day mean nighttime leaf temperature (oC) 
      real(r8), intent (in) :: Kj2Kc                      !ratio:  Kj / Kc
      real(r8), intent (in) :: Wc2Wjb0                    !the baseline ratio of rubisco-limited rate vs light-limited photosynthetic rate (Wc:Wj)
      real(r8), intent (in) :: JmaxCoef                   !coefficient determining the response of electron transport rate to light availability (unitless) and humidity
      real(r8), intent (in) :: Fc                         !the temperature adjustment factor for Vcmax 
      real(r8), intent (in) :: Fj                         !the temperature adjustment factor for Jmax 
      real(r8), intent (in) :: NUEc                       !nitrogen use efficiency for carboxylation 
      real(r8), intent (in) :: NUEj                       !nitrogen use efficiency for electron transport
      real(r8), intent (in) :: NUEcref                    !nitrogen use efficiency for carboxylation under reference climates
      real(r8), intent (in) :: NUEjref                    !nitrogen use efficiency for electron transport under reference climates
      real(r8), intent (in) :: NUEr                       !nitrogen use efficiency for respiration
      real(r8), intent (inout) :: Kc                      !conversion factors from Vc,max to Wc 
      real(r8), intent (inout) :: Kj                      !conversion factor from electron transport rate to Wj 
      real(r8), intent (inout) :: ci                      !inter-cellular CO2 concentration (Pa) 
      real(r8), intent (out) :: Vcmax                     !the maximum carboxyaltion rate (umol/m2/s) 
      real(r8), intent (out) :: Jmax                      !the maximum electron transport rate (umol/m2/s) 
      real(r8), intent (out) :: JmaxL                     !the electron transport rate with maximum daily radiation (umol/m2/s)  
      real(r8), intent (out) :: JmeanL                    !the electron transport rate with mean radiation (umol/m2/s) 
      real(r8), intent (out)  :: Net                      !nitrogen content for electron transport(g N/m2 leaf)
      real(r8), intent (out)  :: Ncb                      !nitrogen content for carboxylation(g N/m2 leaf)
      real(r8), intent (out)  :: Nresp                    !nitrogen content for respiration(g N/m2 leaf)
      real(r8), intent (out)  :: PSN                      !daily photosynthetic rate(g C/day/m2 leaf)
      real(r8), intent (out)  :: RESP                     !daily respiration rate(g C/day/m2 leaf)
      !-------------------------------------------------------------------------------------------------------------------------------
      !intermediate variables
      real(r8) :: A                                       !Gross photosynthetic rate (umol CO2/m2/s)
      real(r8) :: Wc2Wj                                   !ratio: Wc/Wj  
      real(r8) :: ELTRNabsorb                             !absorbed electron rate, umol electron/m2 leaf /s
      real(r8) :: Jmaxb0act                               !base value of Jmax (umol/m2/s) 
      real(r8) :: theta_cj                                !interpolation coefficient
      real(r8) :: theta                                   !light absorption rate (0-1)
      real(r8) :: Vcmaxnight                              !Vcmax during night (umol/m2/s)
      real(r8) :: Wc                                      !rubisco-limited photosynthetic rate (umol/m2/s)
      real(r8) :: Wj                                      !light limited photosynthetic rate (umol/m2/s)
      real(r8) :: NUECHG                                  !the nitrogen use efficiency change under current conidtions compared to reference climate conditions (25oC and 385 ppm )
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      theta_cj    = 0.95_r8
      theta       = 0.292_r8 / (1.0_r8 + 0.076_r8 / (Nlc * Cb))
      ELTRNabsorb = theta  * PARi10
      Jmaxb0act   = Jmaxb0 * FNCa * Fj
      Jmax        = Jmaxb0act + JmaxCoef * ELTRNabsorb
      JmaxL       = theta  * PARimx10 / (sqrt(1.0_r8 + 
     &(theta * PARimx10 / Jmax)**2.0_r8))        
      NUEchg      = (NUEc / NUEcref) * (NUEjref / NUEj)
      Wc2Wj       = Wc2Wjb0 * (NUEchg**0.5_r8)
      Wc2Wj       = min(1.0_r8, Wc2Wj)
      Vcmax       = Wc2Wj * JmaxL * Kj2Kc
      JmeanL      = theta * PARi10 / (sqrt(1.0_r8 +
     &(ELTRNabsorb / Jmax)**2.0_r8))
      
      if(KcKjFlag.eq.0) then      !update the Kc,Kj, anc ci information
       
       ! From SDGVM - for consistency
       ! in LUNA ko, kc, and tau are at the 10-day mean leaf temp
       ! - need to check that units are consistent, in SDGVM these are all in Pa
       kc  = exp(35.8d0   - 80.5d0 /(0.00831d0*(tleafd10+273.15)))
       ko  = exp(9.6d0    - 14.51d0/(0.00831d0*(tleafd10+273.15)))*1000.0d0
       tau = exp(-3.949d0 + 28.99d0/(0.00831d0*(tleafd10+273.15)))
       
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
       CALL ASSIMILATION_CALC(0.0d0,1.0d0,
     &Vcmax,0.0d0,rd,0.0d0,JmeanL,0.0d0,0.0d0,
     &0.0d0,tleafd10,1.0d0,1.0d0,0.0d0,ga,relh10,1,1.0d0,ko,kc,tau,
     &forc_pbot10,O2a10,CO2a10,A,gs,ci,1,2,1,1,.FALSE.,gs_func,ftg0,
     &ftg1,z)
       
      else
       Wc = Kc * Vcmax
       Wj = Kj * JmeanL
       A  = (1.0_r8 - theta_cj) * max(Wc, Wj) + theta_cj * min(Wc, Wj)
      endif
      
      ! Cv converts from umolm-2s-1 to gC m-2 hour-1
      PSN        = Cv * A * hourpd
      Vcmaxnight = VcmxTKattge(tair10, tleafn10) / 
     &VcmxTKattge(tair10, tleafd10) * Vcmax
      RESP       = Cv * 0.015_r8 * (Vcmax * hourpd + Vcmaxnight * (24.0_r8 - hourpd))                 
      Net        = Jmax  / Fj
      Ncb        = Vcmax / Fc
      Nresp      = RESP  / NUEr
      
      end subroutine LUNA_Ninvestments
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      !Calculate the Nitrogen use effieciency dependence on CO2 and leaf temperature
      subroutine NUE(o2a, ci, tgrow, tleaf, NUEj,NUEc,Kj2Kc)
      
      ! this function can be used to calculate NUEref aswell, it just needs to be called with both temps at 25oC and a ci that maintains a ci:ca of 0.7
      
      ! uses module constants of Fc25, Fj25, tfrz, rgas
      
      implicit none
      real(r8), intent (in) :: o2a                        !air O2 partial presuure (Pa)
      real(r8), intent (in) :: ci                         !leaf inter-cellular [CO2] (Pa) - originally incorrectly labelled as PPM
      real(r8), intent (in) :: tgrow                      !10 day growth temperature (oC), 24 hour mean temperature
      real(r8), intent (in) :: tleaf                      !leaf temperature (oC)
      real(r8), intent (out):: NUEj                       !nitrogen use efficiency for electron transport under refernce environmental conditions (25oC and 385 ppm co2)
      real(r8), intent (out):: NUEc                       !nitrogen use efficiency for carboxylation under reference environmental conditions  (25oC and 385 ppm co2)
      real(r8), intent (out):: Kj2Kc                      !the ratio of Kj to Kc 
      !------------------------------------------------
      !intermediate variables
      real(r8) :: Fj                                      !the temperatuer adjust factor for Jmax 
      real(r8) :: Fc                                      !the temperatuer adjust factor for Vcmax 
      real(r8) :: Kc                                      !conversion factor from Vcmax to Wc 
      real(r8) :: Kj                                      !conversion factor from J to W 
      real(r8) :: k_o                                     !Rubsico O2 specifity
      real(r8) :: k_c                                     !Rubsico CO2 specifity
      real(r8) :: awc                                     !second deminator term for rubsico limited carboxylation rate based on Farquhar model
      real(r8) :: c_p                                     !CO2 compenstation point (Pa)
      
      real(r8), parameter :: Fc25 = 294.2_r8              ! Fc25 = 6.22*47.3 #see Rogers (2014) Photosynthesis Research 
      real(r8), parameter :: Fj25 = 1257.0_r8             ! Fj25 = 8.06*156  #see COSTE 2005 and Xu et al 2012
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      Fc  = VcmxTKattge(tgrow, tleaf) * Fc25
      Fj  = JmxTKattge(tgrow, tleaf)  * Fj25
c      k_c = Kc25 * exp((79430.0_r8 / (rgas*1.e-3_r8 * 
c     &(25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
c      k_o = Ko25 * exp((36380.0_r8 / (rgas*1.e-3_r8 *
c     &(25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
c      c_p = Cp25 * exp((37830.0_r8 / (rgas*1.e-3_r8 * 
c     &(25.0_r8 + tfrz))) * (1.0_r8 - (tfrz + 25.0_r8) / (tfrz + tleaf)))
      
      ! From SDGVM - for consistency
      ! in LUNA ko, kc, and tau are at the 10-day mean leaf temp
      ! - need to check that units are consistent, in SDGVM these are all in Pa
      k_c = exp(35.8d0   - 80.5d0 /(0.00831d0*(tleafd10+273.15)))
      k_o = exp(9.6d0    - 14.51d0/(0.00831d0*(tleafd10+273.15)))*1000.0d0
      tau = exp(-3.949d0 + 28.99d0/(0.00831d0*(tleafd10+273.15)))
      c_p =  0.5d0*o2a/tau
     
      awc = k_c * ( 1.0_r8 + o2a/k_o )
      Kj  = max( ci-c_p,0.0_r8 ) / ( 4.0_r8*ci + 8.0_r8*c_p )
      Kc  = max( ci-c_p,0.0_r8 ) / ( ci+awc )
     
      NUEj  = Kj * Fj
      NUEc  = Kc * Fc  
      Kj2Kc = Kj / Kc
     
      end subroutine NUE
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      
      
      
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
      !Calculate the temperature response for respiration, following Bernacchi PCE 2001
      real(r8) function  RespTBernacchi(tleaf)
      implicit none
      real(r8), intent(in):: tleaf  !leaf temperature (oC)
      
      RespTBernacchi= exp(18.72_r8-46.39_r8/(rgas*1.e-6_r8 *
     &(tleaf+tfrz)))
      
      end function RespTBernacchi
      !-------------------------------------------------------------------------------------------------------------------------------------------------       
