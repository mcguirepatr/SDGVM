*----------------------------------------------------------------------*
*                          SUBROUTINE NPPCALC                          *
*                          ******************                          *
*----------------------------------------------------------------------*
      SUBROUTINE NPPCALC(npp_eff,c3,soilc,soiln,minn,soil2g,wtwp,
     &wtfc,rlai,
     &t,rh,ca,oi,rn,qdirect,qdiff,can2a,can2g,canrd,canres,suma,amx,
     &amax,hrs,canga,p,mnth,day,nleaf_sum,fpr,gsm,tleaf_n,tleaf_p,
     &ncalc_type,vcmax_type,
     &leaf_nit,vcmax,jmax,pnlc,enzs,ft,kg,sla,can_clump,tassim,tgs,
     &tci,hw_j,cstype,subd_par,thty_dys,year,lat,swr,cld,
     &read_par,env_vcmax,env_jmax,soilp_map,ga,ftvna,ftvnb,ftjva,ftjvb,
     &ftg0,ftg1,par_loops,s070607,gs_func,ce_light,ce_ci,ce_t,
     &ce_maxlight,ce_ga,ce_rh,ttype,
     &calc_zen,cos_zen,ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ)
*----------------------------------------------------------------------*
      
      IMPLICIT NONE
      
      INCLUDE 'array_dims.inc'
      INCLUDE 'param.inc'

      REAL*8 suma,sumd,rlai,soilc,soil2g,wtfc,nmult,npp_eff
      REAL*8 soiln,y1,y0,x1,x0,mmult,minn(3),nup,t
      REAL*8 tk,up,kc,ko,tau,can_sum,nleaf_sum,cld,gb
      REAL*8 rem,can(12),vm(12),jm(12),jmx(12),oi
      REAL*8 canres,dresp(12),drespt(12),upt(12),nupw
      REAL*8 vmx(12),rh,c1,c2,canga,p,ga,nleaf(12)
      REAL*8 q,qdiff,qdirect,cs,sla,can_clump
      REAL*8 q_sd,qdiff_sd,brent_solver,t_scalar
      REAL*8 qdirect_sd,canrd,cos_zen,k
      REAL*8 can2a,can2g,rn,wtwp,kg,rht,tpav
      REAL*8 tppcv,tpgsv,ca,rd(12),tpaj,tppcj,tpgsj,hrs
      REAL*8 a(12),gs(12),ci(12),tpac4,tppcc4,tpgsc4,xvmax
      REAL*8 a_sd(12),gs_sd(12),ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ
      REAL*8 ci_sd(12),ftvna,ftvnb,ftjva,ftjvb,ftg0,ftg1
      REAl*8 apar,fpr,gsm,ncan,pleaf(12),slaleaf,tleaf,pcan,slaton
      REAL*8 ce_light(30,12),ce_ci(30,12),ce_t(30)
      REAL*8 ce_maxlight(30,12),ce_ga(30,12),ce_rh(30),maxlight(12)
      REAL*8 asunlit,ashade,pcsunlit,pcshade,gssunlit,gsshade
* light-limited assimilation rate (j) and irradiance (q) for sunlit
* and shade
      REAL*8 jsunlit(12),jshade(12),qsunlit(12),qshade(12)
      REAL*8 jsunlit_sd(12),jshade_sd(12)
      REAL*8 qsunlit_sd(12),qshade_sd(12)
* fraction of sunlit and shade
      REAL*8 fsunlit(12),fshade(12),tleaf_n,tleaf_p
      REAL*8 fsunlit_sd(12),fshade_sd(12)
      REAL*8 leaf_nit,vcmax(12),jmax(12),pnlc(12),enzs(12)
      REAL*8 tassim,tgs,tci
      REAL*8 max_daily_pchg, max_dpchg

      REAL*8 ax,amax,amx,lyr,qt,lat,swr,env_vcmax,env_jmax
      INTEGER i,lai,c3,mnth,day,ncalc_type,vcmax_type,ft,oday
      INTEGER hw_j,cstype,subd_par,ii,omnth,par_loops,gs_func
      INTEGER thty_dys,year,no_day,read_par,soilp_map,s070607
      INTEGER ttype,calc_zen,luna_calc_days
      REAL*8  sd_scale(par_loops+1),sd_scale2,nup_rate
      LOGICAL output,gold

      !print*, 'NPPcalc', day

      omnth  = 7
      oday   = 10
      output = .FALSE.
   
      if(hw_j.eq.3) then
        hw_j = 0 
        gold = .TRUE.
      endif 
 
      ! parameter - interval between luna calculations of Vcmax (in days)
      luna_calc_days = 10

      !if(output.and.(mnth.eq.omnth).and.(day.eq.oday)) then
       ! print*, 'NPPCALC'
       ! print*, t,rn,soil2g,wtwp
       ! print*, ''
      !endif

      !print*, vcmax(1)
      !print*, ce_light(30,1),ce_maxlight(30,1)  
      !print*, sum(ce_light(:,1))/30.d0,sum(ce_maxlight(:,1))/30.d0,
      !&vcmax(1),jmax(1)  
      !print*, qdirect, qdirect, qdirect, qdiff
      apar=0.0d0

      rem = rlai - int(rlai)
      lai = int(rlai) + 1
      tk = 273.0d0 + t

      suma = 0.0d0
      sumd = 0.0d0
      if(vcmax_type.ne.9) then
        vcmax(:) = 0.0d0
        jmax(:)  = 0.0d0
      endif
      leaf_nit  = 0.0d0
      nleaf_sum = 0.0d0
      canres    = 0.0d0
      nleaf(:)  = 0.0d0

      kc  = exp(35.8d0 - 80.5d0/(0.00831d0*tk))
      ko  = exp(9.6d0 - 14.51d0/(0.00831d0*tk))*1000.0d0
      tau = exp(-3.949d0 + 28.99d0/(0.00831d0*tk))

      ! changed by Ghislain 08/10/03      IF (rlai.GT.0.1d0) THEN
      q = qdiff + qdirect
      
      !if there are light and leaves - calculate canopy properties 
      !print*, '' 
      IF ((rlai.GT.0.1d0).AND.(q.GT.0.0d0)) THEN
      !print*, 'calculate N and Vcmax'

!      IF (soil2g.GT.wtwp) THEN
!        kg = maxc*((soil2g - wtwp)/(wtfc - wtwp))**p_kgw
!        IF (kg.GT.maxc)  kg = maxc
!      ELSE
!        kg = 0.0d0
!      ENDIF

      ! Nitrogen uptake - occurs on a daily basis and is allocated on a daily basis
      !                 - i.e. no N storage from day to day, this is a bit weird
!      nupw  = p_nu1*(exp(p_nu2*soilc))*kg**p_nu3
!      IF (nupw.LT.0.0d0)   nupw = 0.0d0

!      nmult = soiln*p_nu4
!      IF (nmult.GE.1.0d0)  nmult = 1.0d0

!      nup   = nupw*nmult
!      up    = nup
!      up    = nupw*nmult
!      if (up.LT.0.0d0)     up = 0.0d0
      
      if(ncalc_type.eq.0) then
        up = NUP_RATE(soilc,soiln,p_nu1,p_nu2,p_nu4)*kg**p_nu3
      elseif(ncalc_type.eq.1) then
        ! as above but without soil water limitation effect
        up = NUP_RATE(soilc,soiln,p_nu1,p_nu2,p_nu4)
      endif

      ! canopy N scaling
      k = 0.5d0
      if(cstype.eq.1) k = 0.5d0 * can_clump
      
      ! Total the Beer's Law values for canopy
      can_sum = 0.0d0
      DO i=1,lai-1
        if(cstype.lt.2) then 
          can_sum = can_sum + exp(-k*real(i))
        elseif(cstype.eq.2) then
          ! use light proportion to scale
          can_sum = can_sum + sum(ce_light(:,i))
        else
          PRINT*, 'cstype ',cstype,' undefined. set to a value
     &<3 in <input.dat> or define additional canopy scaling method'
          STOP
        endif
     
      ENDDO
      if(cstype.lt.2) can_sum = can_sum + exp(-k*real(lai))*rem
      if(cstype.eq.2) can_sum = can_sum + sum(ce_light(:,lai))*rem

      !calculate maximum daily change in vcmax for LUNA model after Ali, Xu, et al 2015
      IF(vcmax_type.eq.9) THEN
         max_dpchg = max_daily_pchg(sum(ce_t(:))/30.d0)  
      ENDIF
      !print*, max_dpchg
 
      !start first LAI loop 
      ! - canopy scaling of variables that do not vary with light
      DO i=1,lai
        lyr = real(i)

        !proportion of canopy N in the LAI layer 'lyr'
        ! - this is supposed to be proportional to the light in layer 'lyr' but this is inconsistent with the rad scheme used by SDGVM
        if(cstype.lt.2) then 
          can(i) = exp(-k*(lyr)) / can_sum
        elseif(cstype.eq.2) then
          ! use light proportion to scale
          can(i) = sum(ce_light(:,i)) / can_sum
        endif
        if((i.eq.lai).and.(s070607.eq.0)) can(i) = can(i) * rem

        ! calculate leaf N in canopy layer i
        IF (ncalc_type.le.1) THEN 
           !default topleaf N 
           nleaf(i) = up*can(i)*p_nleaf
        ELSEIF (ncalc_type.eq.2) THEN
           !based on trait regressions/specified top leaf N
           !calculate total canopy n
           ncan     = tleaf_n*can_sum
           nleaf(i) = can(i)*ncan
        ELSE
           PRINT*, 'ncalc_type ',ncalc_type,' undefined. set to a value
     &<3 in <input.dat> or define your own N calculation method'
           STOP
        ENDIF
        ! this scales the bottom layer from a fractional layer to a full layer
        ! - unless total LAI is <1, in which case leaf N is assumed for a full leaf layer 
        ! calculations below this assume a full leaf layer and then scale all fluxes in lowest lai layer by rem 
        if((i.eq.lai).and.(lai.ne.1)) nleaf(i) = nleaf(i) / rem
        !if(i.eq.lai) print*, rlai, can(i), up*p_nleaf, nleaf(i)

        !if soil water not completely limiting calculate photosynthesis
        IF (kg.GT.1.0e-10) THEN

          !calculate Vcmax@25oC (vm) 
          IF(vcmax_type.eq.0) THEN
            !default 070607 SDGVM
            vm(i) = nleaf(i)*p_vm         

          ELSEIF(vcmax_type.eq.1) THEN
            !from Walker et al - N only
            vm(i) = exp(3.712+0.65*log(nleaf(i)))
            !from Walker et al - N & P
            IF(soilp_map.eq.1) THEN
              !assume same change in leaf p through the canopy as leaf n
              pcan = tleaf_p*can_sum
              pleaf(i) = can(i)*pcan
              vm(i) = exp(3.946+0.921*log(nleaf(i))+
     &0.121*log(pleaf(i))+0.282*log(nleaf(i))*log(pleaf(i)))
            ENDIF

          ELSEIF((vcmax_type.eq.2).OR.(vcmax_type.eq.3).or.
     &(vcmax_type.eq.5).or.(vcmax_type.eq.6)) THEN
            !vcmax based on environment 
            vm(i) = env_vcmax * can(i)/can(1)

          ELSEIF(vcmax_type.eq.4) THEN
            !specified as a PFT parameter
            vm(i) = ftvna + ftvnb*nleaf(i)

          ELSEIF(vcmax_type.eq.7) THEN
            ! after Maire et al 2012
            
            !if(i.eq.1) print'(30(f6.1,1x))', ce_light(:,1) *1d6
            !if(i.eq.1) print'(f6.1)',   sum(ce_light(:,1))/30.d0 *1d6
            
            ! calculate Vcmax at the mean temperature of the last month in mol m-2s-1
            vm(i) = brent_solver(0,1d-6,5d-4,
     &oi,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &sum(ce_t(:))/30.d0,0.d0,1,0.d0,0.d0,0.d0,i,
     &sum(ce_ci(:,i))/30.d0,sum(ce_light(:,i))/30.d0,ttype,
     &ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ,sum(ce_t(:))/30.d0)
            
            ! convert back to value at 25oC 
            ! invert temp correction scalar to get values at 25oC
            vm(i) = vm(i) / T_SCALAR(sum(ce_t(:))/30.d0,'v',ttype,
     &ftToptV,ftHaV,ftHdV,sum(ce_t(:))/30.d0) 
            jm(i) = jm(i) / T_SCALAR(sum(ce_t(:))/30.d0,'j',ttype,
     &ftToptJ,ftHaJ,ftHdJ,sum(ce_t(:))/30.d0) 

            ! convert to umol m-2s-1
            vm(i) = vm(i) * 1d6

          ELSEIF(vcmax_type.eq.8) THEN
            !specified as a PFT parameter but jmax specified as the Walker function of Vcmax
            vm(i) = ftvna + ftvnb*nleaf(i)

          ELSEIF(vcmax_type.eq.9) THEN
            ! use LUNA model after Ali, Xu, et al 2015

            vm(i) = vcmax(i)
            jm(i) = jmax(i)

            IF(mod(day,luna_calc_days).eq.0) THEN
            !print*, 'calc LUNA' 
            !print*, luna_calc_days*max_dpchg,vm(i),jm(i)
              call LUNA(i,vm(i),jm(i),PNlc(i),enzs(i),
     &sum(ce_ga(:,i))/30.d0,
     &sla,nleaf(i),sum(ce_light(:,i))/30.d0,sum(ce_maxlight(:,i))/30.d0,
     &ca,oi,hrs,sum(ce_rh(:))/30.d0,sum(ce_t(:))/30.d0,  
     &gs_func,ftg0,ftg1,max_dpchg*luna_calc_days,ttype,ftToptV,ftHaV,
     &ftHdV,ftToptJ,ftHaJ,ftHdJ)  
            ENDIF

            !print*, 'LUNA vcm:', vm(i)
            ! convert to umol m-2s-1
            !vm(i) = vm(i) * 1d6

          ELSE
            PRINT*, 'vcmax_type:',vcmax_type,'undefined. set to a value
     &1-9 in <input.dat> or define your own vcmax calculation method'
            STOP
          ENDIF
          
          !calculate Jmax@25oC (jm) 
          IF((vcmax_type.eq.0).or.(vcmax_type.eq.5).or.
     &(vcmax_type.eq.6)) THEN
            !070607 default - Wullschleger 1993 
            jm(i)  = (29.1d0 + 1.64d0*vm(i))
          ELSEIF((vcmax_type.eq.1).or.(vcmax_type.eq.7).or.
     &(vcmax_type.eq.8)) THEN
            !Walker -  only vcmax
            jm(i)  = exp(1.d0+0.89d0*log(vm(i)))
          ELSEIF((vcmax_type.eq.2).OR.(vcmax_type.eq.3)) THEN
            !from van Bodegom for TERRABITES project
            ! - not a function of vcmax 
            jm(i) = env_jmax * can(i)/can(1)
          ELSEIF(vcmax_type.eq.4) THEN
            !specified as a PFT parameter
            jm(i)  = ftjva + ftjvb*vm(i)
          ELSEIF(vcmax_type.eq.9) THEN
            ! do nothing - jmax already defined above            
          ELSE
            PRINT*, 'vcmax_type ',vcmax_type,' undefined. set to a value
     &<1-9 in <input.dat>'
            STOP
          ENDIF

          !temperature scale Vcmax & Jmax
          !if(vcmax_type.le.6) then
          vmx(i) = vm(i) * T_SCALAR(t,'v',ttype,ftToptV,ftHaV,ftHdV,
     &sum(ce_t(:))/30.d0)
          jmx(i) = jm(i) * T_SCALAR(t,'j',ttype,ftToptJ,ftHaJ,ftHdJ,
     &sum(ce_t(:))/30.d0)

          !leaf age scale Vcmax & Jmax
          vmx(i) = vmx(i)*npp_eff
          jmx(i) = jmx(i)*npp_eff
          
          !water limitation scale Vcmax & Jmax
          vmx(i) = vmx(i)*kg**p_nu3
          jmx(i) = jmx(i)*kg**p_nu3

          !convert from umol to mol
          vmx(i) = vmx(i) * 1d-6
          jmx(i) = jmx(i) * 1d-6

          !dark respiration (daytime) mol m-2 s-1
          ! assumes rd is proportional to soil water lim
          if(s070607.eq.1) then
            rd(:) = vmx(i)*0.02d0
          else
            !SDGVM previously calculated rd but wasn't layer specific
            rd(i) = vmx(i)*p_rd
          endif

          !if subd_par = 0 then calculate light here  
          IF(subd_par.eq.0) then
          !calculate incident light in canopy layer 
            CALL GOUDRIAANSLAW(lyr-0.5d0,rlai,qdirect,qdiff,
     &fsunlit(i),qsunlit(i), fshade(i),qshade(i),can_clump,cos_zen,
     &s070607,gold)

            !if(output.and.(mnth.eq.omnth).and.(day.eq.oday)) then
              !print*, qdirect+qdiff,qsunlit(i)+qshade(i),i
            !endif

            !calculate the sunlit and shaded light limited 
            !electron transport rate
            CALL JCALC(qsunlit(i),jmx(i),hw_j,jsunlit(i))
            CALL JCALC(qshade(i),jmx(i),hw_j,jshade(i))
          !switch 1 endif
          ENDIF

        !water limitation IF
        ELSE

          IF(subd_par.eq.0) then
            CALL GOUDRIAANSLAW(lyr-0.5d0,rlai,qdirect,qdiff,
     &fsunlit(i),qsunlit(i), fshade(i),qshade(i),can_clump,cos_zen,
     &s070607,gold)

            vmx(i)     = 0.0d0
            jsunlit(i) = 0.0d0
            jshade(i)  = 0.0d0
          ENDIF

          !LUNA model has a term that reduces Vcmax during drought and periods of stress 
          IF(vcmax_type.eq.9) THEN
            vm(i) = vcmax(i)
            jm(i) = jmax(i)

            call LUNA_nogrowth(max_dpchg,vm(i),jm(i),enzs(i),lai)  
            !temperature scale Vcmax & Jmax
            !if(vcmax_type.le.6) then
            vmx(i) = vm(i) * T_SCALAR(t,'v',ttype,ftToptV,ftHaV,ftHdV,
     &sum(ce_t(:))/30.d0)
            jmx(i) = jm(i) * T_SCALAR(t,'j',ttype,ftToptJ,ftHaJ,ftHdJ,
     &sum(ce_t(:))/30.d0)

            !leaf age scale Vcmax & Jmax
            vmx(i) = vmx(i)*npp_eff
            jmx(i) = jmx(i)*npp_eff
          
            !water limitation scale Vcmax & Jmax
            vmx(i) = vmx(i)*kg**p_nu3
            jmx(i) = jmx(i)*kg**p_nu3

            !convert from umol to mol
            vmx(i) = vmx(i) * 1d-6
            jmx(i) = jmx(i) * 1d-6

          ENDIF


        !water limitation IF
        ENDIF


        !leaf respiration - nighttime 
        if((vcmax_type.eq.2).or.(vcmax_type.eq.3).or.(vcmax_type.eq.5) 
     &.or.(vcmax_type.eq.6).or.(vcmax_type.eq.7) ) 
     &then
          !these vcmax types assume no leaf N and that day AND night leaf resp = f(vcmax)
          nleaf(i) = 0.d0
          nleaf_sum= 0.d0
          dresp(i) = vmx(i)*p_rd ! * nightr:dayr scalar
        else
          dresp(i)  = nleaf(i)*p_dresp
        endif 
        drespt(i) = dresp(i)         ! I'm not sure what drespt does - nothing I think

        !sum to get canopy N & Respiration
        IF (i.LT.lai) THEN
          nleaf_sum = nleaf_sum + nleaf(i)
          canres    = canres + dresp(i)
        ELSE
          nleaf_sum = nleaf_sum + nleaf(i)*rem
          canres    = canres + dresp(i)*rem
        ENDIF
        
      !initial LAI loop
      ENDDO

      !else there is no light or leaves  
      ELSE
        
        IF(vcmax_type.eq.9) THEN
          !LUNA model retains the final value of vcmax pre-total leaf loss to initialise the following year
          vm(:)     = vcmax(:) 
          jm(:)     = jmax(:)
        ELSE
          vm(:)     = 0.0d0
          jm(:)     = 0.0d0
        ENDIF

        can(:)     = 0.0d0
        upt(:)     = 0.0d0
        nleaf(:)   = 0.0d0
        dresp(:)   = 0.0d0
        drespt(:)  = 0.0d0
        jmx(:)     = 0.0d0
        vmx(:)     = 0.0d0
        jshade(:)  = 0.0d0
        jsunlit(:) = 0.0d0
        fsunlit(:) = 0.0d0
        fshade(:)  = 0.0d0
        qsunlit(:) = 0.0d0
        qshade(:)  = 0.0d0
        
        canres   = 0.0d0
        !vcmax    = 0.0d0
        !jmax     = 0.0d0
        leaf_nit = 0.0d0        
 
      !end light and leaves loop
      ENDIF
      
      !output varibles - vcmax & jmax at 25oC (all layers) & leaf N (top layer only), all per unit leaf area, without water stress scaling
      !print*, vcmax(1),vm(1) 
      vcmax    = vm(:)
      jmax     = jm(:)
      leaf_nit = nleaf(1)
      !print*, vcmax, jmax  
      !print'(12(f8.4,1x))', vmx(:) * 1d6 
      !print'(12(f8.4,1x))', vm(:)  
      

*----------------------------------------------------------------------*
* Assimilation calculations using subroutines ASSVMAX and ASSJ.        *
*----------------------------------------------------------------------*

      !zero layer arrays 
      a(:)  = 0.0d0
      ci(:) = 0.0d0
      gs(:) = 0.0d0

       !if(output.and.(mnth.eq.omnth).and.(day.eq.oday)) then
        ! print*, ''
        ! print*, t,rn,soil2g,wtwp
        ! print*, ''
       !endif
      !if there are leaves
      IF (rlai.GT.0.1d0) THEN
        !leaf boundary layer conductance 
        !- this needs attention, this is not right to use canga divided by lai
        ga = canga/(8.3144d0*tk/p)/rlai

        !subd_par if statement mean daily PAR or downscaled sub-daily PAR 
        ! - if a daily mean PAR is used (SDGVM default) use above calculated PAR values
        IF(subd_par.eq.0) THEN
          if(output.and.(mnth.eq.omnth).and.(day.eq.oday)) then
            print*, 'calculate assimilation based on mean daily PAR'
          endif

          !assimilation LAI loop - fshade, fsunlit, jshade, jsunlit have been calculated above 
          DO i=1,lai
            CALL ASSIMILATION_CALC(fshade(i),fsunlit(i),
     &vmx(i),xvmax,rd(i),jshade(i),jsunlit(i),upt(i),qshade(i),
     &qsunlit(i),t,rn,soil2g,wtwp,ga,rh,C3,kg,ko,kc,tau,p,oi,ca,
     &a(i),gs(i),ci(i),day,mnth,oday,omnth,.FALSE.,gs_func,ftg0,ftg1,i)
          ENDDO

          if(output.and.(mnth.eq.omnth).and.(day.eq.oday)
     &) then
            print*, a(1),rd(1),gs(1),ci(1),vmx(1)
            print*, ''
         endif

        !if switch 1 not equal to 0: use sub-daily PAR loop
        ELSE

          !scaling factor for integration of subdaily variability 
          sd_scale(:) = 1.0d0
          sd_scale(1) = 2.0d0
          sd_scale(par_loops+1) = 2.0d0
          sd_scale2   = 1.0d0 / real(par_loops)

          !zero arrays
          a  = 0.0
          ci = 0.0
          gs = 0.0
          fsunlit = 0.0
          fshade  = 0.0
          qsunlit = 0.0
          qshade  = 0.0
          qdirect = 0.0
          qdiff   = 0.0        

          a_sd  = 0.0
          ci_sd = 0.0
          gs_sd = 0.0
          fsunlit_sd = 0.0
          fshade_sd  = 0.0
          qsunlit_sd = 0.0
          qshade_sd  = 0.0
          qdirect_sd = 0.0 
          qdiff_sd   = 0.0        
        
          DO ii=1,par_loops+1
            !print*, '' 
            !print*, par_loops
            !print*, 'par_loop:',ii
            CALL PFD(lat,no_day(year,mnth,day,thty_dys),hrs,cld,
     &qdirect_sd,qdiff_sd,q_sd,swr,read_par,subd_par,(ii-1),
     &par_loops,calc_zen,cos_zen)
            
            ! at dawn assume diffuse light is 5 umol/m2/s
            ! and cos_zen a fraction above zero
            IF(qdirect_sd+qdiff_sd.lt.5.0d-6) qdiff_sd = 5.0d-6  
            IF(cos_zen.le.0d0)                cos_zen  = cos(3.141/2)

            DO i=1,lai
              !in the average PAR version of the model this is contained within a water limitation if statement 
              CALL GOUDRIAANSLAW(real(i)-0.5d0,rlai,qdirect_sd,
     &qdiff_sd,fsunlit_sd(i),qsunlit_sd(i),fshade_sd(i),
     &qshade_sd(i),can_clump,cos_zen,s070607,gold)

              !if(output.and.(mnth.eq.omnth).and.(day.eq.oday)
      !&.and.(ii.eq.4)) then
                !print*, qdirect_sd+qdiff_sd
                !print*, qsunlit_sd(i)+qshade_sd(i)
      !&qsunlit_sd(i,ii)+qshade_sd(i,ii),i,ii
              !endif

              !calculate jshade and jsunlit
              CALL JCALC(qsunlit_sd(i),jmx(i),hw_j,
     &jsunlit_sd(i))
              CALL JCALC(qshade_sd(i),jmx(i),hw_j,jshade_sd(i))
              !calculate assimilation
              CALL ASSIMILATION_CALC(fshade_sd(i),fsunlit_sd(i),
     &vmx(i),xvmax,rd(i),jshade_sd(i),jsunlit_sd(i),upt(i),
     &qshade_sd(i),qsunlit_sd(i),t,rn,soil2g,wtwp,ga,rh,C3,kg,ko,
     &kc,tau,p,oi,ca,a_sd(i),gs_sd(i),ci_sd(i),day,mnth,oday,
     &omnth,.FALSE.,gs_func,ftg0,ftg1,i)

            !LAI loop
            ENDDO

            !if(output.and.(mnth.eq.omnth).and.(day.eq.oday)
      !&.and.(ii.eq.4)) then
            !  print*, a_sd(1), ii
            !  print*, ''
            !endif

          !numerically integrate canopy values to get daily mean (mean is necessary because per second values are scaled to daytime values in doly)   
            a(:)       =  a(:)       + a_sd(:)       / sd_scale(ii) 
            ci(:)      =  ci(:)      + ci_sd(:)      / sd_scale(ii)  
            gs(:)      =  gs(:)      + gs_sd(:)      / sd_scale(ii)  
            fsunlit(:) =  fsunlit(:) + fsunlit_sd(:) / sd_scale(ii)  
            fshade(:)  =  fshade(:)  + fshade_sd(:)  / sd_scale(ii)  
            qsunlit(:) =  qsunlit(:) + qsunlit_sd(:) / sd_scale(ii)  
            qshade(:)  =  qshade(:)  + qshade_sd(:)  / sd_scale(ii)  
            qdirect    =  qdirect    + qdirect_sd    / sd_scale(ii)  
            qdiff      =  qdiff      + qdiff_sd      / sd_scale(ii)  

            if(ii.eq.1) maxlight(:) = 
     &fsunlit_sd(:)*qsunlit_sd(:)+fshade_sd(:)*qshade_sd(:)

          !sub-daily loop
          ENDDO

          a       =  a        * sd_scale2  
          ci      =  ci       * sd_scale2  
          gs      =  gs       * sd_scale2  
          fsunlit =  fsunlit  * sd_scale2
          fshade  =  fshade   * sd_scale2 
          qsunlit =  qsunlit  * sd_scale2
          qshade  =  qshade   * sd_scale2
          qdirect =  qdirect  * sd_scale2
          qdiff   =  qdiff    * sd_scale2 

        !subd_par if statement mean daily PAR or downscaled sub-daily PAR 
        ENDIF
      
      ! calculate environment variables for Maire & LUNA vcmax calc
      ce_t(1:29)  = ce_t(2:30)
      ce_rh(1:29) = ce_rh(2:30)
      do i=1,lai
        ce_ci(1:29,i)       = ce_ci(2:30,i) 
        ce_ga(1:29,i)       = ce_ga(2:30,i) 
        ce_light(1:29,i)    = ce_light(2:30,i)
        ce_maxlight(1:29,i) = ce_maxlight(2:30,i)
      enddo

      ce_t(30)       = t
      ce_rh(30)      = rh
      ce_ci(30,:)    = ci(:)
      ce_ga(30,:)    = ga
      ce_light(30,:) = fsunlit(:)*qsunlit(:)+fshade(:)*qshade(:)    
      ce_maxlight(30,:) = maxlight 
      !print*, ce_ci(30,:)

      !rlai if
      ELSE
        dresp(:) = 0.0d0
        ci(:)    = ca
      !End rlai if
      ENDIF

*      if(day.eq.160) print*, qdirect + qdiff
 
      !sum canopy layers to get canopy integrated values
      can2a = 0.0d0
      can2g = 0.0d0
      canrd = 0.0d0
      apar  = 0.0d0
      DO i=1,lai
        IF (i.LT.lai) THEN
          can2a = can2a + a(i)
          can2g = can2g + gs(i)
          canrd = canrd + rd(i)
          apar  = apar  + fsunlit(i)*qsunlit(i)+fshade(i)*qshade(i)
        ELSE
          can2a = can2a + a(i)*rem
          can2g = can2g + gs(i)*rem
          canrd = canrd + rd(i)*rem
          apar  = apar + rem*(fsunlit(i)*qsunlit(i)+fshade(i)*qshade(i))
        ENDIF
      ENDDO
      !print*, can2a,canrd,canres,can2g

      !calculate fapar
      fpr    = apar/(qdiff+qdirect)

      !topleaf values
      tassim =  a(1)
      tgs    = gs(1)
      tci    = ci(1)

      ! Cumulate daily assimilation of the last lai layer, 'suma'
      IF (lai.GT.1) THEN
        suma = suma + (rem*a(lai) + (1.0d0 - rem)*a(lai - 1))*
     &3600.0d0*hrs - (rem*dresp(lai) + (1.0d0 - rem)*dresp(lai -
     &1))*3600.0d0*(24.0d0 - hrs)
        sumd = sumd + (rem*dresp(lai) + (1.0d0 - rem)*dresp(lai -
     &1))*3600.0d0*(24.0d0 - hrs)
      ELSE
!            suma = suma + rem*a(lai)*3600.0d0*hrs - rem*dresp(lai)*
!     &3600.0d0*(24.0d0 - hrs)
        sumd = sumd + rem*dresp(lai)*3600.0d0*(24.0d0 - hrs)
        suma = 0.0d0
      ENDIF
        
      if(output.and.(mnth.eq.omnth).and.(day.eq.oday)) then
        print*, '\n integrated daily means for each layer'
        print*,  a(:)
        print*, rd(:) * 1d6
        print*, gs(:)
        print*, ci(:)
      endif
      
      gsm = gs(1)
 
      !convert suma to mols
      suma = suma*1.d-6

      RETURN
      END

*----------------------------------------------------------------------*
*      electron transport rate subroutine                              *
*      extracted from main nppcalc subroutine APW 30 Sept 2013         *
*----------------------------------------------------------------------*
      SUBROUTINE JCALC(q,jmx,fw1984,j)

      IMPLICIT NONE

      real*8  q,jmx,j,I2,qa,qb,qc
      integer fw1984

      IF(fw1984.eq.0) THEN
         !use Harley 1992 J to Jmax realtionship
         j=0.24d0*q/(1.0d0+(0.24d0**2)*(q**2)/(jmx**2))**0.5d0
      ELSE
         !use Farquhar & Wong 1984 J to Jmax realtionship
         I2 = q * 0.36d0    
         qa = 0.7d0
         qb = I2+jmx
         qc = I2*jmx
         j  = (qb - (qb**2 - 4*qa*qc)**0.5d0) / (2*qa)
      ENDIF
      END SUBROUTINE

*----------------------------------------------------------------------*
*      assimilation subroutine                                         *
*      extracted from main nppcalc subroutine APW 30 Sept 2013         *
*      LAI layer arrays have been removed from this subroutine         *
*      this subroutine is now called from within an LAI loop           *
*----------------------------------------------------------------------*
      SUBROUTINE ASSIMILATION_CALC(fshade,fsunlit,vmx,xvmax,rd,
     &jshade,jsunlit,upt,qshade,qsunlit,t,rn,soil2g,wtwp,ga,rh,C3,kg,
     &ko,kc,tau,p,oi,ca,a,gs,ci,day,mnth,oday,omnth,output,gs_func,g0,
     &g1,i)

      IMPLICIT NONE

      !input variable by layer variables
      real*8  fshade,fsunlit,vmx,rd,jshade,jsunlit,upt
      real*8  qshade,qsunlit

      !input variables constant over layers
      real*8   t,rn,soil2g,wtwp,ga,rh,kg,ko,kc,tau,p,oi,ca,g0,g1
      integer  C3,day,mnth,oday,omnth,gs_func,i

      !internal variables
      real*8  tpav,tppcv,cs,tpgsv,brent_solver,gs_leaf 
      real*8  tpaj,tppcj,tpgsj,dv
      real*8  tpac4,tppcc4,tpgsc4
      real*8  ashade,pcshade,gsshade
      real*8  asunlit,pcsunlit,gssunlit

      !output variables
      real*8  a,gs,ci
      real*8  xvmax      !need to make this an array and set it up in canopy properties lai loop
      logical output,brent     !flag to print troubleshooting output to screen

      brent = .false. ! if the brent solver needs reimplementing here 
                      ! then it will be necessary to pass the pft temp scalar parameters to this subroutine  
      
      !calculate VPD
      dv = 0.6108d0*exp(17.269d0*t/(237.3d0 + t))*(1.0d0 - 
     &rh/100.0d0)

      ! print*, kg

      IF ((t.GT.0.0d0).AND.(rn.GT.0.0d0).AND.(soil2g.GT.wtwp).AND.
     &(fshade+fsunlit.GT.0.1d0)) THEN

          a  = 0.0d0
          ci = 0.0d0
          gs = 0.0d0

          !C3/C4 if
          IF (c3.EQ.1) THEN
       !     if(.not.brent) then
              CALL ASSVMAX(tpav,tppcv,cs,tpgsv,oi,ca,vmx,
     &rh,kg,ko,kc,tau,rd,ga,t,p,gs_func,g0,g1,dv,i)
       !     else
       !       tpav  = brent_solver(1,0.d0,50.d0,
      !&oi,ca,vmx,0.d0,rh,kg,rd,ga,t,p,gs_func,g0,g1,dv,i,0.d0,0.d0,ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ,sum(ce_t(:))/30.d0)
       !       cs    = ca - tpav*1.0e-6*p/ga
       !       tpgsv = GS_LEAF(g0,g1,tpav,dv,t,cs,gs_func) * kg
       !       tppcv = cs - tpav*1.0e-6*p/tpgsv
       !     endif

            !light or shade leaves
            IF (fshade.GT.0.01d0) THEN
        !      if(.not.brent) then
                CALL ASSJ(tpaj,tppcj,tpgsj,oi,ca,rh,kg,tau,rd,ga,t,p,
     &jshade,gs_func,g0,g1,dv,i)
        !      else
        !        tpaj  = brent_solver(2,0.d0,40.d0,
      !&oi,ca,vmx,jshade,rh,kg,rd,ga,t,p,gs_func,g0,g1,dv,i,0.d0,0.d0,ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ,sum(ce_t(:))/30.d0)
        !        cs    = ca - tpaj*1.0e-6*p/ga
        !        tpgsj = GS_LEAF(g0,g1,tpaj,dv,t,cs,gs_func) * kg
        !        tppcj = cs - tpaj*1.0e-6*p/tpgsj
        !      endif

              CALL LIMITATIONPROCESS(tpav,tppcv,tpgsv,tpaj,tppcj,
     &tpgsj,ashade,pcshade,gsshade)
              a  = fshade*ashade
              ci = fshade*pcshade
              gs = fshade*gsshade
              !if(i.eq.1) print*, 'C3 sun:',a,gs,ci,ga,cs
            ENDIF

            IF (fsunlit.GT.0.01d0) THEN
        !      if(.not.brent) then
                CALL ASSJ(tpaj,tppcj,tpgsj,oi,ca,rh,kg,tau,rd,ga,t,p,
     &jsunlit,gs_func,g0,g1,dv,i)
        !      else
        !        tpaj  = brent_solver(2,0.d0,40.d0,
      !&oi,ca,vmx,jsunlit,rh,kg,rd,ga,t,p,gs_func,g0,g1,dv,i,0.d0,0.d0,ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ,sum(ce_t(:))/30.d0)
        !        cs    = ca - tpaj*1.0e-6*p/ga
        !        tpgsj = GS_LEAF(g0,g1,tpaj,dv,t,cs,gs_func) * kg
        !        tppcj = cs - tpaj*1.0e-6*p/tpgsj
        !      endif

              CALL LIMITATIONPROCESS(tpav,tppcv,tpgsv,tpaj,tppcj,
     &tpgsj,asunlit,pcsunlit,gssunlit)
              a  = a  + fsunlit*asunlit
              ci = ci + fsunlit*pcsunlit
              gs = gs + fsunlit*gssunlit
              !if(i.eq.1) print*, 'C3 sun:',a,gs,ci,ga,cs
            ENDIF

          !C3/C4 if
          ELSE
            !print*, rd
            !sunlit or shade if 
            IF (fshade.GT.0.01d0) THEN
              CALL ASSC4_colim(ashade,pcshade,gsshade,ca,vmx,rh,kg,t,
     &qshade,p,upt,rd,dv,g0,g1,gs_func,ga)
              a  = fshade*ashade
              ci = fshade*pcshade
              gs = fshade*gsshade
              !print*, 'C4 shade:',a,gs,ci,ga,cs
            ENDIF
            IF (fsunlit.GT.0.01d0) THEN
              CALL ASSC4_colim(asunlit,pcsunlit,gssunlit,ca,vmx,rh,kg,t,
     &qsunlit,p,upt,rd,dv,g0,g1,gs_func,ga)
              a  = a  + fsunlit*asunlit
              ci = ci + fsunlit*pcsunlit
              gs = gs + fsunlit*gssunlit
              !if(i.eq.1) print*, 'C4 sun:',a,gs,ci,ga,cs
            ENDIF

            IF (gs .LT. 0.0d0)  gs = 0.0d0

          !C3/C4 if
          ENDIF


        !assimilation if
        ELSE
          a  = 0.0d0
          gs = g0 
          ci = ca

        !end assimilation if
        ENDIF

      !if(output.and.(mnth.eq.omnth).and.(day.eq.oday)) then
      !  print*, qsunlit+qshade,a,gs,ci
      !  print*, ''
      !endif

      END SUBROUTINE


*----------------------------------------------------------------------*
*     determin which process light, or Rubisco is the limiting         *
*     factor and return the assimilation, partial presure and          *
*     stomatal conductance corresponding to the limiting process       *
*----------------------------------------------------------------------*


      SUBROUTINE LIMITATIONPROCESS(av,pcv,gsv,aj,pcj,gsj,a,pc,gs)
      REAL*8 av,pcv,gsv,aj,pcj,gsj,a,pc,gs

      IF (av.LT.aj) THEN
         a = av
         gs = gsv
         pc = pcv
      ELSE
         a = aj
         gs = gsj
         pc = pcj
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSVMAX                          *
*                          ******************                          *
*                                                                      *
* ASSVMAX computes the assimilation rate by forming a cubic equation   *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. The variables here *
* correspond to the Doly paper, in the main program they are such that *
*        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           *
*     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSVMAX(a,ci,cs,gs,oi,ca,vmx,rh,kg,ko,kc,tau,rd,ga,t,p,
     &gs_func,g0,g1,dv,i)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,ci,gs,oi,ca,vmx,rh,ko,kc,tau,rd,fa0,fa00,faf,gam,step
      REAL*8 p,ga,x,y,dv,t,a0,cs,af,fa,kg,ax,bx,cx,g0,g1!,gsmin
      REAL*8 gs_leaf
      INTEGER i,gs_func

!      gsmin = 0.005d0

!      x = 3.0d0 + 0.2d0*t
!      if (x.LT.0.25d0) x = 0.25d0
!      y = 1.5d0

!      dv = 0.6108d0*exp(17.269d0*t/(237.3d0 + t))*(1.0d0 - 
!     &rh/100.0d0)

      !if (5.gt.30) then
      !do i=0,200
      !a = real(i)/10.0d0
      !cs = ca - a*1.0e-6*p/ga
      !gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
      !&1.54d0*t))*kg
      !ci = cs - a*1.0e-6*p/gs
      !gam = 1.0d0 - 0.5d0*oi/tau/ci
      !fa0 = a - (gam*vmx*ci/(ci + kc*(1.0d0 + oi/ko)) -
      !&rd)*1.0e6
      !print'(7f10.3)',a,fa0,cs,gs,ci,gam,(x*a/(1.0d0 + dv/y)/(cs*10.0d0
      !& - 1.54d0*t))
      !enddo
      !endif

      step = 1.0d0
      a = 0.0d0
      cs = ca - a*1.0e-6*p/ga
!      gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
      gs =  GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg

      ci = cs - a*1.0e-6*p/gs
      gam = 1.0d0 - 0.5d0*oi/tau/ci
      fa0 = a - (gam*vmx*ci/(ci + kc*(1.0d0 + oi/ko)) -
     &rd)*1.0e6

      a0 = a
      fa00 = fa0
      !if(i.eq.1) print*, vmx,gs,ci,rd,gam,kc,ko

50    CONTINUE
        a = a + step
        cs = ca - a*1.0e-6*p/ga
!        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
        gs =  GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
        faf = a - (gam*vmx*ci/(ci + kc*(1.0d0 + oi/ko)) -
     &rd)*1.0e6
      IF (faf.LT.0.0d0) THEN
        fa0 = faf
        a0 = a
        GOTO 50
      ENDIF
      af = a
      !if(i.eq.1) print*, 'a =',0.0,'fa =',fa00,'b =',a,'fb =',faf 

      IF (af.GT.100d0) THEN
        a = 0.0d0
        cs = ca - a*1.0e-6*p/ga
!        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
        gs =  GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
      ELSE
        a = (a0 + af)/2.0d0
        cs = ca - a*1.0e-6*p/ga
!        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
        gs =  GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
        fa = a - (gam*vmx*ci/(ci + kc*(1.0d0 + oi/ko)) -
     &rd)*1.0e6

        bx = ((a+a0)*(faf-fa)/(af-a)-(af+a)*(fa-fa0)/(a-a0))/(a0-af)
        ax = (faf-fa-bx*(af-a))/(af**2-a**2)
        cx = fa-bx*a-ax*a**2

        IF (abs(ax).GT.0.0d0) THEN
          IF (bx**2-4.0d0*ax*cx.GT.0.0d0) THEN
            a = (-bx+(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
        IF (a.GT.af)  a = (-bx-(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
          ELSE
            a = 0.0d0
          ENDIF
        ELSE
          IF (abs(bx).GT.0.0d0) THEN
            a =-cx/bx
          ELSE
            a = 0.0d0
          ENDIF
        ENDIF

      ENDIF

        !if(i.eq.1) print*, 'a =',0.0,'fa =',fa00,'assim=',a
        cs = ca - a*1.0e-6*p/ga
!        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
        gs =  GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
*      print'(6f10.3)',a,fa0,cs,gs,ci,gam
*      print*

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSJ                             *
*                          ***************                             *
*                                                                      *
* ASSVMAX computes the assimilation rate by forming a cubic equation   *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. The variables here *
* correspond to the Doly paper, in the main program they are such that *
*        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           *
*     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSJ(a,ci,gs,oi,ca,rh,kg,tau,rd,ga,t,p,j,gs_func,g0,g1,
     &dv,i)
*----------------------------------------------------------------------*
      IMPLICIT NONE

      REAL*8 a,ci,gs,oi,ca,rh,tau,rd,j,gam,fa0,fa00,faf,a0,af,step
      REAL*8 p,ga,x,y,dv,t,cs,fa,kg
      REAL*8 ax,bx,cx,gs_leaf,g0,g1
      INTEGER gs_func,i

!      gsmin = 0.005d0

!      x = 3.0d0 + 0.2d0*t
!      if (x.LT.0.25d0) x = 0.25d0
!      y = 1.5

!      dv = 0.6108d0*exp(17.269d0*t/(237.3d0 + t))*(1.1d0 - 
!     &rh/100.0d0)

      step = 1.0d0
      a = 0.0d0
      cs = ca - a*1.0e-6*p/ga
!      gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
      gs =  GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
      ci = cs - a*1.0e-6*p/gs
      gam = 1.0d0 - 0.5d0*oi/tau/ci
      fa0 = a - (gam*j*ci/4.0d0/(ci + oi/tau) -
     &rd)*1.0e6

      a0 = a
      fa00 = fa0
      !if(i.eq.1) print*, gam,j,ci,rd

50    CONTINUE
        a = a + step
        cs = ca - a*1.0e-6*p/ga
!        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
        gs =  GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
        faf = a - (gam*j*ci/4.0d0/(ci + oi/tau) -
     &rd)*1.0e6
      IF (faf.LT.0.0d0) THEN
        fa0 = faf
        a0 = a
        GOTO 50
      ENDIF
      af = a
      
      !if(i.eq.1) print*, 'a =',0.0,'fa =',fa00,'b =',a,'fb =',faf 

      IF (af.GT.100d0) THEN
        a = 0.0d0
        cs = ca - a*1.0e-6*p/ga
!        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
        gs =  GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
      ELSE

        a = (a0 + af)/2.0d0
        cs = ca - a*1.0e-6*p/ga
!        gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
        gs =  GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci = cs - a*1.0e-6*p/gs
        gam = 1.0d0 - 0.5d0*oi/tau/ci
        fa = a - (gam*j*ci/4.0d0/(ci + oi/tau) -
     &rd)*1.0e6

        bx = ((a+a0)*(faf-fa)/(af-a)-(af+a)*(fa-fa0)/(a-a0))/(a0-af)
        ax = (faf-fa-bx*(af-a))/(af**2-a**2)
        cx = fa-bx*a-ax*a**2

        IF (abs(ax).GT.0.0d0) THEN
          IF (bx**2-4.0d0*ax*cx.GT.0.0d0) THEN
            a = (-bx+(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
        IF (a.GT.af)  a = (-bx-(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
          ELSE
            a = 0.0d0
          ENDIF
        ELSE
          IF (abs(bx).GT.0.0d0) THEN
            a =-cx/bx
          ELSE
            a = 0.0d0
          ENDIF
        ENDIF

      ENDIF

      !if(i.eq.1) print*, 'a =',0.0,'fa =',fa00,'assim=',a
 
      cs = ca - a*1.0e-6*p/ga
!      gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
      gs =  GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
      ci = cs - a*1.0e-6*p/gs
      gam = 1.0d0 - 0.5d0*oi/tau/ci

      RETURN
      END

*----------------------------------------------------------------------*
      FUNCTION FASSV(a,gstar,km,ca,vmx,kg,rd,ga,t,p,gs_func,g0,g1,dv,i)

      IMPLICIT NONE
      REAL*8 a,ci,gs,ca,vmx,rd,gam
      REAL*8 p,ga,dv,t,a0,cs,kg,g0,g1
      REAL*8 gs_leaf,fassv,gstar,km
      INTEGER gs_func,i

      cs  = ca - a*1.0e-6*p/ga
      gs  = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
      ci  = cs - a*1.0e-6*p/gs
      gam = 1.0d0 - gstar/ci
      fassv = (gam*vmx*ci/(ci + km) - rd)*1.0e6

      !if(i.eq.1) print*, vmx,gs,ci,rd,t,gam,kc,ko

      END
        
        
*----------------------------------------------------------------------*
        
      FUNCTION FASSJ(a,gstar,ca,kg,rd,ga,t,p,j,gs_func,g0,g1,dv,i)

      IMPLICIT NONE
      REAL*8 a,ci,gs,ca,rd,j,gam
      REAL*8 p,ga,dv,t,cs,kg,gstar
      REAL*8 gs_leaf,g0,g1,fassj
      INTEGER gs_func,i

      cs  = ca - a*1.0e-6*p/ga
      gs  = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
      ci  = cs - a*1.0e-6*p/gs
      gam = 1.0d0 - gstar/ci
      fassj = (gam*j*ci/4.0d0/(ci + 2*gstar) - rd)*1.0e6
      END
*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION GS                                 *
*                          ****************                            *
*    Gives stomatal conductance in mol m-2 s-1                         *  
*    GS_WOOD, the default gs model is a form of the Leuning 1995 model *
*    assuming gstar = 1.54t, this is inconsistent with the gstar in the* 
*    assimilations calculations                                        *
*----------------------------------------------------------------------*

      FUNCTION GS_LEAF(g0,g1,a,dv,t,cs,gs_func)

      IMPLICIT NONE

      REAL*8 :: gs_leaf,gs_wood,gs_med
      REAL*8 :: g0,g1,a,dv,t,cs
      INTEGER:: gs_func

      if(gs_func.eq.0) then
        gs_leaf = GS_WOOD(g0,a,dv,t,cs) 
      elseif(gs_func.eq.1) then
        gs_leaf = GS_MED(g0,a,dv,g1,cs) 
      else
        print*, 'No gs function defined. STOPPING.'
        print*, 'gs function:',gs_func
        stop
      endif
      
      END      
*----------------------------------------------------------------------*
*    GS_WOOD, the default gs model is a form of the Leuning 1995 model *
*    assuming gstar = 1.54t, this is inconsistent with the gstar in the* 
*    assimilations calculations                                        *
*----------------------------------------------------------------------*
      
      FUNCTION GS_WOOD(g0,a,dv,t,cs)

      IMPLICIT NONE

      real*8 :: gs_wood,g0,a,dv,t,cs
      real*8 :: x
      
      x = 3.0d0 + 0.2d0*t
      if (x.LT.0.25d0) x = 0.25d0

!      gs = gsmin + (x*a/(1.0d0 + dv/y)/(cs*10.0d0 - 
!     &1.54d0*t))*kg
      gs_wood = g0 + (x*a/(1.0d0 + dv/1.5d0)/(cs*10.0d0 - 
     &1.54d0*t))
      
      END

*----------------------------------------------------------------------*

      FUNCTION GS_MED(g0,a,dv,g1,cs)

      IMPLICIT NONE

      real*8 :: gs_med,g0,a,dv,g1,cs
      
      gs_med = g0 + (g1/(1.0d0 + dv**0.5d0))*(a/(cs/0.1013d0))
      
      END

*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE FARQ_PARS                        *
*                          ******************                          *
*     Gives M-M parameters from Farquhar model of photosynthesis       *  
*     Parameters Kc and Ko given in Pa                                 *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE FARQ_PARS(t,kc,ko,tau,alpha,farq_pars_func)

      IMPLICIT NONE
       
      REAL*8  t,tk,kc,ko,tau,alpha
      INTEGER farq_pars_func

      tk = t + 273.15d0  
      
      alpha = 0.24d0
      kc  = exp(35.8d0 - 80.5d0/(0.00831d0*tk))
      ko  = exp(9.6d0 - 14.51d0/(0.00831d0*tk))*1000.0d0
      tau = exp(-3.949d0 + 28.99d0/(0.00831d0*tk))

      END SUBROUTINE
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION T_SCALAR                           *
*                          *****************                           *
*     Gives vcmax or jmax temperature scalar                           *  
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION T_SCALAR(t,jv,ttype,ftTopt,ftkHa,ftkHd,tmonth)

      IMPLICIT NONE

      REAL*8    :: t,qt,t_scalar,ftTopt,ftkHa,ftkHd,ftHa,ftHd,tmonth
      REAL*8    :: Tsk,Trk,R,deltaS,dS
      INTEGER, intent(in) :: ttype 
      CHARACTER :: jv

      !print*, 'T_SCALAR, ttype:', ttype
      !if(ttype.eq.0) print*, 'here, ttype:', ttype

      !if(0.eq.0) print*, '0 equals 0'
      
      if(ttype.eq.0) then
      !if(0.eq.0) then
        ! SDGVM original
        !print*, 'ttype original'
        qt = 2.3d0
        if (t.gt.30.0d0) t = 30.d0
        t_scalar = qt**(t/10.0d0)/qt**(2.5d0)
   
      else if(ttype.ge.1) then
        ! modified Arrhenius
        
        R = 8.31446        

        ! convert from kJ to J
        ! parameter values are read in as kJ... to save space
        ftHa = ftkHa * 1e3
        ftHd = ftkHd * 1e3

        if(ttype.eq.2) then
          ! employ Kattge&Knorr scaling based on mean temp of previous month 
          
          ftHa = 71513
          ftHd = 49884

          if(jv.eq.'v') then
            deltaS = 668.39d0 - 1.07d0*tmonth
          else
            deltaS = 659.70d0 - 0.75d0*tmonth
          endif
        else 
          ! use fixed delta S based on input PFT parameters
          deltaS = ftHd/(ftTopt+273.15) + ( R*log(ftHa/(ftHd-ftHa)) )
        endif

        Tsk    = t + 273.15
        Trk    = 298.15
  
        ! see Medlyn etal 2002 or Kattge&Knorr 2007  
        t_scalar = exp(ftHa*(Tsk-Trk) / (R*Tsk*Trk)) * ( 
     & (1 + exp((Trk*deltaS-ftHd) / (Trk*R)) ) 
     &/(1 + exp((Tsk*deltaS-ftHd) / (Tsk*R)) ) )  

        !adjust jmax based on kattge&knorr Jmax to Vcmax ratio relationship to temp
        if((jv.eq.'j').and.(ttype.ge.2)) 
     &t_scalar = t_scalar * (2.59-0.035*tmonth)/1.715

      else
        t_scalar = 0.d0
        print*, 'temperature scalar type',ttype,'undefined'
        stop
      endif


      END

*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION VCMAX                              *
*                          ****************                            *
*     Gives vcmax in mol m-2 s-1                                       *  
*     Calculate vcmax from Amax which is in turn calculated from       *
*     N uptake rate, which is a function of soil c and N               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION VCMAX_FROM_AMAX(amx,gs_func,prd,g0,g1,farq_pars_func)

      IMPLICIT NONE
   
      REAL*8 amx,ci,ca,gs,vcmax_from_amax,prd,kc,ko,km,gstar,dv,t,g0,g1
      REAL*8 gs_leaf,tau,alpha
      INTEGER gs_func,farq_pars_func

      ! check units
      dv = 0.75 ! value used in Woodward 1994 Ecophys. of Photosyn.
      ca = 38   ! assumed value used in Woodward 1994
      t  = 25.d0 ! temperature - Woodward 1994 used Amax at optimum temp, for now assumed at 25oC

      gs =  GS_LEAF(g0,g1,amx,dv,t,ca,gs_func)
      ci = ca - amx*0.101325/gs 
     
      CALL FARQ_PARS(t,kc,ko,tau,alpha,farq_pars_func)
      km = kc*(1 + 21000.d0/ko)
      gstar = 0.5d0*21000.d0/tau 

      vcmax_from_amax = amx / ((ci - gstar) / (ci + km) - prd)

      !print'(6f9.4)', vcmax_from_amax,amx,gs,ci,km,gstar
      END
*----------------------------------------------------------------------*
      
*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION MAIRE VCMAX                        *
*                          ********************                        *
*     Gives vcmax AT LEAF TEMPERATURE  in mol m-2 s-1                  *  
*     Calculate vcmax that satisfies the condition that Wc = Wj        *
*     under the past months environmental conditions                   *
*     see Maire et al 2012 PLoS ONE                                    *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION VCMAX_MAIRE(vm,vt,jt,ci,km,gstar,alpha,light,
     &farq_pars_func,i)
  
      IMPLICIT NONE
      REAL*8  :: vcmax_maire,t,ci,oi,light,vm
      REAL*8  :: vt,jt,km,gstar,alpha
      INTEGER :: farq_pars_func,i,ttype

      vcmax_maire = vm*
     &(1.d0+(alpha*light/(jt*(exp(1.d0)*(vm/vt)**0.89d0)))**2.d0)**0.5d0
     &*(4*ci+8*gstar) - alpha*light*(ci+km) 

c      vcmax_maire = alpha*light / 
c     &(1.d0+(alpha*light/(jt*(exp(1.d0)*(vm/vt)**0.89d0)))**2.d0)**0.5d0
c     &*(ci+km)/(4*ci+8*gstar)  
c      if(i.eq.1) print*, 'maire vcmax calc'
c      if(i.eq.1) print'(4f8.4)', vt,jt,km,gstar
c      if(i.eq.1) print'(3f14.8)', vcmax_maire,ci,light
      END  

*----------------------------------------------------------------------*
   
*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION AMAX                               *
*                          ****************                            *
*     Gives Amax from N uptake rate in umol m-2 s-1                    *  
*     See Woodward etal 1995                                           *
*----------------------------------------------------------------------*
      FUNCTION F_AMAX(soilc,soiln,t)

      IMPLICIT NONE
      REAL*8 f_amax,nup_rate,nup_rate1,nup_rate2,nup_tadj
      REAL*8 soilc,soiln,p_nu2,t
      
      p_nu2 = -8e-5

      nup_rate1 = NUP_RATE(soilc,soiln,120.d0,p_nu2,1/600.d0)
      nup_rate2 = NUP_TADJ(nup_rate1,t,soilc)

      !print*, nup_rate1,nup_rate2

      f_amax = 190.d0 * nup_rate2 / (360.d0+nup_rate2)

      END
*----------------------------------------------------------------------*
      
*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION NUP_TADJ                           *
*                          ****************                            *
*     Gives temperature adjusted n uptake rate                         *  
*     See Woodward etal 1995                                           *
*----------------------------------------------------------------------*
      FUNCTION NUP_TADJ(nup_rate,t,soilc)

      IMPLICIT NONE
      REAL*8 nup_rate,soilc,t,nup_tadj
      REAL*8 kt,u1,u2,u3,tk
 
      tk = t + 273.15d0
      u1 = 40.8d0 + 0.01d0*t - 0.002d0*t**2
      u2 = 0.738d0 - 0.002d0*t
      u3 = 97.412d0 - 2.504d0*log(nup_rate)
      if((soilc.gt.13000.d0).and.(t.lt.15.d0)) then
        kt = (1 + (15-t)/30) * (1 + (soilc - 1.3e4)/1e4)
      else
        kt = 1
      endif

      nup_tadj = kt * exp(u1 - u3/(8.31e-3*tk)) / 
     &(1 + exp((u2 * tk - 205.9)/8.31e-3*tk))
    
      END      
*----------------------------------------------------------------------*
      
*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION NUP_RATE                           *
*                          ****************                            *
*     Gives n uptake rate in umol g- 1 m-2 s-1                         *  
*     See Woodward etal 1995 & 1994                                    *
*----------------------------------------------------------------------*
      FUNCTION NUP_RATE(soilc,soiln,p_nu1,p_nu2,p_nu4)

      IMPLICIT NONE
      REAL*8 nup_rate,nupw,nmult,soilc,soiln,p_nu1,p_nu2,p_nu4
      
      ! Nitrogen uptake - occurs on a daily basis and is allocated on a daily basis
      !                 - i.e. no N storage from day to day, this is a bit weird
      nupw  = p_nu1*(exp(p_nu2*soilc))
      IF (nupw.LT.0.0d0) nupw = 0.0d0

      nmult = soiln*p_nu4
      IF (nmult.GE.1.0d0) nmult = 1.0d0

      nup_rate = nupw*nmult
      if (nup_rate.LT.0.0d0) nup_rate = 0.0d0

      END
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSC4                            *
*                          ****************                            *
*                                                                      *
* ASS4 computes the assimilation rate by solving the 
* equations defining assimilation rate, stomatal conductance,          *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate.This is then       
* compared with assimilation limited to light dependant and vmax       *
* dependant assimilation and the minimum one is taken. 
*                                                                      *
*     Collatz etal 1992                                                *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSC4(a,ci,gs,ca,vmax,rh,kg,t,q,p,up,rd,dv,g0,g1,
     &gs_func,ga)
*----------------------------------------------------------------------*
      IMPLICIT NONE

      REAL*8 gs_leaf,c4_jc
      REAL*8 a,pc,gs,pa,vmax,rh,t,q,up,p,vmq,absx,qrub,f,jqq,dv
      REAL*8 amxt,ac4t,kg,rd,ca,cs,ci,jc,je,k,kt
      REAL*8 fa,fa0,fa00,faf,a0,af,step
      REAL*8 ax,bx,cx,g0,g1,ga,maxa
      INTEGER gs_func

      ! light limited rate (Ji)
      absx = 0.8d0
      qrub = 0.11d0
      f    = 0.6d0
      jqq  = absx*qrub*f*q*1.0d6 - rd*1.0d6

      ! RuBisCO limited rate (Je)
      je = vmax*1.0e6 - rd*1.0d6

      ! Determine max possible assimilation rate 
      maxa = jqq
      if(a.gt.je) maxa = je

      ! low CO2 flux rate (Jc)
      k     = 0.70d0
      kt    = k * 2**( (t-25)/10 ) 

      step = 1.0d0
      a    = 0.0d0
      cs   = ca - a*1.0e-6*p/ga
      gs   = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
      ci   = cs - a*1.0e-6*p/gs
      fa0  = a - c4_jc(ci,p,rd,kt)
      a0   = a
      fa00 = fa0

50    CONTINUE
        a   = a + step
        cs  = ca - a*1.0e-6*p/ga
        gs  = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci  = cs - a*1.0e-6*p/gs
        faf = a - c4_jc(ci,p,rd,kt)
        !print*, 'C4 solv loop:',ca,cs,ci,ga,gs,a,faf,p,rd
      IF ((faf.LT.0.0d0).and.(a.lt.maxa)) THEN
        fa0 = faf
        a0  = a
        GOTO 50
      ENDIF
      af = a
      
      IF (af.GT.maxa) THEN
        jc = 1000.0d0
        !print*, 'jc > maxa'
        !cs = ca - a*1.0e-6*p/ga
        !gs = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        !ci = cs - a*1.0e-6*p/gs
      ELSE
        a  = (a0 + af)/2.0d0
        cs = ca - a*1.0e-6*p/ga
        gs = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci = cs - a*1.0e-6*p/gs
        fa = a - c4_jc(ci,p,rd,kt)

        bx = ((a+a0)*(faf-fa)/(af-a)-(af+a)*(fa-fa0)/(a-a0))/(a0-af)
        ax = (faf-fa-bx*(af-a))/(af**2-a**2)
        cx = fa-bx*a-ax*a**2

        IF(abs(ax).GT.0.0d0) THEN
          IF(bx**2-4.0d0*ax*cx.GT.0.0d0) THEN
            a = (-bx+(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
            IF(a.GT.af)  a = (-bx-(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
          ELSE
            a = 0.0d0
          ENDIF
        ELSE
          IF(abs(bx).GT.0.0d0) THEN
            a =-cx/bx
          ELSE
            a = 0.0d0
          ENDIF
        ENDIF
        jc = a
      ENDIF

      ! just take the minimum instead of co-limitation
      a = maxa 
      if(a.gt.jc) a = jc
      !print*, 'C4 processes:',jqq,jc,je

      ! calculate gs and ci from final a
      cs = ca - a*1.0e-6*p/ga
      gs = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
      ci = cs - a*1.0e-6*p/gs

      RETURN
      END
 
*----------------------------------------------------------------------*
      SUBROUTINE ASSC4_colim(a,ci,gs,ca,vmax,rh,kg,t,q,p,up,rd,dv,g0,g1,
     &gs_func,ga)
 
      IMPLICIT NONE

      REAL*8 gs_leaf,c4_colim 
      REAL*8 a,pc,gs,pa,vmax,rh,t,q,up,p,vmq,absx,qrub,f,jqq,dv
      REAL*8 amxt,ac4t,kg,rd,ca,cs,ci,je,k,kt,alpha
      REAL*8 fa,fa0,fa00,faf,a0,af,step
      REAL*8 ax,bx,cx,g0,g1,ga,maxa
      INTEGER gs_func

      ! light limited rate (Ji)
      absx  = 0.8d0
      qrub  = 0.11d0
      f     = 0.6d0
      alpha = 0.04d0
      jqq   = absx*qrub*f*q*1.0d6 - rd*1.0d6

      ! RuBisCO limited rate (Je)
      je = vmax*1.0e6 - rd*1.0d6

      ! Determine max possible assimilation rate 
      maxa = jqq
      if(maxa.gt.je) maxa = je
      !print*, 'C4 min(Je,Ji):',maxa

      ! solve colimited net assimiilation rate 
      k     = 0.70d0
      kt    = k * 2**( (t-25)/10 ) 

      step = 1.0d0
      a    = 0.0d0
      cs   = ca - a*1.0e-6*p/ga
      gs   = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
      ci   = cs - a*1.0e-6*p/gs
      faf  = a - c4_colim(ci,p,rd,alpha,q,vmax,kt)
      a0   = a
      fa0  = faf

50    CONTINUE
        !print*, 'C4 solv loop:',ca,cs,ci,ga,gs,a,faf,vmax,rd
        a   = a + step
        cs  = ca - a*1.0e-6*p/ga
        gs  = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci  = cs - a*1.0e-6*p/gs
        faf = a - c4_colim(ci,p,rd,alpha,q,vmax,kt)
      !IF ((faf.LT.0.0d0).and.(a.le.maxa+1)) THEN
      IF (faf.LT.0.0d0) THEN
        fa0 = faf
        a0  = a
        GOTO 50
      ENDIF
      af = a
      !print*, 'C4 solv done:',ca,cs,ci,ga,gs,a,faf,vmax,rd,t
      
      IF( (af.GT.maxa+1).and.(maxa.gt.0.0d0) ) THEN
        print*, 'C4 solver: a > maxa'
        print*, 'a:',af,'maxa:',maxa,'Ji:',jqq,'Je:',je
        print*, 'FORTRAN STOP'
        stop
        !cs = ca - a*1.0e-6*p/ga
        !gs = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        !ci = cs - a*1.0e-6*p/gs
      ELSE
        a  = (a0 + af)/2.0d0
        cs = ca - a*1.0e-6*p/ga
        gs = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
        ci = cs - a*1.0e-6*p/gs
        fa = a - c4_colim(ci,p,rd,alpha,q,vmax,kt)

        bx = ((a+a0)*(faf-fa)/(af-a)-(af+a)*(fa-fa0)/(a-a0))/(a0-af)
        ax = (faf-fa-bx*(af-a))/(af**2-a**2)
        cx = fa-bx*a-ax*a**2

        IF(abs(ax).GT.0.0d0) THEN
          IF(bx**2-4.0d0*ax*cx.GT.0.0d0) THEN
            a = (-bx+(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
            IF(a.GT.af)  a = (-bx-(bx**2-4.0d0*ax*cx)**0.5)/(2.0d0*ax)
          ELSE
            a = 0.0d0
          ENDIF
        ELSE
          IF(abs(bx).GT.0.0d0) THEN
            a =-cx/bx
          ELSE
            a = 0.0d0
          ENDIF
        ENDIF
      ENDIF

      ! just take the minimum instead of co-limitation
      !a = maxa 
      !if(a.gt.jc) a = jc
      !print*, 'C4 processes:',a,jqq,je
      !print*, ''

      ! calculate gs and ci from final a
      cs = ca - a*1.0e-6*p/ga
      gs = GS_LEAF(g0,g1,a,dv,t,cs,gs_func) * kg
      ci = cs - a*1.0e-6*p/gs

      RETURN
      END
 
*----------------------------------------------------------------------*
*     Low CO2 photosynthesis in C4 plants (Jc)
*     Collatz etal 1992
*----------------------------------------------------------------------*
 
      FUNCTION C4_Jc(ci,p,rd,kt)
      
      IMPLICIT NONE

      REAL*8 :: c4_jc
      REAL*8 :: ci,p,rd,kt

      c4_jc = (kt*ci/p - rd)*1.0d6

      END

*----------------------------------------------------------------------*
*     Colimited photosynthesis in C4 plants 
*     Collatz etal 1992
*----------------------------------------------------------------------*
 
      FUNCTION C4_colim(ci,p,rd,alpha,q,vmax,kt)
      
      IMPLICIT NONE

      REAL*8 :: c4_colim,quad_smallroot
      REAL*8 :: ci,p,rd,alpha,q,vmax,theta,beta,m,b,c,agross,kt
      
      theta = 0.83d0
      beta  = 0.93d0

      ! colimited Je and Ji photosynthesis
      b = vmax + alpha*q
      c = vmax*alpha*q  
      m = quad_smallroot(theta,-b,c)
      !print*, '     C4 colim M:',m,b,c
      
      ! colimited above and Jc photosynthesis
      b = m + kt*ci/p
      c = m*kt*ci/p
      agross = quad_smallroot(beta,-b,c)
      !print*, '     C4 colim A:',agross,b,c
 
      c4_colim = (agross - rd)*1.0d6

      END


*----------------------------------------------------------------------*
      FUNCTION quad_smallroot(a,b,c)

      IMPLICIT NONE

      REAL*8 :: quad_smallroot
      REAL*8 :: a,b,c,q,r1,r2

      if(b.lt.0.0d0) then
        q = 0.50d0*( -b + sqrt(b*b-4.0d0*a*c) )
      else
        q = 0.50d0*( -b - sqrt(b*b-4.0d0*a*c) )
      endif
      r1 = q/a
      r2 = c/q
      !print*,'            roots:',r1,r2

      quad_smallroot = r1
      if(r2.lt.r1)  quad_smallroot = r2

      END

*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

!####Brent solver from wikipedia
      FUNCTION BRENT_SOLVER(func,i1,i2,
     &oi,ca,vmx,j,rh,kg,rd,ga,t,p,gs_func,g0,g1,dv,i,ci,light,ttype,
     &ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ,tmonth) 
 
      IMPLICIT NONE
      
      REAL*8 :: brent_solver,i1,i2,errortol
      REAL*8 :: a,b,c,d,s,fa,fb,fc,fs,tmp,tmp2,fa00
      REAL*8 :: fassv,fassj,vcmax_maire
      REAL*8 :: ftToptV,ftHaV,ftHdV,ftToptJ,ftHaJ,ftHdJ,tmonth
      REAL*8 :: oi,ca,rh,rd,vmx,j,p,ga,dv,t,cs,kg,g0,g1,ci,light
      REAL*8 :: kc,ko,tau,km,gstar,vt,jt,alp,t_scalar
      INTEGER :: func,gs_func,i,n,farq_pars_func,ttype
      LOGICAL :: mflag,done

      ! Error tolerance umol m-2 s-1
      errortol = 1d-3
      done     = .false. 
      fs       = -1.d0 
      farq_pars_func = 1

      ! set parameter values
      CALL FARQ_PARS(t,kc,ko,tau,alp,farq_pars_func)
      km    = kc*(1 + oi/ko)
      gstar = 0.5d0*oi/tau 
     
      !initial guess  
      a  = i1
      b  = i2
c      if (func.eq.0) then
        ! Error tolerance mol m-2 s-1
        errortol = 1d-7
        ! vcmax25 & jmax25 to leaf t scalar 
        vt = T_SCALAR(tmonth,'v',ttype,ftToptV,ftHaV,ftHdV,tmonth)
        jt = T_SCALAR(tmonth,'j',ttype,ftToptJ,ftHaJ,ftHdJ,tmonth)
        fa = VCMAX_MAIRE(a,vt,jt,ci,km,gstar,alp,light,farq_pars_func,i)
        fb = VCMAX_MAIRE(b,vt,jt,ci,km,gstar,alp,light,farq_pars_func,i)
c      elseif (func.eq.1) then
c        fa = a -FASSV(a,gstar,km,ca,vmx,kg,rd,ga,t,p,gs_func,g0,g1,dv,i)
c        fb = b -FASSV(b,gstar,km,ca,vmx,kg,rd,ga,t,p,gs_func,g0,g1,dv,0)
c        if(fa.ge.0.d0) then
c          ! in this case rd is greater than gross a 
c          ! therefore assume a = 0 and anet = rd i
c          b = 0.d0
c          done = .true.
c        endif
c      elseif(func.eq.2) then
c        fa = a - FASSJ(a,gstar,ca,kg,rd,ga,t,p,j,gs_func,g0,g1,dv,i)
c        fb = b - FASSJ(b,gstar,ca,kg,rd,ga,t,p,j,gs_func,g0,g1,dv,0)
c        if(fa.ge.0.d0) then
c          ! in this case rd is greater than gross a 
c          ! therefore assume a = 0 and anet = rd i
c          b = 0.d0
c          done = .true.
c        endif
c      else 
c        print*, 'incorrect solver function number specified'
c      endif
      
      fa00 = fa 
      !if(i.eq.1) print*, i,fa,fb
      !if f(a) f(b) >= 0 then error-exit
      if((.not.done).and.(fa*fb.ge.0).and.(func.gt.0)) then
        print*, 'f(a) and f(b) do not span 0'
        print*, 'a =',a,'fa =',fa,'b =',b,'fb =',fb
        print*, 'solver function:', func 
        stop
      elseif((.not.done).and.(abs(fa).lt.abs(fb))) then
      !if |f(a)| < |f(b)| then swap (a,b) end if
        tmp = a
        a   = b
        b   = tmp
        tmp = fa
        fa  = fb
        fb  = tmp
      endif
 
      c     = a
      fc    = fa
      mflag = .TRUE.
 
      n = 0 
      do 
        if(done.or.(fb.eq.0.d0).or.(fs.eq.0.d0).or.
     &(abs(b-a).lt.errortol)) exit 
    
        if((fa.ne.fc).and.(fb.ne.fc)) then
          ! Inverse quadratic interpolation
          s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc /
     &(fb - fa) / (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb)
        else
          ! Secant Rule
          s = b - fb * (b - a) / (fb - fa)
        endif
    
        tmp2 = (3.d0 * a + b) / 4.d0
        if((.not.(((s.gt.tmp2).and.(s.lt.b)).or.
     &((s.lt.tmp2).and.(s.gt.b)))).or.
     &(mflag.and.(abs(s - b).ge.(abs(b - c)/2.d0))).or.
     &(.not.mflag.and.(abs(s - b).ge.(abs(c - d)/2.d0))).or.
     &(mflag.and.(abs(b - c).lt.errorTol)).or.
     &(.not.mflag.and.(abs(c - d).lt.errorTol))) then
          s = (a + b) / 2.d0
          mflag = .TRUE.
        else
          mflag = .FALSE.
        endif 
    
        !calculates a new value based on function
        !s is the new value 
        !fs = function(s)
    
c        if (func.eq.0) then
        fs = VCMAX_MAIRE(s,vt,jt,ci,km,gstar,alp,light,farq_pars_func,i)
c        elseif (func.eq.1) then
c        fs=s-FASSV(s,gstar,km,ca,vmx,kg,rd,ga,t,p,gs_func,g0,g1,dv,0)
c        elseif(func.eq.2) then
c        fs = s - FASSJ(s,gstar,ca,kg,rd,ga,t,p,j,gs_func,g0,g1,dv,0)
c        endif  
   
        d  = c
        c  = b
        fc = fb
        if((fa * fs).lt.0) then
          b  = s
          fb = fs
        else
          a  = s
          fa = fs
        endif 
    
        !# if |f(a)| < |f(b)| then swap (a,b)
        if(abs(fa).lt.abs(fb)) then
          tmp = a
          a   = b
          b   = tmp
          tmp = fa
          fa  = fb
          fb  = tmp
        endif   

        n = n + 1
      enddo
     
      !if((i.eq.1).and.(func.eq.0)) then
      !  print*, 't =',t,'ci =',ci,'light=',light,'vm =',b*1d6
      !  print*, 'a =',i1,'fa =',fa00,'final=',b,'n =',n,done
      !endif      

      brent_solver = b

      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   THE SUBROUTINES BELOW ARE NOT USED   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSVMAX                          *
*                          ******************                          *
*                                                                      *
* ASSVMAX computes the assimilation rate by forming a cubic equation   *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. The variables here *
* correspond to the Doly paper, in the main program they are such that *
*        po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd                           *
*     = oi,ca,vmx(i),c1,c2,rh,ko,kc,kg,tau,rd.                         *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSVMAX2(a,w,pc,gs,po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd
     &,acheck)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,w,pc,gs,po,pa,vmax,g0,g1,rh,ko,kc,kg,tau,rd,acheck,r2q3
      REAL*8 ax,bx,cx,dx,ex,fx,gx,hx,p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc
      REAL*8 arg1

      pi = atan(1.0d0)*4.0d0
      sc = 1000000.0d0
*----------------------------------------------------------------------*
* 'a2,b2,c2,d2' are such that 'w=(a2a + b2)/(c2a + d2)'.               *
*----------------------------------------------------------------------*
      ax = ko*pa*vmax*(kg*g1*rh - 160.0d0)
      bx = kg*g0*pa**2.0d0*vmax*ko
      cx =-160.0d0*pa*ko + kc*kg*ko*g1*rh + pa*kg*ko*g1*rh +
     & kc*kg*po*g1*rh
      dx = pa**2.0d0*kg*ko*g0 + kc*kg*po*g0*pa + kc*kg*ko*g0*pa

*----------------------------------------------------------------------*
* 'ex,fx,gx,hx' are such that '0.5po/(taupc)-1=(exa+fx)/(gxa+hx)'.     *
*----------------------------------------------------------------------*
      ex =-tau*pa*kg*g1*rh + 160.0d0*tau*pa + 0.5d0*po*kg*g1*rh
      fx =-tau*pa**2*kg*g0 + 0.5d0*po*kg*g0*pa
      gx = tau*pa*(kg*g1*rh - 160.0d0)
      hx = tau*pa**2.0d0*kg*g0

*----------------------------------------------------------------------*
* Cubic coefficients 'p,q,r'.                                          *
*----------------------------------------------------------------------*
      p = sc*rd + (sc*ax*ex + cx*hx + dx*gx)/(cx*gx)
      q = (sc*rd*cx*hx + sc*rd*dx*gx + sc*ax*fx + sc*bx*ex +
     & dx*hx)/(cx*gx)
      r = (sc*rd*dx*hx + sc*bx*fx)/(cx*gx)

*----------------------------------------------------------------------*
* Sove the cubic equaiton to find 'a'.                                 *
*----------------------------------------------------------------------*
      Qt = (p**2 - 3.0d0*q)/9.0d0
      Rt = (2.0d0*p**3 - 9.0d0*p*q + 27.0d0*r)/54.0d0

      r2q3 = Rt**2 - Qt**3

      IF (r2q3.LT.0.0d0) THEN
        arg1 = Rt/Qt**1.5d0
        IF (arg1.GT.1.0d0)  arg1 = 1.0d0
        IF (arg1.LT.-1.0d0)  arg1 = -1.0d0
        th = acos(arg1)
        a1 =-2.0d0*Qt**0.5d0*cos(th/3.0d0) - p/3.0d0
        a2 =-2.0d0*Qt**0.5d0*cos((th + 2.0d0*pi)/3.0d0) - p/3.0d0
        a3 =-2.0d0*Qt**0.5d0*cos((th + 4.0d0*pi)/3.0d0) - p/3.0d0
        a = a1
        IF (a2.GT.a)  a = a2
        IF (a3.GT.a)  a = a3
      ELSE
        at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5d0)**(1.0d0/3.0d0)
        IF (abs(at).GT.1E-6) THEN
          bt = Qt/at
        ELSE
          bt = 0.0d0
        ENDIF
        a1 = at + bt - p/3.0d0
*        a2 =-0.5d0*(at + bt) - p/3.0d0 + i*3.0d0**0.5d0*(at - bt)/2.0d0
*        a3 =-0.5d0*(at + bt) - p/3.0d0 - i*3.0d0**0.5d0*(at - bt)/2.0d0
        a = a1
      ENDIF

*----------------------------------------------------------------------*
* Compute 'gs,pc,w' corresponding to 'a'.                              *
*----------------------------------------------------------------------*
      gs = (g0 + g1*a*rh/pa)*kg
      pc = pa - a*160.0d0/gs
      w = vmax*pc/(pc + kc*(1 + po/ko))
      acheck = (w*(1.0d0 - 0.5*po/(tau*pc)) - rd)*sc


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSJ                             *
*                          ***************                             *
*                                                                      *
* ASSJ computes the assimilation rate by forming a cubic equation      *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. The variables here *
* correspond to the Doly paper, in the main program they are such that *
*        po,pa,j,g0,g1,rh,kg,tau,rd                                    *
*     = oi,ca,j(i),c1,c2,rh,kg,tau,rd                                  *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSJ2(a,w,pc,gs,po,pa,j,g0,g1,rh,kg,tau,rd,acheck)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,w,pc,gs,po,pa,j,g0,g1,rh,kg,tau,rd,acheck,r2q3,arg1
      REAL*8 ax,bx,cx,dx,ex,fx,gx,hx,p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc

      pi = atan(1.0d0)*4.0d0
      sc = 1000000.0d0

*----------------------------------------------------------------------*
* 'a2,b2,c2,d2' are such that 'w=(a2a + b2)/(c2a + d2)'.               *
*----------------------------------------------------------------------*
      ax = j*tau*pa*(kg*g1*rh - 160.0d0)
      bx = j*tau*pa*kg*g0*pa
      cx = 4.0d0*pa*kg*tau*g1*rh - 640.0d0*pa*tau + 4.0d0*po*kg*g1*rh
      dx = 4.0d0*pa**2.0d0*kg*tau*g0 + 4.0d0*po*kg*g0*pa

*----------------------------------------------------------------------*
* 'ex,fx,gx,hx' are such that '0.5po/(taupc)-1=(exa+fx)/(gxa+hx)'.     *
*----------------------------------------------------------------------*
      ex =-tau*pa*kg*g1*rh + 160.0d0*tau*pa + 0.5d0*po*kg*g1*rh
      fx =-tau*pa**2.0d0*kg*g0 + 0.5d0*po*kg*g0*pa
      gx = tau*pa*(kg*g1*rh - 160.0d0)
      hx = tau*pa**2.0d0*kg*g0

*----------------------------------------------------------------------*
* Cubic coefficients 'p,q,r'.                                          *
*----------------------------------------------------------------------*
      p = sc*rd + (sc*ax*ex + cx*hx + dx*gx)/(cx*gx)
      q = (sc*rd*cx*hx + sc*rd*dx*gx + sc*ax*fx + sc*bx*ex +
     & dx*hx)/(cx*gx)
      r = (sc*rd*dx*hx + sc*bx*fx)/(cx*gx)

*----------------------------------------------------------------------*
* Sove the cubic equaiton to find 'a'.                                 *
*----------------------------------------------------------------------*
      Qt = (p**2 - 3.0d0*q)/9.0d0
      Rt = (2.0d0*p**3 - 9.0d0*p*q + 27.0d0*r)/54.0d0

      r2q3 = Rt**2 - Qt**3

      IF (r2q3.LT.0.0d0) THEN
        arg1 = Rt/Qt**1.5d0
        IF (arg1.GT.1.0d0)  arg1 = 1.0d0
        IF (arg1.LT.-1.0d0)  arg1 = -1.0d0
        th = acos(arg1)
        a1 =-2.0d0*Qt**0.5d0*cos(th/3.0d0) - p/3.0d0
        a2 =-2.0d0*Qt**0.5d0*cos((th + 2.0d0*pi)/3.0d0) - p/3.0d0
        a3 =-2.0d0*Qt**0.5d0*cos((th + 4.0d0*pi)/3.0d0) - p/3.0d0
        a = a1
        IF (a2.GT.a)  a = a2
        IF (a3.GT.a)  a = a3
      ELSE
        at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5d0)**(1.0d0/3.0d0)
        IF (abs(at).GT.1E-6) THEN
          bt = Qt/at
        ELSE
          bt = 0.0d0
        ENDIF
        a1 = at + bt - p/3.0d0
*        a2 =-0.5d0*(at + bt) - p/3.0d0 + i*3.0d0**0.5d0*(at - bt)/2.0d0
*        a3 =-0.5d0*(at + bt) - p/3.0d0 - i*3.0d0**0.5d0*(at - bt)/2.0d0
        a = a1
      ENDIF

*----------------------------------------------------------------------*
* Compute 'gs,pc,w' corresponding to 'a'.                              *
*----------------------------------------------------------------------*
      gs = (g0 + g1*a*rh/pa)*kg
      pc = pa - a*160.0d0/gs
      w = j*pc/4.0d0/(pc + po/tau)
      acheck = (w*(1.0d0 - 0.5d0*po/(tau*pc)) - rd)*sc


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE XAMAX                            *
*                          ****************                            *
*                                                                      *
* This subroutine calculates the value of 'amax' at the maximum        *
* iradience as opposed to the average irradiance in the main program.  *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE XAMAX(a,xt,oi,vmx,rd,j,dresp,pc)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,xt,oi,vmx,rd,j,dresp,pc
      REAL*8 t,tk,c1,c2,kc,ko,tau,tpav,tpwv,tpaj,tpwj,sc
*      REAL*8 tppcj,tppcv,tpgsv,tpgsj,acheck

      t = 2.64d0 + 1.14d0*xt
      tk = 273.0d0 + t
      sc = 1000000.0d0

      c1 = 142.4d0 - 4.8d0*t
      IF (c1.GT.80.0d0)  c1 = 80.0d0
      IF (c1.LT.8.0d0)  c1 = 8.0d0
      c2 = 12.7d0 - 0.207d0*t
      IF (c2.GT.10.0d0)  c2 = 10.0d0
      IF (c2.LT.6.9d0)  c2 = 6.9d0
      kc = exp(35.8d0 - 80.5d0/(0.00831d0*tk))
      ko = exp(9.6d0 - 14.51d0/(0.00831d0*tk))*1000.0d0
      tau = exp(-3.949d0 + 28.99d0/(0.00831d0*tk))

*      rht = rh
*      IF (kg*c2*rh.LT.170.0d0)  rht = 170.0d0/(kg*c2)
*      CALL ASSVMAX(tpav,tpwv,tppcv,tpgsv,oi,ca,vmx,c1,c2,
*     &rht,ko,kc,kg,tau,rd,acheck)
*      CALL ASSJ(tpaj,tpwj,tppcj,tpgsj,oi,ca,j,c1,c2,rht,
*     &kg,tau,rd,acheck)

      tpwv = vmx*pc/(pc + kc*(1.0d0 + oi/ko))
      tpav = (tpwv*(1.0d0 - 0.5d0*oi/(tau*pc)) - rd)*sc

      tpwj = j*pc/(4.0d0*(pc + oi/tau))
      tpaj = (tpwj*(1.0d0 - 0.5d0*oi/(tau*pc)) - rd)*sc

      IF (tpav.LT.tpaj) THEN
        A = tpav
      ELSE
        A = tpaj
      ENDIF
      IF (A.LT.-dresp)  A =-dresp


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ASSC4                            *
*                          ****************                            *
*                                                                      *
* ASS4 computes the assimilation rate by forming a cubic equation      *
* from the equations defining assimilation rate, stomatal conductance, *
* CO2 partial pressure (assimilation diffusion equaiton) and the       *
* carboxylation rate. The cubic equation is solved and the largest     *
* real root is assumed to be the assimilation rate. This is then       *
* compared with assimilation limited to light dependant and vmax       *
* dependant assimilation and the minimum one is taken. The variables   *
* here correspond to the Doly paper, in the main program in the        *
* folowing mannor.                                                     *
*        po,pa,j,g0,g1,rh,kg,tau,rd                                    *
*     = oi,ca,j(i),c1,c2,rh,kg,tau,rd                                  *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE ASSC42(a,w,pc,gs,po,pa,vmax,g0,g1,rh,kg,tau,rd,t,qg,
     &acheck,up)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 a,w,pc,gs,po,pa,vmax,g0,g1,rh,kg,tau,rd,t,qg,acheck,up
      REAL*8 p,q,r,qt,rt,th,a1,a2,a3,pi,at,bt,sc,kp,kmq,pg,l,vmq,absx
      REAL*8 qrub,f,a0,b0,b1,c1,d1,e1,f1,b2,c2,d2,b3,c3,d3,wq,wv
      REAL*8 a4,b4,c4,d4,f2,aq,av,a2t,r2q3,arg1

      pi = atan(1.0d0)*4.0d0
      sc = 1000000.0d0

      kp = 0.7d0
      kmq = 2.0d0
      pg = 101325.0d0
      l = 0.000005d0
      vmq = 2.0d0
      absx = 0.8d0
      qrub = 0.11d0
      f = 0.6d0

      vmax = up*360.0d0*vmq**((t - 25.0d0)/10.0d0)/(360.0d0 + up)
      vmax = vmax/((1 + exp(0.3d0*(10.0d0 - t)))*(0.8d0 +
     &exp(0.14d0*(t - 36.0d0))))
      vmax = vmax*1.0e-6*1.25d0

*----------------------------------------------------------------------*
* 'a' and 'b' are such that 'jc=a pc  - b'.                            *
*----------------------------------------------------------------------*
      a0 = kp*kmq**((t - 25.0d0)/10.0d0)/pg
      b0 = l*kmq**((t - 25.0d0)/10.0d0)/pg

*----------------------------------------------------------------------*
* Reduce equations to a cubic in 'a'.                                  *
*----------------------------------------------------------------------*
      a1 =-b0*sc
      b1 = a0*sc
      c1 =-0.5d0*po
      d1 = tau
      e1 =-tau*rd*sc
      f1 = tau

      a2 = g0*kg*(pa**2.0d0)
      b2 = g0*kg*pa
      c2 =-pa*160.0d0 + g1*rh*kg*pa
      d2 = g1*rh*kg

      a3 = a1*b2 + a2*b1
      b3 = a1*d2 + c2*b1
      c3 = c1*b2 + a2*d1
      d3 = c1*d2 + c2*d1

      f2 =-f1

      a4 = f2*c2*d2
      b4 = f2*(a2*d2 + b2*c2) + d3*b3 + e1*c2*d2
      c4 = f2*a2*b2 + c3*b3 + d3*a3 + e1*(a2*d2 + c2*b2)
      d4 = c3*a3 + e1*a2*b2

*----------------------------------------------------------------------*
* Cubic coefficients 'p,q,r'.                                          *
*----------------------------------------------------------------------*
      p = b4/a4
      q = c4/a4
      r = d4/a4

      a2t = a2
*----------------------------------------------------------------------*
* Sove the cubic equaiton to find 'a'.                                 *
*----------------------------------------------------------------------*
      Qt = (p**2 - 3.0d0*q)/9.0d0
      Rt = (2.0d0*p**3 - 9.0d0*p*q + 27.0d0*r)/54.0d0

      r2q3 = Rt**2 - Qt**3

      IF (r2q3.LT.0.0d0) THEN
        arg1 = Rt/Qt**1.5d0
        IF (arg1.GT.1.0d0)  arg1 = 1.0d0
        IF (arg1.LT.-1.0d0)  arg1 = -1.0d0
        th = acos(arg1)
        a1 =-2.0d0*Qt**0.5d0*cos(th/3.0d0) - p/3.0d0
        a2 =-2.0d0*Qt**0.5d0*cos((th + 2.0d0*pi)/3.0d0) - p/3.0d0
        a3 =-2.0d0*Qt**0.5d0*cos((th + 4.0d0*pi)/3.0d0) - p/3.0d0
        a = a1
        IF (a2.GT.a)  a = a2
        IF (a3.GT.a)  a = a3
      ELSE
        at =-Rt/abs(Rt)*(abs(Rt) + r2q3**0.5d0)**(1.0d0/3.0d0)
        IF (abs(at).GT.1E-6) THEN
          bt = Qt/at
        ELSE
          bt = 0.0d0
        ENDIF
        a1 = at + bt - p/3.0d0
*        a2 =-0.5d0*(at + bt) - p/3.0d0 + i*3.0d0**0.5d0*(at - bt)/2.0d0
*        a3 =-0.5d0*(at + bt) - p/3.0d0 - i*3.0d0**0.5d0*(at - bt)/2.0d0
        a = a1
      ENDIF

      gs = (g0 + g1*a*rh/pa)*kg
      pc = pa - a*160.0d0/gs
      w = a0*pc - b0
      acheck = (w*(1.0d0 - 0.5d0*po/(tau*pc)) - rd)*sc

      a2 = a2t

*----------------------------------------------------------------------*
* Compute light dependent assimilation rate.                           *
*----------------------------------------------------------------------*

      wq = absx*qrub*f*qg*1.0d0
      p = tau*c2
      q = tau*a2 + wq*0.5d0*po*sc*d2 - tau*c2*(wq*sc - rd*sc)
      r = wq*0.5d0*po*sc*b2 - tau*a2*(wq*sc - rd*sc)
      aq = (-q + (abs(q)**2.0d0 - 4.0d0*p*r)**0.5d0)/(2.0d0*p)

*----------------------------------------------------------------------*
* Compute 'vmax' dependent assimilation rate.                          *
*----------------------------------------------------------------------*

      wv = vmax*vmq**((t - 25.0d0)/10.0d0)/((1.0d0 + exp(0.3d0*
     &(10.0d0 - t)))*(0.8d0 + exp(0.14d0*(t - 36.0d0))))

      q = tau*a2 + wv*0.5*po*sc*d2 - tau*c2*(wv*sc - rd*sc)
      r = wv*0.5*po*sc*b2 - tau*a2*(wv*sc - rd*sc)
      av = (-q + (abs(q)**2.0d0 - 4.0d0*p*r)**0.5d0)/(2.0d0*p)

*----------------------------------------------------------------------*
* Find limiting assimilation value.                                    *
*----------------------------------------------------------------------*
      IF (aq.LT.a) THEN
        a = aq
        w = wq
        gs = (g0 + g1*a*rh/pa)*kg
        pc = pa - a*160.0d0/gs
      ENDIF

      IF (av.LT.a) THEN
        a = av
        w = wv
        gs = (g0 + g1*a*rh/pa)*kg
        pc = pa - a*160.0d0/gs
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE XVJMAX                           *
*                          *****************                           *
*                                                                      *
* This subroutine calculates the value of 'amax' at the maximum        *
* iradience as opposed to the average irradiance in the main program.  *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE XVJMAX(vmx,jmx,j,xt,xq,sum,nup,oi,dresp)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 vmx,jmx,j,xt,xq,sum,nup,oi
      REAL*8 tk,shapex,maxx,conv,cc,ea,num,denom,up,kc,ko,tau,can,upt,am
      REAL*8 vm,jm,irc,q,t,aa,bb,dresp

      t = 2.64d0 + 1.14d0*xt
      tk = 273.0d0 + t
      q = 2.0d0*xq

      nup = 1.5d0*nup

      shapex = 40.8d0 + 0.01d0*t - t**2*0.002d0
      maxx = 0.738d0 - 0.002d0*t
      IF (nup.LE.0.0d0) THEN
         print *,'Problem in XVJMAX nup=',nup
         stop
      ENDIF
      conv = 97.4115d0 - 2.5039d0*log(nup)

      cc = 36.0d0*exp(-t*0.14d0) + 20.0d0
      ea = 81000.0d0*exp(-t*0.14d0) + 50000.0d0
      num = exp(shapex - conv/(0.00831d0*tk))
      denom = 1.0d0 + exp((maxx*tk - 205.9d0)/(0.00831d0*tk))
      up = num/denom

      kc = exp(35.8d0 - 80.5d0/(0.00831d0*tk))
      ko = exp(9.6d0 - 14.51d0/(0.00831d0*tk))*1000.0d0
      tau = exp(-3.949d0 + 28.99d0/(0.00831d0*tk))
      aa = 24.5d0/(24.5d0 + kc + kc*oi/ko)
      bb = 0.5d0/(24.5d0 + kc + kc*oi/ko)


      can = 1.0d0/sum
      upt = up*can
      am = upt*190.0d0/(360.0d0 + upt)

      dresp = exp(cc-(ea/(8.3144d0*tk)))*upt/50.0d0

      vm =-1.0d0*(am + 0.82d0)/(-1.0d0*aa + bb*oi/tau)
      jm = 29.1d0 + 1.64d0*vm
      vm = vm/1000000.0d0
      jm = jm/1000000.0d0

      IF (t.GT.6.34365d0) THEN
        jmx = jm*(1.0d0 + 0.0409d0*(t-25.0d0)
     &- 1.54E-3*(t - 25.0d0)**2 - 9.42E-5*(t - 25.0d0)**3)
      ELSE
        jmx = jm*0.312633d0
      ENDIF

      IF (t.GT.9.51718d0) THEN
        vmx = vm*(1.0d0 + 0.0505d0*(t - 25.0d0)
     &- 0.248E-3*(t - 25.0d0)**2 - 8.09E-5*(t-25.0d0)**3)
      ELSE
        vmx = vm*0.458928d0
      ENDIF

* Mark we have to discuss about this calculation for irc here...
      irc = q*exp(-0.5d0)
      j = 0.24d0*irc/(1.0d0 + (0.24d0**2)*(irc**2)/(jmx**2))
     &**0.5d0


      RETURN
      END

