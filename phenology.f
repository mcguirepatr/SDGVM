*----------------------------------------------------------------------*
*                          PHENOLOGY                                   *
*                          *********                                   *
*   Unified SDGVM tree and grass phenology and allocation subroutine   * 
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE PHENOLOGY(ftphen,ftagh,ftdth,
     &bbm,bb0,bbmax,bblim,ssm,sss,sslim,
     &nppstore,leafnpp,stemnpp,rootnpp,tran,rlai,lai,lairat,rem,leafls,
     &stemls,rootls,leafmol,respref,soil2g,wtwp,tmem,leafv,stemv,rootv,
     &daysoff,laimax,leaflit,stemlit,rootlit,mnth,day,s_ln,s_sr,s_sn,
     &s_rr,s_rn,bb,ss,bbgs,dsbb,nppstorx,nppstor2,daynpp,
     &maxlai,wtfc,yield,resp,sm_trig,suma,tsumam,stemfr,lmor_sc,
     &chill,dschill,s_lr,leafresp,rootresp,stemresp,
     &phen_cor,s070607,kg,peak_lai)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 nppstore,rootnpp,tran,rlai,lairat,leafmol,respref,soil2g
      REAL*8 wtwp,rem,leaflit,leafv(3600),stemv(1000),rootv(1000)
      REAL*8 ftagh,leafnpp,stemnpp,laiinc,tmem(200),resp,nppstor2,s_ln
      REAL*8 daysoff,laimax,stemnppr,rootnppr,nppstorx,wtfc,s_lr,s_mr
      REAL*8 s_sn,s_rr,s_rn,stemlit,rootlit,daynpp,maxlai,bb0,bbmax,s_sr
      REAL*8 bblim,sslim,bbsum,yield,sm_trig(30),smtrig,xdaynpp,ans
      REAL*8 suma(360),tsuma,tsumam,stemfr,maint,yy,lmor_sc(3600)
      REAL*8 leafresp,rootresp,stemresp,rtemp
      REAL*8 swc_bbthresh, swc_senthresh, minnppstore,lresp,kg,lgrowth
      REAL*8 peak_lai 
      REAL*8 sumrr,sumsr,sumlr,summr,resp_r,resp_s,resp_m,resp_l
      ! not sure what this save line does
      SAVE sumrr,sumsr,sumlr,summr
      INTEGER lai,leafls,stemls,rootls,mnth,day,ij,i,bb,gs,bbgs,sssum
      INTEGER bb2bbmin,bb2bbmax,bbm,ssm,sss,ss,dsbb
      INTEGER ftphen,ftdth,harvest,chill,dschill,phen_cor,s070607
      INCLUDE 'param.inc'


      ! set grass/tree parameter differences
      IF (ftphen.EQ.1) THEN
        ! minimum number of days between bud burst
        bb2bbmin = 285
        ! maximum number of days between bud burst
        bb2bbmax = 390
        ! length of leaf growing season
        gs       = 60
        ! fraction of max available soil water (fc - wp) needed for budburst
        swc_bbthresh  = 0.25d0
        ! fraction of max available soil water (fc - wp) that triggers senescence 
        swc_senthresh = 0.0d0
        ! minimum nppstore required for leaf growth
        minnppstore   = 0.0d0

      ELSEIF (ftphen.EQ.2) THEN
        bb2bbmin = 315
        bb2bbmax = 375
        gs       = 30
        swc_bbthresh  = 0.5d0
        swc_senthresh = 0.5d0
        minnppstore   = 5.0d0
      ENDIF

      ! zero variables 
      leaflit = 0.0d0
      stemlit = 0.0d0
      rootlit = 0.0d0
      yield   = 0.0d0
      resp    = 0.0d0

      ! crop harvest trigger
      harvest = 0
      ! day of senescence trigger
      ss      = 0

      !rem = rlai - int(rlai)
      !lai = int(rlai) + 1

      xdaynpp = daynpp

*----------------------------------------------------------------------*
* Compute C in live root and stem biomass and respiration stores
*----------------------------------------------------------------------*
      stemnppr = 0.0d0
      rootnppr = 0.0d0
      DO ij=1,stemls
        stemnppr = stemnppr + stemv(ij)
      ENDDO
      DO ij=1,rootls
        rootnppr = rootnppr + rootv(ij)
      ENDDO

*----------------------------------------------------------------------*
* Draw C from store for root growth 
*----------------------------------------------------------------------*
      IF ((nppstore.GT.0.0d0).AND.(daynpp.GT.0.0d0)) THEN
        yy = nppstore*p_rootfr
      ELSE
        yy = 0.0d0
      ENDIF

      ! update tracking variables
      nppstore = nppstore - yy
      xdaynpp  = xdaynpp  - yy
      rootnpp  = rootnpp  + yy
      s_rn     = s_rn     + yy

*----------------------------------------------------------------------*
* Age and grow roots in live root biomass and respiration store daily vector
*----------------------------------------------------------------------*
      DO ij=1,rootls-1
        rootv(rootls+1-ij) = rootv(rootls-ij)
      ENDDO
      rootv(1) = yy
      IF (rootv(1).LT.0.0d0) THEN
        print*,'rootv negativ: rootv=',rootv(1),'tran=',tran
      ENDIF

*----------------------------------------------------------------------*
* Root respiration
*----------------------------------------------------------------------*
      resp_r = resp
      DO ij=1,rootls-1
        ! calculate respiration for live roots 
        rtemp     = respref*rootv(ij)
        ! remove respiration from live biomass and respiration store 
        rootv(ij) = rootv(ij) - rtemp
        ! update tracking variables
        resp      = resp + rtemp
        rootresp  = rootresp + rtemp
        s_rr      = s_rr + rtemp
        rootnpp   = rootnpp - rtemp
      ENDDO
      resp_r = resp - resp_r

*----------------------------------------------------------------------*
* Calculate chilling requirement for budburst.                         *
* chill is 0/1 switch indicating 1(0) - chilling requirement (not) met *
* chilling requirement is -100 degree days, between -15 and -5, over 20 days 
* dschill is days since chilling requirement was met
*----------------------------------------------------------------------*
      IF (chill.EQ.0) THEN
        bbsum = 0.0d0
        DO i=1,20
          IF (tmem(i).LT.-5.0d0)  bbsum = bbsum + max(-10.0d0,tmem(i)+
     &5.0d0)
        ENDDO
        IF (bbsum.LT.-100) THEN
          chill   = 1
          dschill = 1
        ENDIF
      ENDIF

      IF (chill.EQ.1) THEN
        dschill = dschill + 1
      ENDIF

      IF (dschill.GT.260) THEN
        chill   = 0
        dschill = 0
      ENDIF

*----------------------------------------------------------------------*
* If soil water is sufficient, calculate growing degree days           *  
* if days since budburst (dsbb) is greater than maximum (bb2bbmax), relax soil water requirement *
*----------------------------------------------------------------------*
      IF (((bb.eq.0).AND.(soil2g.GT.wtwp+swc_bbthresh*(wtfc-wtwp))).OR.
     &((dsbb.GT.bb2bbmax).AND.(soil2g.GT.wtwp+0.1d0*(wtfc-wtwp)))) THEN

        ! calcualte cumulative soil water 
        smtrig = 0.0d0
        DO i=1,30
          smtrig = smtrig + sm_trig(i)
        ENDDO

        IF ((smtrig.GT.30.0d0).OR.(dsbb.GT.bb2bbmax)) THEN
          ! Calculate degree days for temp > bb0 and over bbm days
          bbsum = 0.0d0
          DO i=1,bbm
            IF (tmem(i).GT.bb0)  bbsum = bbsum + min(bbmax,tmem(i)-bb0)
          ENDDO

*----------------------------------------------------------------------*
* Bud burst occurance.                                                 *
* this routine occurs only on the day of budburst                      *
*----------------------------------------------------------------------*
          IF ((real(bbsum).GE.real(bblim)*exp(-0.01d0*real(dschill)))
     &.OR.(dsbb.GT.bb2bbmax)) THEN

            !bb is the day of the year of budburst
            bb = (mnth-1)*30 + day
            bbgs = 0
            dsbb = 0

*----------------------------------------------------------------------*
* LAI control.                                                         *
*----------------------------------------------------------------------*
*   - for grasses proportion of store going into leaf production *
*----------------------------------------------------------------------*
            !IF (ftphen.EQ.1) THEN
              
              IF(s070607.eq.1) THEN
!                nppstorx = nppstore
!              ELSEIF(phen_cor.eq.1) THEN 
!                ! restricts the maximum amount of the npp store to be used for leaf growth to 62.5% (i.e. 50% ends up as leaf mass)
!                nppstorx = 0.5 * 1.25 * nppstore
!                !nppstorx = 0.1 * 1.25 * nppstore
!              ELSE
!                nppstorx = nppstore
!              ENDIF
              
*----------------------------------------------------------------------*
*   - for trees adjust proportion of store available for stem production based on suma. *
*----------------------------------------------------------------------*
            ELSE !IF (ftphen.EQ.2) THEN
  
              !tsuma is the annual C balance of the lowest LAI layer.
              tsuma = 0.0d0
              DO i=1,360
                tsuma = tsuma + suma(i)
              ENDDO
              !leafls is leaf life span in days 
              !p_opt is the 'canopy optimisation correction' currently 1.5, 
              ! - can be thought of as the multipier on leaf C costs to account for roots and stem needed to support those leaves
              !maint is the age-based mean resisdence time of leaves 
              maint = max(1.0d0,(real(leafls)/360.0d0))
              !reduce bottom layer C balance by the annual cost of leaves multiplied by p_opt
              p_opt = 1.0d0
              tsuma = tsuma - leafmol*1.25d0/maint*p_opt
              !tsuma = tsuma -leafmol*1.25d0/maint -leafmol*1.25d0*p_opt

              !restrict the bottom layer C balance to be within +/- the annual cost of leaves * p_opt   
              IF (tsuma.GT.leafmol*1.25d0/maint*p_opt)
     &tsuma = leafmol*1.25d0/maint*p_opt
              IF (tsuma.LT.-leafmol*1.25d0/maint*p_opt)
     &tsuma = -leafmol*1.25d0/maint*p_opt
              
              !stemfr is the target amount of C to allocate to leaves (inc. growth respiration) for the year
              !p_laimem is the 'LAI memory' currently 0.5, effectively smoothing the rate of change of LAI from one year to the next 
              stemfr = stemfr + tsuma*p_laimem
  
              IF (stemfr.LT.10.0d0) stemfr = 10.0d0
  
              !if target leaf allocation is more than 75% of the nppstore allow only 75% of the store for leaf allocation 
              !and reduce the value of lai which is passed to next year by 5%  
              !nppstor2 is the nppstore - leaf C allocation 
              IF (stemfr.LT.0.75d0*nppstore) THEN
                nppstor2 = nppstore - stemfr
              ELSE
                nppstor2 = nppstore*0.25d0
                stemfr = stemfr*0.95
              ENDIF
  
              !nppstorx = nppstore
              nppstorx = nppstore - nppstor2
              peak_lai = nppstorx/leafmol/1.25d0 

            ENDIF
          ENDIF
        ENDIF
      ENDIF

*----------------------------------------------------------------------*
* Compute length of current leaf growing season (bbgs), 
* and stop leaf growth when current growing season (bbgs) = max growing season (gs)                                                      *
*----------------------------------------------------------------------*
      IF (bb.GT.0)  bbgs = bbgs + 1
      IF (bbgs-gs.gt.bb2bbmin) THEN
        bb   = 0
        bbgs = 0
      ENDIF

      ! calculate days since budburst
      IF (dsbb.LT.500) dsbb = dsbb + 1

*----------------------------------------------------------------------*
* Set LAI increase (laiinc) 
*----------------------------------------------------------------------*
      IF ((bb.GT.0).AND.(bbgs.LT.gs).AND.(nppstore.GT.1.0d0)) THEN

        lresp = 1.25d0
        IF((s070607.eq.1).AND.(ftphen.eq.1)) lresp = 1.0d0

        !daily LAI increment
        IF (ftphen.EQ.1) THEN
          IF(s070607.eq.1) THEN
            laiinc = lairat*nppstorx/leafmol/lresp
          ELSE 
            !laiinc = lairat*nppstorx/leafmol/lresp
            laiinc = lairat*kg*nppstorx/leafmol/lresp
          ENDIF
        ELSE IF (ftphen.EQ.2) THEN
          laiinc = lairat*nppstorx/leafmol/lresp
        ENDIF

        IF ((rlai.GT.0).and.(nppstore.LT.minnppstore)) laiinc = 0.0d0
        IF (rlai+laiinc.GT.maxlai)  laiinc = maxlai - rlai
        IF (rlai+laiinc.GT.11.5d0)  laiinc = 11.5d0 - rlai

      ELSE

        laiinc = 0.0d0

      ENDIF

      !print*,mnth,day,bb,bbgs,laiinc,gs,bb2bbmin

*----------------------------------------------------------------------*
* Senescence                                                           *
*----------------------------------------------------------------------*
      IF (rlai.GT.1.0e-6) THEN
        IF (abs(ftagh).LT.1.0E-6) THEN

*----------------------------------------------------------------------*
*   - for trees or grasses 
*----------------------------------------------------------------------*
* Soil water based  
*   - drop all leaves               *
*----------------------------------------------------------------------*
          IF (soil2g.LT.wtwp*0.0d0) THEN
            laiinc = -rlai
            ss = day + (mnth - 1)*30
  
*----------------------------------------------------------------------*
* Temperature based, senescence occurs there are 'sss' days colder  *
* than 'sslim' out of the last 'ssm' days.                             *
*   - drop all leaves               *
*----------------------------------------------------------------------*
          ELSEIF (bbgs.GT.100) THEN
            sssum = 0
            DO i=1,ssm
              IF (tmem(i).LT.sslim)  sssum = sssum + 1
            ENDDO
            IF (sssum.GE.sss) THEN
              laiinc =-rlai
              ss = day + (mnth - 1)*30
*          print*,'senescence temp'
*          print'(100f6.1)',(tmem(i),i=1,ssm)
            ENDIF
          ENDIF

*----------------------------------------------------------------------*
*   - for crops, age or harvest based 
*----------------------------------------------------------------------*
        ELSE
          IF ((bbgs.GT.leafls).OR.(bbgs.GT.ftdth)) THEN
            harvest = 1
            laiinc =-rlai
            ss = day + (mnth - 1)*30
          ENDIF

        ENDIF
      ENDIF

*----------------------------------------------------------------------*
* Pay for new leaves.
*----------------------------------------------------------------------*
      if(ftphen.eq.2) resp_m = resp - resp_m
      IF (laiinc.GT.0.0d0) THEN
        IF (ftphen.EQ.1) THEN
          ! leaf growth
          lgrowth  = laiinc*leafmol
          ! leaf growth respiration
          leafresp = 0.25d0*laiinc*leafmol
        ELSE IF (ftphen.EQ.2) THEN
          lgrowth  = laiinc*leafmol
          leafresp = 0.25d0*laiinc*leafmol
        ENDIF
        ! update store
        nppstore = nppstore - lgrowth - leafresp
        ! update tracking variables
        nppstorx = nppstorx - lgrowth - leafresp 
        xdaynpp  = xdaynpp  - lgrowth - leafresp 
        leafnpp  = leafnpp  + lgrowth 
        resp     = resp     + leafresp 
        s_mr     = leafresp
        s_ln     = s_ln + lgrowth + lresp 
      ENDIF
*      print*,'maint & lai inc ',12.0*0.25d0*laiinc*leafmol,laiinc
      if(ftphen.eq.2) resp_m = resp - resp_m
      !s_mr = 0.25d0*laiinc*leafmol


*----------------------------------------------------------------------*
* Age leaves by one day, kill any which have died of old age,
* adjust by laiinc (+ or -), and then sum to get 'rlai'.
* 'leafls' is an integer variable of leaf lifespan in days.
*----------------------------------------------------------------------*
      CALL LAIALT(leafv,rlai,laiinc,leafls,leaflit)

      ! update lai variables
      rem = rlai - int(rlai)
      lai = int(rlai) + 1

      ! crop harvest
      IF (harvest.EQ.1) THEN
        yield = leaflit*ftagh
        leaflit = (1.0d0 - ftagh)*leaflit
      ELSE
        yield = 0.0d0
      ENDIF

*----------------------------------------------------------------------*
* leaf death not through age mortality 
*----------------------------------------------------------------------*
      rlai = rlai + leaflit
      DO i=1,leafls
        IF (leafv(i).GT.0.0d0) THEN
          IF (leafv(i).GT.1.0e-6) THEN
            !leaflit  = leaflit + leafv(i)*(1 -lmor_sc(i))
            !leafv(i) = leafv(i)*lmor_sc(i)
            leaflit = leaflit + leafv(i)
            leafv(i)=leafv(i)*lmor_sc(i)
            leaflit = leaflit - leafv(i)
          ELSE
            leaflit  = leaflit + leafv(i)
            leafv(i) = 0.0d0
          ENDIF
        ENDIF
      ENDDO
      rlai = rlai - leaflit

*----------------------------------------------------------------------*
* Stem respiration, aging, and NPP                                     *
*----------------------------------------------------------------------*

      resp_s = resp - resp_s
      s_sr   = resp
      DO ij=1,stemls-1
        rtemp     = respref*stemv(ij)
        stemv(ij) = stemv(ij) - rtemp
        resp      = resp + rtemp
        stemresp  = stemresp + rtemp
        stemnpp   = stemnpp  - rtemp
      ENDDO
      s_sr = resp - s_sr
      resp_s = resp - resp_s

      DO ij=1,stemls-1
        stemv(stemls+1-ij) = stemv(stemls-ij)
      ENDDO

      IF ((nppstore.GT.0.0d0).and.(daynpp.GT.0.0d0)) THEN
        yy       = nppstore*p_stemfr
        stemv(1) = yy
        s_sn     = s_sn + yy
        stemnpp  = stemnpp + yy
        nppstore = nppstore - yy
      ELSE
        stemv(1) = 0.0d0
      ENDIF

      IF (stemnppr.GT.0d0) THEN
        IF ((stemv(1).LT.p_stmin)) then
          stemnpp  = stemnpp + p_stmin - stemv(1)
          nppstore = nppstore - p_stmin + stemv(1)
          s_sn     = s_sn + p_stmin - stemv(1)
          stemv(1) = p_stmin
        ENDIF
      ENDIF

*      print*,mnth,day,stemv(1),stemnppr*respref

*----------------------------------------------------------------------*

      IF (rlai.GT.laimax)  laimax = rlai
      IF (.NOT.(rlai.GT.0.0d0))  daysoff = daysoff + 1.0d0

      ! if store is zero or less kill all live pools
      IF (nppstore.LT.0.0d0) THEN
        DO i=1,leafls
          leaflit = leaflit + leafv(i)
          leafv(i) = 0.0d0
        ENDDO 
c added by Ghislain 07/10/03
        stemlit = 0.0d0
        DO i=1,stemls
          stemlit = stemlit + stemv(i)
          stemv(i) = 0.0d0
        ENDDO 
c added by Ghislain 07/10/03
        rootlit = 0.0d0
        DO i=1,rootls
          rootlit = rootlit + rootv(i)
          rootv(i) = 0.0d0
        ENDDO 

        nppstore = 0.0d0
        rlai = 0.0d0

      ELSE
        stemlit = 0.0d0
        rootlit = 0.0d0
      ENDIF

      ans = 0.0d0
      DO i=1,stemls
        ans = ans + stemv(i)
      ENDDO 

      if ((day == 1).and.(mnth == 1)) then
        sumsr=0.0
        sumrr=0.0
        sumlr=0.0
        summr=0.0
      endif

      sumsr = sumsr + 12.0*s_sr
      sumrr = sumrr + 12.0*s_rr
      sumlr = sumlr + 12.0*s_lr
      summr = summr + 12.0*s_mr

      if ((day == -30).and.(mnth == 12)) 
     &print'(4f10.3)',sumsr,sumrr,sumlr,summr

      RETURN
      END

*----------------------------------------------------------------------*
*                          LAIALT                                      *
*----------------------------------------------------------------------*
      SUBROUTINE LAIALT(leafv,rlai,laiinc,leafls,leaflit)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      REAL*8 leafv(3600),rlai,laiinc,sum,leaflit
      INTEGER leafls,i

      leaflit = leafv(leafls)
      DO i=1,leafls-1
        leafv(leafls+1-i) = leafv(leafls-i)
      ENDDO
      leafv(1) = 0.0d0

      IF (laiinc.GT.0.0d0) THEN
        leafv(1) = laiinc
      ELSEIF (laiinc+leaflit.LT.0.0d0) THEN
        sum = laiinc+leaflit
        i = leafls
20      CONTINUE
          sum = sum + leafv(i)
          IF (sum.GT.1.0e-6) then
            leaflit = leaflit + leafv(i) - sum
            leafv(i) = sum
          ELSE
            leaflit = leaflit + leafv(i)
            leafv(i) = 0.0d0
            i = i - 1
          ENDIF
        IF ((.NOT.(sum.GT.0.0d0)).AND.(i.GT.0))  GOTO 20
      ENDIF

      rlai = 0.0d0
      DO i=1,leafls
        rlai = rlai +  leafv(i)
      ENDDO


      RETURN
      END

