*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE COVER2                           *
*                          ****************                            *
*----------------------------------------------------------------------*
      SUBROUTINE COVER2(nft,ftmor,ftppm0,cov,bio,bioleaf,nppstore,
     &npp,nps,tmp,prc,slc,rlc,c3old,c4old,firec,ppm,hgt,
     &fireres,fprob,fprob_prescr,ftprop,ftstmx,stemdp,rootdp,ftsls,
     &ftrls,ilanduse,prescr_fire,nat_map,ic0,burn,harvest,leafdp,
     &flulccc,ftphen,atprop2,atharvest,aggmap_SDGVM_to_aggHyde,
     &ftprop_init,yield,lat,ftprops,debug)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      INTEGER, PARAMETER :: n_at = 7 !number of aggregated (Hyde, functional) types
      INTEGER, PARAMETER :: NS = 16 !number of SDGVM functional types
      LOGICAL :: debug !used to print out more debugging info
      REAL*8 cov(maxage,maxnft),bio(maxage,2,maxnft),bioleaf(maxnft)
      REAL*8 nppstore(maxnft),npp(maxnft),nps(maxnft),tmp(12),prc(12)
      REAL*8 slc(maxnft),rlc(maxnft),firec
      REAL*8 ppm(maxage,maxnft),hgt(maxage,maxnft),ftppm0(maxnft),told
      REAL*8 ftstmx(maxnft),stemdp(1000,maxnft),rootdp(1000,maxnft)
      REAL*8 fprob,ftprop(maxnft),tot_ngcov,ngcov(maxnft)
      REAL*8 gold,c3old,c4old,fri,norm
      REAL*8 grassrc,ic0(8),sumc,leafdp(3600,maxnft)
      REAL*8 loss,sum_cov(maxnft),flulccc(maxnft)
      REAL*8 lossfrac,lossfrac_nowoodh
      REAL*8 cneed, cneed_leaf, cneed_store
      REAL*8 ftpropnew,ftprop3(maxnft)
      REAL*8 sum_cov_test(maxnft)
      REAL*8 atprop2(n_at,n_at),THRESH
      REAL*8 atharvest(n_at),fprob_prescr(12)
      REAL*8 ft2frac(maxnft),at2prop,woodh,totft,loss_nowoodh
      REAL*8 ftprop_init(maxnft),yield(maxnft),lat
      REAL*8 ftprops(maxnft)
      INTEGER ftsls(maxnft),ftrls(maxnft),nft,ftmor(maxnft),year,i,j
      INTEGER ft,fireres,ilanduse,nat_map(8),age,ftphen(maxnft)
      INTEGER ft2,at,at2,aggmap_SDGVM_to_aggHyde(NS)
      INTEGER ft3,at3,mm
      LOGICAL burn,harvest,prescr_fire,fprob_allok
      LOGICAL compute_covchange,change_cover,corrct_ft,corrct_ft2

      IF(debug .EQV. .TRUE.) THEN
        PRINT '(A)','GG4a ftprop '
        PRINT '(16F11.6)',ftprop(1:nft)
        PRINT '(A)','GC4a cov '
        DO j=1,6 !show first six years of cover
          PRINT '(16F11.6)',cov(j,1:nft)
        ENDDO
      ENDIF


      IF (ilanduse.eq.2) THEN
*----------------------------------------------------------------------*
* Add up grass coverage.                                               *
*----------------------------------------------------------------------*
        gold = 0.0d0
        DO ft=3,4
          DO age=1,ftmor(nat_map(ft))
            gold = gold + cov(age,nat_map(ft))
          ENDDO
        ENDDO

*----------------------------------------------------------------------*
* Add up tree coverage.                                                *
*----------------------------------------------------------------------*
        told = 0.0d0
        DO ft=5,8
          DO age=1,ftmor(nat_map(ft))
            told = told + cov(age,nat_map(ft))
          ENDDO
        ENDDO

*        IF (gold.GT.0.0d0) THEN
*          c3old = cov(1,2) + cov(2,2)
*          c4old = cov(1,3) + cov(2,3)
*        ENDIF

      ENDIF

*----------------------------------------------------------------------*
* Reduce ft area due to change in the LULCC database                   *
* and put this as bare ground ready for new growth 'ngrowth'.          *
*----------------------------------------------------------------------*
      flulccc(:) = 0d0
      tot_ngcov = 0.0d0
      ngcov(:) = 0.0d0

      ! HYDE aggregated land cover: 1 primf, 2 primn, 3 secdf, 4 secdn, 5 C3 crop, 6 C4 crop, 7 urban
C The following ordering is the order of ft's in the input.dat file
      !                                    ft       at
      !            aggmap_SDGVM_to_aggHyde(1)     =  0 !BARE 
      !            aggmap_SDGVM_to_aggHyde(2)     =  7 !urban 
      !            aggmap_SDGVM_to_aggHyde(3:4)   =  2 !C3p, C4p 
      !            aggmap_SDGVM_to_aggHyde(5)     =  5 !C3crop 
      !            aggmap_SDGVM_to_aggHyde(6)     =  6 !C4crop 
      !            aggmap_SDGVM_to_aggHyde(7:8)   =  4 !C3s, C4s
      !            aggmap_SDGVM_to_aggHyde(9:12)  =  1 !Ev_Bp, Dc_Bp, Ev_Np, Dc_Np 
      !            aggmap_SDGVM_to_aggHyde(13:16) =  3 !Ev_Bs, Dc_Bs, Ev_Ns, Dc_Ns  
      IF (ilanduse.ne.2 ) THEN
      sum_cov(:)     = 0.0d0
      loss           = 0.0d0
      ftprop_init    = ftprop
      DO ft=2,nft
        !print*, 'G2',ft,ftphen(ft)
        !PRINT*, 'G2a',ft,ftprop(ft) 
        CALL CALC_CORRCT_FT2(ft,lat,corrct_ft)

        DO age=1,ftmor(ft)
          sum_cov(ft) = sum_cov(ft) + cov(age,ft)
        ENDDO
        IF(debug .EQV. .TRUE.) THEN
          print*, 'ft,SUM_COV(ft),ftprop(ft):',ft, sum_cov(ft),
     &             ftprop(ft)*1d-2
        ENDIF
       

        IF ( (ilanduse.GE.4) .AND. (ilanduse.LE.6) ) THEN 
          ! do for both ftphen(ft) == 1 and 2
          compute_covchange = .TRUE. !compute cover change 
        ELSE
          !IF ( ftphen(ft).EQ.2 ) THEN
          !! need to make this tree specific 
          !  compute_covchange = .TRUE.
          !ELSE
          compute_covchange = .FALSE.
          !ENDIF
        ENDIF

           
!        if((sum_cov(ft).le.0d0)) then !PCM
!          sum_cov(ft) = ftprop(ft)
!        endif
!        if((compute_covchange.EQV..TRUE.).AND.(sum_cov(ft).gt.0d0)) then!PCM Site=11 has troubles here
!       we need to handle the sum_cov(ft).eq.0d0 case as well

        if((compute_covchange.EQV..FALSE.)) then
          woodh = 0.0d0
        else 
          if(ftprop(ft).GT.0) then
             CALL ADJUST_FTPROP(nft,n_at,NS,ft,
     &            sum_cov,ftprop,ftprop_init,ftprops,
     &            aggmap_SDGVM_to_aggHyde,ft2frac,
     &            atprop2,corrct_ft,corrct_ft2,lat,woodh,debug)
          end if
        endif !compute_covchange !PCM move this endif here, so that ngcov gets calculated right 

        !IF(ftprop(ft).GE.0.0d0) THEN
        ! loss > 0 if there is a loss in cover; loss < 0 if there is a gain in cover

        loss = sum_cov(ft) - ftprop(ft)*1d-2
        loss_nowoodh = sum_cov(ft) - (ftprop(ft)+woodh)*1d-2

        !ELSE !lose all cover if ftprop(ft).LT.0.0d0
        !  loss = sum_cov(ft)
        !  loss_nowoodh = sum_cov(ft) 
        !  !ftprop(ft) = 0.0d0
        !ENDIF

        IF((loss.GT.0.0d0) .and. (sum_cov(ft).GT.0.0d0)) THEN
            lossfrac         = MIN(loss/sum_cov(ft),1.0d0)
            lossfrac_nowoodh = MIN(loss_nowoodh/sum_cov(ft),1.0d0)
            CALL LULCC_LOSS2(nft,ftmor,cov,ppm,bio,bioleaf,nppstore,hgt,
     &lossfrac,lossfrac_nowoodh,npp,nps,slc,rlc,fireres,flulccc,harvest,
     &leafdp,sum_cov,ft,debug)
        ENDIF

         !ngcov(ft) = MAX(-loss_nowoodh,0.d0) * sum_cov(ft) !new growth only for -loss>0
        IF(loss.LT.0.0d0) THEN
           !ngcov(ft) = MAX(-loss,0.d0) * sum_cov(ft) !new growth only for -loss>0
           ngcov(ft) = -loss !new growth only for -loss>0
           tot_ngcov = tot_ngcov + ngcov(ft) 
        ENDIF
        IF (debug .EQV. .TRUE.) THEN
! ft SUM_COV(ft):           5   9.6228569383682772E-011
!GG0 5 9.990000 0.000000   -9.990000   -9.990000 0.000000 9.990000 0.000000 0.000000
! ft SUM_COV(ft):           3  0.42309430828850630     
!GG0 3 0.620000 0.000000   -0.196906   -0.196906 0.423094 0.196906 0.211759 0.211335
           PRINT '(A I2 F9.6 F9.6 F12.6 F12.6 F9.6 F9.6 F9.6 F9.6)',
     &              'GG0',ft,
     &              ftprop(ft)*1d-2, woodh*1d-2,
     &              loss,loss_nowoodh,
     &              sum_cov(ft), ngcov(ft),
     &              cov(1,ft),cov(2,ft)
        ENDIF
!         ENDIF
!PCM        endif !compute_covchange !PCM move this endif, so that ngcov gets calculated right 

      ENDDO
      ENDIF

      WHERE(ftprop(:).GT.100.0d0)
         ftprop = 100.0
      ENDWHERE

      WHERE(ftprop(:).LT.0.0d0)
         ftprop = 0.0
      ENDWHERE

      IF(debug .EQV. .TRUE.) THEN
        PRINT '(A)','GG4b ftprop '
        PRINT '(16F11.6)',ftprop(1:nft)
        PRINT '(A)','GS4b sum_cov '
        PRINT '(16F11.6)',sum_cov(1:nft)
        PRINT '(A)','GN4b tot_ngcov '
        PRINT '(1F11.6)',tot_ngcov
        PRINT '(A)','GN4b ngcov '
        PRINT '(16F11.6)',ngcov(1:nft)
      ENDIF

      IF(debug .EQV. .TRUE.) THEN
        PRINT '(A)','GF4b prc '
        PRINT '(12F11.6)',prc
        PRINT '(A)','GF4b tmp '
        PRINT '(12F11.6)',tmp
      ENDIF

      IF(prescr_fire .EQV. .TRUE.) THEN
        fprob_allok = .TRUE.
        DO mm=1,12
         IF( fprob_prescr(mm).GT.1.0d0
     & .OR.  fprob_prescr(mm).LT.0.0d0 ) THEN
          fprob_allok = .FALSE.
          PRINT '(A)',
     & 'fprob_prescr is not in bounds. For the month:'
          PRINT '(I4,12F11.6)',mm,fprob_prescr(mm)
          CONTINUE
!          PRINT '(A)','Stopping.'
!          STOP
         ENDIF
        END DO
!PCM: for now, use equal-weighted average of all months 
        IF(fprob_allok) THEN
           fprob = SUM(fprob_prescr(1:12)/12.0d0)
        ELSE
           RETURN
        END IF
        fri   = -1.0d0
      ELSE
*----------------------------------------------------------------------*
* Compute the likelihood of fire in the current year 'fprob'.          *
* 'fri' is the fire index                                              *
*----------------------------------------------------------------------*
        CALL FIRE2(prc,tmp,fri,fprob,burn)
      ENDIF

      IF(debug .EQV. .TRUE.) THEN
        PRINT '(A)','GF4b fri fprob'
        PRINT '(F11.6 F11.6)',fri,fprob
        PRINT '(A)','GF4b burn '
        PRINT *,burn
      END IF

*----------------------------------------------------------------------*
* Take off area burnt by fire together with plants past there sell by  *
* date and put this as bare ground ready for new growth 'ngrowth'.     *
* Also shift cover and biomass arrays one to the right.                *
*----------------------------------------------------------------------*
      CALL NEWGROWTH2(nft,ftmor,cov,ppm,bio,bioleaf,nppstore,hgt,fprob,
     &npp,nps,tot_ngcov,ngcov,slc,rlc,fireres,firec,harvest,leafdp,
     &flulccc,atharvest,n_at,aggmap_SDGVM_to_aggHyde,NS,yield,ft2frac,
     &sum_cov,debug)

*----------------------------------------------------------------------*

      IF (ilanduse.eq.2) THEN
*----------------------------------------------------------------------*
* Compute percentages of decid, ever and grass for the new growth      *
* ('ngrowth') this year.                                               *
*----------------------------------------------------------------------*
*       CALL PERC2(dof,tmin,ftprop,nft)

*----------------------------------------------------------------------*
* Restrict the rate of trees taking over grassland.                    *
*----------------------------------------------------------------------*
       CALL GRASSREC2(nft,ftprop,gold,tot_ngcov,ngcov,grassrc,nat_map,
     &ilanduse)

*----------------------------------------------------------------------*
* Restrict the rate of bare ground reclimation 0.2, means 20% can be   *
* reclaimed every year and set bare ground cover array.                *
*----------------------------------------------------------------------*
*       CALL BAREREC2(ftprop,cov,bpaold,ngcov,barerc)

*----------------------------------------------------------------------*
* Compute c4 c3 grass split.                                           *
*----------------------------------------------------------------------*
*       CALL c3c4_2(ftprop,c3old,c4old,npp,nps)
      ENDIF

*----------------------------------------------------------------------*
* Set cover arrays to adjust to ftprop as best they can.               *
* ftprop contains the total proportion of that cover, not the          *
* proportion of bare land to assign. Calculate ftprop3.                *
* ftprop3 contains the proportion of the new growth ngcov to assign to *
* each ft.                                                             *
*----------------------------------------------------------------------*
      IF(debug .EQV. .TRUE.) THEN
        PRINT '(A)','GN4c tot_ngcov '
        PRINT '(1F11.6)',tot_ngcov
        PRINT '(A)','GN4c ngcov '
        PRINT '(16F11.6)',ngcov(1:nft)
        PRINT '(A)','GC4c cov '
        DO j=1,6 !show first six years of cover
          PRINT '(16F11.6)',cov(j,1:nft)
        ENDDO

        PRINT '(A)','GG4c ftprop '
        PRINT '(16F11.6)',ftprop(1:nft)
      ENDIF
      norm = 0.0d0
      DO ft=1,nft           
        !PRINT*, 'G7',ft,ftprop(ft) 
        IF (ftprop(ft).GT.0.0d0) THEN
          sum_cov(ft) = 0.0d0
*Start at age=2, since cov array has already been shifted in time by one year. *PCM
*Ignoring age=1 means that we aren't including crop ft's in this normalization, as before.
          DO age=2,ftmor(ft)
            sum_cov(ft) = sum_cov(ft) + cov(age,ft)
          ENDDO
*the variable ftprop3(ft) will contain the proportion of new growth land to assign 
* for each ft 
          ftprop3(ft)=ftprop(ft)/100.0d0 - sum_cov(ft)
          IF (ftprop3(ft).LT.0.0d0) ftprop3(ft)=0.0d0 !can't remove cov
          norm = norm + ftprop3(ft)
        ELSE
          ftprop3(ft)=0.0d0
        ENDIF
      END DO

      DO ft=1,nft
        ftprop3(ft) = 100.0d0*ftprop3(ft)/norm
      ENDDO

      IF(debug .EQV. .TRUE.) THEN
        PRINT '(A)','GG4d ftprop '
        PRINT '(16F11.6)',ftprop(1:nft)
        PRINT '(A)','GS4d sum_cov '
        PRINT '(16F11.6)',sum_cov(1:nft)
      ENDIF

      IF(ilanduse .LT. 3) THEN
         cov(1,1) = tot_ngcov*ftprop3(1)/100.0d0
      ELSE
         !cov(1,1) = ngcov(1)
         ftprop(1) = MIN(MAX(100.0d0-SUM(ftprop(2:nft)),0.0d0),100.0d0) 
         cov(1,1) = ftprop(1)/100.0d0
      ENDIF
      !PRINT*, 'G8' 

*----------------------------------------------------------------------*
* Set cover arrays for this years ft proportions, take carbon from     *
* litter to provide nppstore and canopy.                               *
*----------------------------------------------------------------------*
      DO ft=2,nft
        !PRINT*, 'G9',ft,ftprop3(ft) 
        IF (ftprop3(ft).GT.0.0d0) THEN

*----------------------------------------------------------------------*
* If no veg exists then set nppstore to initial value and reinitialise.*
*----------------------------------------------------------------------*
          year = 1
          !PRINT*, 'G10' 
10        CONTINUE
            IF (cov(year,ft).GT.0.0d0) THEN
              goto 20
            ELSE
              year = year + 1
              IF (year.EQ.ftmor(ft)+1) THEN
                nppstore(ft) = ftstmx(ft)
                DO i=1,ftsls(ft)
                  stemdp(i,ft) = 0.0d0
                ENDDO
                DO i=1,ftrls(ft)
                  rootdp(i,ft) = 0.0d0
                ENDDO
                GOTO 20
              ENDIF
            ENDIF
          GOTO 10
20        CONTINUE
*----------------------------------------------------------------------*

          IF(ilanduse .LT. 3) THEN
            cov(1,ft) = ftprop3(ft)*tot_ngcov/100.0d0
          ELSE
            cov(1,ft) = ngcov(ft)
          ENDIF

          IF(cov(1,ft).GT.0.0d0) THEN !only when new growth
            ppm(1,ft) = ftppm0(ft)
            hgt(1,ft) = 0.004d0

            ! take new cohort carbon from stem litter (slc), then soil C
            cneed_leaf   = bioleaf(ft)*cov(1,ft)
            cneed_store  = nppstore(ft)*cov(1,ft)
            cneed        = cneed_leaf + cneed_store
            slc(ft) = slc(ft) - cneed 
            !slc(ft) = slc(ft) - (nppstore(ft) + bioleaf(ft))*cov(1,ft)
            IF (slc(ft).LT.0.0d0) THEN
              sumc = 0.0d0
              DO i=1,8
                 sumc = sumc + ic0(i)
              ENDDO

              ! if soil C is insufficient to support required C do not use 
              IF ((slc(ft)+sumc) .LT. 0.0d0 ) THEN
                nppstore(ft) = nppstore(ft) + slc(ft)*cneed_store/cneed 
                bioleaf(ft)  = bioleaf(ft)  + slc(ft)*cneed_leaf/cneed 
                cov(1,ft)    = cov(1,ft) * (1.0d0 + slc(ft)/cneed)  
              ELSE
                DO i=1,8
                  ic0(i) = ic0(i)*(1.0d0+slc(ft)/sumc)
                ENDDO
              ENDIF

              slc(ft) = 0.0d0
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      !loop added for testing purposes
      DO ft=1,nft
        sum_cov_test(ft) = 0.0
        DO age=1,ftmor(ft)
          sum_cov_test(ft) = sum_cov_test(ft) + cov(age,ft)
        ENDDO
      ENDDO

!      IF (ilanduse.ge.3 .and. ilanduse.le.6) THEN
!       WRITE(*,*) ' BARE       CITY       C3p        C4p        ',
!     &'C3crop     C4crop      C3s        C4s        Ev_Bp       ', 
!     &'Ev_Np      Dc_Bp        Dc_Np      Ev_Bs      Ev_Ns      ', 
!     &'Dc_Bs      Dc_Ns'
!      ELSE
!       WRITE(*,*) ' BARE       CITY       C3         C4         ',
!     &'C3crop     C4crop       Ev_Bl       ', 
!     &'Ev_Nl      Dc_Bl        Dc_Nl       ' 
!      ENDIF
!
      IF(debug .EQV. .TRUE.) THEN
        PRINT '(A)','GG3 COV AGE=1 '
        PRINT '(16F11.6)',cov(1,1:nft)
        PRINT '(A)','GG4 COV '
        PRINT '(16F11.6)',sum_cov_test(1:nft)
!return ftprop for using in gross transitions for next year 
        PRINT '(A)','GG4e ftprop '
        PRINT '(16F11.6)',ftprop(1:nft)
        PRINT '(A)','GH4e ftprop3 '
        PRINT '(16F11.6)',ftprop3(1:nft)
        !PRINT '(A)','GG5 BIOL'
        !PRINT '(16F11.6)',bioleaf(1:nft)
        !PRINT '(A)','GG6 NPPS'
        !PRINT '(16F11.6)',nppstore(1:nft)
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE ADJUST_FTPROP                    *
*                          *****************                           *
*----------------------------------------------------------------------*
      SUBROUTINE ADJUST_FTPROP(nft,n_at,NS,ft,
     &  sum_cov,ftprop,ftprop_init,ftprops,
     &  aggmap_SDGVM_to_aggHyde,ft2frac,atprop2,corrct_ft,
     &  corrct_ft2,lat,woodh,debug)

      INCLUDE 'array_dims.inc'
      LOGICAL :: debug !used to print out more debugging info
      INTEGER nft,ft,n_at,NS
      INTEGER aggmap_SDGVM_to_aggHyde(NS)
      LOGICAL corrct_ft,corrct_ft2
      REAL*8 sum_cov(maxnft)
      REAL*8 ftprop(maxnft),ftprop_init(maxnft),lat,ftprops(maxnft)
      REAL*8 atprop2(n_at,n_at)
      REAL*8 atharvest(n_at)
      REAL*8 ft2frac(maxnft),woodh

      INTEGER at,at2,ft2
      REAL*8 tot_secdn,barefrac_secdn,nonbarefrac

      at = aggmap_SDGVM_to_aggHyde(ft)
      CALL CALC_FT2FRAC2(ftprop_init,aggmap_SDGVM_to_aggHyde,nft,
     &ft2frac,debug)

      if(at.EQ.4) then
        !use ftprops = states vector from net transitions; this allows for transitions from bare secondary non-forest to non-bare secondary non-forest 
        tot_secdn = ftprops(1)+ftprops(7)+ftprops(8)
        if(tot_secdn.GT.0.0d0)then
          barefrac_secdn = ftprops(1)/tot_secdn
        else
          barefrac_secdn = 0.0d0
        endif
        nonbarefrac = 1.0d0-barefrac_secdn
      else
        nonbarefrac = 1.0d0 
      endif

      !if(at.EQ.1 .OR. at.EQ.2) then
      if(at.EQ.1) then
      ! losses to ft from wood harvest in ft 
        woodh = ft2frac(ft)*atharvest(at)
        ftprop(ft) = ftprop(ft) - woodh
      else
        woodh = 0.0d0
      endif

      ! losses to ft from conversion to urban cover (at2.EQ.7)
      ftprop(ft) = ftprop(ft)-ft2frac(ft)*atprop2(at,7)

      ! gains to ft from conversion from urban cover (at2.EQ.7)
      ftprop(ft) = ftprop(ft)+ft2frac(ft)*atprop2(7,at)

      DO ft2=2,nft
        at2 = aggmap_SDGVM_to_aggHyde(ft2)
        CALL CALC_CORRCT_FT2(ft2,lat,corrct_ft2)

         
        if(at2.EQ.4 .AND. at.EQ.2 ) then
        ! gains to ft2 from wood harvest in ft 
        !   for at2.EQ.4 (secdn) from at.EQ.2 (primn)
          !ftprop(ft2) = ftprop(ft2)+ft2frac(ft2)*woodh
          ftprop(ft2) = ftprop(ft2)+ft2frac(ft2-4)*woodh !assumes primary cover fraction
        endif

        if(at2.EQ.3 .AND. at.EQ.1 ) then
        ! gains to ft2 from wood harvest in ft 
        !   for at2.EQ.3 (secdf) from at.EQ.1 (primf)
          !ftprop(ft2) = ftprop(ft2)+ft2frac(ft2)*woodh
          ftprop(ft2) = ftprop(ft2)+ft2frac(ft2-4)*woodh !assumes primary cover fractions
        endif

        ! losses from ft to ft2: !atprop2 additive in %/year
        IF(ft2frac(ft2).GT.0.0d0) THEN !check for pre-existing cover in the ft2
        ! split the losses from each ft by a fraction ft2frac(ft)
        ! split the losses to each ft to each ft2 by a fraction ft2frac(ft2)
          ftprop(ft) = ftprop(ft) -
     &  ft2frac(ft)*ft2frac(ft2)*atprop2(at,at2)
        ELSE IF(corrct_ft2) THEN !check if tropical or temperate
          ! split the losses from each ft by a fraction ft2frac(ft)
          ftprop(ft) = ftprop(ft) -
     &  ft2frac(ft)*             atprop2(at,at2)
        ENDIF

        ! gains to ft from ft2:
        IF(ft2frac(ft).GT.0.0d0) THEN ! check for pre-existing cover in the ft !also check if tropical or temperate 
        ! split the gains from each at2 from each ft2 by a fraction ft2frac(ft2)
        ! split the gains to each at to each ft by a fraction ft2frac(ft)
           ftprop(ft) = ftprop(ft) +
     &  ft2frac(ft)*ft2frac(ft2)*atprop2(at2,at)
        ELSE IF((corrct_ft.EQV..TRUE.).AND.
     &(nonbarefrac.GT.0.0d0)) THEN ! check if tropical or temperate !also, only do this if there is non-bare secdn ground or if it's not secdn
           ftprop(ft) = ftprop(ft) +
     &             ft2frac(ft2)*atprop2(at2,at)
        ENDIF

        IF ((debug.EQV..TRUE.).AND.((at.eq.1).OR.(at2.eq.1))) THEN
            PRINT
     &     '(A I2 I2 I3 I3 F11.6 F11.6 F11.6 F11.6 F11.6 F11.6 F11.6)',
     &              'GHG1',
     &              at2,at,ft2,ft,
     &              ft2frac(ft),ft2frac(ft2),atprop2(at2,at),
     &              ftprop(ft),woodh,sum_cov(ft),atharvest(at)
            PRINT
     &     '(A I2 I2 I3 I3 F11.6 F11.6 F11.6 F11.6 F11.6 F11.6 F11.6)',
     &              'GHL1',
     &              at,at2,ft,ft2,
     &              ft2frac(ft),ft2frac(ft2),atprop2(at,at2),
     &              ftprop(ft),woodh,sum_cov(ft),atharvest(at)
        ENDIF
      ENDDO
!         totft = ftprop(ft) + woodh 
!
!         IF( totft .GT. 1d-1 ) THEN
!            loss = 1d0 - totft*1d-2/sum_cov(ft)
!            loss_nowoodh = 1d0 - (totft-woodh)*1d-2/sum_cov(ft)

         !IF( ftprop(ft) .GT. 1d-1 ) THEN !this discretization leads to carbon imbalances
         !IF( ftprop(ft) .GT. 1d-2 ) THEN !this prevents cov<0
         !IF(sum_cov(ft).GT.0.0d0) THEN
         !IF(ftprop(ft) .GT. 0.0d0) THEN
         ! loss = 1.0d0 - ftprop(ft)*1d-2/sum_cov(ft) 
         ! loss > 0 if there is a loss in cover; loss < 0 if there is a gain in cover
         ! loss_nowoodh = 1.0d0 - (ftprop(ft)+woodh)*1d-2/sum_cov(ft)
         !ELSE
         !   !ftloss_prop = 1d0   !PCM why isn't this ftloss_prop(ft) = 1d0 ??
         !   loss = 1d0   !PCM 
         !   loss_nowoodh = 1d0
         !ENDIF
         !ELSE
         !   loss = 0d0   
         !   loss_nowoodh = 0d0
         !ENDIF


!         print*, 'ftprop is less than sum_cov:',
!     &ft, ftprop(ft)*1d-2, sum_cov(ft), loss 
!         print*, (ftprop(ft)*1d-2) - sum_cov(ft)

!         change_cover = .False. 
!         change_cover = .True. 
!!         IF(sum_cov(ft).GT.0.0d0) THEN
!         IF( ilanduse .LT. 3 ) THEN
!!           IF( ( (ftprop(ft)*1d-2) - sum_cov(ft)) .lt. -5d-3  ) THEN
!           IF( ( ftprop(ft)*1d-2 - sum_cov(ft)) .ne. 0.0  ) THEN
!             change_cover = .True. !settings for net transitions
!           ENDIF
!         ELSE IF( ilanduse .GE. 3 .AND. ilanduse .LE. 6 ) THEN
!           IF( ftprop(ft)*1d-2 .NE. sum_cov(ft) ) THEN 
!             change_cover = .True. !settings for gross transitions
!           ENDIF
!         ENDIF
!!         ENDIF

!         IF( change_cover ) THEN 
         ! NET transitions: IF( ( (ftprop(ft)*1d-2) - sum_cov(ft)) .lt. -5d-3  ) THEN
         ! GROSS transitions IF( ftprop(ft)*1d-2 .LT. sum_cov(ft) ) THEN 

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE GROWTH2
*                          *****************                           *
*----------------------------------------------------------------------*
      SUBROUTINE GROWTH2(nft,ftmor,ftwd,ftxyl,ftpd,ftgr0,ftgrf,cov,bio,
     &nppstore,npp,lai,nps,npr,evp,slc,rlc,sln,rln,stembio,rootbio,ppm,
     &hgt,leaflit,debug)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftwd(maxnft),ftxyl(maxnft),ftpd(maxnft),cov(maxage,maxnft)
      REAL*8 bio(maxage,2,maxnft),npp(maxnft),lai(maxnft),nps(maxnft)
      REAL*8 npr(maxnft),evp(maxnft),nppstore(maxnft),ftgr0(maxnft)
      REAL*8 rootbio,slc(maxnft),rlc(maxnft),sln(maxnft),rln(maxnft)
      REAL*8 stembio,ppm(maxage,maxnft),hgt(maxage,maxnft),ftgrf(maxnft)
      REAL*8 leaflit(maxnft),ftmat(maxnft)
      INTEGER nft,ftmor(maxnft),age,ft,i
      LOGICAL debug

*----------------------------------------------------------------------*
* Initialise litter arrays, and add on leaf litter computed in DOLY.   *
*----------------------------------------------------------------------*
      DO ft=1,nft
        slc(ft) = 0.0d0
        rlc(ft) = 0.0d0
        sln(ft) = 0.0d0
        rln(ft) = 0.0d0
      ENDDO

      DO ft=1,nft
        DO age=1,ftmor(ft)
          slc(ft) = slc(ft) + leaflit(ft)*cov(age,ft)
        ENDDO
      ENDDO

*----------------------------------------------------------------------*
* Compute maturity of plant fts.                                       *
*----------------------------------------------------------------------*
*      CALL VEGMAT2(nft,ftmor,ftwd,ftxyl,ftpd,evp,lai,npp,
*     &     nps,ftmat,maxnft)

      DO ft=1,nft
        ftmat(ft) = real(ftmor(ft))
      ENDDO

*----------------------------------------------------------------------*
* Thin vegetation where npp is not sufficient to maintain sensible     *
* growth rate.                                                         *
*----------------------------------------------------------------------*
      CALL THIN2(nft,ftmor,ftmat,ftwd,ftxyl,ftpd,ftgr0,ftgrf,cov,bio,
     &nppstore,npp,lai,nps,evp,slc,ppm,hgt,debug)

*----------------------------------------------------------------------*
* Add on biomasses for the year.                                       *
*----------------------------------------------------------------------*
      CALL ADDBIO2(nft,ftmor,npp,nps,npr,bio)

*----------------------------------------------------------------------*
* Compute litter arrays.                                               *
*----------------------------------------------------------------------*
      IF(debug .EQV. .TRUE.) THEN
        DO ft=2,nft
           print*,'BEFORE MKLIT: ft,slc,rlc: ',ft,slc(ft),rlc(ft)
        ENDDO
        DO ft=1,nft
           print*,'BEFORE MKLIT: ft,ftgr0,ftmor',ft,ftgr0(ft),ftmor(ft)
        ENDDO
      ENDIF
      CALL MKLIT2(nft,ftgr0,ftmor,ftmat,cov,bio,slc,rlc,sln,rln,npp,nps,
     &npr,debug)
      IF(debug .EQV. .TRUE.) THEN
        DO ft=2,nft
           print*,'AFTER MKLIT: ft,slc,rlc: ',ft,slc(ft),rlc(ft)
        ENDDO
      ENDIF

*----------------------------------------------------------------------*
*Compute leaf root and stem biomasses.                                 *
*----------------------------------------------------------------------*
      stembio = 0.0d0
      rootbio = 0.0d0
      DO ft=1,nft
        DO i=1,ftmor(ft)
          stembio = stembio + bio(i,1,ft)*cov(i,ft)
          rootbio = rootbio + bio(i,2,ft)*cov(i,ft)
        ENDDO
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE THIN2                            *
*                          ***************                             *
*----------------------------------------------------------------------*
      SUBROUTINE THIN2(nft,ftmor,ftmat,ftwd,ftxyl,ftpd,ftgr0,ftgrf,cov,
     &bio,nppstore,npp,lai,nps,evp,slc,ppm,hgt,debug)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftwd(maxnft),ftmat(maxnft),ftxyl(maxnft),ftpd(maxnft)
      REAL*8 cov(maxage,maxnft),bio(maxage,2,maxnft),npp(maxnft)
      REAL*8 lai(maxnft),nps(maxnft),evp(maxnft),nppstore(maxnft)
      REAL*8 storelit(maxnft),slc(maxnft),ppm(maxage,maxnft)
      REAL*8 hgt(maxage,maxnft),ftgr0(maxnft),ftgrf(maxnft),pbionew
      REAL*8 pbioold,no,pbio,ftcov(maxnft),covnew(maxage,maxnft)
      REAL*8 ppmnew(maxage,maxnft),shv,lmv,nv,pi,hwv(maxage),emv,totno
      REAL*8 totcov,scale(maxage),oldbio(maxage,2,maxnft),dimold,g0,gf
      REAL*8 gm,grate(maxage),stlit2,hgtnew(maxage),dimnew,maxhgt,minhgt
      REAL*8 hc1,hc2,tcov1,tcov2,sum,xxx
      INTEGER nft,ftmor(maxnft),ft,i,year,age
      LOGICAL debug

      pi = 3.1415926d0

      hc1 = 0.05d0
      hc2 = 2.0d0

      DO ft=1,nft
        DO year=1,ftmor(ft)
          DO i=1,2
            oldbio(year,i,ft) = bio(year,i,ft)
          ENDDO
        ENDDO
        storelit(ft) = 0.0d0
      ENDDO

*----------------------------------------------------------------------*
* FT loop for thinning and height competition.                         *
*----------------------------------------------------------------------*
      DO ft=1,nft
!        IF(debug)THEN
!          DO year=1,ftmor(ft)
!                PRINT *,'GT0A',ft,year,cov(year,ft),ppm(year,ft)
!          ENDDO
!        ENDIF
        IF (ftgr0(ft).GT.0.0d0) THEN

        IF (npp(ft)*nps(ft).GT.0.0d0) THEN

*----------------------------------------------------------------------*
* Height competition.                                                  *
*----------------------------------------------------------------------*
          minhgt = 10000.0d0
          maxhgt =-10000.0d0
          i = 0
          DO year=1,ftmor(ft)
            IF (hgt(year,ft).gt.0) i = i + 1
            scale(year) = hgt(year,ft)
            scale(year) = real(i)
            IF (scale(year).LT.minhgt)  minhgt = scale(year)
            IF (scale(year).GT.maxhgt)  maxhgt = scale(year)
          ENDDO

          IF (maxhgt-minhgt.LT.0.0d0) THEN
            DO year=1,ftmor(ft)
              scale(year) = (scale(year) - minhgt)/(maxhgt - minhgt)
              scale(year) = (scale(year)*hc1 + 1.0d0 - hc1)**hc2
            ENDDO

            tcov1 = 0.0d0
            tcov2 = 0.0d0
            DO year=1,ftmor(ft)
              covnew(year,ft) = cov(year,ft)*scale(year)
              tcov1 = tcov1 + cov(year,ft)
              tcov2 = tcov2 + covnew(year,ft)
            ENDDO

            DO year=1,ftmor(ft)
!              IF(debug)THEN
!                PRINT *,'GT0B',
!     &ft,year,cov(year,ft),covnew(year,ft),ppm(year,ft)
!              ENDIF
              IF (covnew(year,ft).GT.0.0d0) THEN
              covnew(year,ft) = covnew(year,ft)*tcov1/tcov2
              DO i=1,2
                bio(year,i,ft) = bio(year,i,ft)*cov(year,ft)/
     &covnew(year,ft)
                oldbio(year,i,ft) = oldbio(year,i,ft)*cov(year,ft)/
     &covnew(year,ft)
              ENDDO
              ppm(year,ft) = ppm(year,ft)*cov(year,ft)/covnew(year,ft)
              cov(year,ft) = covnew(year,ft)
              ENDIF
!              IF(debug)THEN
!                PRINT *,'GT1A',
!     &ft,year,cov(year,ft),covnew(year,ft),ppm(year,ft)
!              ENDIF
            ENDDO

          ENDIF
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
* Height competition.                                                  *
*----------------------------------------------------------------------*

          minhgt = 10000.0d0
          maxhgt =-10000.0d0
          i = 0
          DO year=1,ftmor(ft)
            IF (hgt(year,ft).gt.0) i = i + 1
            scale(year) = 1.0d0 - cov(year,ft)
            IF (scale(year).LT.minhgt)  minhgt = scale(year)
            IF (scale(year).GT.maxhgt)  maxhgt = scale(year)
          ENDDO

          IF (maxhgt-minhgt.GT.0.0d0) THEN
            DO year=1,ftmor(ft)
              scale(year) = (scale(year) - minhgt)/(maxhgt - minhgt)
              scale(year) = (scale(year)*hc1 + 1.0d0 - hc1)**hc2
            ENDDO

            tcov1 = 0.0d0
            tcov2 = 0.0d0
            DO year=1,ftmor(ft)
              covnew(year,ft) = cov(year,ft)*scale(year)
              tcov1 = tcov1 + cov(year,ft)
              tcov2 = tcov2 + covnew(year,ft)
            ENDDO

            DO year=1,ftmor(ft)
              IF (covnew(year,ft).GT.0.0d0) THEN
              covnew(year,ft) = covnew(year,ft)*tcov1/tcov2
              DO i=1,2
                bio(year,i,ft) = bio(year,i,ft)*cov(year,ft)/
     &covnew(year,ft)
                oldbio(year,i,ft) = oldbio(year,i,ft)*cov(year,ft)/
     &covnew(year,ft)
              ENDDO
              ppm(year,ft) = ppm(year,ft)*cov(year,ft)/covnew(year,ft)
              cov(year,ft) = covnew(year,ft)
              ENDIF
!              IF(debug)THEN
!                PRINT *,'GT1B',
!     &ft,year,cov(year,ft),covnew(year,ft),ppm(year,ft)
!              ENDIF
            ENDDO

          ENDIF
*----------------------------------------------------------------------*

          g0 = ftgr0(ft)
          gf = ftgrf(ft)
          gm = real(ftmor(ft))/10.0
          DO year=1,ftmor(ft)
            grate(year) = (gf - g0)/gm*real(year - 1) + g0
            IF (real(year).GE.gm)  grate(year) = gf
          ENDDO

          totno = 0.0d0
          totcov = 0.0d0
          DO year=1,ftmor(ft)
            totno = totno + cov(year,ft)*ppm(year,ft)
            totcov = totcov + cov(year,ft)
          ENDDO
          totno = totno/totcov

*----------------------------------------------------------------------*
* Compute cover and ppm to sustain a minimum growth rate, put these    *
* values in covnew and ppmnew.                                         *
*----------------------------------------------------------------------*
          DO year=1,ftmor(ft)

            shv = npp(ft)*nps(ft)/100.0d0/100.0d0
            emv = evp(ft)/3600.0d0/1000.0d0
            lmv = lai(ft)*1.3d0/(lai(ft) + 3.0d0)
            lmv = 1.0d0
            nv = 1.0d0

            IF (ppm(year,ft)*cov(year,ft).GT.0.0d0) THEN
*----------------------------------------------------------------------*
* Calculate the increase in diameter produced by stem NPP 'nps'.       *
*----------------------------------------------------------------------*

              IF (emv.GT.0.0d0) THEN
! hwv = theoretical maximum height (hydrolics)
                hwv(year) = ftpd(ft)*lmv/nv*(ftxyl(ft)*shv/emv/ftwd(ft)/
     &10000.0d0)**0.5d0
              ELSE
                hwv(year) = hgt(year,ft)*0.9d0
              ENDIF

! Calculate new height.
              IF (hwv(year).GT.0.0d0) THEN
                hgtnew(year) = (hwv(year)-hgt(year,ft))/hwv(year)*0.5d0
              ELSE
                hgtnew(year) = 0.0d0
              ENDIF
              IF (hgtnew(year).LT.0.0d0)  hgtnew(year) = 0.0d0
              hgtnew(year) = hgtnew(year) + hgt(year,ft)

! pbio = g/individual
              pbioold = bio(year,1,ft)/ppm(year,ft)
              pbionew = (bio(year,1,ft) + npp(ft)*nps(ft)*(1.0 - 
     &stlit2(year,ftmat(ft)))/100.0d0)/ppm(year,ft)

!              print'(''pbio'',3f12.1)',pbionew,pbioold,pbionew-pbioold

! Old diameter
              IF (hgt(year,ft).GT.0.0d0) THEN
                dimold = 2.0d0*(pbioold/1000000.0d0/hgtnew(year)/
     &pi/ftwd(ft))**0.5d0
              ELSE
                dimold = 0.0d0
              ENDIF

! New diameter
              IF (hgtnew(year).GT.0.0d0) THEN
                dimnew = 2.0d0*(pbionew/1000000.0d0/hgtnew(year)/
     &pi/ftwd(ft))**0.5d0
              ELSE
                dimnew = 0.0d0
              ENDIF

              IF ((dimnew-dimold)/2.0d0.GT.grate(year)) THEN
*----------------------------------------------------------------------*
* No thinning required.                                                *
*----------------------------------------------------------------------*
                ppmnew(year,ft) = ppm(year,ft)
                covnew(year,ft) = cov(year,ft)
              ELSE
*----------------------------------------------------------------------*
* Thinning required.                                                   *
*----------------------------------------------------------------------*
*               print'(i4,3f8.4)',year,dimnew-dimold,dimnew,dimold
                xxx = (grate(year)+dimold/2.0d0)**2.0d0*1000000.0d0*
     &hgtnew(year)*pi*ftwd(ft)
                ppmnew(year,ft) = (bio(year,1,ft) + npp(ft)*nps(ft)*
     &(1.0d0 - stlit2(year,ftmat(ft)))/100.0d0)/xxx
                covnew(year,ft) = cov(year,ft)

              ENDIF

            ELSE

              pbio = 0.0d0
              scale(year) = 0.0d0
              ppmnew(year,ft) = 0.0d0
              covnew(year,ft) = 0.0d0
              hgtnew(year) = 0.0d0

            ENDIF

          ENDDO

*----------------------------------------------------------------------*
* Correct biomass array, and adjust litter for any thinned trees.      *
*----------------------------------------------------------------------*
          DO year=1,ftmor(ft)
            !IF(debug)THEN
            !  PRINT *,'GT2',
!     &ft,year,cov(year,ft),covnew(year,ft),ppm(year,ft)
            !ENDIF
            IF ((ppm(year,ft).GT.0.0d0).AND.(cov(year,ft).GT.0.0d0))
     &THEN
                no = ppm(year,ft)*cov(year,ft) - 
     &ppmnew(year,ft)*covnew(year,ft)

              DO i=1,2
                pbio = bio(year,i,ft)/ppm(year,ft)
                slc(ft) = slc(ft) + pbio*no

                IF (cov(year,ft).GT.0.0d0) THEN
                    bio(year,i,ft) = (oldbio(year,i,ft)*cov(year,ft) - 
     &pbio*no)/covnew(year,ft)
                    bio(year,i,ft) = oldbio(year,i,ft)*ppmnew(year,ft)/
     &ppm(year,ft)
                ENDIF
              ENDDO
              storelit(ft) = storelit(ft) + nppstore(ft)/ppm(year,ft)*no
!              IF(debug)THEN
!                PRINT *,'GT3',ft,year,cov(year,ft), covnew(year,ft)
!              ENDIF
              cov(year,ft) = covnew(year,ft)
              ppm(year,ft) = ppmnew(year,ft)
              hgt(year,ft) = hgtnew(year)
            ELSE
              cov(year,ft) = 0.0d0
              IF(debug)THEN
                PRINT *,'GT4',ft,year,cov(year,ft), covnew(year,ft)
              ENDIF
              ppm(year,ft) = 0.0d0
              hgt(year,ft) = 0.0d0
            ENDIF

          ENDDO
*----------------------------------------------------------------------*

        ELSE
          DO year=1,ftmor(ft)
            ppm(year,ft) = 0.0d0
            hgt(year,ft) = 0.0d0
          ENDDO
        ENDIF

      ENDIF
*----------------------------------------------------------------------*
* End of ft loop.
*----------------------------------------------------------------------*
      ENDDO

*----------------------------------------------------------------------*
* Correct nppstore to account for thinning.                            *
*----------------------------------------------------------------------*
      DO ft=1,nft
        IF (ftgr0(ft).GT.0.0d0) THEN

        IF (storelit(ft).GT.0.0d0) THEN
          ftcov(ft) = 0.0d0
          DO age=1,ftmor(ft)
            ftcov(ft) = ftcov(ft) + cov(age,ft)
          ENDDO
          IF (ftcov(ft).GT.0.0d0) THEN
            slc(ft) = slc(ft) +  storelit(ft)
            nppstore(ft) = (nppstore(ft)*ftcov(ft) - 
     &storelit(ft))/ftcov(ft)
          ENDIF
        ENDIF

        ENDIF
      ENDDO
*----------------------------------------------------------------------*


      RETURN
      END

*----------------------------------------------------------------------*
*                             FUNCTION find2                           *
*                             *************                            *
*----------------------------------------------------------------------*
      FUNCTION find2(tmp,prc)
*----------------------------------------------------------------------*
      REAL*8 find2,tmp(12),prc(12)
      REAL*8 tmplim,prclim
      INTEGER i

      tmplim = -5.0d0
      prclim = 50.0d0

      find2 = 0.0d0
      DO i=1,12
        IF (tmp(i).GT.tmplim) THEN
          IF (prc(i).LT.prclim) THEN
            find2 = find2 + prc(i)/(12.0d0*prclim)
          ELSE
            find2 = find2 + 1/12.0d0
          ENDIF
        ELSE
          find2 = find2 + 1/12.0d0
        ENDIF
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                             SUBROUTINE c3c42_2                       *
*                             ***************                          *
*----------------------------------------------------------------------*
      SUBROUTINE c3c42_2(ftprop,npp,range)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftprop(maxnft),npp(maxnft),range
      REAL*8 grass,nd

      grass = ftprop(2)
      nd = npp(2) - npp(3)
      ftprop(2) = grass*nd/(2.0d0*range) + grass/2.0d0
      IF (ftprop(2).GT.grass)  ftprop(2) = grass
      IF (ftprop(2).LT.0.0d0)  ftprop(2) = 0.0d0
      ftprop(3) = grass - ftprop(2)


      RETURN
      END

*----------------------------------------------------------------------*
*                             SUBROUTINE c3c4_2                        *
*                             ***************                          *
*----------------------------------------------------------------------*
      SUBROUTINE c3c4_2(ftprop,c3old,c4old,npp,nps)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftprop(maxnft),npp(maxnft),nps(maxnft),c3old,c4old
      REAL*8 grass,c3p,c4p,adj

      IF (c3old+c4old.GT.0.0d0) THEN
        c3p = c3old/(c3old + c4old)
        c4p = c4old/(c3old + c4old)
      ELSE
        c3p = 0.5d0
        c4p = 0.5d0
      ENDIF

      IF (npp(2)*nps(2)+npp(3)*nps(3).GT.0.0d0) THEN
        adj = npp(2)*nps(2)/(npp(2)*nps(2) + npp(3)*nps(3)) - 0.5d0
      ELSE
        adj = 0.0d0
      ENDIF

      c3p = c3p + adj
      c4p = c4p - adj

      IF (c3p.LT.0.0d0) THEN
        c3p = 0.0d0
        c4p = 1.0d0
      ENDIF

      IF (c4p.LT.0.0d0) THEN
        c4p = 0.0d0
        c3p = 1.0d0
      ENDIF

      grass = ftprop(2) + ftprop(3)

      ftprop(2) = grass*c3p
      ftprop(3) = grass*c4p

      RETURN
      END

*----------------------------------------------------------------------*
*                             SUBROUTINE GRASSREC2                     *
*                             *******************                      *
*----------------------------------------------------------------------*
      SUBROUTINE GRASSREC2(nft,ftprop,gold,tot_ngcov,ngcov,x,nat_map,
     &ilanduse)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftprop(maxnft),gold,ngcov(maxnft),tot_ngcov
      REAL*8 x,ntcov,ftt,ftpropo(maxnft)
      INTEGER nft,ft,nat_map(8),ilanduse

      ntcov = 0.0d0
      ftt = 0.0d0
      DO ft=5,8
        IF(ilanduse .LT. 3) THEN
          ntcov = ntcov + ftprop(nat_map(ft))*tot_ngcov/100.0d0
        ELSE
          ntcov = ntcov + ftprop(nat_map(ft))*ngcov(ft)/100.0d0
        ENDIF
        ftt = ftt + ftprop(nat_map(ft))
      ENDDO

      IF (ntcov.GT.gold*x) THEN
        DO ft=4,nft
          ftpropo(ft) = ftprop(ft)
          IF(ilanduse .LT. 3) THEN
            ftprop(ft) = 100.0d0*gold*x*ftprop(ft)/(ftt*tot_ngcov)
          ELSE
            ftprop(ft) = 100.0d0*gold*x*ftprop(ft)/(ftt*ngcov(ft))
          ENDIF
          ftprop(2) = ftprop(2) + ftpropo(ft) - ftprop(ft)
        ENDDO

        ntcov = 0.0d0
        DO ft=4,nft
          IF(ilanduse .LT. 3) THEN
            ntcov = ntcov + ftprop(ft)*tot_ngcov/100.0d0
          ELSE
            ntcov = ntcov + ftprop(ft)*ngcov(ft)/100.0d0
          ENDIF
        ENDDO
        IF (abs(ntcov-gold*x).GT.0.000001d0) WRITE(11,
     &'(''Treerec subroutine error'')')
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                             SUBROUTINE BAREREC2                      *
*                             ******************                       *
*----------------------------------------------------------------------*
      SUBROUTINE BAREREC2(ftprop,cov,bpaold,ngcov,x)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftprop(maxnft),bpaold,ngcov,x,cov(maxage,maxnft)
      REAL*8 nbp,cgcov

* Current growth 'cgcov'
      cgcov = 1.0d0 - ngcov

      nbp = ngcov*ftprop(1)/100.0d0

      IF (nbp.LT.bpaold*(1.0d0-x)) THEN
        ngcov = 1.0d0 - bpaold*(1.0d0 - x) - cgcov
        cov(1,1) = bpaold*(1.0d0 - x)
      ELSE
        cov(1,1) = ngcov*ftprop(1)/100.0d0
        ngcov = ngcov*(1.0d0 - ftprop(1)/100.0d0)
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE ADDBIO2                        *
*                            *****************                         *
* Add on biomass to bio array.                                         *
*----------------------------------------------------------------------*
      SUBROUTINE ADDBIO2(nft,ftmor,npp,nps,npr,bio)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 npp(maxnft),nps(maxnft),npr(maxnft),bio(maxage,2,maxnft)
      REAL*8 npps,nppr
      INTEGER nft,ftmor(maxnft),ft,age

      DO ft=1,nft
        npps = nps(ft)*npp(ft)/100.0d0
        nppr = npr(ft)*npp(ft)/100.0d0
        DO age=1,ftmor(ft)
          bio(age,1,ft) = bio(age,1,ft) + npps
          bio(age,2,ft) = bio(age,2,ft) + nppr
        ENDDO
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE PERC2                          *
*                            ***************                           *
*----------------------------------------------------------------------*
      SUBROUTINE PERC2(fttags,dof,tmin,ftprop,nft)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 dof(maxnft),tmin,ftprop(maxnft),dec1
      REAL*8 evbl,evnl,dec,dcnl,drdc,triscale2,sum
      REAL*8 t1eb,t2eb,t3eb,d1eb,d2eb,d3eb
      REAL*8 t1en,t2en,t3en,d1en,d2en,d3en
      REAL*8 t1d,t2d,t3d,d1d,d2d,d3d
      REAL*8 t1dd,t2dd,t3dd,d1dd,d2dd,d3dd
      REAL*8 t1dn,t2dn,t3dn,d1dn,d2dn,d3dn
      REAL*8 t1d1,t2d1,t3d1,d1d1,d2d1,d3d1,ndof
      INTEGER i,ft,ntags,nft
      CHARACTER fttags(100)*1000,st1*1000


      t1eb = 100.0d0   ! 0
      t2eb = -10.0d0   ! 1
      t3eb = -25.0d0   ! 0
      d1eb =  -1.0d0   ! 0
      d2eb =   0.0d0   ! 1
      d3eb = 150.0d0   ! 0

      t1en = -20.0d0   ! 0
      t2en = -40.0d0   ! 1
      t3en = -60.0d0   ! 0
      d1en = 100.0d0   ! 0
      d2en = 180.0d0   ! 1
      d3en = 260.0d0   ! 0

      t1d =  10.0d0    ! 0
      t2d = -15.0d0    ! 1  broadleaf
      t3d = -50.0d0    ! 0
      d1d =-999.0d0    ! 0
      d2d = 100.0d0    ! 1
      d3d = 140.0d0    ! 0

      t1d1 = 100.0d0   ! 0
      t2d1 =  20.0d0   ! 1  broadleaf
      t3d1 =  10.0d0   ! 0
      d1d1 = 150.0d0   ! 0
      d2d1 = 200.0d0   ! 1
      d3d1 = 400.0d0   ! 0

      t1dd =  12.0d0   ! 0
      t2dd = -15.0d0   ! 1  broadleaf
      t3dd = -30.0d0   ! 0
      d1dd =  75.0d0   ! 0
      d2dd = 200.0d0   ! 1
      d3dd = 400.0d0   ! 0

      t1dn = -60.0d0   ! 0
      t2dn = -70.0d0   ! 1
      t3dn = -80.0d0   ! 0
      d1dn = 180.0d0   ! 0
      d2dn = 240.0d0   ! 1
      d3dn = 300.0d0   ! 0

      nft = 1
      st1='Dc_Nl'
*      ndof = dof(ntags(fttags,st1,nft))
      ndof = 100.0d0
      evbl = triscale2(tmin,ndof,t1eb,t2eb,t3eb,d1eb,d2eb,d3eb)
      evnl = triscale2(tmin,ndof,t1en,t2en,t3en,d1en,d2en,d3en)
      dec  = triscale2(tmin,ndof,t1d ,t2d ,t3d ,d1d ,d2d ,d3d )
      drdc = triscale2(tmin,ndof,t1dd,t2dd,t3dd,d1dd,d2dd,d3dd)
      dcnl = triscale2(tmin,ndof,t1dn,t2dn,t3dn,d1dn,d2dn,d3dn)
      dec1 = triscale2(tmin,ndof,t1d1,t2d1,t3d1,d1d1,d2d1,d3d1)

*      IF (tex-1.LT.0.01d0)  dcnl = 0.0d0

      sum = evbl + evnl + dec + drdc + dcnl + dec1
      IF (sum.GT.0.001d0) THEN
        ftprop(1) = 0.0d0
        ftprop(2) = 0.0d0
        ftprop(3) = 0.0d0
        st1='Ev_Bl'
        ftprop(ntags(fttags,st1,nft)) = 100.0d0*evbl/sum
        st1='Ev_Nl'
        ftprop(ntags(fttags,st1,nft)) = 100.0d0*evnl/sum
        st1='Dc_Bl'
        ftprop(ntags(fttags,st1,nft)) = 100.0d0*(dec + drdc + dec1)/sum
        st1='Dc_Nl'
        ftprop(ntags(fttags,st1,nft)) = 100.0d0*dcnl/sum
      ELSE
        DO i=1,nft
          ftprop(i) = 0.0d0
        ENDDO
        ftprop(2) = 100.0d0
      ENDIF

      ftprop(1) = 100.0d0
      DO ft=2,nft
	  ftprop(1) = ftprop(1) - ftprop(ft)
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE MKLIT2                         *
*                            ****************                          *
* Make litter.                                                         *
*----------------------------------------------------------------------*
      SUBROUTINE MKLIT2(nft,ftgr0,ftmor,ftmat,cov,bio,slc,rlc,sln,rln,
     &npp,nps,npr,debug)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftmat(maxnft),cov(maxage,maxnft),bio(maxage,2,maxnft)
      REAL*8 slc(maxnft),rlc(maxnft),sln(maxnft),ftgr0(maxnft)
      REAL*8 rln(maxnft),npp(maxnft),nps(maxnft),npr(maxnft),npps,nppr
      REAL*8 sl,rl,stlit2
      INTEGER nft,ftmor(maxnft),ft,age
      LOGICAL debug

      DO ft=1,nft
        IF(debug .EQV. .TRUE.) THEN
           print*,'IN MKLIT:
     &ft,npp,nps,npr,cov(1,),cov(2,),cov(3,),ftgr0,ftmor',
     &ft,npp(ft),nps(ft),npr(ft),cov(1,ft),cov(2,ft),cov(3,ft),
     &ftgr0(ft),ftmor(ft)
        ENDIF
        IF (ftgr0(ft).GT.0.0d0) THEN
          npps = nps(ft)*npp(ft)/100.0d0
          nppr = npr(ft)*npp(ft)/100.0d0
          DO age=2,ftmor(ft)
            sl = stlit2(age,ftmat(ft))
            rl = 0.5d0
            slc(ft) = slc(ft) + sl*npps*cov(age,ft)
            rlc(ft) = rlc(ft) + rl*nppr*cov(age,ft)
            bio(age,1,ft) = bio(age,1,ft) - sl*npps
            bio(age,2,ft) = bio(age,2,ft) - rl*nppr
          ENDDO
        ENDIF
      ENDDO

      DO ft=1,nft
        sln(ft) = 0.0d0
        rln(ft) = 0.0d0
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
      FUNCTION stlit2(age,mat)
*----------------------------------------------------------------------*
      REAL*8 stlit2,mat,temp
      INTEGER age

      stlit2 = 0.9d0*real(age)/mat + 0.1d0
      IF (stlit2.GT.1.0)  stlit2 = 1.0d0

      IF (age.LT.mat) then
        stlit2 = 0.0d0
      ELSE
        stlit2 = 1.0d0
      ENDIF

      temp = age/mat
      IF (temp.GT.1.0d0) temp = 1.0d0

      stlit2 = temp**0.5d0
      stlit2 = stlit2*0.1d0
*      stlit2 = 0.0d0

      RETURN
      END


*----------------------------------------------------------------------*
*                            SUBROUTINE VEGMAT2                        *
*                            *****************                         *
* Compute the years to reach maturity for each of the fts.             *
*----------------------------------------------------------------------*
      SUBROUTINE VEGMAT2(nft,ftmor,ftwd,ftxyl,ftpd,evp,lai,npp,nps,
     &ftmat)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 ftwd(maxnft),ftxyl(maxnft),ftpd(maxnft),lai(maxnft)
      REAL*8 npp(maxnft),nps(maxnft),evp(maxnft),ftmat(maxnft),emxv,wd
      REAL*8 xyl,pd,lv,fs,shv,pi,lmvt,nvt,hwvt,dvt,massvt,minvt,ppvt
      INTEGER ft,nft,ftmor(maxnft)

      pi = 3.14159d0
      DO ft=2,nft
	  IF (npp(ft).GT.1.0e-6) THEN
          emxv = evp(ft)/3600.0d0/1000.0d0
          lv = lai(ft)
          fs = nps(ft)/100.0d0
          shv = fs*npp(ft)/100.0d0

          lmvt = lv*1.3d0/(lv + 3.0d0)
          nvt = 1.0d0 + 26.0d0*exp(-0.9d0*lv)
          nvt = 1.0d0

          wd = ftwd(ft)
          xyl = ftxyl(ft)
          pd = ftpd(ft)

          hwvt = pd*lmvt*sqrt(xyl*shv/(emxv*wd*10000.0d0))
     &/nvt
          dvt = 0.0028d0*hwvt**1.5d0
          massvt = pi*(dvt/2.0d0)**2.0d0*hwvt*wd*nvt
          minvt = pi*((dvt + 0.001d0/nvt)/2.0d0)**2.0d0*hwvt*wd*nvt
     &- massvt
          ppvt = shv/minvt
          ftmat(ft) = massvt*ppvt/shv
          ftmat(ft) = 1.0d0
	  ELSE
	    ftmat(ft) = 0.0d0
        ENDIF

        IF (ftmor(ft).LT.ftmat(ft)) ftmat(ft) = ftmor(ft)
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE FIRE2                          *
*                            ***************                           *
* Compute the fire return interval 'fri' and the probability of a fire *
* in the current year 'fprob'.                                         *
*----------------------------------------------------------------------*
      SUBROUTINE FIRE2(prc,tmp,fri,fprob,burn)
*----------------------------------------------------------------------*
      REAL*8 fri,fprob,tmp(12),prc(12),prct(12),prco,lim1,lim2,maxfri
      REAL*8 pow,tlim1,tlim2,totp,indexx,tadj,weight
      INTEGER k,no,ind,minx,i
      LOGICAL burn
      INCLUDE 'param.inc'

*----------------------------------------------------------------------*
* Fire model parameters.                                               *
*----------------------------------------------------------------------*
      no = 3
      lim1 = 150.0d0
      lim2 = 50.0d0
      tlim1 =  0.0d0
      tlim2 = -5.0d0
      weight = 0.5d0
      maxfri = 800.0d0
      pow = 3.0d0

*----------------------------------------------------------------------*
* Adjust for temperature.                                              *
*----------------------------------------------------------------------*
      DO i=1,12
        IF (tmp(i).GT.tlim1) THEN
          tadj = 0.0d0
        ELSEIF (tmp(i).LT.tlim2) THEN
          tadj = lim1
        ELSE
          tadj = (tlim1 - tmp(i))/(tlim1 - tlim2)*lim1
        ENDIF
        prct(i) = min(lim1,prc(i) + tadj)
      ENDDO

*----------------------------------------------------------------------*
* Calculate yearly component of index.                                 *
*----------------------------------------------------------------------*
      totp = 0.0d0
      DO k=1,12
        totp = totp + prct(k)/lim1/12.0d0
      ENDDO

*----------------------------------------------------------------------*
* Calculate monthly component of index.                                *
*----------------------------------------------------------------------*
      prco = 0.0d0
      DO k=1,no
        ind = minx(prct)
        prco = prco + min(lim2,prct(ind))/real(no)/lim2
        prct(ind) = lim2
      ENDDO

*----------------------------------------------------------------------*
* Compute fire return interval, and convert to probability.            *
*----------------------------------------------------------------------*
      indexx = weight*(prco) + (1.0d0 - weight)*totp

      fri = indexx**pow*maxfri

      IF (fri.LT.2.0d0)  fri = 2.0d0

      IF (burn) THEN 
        fprob = 1.0d0 
      ELSE 
        fprob = (1.0d0 - exp(-1.0d0/fri))*p_fprob
      ENDIF 


      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE NEWGROWTH2                     *
*                            ********************                      *
* Compute newgrowth and alter cover array accordingly.                 *
*----------------------------------------------------------------------*
      SUBROUTINE NEWGROWTH2(nft,ftmor,cov,ppm,bio,bioleaf,nppstore,hgt,
     &fprob,npp,nps,tot_ngcov,ngcov,slc,rlc,fireres,firec,harvest,
     &leafdp,flulccc,atharvest,n_at,aggmap_SDGVM_to_aggHyde,NS,yield,
     &ft2frac,sum_cov,debug)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 bio(maxage,2,maxnft),cov(maxage,maxnft),ppm(maxage,maxnft)
      REAL*8 hgt(maxage,maxnft),fprob,npp(maxnft),nppstore(maxnft)
      REAL*8 nps(maxnft),ngcov(maxnft),tot_ngcov
      REAL*8 slc(maxnft),rlc(maxnft),bioleaf(maxnft)
      REAL*8 tmor,tmor0,npp0,firec,xfprob,leafdp(3600,maxnft)
      REAL*8 flulccc(maxnft)
      REAL*8 atharvest(n_at),loss_frac,remain_frac
      REAL*8 fprbtm,yield(maxnft),ft2frac(maxnft),dflulccc
      REAL*8 sum_cov(maxnft)
      INTEGER at,aggmap_SDGVM_to_aggHyde(NS)
      INTEGER nft,ftmor(maxnft),ft,age,fireres,n_at,NS
      LOGICAL harvest
      LOGICAL debug !used to print out more debugging info

      xfprob = fprob

      firec = 0.0d0
*----------------------------------------------------------------------*
* Take away veg that has died of old age ie > than ftmor(ft), and      *
* shift cover array on one year.                                       *
*----------------------------------------------------------------------*
      !tot_ngcov = 0.0d0  !This is currently accumulated over function
                          !calls 
      DO ft=2,nft
        tot_ngcov = tot_ngcov + cov(ftmor(ft),ft)
        ngcov(ft) = ngcov(ft) + cov(ftmor(ft),ft)
      
        slc(ft) = slc(ft)  + (bio(ftmor(ft),1,ft) + bioleaf(ft) +
     &nppstore(ft))*cov(ftmor(ft),ft)
        rlc(ft) = rlc(ft)  + bio(ftmor(ft),2,ft)*cov(ftmor(ft),ft)
        IF (debug .EQV. .TRUE.) THEN
                PRINT '(A I3 F9.6 F9.6 I6 F9.6)','GG1B',
     &              ft, tot_ngcov,ngcov(ft),ftmor(ft),cov(ftmor(ft),ft)
        ENDIF
      ENDDO

      IF(debug .EQV. .TRUE.) THEN
        DO ft=2,nft
           print*,'NEWGROWTH STEP 1: ft,slc,rlc: ',ft,slc(ft),rlc(ft)
        ENDDO
      ENDIF

      CALL SHIFT2(ftmor,cov,ppm,bio,hgt,1,nft)

*----------------------------------------------------------------------*
* Take away veg that is burnt or has died through a hard year ie small *
* LAI.                                                                 *
*----------------------------------------------------------------------*
      DO ft=2,nft
*----------------------------------------------------------------------*
* Compute harvest losses for each aggregated type (at)                 *
*----------------------------------------------------------------------*
        at = aggmap_SDGVM_to_aggHyde(ft)
        loss_frac = 0.0d0
        remain_frac = 1.0d0
        IF (at.NE.0) THEN
          IF ((atharvest(at).GT.0.0d0).AND.
     &         ((at.EQ.1).OR.(at.EQ.3))) THEN !only for primf&secdf
            loss_frac  = ft2frac(ft)*atharvest(at)*1.0d-2
            remain_frac = 1.0d0 - loss_frac  
          ENDIF
        ENDIF
*----------------------------------------------------------------------*
* 'tmor' is the mortality rate of the forest based on 'npp'.           *
*----------------------------------------------------------------------*
        npp0 = 0.2d0
        tmor0 = 6.0d0
        IF (npp(ft)/100.0d0.LT.0.1d0) THEN
          tmor = 1.0d0
        ELSEIF (npp(ft)/100.0d0.LT.3.0d0) THEN
          tmor = (0.04d0 - npp0)/3.0d0*(npp(ft)/100.0d0) + npp0
        ELSEIF (npp(ft)/100.0d0.LT.tmor0) THEN
          tmor = 0.04d0/(3.0d0 - tmor0)*(npp(ft)/100.d0) + 
     &3.0d0*0.04d0/(tmor0 - 3.0d0) + 0.04d0
        ELSE
          tmor = 0.0d0
        ENDIF
        IF (tmor.LT.0.01d0)  tmor = 0.01d0

*        IF (npp(ft)*nps(ft)/100.0d0.GT.10.0d0) THEN
        tmor = 0.002d0
*        ELSE
*          tmor = 1.0d0 - npp(ft)*nps(ft)/1000.0d0
*          if (ft.eq.6) print*,'hello ',tmor
*        ENDIF

*----------------------------------------------------------------------*
* kill off trees not able to sustain any growth rate                   *
*----------------------------------------------------------------------*
        IF ((.NOT.(npp(ft)*nps(ft).GT.0.0d0)).OR.
     &(.NOT.(nppstore(ft)).GT.0.0d0)) THEN
          tmor = 1.0d0
*          WRITE(*,*) 'Death ft npp nps store: ',ft,npp(ft),nps(ft),
*     &nppstore(ft)
        ENDIF
*----------------------------------------------------------------------*
        !tmor = 0.0d0 !PCM turn tmor mortality off 

        sum_cov(ft) = 0.d0 
        DO age=2,ftmor(ft) 
          IF (age.LT.fireres) THEN
            fprob = xfprob
          ELSE
            fprob = 0.0d0
          ENDIF
          IF (fireres.LT.0) fprob = real(-fireres)/1000.0d0
!          fprob = 0.0d0 !PCM turn fire off

          fprbtm = fprob-fprob*tmor+tmor
          
          !calculate litter from cover loss  
          tot_ngcov = tot_ngcov + cov(age,ft)*fprbtm
          ngcov(ft) = ngcov(ft) + cov(age,ft)*fprbtm
          
          IF (debug .EQV. .TRUE.) THEN
                PRINT '(A I3 F9.6 F9.6 F9.6 F9.6 F9.6 F9.6)','GG1C',
     &         ft,tot_ngcov,ngcov(ft),
     &         cov(age,ft),fprbtm,
     &         fprob,tmor
          ENDIF
          slc(ft) = slc(ft) + 
     &( bio(age,1,ft) + bioleaf(ft) + nppstore(ft) ) * 
     &(tmor - 0.2d0*fprob*tmor + 0.2d0*fprob) * cov(age,ft)

          rlc(ft) = rlc(ft) + bio(age,2,ft)*fprbtm*cov(age,ft)

          firec   = firec + 
     &( bio(age,1,ft) + bioleaf(ft) + nppstore(ft) ) *
     &(0.8d0*fprob - 0.8d0*fprob*tmor) * cov(age,ft)

          !update cover array
          cov(age,ft) = cov(age,ft)*(1.0d0 - fprob)*(1.0d0 - tmor)
          !update sum_cov for next time
          sum_cov(ft) = sum_cov(ft) + cov(age,ft)
          
        ENDDO

        !split the loop into two, since we have a new cov array with a
        !different sum_cov
        DO age=2,ftmor(ft) 
          ! calculate litter loss and biomass loss from harvest
          ! for a harvest (copice style) the cov array remains unchanged
          IF (harvest) THEN
            !this can act as a coppice type harvest or a fire that leaves the root mass intact and cover intact
            if(age.EQ.2) print*, 'harvest'
            !add harvested/removed leaf biomass to surface soil litter
            slc(ft) = slc(ft) + bioleaf(ft)*cov(age,ft)
            !add harvested/removed wood biomass (including 50% nppstore) to firec losses
            IF(sum_cov(ft).GT.0.0d0)THEN
              flulccc(ft) = flulccc(ft) +
     &bio(age,1,ft)*cov(age,ft)/sum_cov(ft)
            ELSE
              flulccc(ft) = 0.0d0
            ENDIF 
            firec   = firec   + nppstore(ft) * 0.50d0 * cov(age,ft) 
            bio(age,1,ft) =  0.0d0
          ENDIF
          
          ! calculate litter loss and biomass loss from harvest
          ! the cover array has already been adjusted in the COVER
          ! routine for harvest
          IF ((loss_frac.GT.0.0d0).AND.((at.EQ.1).OR.(at.EQ.3))) THEN
            !add harvested/removed leaf biomass to surface soil litter
            slc(ft) = slc(ft) + bioleaf(ft)*cov(age,ft)*loss_frac
            !add harvested/removed wood biomass (including 50% nppstore) to lulccc losses
            IF(sum_cov(ft).GT.0.0d0)THEN
              dflulccc = (bio(age,1,ft) + nppstore(ft) * 0.50d0) *
     &cov(age,ft)*loss_frac/sum_cov(ft)
            ELSE
              dflulccc = 0.0d0
            ENDIF 
            !flulccc(ft) = flulccc(ft) + dflulccc !no double counting
            !add harvested/removed wood biomass (including 50% nppstore) to yield
            yield(ft) = yield(ft) + dflulccc
            bio(age,1,ft) = bio(age,1,ft)*remain_frac
          ENDIF

        ENDDO

        IF(debug) THEN
           PRINT *,
     &'harvest',ft,at,loss_frac,dflulccc,flulccc(ft),yield(ft)
        ENDIF

        IF (harvest) THEN
          !this can act as a coppice type harvest or a fire that leaves the root mass intact and cover intact
          !remove leaf mass - this leaves root mass and 50% nppstore untouched 
          bioleaf(ft)   = 0.0d0
          leafdp(:,ft)  = 0.0d0 
          nppstore(ft)  = nppstore(ft) * 0.50d0
        ENDIF

        IF ((loss_frac.GT.0.0d0).AND.((at.EQ.1).OR.(at.EQ.3))) THEN !only for primf&secdf
          !this can act as a coppice type harvest or a fire that leaves the root mass intact and cover intact
          !remove leaf mass - this leaves root mass and 50% nppstore untouched 
          bioleaf(ft)   = bioleaf(ft)*remain_frac 
          leafdp(:,ft)  = leafdp(:,ft)*remain_frac  
          nppstore(ft)  = nppstore(ft) * (1.0d0 - loss_frac*0.50d0)
        ENDIF

      ENDDO

      fprob = xfprob
      IF(debug .EQV. .TRUE.) THEN
        DO ft=2,nft
           print*,'END NEWGROWTH: ft,slc,rlc: ',ft,slc(ft),rlc(ft)
        ENDDO
      ENDIF

      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE LULCC_LOSS2                    *
*                            ********************                      *
* Compute newgrowth area and alter cover array according to land-use   *
* and land-cover change database                                       *
*----------------------------------------------------------------------*
      SUBROUTINE LULCC_LOSS2(nft,ftmor,cov,ppm,bio,bioleaf,nppstore,hgt,
     &loss,loss_nowoodh,npp,nps,slc,rlc,fireres,flulccc,harvest,leafdp,
     &sum_cov,ft,debug)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 bio(maxage,2,maxnft),cov(maxage,maxnft),ppm(maxage,maxnft)
      REAL*8 hgt(maxage,maxnft),npp(maxnft),nppstore(maxnft)
      REAL*8 loss,loss_nowoodh
      REAL*8 nps(maxnft),dflulccc,sum_cov(maxnft),orig_sum_cov
      REAL*8 slc(maxnft),rlc(maxnft),bioleaf(maxnft)
      REAL*8 tmor,tmor0,npp0,flulccc(maxnft),xfprob,leafdp(3600,maxnft)
      INTEGER nft,ftmor(maxnft),ft,age,fireres
      LOGICAL harvest
      LOGICAL debug !used to print out more debugging info


*----------------------------------------------------------------------*
* kill off pfts that have lost cover according to the landuse database * 
*----------------------------------------------------------------------*

        IF(debug .EQV. .TRUE.) THEN
          PRINT '(A)','GG7 LULCC_LOSS'
        ENDIF
        orig_sum_cov = sum_cov(ft)
        sum_cov(ft) = 0.0d0 
        DO age=1,ftmor(ft)
           rlc(ft) = rlc(ft) + bio(age,2,ft) * 
     &loss * cov(age,ft)

           IF(orig_sum_cov.GT.0.0d0) THEN
             dflulccc = 
     &( bio(age,1,ft) + bioleaf(ft) + nppstore(ft) ) *    
     &loss * cov(age,ft)/orig_sum_cov
           ELSE
             dflulccc = 0.0d0
           ENDIF

           flulccc(ft) = flulccc(ft) + dflulccc 

           !update cover array
           cov(age,ft) = cov(age,ft)*( 1.0d0 - loss )
           !update sum_cov for next time
           sum_cov(ft) = sum_cov(ft) + cov(age,ft)
           IF(debug .EQV. .TRUE.) THEN
             PRINT '(2I5,7F12.7)',age,ft,loss,loss_nowoodh,
     &flulccc(ft),bio(age,1,ft),
     &bioleaf(ft),nppstore(ft),cov(age,ft)
           ENDIF
        ENDDO
        

      RETURN
      END

*----------------------------------------------------------------------*
*                            SUBROUTINE SHIFT2                         *
*                            ****************                          *
* Shift shifts elements of an array to the right and leaves the first  *
* element equal to zero.                                               *
*----------------------------------------------------------------------*
      SUBROUTINE SHIFT2(ftmor,cov,ppm,bio,hgt,ft1,ft2)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 cov(maxage,maxnft),ppm(maxage,maxnft),bio(maxage,2,maxnft)
      REAL*8 hgt(maxage,maxnft)
      INTEGER ftmor(maxnft),ft1,ft2,ft,age,i

      DO ft=ft1,ft2
        DO age=2,ftmor(ft)
          cov(ftmor(ft)-age+2,ft) = cov(ftmor(ft)-age+1,ft)
          ppm(ftmor(ft)-age+2,ft) = ppm(ftmor(ft)-age+1,ft)
          hgt(ftmor(ft)-age+2,ft) = hgt(ftmor(ft)-age+1,ft)
          DO i=1,2
            bio(ftmor(ft)-age+2,i,ft) = bio(ftmor(ft)-age+1,i,ft)
          ENDDO
        ENDDO
        cov(1,ft) = 0.0d0
        ppm(1,ft) = 0.0d0
        hgt(1,ft) = 0.0d0
        bio(1,1,ft) = 0.0d0
        bio(1,2,ft) = 0.0d0
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          NATURAL_VEG2                                *
*                          ***********                                 *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE NATURAL_VEG2(tmp,prc,ftprop,nat_map)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 tmp(12,31),prc(12,31),ftprop(maxnft),min_tmp,av_mnth(12)
      REAL*8 x(6,6),evap,tran,bucket,dep,prc_ind
      REAL*8 triscale,sum,prop(6)
      INTEGER day,mnth,ft,nat_map(8),vdaysoff(5),i,year,daysoff,minx
      LOGICAL leaves
      DATA (x(1,i),i=1,6)/ -15.0, -10.0, 100.0,   -1.0,   0.0, 150.0/!EB
      DATA (x(2,i),i=1,6)/ -70.0, -45.0, -20.0,  -10.0,   0.0, 350.0/!EN
      DATA (x(3,i),i=1,6)/ -50.0, -15.0,  10.0, -999.0, 100.0, 140.0/!DB
      DATA (x(4,i),i=1,6)/  10.0,  20.0, 100.0,  150.0, 200.0, 400.0/!DB
      DATA (x(5,i),i=1,6)/ -30.0, -15.0,  12.0,  75.0,  200.0, 400.0/!DB
      DATA (x(6,i),i=1,6)/ -80.0, -70.0, -60.0, 200.0,  280.0, 300.0/!DN

*----------------------------------------------------------------------*
* Find average monthly min, and convert to an absolute min.            *
*----------------------------------------------------------------------*
      min_tmp = 1000.0d0
      DO mnth=1,12
        av_mnth(mnth) = 0.0d0
        DO day=1,30
          IF (min_tmp.GT.tmp(mnth,day)) min_tmp = tmp(mnth,day)
          av_mnth(mnth) = av_mnth(mnth) + tmp(mnth,day)     
        ENDDO
        av_mnth(mnth) = av_mnth(mnth)/30.0d0     
      ENDDO
      min_tmp = av_mnth(minx(av_mnth))*1.29772d0 - 19.5362d0

*----------------------------------------------------------------------*
* Find days off index. Based on an imaginary bucket.                   *
*----------------------------------------------------------------------*
      dep = 100.0d0
      evap = 0.5d0
      tran = 0.5d0

      DO i=1,5
        vdaysoff(i) = 0
      ENDDO
      leaves = .false.

      bucket = 0.0d0

      DO year=1,10
        daysoff = 0
        DO mnth=1,12
          DO day=1,30
            bucket = bucket + prc(mnth,day)
            IF (leaves) THEN
              bucket = bucket - evap - tran
            ELSE
              bucket = bucket - evap
            ENDIF
            IF (bucket.LT.0.0) THEN
              leaves = .false.
              bucket = 0.0
            ENDIF 
            IF (bucket.gt.dep/2.0) leaves = .true.
            IF (.NOT.(leaves)) daysoff = daysoff + 1
            IF (bucket.GT.dep) bucket = dep
          ENDDO
        ENDDO
        DO i=1,4
          vdaysoff(6-i) = vdaysoff(5-i)
        ENDDO
        vdaysoff(1) = daysoff
      ENDDO
      prc_ind = real(vdaysoff(1)+vdaysoff(2)+vdaysoff(3)+vdaysoff(4))
     &/4.0d0

*----------------------------------------------------------------------*
* Compute fractions of natural vegetation appropriate to seed new      *
* ground.                                                              *
*----------------------------------------------------------------------*
      DO i=1,6
        prop(i) =  triscale(min_tmp,x(i,1),x(i,2),x(i,3))*
     &             triscale(prc_ind,x(i,4),x(i,5),x(i,6))
      ENDDO

*----------------------------------------------------------------------*
* Map the natural functional types to the list given in the input file.*
*----------------------------------------------------------------------*
      ftprop(nat_map(1)) = 0.0d0
      ftprop(nat_map(2)) = 0.0d0
      ftprop(nat_map(3)) = 0.0d0
      ftprop(nat_map(4)) = 0.0d0
      ftprop(nat_map(5)) = prop(1)
      ftprop(nat_map(6)) = prop(2)
      ftprop(nat_map(7)) = prop(3) + prop(4) + prop(5)
      ftprop(nat_map(8)) = prop(6)

*----------------------------------------------------------------------*
* Normalise the fractions.                                             *
*----------------------------------------------------------------------*
      sum = 0.0d0
      DO ft=1,8
        sum = sum + ftprop(nat_map(ft))
      ENDDO
      IF (sum.GT.0.0d0) THEN
        DO ft=1,8
          ftprop(nat_map(ft)) = 100.0d0*ftprop(nat_map(ft))/sum
        ENDDO
      ELSE
        ftprop(nat_map(3)) = 100.0d0
      ENDIF

*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*      DO k=-80,100
*        DO j=0,360
*        min_tmp =  (real(k) + 19.5632d0)/1.29772d0
*        prc_ind = real(j)
*      DO i=1,6
*        ftprop(i) =  triscale(min_tmp,x(i,1),x(i,2),x(i,3))*
*     &               triscale(prc_ind,x(i,4),x(i,5),x(i,6))
*      ENDDO
*      sum=ftprop(1)+ftprop(2)+ftprop(3)+ftprop(4)+ftprop(5)+ftprop(6)
*        IF (sum.LT.1.0e-6) print'(2i4)',k,j
*        ENDDO
*      ENDDO
*      STOP
*----------------------------------------------------------------------*


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                     SUBROUTINE CALC_FT2FRAC2                         *
*                     ***********************                          *
*----------------------------------------------------------------------*
      SUBROUTINE CALC_FT2FRAC2(ftprop_init,aggmap_SDGVM_to_aggHyde,nft,
     &ft2frac,debug)
      INCLUDE 'array_dims.inc'
      INTEGER, PARAMETER :: NS = 16 !number of SDGVM functional types
      REAL*8 ftprop_init(maxnft)
      INTEGER ft2,at,at2,at3,ft3,nft,aggmap_SDGVM_to_aggHyde(NS)
      REAL*8 ft2frac(maxnft),at2prop
      LOGICAL debug

!compute fractions of cover (ft2frac(ft2)) of each ft2 in each aggregated class at2 
      DO ft2=2,nft
        at2 = aggmap_SDGVM_to_aggHyde(ft2)
        at2prop = 0.0d0
        ft2frac(ft2) = 0.0d0 
        DO ft3=2,nft
          at3 = aggmap_SDGVM_to_aggHyde(ft3)
          if(at2.eq.at3) then
            at2prop = at2prop + ftprop_init(ft3)  
          endif
        ENDDO
        if(at2prop.GT.0.0d0) then
          ft2frac(ft2) = ftprop_init(ft2)/at2prop
        endif
        !IF (debug .EQV. .TRUE. ) THEN
        !     PRINT '(A I2 I2 F11.6)', 'GH0',at2,ft2,ft2frac(ft2)
        !ENDIF
      ENDDO

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                     SUBROUTINE CALC_CORRCT_FT2                       *
*                     ***********************                          *
*----------------------------------------------------------------------*
      SUBROUTINE CALC_CORRCT_FT2(ft,lat,corrct_ft)
      INTEGER ft
      REAL*8 lat
      LOGICAL corrct_ft

      corrct_ft = .FALSE.
      IF(ABS(lat).LE.24.5)THEN
       IF((ft.EQ.4).OR.(ft.EQ.8).OR.(ft.EQ.9).OR.(ft.EQ.13))THEN !Tropical
           corrct_ft = .TRUE.
       ELSE IF((ft.EQ.2).OR.(ft.EQ.5).OR.(ft.EQ.6))THEN !URBAN or crop
           corrct_ft = .TRUE.
       ELSE
           corrct_ft = .FALSE.
      ENDIF
       ELSE
         IF((ft.EQ.3).OR.(ft.EQ.7))THEN !C3 grasses
           corrct_ft = .TRUE.
         ELSE IF((ft.GT.9).AND.(ft.LE.12))THEN !Temperate trees (primary)
           corrct_ft = .TRUE.
         ELSE IF((ft.GT.13).AND.(ft.LE.16))THEN !Temperate trees (secondary)
           corrct_ft = .TRUE.
         ELSE IF((ft.EQ.2).OR.(ft.EQ.5).OR.(ft.EQ.6))THEN !URBAN or crop
           corrct_ft = .TRUE.
         ELSE
           corrct_ft = .FALSE.
         ENDIF
       ENDIF
      
      RETURN
      END

