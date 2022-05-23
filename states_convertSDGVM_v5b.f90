! This is a test program written by Patrick McGuire that reads the states.nc file from the Hyde LUC database,
!      and then converts it to the SDGVM format with the SDGVM land cover
!      classes. This program also uses the ESA CCI data as input.
! The Hyde LUC database is described here:  
!   LUH2-GCB2021: Land-Use Harmonization Data Set for GCB2021: years 850-2021
!   title: UofMD LUH2-GCB2021 dataset prepared for Global Carbon Budget
!    references: Hurtt et al. 2020, Chini et al. 2021
! This is derived from a simple UCAR example which reads a NETCDF file : simple_xy_rd.f90
! This is also derived from another test program, states_rd.f90, written by Patrick McGuire.
! And this is also derived from Anthony Walker's program process_LULCC_4SDGVM_v1a1.R, as adapted for JASMIN by Patrick McGuire.
      
!     Note from Patrick McGuire (p.mcguire@reading.ac.uk) 5/15/2022
!     to compile this on the CEDA JASMIN supercomputer:
!       gfortran states_convertSDGVM_v5b.f90 -o states_convertSDGVM_v5b `nf-config --fflags --flibs`

      PROGRAM states_convertSDGVM_v5b
      USE netcdf
      IMPLICIT NONE

      ! half-resolution version computed with: module load jasppy ! on JASMIN
      !                                        cdo gridboxmean,2,2 transitions.nc transitions2b.nc
      !                                        cdo gridboxmean,2,2 states.nc states2b.nc
      CHARACTER (LEN = *), PARAMETER :: fname = "states2b.nc" !half-resolution version of states.nc
      CHARACTER (LEN = *), PARAMETER :: fname_t = "transitions2b.nc" !half-resolution version of transitions.nc
      CHARACTER (LEN = *), PARAMETER :: wdh = &
       "/gws/nopw/j04/nexcs/pmcguire/TRENDYv10/db/LUH2_GCB_2021/"
      CHARACTER (LEN = *), PARAMETER :: wdg = "/gws/nopw/j04/nexcs/pmcguire/sdgvmD/data/land_use/global/ESACCILCP2014/30min/"
      CHARACTER (LEN = *), PARAMETER :: fname_s = "cont_lu"
      CHARACTER (LEN = *), PARAMETER :: esadate = "2009"
      !CHARACTER (LEN = *), PARAMETER :: print_type='unagg'
      !CHARACTER (LEN = *), PARAMETER :: print_type='agg'
      CHARACTER (LEN = *), PARAMETER :: print_type='sdgvm'

      !    DIMENSIONS(sizes): time(1172), lon(1440), lat(720)
      ! INTEGER, PARAMETER :: NX = 1440, NY = 720, NT = 1172, NV = 14, NV2 = 6
      INTEGER, PARAMETER :: NX = 720, NY = 360, NT = 1172, NV = 14, NV2 = 6
      REAL, PARAMETER    :: misval = 1e19
      INTEGER, PARAMETER :: NE = 10, NE2 = 15 
      INTEGER, PARAMETER :: NVT = 118 
      INTEGER, PARAMETER :: DT = 20 !number of years to skip between prints
      !INTEGER, PARAMETER :: ST = 21 !starting year for prints
      INTEGER, PARAMETER :: ST = 1001 !starting year for prints
      REAL :: data_in(NV, NX, NY), dummya(NX, NY)
      REAL :: data_in_old(NV, NX, NY)
      REAL :: data_in_new(NV, NX, NY)
      REAL :: data_in_t(NVT, NX, NY)
      REAL :: data_out(NV2, NX, NY)
      REAL :: data_out_SDGVM(NE, NX, NY)
      LOGICAL :: mask(NX, NY) = .false.
      LOGICAL :: mask_SDGVM(NX, NY) = .false.
      LOGICAL :: esamask(NE2,NX, NY) = .false.
      ! setup ESA arrays
      REAL :: esaarray(NE2,NX,NY) = 0.0
      CHARACTER (LEN = 200) :: fname_s2
      LOGICAL :: compute_next_year = .true.
      REAL, PARAMETER :: NA= -999.0

      ! ESA PFTs
      CHARACTER(LEN=7),PARAMETER :: varname_esa_pfts(NE)=(/'  BARE',' Ev_Bl',' Dc_Bl',' Ev_Nl',' Dc_Nl', &
          ' Shrub','    C3','    C4','C3crop','C4crop'/)

     !  CHARACTER(LEN=7),PARAMETER :: varname(NV)=(/'primf', 'primn', 'secdf', 'secdn', 'urban', &
     !     'c3ann', 'c4ann', 'c3per', 'c4per', 'c3nfx', 'pastr', 'range', &
     !     'secmb', 'secma'/)
     ! using Anthony's order:
      CHARACTER(LEN=7),PARAMETER :: varname(NV)=(/'primf', 'secdf', 'secdn', 'primn', &
          'c3ann', 'c3per', 'c3nfx', 'c4ann', 'c4per', 'pastr','range',               &
          'urban','secmb', 'secma'/)
      CHARACTER(LEN=5) :: from_t, to_t

      CHARACTER(LEN=17),PARAMETER :: varname_t(NVT)=(/ &
           'primf_to_secdn','primf_to_urban', 'primf_to_c3ann','primf_to_c4ann',&
           'primf_to_c3per','primf_to_c4per', 'primf_to_c3nfx','primf_to_pastr',&
           'primf_to_range',&
           'secdf_to_secdn','secdf_to_urban', 'secdf_to_c3ann','secdf_to_c4ann',&
           'secdf_to_c3per','secdf_to_c4per', 'secdf_to_c3nfx','secdf_to_pastr',&
           'secdf_to_range',&
           'secdn_to_secdf','secdn_to_urban','secdn_to_c3ann', 'secdn_to_c4ann',&
           'secdn_to_c3per','secdn_to_c4per','secdn_to_c3nfx', 'secdn_to_pastr',&
           'secdn_to_range',&
           'primn_to_secdf','primn_to_urban','primn_to_c3ann', 'primn_to_c4ann',&
           'primn_to_c3per','primn_to_c4per','primn_to_c3nfx', 'primn_to_pastr',&
           'primn_to_range',&
           'c3ann_to_secdf','c3ann_to_secdn','c3ann_to_urban','c3ann_to_c4ann',&
           'c3ann_to_c3per','c3ann_to_c4per','c3ann_to_c3nfx','c3ann_to_pastr',&
           'c3ann_to_range',&
           'c3per_to_secdf','c3per_to_secdn','c3per_to_urban','c3per_to_c3ann',&
           'c3per_to_c4ann','c3per_to_c4per','c3per_to_c3nfx','c3per_to_pastr',&
           'c3per_to_range',&
           'c3nfx_to_secdf','c3nfx_to_secdn','c3nfx_to_urban','c3nfx_to_c3ann',&
           'c3nfx_to_c4ann','c3nfx_to_c3per','c3nfx_to_c4per','c3nfx_to_pastr',&
           'c3nfx_to_range',&
           'c4ann_to_secdf','c4ann_to_secdn','c4ann_to_urban','c4ann_to_c3ann',&
           'c4ann_to_c3per','c4ann_to_c4per','c4ann_to_c3nfx','c4ann_to_pastr',&
           'c4ann_to_range',&
           'c4per_to_secdf','c4per_to_secdn','c4per_to_urban','c4per_to_c3ann',&
           'c4per_to_c4ann','c4per_to_c3per','c4per_to_c3nfx','c4per_to_pastr',&
           'c4per_to_range',&
           'pastr_to_secdf','pastr_to_secdn','pastr_to_urban','pastr_to_c3ann',&
           'pastr_to_c4ann','pastr_to_c3per','pastr_to_c4per','pastr_to_c3nfx',&
           'pastr_to_range',&
           'range_to_secdf','range_to_secdn','range_to_urban','range_to_c3ann',&
           'range_to_c4ann','range_to_c3per','range_to_c4per','range_to_c3nfx',&
           'range_to_pastr',&
           'urban_to_secdf','urban_to_secdn', 'urban_to_c3ann','urban_to_c4ann',&
           'urban_to_c3per','urban_to_c4per', 'urban_to_c3nfx','urban_to_pastr',&
           'urban_to_range',&
           'primf_harv    ','primn_harv    ',&
           'secmf_harv    ','secyf_harv    ',&
           'secnf_harv    ','primf_bioh    ',&
           'primn_bioh    ','secmf_bioh    ',&
           'secyf_bioh    ','secnf_bioh    ' /)

      ! aggregate Hyde landcover types 
      CHARACTER(LEN=9),PARAMETER :: varname2(NV2)=(/'   forest', ' nforestr', '   c3crop', '   c4crop', &
          '   pastnr', '    urban' /)
      ! primary and secondary forest
      ! primary and secondary non-forest and rangelands
      ! c3 crops
      ! c4 crops
      ! pasture and not rangelands
      ! urban

      ! This will be the netCDF ID for the file and data variable.
      INTEGER :: ncid, varid(NV)
      INTEGER :: ncid_t, varid_t(NVT)

      INTEGER :: num_land

      ! Loop indexes, and error handling.
      INTEGER :: x, y, t, v, v2, i, lon, lat


      ! Open ESA CCILCP 2014 dataset 
      !setwd(wdg)
      DO i=1,NE 
        WRITE (fname_s2, "(A,I0,A)") wdg//fname_s//"-", i,"-"//esadate//".dat" 
        !esav          <- scan(fname_s2)
        !! esamat <- if(i==1) esav else cbind(esamat,esav)
        !esaarray(,,i) <- as.matrix(esav,nrow=lon_ress)
        OPEN(15,FILE=fname_s2,STATUS='old')
        READ(15,*) esaarray(i,:,:)
        CLOSE(15)
      END DO 

      ! convert missing values      
      WHERE (esaarray .EQ. 255) esaarray = NA      
      WHERE (esaarray .NE. NA) esamask = .TRUE.      

      num_land = COUNT( esamask(1,:,:) .EQV. .TRUE.)
      PRINT *, 'ESA num_land=',num_land,'num_tot=',NX*NY 

      ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
      ! the file.
      CALL CHECK( NF90_OPEN(wdh//fname, NF90_NOWRITE, ncid) )

      DO v=1,NV
      ! Get the varid of the data variable, based on its name.
        CALL check( nf90_inq_varid(ncid, varname(v), varid(v)) )
      END DO

      CALL CHECK( NF90_OPEN(wdh//fname_t, NF90_NOWRITE, ncid_t) )

      DO v=1,NVT
      ! Get the varid of the data variable, based on its name.
        CALL CHECK( NF90_INQ_VARID(ncid_t, varname_t(v), varid_t(v)) )
      END DO

      DO t=ST,NT,1 
        data_out(:,:,:) = 0.0
        IF(t > ST) THEN 
          data_in_old(:,:,:)=data_in(:,:,:)
        END IF
      ! Read the data.
        DO v=1,NV
          CALL CHECK( NF90_GET_VAR(ncid, varid(v), data_in(v,:,:), start=[1,1,t], count=[NX,NY,1]) )
      !  print *,'Finished reading data'
          IF( t == ST .AND. v == 1) THEN
            WHERE(data_in(1,:,:) <= 1.00) mask = .TRUE.
            num_land = COUNT( mask .EQV. .TRUE.)
            PRINT *, 'num_land=',num_land,'num_tot=',NX*NY 
          END IF

        ! for now, use regridded states.nc file from CDO, so we skip the next two steps
        ! regrid from 0.25 to 0.5
        !  dummya = reduce_res(data_in(v,:,:))

        ! convert missing values      
        !  where (dummya > misval) mask = .false.      

        END DO


        IF(compute_next_year) THEN
          ! the secondary vegetation doesn't match up unless we
          ! add the transitions to the states for each time step (rather than
          ! accumulating all the transitions in an open loop).
          ! if(t == ST) then 
           data_in_new(:,:,:)=data_in(:,:,:) 
          ! end if

           DO v=1,NVT
            CALL CHECK( NF90_GET_VAR(ncid_t, varid_t(v), data_in_t(v,:,:), start=[1,1,t], count=[NX,NY,1]) )
            from_t = varname_t(v)(1:5)
            to_t   = varname_t(v)(10:14)
            DO v2=1,NV
              IF((v<NVT-5) .AND. (varname(v2) == from_t)) THEN! skip for bioh
               !PRINT *, varname_t(v), varname(v2), from_t, v, v2
               data_in_new(v2,:,:)= data_in_new(v2,:,:) - data_in_t(v,:,:) 
              END IF
              IF((v<NVT-5) .AND. (varname(v2) == to_t)) THEN ! skip for bioh
               !PRINT *, varname_t(v), varname(v2), to_t, v, v2
               data_in_new(v2,:,:)= data_in_new(v2,:,:) + data_in_t(v,:,:) 
              END IF
            END DO
           END DO
        END If

        !aggregate
        DO v=1,NV
        ! all hyde default classes array 
        !dummya1[,,v] <- dummya
          IF(compute_next_year) THEN
            dummya = data_in_new(v,:,:)
          ELSE
            dummya = data_in(v,:,:)
          END IF

        ! aggregate Hyde landcover types 
        ! primary and secondary forest
          IF(v<=2) THEN
            data_out(1,:,:) = data_out(1,:,:) + dummya
        ! primary and secondary non-forest
          ELSE IF(v<=4) THEN
            data_out(2,:,:) = data_out(2,:,:) + dummya
        ! c3 crops
          ELSE IF(v<=7) THEN
            data_out(3,:,:) = data_out(3,:,:) + dummya
        ! c4 crops
          ELSE IF(v<=9) THEN
            data_out(4,:,:) = data_out(4,:,:) + dummya
        ! ! pasture and rangelands
        !   ELSE IF(v<=11) THEN
        !    data_out(5,:,:) = data_out(5,:,:) + dummya
        ! pasture and not rangelands
          ELSE IF(v==10) THEN
            data_out(5,:,:) = data_out(5,:,:) + dummya
        ! rangelands - add to non-forest
          ELSE IF(v==11) THEN
            data_out(2,:,:) = data_out(2,:,:) + dummya
        ! urban
          ELSE IF(v==12) THEN
            data_out(6,:,:) = dummya 
          END IF
        END DO

        ! change to percent      
        data_out = data_out * 100.0      

        esaarray(11:15,:,:) = data_out(1:5,:,:)
     
        ! process, returns an array of PFT, lon, lat  
        !data_out_SDGVM   <- apply(esaarray,c(1,2),join_hyde)
        DO lon=1,NX
         DO lat=1,NY
           data_out_SDGVM(:,lon,lat) = JOIN_HYDE(esaarray(:,lon,lat))
         END DO
        END DO
      
        WHERE(ABS(data_out_SDGVM(1,:,:)) <= 100.00) mask_SDGVM = .TRUE.

        ! convert missing values to SDGVM missing value
        WHERE (data_out_SDGVM == NA) data_out_SDGVM = 255 

        IF( t == ST ) THEN
          num_land = COUNT( mask_SDGVM .EQV. .TRUE.)
          PRINT *, 'SDGVM num_land=',num_land,'num_tot=',NX*NY 
          WRITE(*,FMT='(A5,A2)',ADVANCE='no')'   t','  '
          IF(print_type=='unagg') THEN
            DO v=1,NV
              WRITE(*,FMT='(A9)', ADVANCE='no') varname(v)
            END DO
          ELSE IF(print_type=='agg') THEN
            DO v=1,NV2
               WRITE(*,FMT='(A9)', ADVANCE='no') varname2(v)
            END DO
          ELSE IF(print_type=='sdgvm') THEN
            DO v=1,NE
               WRITE(*,FMT='(A9)', ADVANCE='no') varname_esa_pfts(v)
            END DO
          END IF
          WRITE(*,*) ! Assumes default "ADVANCE='yes'".
        END IF

        IF( MOD(t-ST,DT) == 0 ) THEN 
           IF(compute_next_year) THEN
               WRITE(*,FMT='(I5)', ADVANCE='no') t+849+1 !t=1 is the year 850
           ELSE
               WRITE(*,FMT='(I5)', ADVANCE='no') t+849 !t=1 is the year 850
           END IF
           IF(print_type=='unagg') THEN
             IF(compute_next_year) THEN
               DO v=1,NV
                  WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_in_new(v,:,:),mask)/num_land
               END DO
             ELSE
               DO v=1,NV
                  WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_in(v,:,:),mask)/num_land
               END DO
             END IF
           ELSE IF(print_type=='agg') THEN
             DO v=1,NV2
                WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_out(v,:,:)/100.0,mask)/num_land 
             END DO
           ELSE IF(print_type=='sdgvm') THEN
             DO v=1,NE
                WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_out_SDGVM(v,:,:)/100.0,mask)/num_land 
             END DO
           END IF
           WRITE(*,*) ! Assumes default "ADVANCE='yes'".
        END IF
      END DO

      ! Close the file, freeing all resources.
      CALL CHECK( NF90_CLOSE(ncid) )
      CALL CHECK( NF90_CLOSE(ncid_t) )

      PRINT *,"*** SUCCESS reading file ", wdh//fname, "! "

    CONTAINS
      SUBROUTINE CHECK(status)
        INTEGER, INTENT ( IN) :: status
        
        IF(status /= NF90_NOERR) THEN 
          PRINT *, TRIM(NF90_STRERROR(status))
          STOP "Stopped"
        END IF
      END SUBROUTINE CHECK  

    !from Anthony Walker's regrid.R
    !reduce_res <- function(mat,red=2){
    !  ! reduce resolution by reduction factor
    !  ! red must be an integer - defaults to 2 reducing resolution by half
    !  ! discounts NAs from calculation
    !
    !  !row indexes
    !  ie <- seq(1,dim(mat)[1],red)
    !  io <- seq(2,dim(mat)[1],red)
    !  !column indexes
    !  je <- seq(1,dim(mat)[2],red)
    !  jo <- seq(2,dim(mat)[2],red)
    !
    ! !NA matrix
    !  mat_na <- is.na(mat)
    !  !zero NAs
    !  mat[is.na(mat)] <- 0
    !
    !  !sum rows
    !  y1   <- mat[io,]+mat[ie,]
    !  y1_na <- mat_na[io,]+mat_na[ie,]
    !  !sum cols
    !  y2    <- y1[,jo]+y1[,je]
    !  y2_na <- y1_na[,jo]+y1_na[,je]
    !
    !  n <- red^2 - y2_na
    !
    !  y2/n
    !}

    FUNCTION JOIN_HYDE(v) RESULT( result_join_hyde ) 
      ! takes a 15 element vector, 1:10 - SDGVM PFTs from ESA data, 11:15 HYDE aggregated land cover 
      ! and combines into a 10 element vector of SDGVM PFTs
       
      ! SDGVM PFTs: 1 BARE, 2 Ev_Bl, 3 Dc_Bl, 4 Ev_Nl, 5 Dc_Nl, 6 Shrub, 7 C3, 8 C4, 9 C3crop, 10 C4crop
      ! HYDE aggregated land cover: 11 forest, 12 non-forest, 13 C3 crop, 14 C4 crop, 15 pasture   
      
      ! first forest area and potential forest area from HYDE is used to adjust forest area in ESA
      ! then C3 and C4 croplands in HYDE are used to adjust ESA
      ! C3 and C4 grassland area is then adjusted to account for changes in pasture and natural grassland cover
      ! rangelands in hyde and shrubs in esa are not explicitly included in the calculation of change 
      ! though a proportion of them could be accounted for implicitly if rangelands habour some forest vegetation
      IMPLICIT NONE
      REAL, INTENT(IN) :: v(15)
      REAL, DIMENSION(10) :: result_join_hyde
      REAL :: ov(10), nofcov, nogcov, for_esa, diff, esahyde_cov, hyde_cov, past_esa
      INTEGER :: fid_min,fid_max,f
      LOGICAL :: past_err 
      REAL, PARAMETER :: NA= -999.0
      ! PRINT(v)
      
      !if(is.finite(sum(v(1:10))) & not is.finite(sum(v(11:15))) ) then 
      IF((SUM(v(1:10)) > 0.0) .AND. (SUM(v(1:10)) < 1000.0)  .AND. &
        .NOT.((SUM(v(11:15)) > 0.0) .AND. (SUM(v(11:15)) < 1000.0)) ) THEN 
        ! when ESA has landcover but hyde not
        ! arises from different land masks
        result_join_hyde= v(1:10)
        
      !else if(is.finite(sum(v))) then
      ELSE IF((SUM(v) > 0.0) .AND. (SUM(v) < 1000.0) ) THEN 
        
        ! create output vector 
        ! print(v)
        ov = v(1:10)
        
        ! subscripts for ESA forest PFTs
        fid_min = 2
        fid_max = 6
        
        ! stores for increase in ESA but when either no forest or grass cover in ESA 
        nofcov  = 0.0
        nogcov  = 0.0
        
        ! change forest cover
        ! ESA forest cover
        for_esa = SUM(v(fid_min:fid_max))
        ! If forest cover is greater in hyde
        IF(v(11)>for_esa) THEN 
          diff   = v(11) - for_esa
          IF(for_esa>0) THEN 
            ! increase forest cover proportionally
            DO f=fid_min,fid_max
               ov(f)  = ov(f) + diff * (v(f) / for_esa)
            END DO
          ELSE 
            ! no forest cover in ESA - forest cover stored for post-processing
            ! write(paste('no ESA forest','BARE:',ov[1],'C3 CROP:',ov[9],'C4 CROP:',ov[10]),'1nofor_error.txt',append=T)
            nofcov = diff
          END IF 
        ELSE IF(for_esa>v(11)) THEN 
          ! forest cover greater in ESA
          diff   = for_esa - v(11)
          ! decrease forest cover proportionally
          DO f=fid_min,fid_max
             ov(f)  = ov(f) - diff * (v(f) / for_esa)
          END DO
        END IF 
        
        
        ! increase or decrease crop cover 
        ! C3 crops
        IF(v(9) .NE. v(13))  ov(9)  = v(13)
        
        ! C4 crops
        IF(v(10) .NE. v(14)) ov(10) = v(14)
        
        
        ! increase or decrease grass (pasture) cover 
        ! total cover so far in hyde ESA combined dataset - accounts for forest cover as yet unassigned to a PFT
        esahyde_cov = SUM(ov) + nofcov
        
        ! total cover in hyde
        hyde_cov    = SUM(v(11:15))
        
        ! ESA grasslands (potential pasture)
        past_esa    = v(7) + v(8)
        
        ! if new cover is different from hyde cover
        past_err = .FALSE. 
        IF(esahyde_cov .NE. hyde_cov) THEN 
          diff   = hyde_cov - esahyde_cov
          IF(diff<0.0) THEN 
            IF(ABS(diff)<past_esa) THEN 
              
              ! change C3/C4 grass cover proportionally
              ov(7)  = ov(7) + diff * (v(7) / past_esa)
              ov(8)  = ov(8) + diff * (v(8) / past_esa)
              
            ELSE
              
              ! ESA pasture needs reduction but reduction is more than pasture cover
              ov(7) = 0.0
              ov(8) = 0.0
              
              ! assume remaining reduction from bare ground
              diff = abs(diff) - past_esa
              IF(ov(1)>=diff) THEN
                 ov(1) = ov(1) - diff
              ELSE
                 ov(1) = 0.0
              END IF
            END IF   
            
          ELSE IF(diff>0.0) THEN 
            IF(past_esa>0.0) THEN 
              
              ! change C3/C4 grass cover proportionally
              ov(7)  = ov(7) + diff * (v(7) / past_esa)
              ov(8)  = ov(8) + diff * (v(8) / past_esa)
              
            ELSE 
              
              ! no grass cover in ESA - grass cover stored for post-processing
              ! write(paste('no ESA grass','BARE:',ov[1],'C3 CROP:',ov[9],'C4 CROP:',ov[10]),'1nopast_error.txt',append=T)
              nogcov = diff
            END IF 
          END IF 
       END IF 
        
        ! test for cover sum >100%
        IF( SUM(ov(1:10))-100.0 > 1e-3 ) PRINT *,'Error,cov:',SUM(ov(1:10)),hyde_cov,'. Bare:',v(1) !'1cov_error.txt'
        
        ! assign as yet unassigned forest cover as negative value to DcBl
        IF(nofcov>0) ov(3)  = -nofcov
        ! assign as yet unassigned grass cover as negative value to C3 grass
        IF(nogcov>0) ov(7)  = -nogcov
        
        result_join_hyde = ov(1:10)
       ELSE 
        result_join_hyde = (/ NA,NA,NA,NA,NA,NA,NA,NA,NA,NA /) 
       END IF 
    END FUNCTION JOIN_HYDE

    ! Assigns forest or grass pfts to ESA HYDE according to latitude when no forest or grass PFTs exist in ESA but do in HYDE 
    !function f_lat_assignPFT(j,m,loff) result (outvec) 
    ! ! this function assigns forest and grassland cover to an appropriate PFT in the combined HYDE ESA dataset
    ! ! when there was no forest or grass cover in the ESA data.
    ! ! in this case the forest cover is assigned as a negative value to deciduous braodleaved PFT m[3,] 
    ! ! in this case the grass  cover is assigned as a negative value to C3 grass PFT m[7,]
    ! ! if in temperate latitudes this negative value is simply switched to positive
    ! ! if in tropical latitudes this negative value is switched to positive and assigned to evergreen broadleaved PFT or C4 grass 
    ! 
    ! ! m is a pft x lat matrix
    ! ! loff is the change in subscript by lat - 0 for temperate (DcBL & C3), 1 for tropical (EvBl & C4)
    ! 
    ! implicit none
    ! real, intent(in) :: v(15)
    ! real, intent(out) :: outvec(10)

    ! ! if there was no forest or grass cover in the ESA data 
    ! if(any(m(:,j)<0))
    !   ! switch DcBl and EvBl PFT to allow loff to work for both forest and grass
    !   m(2:3,j) = m(3:2,j)
    !   ! subscripts of negative cover
    !   sub      = which(m[,j]<0)
    !   ! take absolute value of negative covers and assign to appropriate PFT
    !   m(sub+loff(j),j) = abs(m(sub,j))
    !   ! if tropical PFT is appropriate zero negative cover in temperate PFT
    !   if(loff(j) .NE. 0) m(sub,j) = 0.0
    !   ! switch DcBl and EvBl PFT back to original placement in vector
    !   m(2:3,j) = m(3:2,j)
    !   ! return vector
    !   outvec = m(:,j)
    ! else
    !   outvec = m(:,j)
    ! end if 
    !end function f_lat_assignPFT

    END PROGRAM states_convertSDGVM_v5b
