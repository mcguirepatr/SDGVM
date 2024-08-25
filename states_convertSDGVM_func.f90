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
!       gfortran states_convertSDGVM_func_v7a.f90 `nf-config --fflags --flibs`

!versions:
!    v5b: working version that replicates original R code
!    v6a: enhancement to keep track of both primary & secondary tree FTs
!    v7a: enhancement to keep track of gross transitions
!    v7b: enhancement to get wdg,pname,pname_t,SYR, and FYR from input.dat file read in the sdgvm0.f 

      MODULE FUNCTIONS_CLU
      CONTAINS
      SUBROUTINE states_convertSDGVM_func(SYR,FYR,SYR_FIRE,FYR_FIRE, &
           X0,Y0,DXY,get_transitions,compute_next_year,prescr_fire, &
           pname,pname_t,pname_f,wdg,data_out_SDGVM, &
           data_out_AggHydeTransitions,data_out_AggHarvest, &
           data_out_fire,debug )
      USE netcdf
      IMPLICIT NONE

      INTEGER, PARAMETER :: NS = 17
      INTEGER, PARAMETER :: NV2 = 7 
      INTEGER SYR,FYR,SYR_FIRE,FYR_FIRE,X0,Y0,DXY 
      LOGICAL :: compute_next_year,get_transitions,prescr_fire
      REAL*8, DIMENSION(FYR-SYR+1,NS,DXY,DXY) ::  data_out_SDGVM
      REAL*8, DIMENSION(FYR-SYR+1,NV2,NV2,DXY,DXY) ::  &
                         data_out_AggHydeTransitions
      REAL*8, DIMENSION(FYR-SYR+1,NV2,DXY,DXY) ::  data_out_AggHarvest
      REAL*8, DIMENSION(FYR_FIRE-SYR_FIRE+1,DXY,DXY,12) :: data_out_fire

      INTEGER, PARAMETER :: NX = 720, NY = 360
      ! half-resolution version computed with: module load jasppy ! on JASMIN
      !                                        cdo gridboxmean,2,2 transitions.nc transitions2b.nc
      !                                        cdo gridboxmean,2,2 states.nc states2b.nc
      !CHARACTER (LEN = *), PARAMETER :: fname = "states2b.nc" !half-resolution version of states.nc
      !CHARACTER (LEN = *), PARAMETER :: fname_t = "transitions2b.nc" !half-resolution version of transitions.nc
      !CHARACTER (LEN = *), PARAMETER :: wdh = &
      ! "/gws/nopw/j04/nexcs/pmcguire/TRENDYv10/db/LUH2_GCB_2021/"
      !CHARACTER (LEN = *), PARAMETER :: wdg = &
      ! "/gws/nopw/j04/nexcs/pmcguire/sdgvmD/data/land_use/global/ESACCILCP2014/30min/"
      CHARACTER (LEN = *) :: pname !path to half-resolution version of states.nc
      CHARACTER (LEN = *) :: pname_t !path to half-resolution version of transitions.nc
      CHARACTER (LEN = *) :: pname_f !path to 0.5-degree prescribed-fire files 
      CHARACTER (LEN = *) :: wdg
      CHARACTER (LEN = *), PARAMETER :: fname_s = "cont_lu"
      CHARACTER (LEN = *), PARAMETER :: esadate = "2009"
      !CHARACTER (LEN = *), PARAMETER :: print_type='unagg'
      CHARACTER (LEN = *), PARAMETER :: print_type='agg'
      !CHARACTER (LEN = *), PARAMETER :: print_type='esa'
      !CHARACTER (LEN = *), PARAMETER :: print_type='sdgvm'

      INTEGER, PARAMETER :: NT = 1172, NV = 14
      REAL*8, PARAMETER    :: misval = 1e19
      INTEGER, PARAMETER :: NE = 10, NE2 = 16
      INTEGER, PARAMETER :: NVT = 118 
      INTEGER, PARAMETER :: DT = 20 !number of years to skip between prints
      !INTEGER, PARAMETER :: SINDEX = 21 !starting year for prints
      !INTEGER, PARAMETER :: SINDEX = 1001 !starting year for prints ! for 1851
      !INTEGER, PARAMETER :: SINDEX = 850 !starting year for prints ! for 1700
      INTEGER, PARAMETER :: SYR0 = 850 !starting year for NETCDF data 
      INTEGER, PARAMETER :: SYR0_FIRE = 1901 !starting year for NETCDF data for fire
      INTEGER SINDEX,FINDEX
      REAL*8 :: data_in(NV, DXY, DXY), dummya(DXY, DXY)
      REAL*8 :: data_in_old(NV, DXY, DXY)
      REAL*8 :: data_in_new(NV, DXY, DXY)
      REAL*8 :: data_in_t(NVT, DXY, DXY)
      REAL*8 :: data_in_fire(DXY, DXY, 12)
      REAL*8 :: data_out(NV-2, DXY, DXY)
      REAL*8 :: data_t_agg(NV2, NV2, DXY, DXY)
      REAL*8 :: data_out_harvest(NV2, DXY, DXY)
      LOGICAL :: mask(DXY, DXY)
      LOGICAL :: mask_SDGVM(DXY, DXY)
      LOGICAL :: esamask(NE2,DXY, DXY)
      LOGICAL :: debug !used to print out more debugging info
      ! setup ESA arrays
      INTEGER :: esaarray_global(NE2,NX,NY)
      REAL*8 :: esaarray(NE2,DXY,DXY)
      CHARACTER (LEN = 200) :: fname_s2
      INTEGER, PARAMETER :: NA= -999
      INTEGER :: minx,miny,minii,minjj
      INTEGER :: maxx,maxy,maxii,maxjj
      INTEGER :: miny_fire,minjj_fire
      INTEGER :: maxy_fire,maxjj_fire
      ! This will be the netCDF ID for the file and data variable.
      INTEGER :: ncid, varid(NV)
      INTEGER :: ncid_t, varid_t(NVT)
      INTEGER :: ncid_f, varid_f
      INTEGER :: ndims_f, dimids_f(10), dimlen
      CHARACTER (LEN = NF90_MAX_NAME) :: dimname ! Dimension nam  

      INTEGER :: num_land,blank,num_land_hyde,num_land_SDGVM,num_fire

      ! Loop indexes, and error handling.
      INTEGER :: x, y, t, v, v2, v3, i, lon, lat, year_index, shift_year
      INTEGER :: year,jj
      CHARACTER (LEN = 5) :: year_string



      ! processing parameter
      ! latitudinal offset for assigning correct forest or grass PFT when none exist in ESA
      INTEGER :: lat_offset(360) 

      ! ESA PFTs
      CHARACTER(LEN=7),PARAMETER :: varname_esa_pfts(NE)=(/ '  BARE',' Ev_Bl',' Dc_Bl',' Ev_Nl', &
       ' Dc_Nl', ' Shrub','    C3','    C4','C3crop','C4crop' /)

      CHARACTER(LEN=7),PARAMETER :: varname(NV)=(/'primf', 'primn', 'secdf', 'secdn', &
          'c3ann', 'c4ann', 'c3per', 'c4per', 'c3nfx', 'pastr', 'range', 'urban', &
          'secmb', 'secma'/)
     ! Anthony's original order in R code:
     ! CHARACTER(LEN=7),PARAMETER :: varname(NV)=(/'primf', 'secdf', 'secdn', 'primn', &
     !     'c3ann', 'c3per', 'c3nfx', 'c4ann', 'c4per', 'pastr','range',               &
     !     'urban','secmb', 'secma'/)
      CHARACTER(LEN=7),PARAMETER :: varname_sdgvm(NS)=(/'  BARE',' Ev_Bp',' Dc_Bp',' Ev_Np',' Dc_Np', &
          ' Shrup','   C3p','   C4p','C3crop','C4crop',' Ev_Bs',' Dc_Bs',' Ev_Ns',' Dc_Ns',' Shrus',  &
          '   C3s','   C4s'/)
      !CHARACTER(LEN=7),PARAMETER :: varname_sdgvm(NS)=(/'  BARE',   &
      !    '  CITY','   C3p','   C4p','C3crop', &
      !    'C4crop','   C3s','   C4s',' Ev_Bp',' Ev_Np', &
      !    ' Dc_Bp',' Dc_Np',' Ev_Bs',' Ev_Ns',' Dc_Bs', &
      !    ' Dc_Ns'/)


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

      INTEGER aggmap(NV-2) !don't include 'secmb', 'secma'

      CHARACTER(LEN=2),PARAMETER :: varname_f='ba'

      ! Orig: NV2=6
      ! aggregate Hyde landcover types 
      !CHARACTER(LEN=9),PARAMETER :: varname2(NV2)=(/'   forest', ' nforestr', '   c3crop', '   c4crop', &
      !    '   pastnr', '    urban' /)
      ! primary and secondary forest
      ! primary and secondary non-forest and rangelands
      ! c3 crops
      ! c4 crops
      ! pasture and not rangelands
      ! urban

      ! aggmap(1)  =  1   ! primary forest       -> forest
      ! aggmap(2)  =  2   ! primary non-forest   -> nforestr
      ! aggmap(3)  =  1   ! secondary forest     -> forest
      ! aggmap(4)  =  2   ! secondary non-forest -> nforestr
      ! aggmap(5:7)  =  3 ! c3 crops
      ! aggmap(8:9)  =  4 ! c4 crops
      ! !aggmap(10:11) = 5  ! pasture and rangelands
      ! aggmap(10) =  5   ! pasture and not rangelands
      ! aggmap(11) =  2   ! rangelands - add to primary non-forest
      ! aggmap(12) =  6   ! urban

      ! New: NV2=8
      ! HYDE aggregated land cover: 
      !CHARACTER(LEN=9),PARAMETER :: varname2(NV2)=(/'    primf', &
      !    '    primn', '    secdf', '    secdn', &
      !    '   c3crop', '   c4crop', '   pastnr', '    urban' /)

      ! aggmap(1)  =  1   ! primary forest
      ! aggmap(2)  =  2   ! primary non-forest
      ! aggmap(3)  =  3   ! secondary forest
      ! aggmap(4)  =  4   ! secondary non-forest
      ! aggmap(5)  =  5 ! c3 crops
      ! aggmap(6)  =  6 ! c4 crops
      ! aggmap(7)  =  5 ! c3 crops 
      ! aggmap(8)  =  6 ! c4 crops
      ! aggmap(9)  =  5 ! c3 crops 
      ! !aggmap(10:11) = 7  ! pasture and rangelands
      ! aggmap(10) =  7   ! pasture and not rangelands
      ! !aggmap(11) =  2   ! rangelands - add to primary non-forest
      ! aggmap(11) =  4   ! rangelands - add to secondary non-forest
      ! aggmap(12) =  8   ! urban

      ! New2: NV2=7
      ! HYDE aggregated land cover, with pastnr put in secdn: 
      CHARACTER(LEN=9),PARAMETER :: varname2(NV2)=(/'    primf', &
          '    primn', '    secdf', '    secdn', &
          '   c3crop', '   c4crop', '    urban' /)

      aggmap(1)  =  1   ! primary forest
      aggmap(2)  =  2   ! primary non-forest
      aggmap(3)  =  3   ! secondary forest
      aggmap(4)  =  4   ! secondary non-forest
      aggmap(5)  =  5 ! c3 crops
      aggmap(6)  =  6 ! c4 crops
      aggmap(7)  =  5 ! c3 crops 
      aggmap(8)  =  6 ! c4 crops
      aggmap(9)  =  5 ! c3 crops 
      !!aggmap(10:11) = 7  ! pasture and rangelands
      !aggmap(10) =  7   ! pasture and not rangelands
      !!aggmap(11) =  2   ! rangelands - add to primary non-forest
      !aggmap(11) =  4   ! rangelands - add to secondary non-forest
      aggmap(10:11) = 4  ! pasture and rangelands - add to secondary non-forest
      !aggmap(12) =  8   ! urban
      aggmap(12) =  7   ! urban

      ! latitudinal offset for assigning correct forest or grass PFT when none exist in ESA
      lat_offset(1:132)      = 0
      lat_offset(133:133+96) = 1
      lat_offset(133+97:360) = 0

      IF(compute_next_year) THEN
          shift_year = -1 
      ELSE
          shift_year = 0
      ENDIF

      !PCM since currently, there is no wrapping of the land cover at LON=0, we
      !will keep this for now, for this extension, too.
      minx=MAX(X0,1) 
      miny=MAX(Y0,1) 
      maxx=MIN(X0+DXY-1,NX) 
      maxy=MIN(Y0+DXY-1,NY) 
      minii=MAX(DXY/2-X0,1) 
      minjj=MAX(DXY/2-Y0,1) 
      maxii=MIN(NX-X0+1,DXY) 
      maxjj=MIN(NY-Y0+1,DXY) 

      maxy_fire=NY-miny  !the fire y axis is flipped
      miny_fire=NY-maxy
      maxjj_fire=maxjj
      minjj_fire=minjj
      IF(debug) THEN
        PRINT *,'X0,Y0,minx,maxx,miny,maxy,minii,maxii,minjj,maxjj'
        PRINT *,X0,Y0,minx,maxx,miny,maxy,minii,maxii,minjj,maxjj
        PRINT *,'FIRE:miny,maxy,minjj,maxjj'
        PRINT *,miny_fire,maxy_fire,minjj_fire,maxjj_fire
      ENDIF

      ! Open ESA CCILCP 2014 dataset 
      !setwd(wdg)
      DO i=1,NE 
        esaarray(i,:,:) = 255.0
        esaarray_global(i,:,:) = 255
        WRITE (fname_s2, "(A,I0,A)") wdg(1:blank(wdg))//"/"//fname_s// &
                         "-", i,"-"//esadate//".dat" 
        !esav          <- scan(fname_s2)
        !! esamat <- if(i==1) esav else cbind(esamat,esav)
        !esaarray(,,i) <- as.matrix(esav,nrow=lon_ress)
        OPEN(15,FILE=fname_s2,STATUS='old')
        READ(15,*) esaarray_global(i,:,:)
        esaarray(i,minii:maxii,minjj:maxjj) =  &
                     REAL(esaarray_global(i,minx:maxx,miny:maxy),8)
        CLOSE(15)
      END DO 

      ! convert missing values      
      WHERE (esaarray .EQ. 255.0) esaarray = NA      
      WHERE (esaarray .NE. NA)
        esamask = .TRUE.      
      ELSEWHERE
        esamask = .FALSE.      
      END WHERE

      num_land = COUNT( esamask(1,3:4,3:4) .EQV. .TRUE.)
      IF(debug) THEN
        PRINT *, 'ESA num_land=',num_land,'num_tot=',2*2 
!        PRINT *, pname(1:blank(pname)) 
      ENDIF


      ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
      ! the file.
      CALL CHECK( NF90_OPEN(pname(1:blank(pname)), NF90_NOWRITE, ncid) )

      DO v=1,NV
      ! Get the varid of the data variable, based on its name.
        CALL check( nf90_inq_varid(ncid, varname(v), varid(v)) )
      END DO

      CALL CHECK( NF90_OPEN(pname_t(1:blank(pname_t)), NF90_NOWRITE, &
          ncid_t) )

      DO v=1,NVT
      ! Get the varid of the data variable, based on its name.
        CALL CHECK( NF90_INQ_VARID(ncid_t, varname_t(v), varid_t(v)) )
      END DO

      IF(prescr_fire) THEN
        SINDEX = SYR_FIRE - SYR0_FIRE + 1 + shift_year
        FINDEX = FYR_FIRE - SYR0_FIRE + 1 + shift_year
        DO t=SINDEX,FINDEX 
          year_index = t - SINDEX + 1 
          year       = SYR_FIRE + year_index - 1  
          write (year_string,'(I4)') year 
          IF(debug) then
            write(*,*) 'states_*: prescr_fire', year 
          END IF
!      /gws/nopw/j04/nexcs/pmcguire/TRENDYv13/db/burned_area/&
!      &global_monthly_burned_area_fraction_05deg_1901.nc
!      pname_f = 
!      /gws/nopw/j04/nexcs/pmcguire/TRENDYv13/db/burned_area/&
!      &global_monthly_burned_area_fraction_05deg_
          CALL CHECK( NF90_OPEN(pname_f(1:blank(pname_f))//   &
                      trim(year_string)//'.nc',NF90_NOWRITE, ncid_f) )

          ! Get the varid of the data variable, based on its name.
          CALL CHECK( NF90_INQ_VARID(ncid_f, varname_f, varid_f) )

!          CALL CHECK( NF90_INQUIRE_VARIABLE(ncid_f, varid_f,  &
!                      ndims=ndims_f) )
!          CALL CHECK( NF90_INQUIRE_VARIABLE(ncid_f, varid_f,  &
!                      dimids=dimids_f) )
!          WRITE(*,*) ndims_f
!          DO i = 1, ndims_f
!              CALL CHECK( NF90_INQUIRE_DIMENSION(ncid_f, dimids_f(i), &
!                          dimname, dimlen) )
!              WRITE(*,*) dimids_f(i),dimlen,dimname
!          END DO

          data_in_fire(1:4,1:4,1:12) = NA
!          CALL CHECK( NF90_GET_VAR(ncid_f, varid_f,           &
!                  data_in_fire(1:4,1:4,1:12), &
!                  start=[180,180,1], count=[4,4,12]) )

!    Thie fire data has a flipped y axis, so we use miny_fire array values

          CALL CHECK( NF90_GET_VAR(ncid_f, varid_f,           &
                  data_in_fire(minii:maxii,                   &
                  minjj_fire:maxjj_fire,1:12),                &
                  start=[minx,miny_fire,1], count=[maxii-minii+1,  &
                  maxjj_fire-minjj_fire+1,12]) )
          CALL CHECK( NF90_CLOSE(ncid_f) )

          WHERE(ABS(data_in_fire(:,:,1)) <= 1.00 ) !just check for january
            mask = .TRUE.
          ELSEWHERE
            mask = .FALSE.
          END WHERE

          num_fire = COUNT( mask(3:4,3:4) .EQV. .TRUE.)
          IF(num_fire.eq.0) THEN
            data_in_fire = NA 
!            data_out_fire = 255.0
!            data_out_AggHydeTransitions = 255.0 
!            data_out_AggHarvest = 255.0 
!            data_out_SDGVM = 255.0
             
!            RETURN
          END IF

          WHERE (isNAN(data_in_fire(:,:,1:12)))
               data_in_fire = NA
          ENDWHERE
          data_out_fire(year_index,:,:,1:12) = data_in_fire(:,:,1:12) 
        END DO
      END IF
      IF(debug) THEN
       PRINT *, 'num_fire=',num_fire,'num_tot=',2*2
       WRITE(*,*)'in states_*; data_out_fire'
!      DO jj=1,FYR_FIRE-SYR_FIRE+1
       WRITE(*,*)'3,3'
       DO jj=1,4
         WRITE(*,'(12F11.4)')data_out_fire(jj,3,3,1:12)
       END DO
       WRITE(*,*)'3,4'
       DO jj=1,4
         WRITE(*,'(12F11.4)')data_out_fire(jj,3,4,1:12)
       END DO
       WRITE(*,*)'4,3'
       DO jj=1,4
         WRITE(*,'(12F11.4)')data_out_fire(jj,4,3,1:12)
       END DO
       WRITE(*,*)'4,4'
       DO jj=1,4
         WRITE(*,'(12F11.4)')data_out_fire(jj,4,4,1:12)
       END DO
      ENDIF

      !DO t=ST,NT,1 
      SINDEX = SYR - SYR0 + 1 + shift_year
      FINDEX = FYR - SYR0 + 1 + shift_year
      DO t=SINDEX,FINDEX 
      !  write(*,*)'ST1',t,SYR,SYR0,shift_year,FINDEX,FYR
        year_index = t - SINDEX + 1 
        data_out(:,:,:) = 0.0
        data_t_agg(:,:,:,:) = 0.0
        data_out_harvest(:,:,:) = 0.0
        IF(t > SINDEX) THEN 
          data_in_old(:,:,:)=data_in(:,:,:)
        END IF
      ! Read the data.
        DO v=1,NV
          data_in(v,:,:) = NA
          CALL CHECK( NF90_GET_VAR(ncid, varid(v),   &
                  data_in(v,minii:maxii,minjj:maxjj), &
                  start=[minx,miny,t],                &
                  count=[maxii-minii+1,maxjj-minjj+1,1]) )
      !  print *,'Finished reading data'
          IF( t == SINDEX .AND. v == 1) THEN
            WHERE(ABS(data_in(1,:,:)) <= 1.00 )
              mask = .TRUE.
            ELSEWHERE
              mask = .FALSE.
            END WHERE
            num_land_hyde = COUNT( mask(3:4,3:4) .EQV. .TRUE.)
            IF(debug) THEN
              PRINT *, 'num_land_hyde=',num_land_hyde,'num_tot=', &
                 maxii*maxjj,'val=',data_in(1,3,3),'mask=',mask(3,3) 
              PRINT *, '          1     2    3     4'
              PRINT *, '    1 ',data_in(1,:,1)
              PRINT *, '    2 ',data_in(1,:,2)
              PRINT *, '    3 ',data_in(1,:,3)
              PRINT *, '    4 ',data_in(1,:,4)
            ENDIF
            IF(num_land_hyde.eq.0) THEN
              data_out_AggHydeTransitions = 255.0 
              data_out_AggHarvest = 255.0 
              data_out_SDGVM = 255.0
              data_out_fire = 255.0
              RETURN
            END IF
          END IF

        ! for now, use regridded states.nc file from CDO, so we skip the next two steps
        ! regrid from 0.25 to 0.5
        !  dummya = reduce_res(data_in(v,:,:))

        ! convert missing values      
        !  where (dummya > misval) mask = .false.      

        END DO


        WHERE (isNAN(data_in(:,:,:)))
               data_in = NA
        ENDWHERE

        IF(get_transitions .AND. (t .LT. FINDEX)) THEN
          ! the secondary vegetation doesn't match up unless we
          ! add the transitions to the states for each time step (rather than
          ! accumulating all the transitions in an open loop).
          ! if(t == SINDEX) then 
           IF(compute_next_year) THEN
             data_in_new(:,:,:)=data_in(:,:,:) 
           ENDIF
          ! end if
          ! write(*,*)'ST1b tr',t,SYR,SYR0,shift_year,FINDEX,FYR

           DO v=1,NVT-5 ! skip bioh 
            data_in_t(v,:,:) = NA        !data_in_t = transitions matrix element for transition with the name varname_t(v)
            CALL CHECK( NF90_GET_VAR(ncid_t, varid_t(v),  &
                      data_in_t(v,minii:maxii,minjj:maxjj), &
                      start=[minx,miny,t], &
                      count=[maxii-minii+1,maxjj-minjj+1,1]) )
            !write(*,*)v,NVT-5,NVT-10,NVT,varname_t(v)
            IF(v.LE.NVT-10) THEN
               from_t = varname_t(v)(1:5)   !from_t = the state from which the transition is coming
               to_t   = varname_t(v)(10:14) !to_t   = the state to   which the transition is going 
            ELSE
               to_t = 'harv'
               IF(     varname_t(v)(1:10).EQ.'primf_harv') THEN 
                  from_t = 'primf'
               ELSE IF(varname_t(v)(1:10).EQ.'primn_harv') THEN 
                  from_t = 'primn'
               ELSE IF(varname_t(v)(1:10).EQ.'secmf_harv') THEN  !wood-area harvest from mature forest
                  from_t = 'secdf' !don't aggregate young and mature; only count mature
!               ELSE IF(varname_t(v)(1:10).EQ.'secyf_harv') THEN  !wood-area harvest from young forest 
!                  from_t = 'secdf' !don't aggregate young and mature
               ELSE IF(varname_t(v)(1:10).EQ.'secnf_harv') THEN  !wood-area harvest from non-forest 
                  from_t = 'secdn' !note the different spelling from 'secnf_harv'
               ENDIF
            ENDIF

            dummya = data_in_t(v,:,:)
            WHERE (isNAN(dummya(:,:)))
               dummya = NA
            ENDWHERE
            ! there are some gridcells in northern South America that
            ! have very small negative values for the c3nfx_to_c4ann
            ! transition in the year 1700, for example
            WHERE ((dummya(:,:).NE.NA).AND.(dummya(:,:).LT.0.0d0))
               dummya = 0.0  
            ENDWHERE
            
            DO v2=1,NV-2
              IF(varname(v2) == from_t) THEN
               IF(compute_next_year) THEN
                 data_in_new(v2,:,:)= data_in_new(v2,:,:) - dummya 
               ENDIF

               ! aggregate harvests for Hyde landcover types;
               ! this does not combine mature forest with young forest;
               ! only wood harvest from mature forest is counted
               IF(to_t == 'harv') THEN
                 IF(varname_t(v)(1:10).NE.'secyf_harv') THEN
                   data_out_harvest(aggmap(v2),:,:) = &
                        data_out_harvest(aggmap(v2),:,:) + dummya
                   if(t == SINDEX) then
                    IF(debug) THEN
                      PRINT *, varname_t(v), from_t,to_t, v2, &
                               aggmap(v2),100.0*dummya(3,3),  &
                            data_out_harvest(aggmap(v2),3,3)
                    ENDIF
                   end if
                 ENDIF
               ELSE
               ! aggregate transitions for  Hyde landcover types 
                 DO v3=1,NV-2
                  !don't allow aggregated classes to have non-zero transition
                  !probabilities to themselves (i.e. c3per_to_c3ann)
                  IF(varname(v3) == to_t .AND. aggmap(v2).NE.aggmap(v3)) THEN
                    !print 3,3 coords of 4x4 region
                    data_t_agg(aggmap(v2),aggmap(v3),:,:) = data_t_agg(aggmap(v2),aggmap(v3),:,:) + dummya
                    if(t == SINDEX) then
                     IF(debug) THEN
                     IF(num_land_hyde.gt.0)THEN
                      PRINT *, varname_t(v), from_t,varname(v3), v2, v3,aggmap(v2),aggmap(v3), &
                              100.0*SUM(dummya(3:4,3:4),mask(3:4,3:4))/num_land_hyde,  &
                              SUM(data_t_agg(aggmap(v2),aggmap(v3),3:4,3:4),mask(3:4,3:4))/num_land_hyde
                     ELSE
                      PRINT *,'No Valid Hyde Data (in transitions code)'
                     ENDIF
                     ENDIF
                    end if
                  END IF
                 END DO
               END IF
              END IF
              IF(varname(v2) == to_t) THEN
               IF(compute_next_year) THEN
                 data_in_new(v2,:,:)= data_in_new(v2,:,:) + dummya 
               ENDIF
               !don't aggregate transitions twice
               !! aggregate Hyde landcover types 
               !DO v3=1,NV
               ! IF(varname(v3) == from_t) THEN
               !   PRINT *, varname_t(v), varname(v3), to_t, v3, v2, dummya(1,1)
               !   data_t_agg(aggmap(v3),aggmap(v2),:,:) = data_t_agg(aggmap(v3),aggmap(v2),:,:) + dummya
               ! END IF
               !END DO
              END IF
            END DO
           END DO
        ! change to percent      
           WHERE (data_t_agg.GE.0.0d0)
             data_t_agg = data_t_agg * 100.0      
           ENDWHERE
           WHERE (data_t_agg.LT.0.0d0)
             data_t_agg = NA      
           ENDWHERE

           WHERE (data_out_harvest.GE.0.0d0)
             data_out_harvest = data_out_harvest * 100.0      
           ENDWHERE
           WHERE (data_out_harvest.LT.0.0d0)
             data_out_harvest = NA 
           ENDWHERE


           if(debug) then

             write(*,FMT="(A)") 'STSS1'
             IF(num_land_hyde.GT.0)THEN
              write(*,*)'     ','     primf', '     primn', '     secdf', &
                       '     secdn', '    c3crop', '    c4crop', &
                       '     urban', '     harv' 
!          use the 3:4,3:4 coordinates of the 4x4 region for averaging
              DO v3=1,7
               IF(v3.EQ.1)THEN
                 write(*,FMT="(A)",ADVANCE='no') ' primf'
               ELSE IF(v3.EQ.2)THEN
                 write(*,FMT="(A)",ADVANCE='no') ' primn'
               ELSE IF(v3.EQ.3)THEN
                 write(*,FMT="(A)",ADVANCE='no') ' secdf'
               ELSE IF(v3.EQ.4)THEN
                 write(*,FMT="(A)",ADVANCE='no') ' secdn'
               ELSE IF(v3.EQ.5)THEN
                 write(*,FMT="(A)",ADVANCE='no') 'c3crop'
               ELSE IF(v3.EQ.6)THEN
                 write(*,FMT="(A)",ADVANCE='no') 'c4crop'
               ELSE IF(v3.EQ.7)THEN
                 write(*,FMT="(A)",ADVANCE='no') ' urban'
               END IF
               DO v2=1,7
                 write(*,FMT="(E10.3)",ADVANCE='no') SUM(data_t_agg(v3,v2,3:4,3:4), &
                   mask(3:4,3:4))/num_land_hyde
               END DO
               write(*,FMT="(E10.3)",ADVANCE='no') SUM(data_out_harvest(v3,3:4,3:4),&
                  mask(3:4,3:4))/num_land_hyde
               write(*,*) !advance = 'yes'
             END DO
            ELSE
             write(*,*) 'No Valid SDGVM land points'
            ENDIF

           endif
        ELSE
           data_t_agg = NA
           data_out_harvest = NA
        END IF
        !write(*,*)'ST2'
    


        !aggregate
        DO v=1,NV-2
        ! all hyde default classes array 
        !dummya1[,,v] <- dummya
          IF(compute_next_year) THEN
            dummya = data_in_new(v,:,:)
          ELSE
            dummya = data_in(v,:,:)
          END IF

        ! aggregate Hyde landcover types 
          data_out(aggmap(v),:,:) = data_out(aggmap(v),:,:) + dummya
        END DO

        ! change to percent      
        data_out = data_out * 100.0      

        esaarray(11:NE2,:,:) = data_out(1:6,:,:)
        esamask(11:NE2,:,:) = esamask(1:6,:,:)
    
        !IF( (t == SINDEX) .AND. (debug .EQV. .TRUE.) ) THEN
        IF( debug .EQV. .TRUE. ) THEN
          WRITE(*,*)'JH0' ! Assumes default "ADVANCE='yes'".
          IF(num_land.GT.0.0d0)THEN
            DO v=1,NE2
            !WRITE(*,FMT='(F9.4)', ADVANCE='no') esaarray(v,3,3)/100.0 
            WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(esaarray(v,3:4,3:4), & 
               esamask(v,3:4,3:4))/num_land/100.0 
            END DO
            WRITE(*,*) ! Assumes default "ADVANCE='yes'".
          ELSE
            WRITE(*,*) 'no valid data for esaarray' ! Assumes default "ADVANCE='yes'".
          ENDIF

        ENDIF

        ! deal with forest or grass cover where no forest or grass cover existed in ESA
        ! process, returns an array of PFT, lon, lat  
        ! R code: data_out_SDGVM   <- apply(esaarray,c(1,2),join_hyde)
        !DO lon=1,DXY !PCM skip lon,lat from 1,2
        ! DO lat=1,DXY
        DO lon=3,DXY
         DO lat=3,DXY
           data_out_SDGVM(year_index,:,lon,lat) = JOIN_HYDE(esaarray(:,lon,lat))
         END DO
        END DO

        !data_out_SDGVM(year_index,:,3,3) = JOIN_HYDE(esaarray(:,3,3))

      
        WHERE(ABS(data_out_SDGVM(year_index,1,:,:)) <= 100.00)
          mask_SDGVM = .TRUE.
        ELSEWHERE
          mask_SDGVM = .FALSE.
        END WHERE


        ! convert missing values to SDGVM missing value
        WHERE (data_out_SDGVM == NA)
          data_out_SDGVM = 255.0 
        END WHERE 
        num_land_SDGVM = COUNT( mask_SDGVM(3:4,3:4) .EQV. .TRUE.)

        ! deal with forest or grass cover where no forest or grass cover existed in ESA
        CALL F_LAT_ASSIGNPFT(data_out_SDGVM(year_index,:,:,:), &
                           NS,NY,DXY,                          &
                           minii,minjj,maxii,maxjj,lat_offset,Y0,debug)
  
        !IF( (t == SINDEX) .AND. (debug .EQV. .TRUE.) ) THEN
        IF( debug .EQV. .TRUE. ) THEN
          WRITE(*,*)'JH1' ! Assumes default "ADVANCE='yes'".
          IF(num_land_SDGVM.GT.0) THEN
           DO v=1,NS
            WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_out_SDGVM(year_index, &
              v,3:4,3:4),mask_SDGVM(3:4,3:4))/num_land_SDGVM/100.0 
           END DO
           WRITE(*,*) ! Assumes default "ADVANCE='yes'".
          ELSE
           WRITE(*,*) 'no valid SDGVM PFT data' ! Assumes default "ADVANCE='yes'".
          ENDIF
            
        ENDIF

        IF(get_transitions) THEN
          data_out_AggHydeTransitions(year_index,:,:,:,:) = data_t_agg 
          data_out_AggHarvest(year_index,:,:,:) = data_out_harvest
        ! convert missing values to SDGVM missing value
          WHERE (data_out_AggHydeTransitions == NA)
            data_out_AggHydeTransitions = 255.0 
          END WHERE 
          WHERE (data_out_AggHarvest == NA)
            data_out_AggHarvest = 255.0 
          END WHERE 
        ENDIF

        IF( t == SINDEX ) THEN
          IF(debug) THEN
            PRINT *, 'SDGVM num_land=',num_land_SDGVM,'num_tot=',maxii*maxjj
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
              DO v=1,NS
                 WRITE(*,FMT='(A9)', ADVANCE='no') varname_sdgvm(v)
              END DO
            ELSE IF(print_type=='esa') THEN
              DO v=1,NE
                 WRITE(*,FMT='(A9)', ADVANCE='no') varname_esa_pfts(v)
              END DO
            END IF
            WRITE(*,*) ! Assumes default "ADVANCE='yes'".
          END IF
        END IF

        IF( (MOD(t-SINDEX,DT) == 0) .AND. (debug .EQV. .TRUE.) ) THEN 
           IF(compute_next_year) THEN
               WRITE(*,FMT='(I5)', ADVANCE='no') t+849+2 !t=1 is the year 850
           ELSE
               WRITE(*,FMT='(I5)', ADVANCE='no') t+849 !t=1 is the year 850
           END IF
           IF(print_type=='unagg') THEN
             IF(compute_next_year) THEN
               DO v=1,NV
!                  WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_in_new(v,:,:),mask)/num_land
!                  WRITE(*,FMT='(F9.4)', ADVANCE='no') data_in_new(v,3,3)
                  WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_in_new(v,3:4,3:4),mask(3:4,3:4))/num_land
               END DO
             ELSE
               DO v=1,NV
!                  WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_in(v,:,:),mask)/num_land
!                  WRITE(*,FMT='(F9.4)', ADVANCE='no') data_in(v,3,3)
                  WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_in(v,3:4,3:4),mask(3:4,3:4))/num_land
               END DO
             END IF
           ELSE IF(print_type=='agg') THEN
             DO v=1,NV2
!                WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_out(v,:,:)/100.0,mask)/num_land 
!                WRITE(*,FMT='(F9.4)', ADVANCE='no') data_out(v,3,3)/100.0
                WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_out(v,3:4,3:4)/100.0,mask(3:4,3:4))/num_land 
             END DO
           ELSE IF(print_type=='sdgvm') THEN 
             DO v=1,NS
!                WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_out_SDGVM(year_index,v,:,:)/100.0,mask_SDGVM)/num_land 
!                WRITE(*,FMT='(F9.4)', ADVANCE='no') data_out_SDGVM(year_index,v,3,3)/100.0 
                WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_out_SDGVM(year_index,  &
                      v,3:4,3:4)/100.0,mask_SDGVM(3:4,3:4))/num_land_SDGVM 
             END DO
           ELSE IF(print_type=='esa') THEN 
             DO v=1,NE
!                WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_out_SDGVM(year_index,v,:,:)/100.0,mask_SDGVM)/num_land 
!                WRITE(*,FMT='(F9.4)', ADVANCE='no') esaarray(v,3,3)/100.0 
                WRITE(*,FMT='(F9.4)', ADVANCE='no') SUM(data_out_SDGVM(year_index,v,3:4,3:4)/100.0, &
                       esamask(v,3:4,3:4))/num_land 
             END DO
           END IF
           WRITE(*,*) ! Assumes default "ADVANCE='yes'".
        END IF

      END DO

      ! Close the file, freeing all resources.
      CALL CHECK( NF90_CLOSE(ncid) )
      CALL CHECK( NF90_CLOSE(ncid_t) )

     !PRINT *,"*** SUCCESS reading file ", pname(1:blank(pname)), "! "
     !CALL FLUSH()

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
      !ORIG:
      ! takes a 15 element vector, 1:10 - SDGVM PFTs from ESA data, 11:15 HYDE aggregated land cover 
      ! and combines into a 10 element vector of SDGVM PFTs
       
      ! SDGVM PFTs: 1 BARE, 2 Ev_Bl, 3 Dc_Bl, 4 Ev_Nl, 5 Dc_Nl, 6 Shrub, 7 C3, 8 C4, 9 C3crop, 10 C4crop
      ! HYDE aggregated land cover: 11 forest, 12 non-forest, 13 C3 crop, 14 C4 crop, 15 pasture   
      
      ! first forest area and potential forest area from HYDE is used to adjust forest area in ESA
      ! then C3 and C4 croplands in HYDE are used to adjust ESA
      ! C3 and C4 grassland area is then adjusted to account for changes in pasture and natural grassland cover
      ! rangelands in hyde and shrubs in esa are not explicitly included in the calculation of change 
      ! though a proportion of them could be accounted for implicitly if rangelands habour some forest vegetation

      !NEW:
      ! takes a 17 element vector, 1:10 - ESA PFTs, 11:17 HYDE aggregated land cover 
      ! and combines into a 17 element vector of SDGVM PFT transitions

      ! ESA   PFTs: 1 BARE, 2 Ev_Bl, 3 Dc_Bl, 4 Ev_Nl, 5 Dc_Nl, 6 Shrub, 7 C3, 8 C4, 9 C3crop, 10 C4crop,
      ! HYDE aggregated land cover: 11 primf, 12 primn, 13 secdf, 14 secdn, 15 C3 crop, 16 C4 crop, 17 pasture

      ! SDGVM PFTs: 1 BARE, 2 Ev_Bp, 3 Dc_Bp, 4 Ev_Np, 5 Dc_Np, 6 Shrup, 7 C3p, 8 C4p, 9 C3crop, 10 C4crop,
      !                    11 Ev_Bs,12 Dc_Bs,13 Ev_Ns,14 Dc_Ns,15 Shrus, 16 C3s, 17 C4s
      
      ! first forest area and potential forest area from HYDE is used to adjust forest area in ESA
      ! then C3 and C4 croplands in HYDE are used to adjust ESA
      ! C3 and C4 grassland area is then adjusted to account for changes in pasture and natural grassland cover
      ! rangelands in hyde and shrubs in esa are not explicitly included in the calculation of change 
      ! though a proportion of them could be accounted for implicitly if rangelands habour some forest vegetation

      !NEW2:
      ! takes a 17 element vector, 1:10 - ESA PFTs, 11:16 HYDE aggregated land cover 
      ! and combines into a 17 element vector of SDGVM PFT transitions

      ! ESA   PFTs: 1 BARE, 2 Ev_Bl, 3 Dc_Bl, 4 Ev_Nl, 5 Dc_Nl, 6 Shrub, 7 C3, 8 C4, 9 C3crop, 10 C4crop,
      ! HYDE aggregated land cover: 11 primf, 12 primn, 13 secdf, 14 secdn, 15 C3 crop, 16 C4 crop
      !!  example:
      !!            0.7800   0.0100   0.0300   0.0100   0.0000   0.0400   0.0400 0.0000 0.0000   0.0000   
      !!                              0.7420   0.0000   0.0000       0.0000   0.0000   0.0000

      ! SDGVM PFTs: 1 BARE, 2 Ev_Bp, 3 Dc_Bp, 4 Ev_Np, 5 Dc_Np, 6 Shrup, 7 C3p, 8 C4p, 9 C3crop, 10 C4crop,
      !                    11 Ev_Bs,12 Dc_Bs,13 Ev_Ns,14 Dc_Ns,15 Shrus, 16 C3s, 17 C4s
      !!  example:
      !!            0.7800   0.0824   0.2473   0.0824   0.0000   0.3298   0.0000   0.0000  0.0000   0.0000
      !!                     0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000
      
      ! first forest area and potential forest area from HYDE is used to adjust forest area in ESA
      ! then C3 and C4 croplands in HYDE are used to adjust ESA
      ! C3 and C4 grassland area is then adjusted to account for changes in natural grassland cover
      ! pasture and rangelands in hyde and shrubs in esa are not explicitly included in the calculation of change 
      ! though a proportion of them could be accounted for implicitly if rangelands habour some forest vegetation
      IMPLICIT NONE
      INTEGER, PARAMETER :: NE=10, NH=6, NS=17, DF=9, DN=9
      REAL*8, INTENT(IN) :: v(NE+NH)
      REAL*8, DIMENSION(NS) :: result_join_hyde
      REAL*8 :: ov(NS), nofcov, nogcov, for_esa, diff, esahyde_cov, hyde_cov, past_esa
      !ORIG:
      !INTEGER, PARAMETER :: NE=10, NH=5, NS=10 
      !NEW:
      !INTEGER, PARAMETER :: NE=10, NH=7, NS=15, DF=9
      !NE=NUMESA, NH=NUM_Agg_Hyde, NS=NUM_SDGVM,
      ! DF = Difference in Index from first primary tree FT to first secondary ! tree FT 
      !REAL, INTENT(IN) :: v(NE+NH)
      !REAL, DIMENSION(NS) :: result_join_hyde
      !REAL :: ov(NS), nofcov, nogcov, for_esa, diff, esahyde_cov, hyde_cov, past_esa
      REAL :: ovp(NS), ovs(NS)
      INTEGER :: fid_min,fid_max,f
      INTEGER :: nid_min,nid_max,n
      INTEGER :: fsec,nsec
      REAL*8, PARAMETER :: NA= -999.0
      REAL for_hyde, primf_frac_hyde, secdf_frac_hyde 
      REAL past_hyde, primn_frac_hyde, secdn_frac_hyde 
      LOGICAL :: do_pasture
      ! PRINT(v)
      
      ! subscripts for ESA forest PFTs
      fid_min = 2
      fid_max = 6
        
      ! subscripts for ESA non-forest, non-crop PFTs
      nid_min = 7
      nid_max = 8
        
      !if(is.finite(sum(v(1:10))) & not is.finite(sum(v(11:15))) ) then 
      IF((SUM(v(1:NE)) > 0.0) .AND. (SUM(v(1:NE)) < 2000.0)  .AND. &
        .NOT.((SUM(v((NE+1):(NE+NH))) > 0.0) .AND. (SUM(v((NE+1):(NE+NH))) < 2000.0)) ) THEN 
        ! when ESA has landcover but hyde not
        ! arises from different land masks
        IF(NS.EQ.NE) THEN
          result_join_hyde(1:NE) = v(1:NE)
        ELSE IF(NS.GT.NE) THEN
          result_join_hyde(1:NE) = v(1:NE)
          result_join_hyde((NE+1):NS) = 0.0  !PCM assume secondary vegetation = 0, for now
        ENDIF
        
      !else if(is.finite(sum(v))) then
      ELSE IF((SUM(v) > 0.0) .AND. (SUM(v) < 2000.0) ) THEN 
        
        ! create output vector 
        ! print(v)
        ov(1:NE) = v(1:NE)
        ov(NE+1:(NE+NH)) = 0.0d0
        
        ! stores for increase in ESA but when either no forest or grass cover in ESA 
        nofcov  = 0.0
        nogcov  = 0.0
        
        ! change forest cover
        ! ESA forest cover
        for_esa = SUM(v(fid_min:fid_max))
        ! Hyde forest cover (primary + secondary)
        for_hyde = v(11) + v(13)
        IF(for_hyde .GT. 0) THEN
          primf_frac_hyde = v(11)/for_hyde
          secdf_frac_hyde = v(13)/for_hyde
        ELSE
          primf_frac_hyde = 0.0
          secdf_frac_hyde = 0.0 
        ENDIF
        IF(debug)THEN
          write(*,*) 'JI0a',v
          write(*,*) 'JI0b',ov
          write(*,*) 'JI1',for_esa,for_hyde,primf_frac_hyde,secdf_frac_hyde
        ENDIF
        !PCM initialize primary and secondary fractions of output field 
        IF(NS.GT.NE) THEN
          ovp(fid_min:fid_max) = ov(fid_min:fid_max)
          ovs((fid_min+DF):(fid_max+DF)) = ov(fid_min:fid_max)
          ov(fid_min:fid_max) = ovp(fid_min:fid_max)
          ov((fid_min+DF):(fid_max+DF)) = ovs((fid_min+DF):(fid_max+DF))
        ENDIF

        ! If forest cover is greater in hyde
        IF(for_hyde > for_esa) THEN 
          diff   = for_hyde  - for_esa
          IF(debug)THEN
            write(*,*) 'JI2',for_esa,for_hyde,diff
          ENDIF
          IF(for_esa>0) THEN 
            ! increase forest cover proportionally
            DO f=fid_min,fid_max
               ov(f)  = primf_frac_hyde * (ov(f) + diff * (v(f) / for_esa))
               fsec = f + DF !secondary tree FT's
               ov(fsec)  = secdf_frac_hyde * (ov(fsec) + diff * (v(f) / for_esa))
            END DO
          ELSE 
            ! no forest cover in ESA - forest cover stored for post-processing
            !write(paste('no ESA forest','BARE:',ov[1],'C3 CROP:',ov[9],'C4 CROP:',ov[10]),'1nofor_error.txt',append=T)
            print *,'no ESA forest: ','BARE:',ov(1), &
                    'C3 CROP:',ov(9),'C4 CROP:',ov(10)
            nofcov = diff
          END IF 
        ELSE IF(for_esa>for_hyde) THEN 
          ! forest cover greater in ESA
          diff   = for_esa - for_hyde 
          IF(debug)THEN
            write(*,*) 'JI3',for_esa,for_hyde,diff
          ENDIF
          ! decrease forest cover proportionally
          DO f=fid_min,fid_max
             ov(f)  = primf_frac_hyde * (ov(f) - diff * (v(f) / for_esa))
             fsec = f + DF !secondary tree FT's
             ov(fsec)  = secdf_frac_hyde * (ov(fsec) - diff * (v(f) / for_esa))
          END DO
        END IF 
        IF(debug)THEN
          write(*,*) 'JI4b',ov
        ENDIF
        
        
        ! increase or decrease crop cover 
        ! C3 crops
        IF(v(9) .NE. v(15))  ov(9)  = v(15)
        
        ! C4 crops
        IF(v(10) .NE. v(16)) ov(10) = v(16)
        
        ! increase or decrease grass (pasture) cover 
        ! total cover so far in hyde ESA combined dataset - accounts for forest cover as yet unassigned to a PFT
        !PCM don't skip BARE soil anymore ! PCM: skip BARE soil (ov(1)) for now in esahyde_cov calc
        esahyde_cov = SUM(ov(1:(NE+NH))) + nofcov
        
        ! total cover in hyde
        hyde_cov    = SUM(v(11:(NE+NH)))
        
        ! change grassland cover
        ! ESA grasslands (potential pasture)
        past_esa    = v(7) + v(8)
        ! Hyde non-forest cover (primary + secondary + pasture + rangeland)
        past_hyde = v(12) + v(14) 
        IF(past_hyde .GT. 0) THEN
          primn_frac_hyde = v(12)/past_hyde
          secdn_frac_hyde = v(14)/past_hyde
        ELSE
          primn_frac_hyde = 0.0
          secdn_frac_hyde = 0.0 
        ENDIF


        IF(debug)THEN
          write(*,*) 'JJ0a',v
          write(*,*) 'JJ0b',ov
          write(*,*) 'JJ1',past_hyde,primn_frac_hyde,secdn_frac_hyde
        ENDIF
        !PCM initialize primary and secondary fractions of output field 
        IF(NS.GT.NE) THEN
          ovp(nid_min:nid_max) = ov(nid_min:nid_max)
          ovs((nid_min+DN):(nid_max+DN)) = ov(nid_min:nid_max)
          ov(nid_min:nid_max) = ovp(nid_min:nid_max)
          ov((nid_min+DN):(nid_max+DN)) = ovs((nid_min+DN):(nid_max+DN))
        ENDIF

        
        IF(debug)THEN
          write(*,*) 'JJ2',esahyde_cov,hyde_cov
        ENDIF
        ! if new cover is different from hyde cover
        IF(esahyde_cov .NE. hyde_cov) THEN 
          diff   = hyde_cov - esahyde_cov
          IF(debug)THEN
            write(*,*) 'JJ3',diff,past_esa
          ENDIF
          IF(diff<0.0) THEN 
            IF(ABS(diff)<past_esa) THEN 
              IF(debug)THEN
                write(*,*) 'JJ4'
              ENDIF
              
              ! change C3/C4 grass cover proportionally
              DO n=nid_min,nid_max
                ov(n)  = primn_frac_hyde * (ov(n) + diff * (v(n) / past_esa))
                nsec = n + DN !secondary non-forest
                ov(nsec)  = secdn_frac_hyde * (ov(nsec) +  diff * (v(n) / past_esa))
              ENDDO
              
            ELSE
              
              ! ESA pasture needs reduction but reduction is more than pasture cover
              ov(7) = 0.0
              ov(8) = 0.0
              ov(7+DN) = 0.0
              ov(8+DN) = 0.0
              
              ! assume remaining reduction from bare ground
              diff = abs(diff) - past_esa
              IF(debug)THEN
                write(*,*) 'JJ5',diff,past_esa
              ENDIF
              IF(ov(1)>=diff) THEN
                 ov(1) = ov(1) - diff
              ELSE
                 ov(1) = 0.0
              END IF
            END IF   
            
          ELSE IF(diff>0.0) THEN 
            IF(past_esa>0.0) THEN 
             
              IF(debug)THEN  
                write(*,*) 'JJ6',diff,past_esa,v(nid_min),ov(nid_min)
              ENDIF
              ! change C3/C4 grass cover proportionally
              DO n=nid_min,nid_max
                ov(n)  = primn_frac_hyde * (ov(n) + diff * (v(n) / past_esa))
                nsec = n + DN !secondary non-forest
                ov(nsec)  = secdn_frac_hyde * (ov(nsec) +  diff * (v(n) / past_esa))
              ENDDO
              IF(debug)THEN
                write(*,*) 'JJ7',diff,past_esa,v(nid_min),ov(nid_min)
              ENDIF
              
            ELSE 
              
              IF(debug)THEN
                write(*,*) 'JJ8',diff
              ENDIF
              ! no grass cover in ESA - grass cover stored for post-processing
              !write(paste('no ESA grass','BARE:',ov[1],'C3 CROP:',ov[9],'C4 CROP:',ov[10]),'1nopast_error.txt',append=T)
              !PRINT *,'no ESA grass: ','BARE:',ov(1),'C3 CROP:', &
              !            ov(9),'C4 CROP:',ov(10)
              nogcov = diff
            END IF 
          END IF 
        END IF 
        
        ! test for cover sum >100%
        !IF( SUM(ov(1:NS))-100.0 > 1e-3 ) PRINT *,'Error,cov:',SUM(ov(1:NS)),hyde_cov,'. Bare:',v(1) !'1cov_error.txt'
        
        ! assign as yet unassigned forest cover as negative value to DcBp
        !IF(nofcov>0) ov(3)  = -nofcov
        ! assign as yet unassigned grass cover as negative value to C3p
        !IF(nogcov>0) ov(7)  = -nogcov

        IF(nofcov>0) THEN !PCM3 
              ov(2) = 0.0
              ov(3) = -primf_frac_hyde*nofcov !set to Dc_Bp
              ov(4:6) = 0.0
              ov(2+DF) = 0.0
              ov(3+DF) = -secdf_frac_hyde*nofcov !set to Dc_Bs
              ov(4+DF:6+DF) = 0.0
        END IF

        IF(nogcov>0) THEN !PCM3
              ov(7) = -primn_frac_hyde*nogcov !set to C3p 
              ov(8) = 0.0
              ov(7+DN) = -secdn_frac_hyde*nogcov !set to C3s 
              ov(8+DN) = 0.0
        END IF


        
        result_join_hyde = ov(1:NS)

        WHERE (isNAN(result_join_hyde)) !extra check for NaNs
               result_join_hyde = NA
        ENDWHERE

       ELSE 
        result_join_hyde = (/ NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA /) 
       END IF 
      END FUNCTION JOIN_HYDE

      SUBROUTINE F_LAT_ASSIGNPFT(m,NS,NY,DXY,minii,minjj,maxii,maxjj, &
                                  loff,Y0,debug)
      ! Assigns forest or grass pfts to ESA HYDE according to latitude when no forest or grass PFTs exist in ESA but do in HYDE 
  
      ! this function assigns forest and grassland cover to an appropriate PFT in the combined HYDE ESA dataset
      !     when there was no forest or grass cover in the ESA data.
      !     in this case the forest cover is assigned as a negative value to primary deciduous broadleaved (DcBl) PFT m(3,j) 
      !     in this case the forest cover is assigned as a negative value to secondary deciduous broadleaved (DcBs) PFT m(11,j) 
      !     in this case the grass  cover is assigned as a negative value to C3 primary grass (C3p) PFT m(7,j)
      !     in this case the grass  cover is assigned as a negative value to C3 secondary grass (C3p) PFT m(16,j)
      ! if in temperate latitudes this negative value is simply switched to positive
      ! if in tropical latitudes this negative value is switched to positive and assigned to evergreen broadleaved PFT or C4 grass 
      IMPLICIT NONE
     
      INTEGER :: NS !number of PFTs
      INTEGER :: DXY !dimension of array 
      INTEGER :: NY !dimension of latitude loff array 
      ! m is a pft x DXY x DXY matrix
      REAL*8, DIMENSION(NS,DXY,DXY) :: m
      ! loff is the change in subscript by lat - 0 for temperate (DcBp & C3p & DcBs & C3s), 1 for tropical (EvBp & C4p & EvBs & C4s)
      INTEGER, DIMENSION(NY) :: loff 
      INTEGER :: i,j !longitude & latitude dummy indices
      INTEGER :: Y0 !latitude index
      INTEGER :: minii,minjj,maxii,maxjj ! longitude, latitude indices
      INTEGER, DIMENSION(2) :: indx !'permuted DcBp & unpermuted C3p' or 'permuted DcBs & unpermuted C3s'
      INTEGER :: k,kk
      REAL*8 :: tmp3
      LOGICAL :: debug !used to print out more debugging info
    
      
      DO i=minii,maxii !longitude loop 
       DO j=minjj,maxjj !latitude loop
      ! if there was no forest or grass cover in the ESA data 
        IF(any(m(:,i,j)<0)) THEN

          indx = [2,7] !permuted DcBp & unpermuted C3p
          ! switch DcBp and EvBp PFT to allow loff to work for both forest and grass
          tmp3 = m(3,i,j)
          m(3,i,j) =  m(2,i,j) 
          m(2,i,j) = tmp3 

          DO k = 1,2
            kk = indx(k)
            IF(m(kk,i,j)<0) THEN
            ! take absolute value of negative covers and assign to appropriate PFT
              m(kk+loff(Y0+j-1),i,j) = ABS(m(kk,i,j))
            ! if tropical PFT is appropriate, zero negative cover in temperate PFT
              IF(loff(Y0+j-1) .NE. 0) m(kk,i,j) = 0.0
            END IF
          END DO

          ! switch DcBp and EvBp PFT back to original placement in vector
          tmp3 = m(3,i,j)
          m(3,i,j) =  m(2,i,j) 
          m(2,i,j) = tmp3 


          indx = [11,16] !permuted DcBs & unpermuted C3s
          ! switch DcBs and EvBs PFT to allow loff to work for both forest and grass
          tmp3 = m(12,i,j)
          m(12,i,j) =  m(11,i,j) 
          m(11,i,j) = tmp3 

          DO k = 1,2
            kk = indx(k)
            IF(m(kk,i,j)<0) THEN
            ! take absolute value of negative covers and assign to appropriate PFT
              m(kk+loff(Y0+j-1),i,j) = ABS(m(kk,i,j))
            ! if tropical PFT is appropriate, zero negative cover in temperate PFT
              IF(loff(Y0+j-1) .NE. 0) m(kk,i,j) = 0.0
            END IF
          END DO

          ! switch DcBs and EvBs PFT back to original placement in vector
          tmp3 = m(12,i,j)
          m(12,i,j) =  m(11,i,j) 
          m(11,i,j) = tmp3 

        END IF 
       END DO
      END DO

      END SUBROUTINE F_LAT_ASSIGNPFT

     END SUBROUTINE states_convertSDGVM_func
     END MODULE FUNCTIONS_CLU
