*----------------------------------------------------------------------*
*                          SUBROUTINE READCO2                          *
*                          ******************                          *
* This routine reads in co2 values from the file co2 which in found in *
* the stinput directory. If only one record exists in the file then    *
* the co2 value in this record is used throughout the simulation.      *
*----------------------------------------------------------------------*
      SUBROUTINE READCO2(stco2,co2,daily_co2,yr0a,yrfa,year0set,spinl,
     &nyears)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8    co2(maxyrs,12,31),ca,co2const
      INTEGER   norecs,year,const,blank,kode,daily_co2,day,mnth
      INTEGER   yr0a,yrfa,prev_year,spinl,dummy,nyears
      CHARACTER stco2*1000
      LOGICAL   co2spin,year0set

      OPEN(98,FILE=stco2,STATUS='OLD',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Co2 file does not exist.'
        WRITE(*,'('' "'',A,''"'')') stco2(1:blank(stco2))
        STOP
      ENDIF

      norecs = 0
      const  = 0

10    CONTINUE
!      print*, daily_co2

      IF(daily_co2.eq.1) THEN 

        READ(98,*,end=99) year,mnth,day,ca 

        IF ((year.eq.yr0a).and.(mnth.eq.1).and.(day.eq.1)) prev_year =
     &year-1
        !print*, prev_year, norecs      

        IF (year0set) THEN
          IF ((year.GE.yr0a).AND.(year.LE.yrfa)) THEN
            ! standard behaviour
            ! assign years of CO2 data from 1 to nyears in co2 array
            IF(prev_year.EQ.(year-1)) norecs = norecs + 1 
            co2(norecs,mnth,day) = ca 
        
          ELSE IF ((norecs.eq.0).and.(yr0a.LT.year)) THEN
            ! when yr0a is less than start year of CO2 record,
            ! fill co2 array with first day of co2 data up until and 
            ! including the year in which the data series begins
            ! first day instead of whole year is not 100% satisfactory, but unlikely this option will be used
            ! alternative behaviour copuld be to abort
            do dummy=yr0a,year
              norecs = norecs + 1
              co2(norecs,:,:) = ca
            enddo
          ENDIF
        
        ELSE IF(norecs.LT.spinl) THEN
          ! in the case where only a spin is requested and co2const < 0 
          IF(prev_year.EQ.(year-1)) norecs = norecs + 1 
          co2(norecs,mnth,day) = ca 
        ENDIF

        prev_year = year 

      ELSE 

        READ(98,*,end=99) year,ca

        IF (year0set) THEN
          IF ((year.GE.yr0a).AND.(year.LE.yrfa)) THEN
            ! standard behaviour
            ! assign years of CO2 data from 1 to nyears in co2 array
            norecs = norecs + 1
            co2(norecs,:,:) = ca
        
          ELSE IF ((norecs.eq.0).and.(yr0a.LT.year)) THEN
            ! when yr0a is less than start year of CO2 record,
            ! fill co2 array with first year of co2 data up until and 
            ! including the year in which the data series begins
            ! alternative behaviour copuld be to abort
            do dummy=yr0a,year
              norecs = norecs + 1
              co2(norecs,:,:) = ca
            enddo
          ENDIF
        
        ELSE IF(norecs.LT.spinl) THEN
          ! in the case where only a spin is requested and co2const < 0 
          norecs = norecs + 1
          co2(norecs,:,:) = ca
        ENDIF

      ENDIF 

      const = const + 1
      GOTO 10
99    CONTINUE

      !print*, spinl, yr0, yr0a, yrf, yrfa
      CLOSE(98)
      
      ! if only a single value supplied in the file
      IF (const.EQ.1) THEN
     
        DO year=yr0a,yrfa
          co2(year-yr0a+1,:,:) = ca
        ENDDO
     
      ELSE IF ((norecs.LT.nyears).AND.(norecs.LT.(nyears-spinl))) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Only',norecs,'years of data in CO2 file, 
     &but requested simulation is',nyears,
     &'. Need additional years of CO2 data.'
        STOP
     
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                         SUBROUTINE LANDUSE1                          *
*----------------------------------------------------------------------*
      SUBROUTINE LANDUSE1(luse,fire,harvest,yr0a,yrfa,year0set,spinl)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      INTEGER yr0a,yrfa,year,rep,luse(maxyrs),i,use,spinl,yr0a_spin
      LOGICAL fire(maxyrs),harvest(maxyrs),year0set

      luse(:) = 1000
      
      i = 1
      rep = 0
10    CONTINUE
        READ(98,*,end=20) year,use
      
        IF (use.GT.0) rep = use
        IF (year0set) THEN
          IF (year-yr0a+1.GT.0)  luse(year-yr0a+1) = use
        ELSE
          IF (i.EQ.1) yr0a_spin = year
          luse(year-yr0a_spin+1) = use   
        ENDIF

        i = i+1
      GOTO 10
20    CONTINUE

      if(rep.LE.0) then
        print*, 'ERROR:: Land use mapping equal to or below 0'
        print*, 'you must specifiy at least one year with a cover type'
        STOP
      endif

      DO i=1,maxyrs
        If ((luse(i).GT.0).AND.(luse(i).LT.1000)) rep = luse(i)
        If (luse(i).LT.0) fire(i)    = .TRUE.
        If (luse(i).EQ.0) harvest(i) = .TRUE.
        luse(i) = rep
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE EX_CLIM                          *
*                          ******************                          *
*                                                                      *
* Extract climate data from the climate database for the nearest site  *
* to lat,lon, replace lat,lon with the nearest cell from the database. *
*                                                                      *
*            UNIX                DOS                                   *
*                                                                      *
*            ii(4)          PCM  jj(5)   for beginning of binary file  *
*                                        records                       *
*            ii(7)          PCM  jj(8)   for beginning of binary file  *
*                                        records                       *
*            recl = 728          recl = 730     for binary climate     *
*            recl = 577          recl = 578     for text map           *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE EX_CLIM(stinput,lat,lon,xlatf,xlatres,xlatresn,xlon0,
     &xlonres,xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,xwndv,isite,
     &year0,yearf,siteno,du,swrv,read_par)
*----------------------------------------------------------------------*
      REAL*8 lat,lon,xlon0,xlatf,xlatres,xlonres,ans(12)
      REAL*8  xtmpv(500,12,31),xhumv(500,12,31),xprcv(500,12,31) !PCM
      REAL*8  swrv(500,12,31),xwndv(500,12,31) !PCM
      INTEGER year,year0,yearf,nrec,ncol,ans2(1000),siteno,i,du
      INTEGER nyears,yr0,mnth,day,isite,yrf,blank,recl1,recl2
      INTEGER no_days
      INTEGER recl3                      !PCM
      INTEGER xlatresn,xlonresn,fno,read_par
      !CHARACTER ii(4),jj(5),fname4*1000 !PCM
      CHARACTER ii*7,jj*8,fname4*1000,fname5*1000 !PCM
      CHARACTER fname1*1000,fname2*1000,fname3*1000,stinput*1000
      CHARACTER num*15 !PCM
      INTEGER*2 tmpv(500,12,31),humv(500,12,31),prcv(500,12,31) !PCM
      INTEGER*2 wndv(500,12,31) !PCM
      REAL*8 TMP_MULT,PRC_MULT,HUM_MULT,PRC_MULT1,WND_MULT

      IF (du.eq.1) THEN
        !recl1 = 730  !PCM
        recl1 = 7 + 4*360 !PCM 
        recl2 = 6*xlonresn
        recl3 = 7 + 7*360 !PCM
      ELSE
        !recl1 = 728  !PCM
        recl1 = 7 + 4*360 + 1 !PCM
        recl2 = 6*xlonresn + 1
        recl3 = 7 +  7*360 + 1 !PCM
      ENDIF

      nyears = yearf - year0 + 1

      nrec = int((xlatf - lat)/xlatres + 1.0d0)
      ncol = int((lon - xlon0)/xlonres + 1.0d0)
C      write(*,*) 'nrec,ncol,xlatf,lat,xlatres,lon,xlon0,xlonres'
C      write(*,*) nrec,ncol,xlatf,lat,xlatres,lon,xlon0,xlonres

      IF ((nrec.LE.xlatresn).AND.(ncol.LE.xlonresn)) THEN

      lat = xlatf - xlatres/2.0 - real(nrec - 1)*xlatres
      lon = xlon0 + xlonres/2.0 + real(ncol - 1)*xlonres

      fno = 90
C      write(*,*) 'nrec,ncol,xlatf,lat,xlatres,lon,xlon0,xlonres in IF'
C      write(*,*) nrec,ncol,xlatf,lat,xlatres,lon,xlon0,xlonres

      OPEN(fno+1,file=stinput(1:blank(stinput))//'/maskmap.dat',
     &ACCESS='DIRECT',RECL=recl2,FORM='formatted',STATUS='OLD')

      !READ(fno+1,'(96i6)',REC=nrec) (ans2(i),i=1,xlonresn) !PCM
      READ(fno+1,'(305i6)',REC=nrec) (ans2(i),i=1,xlonresn) !PCM
      CLOSE(fno+1)
      siteno = ans2(ncol)
C      write(*,*) 'siteno after maskmap.dat read'
C      write(*,*) siteno

      isite = 0
      IF (siteno.GT.0) THEN

        isite = 1

C PCM        WRITE(num,'(i3.3)') (siteno-1)/100
C       PCM I kept the trailing zeroes, and I had no leading zeroes 
C       PCM added trailing '.dat'
        WRITE(num,'(i0,a4)')  100*((siteno-1)/100),'.dat'
        WRITE(fname1,'(100a)') (stinput(i:i),i=1,blank(stinput)),
     &'/tmp_',num
        WRITE(fname2,'(100a)') (stinput(i:i),i=1,blank(stinput)),
     &'/hum_',num
        WRITE(fname3,'(100a)') (stinput(i:i),i=1,blank(stinput)),
     &'/prc_',num
        WRITE(fname4,'(100a)') (stinput(i:i),i=1,blank(stinput)),
     &'/swr_',num
        WRITE(fname5,'(100a)') (stinput(i:i),i=1,blank(stinput)),
     &'/wnd_',num

C        write(*,*) 'siteno before mod'
C        write(*,*) siteno
        siteno = mod(siteno-1,100) + 1

        OPEN(fno+1,file=fname1,access='direct',recl=recl1,
     &form='formatted',status='old')                         !PCM
        OPEN(fno+2,file=fname2,access='direct',recl=recl1,
     &form='formatted',status='old')                         !PCM
        OPEN(fno+3,file=fname3,access='direct',recl=recl1,
     &form='formatted',status='old')                         !PCM
C         OPEN(fno+4,file=fname4,access='direct',recl=recl1, !PCM
C     &form='unformatted',status='old')                      !PCM
        OPEN(fno+4,file=fname4,access='direct',recl=recl3,   !PCM
     &form='formatted',status='old')                         !PCM
        OPEN(fno+5,file=fname5,access='direct',recl=recl3,   !PCM
     &form='formatted',status='old')                         !PCM

        IF (du.eq.1) THEN
        DO year=yr0,yrf
            READ(fno+1,1001,                     !PCM
     & REC=(siteno-1)*nyears+year-year0+1) jj,
     &((tmpv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
            READ(fno+2,1001,                     !PCM
     & REC=(siteno-1)*nyears+year-year0+1) jj,
     &((humv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
            READ(fno+3,1001,                     !PCM
     & REC=(siteno-1)*nyears+year-year0+1) jj,
     &((prcv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
            if(read_par.eq.1) 
     & READ(fno+4,1002,                          !PCM
     & REC=(siteno-1)*nyears+year-year0+1)
     &jj, ((swrv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
            READ(fno+5,1001,                     !PCM
     & REC=(siteno-1)*nyears+year-year0+1) jj,
     &((wndv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
C            DO mnth=1,12
C              DO day=1,30
C                prcv(year-yr0+1,mnth,day) = 
C     &int(real(prcv(year-yr0+1,mnth,day))/100.0 )   !PCM
C     &int(real(prcv(year-yr0+1,mnth,day))/10.0 + 0.5) !PCM
C              ENDDO
C            ENDDO
          ENDDO
        ELSE
          DO year=yr0,yrf
C            write(*,*) 'year,siteno,year0,nyears,yr0,yrf'
C            write(*,*) year,siteno,year0,nyears,yr0,yrf
C            write(*,*) '(siteno-1)*nyears+year-year0+1'
C            write(*,*) (siteno-1)*nyears+year-year0+1
            READ(fno+1,1001, !PCM
     & REC=(siteno-1)*nyears+year-year0+1) ii,
     &((tmpv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
1001        FORMAT(A7,360I4)                
1002        FORMAT(A7,360F7.2)                
            READ(fno+2,1001,                   !PCM
     & REC=(siteno-1)*nyears+year-year0+1) ii,
     &((humv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
            READ(fno+3,1001,                   !PCM
     & REC=(siteno-1)*nyears+year-year0+1) ii,
     &((prcv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
            if(read_par.eq.1)
     &READ(fno+4,1002,                       !PCM
     & REC=(siteno-1)*nyears+year-year0+1)
     &ii, ((swrv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
            READ(fno+5,1001,                   !PCM
     & REC=(siteno-1)*nyears+year-year0+1) ii,
     &((wndv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
C PCM            DO mnth=1,12
C PCM              DO day=1,30
C PCM                prcv(year-yr0+1,mnth,day) = 
C PCM  &int(real(prcv(year-yr0+1,mnth,day))/10.0 + 0.5) 
C PCM              ENDDO
C PCM            ENDDO
          ENDDO
        ENDIF
        CLOSE(fno+1)
        CLOSE(fno+2)
        CLOSE(fno+3)
        CLOSE(fno+5)

      ENDIF

      ELSE
        siteno = 0
      ENDIF

      do mnth=1,12
        ans(mnth) = 0.0d0
        do day=1,30
          ans(mnth) = ans(mnth) + real(tmpv(1,mnth,day))/100.0d0
        enddo
      enddo

C   PCM added the following triple DO loop
C     PCM : Each of these 4 variables
C           are required to be passed to sdgvm0 in units of
C           0.01 oC, 0.1 mm, 0.01 % hum, 1 mm/s
C           They are read as:
C           0.1 oC, 0.01 mm, 0.1 % hum, 1 m/s
C     The following scalars convert from read units to sdgvm0 expected units
      TMP_MULT = 10.0                         
      HUM_MULT = 10.0
      PRC_MULT = 1/10.0 
      WND_MULT = 0.001 
      DO year=year0,yearf
       DO mnth=1,12
        DO day=1,no_days(year,mnth,0)
         IF ((year.LE.yrf).AND.(year.GE.yr0)) THEN
          xtmpv(year-yr0+1,mnth,day)= tmpv(year-yr0+1,mnth,day)*TMP_MULT
          xprcv(year-yr0+1,mnth,day)= prcv(year-yr0+1,mnth,day)*PRC_MULT
          xhumv(year-yr0+1,mnth,day)= humv(year-yr0+1,mnth,day)*HUM_MULT
          xwndv(year-yr0+1,mnth,day)= wndv(year-yr0+1,mnth,day)*HUM_MULT
         ENDIF
        ENDDO
      ENDDO
      ENDDO

C PCM2      WRITE(*,*) '00000000'
C PCM2      WRITE(*,*) 'TEMPERATURE (deg C)'
C PCM2      DO mnth=1,12
C PCM2        WRITE(*,'(30F5.1)') xtmpv(1,mnth,1:30)/TMP_MULT/TMP_MULT
C PCM2      ENDDO
C PCM2      WRITE(*,*) 'PRECIP/DAY'
C PCM2      DO mnth=1,12
C PCM2        WRITE(*,'(30F5.2)') xprcv(1,mnth,1:30)/PRC_MULT
C PCM2      ENDDO
C PCM2      WRITE(*,*) 'HUMIDITY (%)'
C PCM2      DO mnth=1,12
C PCM2        WRITE(*,'(30F5.1)') xhumv(1,mnth,1:30)/HUM_MULT/HUM_MULT 
C PCM2      ENDDO
C PCM2      WRITE(*,*) '111111111'

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE EX_CLIM_SITE                     *
*                          ***********************                     *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE EX_CLIM_SITE(stinput,yr0,yrf,tmpv,humv,prcv,wndv,year0,
     &yearf,swrv,read_par)
*----------------------------------------------------------------------*
      REAL*8 tmp,prc,hum,swr,swrv(500,12,31),wnd
      REAL*8 tmpv(500,12,31),humv(500,12,31),prcv(500,12,31)
      REAL*8 wndv(500,12,31)
      INTEGER year,year0,yearf,yr0,mnth,day,yrf,iyear,imnth,iday
      INTEGER blank,no_days,read_par
      CHARACTER stinput*1000

      OPEN(91,file=stinput(1:blank(stinput))//'/site.dat')

      DO year=year0,yearf
        DO mnth=1,12
          DO day=1,no_days(year,mnth,0)
            IF(read_par.eq.1) THEN
              READ(91,*) iyear,imnth,iday,tmp,prc,hum,wnd,swr
            ELSE
              READ(91,*) iyear,imnth,iday,tmp,prc,hum,wnd
            ENDIF

            IF ((iday.NE.day).OR.(imnth.NE.mnth).OR.(iyear.NE.year)) 
     &THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,*) 'Error in climate data file',year,mnth,day
              STOP
            ENDIF

            IF ((year.LE.yrf).AND.(year.GE.yr0)) THEN
              tmpv(year-yr0+1,mnth,day) = tmp*100.0d0
              prcv(year-yr0+1,mnth,day) = prc*10.0d0
              humv(year-yr0+1,mnth,day) = hum*100.0d0
              wndv(year-yr0+1,mnth,day) = wnd/1000.0d0
              IF(read_par.eq.1) swrv(year-yr0+1,mnth,day) = swr
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CLOSE(91)


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE EX_CLIM_SITE_MONTH               *
*                          *****************************               *
*----------------------------------------------------------------------*
      SUBROUTINE EX_CLIM_SITE_MONTH(stinput,yr0,yrf,tmpv,humv,prcv,wndv,
     &cldv,year0,yearf)
*----------------------------------------------------------------------*
      REAL*8 tmp(12),prc(12),hum(12),cld(12),swr(12),wnd(12)
      REAL*8 swrv(500,12,31),cldv(500,12),wndv(500,12,31)
      REAL*8 tmpv(500,12,31),humv(500,12,31),prcv(500,12,31)
      INTEGER year,year0,yearf,yr0,mnth,day,yrf,iyear,imnth,iday
      INTEGER blank,no_days,fno,cld_default,kode,read_par
      LOGICAL cloud
      CHARACTER stinput*1000

      cld_default = 50

      fno = 90

      OPEN(fno+1,file=stinput(1:blank(stinput))//'/tmp.dat',
     &status='old',iostat=kode)
      OPEN(fno+2,file=stinput(1:blank(stinput))//'/prc.dat',
     &status='old',iostat=kode)
      OPEN(fno+3,file=stinput(1:blank(stinput))//'/hum.dat',
     &status='old',iostat=kode)
      INQUIRE(file=stinput(1:blank(stinput))//'/cld.dat',
     &exist=cloud)
      IF (cloud) THEN
        OPEN(fno+4,file=stinput(1:blank(stinput))//'/cld.dat',
     &status='old',iostat=kode)
      ENDIF
      IF (read_par.eq.1) THEN
        OPEN(fno+5,file=stinput(1:blank(stinput))//'/swr.dat',
     &status='old',iostat=kode)
      ENDIF
      OPEN(fno+6,file=stinput(1:blank(stinput))//'/wnd.dat',
     &status='old',iostat=kode)

      DO year=year0,yearf
        READ(fno+1,*) iyear,tmp
        READ(fno+2,*) iyear,prc
        READ(fno+3,*) iyear,hum
        IF (cloud)     READ(fno+4,*) iyear,cld
        IF (read_par.eq.1)  READ(fno+5,*) iyear,cld
        READ(fno+6,*) iyear,wnd
        IF (iyear.NE.year) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'Error in climate data file',year,mnth,day
          STOP
        ENDIF

        IF ((year.LE.yrf).AND.(year.GE.yr0)) THEN
          DO mnth=1,12
            tmpv(year-yr0+1,mnth,1) = tmp(mnth)
            prcv(year-yr0+1,mnth,1) = prc(mnth)
            humv(year-yr0+1,mnth,1) = hum(mnth)
            IF (cloud) THEN
              cldv(year-yr0+1,mnth) = cld(mnth)
            ELSE
              cldv(year-yr0+1,mnth) = cld_default
            ENDIF
            IF (read_par.eq.1) swrv(year-yr0+1,mnth,1) = swr(mnth)
            wndv(year-yr0+1,mnth,1) = wnd(mnth)
          ENDDO
        ENDIF

      ENDDO

      CLOSE(fno+1)
      CLOSE(fno+2)
      CLOSE(fno+3)
      CLOSE(fno+4)
      CLOSE(fno+5)
      CLOSE(fno+6)

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE EX_SOIL                          *
*                          ******************                          *
*                                                                      *
* Extract % sand % silt bulk density and depth from soils database.    *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE EX_SOIL(fname1,lat,lon,sol_chr2,du,l_soil)
*----------------------------------------------------------------------*
      REAL*8 lat,lon,lon0,latf,latr,lonr,xlat,xlon,sol_chr2(10)
      INTEGER blank,row,col,recn,recl1,du,latn,lonn,i,ii,jj,kode
      CHARACTER fname1*1000
      INTEGER indx1(4,4),indx2(4,4),indx3(4,4),indx4(4,4),indx5(4,4)
      INTEGER indx6(4,4),indx7(4,4),indx8(4,4),indx9(4,4),indx10(4,4)
      REAL*8 xx1(4,4),xx2(4,4),xx3(4,4),xx4(4,4)
      REAL*8 xx5(4,4),xx6(4,4),xx7(4,4),xx8(4,4),xx9(4,4),xx10(4,4)
      REAL*8 ynorm,xnorm,rrow,rcol
      REAL*8 ans
      LOGICAL l_soil(20)

      IF (du.eq.1) THEN
        recl1 = 16+10*10
      ELSE
        recl1 = 16+10*10+1
      ENDIF

      OPEN(99,FILE=fname1(1:blank(fname1))//'/readme.dat',STATUS='old',
     &iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Soil data file does not exist.'
        WRITE(*,'('' "'',A,''/readme.dat"'')') fname1(1:blank(fname1))
        STOP
      ENDIF

      READ(99,*)
      READ(99,*) latf,lon0
      READ(99,*)
      READ(99,*) latr,lonr
      READ(99,*)
      READ(99,*) latn,lonn
      CLOSE(99)

*----------------------------------------------------------------------*
* Find the real row col corresponding to lat and lon.                  *
*----------------------------------------------------------------------*
      rrow = (latf - lat)/latr
      rcol = (lon - lon0)/lonr


      ynorm = rrow - real(int(rrow))
      xnorm = rcol - real(int(rcol))
*----------------------------------------------------------------------*

      OPEN(99,FILE=fname1(1:blank(fname1))//'/data.dat',STATUS='old',
     &FORM='formatted',ACCESS='direct',RECL=recl1,iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Soil data-base.'
        WRITE(*,*) 'File does not exist:',fname1(1:blank(fname1)),
     &'/data.dat'
        WRITE(*,*) 'Or record length missmatch, 16 + 9*10.'
        STOP
      ENDIF

      DO ii=1,4
        DO jj=1,4
          row = int(rrow)+jj-1
          col = int(rcol)+ii-1
          IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.
     &(col.LE.lonn)) THEN
            recn = (row-1)*lonn + col
            READ(99,'(f7.3,f9.3,10f10.4)',rec=recn)
     & xlat,xlon,(sol_chr2(i),i=1,10)
            xx1(ii,jj) = sol_chr2(1)
            xx2(ii,jj) = sol_chr2(2)
            xx3(ii,jj) = sol_chr2(3)
            xx4(ii,jj) = sol_chr2(4)
            xx5(ii,jj) = sol_chr2(5)
            xx6(ii,jj) = sol_chr2(6)
            xx7(ii,jj) = sol_chr2(7)
            xx8(ii,jj) = sol_chr2(8)
            xx9(ii,jj) = sol_chr2(9)
            xx10(ii,jj) = sol_chr2(10)
            IF (sol_chr2(1).LT.0.0d0) THEN
              indx1(ii,jj) = 0
            ELSE
              indx1(ii,jj) = 1
            ENDIF
            IF (sol_chr2(2).LT.0.0d0) THEN
              indx2(ii,jj) = 0
            ELSE
              indx2(ii,jj) = 1
            ENDIF
            IF (sol_chr2(3).LT.0.0d0) THEN
              indx3(ii,jj) = 0
            ELSE
              indx3(ii,jj) = 1
            ENDIF
            IF (sol_chr2(4).LT.0.0d0) THEN
              indx4(ii,jj) = 0
            ELSE
              indx4(ii,jj) = 1
            ENDIF
            IF (sol_chr2(5).LT.0.0d0) THEN
              indx5(ii,jj) = 0
            ELSE
              indx5(ii,jj) = 1
            ENDIF
            IF (sol_chr2(6).LT.0.0d0) THEN
              indx6(ii,jj) = 0
            ELSE
              indx6(ii,jj) = 1
            ENDIF
            IF (sol_chr2(7).LT.0.0d0) THEN
              indx7(ii,jj) = 0
            ELSE
              indx7(ii,jj) = 1
            ENDIF
            IF (sol_chr2(8).LT.0.0d0) THEN
              indx8(ii,jj) = 0
            ELSE
              indx8(ii,jj) = 1
            ENDIF
            IF (sol_chr2(9).LT.0.0d0) THEN
              indx9(ii,jj) = 0
            ELSE
              indx9(ii,jj) = 1
            ENDIF
            IF (sol_chr2(10).LT.0.0d0) THEN
              indx10(ii,jj) = 0
            ELSE
              indx10(ii,jj) = 1
            ENDIF
          ELSE
            indx1(ii,jj) = -1
            indx2(ii,jj) = -1
            indx3(ii,jj) = -1
            indx4(ii,jj) = -1
            indx5(ii,jj) = -1
            indx6(ii,jj) = -1
            indx7(ii,jj) = -1
            indx8(ii,jj) = -1
            indx9(ii,jj) = -1
            indx10(ii,jj) = -1
          ENDIF
        ENDDO
      ENDDO

      CLOSE(99)

      IF ((indx1(2,2).EQ.1).OR.(indx1(2,3).EQ.1).OR.(indx1(3,2).EQ.1).OR
     &.(indx1(3,3).EQ.1)) THEN
        CALL BI_LIN(xx1,indx1,xnorm,ynorm,ans)
        sol_chr2(1) = ans
        l_soil(1) = .TRUE.
      ELSE
        l_soil(1) = .FALSE.
      ENDIF
      IF ((indx2(2,2).EQ.1).OR.(indx2(2,3).EQ.1).OR.(indx2(3,2).EQ.1).OR
     &.(indx2(3,3).EQ.1)) THEN
        CALL BI_LIN(xx2,indx2,xnorm,ynorm,ans)
        sol_chr2(2) = ans
        l_soil(2) = .TRUE.
      ELSE
        l_soil(2) = .FALSE.
      ENDIF
      IF ((indx3(2,2).EQ.1).OR.(indx3(2,3).EQ.1).OR.(indx3(3,2).EQ.1).OR
     &.(indx3(3,3).EQ.1)) THEN
        CALL BI_LIN(xx3,indx3,xnorm,ynorm,ans)
        sol_chr2(3) = ans
        l_soil(3) = .TRUE.
      ELSE
        l_soil(3) = .FALSE.
      ENDIF
      IF ((indx4(2,2).EQ.1).OR.(indx4(2,3).EQ.1).OR.(indx4(3,2).EQ.1).OR
     &.(indx4(3,3).EQ.1)) THEN
        CALL BI_LIN(xx4,indx4,xnorm,ynorm,ans)
        sol_chr2(4) = ans
        l_soil(4) = .TRUE.
      ELSE
        l_soil(4) = .FALSE.
      ENDIF
      IF ((indx5(2,2).EQ.1).OR.(indx5(2,3).EQ.1).OR.(indx5(3,2).EQ.1).OR
     &.(indx5(3,3).EQ.1)) THEN
        CALL BI_LIN(xx5,indx5,xnorm,ynorm,ans)
        sol_chr2(5) = ans
        l_soil(5) = .TRUE.
      ELSE
        l_soil(5) = .FALSE.
      ENDIF
      IF ((indx6(2,2).EQ.1).OR.(indx6(2,3).EQ.1).OR.(indx6(3,2).EQ.1).OR
     &.(indx6(3,3).EQ.1)) THEN
        CALL BI_LIN(xx6,indx6,xnorm,ynorm,ans)
        sol_chr2(6) = ans
        l_soil(6) = .TRUE.
      ELSE
        l_soil(6) = .FALSE.
      ENDIF
      IF ((indx7(2,2).EQ.1).OR.(indx7(2,3).EQ.1).OR.(indx7(3,2).EQ.1).OR
     &.(indx7(3,3).EQ.1)) THEN
        CALL BI_LIN(xx7,indx7,xnorm,ynorm,ans)
        sol_chr2(7) = ans
        l_soil(7) = .TRUE.
      ELSE
        l_soil(7) = .FALSE.
      ENDIF
      IF ((indx8(2,2).EQ.1).OR.(indx8(2,3).EQ.1).OR.(indx8(3,2).EQ.1).OR
     &.(indx8(3,3).EQ.1)) THEN
        CALL BI_LIN(xx8,indx8,xnorm,ynorm,ans)
        sol_chr2(8) = ans
        l_soil(8) = .TRUE.
      ELSE
        l_soil(8) = .FALSE.
      ENDIF
      IF ((indx9(2,2).EQ.1).OR.(indx9(2,3).EQ.1).OR.(indx9(3,2).EQ.1).OR
     &.(indx9(3,3).EQ.1)) THEN
        CALL BI_LIN(xx9,indx9,xnorm,ynorm,ans)
        sol_chr2(9) = ans
        l_soil(9) = .TRUE.
      ELSE
        l_soil(9) = .FALSE.
      ENDIF
      IF ((indx10(2,2).EQ.1).OR.(indx10(2,3).EQ.1).OR.(indx10(3,2).EQ.1)
     &.OR.(indx10(3,3).EQ.1)) THEN
        CALL BI_LIN(xx10,indx10,xnorm,ynorm,ans)
        sol_chr2(10) = ans
        l_soil(10) = .TRUE.
      ELSE
        l_soil(10) = .FALSE.
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE EX_CLUMP                         *
*                          ******************                          *
*                                                                      *
* Extract clumpin index from xxxx database.                            *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE EX_CLUMP(fname1,lat,lon,var,du)
*----------------------------------------------------------------------*
      REAL*8 lat,lon,lon0,latf,latr,lonr,xlat,xlon,var,rvar
      INTEGER blank,row,col,recn,recl1,du,latn,lonn,i,ii,jj,kode
      CHARACTER fname1*1000
      INTEGER indx1(4,4)
      REAL*8 xx1(4,4)
      REAL*8 ynorm,xnorm,rrow,rcol
      REAL*8 ans
      LOGICAL l_soil

      IF (du.eq.1) THEN
        recl1 = 6
      ELSE
        recl1 = 6+1
      ENDIF

      OPEN(99,FILE=fname1(1:blank(fname1))//'/readme.dat',STATUS='old',
     &iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Canopy clumping data file does not exist.'
        WRITE(*,'('' "'',A,''/readme.dat"'')') fname1(1:blank(fname1))
        STOP
      ENDIF

      READ(99,*)
      READ(99,*)
      READ(99,*) latf,lon0
      READ(99,*)
      READ(99,*) latr,lonr
      READ(99,*)
      READ(99,*) latn,lonn
      CLOSE(99)

*----------------------------------------------------------------------*
* Find the real row col corresponding to lat and lon.                  *
* PCM: offset this by 2 gridcells, since we have a box of 4 x 4        *
*----------------------------------------------------------------------*
      rrow = 1.0 + (latf - lat)/latr - 2.0
      rcol = 1.0 + (lon - lon0)/lonr - 2.0


      ynorm = rrow - real(int(rrow))
      xnorm = rcol - real(int(rcol))
*----------------------------------------------------------------------*

      OPEN(99,FILE=fname1(1:blank(fname1))//'/clumpingFLIP.dat',
     &STATUS='old',FORM='formatted',ACCESS='direct',RECL=recl1,
     &iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Canopy clumping data-base.'
        WRITE(*,*) 'File does not exist:',fname1(1:blank(fname1)),
     &'/clumpingFLIP.dat'
        WRITE(*,*) 'Or record length missmatch, 6.'
        STOP
      ENDIF

      DO ii=1,4
        DO jj=1,4
          row = int(rrow)+jj-1
          col = int(rcol)+ii-1
          IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.
     &(col.LE.lonn)) THEN
            recn = (row-1)*lonn + col
            READ(99,'(f5.3)',rec=recn) rvar
            xx1(ii,jj) = rvar
            IF (rvar.LT.0.0d0) THEN
              indx1(ii,jj) = 0
            ELSE
              indx1(ii,jj) = 1
            ENDIF
          ELSE
            indx1(ii,jj) = -1
          ENDIF
        ENDDO
      ENDDO

      CLOSE(99)

      IF ((indx1(2,2).EQ.1).OR.(indx1(2,3).EQ.1).OR.(indx1(3,2).EQ.1).OR
     &.(indx1(3,3).EQ.1)) THEN
        CALL BI_LIN(xx1,indx1,xnorm,ynorm,ans)
        var = ans
        l_soil = .TRUE.
      ELSE
        PRINT*, 'No clumping value in any central four gridcells of the'
     & //' clumping database' 
        l_soil = .FALSE.
      ENDIF

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE EX_LU                            *
*                          ****************                            *
*                                                                      *
* Extract land use from land use soils database for the nearest        *
* site to lat,lon, replace lat,lon with the nearest cell from the      *
* database.                                                            *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE EX_LU(fname1,lat,lon,luse,yr0,yrf,du)
*----------------------------------------------------------------------*
      INCLUDE 'array_dims.inc'
      REAL*8 lat,lon,lon0,latf,latr,lonr,xlat,xlon
      INTEGER i,n_fields4000,n,j,du,latn,lonn,blank,row,col,recn,kode
      INTEGER luse(maxyrs),yr0,yrf,rep,years(1000),lu(1000),nrecl
      CHARACTER fname1*1000,st1*4000

      OPEN(99,FILE=fname1(1:blank(fname1))//'/readme.dat')
      READ(99,*)
      READ(99,*) latf,lon0
      READ(99,*)
      READ(99,*) latr,lonr
      READ(99,*)
      READ(99,*) latn,lonn
      READ(99,*)
      READ(99,'(A)') st1
      CLOSE(99)

      n = n_fields4000(st1)
      CALL ST2ARR4000(st1,years,4000,n)

      IF (du.eq.1) THEN
        nrecl = 16+n*3
      ELSE
        nrecl = 16+n*3+1
      ENDIF

      lu(1) = 0
      OPEN(99,FILE=fname1(1:blank(fname1))//'/landuse.dat',STATUS='old',
     &FORM='formatted',ACCESS='direct',RECL=nrecl,iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'File does not exist:',fname1(1:blank(fname1)),'/lan
     &duse.dat'
        STOP
      ENDIF


      row = int((latf - lat)/latr + 1.0d0)
      col = int((lon - lon0)/lonr + 1.0d0)

      recn = (row-1)*lonn + col
      IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.(col.LE.lonn))
     &READ(99,'(f7.3,f9.3,1000i3)',REC=recn) xlat,xlon,(lu(i),i=1,n)

      IF (lu(1).eq.0) THEN
        row = int((latf - (lat + latr))/latr + 1.0d0)
        col = int(((lon) - lon0)/lonr + 1.0d0)
        recn = (row-1)*lonn + col
      IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.(col.LE.lonn))
     &READ(99,'(f7.3,f9.3,1000i3)',
     &REC=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
      ENDIF

      IF (lu(1).eq.0) THEN
        row = int((latf - (lat - latr))/latr + 1.0d0)
        col = int(((lon) - lon0)/lonr + 1.0d0)
        recn = (row-1)*lonn + col
      IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.(col.LE.lonn))
     &READ(99,'(f7.3,f9.3,1000i3)',
     &REC=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
      ENDIF

      IF (lu(1).eq.0) THEN
        row = int((latf - (lat))/latr + 1.0d0)
        col = int(((lon + lonr) - lon0)/lonr + 1.0d0)
        recn = (row-1)*lonn + col
      IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.(col.LE.lonn))
     &READ(99,'(f7.3,f9.3,1000i3)',
     &REC=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
      ENDIF

      IF (lu(1).eq.0) THEN
        row = int((latf - (lat))/latr + 1.0d0)
        col = int(((lon - lonr) - lon0)/lonr + 1.0d0)
        recn = (row-1)*lonn + col
      IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.(col.LE.lonn))
     &READ(99,'(f7.3,f9.3,1000i3)',
     &REC=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
      ENDIF

      IF (lu(1).eq.0) THEN
        row = int((latf - (lat + latr))/latr + 1.0d0)
        col = int(((lon + lonr) - lon0)/lonr + 1.0d0)
        recn = (row-1)*lonn + col
      IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.(col.LE.lonn))
     &READ(99,'(f7.3,f9.3,1000i3)',
     &REC=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
      ENDIF

      IF (lu(1).eq.0) THEN
        row = int((latf - (lat - latr))/latr + 1.0d0)
        col = int(((lon + lonr) - lon0)/lonr + 1.0d0)
        recn = (row-1)*lonn + col
      IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.(col.LE.lonn))
     &READ(99,'(f7.3,f9.3,1000i3)',
     &REC=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
      ENDIF

      IF (lu(1).eq.0) THEN
        row = int((latf - (lat + latr))/latr + 1.0d0)
        col = int(((lon - lonr) - lon0)/lonr + 1.0d0)
        recn = (row-1)*lonn + col
      IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.(col.LE.lonn))
     &READ(99,'(f7.3,f9.3,1000i3)',
     &REC=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
      ENDIF

      IF (lu(1).eq.0) THEN
        row = int((latf - (lat - latr))/latr + 1.0d0)
        col = int(((lon - lonr) - lon0)/lonr + 1.0d0)
        recn = (row-1)*lonn + col
      IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.(col.LE.lonn))
     &READ(99,'(f7.3,f9.3,1000i3)',
     &REC=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
      ENDIF

      CLOSE(99)

      rep = lu(1)
      j = 1
      DO i=1,n
        IF (yr0.GE.years(j+1)) THEN
          rep = lu(i)
          j = i
        ENDIF
      ENDDO

      DO i=1,yrf-yr0+1
        IF ((i+yr0-1.GE.years(j+1)).AND.(j.LT.n)) THEN
          j = j+1
          rep = lu(j)
        ENDIF
        luse(i) = rep
      ENDDO


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE EX_CLU                           *
*                          *****************                           *
*                                                                      *
* Extract land use from land use soils database for the nearest        *
* site to lat,lon, replace lat,lon with the nearest cell from the      *
* database.                                                            *
* The data contains the percent (in byte format) for each pft from 2   *
* to nft                                                               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE EX_CLU(fname1,lat,lon,nft,lutab,cluse,du,l_lu,
     &yr0a,yrfa,year0set,spinl,ilanduse,SYR,NYR,NYR_FILE,
     &SYR_FIRE,NYR_FIRE,prescr_fire,lutab2,pname,pname_t,pname_f,wdg,
     &cluse2,cluseh,fprob_prescrh,debug)
*----------------------------------------------------------------------*
      USE FUNCTIONS_CLU
      INCLUDE 'array_dims.inc'
      INTEGER, PARAMETER :: NX = 720, NY = 360
      INTEGER, PARAMETER :: NS = 17
      INTEGER, PARAMETER :: maxn_at = 7 !max number of aggregrated (functional) types
      INTEGER :: NYR,SYR,NYR_FILE 
      INTEGER :: NYR_FIRE,SYR_FIRE 
      REAL*8 lat,lon,lon0,latf,latr,lonr,classprop(255)
      REAL*8 cluse(maxnft,maxyrs),lutab(255,100),ans
      REAL*8 ftprop(maxnft),rrow,rcol,xx(4,4),xnorm,ynorm,co2const
      REAL*8 cluse2(maxn_at,maxn_at,maxyrs)
      REAL*8 cluseh(maxn_at,maxyrs)
      REAL*8 fprob_prescrh(maxyrs,12),fprob_prescr(12)
      REAL*8 SDGVM_LUC(NYR, NS, 4, 4),xf,x2,xf2,xx2(4,4,maxn_at)
      REAL*8 SDGVM_LUC2(NYR, maxn_at, maxn_at, 4, 4)
      REAL*8 SDGVM_LUC2_HARVEST(NYR, maxn_at, 4, 4)
      REAL*8 SDGVM_LUC_FIRE(NYR_FIRE, 4, 4, 12)
      REAL*8 agclassprop2(255,255), atprop2(maxn_at,maxn_at)
      REAL*8 agclassprop(255),atharvest(maxn_at)
      INTEGER i,n_fields,n,j,du,latn,lonn,blank,row,col,recn,k,x,nft,ift
      INTEGER ii,jj,stcmp,indx(4,4),years(1000),nrecl,yr0a,yrfa
      INTEGER classes(1000),nclasses,kode,spinl,yr_offset,yr_offset_fire
      INTEGER ij,ij1,j1,num_land,years_fire(1000)
      INTEGER ilanduse,k2,k3,agclasses(1000),indx2(4,4,maxn_at)
      INTEGER iat2,iat3
      INTEGER n_fire,j_fire,j1_fire,mm,years_NYR_FILE
      REAL*8 lutab2(255,100) 
      CHARACTER fname1*1000,st1*1000,st2*1000,in2st*1000,st3*1000
      CHARACTER pname*1000,pname_t*1000,wdg*1000,pname_f*1000
      CHARACTER st4*4000
      INTEGER n_fields4000
      LOGICAL l_lu,year0set,prescr_fire
      LOGICAL compute_next_year,get_transitions,mindx(4,4)
      LOGICAL do_fire_interp
      LOGICAL :: debug !used to print out more debugging info


      l_lu = .TRUE.
      IF(ilanduse.GE.3 .AND. ilanduse.LE.6 ) THEN !PCM Use states2b.nc or transitions2b.nc file: half-res
        !PCM hardwire these numbers for now
        ! APW: need to add a readmeLUH2.dat here
        latf = 89.75
        lon0 = -179.75
        latr = 0.5
        lonr = 0.5
        latn = 360
        lonn = 720
        n    = NYR 
        years(1:n) = (/(i, i=SYR,SYR+n-1)/)
        years_NYR_FILE = SYR + NYR_FILE - 1
!        write(*,*)'SYR,years(1)=',SYR,years(1)
!        if(prescr_fire .EQV. .TRUE.) THEN
        n_fire = NYR_FIRE !1901-2020=120
        !write(*,*)'NYR_FIRE=',n_fire
        years_fire(1:n_fire) = (/(i, i=SYR_FIRE,SYR_FIRE+n_fire-1)/)
!        END IF
        nclasses   = NS 
        classes(1:nclasses) = (/(i, i=1,nclasses)/)
        agclasses(1:maxn_at) = (/(i, i=1,maxn_at)/)

      ELSEIF(ilanduse.EQ.0) THEN !PCM original method of using SDGVM land use 
*----------------------------------------------------------------------*
* Read in the readme file 'readme.dat'.                                *
*----------------------------------------------------------------------*
        OPEN(99,FILE=fname1(1:blank(fname1))//'/readme.dat',
     &   status='old',iostat=kode)
        IF (kode.NE.0) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'Land use file does not exist.'
          WRITE(*,'('' "'',A,''/readme.dat"'')') fname1(1:blank(fname1))
          STOP
        ENDIF

        READ(99,*) st1
        st2='CONTINUOUS'
        IF (stcmp(st1,st2).EQ.0) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'landuse is not a continuous field ?'
          WRITE(*,*) 'readme.dat should begin with CONTINUOUS'
          STOP
        ENDIF
        READ(99,*)
        READ(99,*) latf,lon0
        READ(99,*)
        READ(99,*) latr,lonr
        READ(99,*)
        READ(99,*) latn,lonn
        READ(99,*)
        READ(99,'(A)') st4
        n = n_fields4000(st4)
        CALL ST2ARR4000(st4,years,1000,n)
        READ(99,*)
        READ(99,'(A)') st1
        CLOSE(99)
        nclasses = n_fields(st1)
        CALL ST2ARR(st1,classes,1000,nclasses)
*----------------------------------------------------------------------*

        IF (du.eq.1) THEN
          nrecl = 3
        ELSE
          nrecl = 4
        ENDIF
      ENDIF !PCM 

      yr_offset = 0 
      IF ((n.GT.1).AND.(yr0a.LT.years(1))) THEN
          WRITE(*,*) 'Cannot start running in year',yr0a,
     &      ' since landuse map begins in ',years(1)
          WRITE(*,*) 'The first year of land use data
     &will be used and the first year of CO2 data will be used in',
     & yr0a,'. This has the potential to uncouple the land use year
     &from the CO2 year if their datasets start in different years.' 
          WRITE(*,*) 'If this message is telling you that
     &the simulation was requested to start in year 1, then
     &you have requested a spin only with co2const < 0.'
          WRITE(*,*) 'Alternatively you have requested a spin and a run
     &proper and land use and CO2 will be out of sync with climate.'  
          yr_offset = years(1) - yr0a 
          !STOP
      ENDIF


      yr_offset_fire = 0 
!      IF ((n_fire.GT.1).AND.(yr0a.LT.years_fire(1))) THEN
!          WRITE(*,*) 'Cannot start running in year',yr0a,
!     &      ' since prescribed fire map begins in ',years_fire(1)
!          WRITE(*,*) 'The first year of prescribed fire data
!     &will be used and the first year of CO2 data will be used in',
!     & yr0a,'. This has the potential to uncouple the prescribed-fire
!     & year from the CO2 year if their datasets start in different
!     & years.' 
!          WRITE(*,*) 'If this message is telling you that
!     &the simulation was requested to start in year 1, then
!     &you have requested a spin only with co2const < 0.'
!          WRITE(*,*) 'Alternatively you have requested a spin and a run
!     &proper and prescribed fire and CO2 will be out of sync
!     &with climate.'  
!          yr_offset_fire = years_fire(1) - yr0a 
!          !STOP
!      ENDIF

c     look for the first year
      ! Currently this selects the starting landuse of the simulation to be the landuse specified in the 
      ! dataset in the exact start year of the simulation, or the first year after that in the landuse dataset 
      ! if there is no landcover in the exact start year. Thus all landcover information in 
      ! the years prior to the simulation start year are ignored and the land-use from the 
      ! beginning year of the simulation to the first year of land-cover data is the land-use 
      ! specified in the first land-use data year within the range of simulation years.
      ! This is the case unless a spin only is selected with co2const < 0.
      ! This gives the slightly strange behaviour when a spin or spin and run proper are selected, and
      ! co2const > 0, that landuse data is taken only from the range of run proper years even if land use 
      ! data are available for the effective spin years. 
      ! This also applies to land-use specified in the input.dat file.
   
      j=1

   10 CONTINUE
      IF ((j.LT.n).AND.(years(j)-yr_offset.LT.yr0a)) THEN
         j = j + 1
         GOTO 10
      ENDIF
      j1 = j

      j_fire=1
   20 CONTINUE
      IF(prescr_fire .EQV. .TRUE.) THEN
        IF ((j_fire.LT.n_fire)
     & .AND.(years_fire(j_fire)-yr_offset_fire.LT.yr0a)) THEN
         j_fire = j_fire + 1
         GOTO 20 
        ENDIF
      END IF
      j1_fire = j_fire

!      print*, 'EX_CLU:before_loop'
!      print*, 'j,n:', j,n
!      print*, 'start year:', years(j)
!      print*, 'year 2:', years(j+1)
!      print*, 'year n:', years(n)
!      print*, 'j_fire,n_fire:', j_fire,n_fire
!      print*, 'fire start year:', years_fire(j_fire)
!      print*, 'fire year 2:', years_fire(j_fire+1)
!      print*, 'fire year n_fire:', years_fire(n_fire)
*----------------------------------------------------------------------*
* Find the real row col corresponding to lat and lon.                  *
* PCM: offset this by 2 gridcells, since we have a box of 4 x 4        *
*----------------------------------------------------------------------*
      rrow = 1.0 + (latf - lat)/latr - 2.0
      rcol = 1.0 + (lon - lon0)/lonr - 2.0

      ynorm = rrow - real(int(rrow))
      xnorm = rcol - real(int(rcol))
*----------------------------------------------------------------------*
      IF(debug .EQV. .TRUE.) THEN
        PRINT *,rrow,rcol,lat,lon
      ENDIF

      IF(ilanduse.GE.3 .and. ilanduse.LE.6 ) THEN
CPCM Use states2b.nc (compute_next_year==.false.) or transitions2b.nc file (compute_next_year==.true.) 
CPCM get_transitions==.true. : compute transitions
        IF(ilanduse.EQ.3 .or. ilanduse.EQ.5 ) THEN
          compute_next_year = .false.
          get_transitions = .true. !PCM for now
        ELSE IF(ilanduse.EQ.4 .or. ilanduse.EQ.6 ) THEN
          !compute_next_year = .true. !for testing: add transitions to states; PCM
          compute_next_year = .false. 
          get_transitions = .true.
        ENDIF
        ! get the 4 neighboring grid cells for all nclasses for the years range from states2b.nc in the SDGVM_LUC variable
        ! get the 4 neighboring grid cells for all maxn_at*maxn_at for the years range from transitions2b.nc in the SDGVM_LUC2 variable
        !write(*,*)'BEFORE STATES, years(1)=',years(1)
        CALL states_convertSDGVM_func(
     &     years(1),years(n),years_NYR_FILE,
     &     years_fire(1),years_fire(n_fire),
     &     INT(rcol),INT(rrow),4, ! with the definition of rcol, rcol+2 starts at 1 for lon==lon0
     &     get_transitions,compute_next_year,prescr_fire, 
     &     pname,pname_t,pname_f,wdg,
     &     SDGVM_LUC, SDGVM_LUC2, SDGVM_LUC2_HARVEST,SDGVM_LUC_FIRE,
     &     debug)  
       ! WRITE(*,*)'EX_CLU: after states_convertSDGVM_func'
       ! write(*,FMT="(A)") "EX_CLU:SDGVM_LUC2_FIRE"
       ! DO jj=1,n_fire
       !   write(*,FMT="(12F7.4)") SDGVM_LUC_FIRE(jj,3,3,:)
       ! END DO
        IF(debug .EQV. .TRUE.) THEN
          WRITE(*,*)
     &'    t   BARE     Ev_Bp    Dc_Bp    Ev_Np    Dc_Np    ',
     &'Shrup    C3p      C4p      C3crop   C4crop   ',
     &'Ev_Bs    Dc_Bs    Ev_Ns    Dc_Ns    Shrus    C3s    C4s'
        ENDIF
  
!        IF (ANY(ABS(SDGVM_LUC(:,:,3,3)).GT.200.0)) THEN
        IF (ALL(ABS(SDGVM_LUC(:,:,3:4,3:4)).GT.200.0)) THEN
!          WRITE(*,*)'SDGVM_LUC: l_lu=.FALSE.'
          l_lu = .FALSE.
          RETURN
        ENDIF
  
!       check for all but the last year, since it may not be valid for the
!       last year
!        IF (ANY(ABS(SDGVM_LUC2(1:NYR-1,:,:,3,3)).GT.200.0)) THEN
        IF (NYR.GT.1) THEN
         IF (ALL(ABS(SDGVM_LUC2(1:NYR-1,:,:,3:4,3:4)).GT.200.0)) THEN 
!          WRITE(*,*)'SDGVM_LUC2: l_lu=.FALSE.'
          l_lu = .FALSE.
          RETURN
         ENDIF
        ENDIF
  
        IF ((prescr_fire .EQV. .TRUE.) .AND.
     &    (ALL(ABS(SDGVM_LUC_FIRE(:,3:4,3:4,1:12)).GT.200.0)) ) THEN
!          WRITE(*,*)'SDGVM_LUC_FIRE: l_lu=.FALSE.'
          l_lu = .FALSE.
          RETURN
        ENDIF
  
      ELSE
        SDGVM_LUC=0.0
        SDGVM_LUC2=0.0
        SDGVM_LUC2_HARVEST=0.0
        SDGVM_LUC_FIRE=0.0
        IF(debug .EQV. .TRUE.) THEN
         WRITE(*,*)
     &'    t   BARE     Ev_Bl    Dc_Bl    Ev_Nl    Dc_Nl    ',
     &'Shrub    C3       C4       C3crop   C4crop   '
         ENDIF
      ENDIF

      IF(debug .EQV. .TRUE.) THEN
        write(*,FMT="(A,7E10.3)") 'D harv',SDGVM_LUC2_HARVEST(1,:,3,3)
      ENDIF

      IF (prescr_fire .EQV. .TRUE.) THEN !PCM
        DO i=1,yrfa-yr0a+1
!       DO i=1,yrfa-yr0a+1
!!        print*, yr0a,yrfa 
!!        print*, j, years(j), i, years(j) - yr_offset,i+yr0a-1
!        IF ((i.EQ.1).OR.((i+yr0a-1).EQ.years(j)-yr_offset)) THEN
!          WRITE(*,*)'j_fire,n_fire',
!     &                 j_fire,n_fire
          IF(j_fire.LE.n_fire) THEN
              j_fire=j_fire+1
          ELSE
              j_fire=2 !first time in loop j_fire-1 = 1
          END IF

!          WRITE(*,*)'EX_CLU:before_AVE4,j_fire=',j_fire
          fprob_prescr(1:12) = 0.0d0
          DO mm=1,12
              DO ii=1,4
                DO jj=1,4
                  row = int(rrow)+jj-1
                  col = int(rcol)+ii-1
                  !PCM: currently, there is no wrapping at longitude of date-line
                  !for the 4x4 interpolation
                  IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.
     &               (col.LE.lonn)) THEN
                    xf = SDGVM_LUC_FIRE(j_fire-1, ii, jj, mm) !for j_fire=1, years(j)=1901
!                    xf = SDGVM_LUC_FIRE(j_fire-1, 3, 3, mm) !for j_fire=1, years(j)=1901
                    xx(ii,jj) = xf 
                    x = INT(xf) !need this for indx and mindx masking, below
                    IF (x.LT.2 .AND. x.GE.0) THEN
                      indx(ii,jj) = 1
                      mindx(ii,jj) = .true. 
                    ELSE
                      indx(ii,jj) = 0
                      mindx(ii,jj) = .false. 
                    ENDIF
                  ELSE
                    indx(ii,jj) = -1
                    mindx(ii,jj) = .false. 
                  ENDIF
                ENDDO
              ENDDO

              num_land = COUNT( mindx .EQV. .true.)
              !print*, num_land,mm
              
              CALL AVE4(xx,indx,xnorm,ynorm,ans)

              x = int(ans+0.5d0)

              fprob_prescr(mm) = ans !should be between 0 & 1
              !fprob_prescr(mm) = xf !should be between 0 & 1
          END DO
!          WRITE(*,'(A)')'EX_CLU:after_AVE4'
!          WRITE(*,'(12F7.4)')fprob_prescr
          fprob_prescrh(j_fire-1,1:12) = fprob_prescr(1:12)

!        ENDIF ! finished reading
        ENDDO !year loop #1

      END IF ! end of prescr_fire 


      ! read landcover dataset
      DO i=1,yrfa-yr0a+1
!        print*, yr0a,yrfa 
!        print*, j, years(j), i, years(j) - yr_offset,i+yr0a-1
        IF ((i.EQ.1).OR.((i+yr0a-1).EQ.years(j)-yr_offset)) THEN

          st2=in2st(years(j))
          CALL STRIPB(st2)
          !print*, st2(1:4) 
          j=j+1
          

          IF((MOD(years(j-1)-years(1),20)==0).AND.
     &(debug.EQV..TRUE.))THEN 
            WRITE(*,FMT='(A1)', ADVANCE='no') 'D' 

            IF(compute_next_year) THEN
               WRITE(*,FMT='(I5)', ADVANCE='no') years(j-1)+1
            ELSE
               WRITE(*,FMT='(I5)', ADVANCE='no') years(j-1) 
            END IF
          ENDIF


          DO k=1,nclasses
            classprop(classes(k)) = 0
            st3=in2st(classes(k))
            CALL STRIPB(st3)
            IF(ilanduse.EQ.0) THEN !PCM
C PCM: st2 is the year st3 is the class, ranging from 1 to NS (NS=10)
              OPEN(99,FILE=fname1(1:blank(fname1))//'/cont_lu-'//
     &st3(1:blank(st3))//'-'//st2(1:4)//'.dat',STATUS='old',
     &FORM='formatted',ACCESS='direct',RECL=nrecl,iostat=kode)
              IF (kode.NE.0) THEN
                WRITE(*,'('' PROGRAM TERMINATED'')')
                WRITE(*,*) 'Land Use data-base.'
                WRITE(*,*) 'File does not exist:'
                WRITE(*,*) fname1(1:blank(fname1)),
     &'/cont_lu-',st3(1:blank(st3)),'-',st2(1:4),'.dat'
                STOP
              ENDIF
            ENDIF !PCM 

            DO ii=1,4
                DO jj=1,4
                  row = int(rrow)+jj-1
                  col = int(rcol)+ii-1
                  !PCM: currently, there is no wrapping at longitude of date-line
                  !for the 4x4 interpolation
                  IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.
     &(col.LE.lonn)) THEN
                    recn = (row-1)*lonn + col
                    IF(ilanduse.GE.3 .AND. ilanduse.LE.6 ) THEN !PCM
                      xf = SDGVM_LUC(j-1, classes(k), ii, jj) !for j=1, years(j)=1700
                      xx(ii,jj) = xf 
                      x = INT(xf) !need this for indx and mindx masking, below
                    ELSEIF(ilanduse.EQ.0) THEN !PCM
                      READ(99,'(i3)',REC=recn) x
                      xx(ii,jj) = real(x) !convert from integer
                    ENDIF
                    IF (x.LT.200) THEN
                      indx(ii,jj) = 1
                      mindx(ii,jj) = .true. 
                    ELSE
                      indx(ii,jj) = 0
                      mindx(ii,jj) = .false. 
                    ENDIF
                  ELSE
                    indx(ii,jj) = -1
                    mindx(ii,jj) = .false. 
                  ENDIF
                ENDDO
            ENDDO

            num_land = COUNT( mindx .EQV. .true.)

C            IF( MOD(years(j-1)-years(1),20) == 0 ) THEN 
C              if(num_land .GT. 0) THEN
C                WRITE(*,FMT='(F9.4)',ADVANCE='no') 
C     &                 SUM(xx,mindx)/100.0/num_land 
CC                WRITE(*,FMT='(F9.4)',ADVANCE='no') xx(1,1)/100.0 
C              ELSE
CC                WRITE(*,FMT='(F9.4)', ADVANCE='no') -1.00 
C              ENDIF
C            ENDIF


            CALL AVE4(xx,indx,xnorm,ynorm,ans)

            x = int(ans+0.5d0)

            classprop(classes(k)) = ans

            IF((MOD(years(j-1)-years(1),20)==0)
     &.AND.(debug.EQV..TRUE.)) THEN 
              if(num_land .GT. 0) THEN
                WRITE(*,FMT='(F9.4)',ADVANCE='no')  
     &ans 
              ELSE
                WRITE(*,FMT='(F9.4)', ADVANCE='no') -1.00 
              ENDIF
            ENDIF

            IF(ilanduse.EQ.0) THEN !PCM
              CLOSE(99)
            ENDIF !PCM 


          ENDDO ! end of loop over the classes

          IF((MOD(years(j-1)-years(1),20)==0)
     &.AND.(debug.EQV..TRUE.))THEN 
              WRITE(*,*) ! Assumes default "ADVANCE='yes'".
          ENDIF


          IF(ilanduse.GE.3 .AND. ilanduse.LE.6 ) THEN !PCM
            DO k3=1,maxn_at
              DO k2=1,maxn_at
                agclassprop2(agclasses(k3),agclasses(k2)) = 0
              ENDDO

              DO ii=1,4
                DO jj=1,4
                  row = int(rrow)+jj-1
                  col = int(rcol)+ii-1
                  !PCM: currently, there is no wrapping at longitude of date-line
                  !for the 4x4 interpolation
                  IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.
     &(col.LE.lonn)) THEN
                    DO k2=1,maxn_at
                       xf2 = SDGVM_LUC2(j-1,
     &agclasses(k3),agclasses(k2),ii, jj) !for j=1, years(j)=1700
                       xx2(ii,jj,k2) = xf2 
                       IF (xf2.LT.200) THEN
                         indx2(ii,jj,k2) = 1
                       ELSE
                         indx2(ii,jj,k2) = 0
                       ENDIF
                       !IF(debug)THEN
                       !  WRITE(*,*) 'LUC2',ii,jj,xf2,indx2(ii,jj,k2)
                       !ENDIF
                    ENDDO
                  ELSE
                    indx2(ii,jj,:) = -1
                  ENDIF
                ENDDO
              ENDDO


              DO k2=1,maxn_at

C               IF( MOD(years(j-1)-years(1),20) == 0 ) THEN 
C
C                 if(num_land .GT. 0) THEN
C                   WRITE(*,FMT='(F9.4)',ADVANCE='no') 
C     &                 SUM(xx2(:,:,k2),mindx)/100.0/num_land 
C                 ELSE
C                   WRITE(*,FMT='(F9.4)', ADVANCE='no') -1.00 
C                 ENDIF
C               ENDIF

               CALL AVE4(xx2(:,:,k2),indx2(:,:,k2),xnorm,ynorm,ans)
               x2 = int(ans+0.5d0)
               agclassprop2(agclasses(k3),agclasses(k2)) = ans
C               print *,'DD1',k3,k2,agclasses(k3),agclasses(k2),
C     &agclassprop2(agclasses(k3),agclasses(k2))
              ENDDO

C            IF( MOD(years(j-1)-years(1),20) == 0 ) THEN 
C              WRITE(*,*) ! Assumes default "ADVANCE='yes'".
C            ENDIF
            ENDDO ! end of k3 loop over the agclasses
          END IF

          IF(ilanduse.GE.3 .AND. ilanduse.LE.6 ) THEN !PCM
            DO k3=1,maxn_at
              agclassprop(agclasses(k3)) = 0

              DO ii=1,4
                DO jj=1,4
                  row = int(rrow)+jj-1
                  col = int(rcol)+ii-1
                  !PCM: currently, there is no wrapping at longitude of date-line
                  !for the 4x4 interpolation
                  IF ((row.GE.1).AND.(row.LE.latn).AND.(col.GE.1).AND.
     &               (col.LE.lonn)) THEN
                     xf = SDGVM_LUC2_HARVEST(j-1,
     &                 agclasses(k3),ii, jj) !for j=1, years(j)=1700
                     xx(ii,jj) = xf
                     IF (xf.LT.200) THEN
                       indx(ii,jj) = 1
                     ELSE
                       indx(ii,jj) = 0
                     ENDIF
                  ELSE
                    indx(ii,jj) = -1
                  ENDIF
                ENDDO
              ENDDO


              CALL AVE4(xx(:,:),indx(:,:),xnorm,ynorm,ans)
              x = int(ans+0.5d0)
              agclassprop(agclasses(k3)) = ans

            ENDDO ! end of k3 loop over the agclasses
          END IF


c
c Now calculate the ftprop.
c
C PCM classprop is the percentage of each class in that gridcell, after
C     interpolation in the BI_LIN step above 
C PCM lutab(classes(k),ift) is the percentage of the SDGVM
C     land-cover-class that is of the labelled SDGVM ift
C PCM agclassprop2(cl3,cl2) is the percentage of each aggregated Hyde class transitioning
C     from agclass cl3 to agclass cl2 in that gridcell, after
C interpolation in the BI_LIN step above 
          DO ift=2,nft
            ftprop(ift)=0.0d0
            DO k=1,nclasses
              ftprop(ift)=ftprop(ift)+lutab(classes(k),ift)*
     &classprop(classes(k))/100.0d0
              !print*,'AA',ift,k,x,lutab(classes(k),ift),classes(k)
            ENDDO
          ENDDO

          IF(ilanduse.GE.3 .AND. ilanduse.LE.6 ) THEN !PCM
           atprop2(:,:) = agclassprop2(1:maxn_at,1:maxn_at) !PCM: kluge: assumes lutab2(iat3,iat3) = 100.0
           atharvest(:) = agclassprop(1:maxn_at) !PCM: kluge: assumes lutab2(iat3,iat3) = 100.0
!PCM: try simplified version above, first; comment out these lines
!           DO iat3=1,maxn_at
!           DO iat2=1,maxn_at
!            atprop2(iat3,iat2)=0.0d0
!            DO k3=1,maxn_at
!             DO k2=1,maxn_at
!              atprop2(iat3,iat2)=atprop2(iat3,iat2)+
!     &lutab2(agclasses(k3),iat3)*
!     &agclassprop2(agclasses(k3),agclasses(k2))/100.0d0
!             print *,'DD2',iat3,iat2,k3,k2,lutab2(agclasses(k3),iat3),
!     &agclasses(k3),
!     &agclassprop2(agclasses(k3),agclasses(k2)),atprop2(iat3,iat2)
!             ENDDO
!            ENDDO
!           ENDDO
!           ENDDO
          ENDIF

c
c Calculate the bare soil.
c
          IF ((ftprop(2).LE.100).AND.(ftprop(2).GE.0)) THEN
            ftprop(1)=100
            DO ift=2,nft
              ftprop(1)=ftprop(1)-ftprop(ift)
            ENDDO
          ENDIF

        ENDIF ! finished reading

        ! previous SDGVM method finished here, assumes step changes in land-use
        DO ift=1,nft
          cluse(ift,i) = ftprop(ift)
        ENDDO

        IF(ilanduse.GE.3 .and. ilanduse.LE.6 ) THEN
!         IF( MOD(yr0a+i-2,20) == 0 ) THEN 
!          WRITE(*,*)
!     &'   t    BARE    CITY     C3p     C4p  C3crop  C4crop',
!     &'     C3s     C4s   Ev_Bp   Ev_Np   Dc_Bp   Dc_Np',
!     &'   Ev_Bs   Ev_Ns   Dc_Bs   Dc_Ns'
!         ENDIF
!         write(*,'(I4,16F8.2)') yr0a+i-1, cluse(1:16,i)  
        ELSE
         IF(debug .EQV. .True.) THEN
          IF( MOD(yr0a+i-2,20) == 0 ) THEN 
           WRITE(*,*)
     &'   t    BARE    CITY      C3      C4  C3crop  C4crop',
     &'   Ev_Bl   Ev_Nl   Dc_Bl   Dc_Nl'
          ENDIF
          write(*,'(I4,10F8.2)') yr0a+i-1, cluse(1:10,i)  
         ENDIF
        ENDIF

        IF(ilanduse.GE.3 .AND. ilanduse.LE.6 ) THEN !PCM
         DO iat3=1,maxn_at
          cluseh(iat3,i)       = atharvest(iat3)
          DO iat2=1,maxn_at
           cluse2(iat3,iat2,i) = atprop2(iat3,iat2)
          ENDDO
         ENDDO
        ENDIF
      ENDDO !year loop #2

!      IF(prescr_fire .EQV. .TRUE.) THEN
!          WRITE(*,'(A,I3)')'EX_CLU:After year loop, fprob_prescrh'
!          DO jj=1,n_fire
!            WRITE(*,'(12F7.4)')fprob_prescrh(jj,1:12)
!          END DO
!      END IF


      ! TRENDY SDGVM method 
      !- assumes linear interpolation of land-use between years
      !specified in input dataset
      IF ((n-j1+1).ne.1) THEN
        IF(debug .EQV. .True.) THEN
          print*, 'Linear interpolation of dynamic Land-Cover fractions'
          print*, n,j1,n-j1+1
        ENDIF
        DO i=1,years(n)-yr0a-yr_offset
          IF ( (i.EQ.1).OR.((i+yr0a-1).EQ.(years(j1)-yr_offset)) ) THEN
!            print*, j, years(j1), years(j1+1)
            ij  = i
            ij1 = ij + years(j1+1) - years(j1)  
            j1  = j1 + 1
!            print*, i,j1,ij,ij1
          ELSE
!            print*, i,j1,ij,ij1,'INTERPOLATING in TIME'
!       for S3 experiment, it never reaches here
            DO ift=1,nft
              cluse(ift,i) = cluse(ift,ij) + 
     &  ( (real(i)-real(ij))/(real(ij1)-real(ij)) * 
     &  (cluse(ift,ij1) - cluse(ift,ij)) )
      !       if(ift.eq.10) then
      !        print*, cluse(ift,ij), (cluse(ift,ij1) - cluse(ift,ij)),
      !&cluse(ift,i) 
      !      endif 
            ENDDO
          !write(*,'(I4,12F8.2)') yr0a+i-1, cluse(1:12,i)  

            DO iat3=1,maxn_at
             DO iat2=1,maxn_at
              cluse2(iat3,iat2,i) = cluse2(iat3,iat2,ij) + 
     &  ( (real(i)-real(ij))/(real(ij1)-real(ij)) * 
     &  (cluse2(iat3,iat2,ij1) - cluse2(iat3,iat2,ij)) )
             ENDDO
            ENDDO

            DO iat3=1,maxn_at
              cluseh(iat3,i) = cluseh(iat3,ij) + 
     &  ( (real(i)-real(ij))/(real(ij1)-real(ij)) * 
     &  (cluseh(iat3,ij1) - cluseh(iat3,ij)) )
            ENDDO

          ENDIF
        ENDDO
      ENDIF

      do_fire_interp = .false. !this doesn't work now, so turn it off
      IF (do_fire_interp .and. prescr_fire) THEN
      ! TRENDY SDGVM method 
      !- assumes linear interpolation of land-use between years
      !specified in input dataset
       IF ((n_fire-j1_fire+1).ne.1) THEN
        IF(debug .EQV. .True.) THEN
          print*, 'Linear interpolation of dynamic Land-Cover fractions'
          print*, n_fire,j1_fire,n_fire-j1_fire+1
        ENDIF
        DO i=1,years_fire(n_fire)-yr0a-yr_offset
          IF ( (i.EQ.1).OR.((i+yr0a-1).EQ.    
     &          (years_fire(j1_fire)-yr_offset)) ) THEN
!            print*, j_fire, years(j1_fire), years(j1_fire+1)
            ij  = i
            ij1 = ij + years_fire(j1_fire+1) - years_fire(j1_fire)  
            j1_fire  = j1_fire + 1
!            print*, i,j1_fire,ij,ij1
          ELSE
!              print*, i,j1_fire,ij,ij1,'INTERPOLATING in TIME'
              fprob_prescrh(i,1:12) = fprob_prescrh(ij,1:12) + 
     &  ( (real(i)-real(ij))/(real(ij1)-real(ij)) * 
     &  (fprob_prescrh(ij1,1:12) - fprob_prescrh(ij,1:12)) )
          ENDIF
        ENDDO
       ENDIF
      ENDIF

!      IF(prescr_fire .EQV. .TRUE.) THEN
!          WRITE(*,'(A,I3)')'EX_CLU:After TRENDY SDGVM method:fprob'
!          DO jj=1,n_fire
!            WRITE(*,'(12F7.4)')fprob_prescrh(jj,1:12)
!          END DO
!      END IF
!      WRITE(*,'(A,I3)')'EX_CLU:After TRENDY SDGVM method:cluseh(1)'
!      DO jj=1,n
!        WRITE(*,'(F7.4)')cluseh(1,jj)
!      END DO
!      WRITE(*,'(A,I3)')'EX_CLU:After TRENDY SDGVM method:cluseh(2)'
!      DO jj=1,n
!        WRITE(*,'(F7.4)')cluseh(2,jj)
!      END DO
!      WRITE(*,'(A,I3)')'EX_CLU:After TRENDY SDGVM method:cluseh(3)'
!      DO jj=1,n
!        WRITE(*,'(F7.4)')cluseh(3,jj)
!      END DO



!!old method:
!!PCM      IF ((indx(2,2).EQ.1).OR.(indx(2,3).EQ.1).OR.(indx(3,2).EQ.1).OR.
!!PCM     &  (indx(3,3).EQ.1)) THEN
!!PCM4     IF (indx(3,3).EQ.1) THEN
!PCM new method:
!comment this out, since this indx was computed for various time slices 
!      IF ((indx(3,3).EQ.1).OR.(indx(3,4).EQ.1).OR.(indx(4,3).EQ.1).OR.
!     &  (indx(4,4).EQ.1)) THEN
!        l_lu = .TRUE.
!      ELSE
!        l_lu = .FALSE.
!      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
      SUBROUTINE lorc(du,lat,lon,latdel,londel,stmask,xx)
*----------------------------------------------------------------------*
      IMPLICIT NONE
      LOGICAL xx
      REAL*8 lat,lon,latdel,londel,del,latf,lon0
      CHARACTER outc(6200),stmask*1000
      INTEGER i,j,k,x(7),ians(43300),sum1,n,col,row,check,nrecl,du,blank
      INTEGER kode,ii

      IF (du.eq.1) THEN
        nrecl = 6172
      ELSE
        nrecl = 6172 + 1
      ENDIF

      OPEN(99,file=stmask(1:blank(stmask))//'/land_mask.dat',
     &form='formatted',recl=nrecl,access='direct',status='old',
     &iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Land sea mask.'
        WRITE(*,*) 'Either the file doesn''t exist or there is a record 
     & length miss match.'
        WRITE(*,*) 'Check that the correct DOS|UNIX switch is being use
     &d in the input file.' 
        WRITE(*,*) 'Land sea file :',stmask(1:blank(stmask))
        STOP
      ENDIF

      del = 1.0d0/60.0d0/2.0d0
      latf = 90.0d0 - del/2.0d0
      lon0 =-180.0d0 + del/2.0d0

      col = int((lon - lon0)/del + 0.5d0)
      row = int((latf - lat)/del + 0.5d0)

      n = min((latdel/del-1.0d0)/2.0d0,(londel/del-1.0d0)/2.0d0)

*----------------------------------------------------------------------*
* Check nearest pixel for land.                                        *
*----------------------------------------------------------------------*
      sum1 = 0
      READ(99,'(6172a)',rec=row) (outc(j),j=1,6172)
      DO j=1,6172
        CALL base72i(outc(j),x)
        DO k=1,7
          ians(7*(j-1)+k) = x(k)
        ENDDO
      ENDDO
      IF (ians(col).GT.0) sum1 = sum1 + 1
      check = 1

*----------------------------------------------------------------------*

      IF (n.GT.0) THEN
*----------------------------------------------------------------------*
* Check outward diagonals for land.                                    *
*----------------------------------------------------------------------*
        DO ii=1,n/4
          i = 4*ii
          READ(99,'(6172a)',rec=row+i) (outc(j),j=1,6172)
          DO j=1,6172
            CALL base72i(outc(j),x)
            DO k=1,7
              ians(7*(j-1)+k) = x(k)
            ENDDO
          ENDDO
          IF (ians(col+i).GT.0) sum1 = sum1 + 1
          IF (ians(col-i).GT.0) sum1 = sum1 + 1

          READ(99,'(6172a)',rec=row-i) (outc(j),j=1,6172)
          DO j=1,6172
            CALL base72i(outc(j),x)
            DO k=1,7
              ians(7*(j-1)+k) = x(k)
            ENDDO
          ENDDO
          IF (ians(col+i).GT.0) sum1 = sum1 + 1
          IF (ians(col-i).GT.0) sum1 = sum1 + 1
          check = check + 4
        ENDDO
*----------------------------------------------------------------------*
      ENDIF

      CLOSE(99)

      IF (real(sum1).GE.(check+1)/2) then
        xx = .true.
      ELSE
        xx = .false.
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE BASE72I                          *
*                          ******************                          *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE base72i(c,x)
*----------------------------------------------------------------------*
      CHARACTER c
      INTEGER x(7),i

      i=ichar(c)-100
      x(1)=i/64
      i=i-x(1)*64
      x(2)=i/32
      i=i-x(2)*32
      x(3)=i/16
      i=i-x(3)*16
      x(4)=i/8
      i=i-x(4)*8
      x(5)=i/4
      i=i-x(5)*4
      x(6)=i/2
      i=i-x(6)*2
      x(7)=i


      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          FUNCTION n7                                 *
*                          ***********                                 *
*                                                                      *
* Returns the value of a seven digit binary number given as a seven    *
* dimensional array of 0's and 1's                                     *
*                                                                      *
*----------------------------------------------------------------------*
      FUNCTION n7(x)
*----------------------------------------------------------------------*
      INTEGER n7,x(7)

      n7 = 64*x(1)+32*x(2)+16*x(3)+8*x(4)+4*x(5)+2*x(6)+x(7)
 

      RETURN
      END

*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE BI_LIN                           *
*                          *****************                           *
*                                                                      *
* Performs bilinear interpolation between four points, the normalised  *
* distances from the point (1,1) are given by 'xnorm' and 'ynorm'.     *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE BI_LIN(xx,indx,xnorm,ynorm,ans)
*----------------------------------------------------------------------*
      REAL*8 xx(4,4),xnorm,ynorm,ans,av
      INTEGER indx(4,4),iav,ii,jj
      LOGICAL use_bilinear

!      use_bilinear = .False. !use nearest-neighbor sampling instead
      use_bilinear = .TRUE. !use nearest-neighbor sampling instead

*----------------------------------------------------------------------*
* Fill in averages if necessary.                                       *
*----------------------------------------------------------------------*
      DO ii=2,3
        DO jj=2,3
          IF (indx(ii,jj).NE.1) THEN
            av = 0.0d0
            iav = 0
            IF (indx(ii+1,jj).EQ.1) THEN
              av = av + xx(ii+1,jj)
              iav = iav + 1
            ENDIF 
            IF (indx(ii-1,jj).EQ.1) THEN
              av = av + xx(ii-1,jj)
              iav = iav + 1
            ENDIF 
            IF (indx(ii,jj+1).EQ.1) THEN
              av = av + xx(ii,jj+1)
              iav = iav + 1
            ENDIF 
            IF (indx(ii,jj-1).EQ.1) THEN
              av = av + xx(ii,jj-1)
              iav = iav + 1
            ENDIF
            IF (indx(ii+1,jj+1).EQ.1) THEN
              av = av + xx(ii+1,jj+1)
              iav = iav + 1
            ENDIF
            IF (indx(ii-1,jj-1).EQ.1) THEN
              av = av + xx(ii-1,jj-1)
              iav = iav + 1
            ENDIF
            IF (indx(ii+1,jj-1).EQ.1) THEN
              av = av + xx(ii+1,jj-1)
              iav = iav + 1
            ENDIF
            IF (indx(ii-1,jj+1).EQ.1) THEN
              av = av + xx(ii-1,jj+1)
              iav = iav + 1
            ENDIF
            IF (iav.GT.0) THEN
              xx(ii,jj) = av/real(iav)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      IF( use_bilinear .EQV. .TRUE. ) THEN 
*----------------------------------------------------------------------*
* Bilinear interpolation.                                              *
*----------------------------------------------------------------------*
        ans = xx(2,2)*(1.0d0-xnorm)*(1.0d0-ynorm) + 
     &        xx(3,2)*xnorm*(1.0d0-ynorm) + 
     &        xx(2,3)*(1.0d0-xnorm)*ynorm + 
     &        xx(3,3)*xnorm*ynorm

      ELSE
*----------------------------------------------------------------------*
* Nearest grid cell.                                                   *
*----------------------------------------------------------------------*
        !ans = xx(int(xnorm+2.5d0),int(ynorm+2.5d0))
        ans = xx(3,3)
*----------------------------------------------------------------------*
      ENDIF


      RETURN
      END


*----------------------------------------------------------------------*
*                                                                      *
*                          SUBROUTINE AVE4                             *
*                          *****************                           *
*                                                                      *
* Performs averaging of up to four points,                             *
* the normalised distances from the point (1,1) are given by 'xnorm'   *
* and 'ynorm'.                                                         *
*                                                                      *
*----------------------------------------------------------------------*
      SUBROUTINE AVE4(xx,indx,xnorm,ynorm,ans)
*----------------------------------------------------------------------*
      REAL*8 xx(4,4),xnorm,ynorm,ans,av,av4
      INTEGER indx(4,4),iav,ii,jj,iav4
      LOGICAL use_bilinear

      use_bilinear = .False. !use nearest-neighbor sampling instead

*----------------------------------------------------------------------*
* Fill in averages if necessary.                                       *
*----------------------------------------------------------------------*
!      DO ii=2,3
!        DO jj=2,3
!          IF (indx(ii,jj).NE.1) THEN
!            av = 0.0d0
!            iav = 0
!            IF (indx(ii+1,jj).EQ.1) THEN
!              av = av + xx(ii+1,jj)
!              iav = iav + 1
!            ENDIF 
!            IF (indx(ii-1,jj).EQ.1) THEN
!              av = av + xx(ii-1,jj)
!              iav = iav + 1
!            ENDIF 
!            IF (indx(ii,jj+1).EQ.1) THEN
!              av = av + xx(ii,jj+1)
!              iav = iav + 1
!            ENDIF 
!            IF (indx(ii,jj-1).EQ.1) THEN
!              av = av + xx(ii,jj-1)
!              iav = iav + 1
!            ENDIF
!            IF (indx(ii+1,jj+1).EQ.1) THEN
!              av = av + xx(ii+1,jj+1)
!              iav = iav + 1
!            ENDIF
!            IF (indx(ii-1,jj-1).EQ.1) THEN
!              av = av + xx(ii-1,jj-1)
!              iav = iav + 1
!            ENDIF
!            IF (indx(ii+1,jj-1).EQ.1) THEN
!              av = av + xx(ii+1,jj-1)
!              iav = iav + 1
!            ENDIF
!            IF (indx(ii-1,jj+1).EQ.1) THEN
!              av = av + xx(ii-1,jj+1)
!              iav = iav + 1
!            ENDIF
!            IF (iav.GT.0) THEN
!              xx(ii,jj) = av/real(iav)
!            ENDIF
!          ENDIF
!        ENDDO
!      ENDDO
!
      IF( use_bilinear .EQV. .TRUE. ) THEN 
*----------------------------------------------------------------------*
* Bilinear interpolation.                                              *
*----------------------------------------------------------------------*
        ans = xx(2,2)*(1.0d0-xnorm)*(1.0d0-ynorm) + 
     &        xx(3,2)*xnorm*(1.0d0-ynorm) + 
     &        xx(2,3)*(1.0d0-xnorm)*ynorm + 
     &        xx(3,3)*xnorm*ynorm

      ELSE
*----------------------------------------------------------------------*
* Nearest grid cell.                                                   *
*----------------------------------------------------------------------*
        !ans = xx(int(xnorm+2.5d0),int(ynorm+2.5d0))
        !ans = xx(3,3)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* AVE of up to 4 grid cells                                            *
*----------------------------------------------------------------------*
         av4 = 0.0d0
         iav4 = 0
         IF(indx(3,3).EQ.1) THEN
            av4 = av4 + xx(3,3)
            iav4 = iav4 + 1
         ENDIF
         IF(indx(3,4).EQ.1) THEN
            av4 = av4 + xx(3,4)
            iav4 = iav4 + 1
         ENDIF
         IF(indx(4,3).EQ.1) THEN
            av4 = av4 + xx(4,3)
            iav4 = iav4 + 1
         ENDIF
         IF(indx(4,4).EQ.1) THEN
            av4 = av4 + xx(4,4)
            iav4 = iav4 + 1
         ENDIF
         IF (iav4.GT.0) THEN
            ans = av4/real(iav4)
         ELSE
            ans = -999.0
         ENDIF
      ENDIF


      RETURN
      END


