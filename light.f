*----------------------------------------------------------------------*
*                          FUNCTIONS                                   *
*----------------------------------------------------------------------*
*                          DAYLENGTH                                   *
*----------------------------------------------------------------------*
      FUNCTION dayl(lat,day)
*----------------------------------------------------------------------*
      REAL*8 dayl,lat,del,has,tem,conv!,dawn_angle
      INTEGER day
      PARAMETER (conv = 1.74532925E-2)

      !declination angle
      del=-23.4d0*cos(conv*360.0d0*(day + 10.0d0)/365.0d0)
      !tem is the acos of the dawn angle
      tem=-tan(lat*conv)*tan(del*conv)
      !CALL sunrise_angle(day,lat,dawn_angle,tem)

      IF (abs(tem).LT.1.0d0) THEN
        ! has is the hour angle of sunrise in degrees
        has  = acos(tem)/conv
        !CALL sunrise_angle(day,lat,dawn_angle)
        !has  = dawn_angle/conv
        !calculate daylength, 
        ! - 2has is daylength in degrees, 15 is the angular velocity of Earth deg/hr
        dayl = 2.0d0*has/15.0d0
      ELSEIF (tem.GT.0.0d0) THEN
        dayl =  0.0d0
      ELSE
        dayl = 24.0d0
      ENDIF

      !print*, day,dayl,tem,dawn_angle,has

      RETURN
      END


*----------------------------------------------------------------------*
*                          DAWN ANGLE                                  *
*----------------------------------------------------------------------*
      !SUBROUTINE sunrise_angle(day,lat,angle,tem)

      !REAL*8    lat,del,angle,tem,conv
      !INTEGER   day
      !PARAMETER (conv = 1.74532925E-2)

      ! calculates the hour angle at dawn 
      ! - in radians
 
      ! declination angle - in degrees
      !del   = -23.4d0*cos(conv*360.0d0*(day + 10.0d0)/365.0d0)
      ! tem is the cos of the dawn angle
      !tem   = -tan(lat*conv)*tan(del*conv)
      !angle = acos(tem)

      !END SUBROUTINE

*----------------------------------------------------------------------*
*                          PHOTON FLUX DENSITY                         *
*----------------------------------------------------------------------*
      SUBROUTINE pfd(lat,day,hrs,cloud,direct,diffuse,total,swr,
     & read_par,subd_par,t,total_t,calc_zen,cos_zen)
*----------------------------------------------------------------------*
      REAL*8    lat,hrs,del,rlat,toa,pi,conv,cloud,diffprop,alpha,beta
      REAL*8    sigma,dawn_angle,coscst,direct,diffuse,total,clearness
      REAL*8    swr,ts,h,toa_h,solar_const,tem,dJ_R,dJ_K,cos_zen
      INTEGER   day,read_par,subd_par,t,total_t,calc_zen
      PARAMETER (conv = 1.74532925E-2,pi = 3.1415927)
      PARAMETER (solar_const = 1370.0d0)

      !print*, hrs
      IF (hrs.GT.1E-6) THEN
        ! declination angle
        del  = -23.4d0*cos(conv*360.0d0*(day+10.0d0)/365.0d0)
        del  = del*conv
        rlat = lat*conv

        !calculate dawn angle
        dawn_angle = (hrs/2.0d0)*(2.0d0*pi/24.0d0)
        !CALL sunrise_angle(day,lat,dawn_angle,tem)
        !print*, dawn_angle

        !calculate daytime mean cos of the zenith angle  
        cos_zen = 1/dawn_angle*( 
     &    sin(rlat)*sin(del)*dawn_angle + 
     &    cos(rlat)*cos(del)*sin(dawn_angle) )

        !top of atmosphere irradiance
        ! - watts/m2 average over the day 
        ! - between sun rise and sun set
        !toa = solar_const * cos_zen ! daytime mean irradiance (W/m2)
        toa = solar_const *         ! solar constant (W/m2)
     &    24.0d0/(pi*hrs)*(         ! average over the day (daylength)
     &    sin(rlat)*sin(del)*dawn_angle + 
     &    cos(rlat)*cos(del)*sin(dawn_angle) )

c     calculate the mean daylight top of canopy (total) shortwave radiation
        IF(read_par.eq.1) THEN
          !read from the input dataset as mean SWR over 24hrs
          ! - scale to mean over daylight hours only
          total = 24.0d0*swr/hrs   
        ELSE
          !calculate based on top of atmosphere value and cloud cover
          ! - also calculated from Forest ETP 
          coscst = 1.0d0-1.3614d0*cos(rlat)
          alpha  = 18.0d0 - 64.884d0*coscst
          alpha  = alpha*4.1842d0*1.0d4     ! convert from cal/cm to J/m2
          alpha  = alpha/(24.0d0*3600.0d0)  ! convert into w/m2  
          beta   = 0.682d0 - 0.3183d0*coscst
          sigma  = 0.02d0*log(max(cloud,0.001d0))+0.03259d0

          total  = toa*(beta-sigma*cloud*10.0d0)-alpha
        ENDIF

        IF(subd_par.eq.1) THEN
          !sub-daily variation in SWR 

          !calculate using hour angle (h) 
          ! - scale by (instantaneous insolation at h)/(mean daylight insolation) law of cosines
          ! - when t = 0, h is 0 and is solar noon
          ! - when t = total_t, h is the dawn angle i.e. sunrise
          h       = dawn_angle * real(t)/real(total_t)
          ! calculate cos of the zenith angle
          cos_zen = sin(rlat)*sin(del) + cos(rlat)*cos(del)*cos(h)
          !calculate instantaneous insolation at hour angle h
          toa_h   = solar_const * cos_zen 
          !toa_h   = solar_const * ( sin(rlat)*sin(del) + 
      !& cos(rlat)*cos(del)*cos(h) )
!          if(day.eq.160) print*, h,toa,toa_h,toa_h/toa
          !scale SWR (total) from daytime mean to instantaneous 
          ! - scale by the ratio of toa instantaneous insolation to daytime mean
          if(toa.gt.0.0) then
            total = total * toa_h/toa
          else
            total = 0.0
          endif
!          if(day.eq.160) print*, total 
          !reset top of atmosphere to instantaneous value
          toa = toa_h

        ENDIF

        IF (total.LT.0.0d0) total=0.0d0
        if(toa.le.0.0) then
          clearness = 0.0
        else
          clearness = total/toa
        endif

        !calculate diffuse irradiance (from Spitters etal 1986)
        IF(subd_par.eq.1) THEN
          ! hourly data
          dJ_R = 0.847d0 - 1.61d0*cos_zen + 1.04d0*cos_zen**2
          dJ_K = (1.47d0-dJ_R)/1.66d0
          IF(clearness.LT.0.22) THEN
            diffprop = 1.0
          ELSE IF(clearness.LT.0.35) THEN
            diffprop = 1.0-6.4*(clearness-0.22)**2
          ELSE IF(clearness.LT.dJ_K) THEN
            diffprop = 1.47-1.66*clearness
          ELSE
            diffprop = dJ_R
          ENDIF
        ELSE
          ! daily data
          IF(clearness.LT.0.07) THEN
            diffprop = 1.0
          ELSE IF(clearness.LT.0.35) THEN
            diffprop = 1.0-2.3*(clearness-0.07)**2
          ELSE IF(clearness.LT.0.75) THEN
            diffprop = 1.33-1.46*clearness
          ELSE
            diffprop = 0.23d0
          ENDIF
        ENDIF

        diffuse = total*diffprop

c PAR is about 48% of the total irradiance. A better formula could be 
c found.

        total   = 0.48d0*total
        diffuse = 0.48d0*diffuse ! not sure this is really right
                                 ! ... (Ghislain)

c convert to mol/m2/sec from W/m2

        total   = 4.6d-6*total
        diffuse = 4.6d-6*diffuse
        direct  = total - diffuse

      ELSE

        total   = 0.0d0
        diffuse = 0.0d0
        direct  = 0.0d0
      ENDIF

      ! set zenith angle to zero if requested (originaly SDGVM default)
      if(calc_zen.eq.0) cos_zen = 1d0
      
      !print*, hrs, swr, total!, direct, diffuse

      RETURN
      END







*----------------------------------------------------------------------*
*                          PHOTON FLUX DENSITY                         *
*----------------------------------------------------------------------*
      FUNCTION pfd_without_cloud(lat,day,hrs)
*----------------------------------------------------------------------*
      REAL*8 pfd_without_cloud
      REAL*8 lat,hrs,del,rlat,r,sinbt,pi,conv,cloud,st,sdiff,trans
      INTEGER day
      PARAMETER(conv = 1.74532925E-2, pi = 3.1415927)

      cloud = 0.50d0

      IF (hrs.GT.1E-6) THEN
        del  = -23.4d0*cos(conv*360.0d0*(day+10.0d0)/365.0d0)
        del  = del*conv
        rlat = lat*conv

        trans = 0.3d0 + (0.35d0*(1.0d0 - cloud))

* solar radiation in watts/m2 at midday
        sinbt = sin(rlat)*sin(del) + cos(rlat)*cos(del)

        IF (sinbt.GT.1E-6) THEN
          r = 1370.0d0*trans**(1.0d0/sinbt)*sinbt

          sdiff = 0.3d0*(1.0d0 - trans**(1.0d0/sinbt))*1370.0d0*sinbt

          st = r + sdiff

          pfd_without_cloud = 
     &         2.0d0*st/pi*(0.5d0 + cos(rlat)/2.0d0)*2.1d0/1.0e+6

        ELSE
         pfd_without_cloud = 0.0d0
        ENDIF
      ELSE
        pfd_without_cloud = 0.0d0
      ENDIF


      RETURN
      END


*----------------------------------------------------------------------*
*                          PHOTON FLUX DENSITY                         *
*----------------------------------------------------------------------*
      FUNCTION pfd2(lat,day,hrs)
*----------------------------------------------------------------------*
      REAL*8 pfd2,lat,hrs,del,rlat,r,sinbt,pi,conv,rad
      INTEGER day
      PARAMETER(conv = 1.74532925E-2, pi = 3.1415927)

      IF (hrs.GT.1E-6) THEN
	del  = -23.4d0*cos(conv*360.0d0*(day+10.0d0)/365.0d0)
	del  = del*conv
	rlat = lat*conv

* solar radiation in watts/m2 at midday
        sinbt = sin(rlat)*sin(del) + cos(rlat)*cos(del)
        r = 1360.0d0*0.7d0**1.32d0*sinbt

* convert to Langleys/day
        r = r/0.485d0

* account for cloud cover (assume 0.0)
        rad = r*(-0.21843d0*20.0d0 + 58.94408d0)/100.0d0

* solve for daily maximum (mol/m2/sec)
        pfd2 = (rad*41868.0d0*pi)/(hrs*2.0d0*3600.0d0)
        pfd2 = pfd2*0.5d0*4.255E-6
      ELSE
	pfd2 = 0.0d0
      ENDIF


      RETURN
      END

*----------------------------------------------------------------------*
*                          PHOTON FLUX DENSITY                         *
*                                                                      *
* From Total Downward Shortwave at the Surface (TDSS) - for use with   *
* GCM data.                                                            *
*----------------------------------------------------------------------*
      FUNCTION pfds(tdss,hrs)
*----------------------------------------------------------------------*
      REAL*8 pfds,tdss,hrs,pi
      PARAMETER (pi = 3.1415927)

      IF (hrs.GT.1E-6) THEN
*   solve for daily maximum (mol/m2/sec)
        pfds = (tdss*pi*24.0d0)/(hrs*2.0d0)*4.255E-6
*        pfds = pfds*0.5d0*4.255E-6
      ELSE
	pfds = 0.0
      ENDIF


      RETURN
      END
