
      SUBROUTINE GOUDRIAANSLAW_old(lyr,lai,beamrad,diffrad,
     &fsunlit,qsunlit,fshade,qshade,can_clump,cos_zen)
      IMPLICIT NONE
      
      REAL*8 lyr,lai,beamrad,diffrad
      REAL*8 fsunlit,qsunlit,fshade,qshade
      
      REAL*8 soilalbedo,leafscattering,canopyalbedo
      REAL*8 albedobeam,albedodiff
      REAL*8 m,can_clump,cos_zen
      REAL*8 kbeam,kdiff,kbeamstar

c     value from Wang
      soilalbedo=0.15d0
c     leaf transmittance and reflectance are equal. Then albedo is the double
      leafscattering=2.0d0*0.075d0

c     extinction coefficent of the direct beam
      kbeam=0.5d0  * can_clump / cos_zen
c     extinction coefficient of diffuse radiations
      kdiff=0.66d0 
c     I found somewhere an approximate relation: kbeam=0.75*kdiff

      m=sqrt(1.0d0-leafscattering)
      kbeamstar=m*kbeam

      canopyalbedo=2.0d0*kbeam/(kbeam+kdiff)*(1.0d0-m)/(1.0d0+m)
      albedobeam=canopyalbedo+(soilalbedo-canopyalbedo)
     &     *exp(-kbeam*lai)
c     ... hum, I don't know if albedodiff is really different or not... so equal
      albedodiff=albedobeam

      qshade=(1.0d0-albedodiff)*diffrad*kdiff*exp(-kdiff*lyr)
     & +(1.0d0-albedobeam)*beamrad*kbeam*exp(-kbeam*lyr)
     & -(1.0d0-leafscattering)*beamrad*kbeamstar*exp(-kbeamstar*lyr)
      IF (qshade.LT.0.0d0) qshade=0.0d0

      qsunlit=(1.0d0-leafscattering)*kbeamstar*beamrad+qshade


      fsunlit=exp(-kbeam*lyr)
      fshade=1-fsunlit
      

      END


      SUBROUTINE GOUDRIAANSLAW(lyr,lai,beamrad,diffrad,
     &fsunlit,qsunlit,fshade,qshade,can_clump,cos_zen,s070607,gold)
      
      IMPLICIT NONE
      
      REAL*8 lyr,lai,beamrad,diffrad
      REAL*8 fsunlit,qsunlit,fshade,qshade
      REAL*8 rhosoil,leafscattering,rhoh
      REAL*8 rhobeam,rhodiff,rhobeam_can,rhodiff_can
      REAL*8 m,can_clump,cos_zen,Gbeam,canopyalbedo
      REAL*8 kbeam,kdiff,kbeamprime,kdiffprime
      INTEGER s070607
      LOGICAL gold

      if(.not.((s070607.eq.1).or.(gold))) then

        ! Goudriaan's Law
        ! as described in Walker ... from Spitters 1986 and Wang 2003 Func.Plant Biol.     
        ! all below values are for visible wavelengths 

        ! extinction coefficents of direct & diffuse radiation 
        ! these account for both canopy clumping and solar zenith angle, also inform by Bodin & Franklin 2012 GMD  
        ! transmittance equals reflectance
        leafscattering = 2d0*0.075d0
        m          = sqrt(1.0d0-leafscattering)
        Gbeam      = 0.5d0 * can_clump
        kbeam      = Gbeam / cos_zen
        kbeamprime = m*kbeam
        ! kdiff is Gbeam * can_clump / cos_zen integrated over zenith angle 0 to pi/2  
        kdiff      = 0.8d0 * can_clump 
        kdiffprime = m*kdiff

        ! calculate albedos
        rhosoil     = 0.15d0
        rhoh        = (1d0-m)/(1d0+m)
        rhobeam_can = rhoh*2d0*kbeam/(kbeam+kdiff)
        rhodiff_can = 4*Gbeam*rhoh*( 
     & Gbeam*(log(Gbeam)-log(Gbeam+kdiff))/kdiff**2 
     & + 1/kdiff)
        rhobeam     = rhobeam_can + (rhosoil-rhobeam_can)
     & *exp(-2d0*kbeamprime*lai)
        rhodiff     = rhodiff_can + (rhosoil-rhodiff_can)
     & *exp(-2d0*kdiffprime*lai)

       ! calculate absorbed direct & diffuse radiation  
       qshade = (1d0-rhodiff)*diffrad*kdiffprime*exp(-kdiffprime*lyr)
     & + (1d0-rhobeam)*beamrad*kbeamprime*exp(-kbeamprime*lyr)
     & - (1d0-leafscattering)*beamrad*kbeam*exp(-kbeam*lyr)

        IF (qshade.LT.0d0) qshade = 0d0
        qsunlit = (1d0-leafscattering)*kbeam*beamrad + qshade

        ! calculate fraction sunlit vs shaded leaves
        fsunlit = exp(-kbeam*lyr)
        fshade  = 1d0 - fsunlit
      
      else 

c       value from Wang
        rhosoil=0.15d0
c       leaf transmittance and reflectance are equal. Then albedo is the double
        leafscattering=2.0d0*0.075d0

c       extinction coefficent of the direct beam
        kbeam=0.5d0  * can_clump / cos_zen
c       extinction coefficient of diffuse radiations
        kdiffprime=0.66d0 
c       I found somewhere an approximate relation: kbeam=0.75*kdiff

        m=sqrt(1.0d0-leafscattering)
        kbeamprime=m*kbeam

        rhobeam_can=2.0d0*kbeam/(kbeam+kdiffprime)*(1.0d0-m)/(1.0d0+m)
        rhobeam=rhobeam_can+(rhosoil-rhobeam_can)
     &   *exp(-kbeam*lai)
c       ... hum, I don't know if albedodiff is really different or not... so equal
        rhodiff=rhobeam

       qshade=(1.0d0-rhodiff)*diffrad*kdiffprime*exp(-kdiffprime*lyr)
     & +(1.0d0-rhobeam)*beamrad*kbeam*exp(-kbeam*lyr)
     & -(1.0d0-leafscattering)*beamrad*kbeamprime*exp(-kbeamprime*lyr)
        IF (qshade.LT.0.0d0) qshade=0.0d0

        qsunlit=(1.0d0-leafscattering)*kbeamprime*beamrad+qshade


        fsunlit=exp(-kbeam*lyr)
        fshade=1-fsunlit

      endif

      END
      
      
      
      

      SUBROUTINE BEERSLAW(lyr,beamrad,diffrad,
     &     fsunlit,qsunlit, fshade,qshade)
      IMPLICIT NONE
      REAL*8 lyr,beamrad,diffrad
      REAL*8 fsunlit,qsunlit,fshade,qshade
      REAL*8 kbeam,kdiff

c     extinction coefficent of the direct beam
      kbeam=0.5d0
c     extinction coefficient of diffuse radiations
      kdiff=0.66d0     
      
      fshade=1.0d0
      fsunlit=0

      qsunlit=0.0d0
      qshade=beamrad*kbeam*exp(-kbeam*lyr)+diffrad*kdiff*exp(-kdiff*lyr)

      END


      SUBROUTINE SUNFLECT(lyr,beamrad,diffrad,
     &     fsunlit,qsunlit, fshade,qshade)
      IMPLICIT NONE
      REAL*8 lyr,beamrad,diffrad
      REAL*8 fsunlit,qsunlit,fshade,qshade
      REAL*8 kbeam,kdiff

c     extinction coefficent of the direct beam
      kbeam=0.5d0
c     extinction coefficient of diffuse radiations
      kdiff=0.66d0     
      
      fsunlit=0.5d0
      fshade=0.5d0

      qsunlit=beamrad*kbeam*exp(-kbeam*lyr)+
     &     diffrad*kdiff*exp(-kdiff*lyr)
      qshade=qsunlit


      END
