*----------------------------------------------------------------------*
*                       SDGVM0 VERSION 140928                          *
*----------------------------------------------------------------------*
      IMPLICIT NONE

      CHARACTER stver*1000
      PARAMETER(stver = 'VERSION 140928')

      INCLUDE 'array_dims.inc'
      INCLUDE 'param.inc'

      REAL*8 defaulttopsl
      PARAMETER(defaulttopsl = 5.0d0)
      INTEGER, PARAMETER :: NX = 720, NY = 360
      INTEGER, PARAMETER :: NE = 10, maxn_at = 7, NS = 16

      REAL*8 lat,lon,dep,ca(12,31),npp(maxnft),lai(maxnft),evp(maxnft)
      REAL*8 gpp(maxnft),sresp(maxnft),evt(maxnft),soilt,grassrc,resp
      REAL*8 tmp(12,31),prc(12,31),hum(12,31),cld(12),latdel,londel
      REAL*8 leafper,stemper,rootper,avnpp,avgpp,avlai,avdof,co20,co2f
      REAL*8 avrof,infix,nppsold,avnppst,co2const,ic0(8),in0(8),iminn(3)
      REAL*8 avnleaf,avvcmax,avjmax,avleaf_nit,avsla,sl,hrs
      REAL*8 swcold,swcnew,is1,is2,is3,is4,isn,ilsn,tc0(8),tn0(8),sum1
      REAL*8 tminn(3),ts1,ts2,ts3,ts4,tsn,tlsn,oscale,yield(maxnft)
      REAL*8 c0(8,maxnft),n0(8,maxnft),minn(3,maxnft),c0v(8),n0v(8)
      REAL*8 minnv(3),s1(maxnft),s2(maxnft),s3(maxnft),s4(maxnft)
      REAL*8 sn(maxnft),lsn(maxnft),tmin,prctot,npl(maxnft),wtswc,c3old
      REAL*8 nps(maxnft),npr(maxnft),dof(maxnft),rof(maxnft),mnthprc(12)
      REAL*8 co2(maxyrs,12,31),maxcov,maxbio,ftprop(maxnft),tmem(200)
      REAL*8 cov(maxage,maxnft),bio(maxage,2,maxnft),ppm(maxage,maxnft)
      REAL*8 hgt(maxage,maxnft),wdt(maxage,maxnft),leafdp(3600,maxnft)
      REAL*8 stemdp(1000,maxnft),rootdp(1000,maxnft),leafv(3600),c4old
      REAL*8 stemv(1000),rootv(1000),yeartmp,yearprc,yearhum,barerc
      REAL*8 laimax(maxnft),laimnth(12,maxnft),ht(maxnft),avtrn,yrpet
      REAL*8 trn(maxnft),slc(maxnft),rlc(maxnft),sln(maxnft),rln(maxnft)
      REAL*8 mnthtmp(12),firec,dslc,drlc,dsln,drln,lch(maxnft),avlch
      REAL*8 ftgr0(maxnft),ftgrf(maxnft),mnthhum(12),ccheck,iadj,jadj
      REAL*8 bioo(maxnft),covo(maxnft),avevt,avsresp,ftppm0(maxnft)
      REAL*8 tsoilc,tsoiln,soilc(maxnft),soiln(maxnft),isoilc,isoiln
      REAL*8 sumbio,ans1,ftstmx(maxnft),leaflit(maxnft),stemlit(maxnft)
      REAL*8 rootlit(maxnft),ftwd(maxnft),ftxyl(maxnft),ftpd(maxnft)
      REAL*8 ftsla(maxnft),ftcov(maxnft),lon0,lonf,ftrat(maxnft),kd,kx
      REAL*8 input_ftsla(maxnft)
      REAL*8 ftvna(maxnft),ftvnb(maxnft),ftjva(maxnft),ftjvb(maxnft)
      REAL*8 ftg0(maxnft),ftg1(maxnft),amax(maxnft),vcmax_from_amax
      REAL*8 stembio,rootbio,sum,solcoo,biotoo,lutab(255,100),awl(4)
      REAL*8 nppstore(maxnft),ftlmor(maxnft),bioleaf(maxnft),nci(4)
      REAL*8 nppstoreold(maxnft),xlatf,xlatres,xlon0,xlonres,lat0,latf
      REAL*8 lat_lon(2000000,2),nppstorx(maxnft),avppm,avpet,fl(14)
      REAL*8 fpet(maxnft),SETARANDOM,temp_lat,temp_lon,f0(14,maxnft)
      REAL*8 nppstor2(maxnft),maxlai(maxnft),ftbb0(maxnft),wtfc,wtwp
      REAL*8 ftbbmax(maxnft),ftbblim(maxnft),ftsslim(maxnft),nleaf,cal
      REAL*8 flow1(maxnft),flow2(maxnft),h2o,adp(4),sfc(4),sw(4),sswc(4)
      REAL*8 leafnpp(maxnft),stemnpp(maxnft),rootnpp(maxnft),srespm,lchm
      REAL*8 evapm(12,maxnft),tranm(12,maxnft),roffm(12,maxnft),avyield
      REAL*8 photm(12,maxnft),avmnpet(maxnft),avmnt,avmnppt,tc,ts,tsi
      REAL*8 bulk,nupc,nfix,daygpp,evap,tran,roff,pet,yrtran,yrevap
      REAL*8 yrroff,daily_out(douts,maxnft,12,31),ans(12,31),interc,evbs
      REAL*8 soil_chr(10),leafold,stemold,rootold,ftagh(maxnft)
      REAL*8 cluse(maxnft,maxyrs),soil_chr2(10),fprob,ftmix(maxnft)
      REAL*8 oldlai,qdirect,qdiff,sumadp(360,maxnft),petm(12,maxnft)
      REAL*8 sm_trig(30,maxnft),smtrig(30),topsl,suma(360),dayra
      REAL*8 tsumam(maxnft),stemfr(maxnft),lmor_sc(3600,maxnft)
      REAL*8 wilt,field,sat,orgc,lflitold,fpr,daynpp,canga,ga
      REAL*8 soilCtoN,leaf_nit,vcmax(12,maxnft),jmax(12,maxnft)
      REAL*8 tassim,tgs,tci,nup_rate,pnlc(12,maxnft),enzs(12,maxnft)
      REAL*8 vtleaf_n(12,30),qdirsum,qdifsum,mnthp_days(12)             
      REAL*8 gsm(maxnft),tleaf_n,tleaf_p,respsum,cbal,leafg,f_amax
      REAL*8 soilp,soilcn_init,masoilcnv(3),soilcn,masoilcn,masoilcnr
      REAL*8 soilpr,soilp_init,kg(maxnft),kg_beta,ftcan_clump(maxnft)
      REAL*8 map_clump,w_scalar,t_scalar,leafv_sum,stemv_sum,rootv_sum
      REAL*8 aprc_dryqv(10),aprc_dryq,yearprcdryq,prc_week(52),prcq(52)
      REAL*8 matvar,aprc_rel,a2,b2,flulccc
      REAL*8 jmax_int(maxnft),jmax_int_er(maxnft)
      REAL*8 jmax_ci_low,jmax_ci_high
      REAL*8 ftToptV(maxnft),ftHaV(maxnft),ftHdV(maxnft)
      REAL*8 ftToptJ(maxnft),ftHaJ(maxnft),ftHdJ(maxnft)
      REAL*8 jmax_slope(maxnft),jmax_slope_er(maxnft)
      REAL*8 sla_ci_low,sla_ci_high
      REAL*8 sla_int(maxnft),sla_int_er(maxnft)
      REAL*8 sla_slope(maxnft),sla_slope_er(maxnft)
      REAL*8 atprop2(maxn_at,maxn_at)
      REAL*8 atharvest(maxn_at)
      REAL*8 cluse2(maxn_at,maxn_at,maxyrs)
      REAL*8 cluseh(maxn_at,maxyrs)
      REAL*8 sum_cov_test(maxnft)


      INTEGER read_clump,hw_j,cstype,calc_zen,phen_cor,pft_nflds
      INTEGER mswitch,subd_par,switch3,no_slw_lim
      INTEGER soilcn_map,soilp_map,vcmax_type,ncalc_type,read_par,ttype
      INTEGER sites,cycle,yr0,yrf,snp_no,snpshts(1000),dschill(maxnft)
      INTEGER ftbbm(maxnft),ftssm(maxnft),ftsss(maxnft),n_fields
      INTEGER snp_year,ftlls(maxnft),ftsls(maxnft),ftrls(maxnft),day,d
      INTEGER isite,ntags,du,ii,otagsn(douts),otagsnft(douts),sit_grd
      INTEGER ftmor(maxnft),ftc3(maxnft),nft,site0,sitef,nat_map(8)
      INTEGER ilanduse,siteno,iofn,iofnft,iofngft,recl1
      INTEGER icontinuouslanduse,ftphen(maxnft),ftdth(maxnft),kode
      INTEGER i,j,k,l,m,ft,s,w,w1,f
      INTEGER at2,at
      INTEGER blank,site,year,age,bioind,xyearf
      INTEGER mnth,no_days,fireres,xyear0,per,omav(douts),year_out
      INTEGER dolydo(maxnft),luse(maxyrs),fno,bb(maxnft),bbgs(maxnft)
      INTEGER iyear,oymd,oymdft,dsbb(maxnft),chill(maxnft),persum
      INTEGER stcmp,iargc,yearind(maxyrs),idum,outyears,thty_dys
      INTEGER xlatresn,xlonresn,day_mnth,yearv(maxyrs),nyears,narg
      INTEGER met_seqv(maxyrs),metyear,met_yearv(maxyrs)
      INTEGER yr0ms,yrfms,yr0m,yrfm,iyear_adj,yr0a,yrfa
      INTEGER seed1,seed2,seed3,spinl,yr0s,yr0p,yrfp,xseed1,site_dat
      INTEGER ibox,jbox,last_blank,site_out,country_id,outyears1,fti
      INTEGER outyears2,budo(maxnft),seno(maxnft),ss(maxnft),clim_type
      INTEGER check_ft_grow,no_countries,n_param_0,n_param_f,n_param
      REAL*8 xtmpv(500,12,31),xprcv(500,12,31),xhumv(500,12,31)
      REAL*8 xcldv(500,12),xswrv(500,12,31),swr(12,31),mnthswr(12)
      REAL*8 yearswr,mapv(10),maswrv(30),maprc,maswr,maprcr,maswrr
      REAL*8 maprc_init,maswr_init,leafresp,rootresp,stemresp
      REAL*8 matmpv(10),matmp_maxv(10),yeartmp_max,gi,wi,covind
      REAL*8 matmp_minv(10),yeartmp_min,map_daysv(10),yearp_days
      REAL*8 mahumv(10),masoilcv(10),yearsoilc
      REAL*8 masoilwv(10),yearsoilw,tabglitterc,tblgc,tbioleaf
      REAL*8 matmp,matmp_max,matmp_min,map_days
      REAL*8 mahum,masoilc,masoilw,eco2
      REAL*8 matmp_init,matmp_max_init,matmp_min_init,map_days_init
      REAL*8 ce_light(30,12,maxnft),ce_ci(30,12,maxnft),ce_t(30)
      REAL*8 ce_maxlight(30,12,maxnft),ce_ga(30,12,maxnft),ce_rh(30)
      REAL*8 mahum_init,masoilc_init,masoilw_init
      REAL*8 env_vcmax(maxnft),env_jmax(maxnft)
      REAL*8 env_vcmax_min(maxnft),env_jmax_min(maxnft)
      REAL*8 env_vcmax_max(maxnft),env_jmax_max(maxnft)
      REAL*8 env_sla_max(maxnft),env_sla_min(maxnft)
      REAL*8 trait_env_C4_vcmax,trait_env_C4_pepc 
      INTEGER env_vcmax_bounds(maxnft),env_jmax_bounds(maxnft)
      INTEGER env_sla_bounds(maxnft),daily_co2,par_loops,s070607
      REAL*8 can2g,tslc,trlc,active_cov

      CHARACTER st1*1000,st2*1000,st3*1000,st4*1000
      CHARACTER st5*1000,otags(douts)*1000,ofmt(200)*100,in2st*10000
      CHARACTER stinput*1000,stoutput*1000,stinit*1000,stco2*1000
      CHARACTER stmask*1000,country_name*1000,countries(100)*20
      CHARACTER sttxdp*1000,stlu*1000,ststats*1000,buff1*80
      CHARACTER stpname*1000,stpname_t*1000,stwdg*1000
      CHARACTER param_file*1000,date*8,time*10,fttags(maxnft)*1000

      LOGICAL initise,initiseo,speedc,crand,xspeedc,withcloudcover
      LOGICAL l_clim,l_lu,l_soil(20),l_stats,l_regional,l_countries
      LOGICAL out_cov,out_bio,out_bud,out_sen,l_b_and_c,check_c
      LOGICAL land_check,l_parameter,SDGVM_070607,SDGVM_140129
      LOGICAL fire(maxyrs),harvest(maxyrs),met_seq,goudriaan_old
      LOGICAL year0set
      LOGICAL debug

*----------------------------------------------------------------------*
      REAL*8 zs1(maxnft),zs2(maxnft),zs3(maxnft),zs4(maxnft)
      REAL*8 zsn(maxnft),zlsn(maxnft),zic0(8),zin0(8),ziminn(3)
      REAL*8 zdslc,zdrlc,zdsln,zdrln,znpp(maxnft),zlai(maxnft)
      REAL*8 zlat,zlon,zcov(maxage,maxnft),zbio(maxage,2,maxnft)
      REAL*8 zppm(maxage,maxnft),zhgt(maxage,maxnft)
      REAL*8 zwdt(maxage,maxnft),zleafdp(3600,maxnft)
      REAL*8 zstemdp(1000,maxnft),zrootdp(1000,maxnft)
      REAL*8 znps(maxnft),znpl(maxnft),zevp(maxnft),zdof(maxnft)
      REAL*8 zslc(maxnft),zrlc(maxnft),zsln(maxnft),zrln(maxnft)
      REAL*8 znppstore(maxnft),znppstorx(maxnft),znppstor2(maxnft)
      REAL*8 zmaxlai(maxnft),ztmem(200),zsm_trig(30,maxnft),av_hgt_all
      REAL*8 zsumadp(360,maxnft),zstemfr(maxnft),av_hgt(maxnft),max_hgt
      INTEGER zbb(maxnft),zbbgs(maxnft),zdsbb(maxnft)
      INTEGER hi,xi,gs_func
      INTEGER PHASE !PCM
      INTEGER SYR,NYR !PCM
      INTEGER aggmap_SDGVM_to_aggHyde(NS) !PCM
      REAL*8 lutab2(255,100) !PCM
      LOGICAL closed_loop_ft !PCM
      REAL*8 ftprop1(maxnft)


*----------------------------------------------------------------------*

      debug = .FALSE.

      IF (IARGC().GT.0) THEN
        CALL GETARG(1,buff1)
        narg = 1
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) ' Input file must be given as an argument.'
        STOP
      ENDIF

      no_countries = 0
      zlat = -1000.0d0
      zlon = -1000.0d0

*----------------------------------------------------------------------*
* Read input common parameters.                                        *
*----------------------------------------------------------------------*
      OPEN(98,FILE='param.dat',STATUS='OLD',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) ' File does not exist: "param.dat"'
        STOP
      ENDIF

      READ(98,*)
      READ(98,*) st1
      CALL STRIPB(st1)
      st2 = stver
      CALL STRIPBS(st2,st3)
      CALL STRIPB(st2)
      st1=st1(1:blank(st1))
      st2=st2(1:blank(st2))
      IF (stcmp(st1,st2).NE.1) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Version number mismatch between ''sdgvm0.f'' and ''p
     &aram.dat''.'
        STOP
      ENDIF

      READ(98,*)
      READ(98,*) p_sla        ! 4
      READ(98,*)
      READ(98,*) p_stemfr     ! 6
      READ(98,*)
      READ(98,*) p_rootfr     ! 8
      READ(98,*)
      READ(98,*) p_opt        !10
      READ(98,*)
      READ(98,*) p_stmin      !12
      READ(98,*)
      READ(98,*) p_laimem     !14
      READ(98,*)
      READ(98,*) p_resp       !16
      READ(98,*)
      READ(98,*) p_kscale     !18
      READ(98,*)
      READ(98,*) p_nu1,p_nu2,p_nu3,p_nu4
      READ(98,*)
      READ(98,*) p_nleaf
      READ(98,*)
      READ(98,*) p_dresp
      READ(98,*)
      READ(98,*) p_vm
      READ(98,*)
      READ(98,*) p_kgw
      READ(98,*)
      READ(98,*) p_v1,p_v2,p_v3
      READ(98,*)
      READ(98,*) p_j1,p_j2,p_j3
      READ(98,*)
      READ(98,*) p_pet
      READ(98,*)
      READ(98,*) p_bs
      READ(98,*)
      READ(98,*) p_et
      READ(98,*)
      READ(98,*) p_roff
      READ(98,*)
      READ(98,*) p_roff2
      READ(98,*)
      READ(98,*) p_fprob
      READ(98,*)
      READ(98,*) p_city_dep
      READ(98,*)
      READ(98,*) par_loops
      READ(98,*)
      READ(98,*) p_rd
      CLOSE(98)


*----------------------------------------------------------------------*
* Read 'misc_params.dat'.                                              *
*----------------------------------------------------------------------*
      OPEN(98,FILE='misc_params.dat',STATUS='OLD',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,*) ' File does not exist: "misc_params.dat"'
        WRITE(*,*) ' Using screen output options: 1 0. '
        WRITE(*,*) ' Using countries, not regions. '
        site_out = 1
        year_out = 0
        l_regional = .FALSE.
      ELSE
        READ(98,*)
        READ(98,*) site_out,year_out
        READ(98,*)
        READ(98,*) l_regional
      ENDIF
      CLOSE(98)

*----------------------------------------------------------------------*
* Read input data.                                                     *
*----------------------------------------------------------------------*
      OPEN(98,FILE=buff1,STATUS='OLD',iostat=kode)
      IF (kode.NE.0) THEN
        DO i=1,80
          st1(i:i) = buff1(i:i)
        ENDDO
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'SDGVM input file does not exist.'
        WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
        STOP
      ENDIF

      READ(98,'(A)') st1
      st2 = 'DOS'
      st3 = 'UNIX'

      IF (stcmp(st1,st2).EQ.1) THEN
        du = 1
      ELSEIF (stcmp(st1,st3).EQ.1) THEN
        du = 0
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'First line of input file must read DOS or UNIX.'
        STOP
      ENDIF

      READ(98,'(A)') st1
      ii = n_fields(st1) 
C      WRITE(*,*) '0000' 
C      WRITE(*,*) st1
C      WRITE(*,*) '0000' 
C      WRITE(*,*) stinput
      CALL STRIPBS(st1,stinput)
      st1 = stinput
      
      met_seq = .FALSE.
      IF(ii.gt.1) THEN
        st3 = 'seq'
C        WRITE(*,*) '1111' 
C        WRITE(*,*) st1
C        WRITE(*,*) 'aaaa' 
        CALL STRIPBS(st1,st2)
C        WRITE(*,*) st2
C        WRITE(*,*) 'bbbb' 
        IF (stcmp(st2,st3).EQ.1) THEN
          met_seq = .TRUE.
        ELSE
          WRITE(*,*) 'If it exists, second field of second line in 
     &input.dat must read "seq"' 
           STOP 
        ENDIF
      ENDIF

      OPEN(99,FILE=stinput(1:blank(stinput))//'/readme.dat',
     &STATUS='OLD',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Climate data file does not exist:'
        WRITE(*,'('' "'',A,''/readme.dat"'')') stinput(1:blank(stinput))
        STOP
      ENDIF
      READ(99,'(A)') st1

*----------------------------------------------------------------------*
* Determine whether climate is daily or monthly, from first line of    *
* readme file.                                                         *
*----------------------------------------------------------------------*
      st2 = 'DAILY'
      st3 = 'MONTHLY'
      st4 = 'SITED'
      st5 = 'SITEM'
      sit_grd = 0
      clim_type = 0
      IF (stcmp(st1,st2).EQ.1) THEN
        clim_type = 1
        day_mnth = 1
        thty_dys = 1
        READ(99,*)
        READ(99,*) xlatf,xlon0
        READ(99,*)
        READ(99,*) xlatres,xlonres
        READ(99,*)
        READ(99,*) xlatresn,xlonresn
        READ(99,*)
        READ(99,*) xyear0,xyearf
      ELSEIF (stcmp(st1,st3).EQ.1) THEN
        clim_type = 2
        day_mnth = 0
        thty_dys = 1
        READ(99,*)
        READ(99,*) xlatf,xlon0
        READ(99,*)
        READ(99,*) xlatres,xlonres
        READ(99,*)
        READ(99,*) xlatresn,xlonresn
        READ(99,*)
        READ(99,*) xyear0,xyearf
      ELSEIF (stcmp(st1,st4).EQ.1) THEN
        clim_type = 3
        day_mnth = 1
        thty_dys = 0
        sit_grd = 1
        READ(99,*)
        READ(99,*) xlatf,xlon0
        READ(99,*)
        READ(99,*) xyear0,xyearf
      ELSEIF (stcmp(st1,st5).EQ.1) THEN
        clim_type = 4
        day_mnth = 0
        thty_dys = 1
        sit_grd = 1
        READ(99,*)
        READ(99,*) xlatf,xlon0
        READ(99,*)
        READ(99,*) xyear0,xyearf
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*)
     &'First line of climate readme.dat file must read DAILY MONTHLY or 
     &SITE.'
        STOP
      ENDIF
      CLOSE(99)
*----------------------------------------------------------------------*

      !APW - not sure what this does
      st1 = stinput

      READ(98,'(A)') ststats
      CALL STRIPB(ststats)

      READ(98,'(A)') sttxdp  !soil data
      CALL STRIPB(sttxdp)

*----------------------------------------------------------------------*
* Read in type of landuse: 0 = defined by map; 1 = defined explicitly  *
* in the input file; 2 = natural vegetation based on average monthly   *
* temperatures.                                                        *
*----------------------------------------------------------------------*
      READ(98,'(1000a)') st1

      ii = n_fields(st1)
      IF (ii.EQ.3) THEN
        !read in ilanduse 
        CALL STRIPBN(st1,i)
        IF (i.gt.-1)  ilanduse  = i
        !read start year of extraction from LUC database 
        CALL STRIPBN(st1,i)
        IF (i.gt.-1)  SYR   = i
        !read number of years of extraction from LUC database 
        CALL STRIPBN(st1,i)
        IF (i.gt.-1)  NYR   = i
        IF (ilanduse.EQ.5. .OR. ilanduse.EQ.6) NYR = 1 !override number from input.dat if set that way accidently
      ELSE IF (ii.EQ.1) THEN
        READ(st1,*) ilanduse
        SYR = -1 !SYR and NYR not used and not defined here
        NYR = -1
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'ilanduse: either 1 or 3 arguments required'
        STOP
      ENDIF

      fire(:)    = .FALSE.
      harvest(:) = .FALSE.
      IF (ilanduse.EQ.1) THEN
* Use landuse defined in input file.
        CALL LANDUSE1(luse,fire,harvest,yr0a,yrfa,year0set,spinl)
      ENDIF
      IF ((ilanduse.LT.0).OR.(ilanduse.GT.6)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'No landuse defined'
        WRITE(*,*) '0:= Defined from a map.'
        WRITE(*,*) '1:= Defined in the input file.'
        WRITE(*,*) '2:= Natural vegetation.'
        WRITE(*,*) '3:= Defined from the maps in states2b.nc file.'
        WRITE(*,*) '4:= Defined from the maps in transitions2b.nc file.'
        WRITE(*,*)
     & '5:= Defined from the 1st-yr maps in states2b.nc file.'
        WRITE(*,*)
     & '6:= Defined from the 1st-yr maps in transitions2b.nc file.'
        STOP
      ENDIF

      READ(98,'(A)') st1 
      CALL STRIPBS(st1,stpname)

      READ(98,'(A)') st1 
      CALL STRIPBS(st1,stpname_t)

      READ(98,'(A)') st1 
      CALL STRIPBS(st1,stwdg)
      !stwdg=stwdg(1:blank(stwdg))

      READ(98,'(A)') st1
      CALL STRIPBS(st1,stlu)

      READ(98,'(A)') stco2
      CALL STRIPB(stco2)

      READ(98,*) co2const

      READ(98,'(A)') stmask
      CALL STRIPB(stmask)

*----------------------------------------------------------------------*
* read input switches                                                  *
*----------------------------------------------------------------------*

      !line 9 in <input.dat> determines SDVGM version
      !0 is current version 1 is old 070607 version
      READ(98,'(1000a)') st1 !      
      ii = n_fields(st1) !
      IF (ii.EQ.1) THEN !
        CALL STRIPBN(st1,i) !
        if((i.lt.0).or.(i.gt.2)) then
          WRITE(*,'('' PROGRAM TERMINATED'')') !
          WRITE(*,*) 'Line 9 must be 2/1/0' !
          WRITE(*,'('' "'',A,''"'')') st1(1:30) !
          STOP !        
        else
          IF (i.eq.1)then
            SDGVM_070607 = .TRUE. 
            SDGVM_140129 = .FALSE.
          elseif(i.eq.0) then
            SDGVM_070607 = .FALSE.
            SDGVM_140129 = .TRUE.
          else
            SDGVM_070607 = .FALSE.
            SDGVM_140129 = .FALSE.
          endif
        endif        
      ELSE !
        WRITE(*,'('' PROGRAM TERMINATED'')') !
        WRITE(*,*) 'Line 9 must contain 1 field' !
        WRITE(*,'('' "'',A,''"'')') st1(1:30) !
        STOP !        
      ENDIF !

      !read next line of input switches. ln 10 in input.dat
      READ(98,'(1000a)') st1 
      ii = n_fields(st1) 
      IF (ii.EQ.6) THEN 
        CALL STRIPBN(st1,i)
        !read in daily co2 
        IF (i.gt.-1)  daily_co2  = i 
        CALL STRIPBN(st1,i) 
        !read PAR data 
        IF (i.gt.-1)  read_par   = i 
        CALL STRIPBN(st1,i)
        !use sub-daily PAR scaling
        IF (i.gt.-1)  subd_par   = i 
        CALL STRIPBN(st1,i) 
        !canopy clumping index: 0-do not use; 1-from pft values; 2-from a map
        IF (i.gt.-1)  read_clump = i 
        CALL STRIPBN(st1,i) 
        !calculate solar zenith angle and use in canopy light interception 
        IF (i.gt.-1)  calc_zen   = i 
        CALL STRIPBN(st1,i) 
        !no soil water limitation
        IF (i.gt.-1)  no_slw_lim = i 
      ELSE !
        WRITE(*,'('' PROGRAM TERMINATED'')') !
        WRITE(*,*) 'Line 10 must contain 6 fields' !
        WRITE(*,'('' "'',A,''"'')') st1(1:30) !
        STOP !       
      ENDIF !

      if(subd_par.eq.1) then
        print*, ''
        print*, 'SDGVM running with ',par_loops,
     &'sub-daily photosynthesis time-points'
        print*, ''
      else
        print*, ''
        print*, 'SDGVM running with ','1',
     &'sub-daily photosynthesis time-points'
        print*, ''
      endif

      !read next line of input switches. ln 11 in input.dat
      READ(98,'(1000a)') st1 
      ii = n_fields(st1) 
      IF (ii.EQ.5) THEN 
        CALL STRIPBN(st1,i)
        !select canopy nitrogen calculation method
        cstype = 0
        IF (i.ge.10) THEN
          cstype = int(real(i)/10.0)
          i      = i - 10*cstype
        ENDIF
        IF (i.gt.-1)  ncalc_type = i 
        CALL STRIPBN(st1,i) 
        !select vcmax parameterisation
        ttype = 0
        IF (i.ge.10) THEN
          ttype = int(real(i)/10.0)
          i     = i - 10*ttype
        ENDIF
        IF (i.gt.-1)  vcmax_type = i 
        CALL STRIPBN(st1,i) 
        !use soil P from a map to calculate Vcmax (Vcmax switch must also be 1 for P to be used in Vcmax calc) 
        !- even if this variable is 0, it is expected to be read from the soils database 
        IF (i.gt.-1)  soilp_map = i 
        CALL STRIPBN(st1,i) 
        !switch to run with 070607 routines that have been changed in the new version but don't have individual switches (rd, swlim, et multiplier etc)
        IF (i.gt.-1)  s070607   = i 
        CALL STRIPBN(st1,i) 
        !switch stomatal conductance function
        IF (i.gt.-1)  gs_func = i 
        !CALL STRIPBN(st1,i) 
        !switch electron transport function
        !IF (i.gt.-1)  hw_j = i 
      ELSE !
        WRITE(*,'('' PROGRAM TERMINATED'')') !
        WRITE(*,*) 'Line 11 must contain 5 fields' !
        WRITE(*,'('' "'',A,''"'')') st1(1:30) !
        STOP !
      ENDIF !

      
      !the below switches are hard coded as they are unlikely to need changing
      ! - they are only changed if the old (070607) version of the model is used  

      !use soil C:N ratio from a map 
      !- even if this variable is 0, it is expected to be read from the soils database
      ! - 0 is the default for this switch, use this switch to implement a routine that uses soil C:N read from the soils database 
      soilcn_map = 0  

      !grass lai can only take 50% of stored C and leaf growth subject to growth respiration
      phen_cor   = 1 

      !switch electron transport function, 0 - Harley 1992, 1 - Farquhar
      !& Wong 1984
      hw_j       = 0

      !this is commented out for ease of adding new swithces to the input.dat for model development
      !read next line of input switches. ln 12 in input.dat
!      READ(98,'(1000a)') st1 
!      ii = n_fields(st1) 
!      IF (ii.EQ.4) THEN 
!      CALL STRIPBN(st1,i) 
!      !master switch for these switches
!      IF (i.gt.-1)  mswitch = i 
!      CALL STRIPBN(st1,i) 
!      !switch 1 
!      IF (i.gt.-1)  switch1 = i 
!      !switch 2
!      IF (i.gt.-1)  switch2 = i 
!      !switch 3
!      IF (i.gt.-1)  switch3 = i 
!      ELSE !
!        WRITE(*,'('' PROGRAM TERMINATED'')') !
!        WRITE(*,*) 'Line 12 must contain 4 fields' !
!        WRITE(*,'('' "'',A,''"'')') st1(1:30) !
!        STOP !   
!      ENDIF !

      if(gs_func.eq.3) goudriaan_old = .TRUE.

      !set default configurations for standard versions
      IF(SDGVM_070607) THEN
        daily_co2  = 0 
        read_par   = 0
        subd_par   = 0 
        read_clump = 0
        calc_zen   = 0
        no_slw_lim = 0

        cstype     = 0  
        ncalc_type = 0
        ttype      = 0
        vcmax_type = 0
        soilp_map  = 0
        s070607    = 1
        gs_func    = 0

        soilcn_map = 0
        phen_cor   = 0 
        hw_j       = 0
      ENDIF

      IF(SDGVM_140129) THEN
        daily_co2  = 0 
        read_par   = 1
        subd_par   = 1
        read_clump = 0
        calc_zen   = 1 
        no_slw_lim = 0

        cstype     = 0
        ncalc_type = 1
        ttype      = 0
        vcmax_type = 8 ! Kattge 2009 oxisol for EvBl, Jmax~Vcmax from Walker 2014
        soilp_map  = 0 
        s070607    = 0
        gs_func    = 0

        soilcn_map = 0
        phen_cor   = 1
        hw_j       = 0
      ENDIF

      if(goudriaan_old) hw_j = 3

      !parameters for 070607 version
      if(s070607.eq.1) then
        p_et  = 0.7
        p_pet = 0.7
      endif

      READ(98,*)

      !read input and output directories from input.dat
      READ(98,'(A)') stinit
      CALL STRIPB(stinit)
      st2 = 'ARGUMENT'
      IF (stcmp(stinit,st2).EQ.1) THEN
        narg = narg + 1
        CALL GETARG(narg,buff1)
        DO i=1,80
          stinit(i:i) = buff1(i:i)
        ENDDO
      ENDIF

      READ(98,'(A)') stoutput
      CALL STRIPB(stoutput)
      st2 = 'ARGUMENT'
      IF (stcmp(stoutput,st2).EQ.1) THEN
        narg = narg + 1
        CALL GETARG(narg,buff1)
        DO i=1,80
          stoutput(i:i) = buff1(i:i)
        ENDDO
      ENDIF
      READ(98,*)

*----------------------------------------------------------------------*
* Open diagnostics file.                                               *
*----------------------------------------------------------------------*
      st1 = stoutput
      OPEN(11,FILE=st1(1:blank(st1))//'/diag.dat',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'SDGVM output directory does not exist.'
        WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
        STOP
      ENDIF
      WRITE(11,'(''ERRORS'')')

*----------------------------------------------------------------------*
* Read run type switches
*----------------------------------------------------------------------*
      READ(98,'(1000a)') st1
      ii = n_fields(st1)
      IF ((ii.EQ.3).OR.(ii.EQ.4)) THEN
        CALL STRIPBN(st1,i)
        initise = .true.
        IF (i.EQ.0)  initise = .false.
        CALL STRIPBN(st1,j)
        initiseo = .true.
        IF (j.EQ.0)  initiseo = .false.
        CALL STRIPBN(st1,k)
        speedc = .true.
        IF (k.EQ.0)  speedc = .false.
        IF (ii.EQ.4) THEN
          CALL STRIPBN(st1,xseed1)
        ELSE
          xseed1 = 1
        ENDIF
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Line 13 must contain 3 or 4 arguments'
        WRITE(*,'('' "'',A,''"'')') st1(1:30)
        STOP
      ENDIF
      xspeedc = speedc

      READ(98,'(1000a)') st1
      CALL STRIPB(st1)
      i = n_fields(st1)
      IF (i.EQ.7) THEN
        CALL STRIPBN(st1,spinl)
        CALL STRIPBN(st1,yr0s)
        CALL STRIPBN(st1,cycle)
        CALL STRIPBN(st1,j)
        crand = .TRUE.
        IF (j.EQ.0)  crand = .FALSE.
        CALL STRIPBN(st1,yr0p)
        CALL STRIPBN(st1,yrfp)
        CALL STRIPBN(st1,outyears)
      ELSEIF (i.EQ.6) THEN
        CALL STRIPBN(st1,spinl)
        CALL STRIPBN(st1,yr0s)
        CALL STRIPBN(st1,cycle)
        CALL STRIPBN(st1,j)
        crand = .TRUE.
        IF (j.EQ.0)  crand = .FALSE.
        CALL STRIPBN(st1,yr0p)
        CALL STRIPBN(st1,yrfp)
        outyears = yrfp - yr0p + 1
      ELSEIF (i.EQ.5) THEN
        CALL STRIPBN(st1,spinl)
        CALL STRIPBN(st1,yr0s)
        CALL STRIPBN(st1,cycle)
        CALL STRIPBN(st1,j)
        crand = .TRUE.
        IF (j.EQ.0)  crand = .FALSE.
        yr0p = yr0s+1
        yrfp = yr0s
        CALL STRIPBN(st1,outyears)
      ELSEIF (i.EQ.4) THEN
        CALL STRIPBN(st1,spinl)
        CALL STRIPBN(st1,yr0s)
        CALL STRIPBN(st1,cycle)
        CALL STRIPBN(st1,j)
        crand = .TRUE.
        IF (j.EQ.0)  crand = .FALSE.
        yr0p = yr0s+1
        yrfp = yr0s
        outyears = cycle + 1
      ELSEIF (i.EQ.3) THEN
        spinl = 0
        cycle = maxyrs
        j = 0
        crand = .FALSE.
        CALL STRIPBN(st1,yr0p)
        CALL STRIPBN(st1,yrfp)
        CALL STRIPBN(st1,outyears)
        yr0s = yr0p
      ELSEIF (i.EQ.2) THEN
        spinl = 0
        cycle = maxyrs
        j = 0
        crand = .FALSE.
        CALL STRIPBN(st1,yr0p)
        CALL STRIPBN(st1,yrfp)
        yr0s = yr0p
        outyears = yrfp - yr0p + 1
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Line 14 must contain 2-7 fields'
        STOP
      ENDIF
      nyears = spinl+yrfp-yr0p+1
      IF (nyears.GT.maxyrs) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,'('' Trying to simulate '',i4,'' years, maximum allowabl
     &e is '',i4,''.'')') nyears,maxyrs
        WRITE(*,*) 'Either reduce the length of the simulation, or incr
     &ease "maxyrs"'
        WRITE(*,*) 'this is set in array_param.inc, you must re-compile 
     &after altering'
      WRITE(*,*) 'this file.'
        STOP
      ENDIF
      IF ((j.NE.0).AND.(j.NE.1)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,'('' Fourth field of line 14 in the input file must be e
     &ither 0 or 1.'')')
        WRITE(*,'('' Currently set to '',i5,''.'')') j
        STOP
      ENDIF
      outyears = min(outyears,nyears)


*----------------------------------------------------------------------*
* Read in a met sequence if specified   -------------------------------*
*----------------------------------------------------------------------*
      yr0ms = yr0p
      yrfms = yrfp 
      IF(met_seq) THEN
        OPEN(99,FILE='met_seq.dat',iostat=kode)
        IF (kode.NE.0) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'non-sequential met sequence selected but file 
     &does not exist.'
          WRITE(*,'('' "'',A,''"'')') 'met_seq.dat'
          STOP
        ENDIF
        DO i=1,(yrfp-yr0p+1)
          READ(99,*) met_seqv(i)
          IF(i.eq.1) yr0ms = met_seqv(i)
          if(i.gt.1) yr0ms = min(yr0ms,met_seqv(i))
          IF(i.eq.1) yrfms = met_seqv(i)
          if(i.gt.1) yrfms = max(yrfms,met_seqv(i))
        ENDDO
        CLOSE(99)
      ENDIF

*----------------------------------------------------------------------*
* Set 'yr0' and 'yrf' which are the years of actual climate required   *
* for the run. And check that the climate exists in the climate        *
* database                                                             *
*----------------------------------------------------------------------*
      yr0 = min(yr0s,yr0p)
      yrf = max(yr0s+min(spinl,cycle)-1,yrfp)

      ! allow a randomised (or otherwise) sequence of years to be used for 
      ! climate data during the run proper (transient run)
      yr0m = yr0
      yrfm = yrf
      IF(met_seq) yr0m = min(yr0s,yr0ms)
      IF(met_seq) yrfm = max(yr0s+min(spinl,cycle)-1,yrfms)
      IF ((yr0m.LT.xyear0).OR.(yrfm.GT.xyearf)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,'('' Trying to use '',i4,''-'',i4,'' climate.'')') 
     &yr0m,yrfm 
        WRITE(*,'('' Climate database runs from '',i4,''-'',i4,''.'')') 
     &xyear0,xyearf
        STOP
      ENDIF

      ! allow more flexibility in how climate, land-use, and CO2 data are combined.
      ! Specifically, when co2const < 0 allow climate to cycle through spin, 
      ! but use dynamic land-use and co2   
      yr0a = yr0
      yrfa = yrf
      year0set = .TRUE.
C PCM changed the following line for TRENDY S4-S6
      if( (co2const.lt.0.0) .and. (spinl.gt.0) ) then !For TRENDY S1-S3 
C PCM      if( spinl.gt.0 ) then  !For TRENDY S4-S6
        if (spinl.lt.nyears) then
          yr0a = yr0 - spinl
        else
          year0set = .FALSE.
C PCM changed the following two lines for TRENDY S4-S6
          yr0a = 1     !For TRENDY S1-S3
          yrfa = spinl !For TRENDY S1-S3 
C          yr0a = 1700 !For TRENDY S4-S6 
C          yrfa = 2018 !For TRENDY S4-S6 
        endif
      endif

*----------------------------------------------------------------------*
* Set up year climate sequence.                                        *
*----------------------------------------------------------------------*
      DO i=1,maxyrs
        yearind(i) = i
      ENDDO
      idum = 1
      PHASE = 0 !PCM 0-year offset for TRENDY 2018 S1-S3
C      PHASE = 1 !PCM 1-year offset for TRENDY 2018 S4-S6
      DO i=1,nyears
        IF (i.LE.spinl) THEN
          IF (crand) THEN
C PCM            IF ((mod(i-1,cycle).EQ.0).AND.(i.GT.2))
C PCM     & CALL RANDOMV(yearind,1,cycle,idum)
            IF ((mod(i-1,cycle).EQ.0)) 
     & CALL RANDOMV(yearind,1,cycle,idum)
            yearv(i) = yearind(mod(i-1,cycle)+1) + yr0s - 1
          ELSE
C PCM  changed the following line for TRENDY 2018 S4-S6
            yearv(i) = mod(i-1,cycle) + yr0s !For TRENDY S1-S3
C PCM       yearv(i) = mod(i-1+PHASE,cycle) + yr0s !For TRENDY S4-S6
          ENDIF
          met_yearv(i) = yearv(i)
        ELSE
          yearv(i) = i - spinl + yr0p - 1
          IF(met_seq) THEN
            met_yearv(i) = met_seqv(i-spinl)
          ELSE
            met_yearv(i) = yearv(i)
          ENDIF
        ENDIF
      ENDDO

*----------------------------------------------------------------------*
* Open monthly/daily PIXEL output files if required.                   *
*----------------------------------------------------------------------*
      CALL OUTPUT_OPTIONS(otags,omav,ofmt)

      READ(98,'(A)') st1
      CALL STRIPB(st1)
      st2 = 'PIXEL'
      IF (stcmp(st1,st2).EQ.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'The first field of line 15 of the input file must
     & read ''PIXEL''.'
        WRITE(*,*) st1(1:30)
        WRITE(*,*) 'Output variable options:'
        WRITE(*,'(1x,20a4)') (otags(i),i=1,15)
        WRITE(*,'(1x,20a4)') (otags(i),i=16,douts)
        STOP
      ENDIF

      CALL SET_PIXEL_OUT(st1,st2,outyears1,otagsn,otags,
     &oymd)
      outyears1 = min(outyears1,nyears)

      st1 = stoutput
      iofn = 400
      IF (outyears1.GT.0) THEN
      IF ((oymd.EQ.1).OR.(oymd.EQ.3)) THEN
        DO i=1,douts
          IF (otagsn(i).EQ.1) THEN
            iofn = iofn + 1
            OPEN(iofn,file=st1(1:blank(st1))//'/monthly_'//otags(i)(1:3)
     &//'.dat')
          ENDIF
        ENDDO
      ENDIF
      IF ((oymd.EQ.2).OR.(oymd.EQ.3)) THEN
        DO i=1,douts
          IF (otagsn(i).EQ.1) THEN
            iofn = iofn + 1
            OPEN(iofn,file=st1(1:blank(st1))//'/daily_'//otags(i)(1:3)
     &//'.dat')
          ENDIF
        ENDDO
      ENDIF
      ENDIF
      iofnft = iofn

*----------------------------------------------------------------------*
* Determine whether daily or monthly subpixel outputs are required.    *
*----------------------------------------------------------------------*

      READ(98,'(A)') st1
      CALL STRIPB(st1)
      st2 = 'SUBPIXEL'
      IF (stcmp(st1,st2).EQ.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'The first field of line 16 of the input file must
     & read ''SUBPIXEL''.'
        WRITE(*,'('' "'',A,''"'')') st1(1:30)
        WRITE(*,*) 'Output variable options:'
        WRITE(*,'(1x,20a4)') (otags(i),i=1,15)
        WRITE(*,'(1x,20a4)') (otags(i),i=16,douts),
     &'cov ','bio ','bud ','sen '
        STOP
      ENDIF

      CALL SET_SUBPIXEL_OUT(st1,st2,outyears2,otagsnft,
     &otags,oymdft,out_cov,out_bio,out_bud,out_sen)
      outyears2 = min(outyears2,nyears)

*----------------------------------------------------------------------*
* Read in snapshot years.                                              *
*----------------------------------------------------------------------*
      READ(98,'(A)') st1
      CALL STRIPB(st1)
      st2 = 'SNAPSHOTS'
      IF (stcmp(st1,st2).EQ.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'First field of line 17 must read SNAPSHOTS'
        STOP
      ENDIF
      CALL STRIPBS(st1,st2)
      CALL STRIPB(st1)
      CALL ST2ARR(st1,snpshts,100,snp_no)
      READ (98,*)

*----------------------------------------------------------------------*
* Read in compulsory functional types.                                 *
*----------------------------------------------------------------------*
      fttags(1) = 'BARE'
      fttags(2) = 'CITY'
      READ(98,'(A)') st1
      CALL STRIPB(st1)
      st2 = 'BARE'
      IF ((stcmp(st2,st1).EQ.0).OR.(n_fields(st1).NE.2)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Line 19 must read "BARE" forllowed by a number (0-1)
     &.'
        STOP
      ENDIF
      CALL STRIPBS(st1,st2)
      READ(st1,*) ftmix(1)
      READ(98,'(A)') st1
      CALL STRIPB(st1)
      st2 = 'CITY'
      IF ((stcmp(st2,st1).EQ.0).OR.(n_fields(st1).NE.2)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Line 19 must read "BARE" forllowed by a number (0-1)
     &.'
        STOP
      ENDIF
      CALL STRIPBS(st1,st2)
      READ(st1,*) ftmix(2)
      READ (98,*)

* Initialise redundant parameterisation 
      ft = 1
      ftc3(ft) = 0
      ftphen(ft) = 0
      ftagh(ft) = 0.0d0
      ftdth(ft) = 0
      ftmor(ft) = 0
      ftwd(ft) = 0.0d0
      ftxyl(ft) = 0.0d0
      ftpd(ft) = 0.0d0
      ftsla(ft) = 0.0d0
      ftlls(ft) = 0
      ftsls(ft) = 0
      ftrls(ft) = 0
      ftlmor(ft) = 0.0d0
      ftrat(ft) = 0.0d0
      ftbbm(ft) = 0
      chill(ft) = 0
      dschill(ft) = 0
      ftbb0(ft) = 0.0d0
      ftbbmax(ft) = 0.0d0
      ftbblim(ft) = 0.0d0
      ftssm(ft) = 0
      ftsss(ft) = 0
      ftsslim(ft) = 0.0d0
      ftstmx(ft) = 0.0d0
      ftgr0(ft) = 0.0d0
      ftgrf(ft) = 0.0d0
      ftppm0(ft) = 0.0d0
      ftcan_clump(ft) = 0.0d0
      ftvna(ft)  = 0.0d0
      ftvnb(ft)  = 0.0d0
      ftjva(ft)  = 0.0d0
      ftjvb(ft)  = 0.0d0
      ftg0(ft)   = 0.0d0
      ftg1(ft)   = 0.0d0
      ftToptV(ft) = 0.0d0
      ftHaV(ft)   = 0.0d0
      ftHdV(ft)   = 0.0d0
      ftToptJ(ft) = 0.0d0
      ftHaJ(ft)   = 0.0d0
      ftHdJ(ft)   = 0.0d0
      vcmax(:,ft) = 0.0d0
      jmax(:,ft)  = 0.0d0
      pnlc(:,ft)  = 0.0d0
      enzs(:,ft)  = 0.0d0

      ft = 2
      ftc3(ft) = 0
      ftphen(ft) = 0
      ftagh(ft) = 0.0d0
      ftdth(ft) = 0
      ftmor(ft) = 0
      ftwd(ft) = 0.0d0
      ftxyl(ft) = 0.0d0
      ftpd(ft) = 0.0d0
      ftsla(ft) = 0.0d0
      ftlls(ft) = 0
      ftsls(ft) = 0
      ftrls(ft) = 0
      ftlmor(ft) = 0.0d0
      ftrat(ft) = 0.0d0
      ftbbm(ft) = 0
      chill(ft) = 0
      dschill(ft) = 0
      ftbb0(ft) = 0.0d0
      ftbbmax(ft) = 0.0d0
      ftbblim(ft) = 0.0d0
      ftssm(ft) = 0
      ftsss(ft) = 0
      ftsslim(ft) = 0.0d0
      ftstmx(ft) = 0.0d0
      ftgr0(ft) = 0.0d0
      ftgrf(ft) = 0.0d0
      ftppm0(ft) = 0.0d0
      ftcan_clump(ft) = 0.0d0
      ftvna(ft)   = 0.0d0
      ftvnb(ft)   = 0.0d0
      ftjva(ft)   = 0.0d0
      ftjvb(ft)   = 0.0d0
      ftg0(ft)    = 0.0d0
      ftg1(ft)    = 0.0d0
      ftToptV(ft) = 0.0d0
      ftHaV(ft)   = 0.0d0
      ftHdV(ft)   = 0.0d0
      ftToptJ(ft) = 0.0d0
      ftHaJ(ft)   = 0.0d0
      ftHdJ(ft)   = 0.0d0
      vcmax(:,ft) = 0.0d0
      jmax(:,ft)  = 0.0d0
      pnlc(:,ft)  = 0.0d0
      enzs(:,ft)  = 0.0d0

*----------------------------------------------------------------------*
* Read in functional type parameterisation.                            *
*----------------------------------------------------------------------*
      pft_nflds = 34
      IF(ttype.ge.1) pft_nflds = 40

      READ(98,'(A)') st1
      CALL STRIPB(st1)

      IF (n_fields(st1).EQ.1) THEN
        st2 = 'ARGUMENT'

        IF (stcmp(st2,st1).EQ.1) THEN
*----------------------------------------------------------------------*
* From a file.                                                         *
*----------------------------------------------------------------------*
          narg = narg + 1
          CALL GETARG(narg,buff1)
          DO i=1,80
            st1(i:i) = buff1(i:i)
          ENDDO
        ENDIF

        OPEN(99,FILE=st1,STATUS='OLD',iostat=kode)
        IF (kode.NE.0) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'Functional type parameterisation file does not ex
     &ist.'
          WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
          STOP
        ENDIF

        ft = 2

40      CONTINUE
        READ(99,'(1000a)',end=30) st1
        ft = ft + 1

        IF (ft.GT.maxnft) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,'('' Maximum number of fts is'',i4,'', currently m
     &ore than this are defined.'')') maxnft
          WRITE(*,*) 'Either increase "maxnft" defined in "array_dim
     &s.dat" or decrease'
          WRITE(*,*) 'the number of ft parameterisations: this requi
     &res re-compilation.'
          STOP
        ENDIF

        IF (n_fields(st1).NE.pft_nflds) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'The ft parameterisation must contain',
     &pft_nflds,'fields.'
          WRITE(*,'(1x,A,'' has '',i3)') st1(1:blank(st1)),
     &n_fields(st1)
          STOP
        ENDIF

        IF(ttype.eq.0) THEN
          READ(st1,*) fttags(ft),ftmix(ft),ftc3(ft),ftphen(ft),
     &ftagh(ft),ftdth(ft),ftmor(ft),ftwd(ft),ftxyl(ft),ftpd(ft),
     &ftsla(ft),ftlls(ft),ftsls(ft),ftrls(ft),ftlmor(ft),ftrat(ft),
     &ftbbm(ft),ftbb0(ft),ftbbmax(ft),ftbblim(ft),ftssm(ft),ftsss(ft),
     &ftsslim(ft),ftstmx(ft),ftgr0(ft),ftgrf(ft),ftppm0(ft),
     &ftcan_clump(ft),ftvna(ft),ftvnb(ft),ftjva(ft),ftjvb(ft),ftg0(ft),
     &ftg1(ft)
      ftToptV(ft) = 0.0d0
      ftHaV(ft)   = 0.0d0
      ftHdV(ft)   = 0.0d0
      ftToptJ(ft) = 0.0d0
      ftHaJ(ft)   = 0.0d0
      ftHdJ(ft)   = 0.0d0
        ELSEIF(ttype.ge.1) THEN
          READ(st1,*) fttags(ft),ftmix(ft),ftc3(ft),ftphen(ft),
     &ftagh(ft),ftdth(ft),ftmor(ft),ftwd(ft),ftxyl(ft),ftpd(ft),
     &ftsla(ft),ftlls(ft),ftsls(ft),ftrls(ft),ftlmor(ft),ftrat(ft),
     &ftbbm(ft),ftbb0(ft),ftbbmax(ft),ftbblim(ft),ftssm(ft),ftsss(ft),
     &ftsslim(ft),ftstmx(ft),ftgr0(ft),ftgrf(ft),ftppm0(ft),
     &ftcan_clump(ft),ftvna(ft),ftvnb(ft),ftjva(ft),ftjvb(ft),ftg0(ft),
     &ftg1(ft),ftToptV(ft),ftHaV(ft),ftHdV(ft),ftToptJ(ft),ftHaJ(ft),
     &ftHdJ(ft)
        ENDIF

        IF (ftmor(ft).GT.maxyrs) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,'('' Mortality of '',A,'' is '',i4,'', maximum allow
     &able is '',i4,''.'')')
          STOP
        ENDIF

        IF (ftsla(ft).LT.0.0d0) THEN
          ftsla(ft) = 10.0d0**(2.35d0 -
     &0.39d0*log10(real(ftlls(ft))/30.0d0))*2.0d0/10000.0d0
*      ftsla(ft) = 10.0d0**(2.43d0-0.46d0*log10(real(ftlls(ft))/30.0d0))
*     &*2.0/10000.0
        ENDIF

        ftsla(ft) = ftsla(ft)/p_sla

        IF (ftlls(ft).LT.0.0d0) THEN
          ftlls(ft) = int(10.0d0**((2.35d0 -
     &log10(ftsla(ft)*10000.0d0/2.0d0))/0.39d0)*30.0d0+0.5d0)
        ENDIF

        IF (ftmor(ft).GT.maxage) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,'('' Maximum age of "'',A,''" is '',i4,'', maximum a
     &llowable age is '',i4,''.'')') st1(1:blank(st1)),ftmor(ft),maxage
          WRITE(*,*) 'Either increase "maxage" defined in "array_dims"
     &'
          WRITE(*,*) 'or decrease mortality in the parameterisation.'
          STOP
        ENDIF

        GOTO 40
30      CONTINUE
        nft = ft
        READ (98,*)

      ELSE
*----------------------------------------------------------------------*
* From the input file.                                                 *
*----------------------------------------------------------------------*
        ft = 2
98      CONTINUE

        IF (ichar(st1(1:1)).NE.32) THEN
            ft = ft + 1

            IF (ft.GT.maxnft) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,'('' Maximum number of fts is'',i4,'', currently m
     &ore than this are defined.'')') maxnft
              WRITE(*,*) 'Either increase "maxnft" defined in "array_dim
     &s.dat" or decrease'
              WRITE(*,*) 'the number of ft parameterisations: this requi
     &res re-compilation.'
              STOP
            ENDIF

            IF (n_fields(st1).NE.pft_nflds) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,*) 'The ft parameterisation must contain',
     &pft_nflds,'fields.'
              WRITE(*,'(1x,A,'' has '',i3)') st1(1:blank(st1)),n_fields(
     &st1)
              STOP
            ENDIF

            IF(ttype.eq.0) THEN
              READ(st1,*) fttags(ft),ftmix(ft),ftc3(ft),ftphen(ft),
     &ftagh(ft),ftdth(ft),ftmor(ft),ftwd(ft),ftxyl(ft),ftpd(ft),
     &ftsla(ft),ftlls(ft),ftsls(ft),ftrls(ft),ftlmor(ft),ftrat(ft),
     &ftbbm(ft),ftbb0(ft),ftbbmax(ft),ftbblim(ft),ftssm(ft),ftsss(ft),
     &ftsslim(ft),ftstmx(ft),ftgr0(ft),ftgrf(ft),ftppm0(ft),
     &ftcan_clump(ft),ftvna(ft),ftvnb(ft),ftjva(ft),ftjvb(ft),ftg0(ft),
     &ftg1(ft)
      ftToptV(ft) = 0.0d0
      ftHaV(ft)   = 0.0d0
      ftHdV(ft)   = 0.0d0
      ftToptJ(ft) = 0.0d0
      ftHaJ(ft)   = 0.0d0
      ftHdJ(ft)   = 0.0d0
            ELSEIF(ttype.ge.1) THEN
              READ(st1,*) fttags(ft),ftmix(ft),ftc3(ft),ftphen(ft),
     &ftagh(ft),ftdth(ft),ftmor(ft),ftwd(ft),ftxyl(ft),ftpd(ft),
     &ftsla(ft),ftlls(ft),ftsls(ft),ftrls(ft),ftlmor(ft),ftrat(ft),
     &ftbbm(ft),ftbb0(ft),ftbbmax(ft),ftbblim(ft),ftssm(ft),ftsss(ft),
     &ftsslim(ft),ftstmx(ft),ftgr0(ft),ftgrf(ft),ftppm0(ft),
     &ftcan_clump(ft),ftvna(ft),ftvnb(ft),ftjva(ft),ftjvb(ft),ftg0(ft),
     &ftg1(ft),ftToptV(ft),ftHaV(ft),ftHdV(ft),ftToptJ(ft),ftHaJ(ft),
     &ftHdJ(ft)
            ENDIF

            IF (ftsla(ft).LT.0.0d0) THEN
              ftsla(ft) = 10.0d0**(2.35d0 -
     &0.39d0*log10(real(ftlls(ft))/30.0d0))*2.0d0/10000.0d0
*      ftsla(ft) = 10.0d0**(2.43d0-0.46d0*log10(real(ftlls(ft))/30.0d0))
*     &*2.0/10000.0
            ENDIF

            ftsla(ft) = ftsla(ft)/p_sla

            IF (ftlls(ft).LT.0.0d0) THEN
              ftlls(ft) = int(10.0d0**((2.35d0 -
     &log10(ftsla(ft)*10000.0d0/2.0d0))/0.39d0)*30.0d0+0.5d0)
            ENDIF

            IF (ftmor(ft).GT.maxage) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,'('' Maximum age of "'',A,''" is '',i4,'', maximum
     & allowable age is '',i4,''.'')') st1(1:blank(st1)),ftmor(ft),maxag
     &e
              WRITE(*,*) 'Either increase "maxage" defined in "array_dim
     &s.dat"'
              WRITE(*,*) 'or decrease mortality in the parameterisation
     &.'
              STOP
            ENDIF

            READ(98,'(1000a)') st1
            CALL STRIPB(st1)
          GOTO 98        
        ENDIF
        nft = ft
      ENDIF

      ! preverve inpyut sla value to overwrite trait environment
      ! relationships below when needed 
      input_ftsla(:) = ftsla(:)

*----------------------------------------------------------------------*
* Create leaf mortality scales values.                                 *
*----------------------------------------------------------------------*
      lmor_sc(1,1)=0.0d0
      lmor_sc(2,1)=0.0d0
      DO ft=3,nft
        DO i=1,ftlls(ft)
          lmor_sc(i,ft)=(real(ftlls(ft)-i)/real(ftlls(ft)))**ftlmor(ft)
        ENDDO
      ENDDO

*----------------------------------------------------------------------*
* Open site_info file and check output directory exists.               *
*----------------------------------------------------------------------*
      st1 = stoutput
      OPEN(12,FILE=st1(1:blank(st1))//'/site_info.dat',iostat=kode)
      IF (kode.NE.0) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Output directory does not exist:'
        WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
        STOP
      ENDIF
      WRITE(12,'(''COUNTRY           ID     lat      lon     co20  co2f 
     &   sand   silt   bulk   orgc     wp     fc    swc    dep    sC:N 
     &      ''
     &,50A13)') (fttags(i),i=1,nft)

*----------------------------------------------------------------------*
* Open monthly/daily ft output files.                                  *
*----------------------------------------------------------------------*
      IF (outyears2.GT.0) THEN
      iofn = iofnft
      st1 = stoutput
      IF ((oymdft.EQ.1).OR.(oymdft.EQ.3)) THEN
        DO i=1,douts
          IF (otagsnft(i).EQ.1) THEN
            DO ft=1,nft
              iofn = iofn + 1
              OPEN(iofn,file=st1(1:blank(st1))//'/monthly_'//otags(i)(1:
     &3)//'_'//fttags(ft)(1:blank(fttags(ft)))//'.dat')
           ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF ((oymdft.EQ.2).OR.(oymdft.EQ.3)) THEN
        DO i=1,douts
          IF (otagsnft(i).EQ.1) THEN
            DO ft=1,nft
              iofn = iofn + 1
              OPEN(iofn,file=st1(1:blank(st1))//'/daily_'//otags(i)(1:3)
     &//'_'//fttags(ft)(1:blank(fttags(ft)))//'.dat')
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      ENDIF

      DO i=1,nft


        chill(i) = 0
        dschill(i) = 0
      ENDDO

*----------------------------------------------------------------------*
* Read in landuse index mapping.                                       *
*----------------------------------------------------------------------*
      DO j=1,255
        DO k=1,nft
          lutab(j,k) = 0.0d0
        ENDDO
      ENDDO 

      ii = 0
96    CONTINUE
        READ(98,'(1000a)') st1
        ii = ii + 1
        CALL STRIPB(st1)
c     read table of conversion from class to ft's proportion
        IF (ichar(st1(1:1)).NE.32) THEN
          i = n_fields(st1)
          CALL STRIPBN(st1,j)
          persum = 0
          DO k=1,(i-1)/2
            CALL STRIPBN(st1,per)
            persum = persum + per
            CALL STRIPBS(st1,st2)
            l = ntags(fttags,st2)
            IF (l.EQ.-9999) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,*) 'Error in a tag name in the land use mapping.'
              WRITE(*,'('' Line number'',i3,'', within the mapping list.
     &'')') ii
              WRITE(*,'('' Category number'',i3,'', tag number'',i2,''.'
     &')') j,k
              WRITE(*,'('' "'',A,''"'')') st2(1:blank(st2))
              STOP
            ENDIF
            lutab(j,l) = real(per)
          ENDDO
          lutab(j,1) = 100.0d0 - real(persum)

        GOTO 96        
      ENDIF
*----------------------------------------------------------------------*

      READ(98,*) grassrc, barerc, fireres
      READ(98,*)

*----------------------------------------------------------------------*
* Read in sites.                                                       *
*----------------------------------------------------------------------*
      land_check = .true.

      latdel = 0.0d0
      londel = 0.0d0

      IF (sit_grd.EQ.1) THEN
*----------------------------------------------------------------------*
* Single site defined in climate file.                                 *
*----------------------------------------------------------------------*
        READ(98,'(1000a)') st1
        CALL STRIPB(st1)

        IF (n_fields(st1).EQ.1) THEN
          READ(st1,*) sites
        ELSE
          sites = 1
        ENDIF
        DO site=1,sites
          lat_lon(site,1) = xlatf
          lat_lon(site,2) = xlon0
        ENDDO
        READ(98,*)
      ELSE

        READ(98,'(1000a)') st1
        CALL STRIPB(st1)

        IF (n_fields(st1).EQ.1) THEN
          land_check = .false.

          IF (du.eq.1) THEN
            recl1 = 16
          ELSE
            recl1 = 17
          ENDIF

          st2 = 'ARGUMENT'
          IF (stcmp(st2,st1).EQ.1) THEN

            narg = narg + 1
            CALL GETARG(narg,buff1)
            DO i=1,80
              st1(i:i) = buff1(i:i)
            ENDDO

            narg = narg + 1
            CALL GETARG(narg,buff1)
            DO i=1,80
              st2(i:i) = buff1(i:i)
            ENDDO
            CALL STRIPBN(st2,site0)

            narg = narg + 1
            CALL GETARG(narg,buff1)
            DO i=1,80
              st2(i:i) = buff1(i:i)
            ENDDO
            CALL STRIPBN(st2,sitef)

            OPEN(99,FILE=st1,STATUS='OLD',FORM='FORMATTED',
     &ACCESS='DIRECT',RECL=recl1,iostat=kode)
            IF (kode.NE.0) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,*) 'File containing list of sites does not exist.
     &'
              WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
              STOP
            ENDIF

            DO site=site0,sitef
              READ(99,'(F7.3,F9.3)',REC=site) 
     &lat_lon(site-site0+1,1),lat_lon(site-site0+1,2)
            ENDDO
            sites = sitef - site0 + 1

c CLOSE added by Ghislain 15/12/03
            CLOSE(99)

            READ(98,*)

          ELSE

            OPEN(99,FILE=st1,STATUS='OLD',iostat=kode)
            IF (kode.NE.0) THEN
              WRITE(*,'('' PROGRAM TERMINATED'')')
              WRITE(*,*) 'File containing list of sites does not exist.
     &'
              WRITE(*,'('' "'',A,''"'')') st1(1:blank(st1))
              STOP
            ENDIF

* Check for countries to be run.
            READ(98,'(1000a)') st1
            CALL STRIPB(st1)

            IF (n_fields(st1).GT.0) THEN
              no_countries = n_fields(st1)
              READ(st1,*) (countries(i),i=1,no_countries)
              READ(98,*)
            ENDIF

            sites = 0
50          READ(99,*,END=60) lat_lon(sites+1,1),lat_lon(sites+1,2)
              sites = sites + 1
            GOTO 50
60          CONTINUE

c CLOSE added by Ghislain 15/12/03
            CLOSE(99)
          ENDIF

        ELSEIF ((n_fields(st1).EQ.2).OR.(n_fields(st1).EQ.3)) THEN

*----------------------------------------------------------------------*
* List of sites in the input file.                                     *
*----------------------------------------------------------------------*
          sites = 1
          READ(st1,*) lat_lon(sites,1),lat_lon(sites,2)
97        CONTINUE
            READ(98,'(1000a)') st1
            CALL STRIPB(st1)
            IF (ichar(st1(1:1)).NE.32) THEN
              sites = sites + 1
              READ(st1,*) lat_lon(sites,1),lat_lon(sites,2)
            GOTO 97
          ENDIF

        ELSEIF ((n_fields(st1).EQ.4).OR.(n_fields(st1).EQ.5)) THEN

*----------------------------------------------------------------------*
* Box of sites in input file.                                          *
*----------------------------------------------------------------------*
          READ(st1,*) lat0,latf,lon0,lonf
          IF ((lat0.GT.latf).OR.(lon0.GT.lonf)) THEN
            WRITE(*,'('' PROGRAM TERMINATED'')')
            WRITE(*,*) 'lat0 and lon0 must be less than latf and lonf.'
          ENDIF
          IF (lat0.LT.-90.0d0)  lat0 = -90.0d0
          IF (latf.GT. 90.0d0)  latf =  90.0d0
          IF (lon0.LT.-180.0d0)  lon0 = -180.0d0
          IF (lonf.GT. 180.0d0)  lonf =  180.0d0
          sites = 0
          DO i=1,int((latf-lat0)/xlatres+0.5d0)
            DO j=1,int((lonf-lon0)/xlonres+0.5d0)
              sites = sites + 1
              lat_lon(sites,1) = lat0 + (real(i-1)+0.5d0)*xlatres
              lat_lon(sites,2) = lon0 + (real(j-1)+0.5d0)*xlonres
            ENDDO
          ENDDO
          READ(98,*)
          latdel = xlatres
          londel = xlonres

        ELSEIF (n_fields(st1).GE.6) THEN

*----------------------------------------------------------------------*
* Box of sites in input file.                                          *
*----------------------------------------------------------------------*
          READ(st1,*) lat0,latf,lon0,lonf,latdel,londel
          IF ((lat0.GT.latf).OR.(lon0.GT.lonf)) THEN
            WRITE(*,'('' PROGRAM TERMINATED'')')
            WRITE(*,*) 'lat0 and lon0 must be less than latf and lonf.'
          ENDIF
          IF (lat0.LT.-90.0d0)  lat0 = -90.0d0
          IF (latf.GT. 90.0d0)  latf =  90.0d0
          IF (lon0.LT.-180.0d0)  lon0 = -180.0d0
          IF (lonf.GT. 180.0d0)  lonf =  180.0d0
          IF (latdel.LT.1.0d0/120.0d0) latdel = 1.0d0/120.0d0
          IF (londel.LT.1.0d0/120.0d0) londel = 1.0d0/120.0d0
          sites = 0
          IF (latf.GT.lat0+latdel/2.0d0) THEN
            ibox = int((latf-lat0)/latdel+0.5d0)
            iadj = latdel/2.0d0
          ELSE
            ibox = 1
            iadj = 0.0d0
          ENDIF
          IF (lonf.GT.lon0+londel/2.0d0) THEN
            jbox = int((lonf-lon0)/londel+0.5d0)
            jadj = londel/2.0d0
          ELSE
            jbox = 1
            jadj = 0.0d0
          ENDIF
          DO i=1,ibox
            DO j=1,jbox
              sites = sites + 1
              lat_lon(sites,1) = lat0 + real(i-1)*latdel + iadj
              lat_lon(sites,2) = lon0 + real(j-1)*londel + jadj
            ENDDO
          ENDDO
          READ(98,*)
        ELSE
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'Error in the number of fields for site input'
          STOP
        ENDIF

      ENDIF

*---------------------------------------------------------------------*
* Check sites against land mask and disregard when no land.           *
*---------------------------------------------------------------------*
      st2 = stoutput
      IF (land_check) THEN
        CALL LAND_SITE_CHECK(st1,st2,sites,lat_lon,latdel,londel,du,
     &stmask,no_countries)
      ENDIF

*----------------------------------------------------------------------*
* Read in soil characteristics. Defaults used when zero                *
*----------------------------------------------------------------------*
      READ(98,'(1000a)') st1
      IF ((n_fields(st1).LT.8).OR.(n_fields(st1).GT.9)) THEN
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Eight or nine fields must be present for the soil 
     &characteristics: sand %,'
        WRITE(*,*) 'silt %, bulk density and Soil Depth (cm). Set to 
     &zero if default values'
        WRITE(*,*) 'are required.'
        STOP
      ENDIF
      IF (n_fields(st1).EQ.8) THEN
        READ(st1,*) (soil_chr(i),i=1,8)
        topsl = defaulttopsl
      ELSE
        READ(st1,*) (soil_chr(i),i=1,8),topsl
        IF (topsl.LT.1.0e-6) topsl = defaulttopsl
      ENDIF
      IF ((soil_chr(5).LT.0.0d0).OR.(soil_chr(6).LT.0.0d0).OR.
     &(soil_chr(7).LT.0.0d0)) THEN
        l_b_and_c = .TRUE.
      ELSE
        l_b_and_c = .FALSE.
      ENDIF
      READ(98,'(1000a)') st1

      READ(98,'(1000a)') st1

*----------------------------------------------------------------------*
* Check for parameter adjustment
*----------------------------------------------------------------------*
      st2 = 'ARGUMENT'
      IF (stcmp(st2,st1).EQ.1) THEN

        l_parameter = .true.

        narg = narg + 1
        CALL GETARG(narg,buff1)
        DO i=1,80
          param_file(i:i) = buff1(i:i)
        ENDDO

        narg = narg + 1
        CALL GETARG(narg,buff1)
        DO i=1,80
          st2(i:i) = buff1(i:i)
        ENDDO
        CALL STRIPBN(st2,n_param_0)

        narg = narg + 1
        CALL GETARG(narg,buff1)
        DO i=1,80
          st2(i:i) = buff1(i:i)
        ENDDO
        CALL STRIPBN(st2,n_param_f)

        sites = n_param_f - n_param_0 + 1
        n_param = n_param_0 - 1
        READ(98,*)

      ELSE
        l_parameter = .false.
      ENDIF


      CLOSE(98)
*----------------------------------------------------------------------*
* End read of input file 
*----------------------------------------------------------------------*


*----------------------------------------------------------------------*
* Set up parameters for hard wired funcitonal types 'BARE' and 'CITY'. *
*----------------------------------------------------------------------*
      ftmor(1) = 1
      ftlls(1) = 0
      ftsls(1) = 0
      ftrls(1) = 0
      ftxyl(1) = 0.0d0
      ftpd(1) = 0.0d0
      ftstmx(1) = 0.0d0
      nppstore(1) = 0.0d0
      nppstorx(1) = 0.0d0
      nppstor2(1) = 0.0d0
      tsumam(1) = 0.0d0

      ftmor(2) = 1
      ftlls(2) = 0
      ftsls(2) = 0
      ftrls(2) = 0
      ftxyl(2) = 0.0d0
      ftpd(2) = 0.0d0
      ftstmx(2) = 0.0d0
      nppstore(2) = 0.0d0
      nppstorx(2) = 0.0d0
      nppstor2(2) = 0.0d0
      tsumam(2) = 0.0d0

*----------------------------------------------------------------------*
      srespm = 0.0d0

      co2(1,:,:) = 0.0d0

      DO ft=1,nft
        ftxyl(ft) =  ftxyl(ft)*1.0e-9
        ftpd(ft)  =  ftpd(ft)*1.0e+3
      ENDDO

      st1 = stoutput
      st2 = stinput
      st3 = stinit

*----------------------------------------------------------------------*
* Open snapshot files.                                                 *
*----------------------------------------------------------------------*
      IF (snp_no.GT.0) THEN
        DO i=1,snp_no
          fno = 100 + (i-1)*4
          st4=in2st(snpshts(i))
          CALL STRIPB(st4)
          OPEN(fno+1,FILE=st1(1:blank(st1))//
     &'/initbio_'//st4(1:blank(st4))//'.dat')
          OPEN(fno+2,FILE=st1(1:blank(st1))//
     &'/initcov_'//st4(1:blank(st4))//'.dat')
          OPEN(fno+3,FILE=st1(1:blank(st1))//
     &'/initppm_'//st4(1:blank(st4))//'.dat')
          OPEN(fno+4,FILE=st1(1:blank(st1))//
     &'/inithgt_'//st4(1:blank(st4))//'.dat')
        ENDDO
      ENDIF

*----------------------------------------------------------------------*
* Open default output files.                                           *
*----------------------------------------------------------------------*
      OPEN(21,FILE=st1(1:blank(st1))//'/lai.dat')
      OPEN(22,FILE=st1(1:blank(st1))//'/npp.dat')
      OPEN(23,FILE=st1(1:blank(st1))//'/scn.dat')
      OPEN(24,FILE=st1(1:blank(st1))//'/snn.dat')
      OPEN(25,FILE=st1(1:blank(st1))//'/nep.dat')
      OPEN(26,FILE=st1(1:blank(st1))//'/swc.dat')
      OPEN(27,FILE=st1(1:blank(st1))//'/biot.dat')
      OPEN(28,FILE=st1(1:blank(st1))//'/bioind.dat')
      OPEN(29,FILE=st1(1:blank(st1))//'/co2.dat')
      OPEN(30,FILE=st1(1:blank(st1))//'/dof.dat')
      OPEN(31,FILE=st1(1:blank(st1))//'/rof.dat')
      OPEN(32,FILE=st1(1:blank(st1))//'/fcn.dat')
      OPEN(33,FILE=st1(1:blank(st1))//'/nppstore.dat')
      OPEN(34,FILE=st1(1:blank(st1))//'/stembio.dat')
      OPEN(35,FILE=st1(1:blank(st1))//'/rootbio.dat')
      OPEN(36,FILE=st1(1:blank(st1))//'/leafper.dat')
      OPEN(37,FILE=st1(1:blank(st1))//'/stemper.dat')
      OPEN(38,FILE=st1(1:blank(st1))//'/rootper.dat')
      OPEN(39,FILE=st1(1:blank(st1))//'/sresp.dat')
      OPEN(40,FILE=st1(1:blank(st1))//'/evt.dat')
      OPEN(41,FILE=st1(1:blank(st1))//'/gpp.dat')
      OPEN(42,FILE=st1(1:blank(st1))//'/lch.dat')
      OPEN(43,FILE=st1(1:blank(st1))//'/prc.dat')
      OPEN(44,FILE=st1(1:blank(st1))//'/nbp.dat')
      OPEN(45,FILE=st1(1:blank(st1))//'/trn.dat')
      OPEN(46,FILE=st1(1:blank(st1))//'/fab.dat')
      OPEN(47,FILE=st1(1:blank(st1))//'/tmp.dat')
      OPEN(48,FILE=st1(1:blank(st1))//'/hum.dat')
      OPEN(49,FILE=st1(1:blank(st1))//'/ppm.dat')
      OPEN(50,FILE=st1(1:blank(st1))//'/pet.dat')
      OPEN(51,FILE=st1(1:blank(st1))//'/surface_lit.dat')
      OPEN(52,FILE=st1(1:blank(st1))//'/subsurface_lit.dat')
      OPEN(53,FILE=st1(1:blank(st1))//'/tleaf_n.dat')
      OPEN(54,FILE=st1(1:blank(st1))//'/tleaf_p.dat')
      OPEN(55,FILE=st1(1:blank(st1))//'/sla.dat')
      OPEN(56,FILE=st1(1:blank(st1))//'/max_hgt.dat')
      OPEN(57,FILE=st1(1:blank(st1))//'/vcmax_type.dat')
      OPEN(58,FILE=st1(1:blank(st1))//'/anlfn.dat')
      OPEN(59,FILE=st1(1:blank(st1))//'/antlfn.dat')
      OPEN(60,FILE=st1(1:blank(st1))//'/anvcmax.dat')
      OPEN(61,FILE=st1(1:blank(st1))//'/anjmax.dat')
      OPEN(62,FILE=st1(1:blank(st1))//'/presp.dat')
      OPEN(63,FILE=st1(1:blank(st1))//'/swr.dat')
      OPEN(64,FILE=st1(1:blank(st1))//'/mgresp.dat')
      OPEN(65,FILE=st1(1:blank(st1))//'/soilcn.dat')
      OPEN(66,FILE=st1(1:blank(st1))//'/kg_beta.dat')
      OPEN(67,FILE=st1(1:blank(st1))//'/can_clump.dat')
      OPEN(68,FILE=st1(1:blank(st1))//'/qdir.dat')
      OPEN(69,FILE=st1(1:blank(st1))//'/qtotal.dat')
      OPEN(601,FILE=st1(1:blank(st1))//'/field_capacity.dat')
      OPEN(602,FILE=st1(1:blank(st1))//'/wilting_point.dat')
      OPEN(603,FILE=st1(1:blank(st1))//'/abg_litter.dat')
      OPEN(604,FILE=st1(1:blank(st1))//'/blg_c.dat')
      OPEN(605,FILE=st1(1:blank(st1))//'/leafc.dat')
      OPEN(606,FILE=st1(1:blank(st1))//'/lulccc.dat')

*----------------------------------------------------------------------*
* Open optional yearly cover, biomass, budburst and senescence files.  *
*----------------------------------------------------------------------*
      iofn = 200
      iofngft = iofn
      IF (outyears2.GT.0) THEN
      DO ft=1,nft
        IF (out_cov) THEN
          iofn = iofn + 1
          OPEN(iofn,FILE=st1(1:blank(st1))//'/cov_'//
     &fttags(ft)(1:blank(fttags(ft)))//'.dat')
        ENDIF
        IF (out_bio) THEN
          iofn = iofn + 1
          OPEN(iofn,FILE=st1(1:blank(st1))//'/bio_'//
     &fttags(ft)(1:blank(fttags(ft)))//'.dat')
        ENDIF
        IF (out_bud) THEN
          iofn = iofn + 1
          OPEN(iofn,FILE=st1(1:blank(st1))//'/bud_'//
     &fttags(ft)(1:blank(fttags(ft)))//'.dat')
        ENDIF
        IF (out_sen) THEN
          iofn = iofn + 1

          OPEN(iofn,FILE=st1(1:blank(st1))//'/sen_'//
     &fttags(ft)(1:blank(fttags(ft)))//'.dat')
        ENDIF
      ENDDO
      ENDIF

*----------------------------------------------------------------------*
* Open input files for state vector.                                   *
*----------------------------------------------------------------------*
      IF (.NOT.initise) THEN
        OPEN(70,FILE=st3(1:blank(st3))//'/init.dat',status='old',
     &iostat=kode)
        IF (kode.NE.0) THEN
          WRITE(*,'('' PROGRAM TERMINATED'')')
          WRITE(*,*) 'Initialisation directory/file does not exist.'
          WRITE(*,'('' "'',A,''"'')') st3(1:blank(st3))
          STOP
        ENDIF
        OPEN(71,FILE=st3(1:blank(st3))//'/initbio.dat')
        OPEN(72,FILE=st3(1:blank(st3))//'/initcov.dat')
        OPEN(73,FILE=st3(1:blank(st3))//'/initppm.dat')
        OPEN(74,FILE=st3(1:blank(st3))//'/inithgt.dat')
        OPEN(75,FILE=st3(1:blank(st3))//'/initwdt.dat')
        OPEN(76,FILE=st3(1:blank(st3))//'/initleaf.dat')
        OPEN(77,FILE=st3(1:blank(st3))//'/initstem.dat')
        OPEN(78,FILE=st3(1:blank(st3))//'/initroot.dat')
        OPEN(79,FILE=st3(1:blank(st3))//'/initmisc.dat')
      ENDIF

*----------------------------------------------------------------------*
* Open output files for state vector.                                  *
*----------------------------------------------------------------------*
      OPEN(80,FILE=st1(1:blank(st1))//'/init.dat')
      IF (initiseo) THEN
        OPEN(81,FILE=st1(1:blank(st1))//'/initbio.dat')
        OPEN(82,FILE=st1(1:blank(st1))//'/initcov.dat')
        OPEN(83,FILE=st1(1:blank(st1))//'/initppm.dat')
        OPEN(84,FILE=st1(1:blank(st1))//'/inithgt.dat')
        OPEN(85,FILE=st1(1:blank(st1))//'/initwdt.dat')
        OPEN(86,FILE=st1(1:blank(st1))//'/initleaf.dat')
        OPEN(87,FILE=st1(1:blank(st1))//'/initstem.dat')
        OPEN(88,FILE=st1(1:blank(st1))//'/initroot.dat')
        OPEN(89,FILE=st1(1:blank(st1))//'/initmisc.dat')
      ENDIF

*----------------------------------------------------------------------*
*  Read co2 file.                                                      *
*----------------------------------------------------------------------*
      !print*, daily_co2,yr0,yrf
      CALL READCO2 (stco2,co2,daily_co2,yr0a,yrfa,year0set,spinl,
     &nyears)
*----------------------------------------------------------------------*

      site_dat = 0

*----------------------------------------------------------------------*
*  site loop                              
*----------------------------------------------------------------------*
!PCM4 for running a single site in the main run from a gridded spinup;
!PCM4 (other changes below)
!      DO site=10,10
      DO site=1,sites
*----------------------------------------------------------------------*
* closed_loop_ft:                                                      *
*    standard setting is .FALSE.                                       *
*    But, this is set to .TRUE. after the first initialization of      * 
*    ftprop with cluse, so that the gross transitions will update the  *
*    ftprop(ft) values each year instead of using the values from the  *
*    cluse table from the net land use transitions                     *
*----------------------------------------------------------------------*
        closed_loop_ft = .FALSE.

        speedc = xspeedc
        swcnew = 0.0d0

        IF (l_parameter) THEN
          lat = lat_lon(1,1)
          lon = lat_lon(1,2)
        ELSE
          lat = lat_lon(site,1)
          lon = lat_lon(site,2)
        ENDIF

        WRITE(*,*)site,site0+site-1,lat,lon

        IF (abs(xseed1).EQ.0) THEN
          IF (site.EQ.1) THEN
            seed1 = int(SETARANDOM()*10000.0d0+.5)
            seed2 = 2*seed1
            seed3 = 3*seed1
          ENDIF
        ELSEIF (xseed1.LT.0.0) THEN
          IF (site.EQ.1) THEN
            seed1 = -xseed1
            seed2 = 2*seed1
            seed3 = 3*seed1
          ENDIF
        ELSE
          seed1 = xseed1
          seed2 = 2*seed1
          seed3 = 3*seed1
        ENDIF

      l_countries = .TRUE.
      IF (no_countries.GT.0) THEN
        l_countries = .FALSE.
        CALL COUNTRY(stmask,lat,lon,country_name,country_id,.true.)
        DO i=1,no_countries
          IF (blank(country_name).EQ.blank(countries(i))) THEN
            check_c = .TRUE.
            DO j=1,blank(country_name)
              IF (country_name(j:j).NE.countries(i)(j:j)) 
     &check_c = .FALSE.
            ENDDO
          ELSE
            check_c = .FALSE.
          ENDIF
          IF (check_c) l_countries = .TRUE.
        ENDDO
        CALL COUNTRY(stmask,lat,lon,country_name,country_id,.false.)
        DO i=1,no_countries
          IF (blank(country_name).EQ.blank(countries(i))) THEN
            check_c = .TRUE.
            DO j=1,blank(country_name)
              IF (country_name(j:j).NE.countries(i)(j:j)) 
     &check_c = .FALSE.
            ENDDO
          ELSE
            check_c = .FALSE.
          ENDIF
          IF (check_c) l_countries = .TRUE.
        ENDDO
      ENDIF

      IF (l_countries) THEN
*----------------------------------------------------------------------*
* Read in climate.                                                     *
*----------------------------------------------------------------------*
      !print*, 'start clim'
      IF (clim_type.EQ.1) THEN
*----------------------------------------------------------------------*
* DAILY Gridded.                                                       *
*----------------------------------------------------------------------*
        CALL EX_CLIM(st2,lat,lon,xlatf,xlatres,xlatresn,xlon0,xlonres,
     &xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,isite,xyear0,xyearf,
     &siteno,du,xswrv,read_par)
        if(siteno.NE.0) THEN !PCM
          l_clim = .TRUE.    !PCM
          l_stats = .TRUE.   !PCM: we don't have or need this stats data,
                             !PCM  but this flag needs to be set
         ELSE                !PCM
          l_clim = .FALSE.   !PCM
          l_stats = .FALSE.   !PCM
         ENDIF
C PCM2        WRITE(*,*) 'aaaaaaaaa 1st year temperature'
C PCM2        DO mnth=1,12
C PCM2          WRITE(*,'(30f5.0)') xtmpv(1,mnth,1:30) 
C PCM2        ENDDO
C PCM2        WRITE(*,*) 'bbbbbbbbb'
        withcloudcover=.FALSE.
      ELSEIF (clim_type.EQ.2) THEN
*----------------------------------------------------------------------*
* MONTHLY Gridded.                                                     *
*----------------------------------------------------------------------*
        CALL EX_CLIM_WEATHER_GENERATOR(st2,ststats,lat,lon,xlatf,
     &xlatres,xlatresn,xlon0,xlonres,xlonresn,yr0,yrf,xtmpv,xhumv,xprcv,
     &xcldv,isite,xyear0,xyearf,du,seed1,seed2,seed3,l_clim,l_stats,
     &xswrv,read_par)
        withcloudcover=.TRUE.
      ELSEIF (clim_type.EQ.3) THEN
*----------------------------------------------------------------------*
* DAILY Site                                                           *
*----------------------------------------------------------------------*
        CALL EX_CLIM_SITE(st2,yr0,yrf,xtmpv,xhumv,xprcv,xyear0,xyearf,
     &xswrv,read_par)
        withcloudcover=.FALSE.
        siteno = 1
        l_clim = .TRUE.
        l_stats = .TRUE.
      ELSEIF (clim_type.EQ.4) THEN
*----------------------------------------------------------------------*
* MONTHLY Site                                                           *
*----------------------------------------------------------------------*
        CALL EX_CLIM_SITE_MONTH(st2,yr0,yrf,xtmpv,xhumv,xprcv,xcldv,
     &xyear0,xyearf,xswrv,read_par)
        withcloudcover=.FALSE.
        siteno = 1
        CALL GENERATE_MONTHLY(ststats,yr0,yrf,xtmpv,xhumv,xprcv,
     &seed1,seed2,seed3,xswrv,read_par)
        l_clim = .TRUE.
        l_stats = .TRUE.
      ELSE
        WRITE(*,'('' PROGRAM TERMINATED'')')
        WRITE(*,*) 'Error defining climate to read'
        STOP
      ENDIF
      !print*, 'end clim'

*----------------------------------------------------------------------*
* Read in soil parameters.                                             *
*----------------------------------------------------------------------*
! must be read from a file that contains lat lon and 10 soil parameter fields (i.e. 2 more than the original sdgvm) 
      IF (l_clim) THEN
        CALL SOIL_PARAMETERS(ts,tsi,tc,bulk,orgc,wilt,field,sat,dep,
     &sttxdp,lat,lon,soil_chr,soil_chr2,du,l_soil,soilCtoN,soilp)
      !print*, field, wilt 
      !print*, 'end soil'
      ENDIF

*----------------------------------------------------------------------*
* Read in canopy clumping index from a map                             *
*----------------------------------------------------------------------*
      IF(read_clump.eq.0)THEN
      !if no clumping set canopy clumping index to 1
        ftcan_clump(:) = 1
      ELSEIF(read_clump.eq.2)THEN
        CALL EX_CLUMP(stlu,lat,lon,map_clump,du)
        ftcan_clump(:) = map_clump
      ENDIF

*----------------------------------------------------------------------*
* IF in parameter adjustment mode read in parameters and perform       *
* necessary adjustments.                                               *
*----------------------------------------------------------------------*
      IF (l_parameter) THEN
        n_param = n_param + 1
        CALL PARAMETER_ADJUSTMENT(param_file,n_param)
        ts = p_sand
        tsi = p_silt
        tc = 100.0 - ts - tsi
        bulk = p_bulk
        orgc = p_orgc
        dep = p_dep
        topsl = p_topsl
      ENDIF
      !print*, 'end param'
*----------------------------------------------------------------------*
* Read in landuse/cover.                                               *
*----------------------------------------------------------------------*
c     temporaire !
      icontinuouslanduse=1
      !print*, luse(:)
      IF (ilanduse.EQ.0 .OR. (ilanduse.GE.3 .AND. ilanduse.LE.6)) THEN
        IF (icontinuouslanduse.EQ.0) THEN
          CALL EX_LU(stlu,lat,lon,luse,yr0,yrf,du)
c     create the continuous land use (cluse)
          DO year=yr0,yrf
            DO ft=1,nft
              cluse(ft,year-yr0+1) = lutab(luse(year-yr0+1),ft)
              !print*, year,yr0
              !print*, luse(year-yr0+1)
              !print*, lutab(luse(year-yr0+1),ft)
              !print*, cluse(ft,year-yr0+1) 
            ENDDO
          ENDDO
        ELSE            
C          IF(ilanduse.EQ.3 .OR. ilanduse.EQ.4) THEN
C           SYR = 1700
C            NYR = 322
C         ELSE IF(ilanduse.EQ.5 .OR. ilanduse.EQ.6) THEN !preindustrial, constant
C           SYR = 1700
C           NYR = 1 
C         ELSE IF(ilanduse.EQ.0) THEN !SYR and NYR not used and not defined here
C           SYR = -1 
C           NYR = -1 
C         ENDIF
          lutab2(:,:)                      =  0.0
          IF(ilanduse.GE.3 .AND. ilanduse.LE.6) THEN
C The following order is order of land types in the input.dat file
            !aggmap_SDGVM_to_aggHyde(1)     =  0 
            !aggmap_SDGVM_to_aggHyde(2:6)   =  1 
            !aggmap_SDGVM_to_aggHyde(7:8)   =  2 
            !aggmap_SDGVM_to_aggHyde(9)     =  5 
            !aggmap_SDGVM_to_aggHyde(10)    =  6 
            !aggmap_SDGVM_to_aggHyde(11:15) =  3 
            !aggmap_SDGVM_to_aggHyde(16:17) =  4 
C The following ordering is the order of ft's in the input.dat file
            aggmap_SDGVM_to_aggHyde(1)     =  0 
            aggmap_SDGVM_to_aggHyde(2)     =  7 
            aggmap_SDGVM_to_aggHyde(3:4)   =  2 
            aggmap_SDGVM_to_aggHyde(5)     =  5 
            aggmap_SDGVM_to_aggHyde(6)     =  6 
            aggmap_SDGVM_to_aggHyde(7:8)   =  4 
            aggmap_SDGVM_to_aggHyde(9:12)  =  1 
            aggmap_SDGVM_to_aggHyde(13:16) =  3 
            DO at=1,maxn_at
               lutab2(at,at)               = 100.0 !PCM: Kluge: assumes 1<->1 mapping of 'at' to classes for now 
            ENDDO
          ELSE
            aggmap_SDGVM_to_aggHyde(:)     =  0
          END IF

          CALL EX_CLU(stlu,lat,lon,nft,lutab,cluse,du,l_lu,
     &yr0a,yrfa,year0set,spinl,ilanduse,SYR,NYR,lutab2,
     &stpname,stpname_t,stwdg,cluse2,cluseh,debug)

      !loop added for testing purposes
          DO ft=1,nft
            sum_cov_test(ft) = 0.0
            DO age=1,ftmor(ft)
              sum_cov_test(ft) = sum_cov_test(ft)
     &                           + cluse(ft,age)/100.0d0/real(ftmor(ft))
            ENDDO
          ENDDO

          write(*,FMT="(A)") 'SS1'
          IF (ilanduse.GE.3 .AND. ilanduse.LE.6) THEN
            IF(debug .EQV. .TRUE.) THEN
             WRITE(*,*) 'BARE       CITY      C3p        C4p        ',
     &    'C3crop     ','C4crop     C3s       C4s        Ev_Bp      ',
     &    'Ev_Np      ','Dc_Bp      Dc_Np     Ev_Bs      Ev_Ns      ',
     &    'Dc_Bs      ','Dc_Ns      '

             write(*,FMT="(16F11.7)") cluse(1:nft,1)
             !write(*,FMT="(16F7.3)") sum_cov_test(1:nft)
             write(*,FMT="(A)") 'SS2'
             write(*,FMT="(8A14)")
     &      '     ','     primf', '     primn', '     secdf',
     &              '     secdn', '    c3crop', '    c4crop',
     &              '     urban', '      harv' 

             write(*,FMT="(A,8E14.7)") ' primf',cluse2(1,:,1),
     &           cluseh(1,1)
             write(*,FMT="(A,8E14.7)") ' primn',cluse2(2,:,1),
     &           cluseh(2,1)
             write(*,FMT="(A,8E14.7)") ' secdf',cluse2(3,:,1),
     &           cluseh(3,1)
             write(*,FMT="(A,8E14.7)") ' secdn',cluse2(4,:,1),
     &           cluseh(4,1)
             write(*,FMT="(A,8E14.7)") 'c3crop',cluse2(5,:,1),
     &           cluseh(5,1)
             write(*,FMT="(A,8E14.7)") 'c4crop',cluse2(6,:,1),
     &           cluseh(6,1)
             write(*,FMT="(A,8E14.7)") ' urban',cluse2(7,:,1),
     &           cluseh(7,1)
           ENDIF
          ELSE
!             WRITE(*,*) ' BARE   CITY   C3     C4     C3crop ',
!     &'C4crop Ev_Bl  Ev_Nl  Dc_Bl ',
!     &'Dc_Nl '
             !write(*,FMT="(10F7.3)") cluse(1:nft,1)
             !write(*,FMT="(10F7.3)") sum_cov_test(1:nft)
          ENDIF
        ENDIF
      ELSEIF (ilanduse.EQ.1) THEN
        l_lu = .TRUE.
        !yr0a = yr0
        !yrfa = yrf
        !if( (co2const.lt.0.0) .and. (spinl.gt.0) ) then 
        !  if (spinl.lt.nyears) then
        !    yr0a = yr0 - spinl
        !  else
        !    !yr0a = years(1)
        !    yrfa = yr0a + spinl - 1
        !  endif
        !endif
        DO year=yr0a,yrfa
          DO ft=1,nft
            cluse(ft,year-yr0a+1) = lutab(luse(year-yr0a+1),ft)
          ENDDO
        ENDDO
      ELSEIF (ilanduse.EQ.2) THEN
        WRITE(*,*) 'Checking natural vegetation types exist:'
        WRITE(*,*) 'BARE CITY C3 C4 Ev_Bl Ev_Nl Dc_Bl Dc_Nl.'
        l_lu = .TRUE.
        st2 = 'BARE'
        nat_map(1) = ntags(fttags,st2)
        st2 = 'CITY'
        nat_map(2) = ntags(fttags,st2)
        st2 = 'C3'
        nat_map(3) = ntags(fttags,st2)
        st2 = 'C4'
        nat_map(4) = ntags(fttags,st2)
        st2 = 'Ev_Bl'
        nat_map(5) = ntags(fttags,st2)
        st2 = 'Ev_Nl'
        nat_map(6) = ntags(fttags,st2)
        st2 = 'Dc_Bl'
        nat_map(7) = ntags(fttags,st2)
        st2 = 'Dc_Nl'
        nat_map(8) = ntags(fttags,st2)
        DO year=yr0a,yrfa
          DO ft=1,nft
            cluse(ft,year-yr0+1) = 0.0d0
          ENDDO
        ENDDO
        cluse(1,1) = 100.0d0
      ENDIF
      !print*, 'end landuse, ', site
      !print*, l_clim,l_stats,l_soil(1),l_soil(3),l_soil(8),l_lu, 
      !&l_countries

*----------------------------------------------------------------------*
      ENDIF
*----------------------------------------------------------------------*
* End of countries check.                                              *
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
* Start of climate data exists 'if' statement.                         *
*----------------------------------------------------------------------*
        IF((l_clim).AND.(l_stats).AND.(l_soil(1)).AND.(l_soil(3)).AND.
     &(l_soil(8)).AND.(l_lu).AND.(l_countries)) THEN

        site_dat = site_dat + 1
        IF (mod(site_dat,max(site_out,1)).EQ.min(1,site_out)-1) 
c     &WRITE(*,'( '' Site no. '',i5,'', Lat ='',f7.3,'', Lon ='',f9.3,
c     &'' CO2 = '',f9.3)') 
c     &site_dat,lat,lon,ca
     &WRITE(*,'( '' Site no. '',i5,'', Lat ='',f7.3,'', Lon ='',f9.3)') 
     &site_dat,lat,lon

*----------------------------------------------------------------------*
* Write lat lon for general output files.                              *
*----------------------------------------------------------------------*
        DO i=21,69
          WRITE(i,'(f7.3,f9.3,$)') lat,lon
        ENDDO
        DO i=601,606
          WRITE(i,'(f7.3,f9.3,$)') lat,lon
        ENDDO

*----------------------------------------------------------------------*
* Write lat lon for optional yearly ft files                           *
*----------------------------------------------------------------------*
        iofn = iofngft
        IF (outyears2.GT.0) THEN
        DO ft=1,nft
          IF (out_cov) THEN
            iofn = iofn + 1
            WRITE(iofn,'(f7.3,f9.3,$)') lat,lon
          ENDIF
          IF (out_bio) THEN
            iofn = iofn + 1
            WRITE(iofn,'(f7.3,f9.3,$)') lat,lon
          ENDIF
          IF (out_bud) THEN
            iofn = iofn + 1
            WRITE(iofn,'(f7.3,f9.3,$)') lat,lon
          ENDIF
          IF (out_sen) THEN
            iofn = iofn + 1
            WRITE(iofn,'(f7.3,f9.3,$)') lat,lon
          ENDIF
        ENDDO
        ENDIF

*----------------------------------------------------------------------*
* Write lat,lon for monthly/daily output files.                        *
*----------------------------------------------------------------------*
        iofn = 400
        IF (outyears1.GT.0) THEN
        IF ((oymd.EQ.1).OR.(oymd.EQ.3)) THEN
          DO i=1,douts
            IF (otagsn(i).EQ.1) THEN
              iofn = iofn + 1
              WRITE(iofn,'(f7.3,f9.3)') lat,lon
            ENDIF
          ENDDO
        ENDIF
        IF ((oymd.EQ.2).OR.(oymd.EQ.3)) THEN
          DO i=1,douts
            IF (otagsn(i).EQ.1) THEN
              iofn = iofn + 1
              WRITE(iofn,'(f7.3,f9.3)') lat,lon
            ENDIF
          ENDDO
        ENDIF
        ENDIF

*----------------------------------------------------------------------*
* Write lat,lon for monthly/daily ft output files.                     *
*----------------------------------------------------------------------*
        IF (outyears2.GT.0) THEN
        IF ((oymdft.EQ.1).OR.(oymdft.EQ.3)) THEN
          DO i=1,douts
            IF (otagsnft(i).EQ.1) THEN
              DO ft=1,nft
                iofn = iofn + 1
                WRITE(iofn,'(f7.3,f9.3)') lat,lon
             ENDDO
            ENDIF
          ENDDO
        ENDIF
        IF ((oymdft.EQ.2).OR.(oymdft.EQ.3)) THEN
          DO i=1,douts
            IF (otagsnft(i).EQ.1) THEN
              DO ft=1,nft
                iofn = iofn + 1
                WRITE(iofn,'(f7.3,f9.3)') lat,lon
              ENDDO
            ENDIF
          ENDDO
        ENDIF
        ENDIF

*----------------------------------------------------------------------*
* Initialisation of water parameters
*----------------------------------------------------------------------*
      
      CALL WSPARAM(l_b_and_c,wilt,field,sat,nupc,awl,kd,kx,nci,infix,
     &adp,sfc,sw,sswc,ts,tc,tsi,bulk,orgc,dep,topsl,wtfc,wtwp,
     &wtswc,l_parameter)

*----------------------------------------------------------------------*
* Write site info to 'site_info.dat'.                                  *
*----------------------------------------------------------------------*
      CALL CO2_0_F(co20,co2f,yearv,yr0,spinl,co2,co2const,nyears)
      
      CALL COUNTRY(stmask,lat,lon,country_name,country_id,l_regional)
      
      WRITE(12,'(a15,i6,f9.3,f9.3,1x,2f6.1,1x,2f7.1,f7.3,f7.2,1x,
     &3f7.3,f7.1,1x,f7.3,1x,200(2f6.1,1x))')
     & country_name,country_id,lat,lon,co20,co2f,ts,tsi,bulk,orgc,
     &wtwp,wtfc,wtswc,dep,soilCtoN,
     &(cluse(ft,yearv(1)-yr0+1),cluse(ft,yearv(nyears)-yr0+1),ft=1,nft)

*----------------------------------------------------------------------*
* Initialise biomass and cover arrays.                                 *
*----------------------------------------------------------------------*
        IF (initise) THEN
          DO ft=1,nft
            IF (cluse(ft,1).GT.0.0d0) THEN
              DO age=1,ftmor(ft)
                cov(age,ft) = cluse(ft,1)/100.0d0/real(ftmor(ft))
                bio(age,1,ft) = 1000.0d0
                bio(age,2,ft) = 100.0d0
              ENDDO
              bio(1,1,1) = 0.0d0
              bio(1,2,1) = 0.0d0
              DO age=1,ftmor(ft)
                ppm(age,ft) = 0.001d0
                hgt(age,ft) = 5.0d0
                wdt(age,ft) = 0.5d0
              ENDDO
              nppstore(ft) = ftstmx(ft)
              nppstorx(ft) = ftstmx(ft)
              nppstor2(ft) = ftstmx(ft)
            ELSE
              DO age=1,ftmor(ft)
                cov(age,ft) = 0.0d0
                bio(age,1,ft) = 0.0d0
                bio(age,2,ft) = 0.0d0
              ENDDO
              bio(1,1,1) = 0.0d0
              bio(1,2,1) = 0.0d0
              DO age=1,ftmor(ft)
                ppm(age,ft) = 0.0d0
                hgt(age,ft) = 0.0d0
                wdt(age,ft) = 0.0d0
              ENDDO
              nppstore(ft) = 0.0d0
              nppstorx(ft) = 0.0d0
              nppstor2(ft) = 0.0d0
            ENDIF
            nppstore(2) = ftstmx(2)
            nppstore(3) = ftstmx(3)
            nppstorx(2) = ftstmx(2)
            nppstorx(3) = ftstmx(3)
            nppstor2(2) = ftstmx(2)
            nppstor2(3) = ftstmx(3)

            slc(ft) = 0.0d0
            rlc(ft) = 0.0d0
            sln(ft) = 0.0d0
            rln(ft) = 0.0d0

            DO day=1,ftlls(ft)
              leafdp(day,ft) = 0.0d0
            ENDDO
            DO day=1,ftsls(ft)
              stemdp(day,ft) = 0.0d0
            ENDDO
            DO day=1,ftrls(ft)
              rootdp(day,ft) = 0.0d0
            ENDDO
            DO day=1,360
              sumadp(day,ft) = 0.0d0
            ENDDO
            lai(ft) = 0.0d0
            npp(ft) = 400.0d0
            nps(ft) = 40.0d0
            npl(ft) = 40.0d0
            evp(ft) = 1.0d0
            dof(ft) = 0.0
            maxlai(ft) = 2.0d0

            bb(ft) = 0
            ss(ft) = 0
            bbgs(ft) = 0
            dsbb(ft) = 0

            DO day=1,30
              sm_trig(day,ft) = 0.0d0
            ENDDO

            tsumam(ft) = 0.0d0
            stemfr(ft) =10.0d0
          ENDDO ! ft

          is1 = 5.0d0
          is2 = 50.0d0
          is3 = 50.0d0

          is4 = 50.0d0
          isn = 0.0d0
          ilsn = 0.0d0
          DO i=1,8
            ic0(i) = 1000.0d0
            in0(i) = 100.0d0
          ENDDO
          ic0(7) = 5000.0d0
          ic0(8) =  5000.0d0
          iminn(1) = 0.0d0
          iminn(2) = 0.0d0
          iminn(3) = 0.0d0
          dslc = 0.0d0
          drlc = 0.0d0
          dsln = 0.0d0
          drln = 0.0d0
          DO day=1,200
            tmem(day) = real(xtmpv(1,(day-1)/30+1,mod(day-1,30)+1))
     &/100.0d0
          ENDDO
          DO ft=1,nft
            s1(ft) = is1
            s2(ft) = is2
            s3(ft) = is3
            s4(ft) = is4
            sn(ft) = isn
            lsn(ft) = ilsn
          ENDDO
        
          !simple initiations of mean climate and soil variables
          !only used if ncalc_type = >2 or vcmax_type = >2  
          !XXXX
          maprc_init     = 1300.0d0
          maswr_init     = 150.0d0
          soilcn_init    = 15.0d0 
          soilp_init     = 569.20d0 !median value from Yang et al 2012 global soil P map converted to ug P g-1 soil using bulk density data from islscp2
          matmp_init     = 10.0d0
          matmp_max_init = 15.0d0
          matmp_min_init = 5.0d0
          map_days_init  = 5.0d0
          mahum_init     = 75.0d0
          masoilc_init   = 20.0d0
          masoilw_init   = 0.20d0

          ce_t           = 15.d0
          ce_rh          = 75.d0 
          ce_ci          = 28.d0
          ce_ga          = 0.1d0
          ce_light       = 5.0e-4 
          ce_maxlight    = 1.0e-3 

          DO ft=3,nft
            vcmax(:,ft)  = 10.d0 
            jmax(:,ft)   = 20.d0
            pnlc(:,ft)   = 0.1d0
            enzs(:,ft)   = 1.d0
          ENDDO

        ELSE
      
!PCM4 Comment this out for running a single site in the main run from a gridded spinup
         IF ((abs(lat-zlat).GT.0.001).OR.
     &(abs(lon-zlon).GT.0.001)) THEN

!PCM4 for running a single site in the main run from a gridded spinup
!         DO WHILE ((abs(lat-zlat).GT.0.001).OR.
!     &(abs(lon-zlon).GT.0.001))

          READ(70,*) zlat,zlon
          IF ((abs(lat-zlat).GT.0.001).OR.
     &(abs(lon-zlon).GT.0.001)) THEN
            WRITE(*,'("Error lat and lon dont match.")')
            WRITE(*,*) lat,lon
            WRITE(*,*) zlat,zlon
!PCM4 comment the STOP out for running a single site in the main run from a gridded spinup
            STOP
          ENDIF
          READ(70,*) (zs1(ft),ft=1,nft)
          READ(70,*) (zs2(ft),ft=1,nft)
          READ(70,*) (zs3(ft),ft=1,nft)
          READ(70,*) (zs4(ft),ft=1,nft)
          READ(70,*) (zsn(ft),ft=1,nft)
          READ(70,*) (zlsn(ft),ft=1,nft)
          READ(70,*) (zic0(i),i=1,8)
          READ(70,*) (zin0(i),i=1,8)
          READ(70,*) (ziminn(i),i=1,3)
          READ(70,*) zdslc
          READ(70,*) zdrlc
          READ(70,*) zdsln
          READ(70,*) zdrln
          READ(70,*) zlai(4)
          READ(70,*) znpp(4)

          READ(72,*) zlat,zlon,((zcov(i,j),i=1,ftmor(j)),j=2,nft)
          READ(71,*) 
     &zlat,zlon,(((zbio(i,j,k),j=1,2),i=1,ftmor(k)),k=2,nft)
          READ(73,*) zlat,zlon,((zppm(i,j),i=1,ftmor(j)),j=2,nft)
          READ(74,*) zlat,zlon,((zhgt(i,j),i=1,ftmor(j)),j=2,nft)
          READ(75,*) zlat,zlon,((zwdt(i,j),i=1,ftmor(j)),j=2,nft)
          READ(76,*) zlat,zlon,((zleafdp(i,j),i=1,ftlls(j)),j=2,nft)
          READ(77,*) zlat,zlon,((zstemdp(i,j),i=1,ftsls(j)),j=2,nft)
          READ(78,*) zlat,zlon,((zrootdp(i,j),i=1,ftrls(j)),j=2,nft)

          READ(79,*) zlat,zlon
          READ(79,*) (zlai(i),i=1,nft)
          READ(79,*) (znpp(i),i=1,nft)
          READ(79,*) (znps(i),i=1,nft)
          READ(79,*) (znpl(i),i=1,nft)
          READ(79,*) (zevp(i),i=1,nft)
          READ(79,*) (zdof(i),i=1,nft)
          READ(79,*) (zslc(i),i=1,nft)
          READ(79,*) (zrlc(i),i=1,nft)
          READ(79,*) (zsln(i),i=1,nft)
          READ(79,*) (zrln(i),i=1,nft)
          READ(79,*) (znppstore(i),i=1,nft)
          READ(79,*) (znppstorx(i),i=1,nft)
          READ(79,*) (znppstor2(i),i=1,nft)
          READ(79,*) (zbb(i),i=1,nft)
          READ(79,*) (zbbgs(i),i=1,nft)
          READ(79,*) (zdsbb(i),i=1,nft)
          READ(79,*) (zmaxlai(i),i=1,nft)
          READ(79,*) (ztmem(i),i=1,200)
          READ(79,*) maprc_init,maswr_init,soilcn_init,soilp_init,
     &matmp_init,matmp_max_init,matmp_min_init,map_days_init,mahum_init,
     &masoilc_init,masoilw_init

          DO j=1,nft
            READ(79,*) (zsm_trig(i,j),i=1,30)
          ENDDO
          DO ft=3,nft
            READ(79,*) (zsumadp(i,ft),i=1,360)
          ENDDO
          READ(79,*) (zstemfr(ft),ft=1,nft)

          READ(79,*) ce_t 
          READ(79,*) ce_rh       
          DO ft=3,nft
            READ(79,*) (ce_ci(:,i,ft),      i=1,12)
            READ(79,*) (ce_ga(:,i,ft),      i=1,12)       
            READ(79,*) (ce_light(:,i,ft),   i=1,12)
            READ(79,*) (ce_maxlight(:,i,ft),i=1,12)
            READ(79,*) (vcmax(i,ft),i=1,12)
            READ(79,*) (jmax(i,ft),i=1,12)
            READ(79,*) (pnlc(i,ft),i=1,12)
            READ(79,*) (enzs(i,ft),i=1,12)
          ENDDO
!PCM4 comment the ENDIF out for running a single site in the main run from a gridded spinup
          ENDIF
!PCM4 uncomment the ENDDO for running a single site in the main run from a gridded spinup
!          ENDDO !PCM4

          !print*, vcmax(1,:)

          DO ft=1,nft
            s1(ft) = zs1(ft)
            s2(ft) = zs2(ft)
            s3(ft) = zs3(ft)
            s4(ft) = zs4(ft)
            sn(ft) = zsn(ft)
            lsn(ft) = zlsn(ft)
            chill(ft) = 0
            dschill(ft) = 0
          ENDDO
          DO i=1,8
            ic0(i) = zic0(i)
            in0(i) = zin0(i)
          ENDDO
          f0(:,:) = -9999.0d0 ! - could add f0 to the system state files

          DO i=1,3
            iminn(i) = ziminn(i)
          ENDDO
          dslc = zdslc
          drlc = zdrlc
          dsln = zdsln
          drln = zdrln
          lai(4) = zlai(4)
          npp(4) = znpp(4)

          DO j=2,nft
            DO i=1,ftmor(j)
              cov(i,j) = zcov(i,j)
              bio(i,1,j) = zbio(i,1,j)
              bio(i,2,j) = zbio(i,2,j)
              ppm(i,j) = zppm(i,j)
              hgt(i,j) = zhgt(i,j)
              wdt(i,j) = zwdt(i,j)
            ENDDO
          ENDDO
          DO j=2,nft
            DO i=1,ftlls(j)
              leafdp(i,j) = zleafdp(i,j)
            ENDDO
          ENDDO
          DO j=2,nft
            DO i=1,ftsls(j)
              stemdp(i,j) = zstemdp(i,j)
            ENDDO
          ENDDO
          DO j=2,nft
            DO i=1,ftrls(j)
              rootdp(i,j) = zrootdp(i,j)
            ENDDO
          ENDDO

          DO i=1,nft
            lai(i) = zlai(i)
            npp(i) = znpp(i)
            nps(i) = znps(i)
            npl(i) = znpl(i)
            evp(i) = zevp(i)
            dof(i) = zdof(i)
            slc(i) = zslc(i)
            rlc(i) = zrlc(i)
            sln(i) = zsln(i)
            rln(i) = zrln(i)
            nppstore(i) = znppstore(i)
            nppstorx(i) = znppstorx(i)
            nppstor2(i) = znppstor2(i)
            bb(i) = zbb(i)
            bbgs(i) = zbbgs(i)
            dsbb(i) = zdsbb(i)
            maxlai(i) = zmaxlai(i)
          ENDDO

          DO i=1,200
            tmem = ztmem(i)
          ENDDO
          DO j=1,nft
            DO i=1,30
              sm_trig(i,j) = zsm_trig(i,j)
            ENDDO
          ENDDO
          DO ft=3,nft
            DO i=1,360
              sumadp(i,ft) = zsumadp(i,ft)
            ENDDO
          ENDDO
          DO ft=1,nft
            stemfr(ft) = zstemfr(ft)
          ENDDO
      
        ENDIF
      

*----------------------------------------------------------------------*
* Ensure cover array sums to 1.                                        *
*----------------------------------------------------------------------*
        bio(1,1,1)=0.0
        bio(1,2,1)=0.0
        hgt(1,1)=0.0
        wdt(1,1)=0.0
        DO ft=2,nft
          leaflit(ft) = 0.0d0
          yield(ft) = 0.0d0
          stemlit(ft) = 0.0d0
          rootlit(ft) = 0.0d0
          ftcov(ft) = 0.0d0
          DO age=1,ftmor(ft)
*            IF (ppm(age,ft).LT.1.0e-16) then
*              cov(age,ft)=0.0
*              bio(age,1,ft)=0.0
*              bio(age,2,ft)=0.0
*              hgt(age,ft)=0.0
*              wdt(age,ft)=0.0
*            ENDIF
            ftcov(ft) = ftcov(ft) + cov(age,ft)
          ENDDO
        ENDDO

        sum = 0.0d0
        DO ft=2,nft
          sum = sum + ftcov(ft)
        ENDDO

        IF (sum.LT.1.0d0) THEN
          cov(1,1) = 1.0d0 - sum
          ftcov(1) = cov(1,1)
        ELSE
          cov(1,1) = 0.0d0
          ftcov(1) = 0.0d0
          DO ft=2,nft
            ftcov(ft) = ftcov(ft)/sum
            DO age=1,ftmor(ft)
              cov(age,ft)=cov(age,ft)/sum
            ENDDO
          ENDDO
        ENDIF
*----------------------------------------------------------------------*

        biotoo = 0.0d0
        DO ft=1,nft
          DO age=1,ftmor(ft)
            biotoo = biotoo + (bio(age,1,ft) + bio(age,2,ft) 
     &)*cov(age,ft) 
          ENDDO
        ENDDO
        solcoo = 0.0d0
        DO i=1,8
          solcoo = solcoo + ic0(i)
        ENDDO

        c3old = cov(1,2)
        c4old = cov(1,3)

        snp_year = 1

!        soilt = 10.0d0

*----------------------------------------------------------------------*
*                               Year Loop                              *
*----------------------------------------------------------------------*
      DO iyear=1,nyears

        year    = yearv(iyear)
        metyear = met_yearv(iyear)
        IF(debug .EQV. .TRUE.) THEN
          PRINT*, 'YEAR=',iyear,year,metyear
          write(*,FMT="(A)") 'UU1'
          WRITE(*,*) 'BARE       CITY      C3p        C4p        ',
     &    'C3crop     ','C4crop     C3s       C4s        Ev_Bp      ',
     &    'Ev_Np      ','Dc_Bp      Dc_Np     Ev_Bs      Ev_Ns      ',
     &    'Dc_Bs      ','Dc_Ns      '

          write(*,FMT="(16F11.7)") cluse(1:nft,iyear)
          !write(*,FMT="(16F7.3)") sum_cov_test(1:nft)
          write(*,FMT="(A)") 'UU2'
          write(*,FMT="(8A14)")
     &      '     ','     primf', '     primn', '     secdf',
     &              '     secdn', '    c3crop', '    c4crop',
     &              '     urban', '      harv' 

          write(*,FMT="(A,8E14.7)") ' primf',cluse2(1,:,iyear),
     &           cluseh(1,iyear)
          write(*,FMT="(A,8E14.7)") ' primn',cluse2(2,:,iyear),
     &           cluseh(2,iyear)
          write(*,FMT="(A,8E14.7)") ' secdf',cluse2(3,:,iyear),
     &           cluseh(3,iyear)
          write(*,FMT="(A,8E14.7)") ' secdn',cluse2(4,:,iyear),
     &           cluseh(4,iyear)
          write(*,FMT="(A,8E14.7)") 'c3crop',cluse2(5,:,iyear),
     &           cluseh(5,iyear)
          write(*,FMT="(A,8E14.7)") 'c4crop',cluse2(6,:,iyear),
     &           cluseh(6,iyear)
          write(*,FMT="(A,8E14.7)") ' urban',cluse2(7,:,iyear),
     &           cluseh(7,iyear)
        ENDIF

        !APW - not sure what this does
        DO ft=1,nft
          laimax(ft) = 4.601d0 
        ENDDO

        !initialise daily output array with 0
        daily_out(:,:,:,:) = 0.0d0

        nfix = infix

*----------------------------------------------------------------------*
* Set CO2.                                                             *
*----------------------------------------------------------------------*

        ! if there is a spin up but iyear is in the run proper phase
        ! set CO2 to the year 
        iyear_adj = 0
        IF ((spinl.gt.0).AND.(iyear.GT.spinl)) THEN
          speedc = .FALSE.
          if(co2const.gt.0) iyear_adj = spinl
          !ca(:,:) = co2(year-yr0+1,:,:) !
          ca(:,:) = co2(iyear-iyear_adj,:,:) !
        ELSE
          IF (co2const.GT.0.0d0) THEN
            ca(:,:) = co2const
          ELSE
            ca(:,:) = co2(iyear,:,:) !
          ENDIF
        ENDIF

        !print*, year, metyear, iyear, ca(1,1)

        IF (mod(iyear,max(year_out,1)).EQ.min(1,year_out)-1) THEN
           WRITE(*,'('' Year no.'',2i5,'', ca = '',2f6.2)') iyear, !
     & year,ca(1,1),ca(12,31) !
        ENDIF

*----------------------------------------------------------------------*
* Set 'tmp' 'hum' 'prc', calculate:  annual stats 
*----------------------------------------------------------------------*
        DO mnth=1,12
          DO day=1,no_days(year,mnth,thty_dys)
            tmp(mnth,day) = 
     &real(xtmpv(metyear-yr0m+1,mnth,day))/100.0d0
            prc(mnth,day) = 
     &real(xprcv(metyear-yr0m+1,mnth,day))/10.0d0
            hum(mnth,day) = 
     &real(xhumv(metyear-yr0m+1,mnth,day))/100.0d0
            swr(mnth,day) = 
     &xswrv(metyear-yr0m+1,mnth,day)       
            IF (withcloudcover) THEN
               cld(mnth) = 
     &real(xcldv(metyear-yr0m+1,mnth))/1000.0d0
            ELSE
               cld(mnth) = 0.5d0
               if(read_par.eq.1) cld(mnth) = 1.0d0
               !print*, cld(1)
            ENDIF
          ENDDO
        ENDDO
C PCM2        WRITE(*,*)'metyear,yr0m' 
C PCM2        WRITE(*,*) metyear,yr0m 
C PCM2        DO mnth=1,12
C PCM2          WRITE(*,'(30f5.0)') xtmpv(metyear-yr0m+1,mnth,1:30) 
C PCM2        ENDDO
C PCM2        DO mnth=1,12
C PCM2          WRITE(*,'(30f5.1)') tmp(mnth,1:30) 
C PCM2        ENDDO

        !tmin       = 100.0d0
        yeartmp    = 0.0d0
        yearprc    = 0.0d0
        yearhum    = 0.0d0
        yearswr    = 0.0d0
        prc_week(:)= 0.0d0
        prcq(:)    = 0.0d0
        !yearp_days = 0.0d0
        !yearsoilc  = 0.0d0
        w = 1
        d = 1
        DO mnth=1,12
          mnthprc(mnth) = 0.0d0
          mnthtmp(mnth) = 0.0d0
          mnthhum(mnth) = 0.0d0
          mnthswr(mnth) = 0.0d0
          mnthp_days(mnth) = 0.0d0 
          ! calculate monthly values
          DO day=1,no_days(year,mnth,thty_dys)
            ! restrict humidity - presumably for stability
            IF (hum(mnth,day).LT.30.0d0)  hum(mnth,day) = 30.0d0
            IF (hum(mnth,day).GT.95.0d0)  hum(mnth,day) = 95.0d0
            mnthprc(mnth) = mnthprc(mnth) + 
     &prc(mnth,day)
!            IF(prc(mnth,day).gt.0.10d0) mnthp_days(mnth) = 
!     &mnthp_days(mnth) + 1
            mnthtmp(mnth) = mnthtmp(mnth) + 
     &tmp(mnth,day)/no_days(year,mnth,thty_dys)
            mnthhum(mnth) = mnthhum(mnth) + 
     &hum(mnth,day)/no_days(year,mnth,thty_dys)
            mnthswr(mnth) = mnthswr(mnth) + 
     &swr(mnth,day)/no_days(year,mnth,thty_dys)
            ! calculate weekly precip
            if(mod(d,7).eq.0) w = w + 1 
            d = d + 1
            if(w.ge.53) w = 52 
            prc_week(w)   = prc_week(w) + prc(mnth,day) 
          ENDDO ! day
          !IF (mnthtmp(mnth).LT.tmin)  tmin = mnthtmp(mnth)
          yeartmp = yeartmp + mnthtmp(mnth)/12.0d0
          yearprc = yearprc + mnthprc(mnth)
          yearhum = yearhum + mnthhum(mnth)/12.0d0
          yearswr = yearswr + mnthswr(mnth)/12.0d0
          !yearp_days  = yearp_days + mnthp_days(mnth)/12.0d0
        ENDDO ! month
        ! not sure what this does
        !tmin = tmin*1.29772d0 - 19.5362d0

        yeartmp_max = maxval(mnthtmp(:))
        yeartmp_min = minval(mnthtmp(:))
        ! calculate precip of the driest quarter
        do w=1,52
          s = w - 6 
          f = w + 6
          do w1=s,f
            if(w1.lt.1) then
              prcq(w) = prcq(w) + prc_week(52+w1)
            elseif(w1.gt.52) then
              prcq(w) = prcq(w) + prc_week(mod(w1,52))
            else 
              prcq(w) = prcq(w) + prc_week(w1)
            endif
          enddo
        enddo
        yearprcdryq = minval(prcq(:)) 
        ! use end of year value for soil C for trait regressions 
!        DO i=1,nft
!          yearsoilc = yearsoilc + slc(i) + rlc(i)
!        ENDDO
!        yearsoilc   = ic0(1) + ic0(2) + ic0(3) + ic0(4) + ic0(5) +
!     &ic0(6) + ic0(7) + ic0(8)
!        !assume top 30 cm contains 50% of total soil C and convert to kg
!        yearsoilc   = 0.50d0 * yearsoilc/1000.0d0

          IF(soilp_map.eq.1) THEN
            IF(iyear.eq.1) THEN
              soilpr = soilp_init
            ELSE
              soilpr = soilp
            ENDIF

            !replace missing values with mean
            if(soilpr.lt.0.0) soilpr = 569.2
            !restrict soil P to bounds used in regression
            if(soilpr.lt.21) then 
              soilpr = 21
            elseif(soilpr.gt.832) then 
              soilpr = 832
            endif
            tleaf_p  = 10**(0.352*log10(soilpr)-1.758)   
          ELSE
            tleaf_p   = 0.0d0
          ENDIF

* calculate annual mean climate variables for trait relationships, Ordonez, TERRABITES etc
        !IF((vcmax_type.ge.2).OR.(ncalc_type.ge.2)) THEN
          !initialise arrays
          IF(iyear.eq.1) THEN

            ! for TERRABITES need to calculate:
            ! mean an temp (oC)
            ! mean temp of max temp month (oC)
            ! mean temp of min temp month (oC)
            ! tvar tmax - tmin (oC)
            ! mean an prc (mm)
            ! mean an rh (%)
            ! mean annual SW radiation (wm2)
            ! driest quarter of the year (mm), to the nearest week 

            ! these initialisations are only necessary if
            ! averaging periods >1 year are used for the met data  
!            matmpv(:)     = matmp_init
!            matmp_maxv(:) = matmp_max_init
!            matmp_minv(:) = matmp_min_init
!            mapv(:)       = maprc_init
!            mahumv(:)     = mahum_init
!            maswrv(:)     = maswr_init
!            aprc_dryqv(:) = aprc_dryq_init   
!            map_daysv(:)  = map_days_init
!            masoilcv(:)   = masoilc_init
!            masoilwv(:)   = masoilw_init
            matmpv(:)     = yeartmp
            matmp_maxv(:) = yeartmp_max
            matmp_minv(:) = yeartmp_min
            mapv(:)       = yearprc
            mahumv(:)     = yearhum
            maswrv(:)     = yearswr
            aprc_dryqv(:) = yearprcdryq
 

            IF(soilcn_map.eq.1) THEN
              masoilcnv(:) = soilcn_init
            ELSE
              masoilcnv(:) = soilCtoN
            ENDIF

          ELSE

            !age arrays
            DO i=0,8 
              matmpv(10-i)     = matmpv(10-1-i)
              matmp_maxv(10-i) = matmp_maxv(10-1-i)
              matmp_minv(10-i) = matmp_minv(10-1-i)
              mapv(10-i)       = mapv(10-1-i) 
              mahumv(10-i)     = mahumv(10-1-i)
              maswrv(10-i)     = maswrv(10-1-i)
              aprc_dryqv(10-i) = aprc_dryqv(10-1-i)
            !  map_daysv(10-i)  = map_daysv(10-1-i)
            !  masoilcv(10-i)   = masoilcv(10-1-i)
            !  masoilwv(10-i)   = masoilwv(10-1-i)
            ENDDO
            DO i=1,2
              masoilcnv(i+1) = masoilcnv(i)
            ENDDO
          
            !add this years value to the vectori
            matmpv(1)     = yeartmp
            matmp_maxv(1) = yeartmp_max
            matmp_minv(1) = yeartmp_min
            mapv(1)       = yearprc
            mahumv(1)     = yearhum
            maswrv(1)     = yearswr
            aprc_dryqv(1) = yearprcdryq
            ! for the TERRABITES simulations only the values in the current year are used
            ! i.e. these values are not means over 10 years as suggested by the vectors 
            !matmpv(:)     = yeartmp
            !matmp_maxv(:) = yeartmp_max
            !matmp_minv(:) = yeartmp_min
            !mapv(:)       = yearprc
            !mahumv(:)     = yearhum
            !maswrv(:)     = yearswr
            !aprc_dryqv(:) = yearprcdryq
            !map_daysv(:)  = yearp_days
            !masoilcv(:)   = yearsoilc
            !masoilwv(:)   = yearsoilw

            masoilcnv(1)= soilcn

          ENDIF

          !calulate mean values                
          matmp     = 0.0
          matmp_max = 0.0
          matmp_min = 0.0
          maprc     = 0.0
          mahum     = 0.0
          maswr     = 0.0
          aprc_dryq = 0.0
          masoilcn  = 0.0
          do i=1,10
            matmp     = matmp + matmpv(i)    
            matmp_max = matmp_max + matmp_maxv(i)
            matmp_min = matmp_min + matmp_minv(i)
            maprc     = maprc + mapv(i)
            mahum     = mahum + mahumv(i)
            maswr     = maswr + maswrv(i)
            aprc_dryq = aprc_dryq + aprc_dryqv(i)
            !map_days  = map_days + map_daysv(i)
            !masoilc   = masoilc + masoilcv(i)
            !masoilw   = masoilw + masoilwv(i)
          enddo
          do i=1,3
            masoilcn = masoilcn + masoilcnv(i)
          enddo 

          matmp     = matmp/size(matmpv)    
          matmp_max = matmp_max/size(matmp_maxv)
          matmp_min = matmp_min/size(matmp_minv)
          maprc     = maprc/size(mapv)
          mahum     = mahum/size(mahumv)
          maswr     = maswr/size(maswrv)
          aprc_dryq = aprc_dryq/size(aprc_dryqv)
          !map_days  = map_days/size(map_daysv)
          !masoilc   = masoilc/size(masoilcv)
          !masoilw   = masoilw/size(masoilwv)

          masoilcn = masoilcn/size(masoilcnv)

          matvar   = matmp_max - matmp_min
          aprc_rel = aprc_dryq / maprc

          
*          print*, fttags(:)
*          print*, matmp     
*          print*, matmp_max 
*          print*, matmp_min 
*          print*, maprc     
*          print*, map_days  
*          print*, mahum    
*          print*, ' '

      !trait regressions
          !from correspondence with Peter van Bodegom 
          IF(ncalc_type.eq.2) THEN
            !fix climate variables if they are outside range of values used
            ! - in Ordonez trait~climate regressions
            ! write(*,*) maprc, maswr
            ! write(*,*) soilp, soilpr
            maprcr   = maprc
            maswrr   = maswr
            masoilcnr= masoilcn

            if(maprc.lt.300)  maprcr = 300
            if(maswr.gt.192)  maswrr = 192
            if(masoilcn.lt.6.6) then 
              masoilcnr = 6.6
            elseif(masoilcn.gt.55) then 
              masoilcnr = 55
            endif
            !topleaf N
            tleaf_n   = 10**(-14.0319+3.9991*log10(maprcr)+7.9423*log10
     &(masoilcnr)+0.0512*maswrr-1.8651*log10(maprcr)*log10(masoilcnr)
     &-0.0129*log10(maprcr)*maswrr-0.0135*log10(masoilcnr)*maswrr)
          !SLA in dry weight units 
c          tleaf_sla = 10**(7.6236-1.9505*log10(maprcr)-5.4218
c     &*log10(masoilcnr)+1.5910*log10(maprcr)*log10(masoilcnr))

          !convert m2/kg to m2/g
c          tleaf_sla = tleaf_sla/1000  
          ENDIF

          ! Calculate Vcmax, Jmax and SLA for the TERRABITES sims
          env_vcmax(:) = 0.0d0
          env_jmax(:)  = 0.0d0
          
          IF(vcmax_type.ne.3) THEN
          !van Bodegom and Verheijin trait regressions

          !C3 grass/forb
          env_vcmax(3)   = 195.98d0 - 305.25d0*aprc_rel - 
     &5.31d0*matmp_min
          env_jmax(3)    = - 269.21d0 - 0.197d0*maprc - 7.06d0*matvar
          ftsla(3)       = 20.44 + 0.59d0*matmp - 0.018*aprc_dryq +
     &30.76*aprc_rel - 0.089*maswr + 0.00043*maprc*matvar
 
          !C4 grass/forb
          env_vcmax(4)   = 61.02d0 - 0.549d0*mahum
          env_jmax(4)    = 688.03d0 - 8.32d0*mahum
          ftsla(4)       = - 66.94d0 + 1.19d0*mahum + 4.53d0*matmp -
     &5.10*matmp_min + 0.051*mahum*matvar + 0.16*mahum*matmp_min -
     &0.158*mahum*matmp
 
          !Ev_Bl
          env_vcmax(7)   = 54.90d0 + 2.76d0*matvar - 335.63d0*aprc_rel
          env_jmax(7)    = 148.25d0 - 2.91d0*matmp_min
          ftsla(7)    = 14.45d0 + 0.272d0*matmp + 0.0023*maprc -
     &0.0072*aprc_dryq - 0.063*maswr
 
          !Ev_Nl
          env_vcmax(8)   = -199.41d0 - 6.49d0*mahum - 371.52d0*aprc_rel 
     &+ 79.51d0*matmp_max + 48.59d0*matvar - 1.74d0*maswr - 
     &3.53*matmp_max*matvar
          env_jmax(8)    = 795.67d0 - 0.22d0*maprc + 19.09d0*matmp_min -
     &0.78*aprc_dryq - 1.67*maswr
          ftsla(8)       = -10.37 + 0.075d0*mahum - 28.57d0*aprc_rel +
     &0.31d0*matvar - 1.05d0*matmp_max + 0.96d0*matmp - 0.0047*maprc +
     &0.021*aprc_dryq
 
          !Dc_Bl
          env_vcmax(9)   = - 31.55d0 + 287.41d0*aprc_rel + 3.51d0*matvar
     &- 30.11d0*aprc_rel*matvar
          env_jmax(9)    = 74.22d0 + 2.05d0*matvar 
          ftsla(9)       = 28.73d0 + 3.11d0*matmp - 1.62d0*matmp_max -
     &1.54*matmp_min - 0.067d0*maswr
 
          !Dc_Nl - vcmax and jmax same as Ev_Nl
          env_vcmax(10)  = env_vcmax(8)
          env_jmax(10)   = env_jmax(8)
          ftsla(10)      = -50.60d0 + 0.76d0*mahum + 409.75*aprc_rel -
     &5.20d0*mahum*aprc_rel
 
          ELSEIF(vcmax_type.eq.3) THEN
          ! Calculate Vcmax, Jmax and SLA for the TERRABITES sims
          ! - mean values
          
          !C3 grass/forb
          env_vcmax(3)   = 62.220d0
          env_jmax(3)    = 121.25d0
          ftsla(3)       = 20.15d0
 
          !C4 grass/forb
          env_vcmax(4)   = 27.270d0
          env_jmax(4)    = 176.44d0
          ftsla(4)       = 18.50d0
          
          !Ev_Bl
          env_vcmax(7)   = 35.5d0 
          env_jmax(7)    = 84.2d0
          ftsla(7)       = 9.350d0 
 
          !Ev_Nl
          env_vcmax(8)   = 71.0d0
          env_jmax(8)    = 148.2d0
          ftsla(8)       = 4.520d0
 
          !Dc_Bl
          env_vcmax(9)   = 55.0d0
          env_jmax(9)    = 92.870d0
          ftsla(9)       = 14.48d0
 
          !Dc_Nl
          env_vcmax(10)  = env_vcmax(8)
          env_jmax(10)   = env_jmax(8)
          ftsla(10)      = 9.46d0
          ENDIF
       
c        ELSE
c          !satisfy initialisation output when trait regressions are not used 
c          matmp     = matmp_init 
c          matmp_max = matmp_max_init
c          matmp_min = matmp_min_init
c          maprc     = maprc_init
c          mahum     = mahum_init
c          maswr     = maswr_init
c          !map_days  = map_days_init
c          !masoilc   = masoilc_init
c          !masoilw   = masoilw_init
c        
c          masoilcn = soilcn_init
c          soilp    = soilp_init
c        !leaf trait loop  
c        ENDIF         


        IF(vcmax_type.ne.3) THEN
          !restrict to 95 %iles of observed trait data
          
          !print*, env_vcmax(:)
          !print*, env_jmax(:)
          
          !bounds for TRY traits
          env_vcmax_min(1:nft) = (/ 0.0,0.0,25.34,21.72,25.34,21.72,
     &16.3,19.9,19.43,19.9 /)
          env_vcmax_max(1:nft) = (/ 0.0,0.0,116.08,46.33,116.08,46.33,
     &93.99,178.0,130.0,180.0 /)
          env_jmax_min(1:nft)  = (/ 0.0,0.0,44.87,92.3,44.87,92.3,
     &35.8,57.43,41.9,57.43 /)
          env_jmax_max(1:nft)  = (/ 0.0,0.0,230.79,465.4,230.79,465.4,
     &164.7,333.2,206.5,333.2 /)
          env_sla_min(1:nft)   = (/ 0.0,0.0,5.42,5.71,5.42,5.71,
     &4.37,2.4,7.51,4.52 /)
          env_sla_max(1:nft)   = (/ 0.0,0.0,43.34,34.34,43.34,34.34,
     &21.3,9.84,39.85,13.93 /) 

          ! trait regression parameters
          jmax_int(1:nft)      = (/ 0.0,0.0,35.46,-236.96,35.46,-236.96,
     &29.34,30.62,25.54,30.62 /)
          jmax_int_er(1:nft)   = (/ 0.0,0.0,83.95,0.0062,83.95,0.0062,
     &16.03,35.61,15.55,35.61 /)
          jmax_slope(1:nft)     = (/ 0.0,0.0,1.51,15.15,1.51,15.15,
     &0.96,1.88,1.32,1.88 /)
          jmax_slope_er(1:nft)  = (/ 0.0,0.0,1.28,0.0002,1.28,0.0002,
     &0.23,0.467,0.27,0.467 /)
          sla_int(1:nft)        = (/ 0.0,0.0,0.0,45.26,0.0,45.26,
     &0.0,0.0,25.7,0.0 /)
          sla_int_er(1:nft)     = (/ 0.0,0.0,0.0,9.58,0.0,9.58,
     &0.0,0.0,3.84,0.0 /)
          sla_slope(1:nft)      = (/ 0.0,0.0,0.0,-0.97,0.0,-0.97,
     &0.0,0.0,-0.22,0.0 /)
          sla_slope_er(1:nft)   = (/ 0.0,0.0,0.0,0.30,0.0,0.30,
     &0.0,0.0,0.076,0.0 /)

          !zero bounds flag arrays
          env_vcmax_bounds(:) = 0
          env_jmax_bounds(:)  = 0
          env_sla_bounds(:)   = 0
          
          DO ft=1,nft
            if(env_vcmax(ft).lt.env_vcmax_min(ft)) then
              env_vcmax(ft) = env_vcmax_min(ft) 
              env_vcmax_bounds(ft) = 1
            endif
            if(env_vcmax(ft).gt.env_vcmax_max(ft)) then
              env_vcmax(ft) = env_vcmax_max(ft) 
              env_vcmax_bounds(ft) = 1
            endif
            if(env_jmax(ft).lt.env_jmax_min(ft)) then
              env_jmax(ft) = env_jmax_min(ft) 
              env_jmax_bounds(ft) = 1
            endif
            if(env_jmax(ft).gt.env_jmax_max(ft)) then
              env_jmax(ft) = env_jmax_max(ft) 
              env_jmax_bounds(ft) = 1
            endif
            if(ftsla(ft).lt.env_sla_min(ft)) then
              ftsla(ft) = env_sla_min(ft) 
              env_sla_bounds(ft) = 1
            endif
            if(ftsla(ft).gt.env_sla_max(ft)) then
              ftsla(ft) = env_sla_max(ft) 
              env_sla_bounds(ft) = 1
            endif

            ! check that values of vcmax, jmax and sla do not violate observed covariance between these variables. 
            !   Regression CI95 lower range:
            jmax_ci_low  = 
     &jmax_int(ft) - jmax_int_er(ft) + 
     &((jmax_slope(ft)-jmax_slope_er(ft)) * env_vcmax(ft))
            !  Regression CI95 upper range:
            jmax_ci_high  = 
     &jmax_int(ft) + jmax_int_er(ft) + 
     &((jmax_slope(ft)+jmax_slope_er(ft)) * env_vcmax(ft))

            if(env_jmax(ft).lt.jmax_ci_low) then
              !first determine slope and intercept for the perpendicular line:
              ! Jmax=a2*Vmax+b2
              a2 = (-1/(jmax_slope(ft)-jmax_slope_er(ft))) !slope perpendiculair line
              b2 = env_jmax(ft)-(a2*env_vcmax(ft))         ! intercept perpendiculair line
              !search values when slopes cross:
              !(a2*Vmax2)+b2=(a*Vmax)+intercept_lower
              !(a2-slope_lower)*Vmax=intercept_lower-b2
              env_vcmax(ft) = (jmax_int(ft)-jmax_int_er(ft) - b2) /
     &(a2 - jmax_slope(ft)-jmax_slope_er(ft)) 
              env_jmax(ft)  = jmax_int(ft)-jmax_int_er(ft) +
     &((jmax_slope(ft)-jmax_slope_er(ft))*env_vcmax(ft))

              !ensure that the interpolated values still match the observed ranges
              if(env_vcmax(ft).lt.env_vcmax_min(ft)) then
                env_vcmax(ft) = env_vcmax_min(ft) 
              endif
              if(env_jmax(ft).lt.env_jmax_min(ft)) then
                env_jmax(ft) = env_jmax_min(ft) 
              endif

            elseif(env_jmax(ft).gt.jmax_ci_high) then
              a2 = (-1/(jmax_slope(ft)+jmax_slope_er(ft))) !slope perpendiculair line
              b2 = env_jmax(ft)-(a2*env_vcmax(ft))         ! intercept perpendiculair line
              env_vcmax(ft) = (jmax_int(ft)+jmax_int_er(ft) - b2) /
     &(a2 - jmax_slope(ft)+jmax_slope_er(ft)) 
              env_jmax(ft)  = jmax_int(ft)+jmax_int_er(ft) +
     &((jmax_slope(ft)+jmax_slope_er(ft))*env_vcmax(ft))

              !ensure that the interpolated values still match the observed ranges
              if(env_vcmax(ft).gt.env_vcmax_max(ft)) env_vcmax(ft) = 
     &env_vcmax_max(ft) 
              if(env_jmax(ft).gt.env_jmax_max(ft))   env_jmax(ft) = 
     &env_jmax_max(ft) 
            endif

            !evaluate whether Vmax matches the trade-off with SLA and adjust SLA if necessary
            if((ft.eq.4).or.(ft.eq.9))then
              !Regression CI95 lower range:
              sla_ci_low  = sla_int(ft)-sla_int_er(ft) +
     &((sla_slope(ft)-sla_slope_er(ft))*env_vcmax(ft))
              !Regression CI95 upper range:
              sla_ci_high = sla_int(ft)+sla_int_er(ft) +
     &((sla_slope(ft)+sla_slope_er(ft))*env_vcmax(ft))

              if(ftsla(ft).lt.sla_ci_low) then
                ftsla(ft) = sla_ci_low
                if(ftsla(ft).lt.env_sla_min(ft)) ftsla(ft) = 
     &env_sla_min(ft) 
              elseif(ftsla(ft).lt.sla_ci_high) then
                ftsla(ft) = sla_ci_high
                if(ftsla(ft).gt.env_sla_max(ft)) ftsla(ft) = 
     &env_sla_max(ft) 
              endif
            endif
          enddo
        ENDIF

        IF((vcmax_type.ge.2).and.(vcmax_type.le.3)) THEN
          ! C3/C4 crop values same as C3/C4 grass
          env_vcmax(5) = env_vcmax(3)
          env_jmax(5)  = env_jmax(3)
          ftsla(5)     = ftsla(3)
          env_vcmax(6) = env_vcmax(4)
          env_jmax(6)  = env_jmax(4)
          ftsla(6)     = ftsla(4)

          ! correct for CO2 acclimation
          eco2 = 0.0d0
          DO mnth=1,12
            DO day=1,no_days(year,mnth,thty_dys)
              eco2 = eco2 + ca(mnth,day)/0.10130d0
            ENDDO
          ENDDO
          IF(thty_dys.eq.1) THEN
            eco2 = eco2 / 360.0 
          ELSE
            eco2 = eco2 / 365.0  
          ENDIF
      
          !trees
          env_vcmax(7:10) = (-2.570d-4*eco2+1.090d0)/0.9990d0*
     &env_vcmax(7:10)
          env_jmax(7:10)  = (-3.430d-4*eco2+1.120d0)/0.99860d0*
     &env_jmax(7:10)
        
          !C3 grass/crop
          env_vcmax(3) = (-7.570d-4*eco2+1.274d0)/1.0060d0*
     &env_vcmax(3)
          env_jmax(3)  = (-3.80d-4*eco2+1.138d0)/1.00350d0*
     &env_jmax(3)
          env_vcmax(5) = (-7.570d-4*eco2+1.274d0)/1.0060d0*
     &env_vcmax(5)
          env_jmax(5)  = (-3.80d-4*eco2+1.138d0)/1.12250d0*
     &env_jmax(5)
       
          ftsla(:)     = ftsla(:) /
     &( exp(-0.2510d0+8.0d-4*eco2 - 2.90d-5*(eco2/10)**2) / 
     &  exp(-0.2510d0+8.0d-4*35.39- 2.90d-5*35.39**2 )   )  

          ftsla(4)     = ftsla(4) /
     &( exp(-0.06510d0+3.380d-4*eco2 - 3.90d-5*(eco2/10)**2) / 
     &  exp(-0.06510d0+3.380d-4*35.39- 3.90d-5*35.39**2 )   )   

          !print*, env_vcmax(:)
          !print*, env_jmax(:)
          !convert mm2/mg to m2/g
          ftsla(:)  =  ftsla(:) / 1000.0d0

          !print*, env_vcmax(:)
          !print*, env_jmax(:)
        ENDIF

        ! extract vcmax & pepc values for C4
        trait_env_C4_vcmax = env_vcmax(4)
        trait_env_C4_pepc  = env_jmax(4)

        !overwrite trait environment relationship of SLA when not required 
        IF ((vcmax_type.ne.2).and.(vcmax_type.ne.3)) THEN
          ftsla(:) = input_ftsla(:)
        ENDIF



*----------------------------------------------------------------------*
* Carbon at the start of the year.                                     *
*----------------------------------------------------------------------*
        DO ft=1,nft
          bioleaf(ft) = 0.0d0
          if(s070607.eq.1) then
            DO day=1,ftlls(ft)
              bioleaf(ft) = bioleaf(ft) + 
     &leafdp(day,ft)*12.0d0/ftsla(ft)/18.0d0
            ENDDO
          else
            DO day=1,ftlls(ft)
              bioleaf(ft) = bioleaf(ft) + 
     &leafdp(day,ft)/(ftsla(ft)/0.480d0)
            ENDDO
          endif
        ENDDO

        ans1 = 0.0d0
        DO ft=1,nft
          DO i=1,ftmor(ft)
            ! adding stem_carbon=bio(i,1,ft) & root_carbon(i,2,ft) & leaf_carbon & nppstore
            ans1 = ans1 + (bio(i,1,ft) + bio(i,2,ft) + bioleaf(ft) +
     &nppstore(ft))*cov(i,ft)
          ENDDO
          ans1 = ans1 + slc(ft) + rlc(ft) ! adding stem_litter_carbon(ft) & root_litter_carbon(ft)
        ENDDO
        IF(debug .EQV. .TRUE.) THEN
          WRITE(*,'(''Icheck0'',3f13.6)') ccheck
          WRITE(*,'(''Ians1'',3f12.6)') ans1 
          WRITE(*,'(''Itc0(1)'',3f12.6)') ic0(1) 
          WRITE(*,'(''Itc0(2)'',3f12.6)') ic0(2) 
          WRITE(*,'(''Itc0(3)'',3f12.6)') ic0(3) 
          WRITE(*,'(''Itc0(4)'',3f12.6)') ic0(4) 
          WRITE(*,'(''Itc0(5)'',3f12.6)') ic0(5) 
          WRITE(*,'(''Itc0(6)'',3f12.6)') ic0(6) 
          WRITE(*,'(''Itc0(7)'',3f12.6)') ic0(7) 
          WRITE(*,'(''Itc0(8)'',3f12.6)') ic0(8) 
        ENDIF

        ccheck = ans1 + ic0(1) + 
     &ic0(2) + ic0(3) + ic0(4) + ic0(5) + ic0(6) + ic0(7) + ic0(8)

        IF(debug .EQV. .TRUE.) THEN
          WRITE(*,'(''Icheck1'',3f13.6)') ccheck

          PRINT '(A)','SS3a ftprop (from states;from netTransitions) 
     &(previous time step)'
          PRINT '(16F11.7)',ftprop(1:nft)
        ENDIF
*----------------------------------------------------------------------*
* Set land use through ftprop.                                         *
*----------------------------------------------------------------------*
        IF (ilanduse.EQ.2) THEN
          CALL NATURAL_VEG(tmp,prc,ftprop,nat_map)
        ELSE ! 0 and 1 and 3 to 6
          ftprop(1) = 100.0d0
          DO ft=2,nft
            IF (check_ft_grow(tmp,ftbbm(ft),ftbb0(ft),ftbbmax(ft), 
     &ftbblim(ft),chill(ft),dschill(ft)).EQ.1) THEN
              !ftprop(ft) = cluse(ft,year-yr0+1)
              !if((co2const.gt.0.0).and.(spinl.gt.0).and.
            !&(iyear.le.spinl)) then
              !  ftprop(ft) = cluse(ft,1)
              !else
              !  ftprop(ft) = cluse(ft,iyear-iyear_adj)
              !endif  

              ! logic below is identical to the co2 logic
              if((spinl.gt.0).and.(iyear.gt.spinl)) then
                ftprop(ft) = cluse(ft,iyear-iyear_adj)
              else
C PCM temporarily changed the following line for S4v8-S6v8 TRENDY
                if(co2const.gt.0.0) then !For TRENDY S1-S3
C    if((co2const.gt.0.0).and.(spinl.lt.nyears)) then !For TRENDY S4-S6
                  ftprop(ft) = cluse(ft,1)
                else
            !print*, ft, iyear, iyear_adj, cluse(ft,iyear-iyear_adj)
                  ftprop(ft) = cluse(ft,iyear-iyear_adj)
                endif  
              endif  

              ftprop(1)  = ftprop(1) - ftprop(ft)
            ELSE
              ftprop(ft) = 0.0d0
            ENDIF

            IF (ilanduse.GE.3 .AND. ilanduse.LE.6) THEN
            at = aggmap_SDGVM_to_aggHyde(ft)
            DO at2=1,maxn_at
              ! logic below is identical to the co2 logic
              if((spinl.gt.0).and.(iyear.gt.spinl)) then
                atprop2(at,at2) = cluse2(at,at2,iyear-iyear_adj)
              else
C PCM temporarily changed the following line for S4v8-S6v8 TRENDY
                if(co2const.gt.0.0) then !For TRENDY S1-S3
C    if((co2const.gt.0.0).and.(spinl.lt.nyears)) then !For TRENDY S4-S6
                  atprop2(at,at2) = cluse2(at,at2,1)
                else
                  atprop2(at,at2) = cluse2(at,at2,iyear-iyear_adj)
                endif  
              endif  
            ENDDO
            ENDIF

            IF (ilanduse.GE.3 .AND. ilanduse.LE.6) THEN
            at = aggmap_SDGVM_to_aggHyde(ft)
              ! logic below is identical to the co2 logic
              if((spinl.gt.0).and.(iyear.gt.spinl)) then
                atharvest(at) = cluseh(at,iyear-iyear_adj)
              else
C PCM temporarily changed the following line for S4v8-S6v8 TRENDY
                if(co2const.gt.0.0) then !For TRENDY S1-S3
C    if((co2const.gt.0.0).and.(spinl.lt.nyears)) then !For TRENDY S4-S6
                  atharvest(at) = cluseh(at,1)
                else
                  atharvest(at) = cluseh(at,iyear-iyear_adj)
                endif  
              endif  
            ENDIF


          ENDDO
          IF (ftprop(1).LT.0.0d0) THEN
            DO ft=2,nft
              ftprop(ft) = ftprop(ft)*100.0d0/(100.0 - ftprop(1))
            ENDDO
            ftprop(1) = 0.0d0
          ENDIF
        ENDIF

        IF(debug .EQV. .TRUE.) THEN
          PRINT '(A)','SS3b ftprop (from states; from net transitions) '
          PRINT '(16F11.7)',ftprop(1:nft)
        ENDIF

        IF (closed_loop_ft .EQV. .FALSE.) THEN
          ftprop1 = ftprop !net transitions or 1st year of gross transitions 
        ENDIF

        IF(debug .EQV. .TRUE.) THEN
          PRINT '(A)','SS3c ftprop1 (before COVER( ) routine) '
          PRINT '(16F11.7)',ftprop1(1:nft)
        ENDIF
*----------------------------------------------------------------------*
        CALL COVER(nft,ftmor,ftppm0,cov,bio,bioleaf,nppstore,
     &npp,nps,mnthtmp,mnthprc,slc,rlc,c3old,c4old,firec,ppm,hgt,fireres,
     &fprob,ftprop1,ftstmx,stemdp,rootdp,ftsls,ftrls,ilanduse,nat_map,
     &ic0,fire(iyear),harvest(iyear),leafdp,flulccc,ftphen,atprop2,
     &atharvest,aggmap_SDGVM_to_aggHyde,debug)

        IF(debug .EQV. .TRUE.) THEN
          PRINT '(A)','SS3d ftprop1 (after  COVER( ) routine) '
          PRINT '(16F11.7)',ftprop1(1:nft)
        ENDIF

        CALL MKDLIT(nft,ftmor,ftcov,dslc,drlc,dsln,drln,cov,slc,rlc,sln,
     &rln)

        IF(debug .EQV. .TRUE.) THEN
          !PRINT *,'SS3',nppstore(1:nft)
        ENDIF

        IF (ilanduse.GE.3 .AND. ilanduse.LE.6) THEN !turn on after 1st year
           closed_loop_ft = .TRUE.
        ENDIF

        DO i=1,8
          tc0(i) = 0.0d0
          tn0(i) = 0.0d0
        ENDDO

        DO i=1,3
          tminn(i) = 0.0d0
        ENDDO
        ts1 = 0.0d0
        ts2 = 0.0d0
        ts3 = 0.0d0
        ts4 = 0.0d0
        tsn = 0.0d0
        tlsn = 0.0d0

*----------------------------------------------------------------------*
* Find which DOLY runs are required.                                   *
*----------------------------------------------------------------------*
        DO ft=1,nft
          dolydo(ft) = 0
          DO age=1,ftmor(ft)
            IF (cov(age,ft).GT.0.0d0) THEN
              dolydo(ft) = 1
              GOTO 129
            ENDIF
          ENDDO

129       CONTINUE
        ENDDO

*----------------------------------------------------------------------*
* Initialisations that were in doly at the beginning of the year       *
*----------------------------------------------------------------------*

        iminn(3) = iminn(1) + iminn(2)

        isoilc = 0.0d0
        isoiln = 0.0d0
        DO i=1,8
           isoilc = isoilc + ic0(i)
           isoiln = isoiln + in0(i)
        ENDDO
        isoiln = isoiln + iminn(3)

c initialisse for all the ft, even for those doly run is not required
        DO ft=1,nft
          leaflit(ft) = 0.0d0
          yield(ft) = 0.0d0
          stemlit(ft) = 0.0d0
          rootlit(ft) = 0.0d0

          leafnpp(ft) = 0.0d0
          stemnpp(ft) = 0.0d0
          rootnpp(ft) = 0.0d0

          dof(ft) = 0.0d0

          lai(ft) = 0.0d0
          DO day=1,ftlls(ft)
            lai(ft) = lai(ft) + leafdp(day,ft)
          ENDDO

c     water
*           s1(ft) = is1
*           s2(ft) = is2
*           s3(ft) = is3
*           s4(ft) = is4
*           sn(ft) = isn
*           lsn(ft) = ilsn
          DO i=1,8
            c0(i,ft) = ic0(i)
            n0(i,ft) = in0(i)
          ENDDO

          DO i=1,3
            minn(i,ft) = iminn(i)
          ENDDO
          soilc(ft) = isoilc
          soiln(ft) = isoiln

          evp(ft) = 0.0d0

          sresp(ft) = 0.0d0
          lch(ft) = 0.0d0

          ! calculate Amax as a function of N uptake
          ! then calc Vcmax from Amax (SDGVM original method)
          amax(ft) = 0.0d0
          IF(vcmax_type.eq.5) THEN
          !calculation of Vcmax from Amax and N uptake rate
          ! see Woodward 1994 Ecophys. of Photosyn. & ABR
          ! Amax (assumed topleaf) is calculated from N uptake rate which is a function of soil C & N, as described in the above references
          ! Vcmax is then calculated from Amax assuming carboxylation limited photosyn, VDP 0.75 kPa, ca 380 ppmv 
            !the below amax eq is from WS 1994 ABR fig 9
            amax(ft) = 50 * 0.999927**soilc(ft) * 
     &min(1.d0,0.001667*soiln(ft))
c            print*, amax(ft)
            env_vcmax(ft) = VCMAX_FROM_AMAX(amax(ft),gs_func,p_rd,
     &ftg0(ft),ftg1(ft),0)
c            print*, ft, env_vcmax(ft)
          ELSEIF(vcmax_type.eq.6) THEN
          !calculation of Vcmax from Amax and N uptake rate
          ! see Woodward et al 1995 GBC
          ! Amax (assumed topleaf) is calculated from N uptake rate which is a function of soil C & N and temperature, as described in the above reference
          ! Vcmax is then calculated from Amax assuming carboxylation limited photosyn, VDP 0.75 kPa, ca 380 ppmv 
          ! mean annual temprature used but could use growing season or other such temp 
          
            amax(ft) = f_amax(soilc(ft),soiln(ft),matmp)
            env_vcmax(ft) = VCMAX_FROM_AMAX(amax(ft),gs_func,p_rd,
     &ftg0(ft),ftg1(ft),0)
            !print*, ft,soilc(ft),soiln(ft),matmp,amax(ft),env_vcmax(ft)
          ENDIF

          !In case where Amax Vcmax method used, set C4 grass/forb vcmax and PEPC (jmax hijacked) using
          !trait environment relationships or mean
          env_vcmax(4)   = trait_env_C4_vcmax
          env_jmax(4)    = trait_env_C4_pepc
          ! set C4 crops the same
          env_vcmax(6) = env_vcmax(4)
          env_jmax(6)  = env_jmax(4)
          
          
          
* Height (m)
cccn           ht(ft) = 0.807d0*(laimax(ft)**2.13655d0)
          ht(ft) = 0.807d0*(5.0d0**2.13655d0)
          IF (ht(ft).GT.50.0d0)  ht(ft)=50.0d0
          laimax(ft)=0.0d0   
cThis is to do what doly did... but 
cit's strange. Initialisation to 0.0d0 is ok, but the previous 
c one laimax=4.6d0 above seems not good.

          

c     Convert nppstore to mols
          nppstore(ft) = nppstore(ft)/12.0d0
          nppstorx(ft) = nppstorx(ft)/12.0d0
          nppstor2(ft) = nppstor2(ft)/12.0d0
          nppstoreold(ft) = nppstore(ft)

          budo(ft) = 0
          seno(ft) = 0

        ENDDO ! ft's

c      print*,'start monthly loop'

*----------------------------------------------------------------------*
* START OF MONTH LOOP                                                  *
*----------------------------------------------------------------------*
      DO mnth=1,12

c     monthly initialisations
        DO ft=1,nft
           flow1(ft) = 0.0d0
           flow2(ft) = 0.0d0

           evapm(mnth,ft) = 0.0d0
           tranm(mnth,ft) = 0.0d0
           roffm(mnth,ft) = 0.0d0
           petm(mnth,ft) = 0.0d0

           photm(mnth,ft) = 0.0d0

           laimnth(mnth,ft) = 0.0d0
           avmnpet(ft) = 0.0d0

           !traceability analysis soil pool outputs 
           ! are flows consistent due to the mixing of resources, probably not
           daily_out(48,ft,mnth,:) = c0(1,ft)
           daily_out(49,ft,mnth,:) = c0(2,ft)
           daily_out(50,ft,mnth,:) = c0(3,ft)
           daily_out(51,ft,mnth,:) = c0(4,ft)
           daily_out(52,ft,mnth,:) = c0(5,ft)
           daily_out(53,ft,mnth,:) = c0(6,ft)
           daily_out(54,ft,mnth,:) = c0(7,ft)
           daily_out(55,ft,mnth,:) = c0(8,ft)
           daily_out(56,ft,mnth,1) = f0(1,ft)
           daily_out(57,ft,mnth,1) = f0(2,ft)
           daily_out(58,ft,mnth,1) = f0(3,ft)
           daily_out(59,ft,mnth,1) = f0(4,ft)
           daily_out(60,ft,mnth,1) = f0(5,ft)
           daily_out(61,ft,mnth,1) = f0(6,ft)
           daily_out(62,ft,mnth,1) = f0(7,ft)
           daily_out(63,ft,mnth,1) = f0(8,ft)
           daily_out(64,ft,mnth,1) = f0(9,ft)
           daily_out(65,ft,mnth,1) = f0(10,ft)
           daily_out(66,ft,mnth,1) = f0(11,ft)
           daily_out(67,ft,mnth,1) = f0(12,ft)
           daily_out(68,ft,mnth,1) = f0(13,ft)
           daily_out(69,ft,mnth,1) = f0(14,ft)
           daily_out(70,ft,mnth,1) = f0(12,ft) * cal

        ENDDO

        avmnppt = 0.0d0
        avmnt = 0.0d0
*----------------------------------------------------------------------*
* DAILY LOOP.                                                          *
*----------------------------------------------------------------------*
        DO day=1,no_days(year,mnth,thty_dys)
          fpr=0.0d0

          avmnppt = avmnppt + prc(mnth,day)
          avmnt = avmnt + tmp(mnth,day)/no_days(year,mnth,thty_dys)

*----------------------------------------------------------------------*
* Updata daily climate memory.                                         *
*----------------------------------------------------------------------*
          DO i=1,199
            tmem(201-i) = tmem(200-i)
          ENDDO
          tmem(1) = tmp(mnth,day)

*----------------------------------------------------------------------*
* Mix water resources.                                                 *
*----------------------------------------------------------------------*
          CALL MIX_WATER(s1,s2,s3,s4,sn,lsn,ftcov,ftmix,nft)
*----------------------------------------------------------------------*

          DO ft=1,nft
          IF (dolydo(ft).EQ.1) THEN

            DO d=1,ftlls(ft)
              leafv(d)=leafdp(d,ft)
            ENDDO
            DO d=1,ftsls(ft)
              stemv(d)=stemdp(d,ft)
            ENDDO
            DO d=1,ftrls(ft)
              rootv(d)=rootdp(d,ft)
            ENDDO

            DO d=1,30
              smtrig(d) = sm_trig(d,ft)
            ENDDO

            DO d=1,360
              suma(d) = sumadp(d,ft)
            ENDDO

*This is not really usefull since every ft's have the same minn 
*at the moment... but this may (should) be changed in the future
            DO i=1,3
              minnv(i)=minn(i,ft)
            ENDDO
            
            daygpp = 0.0d0
            dayra  = 0.0d0
            evap   = 0.0d0
            tran   = 0.0d0
            roff   = 0.0d0
            pet    = 0.0d0
            fpr    = 0.0d0

            leafold  = leafnpp(ft)
            stemold  = stemnpp(ft)
            rootold  = rootnpp(ft)
            nppsold  = nppstore(ft)
            leafresp = 0.0d0
            rootresp = 0.0d0
            stemresp = 0.0d0
            resp     = 0.0d0 
             
*----------------------------------------------------------------------*
* nppstore mols
* leafnpp  mols
* stemnpp  mols
* rootnpp  mols
* leaflit  mols
* stemlit  mols
* rootlit  mols
*----------------------------------------------------------------------*
            oldlai   = lai(ft)
            nleaf    = 0.0d0
            leaf_nit = 0.0d0
            !jmax     = 0.0d0
            !vcmax    = 0.0d0

            !print*, vcmax(1,ft)

            lflitold = leaflit(ft)
C PCM2      print*,'d2 ','tmp,prc,hum,cld,ft,soilc(ft),s1(ft),year,mnth,day'
C PCM2      print*,'dol ',tmp(mnth,day),prc(mnth,day),hum(mnth,day),cld(mnth), 
C PCM2     &ft,soilc(ft),s1(ft),year,mnth,day
!      stop
!            soilt = 0.97d0*soilt + 0.03d0*tmp(mnth,day)
!      write(*,*) mnth,day,tleaf_n
*            print*,'daily soilc 1:3',soilc(1:3)

            CALL DOLYDAY(ftsla(ft),ftc3(ft),ftphen(ft),ftagh(ft),
     &ftdth(ft),ftlls(ft),ftsls(ft),ftrls(ft),ftbbm(ft),
     &ftbb0(ft),ftbbmax(ft),ftbblim(ft),ftssm(ft),ftsss(ft),ftsslim(ft),
     &ftrat(ft),lat,dep,tmp(mnth,day),prc(mnth,day),hum(mnth,day),
     &cld(mnth),ca(mnth,day),soilc(ft),soiln(ft),minnv,s1(ft),s2(ft),
     &s3(ft),s4(ft),sn(ft),lsn(ft),adp,sfc,sw,sswc,awl,kd,kx,daygpp,
     &dayra,lai(ft),nppstore(ft),nppstorx(ft),nppstor2(ft),evp(ft),
     &dof(ft),evap,tran,roff,interc,evbs,flow1(ft),flow2(ft),year,mnth,
     &day,pet,laimax(ft),ht(ft),leafv,stemv,rootv,leaflit(ft),
     &stemlit(ft),rootlit(ft),bb(ft),ss(ft),bbgs(ft),dsbb(ft),
     &tmem,maxlai(ft),thty_dys,wtfc,wtwp,leafnpp(ft),stemnpp(ft),
     &rootnpp(ft),yield(ft),ft,resp,smtrig,qdirect,qdiff,suma,
     &tsumam(ft),stemfr(ft),lmor_sc(:,ft),nleaf,chill(ft),dschill(ft),
     &fpr,gsm(ft),swr(mnth,day),tleaf_n,tleaf_p,
     &ncalc_type,leaf_nit,vcmax(:,ft),jmax(:,ft),pnlc(:,ft),enzs(:,ft),
     &vcmax_type,
     &leafresp,rootresp,stemresp,daynpp,read_par,kg(ft),
     &ftcan_clump(ft),tassim,tgs,tci,hw_j,cstype,phen_cor,subd_par,
     &env_vcmax(ft),env_jmax(ft),soilp_map,can2g,canga,ga,
     &ftvna(ft),ftvnb(ft),ftjva(ft),ftjvb(ft),ftg0(ft),ftg1(ft),
     &no_slw_lim,par_loops,s070607,gs_func,
     &ce_light(:,:,ft),ce_ci(:,:,ft),ce_t,
     &ce_maxlight(:,:,ft),ce_ga(:,:,ft),ce_rh,
     &sl,hrs,ttype,calc_zen,iyear,
     &ftToptV(ft),ftHaV(ft),ftHdV(ft),ftToptJ(ft),ftHaJ(ft),ftHdJ(ft))
      
C PCM2      print*,'d2 ','tmp,prc,hum,cld,ft,soilc(ft),s1(ft),year,mnth,day'
C PCM2      print*,'d2 ',tmp(mnth,day),prc(mnth,day),hum(mnth,day),cld(mnth), 
C PCM2     &ft,soilc(ft),s1(ft),year,mnth,day
C PCM2      stop

!            write(*,*) mnth,day,tleaf_n
c      cbal = -1*(daygpp) + dayra + leafresp + rootresp + stemresp +
c     &         leafnpp(ft) + rootnpp(ft) + stemnpp(ft) + 
c     &         leafold - stemold - rootold + 
c     &         nppstore(ft) - nppsold
c      cbal = (cbal**2)**0.5
c      if(cbal.gt.0.000001) then
c       write(*,'(i3,15(1x,f8.4))') day, cbal,daygpp,daynpp,dayra,
c     &         leafresp, rootresp, stemresp, leafnpp(ft),  
c     &         rootnpp(ft), stemnpp(ft),leafold, stemold,  
c     &         rootold, nppstore(ft), nppsold
c      endif

*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
* Set daily memories for output.                                       *
*----------------------------------------------------------------------*

            leafv_sum = 0.0
            stemv_sum = 0.0
            rootv_sum = 0.0

            if((fttags(ft).eq.'BARE').or.fttags(ft).eq.'CITY') then
              leafg     = 0.0
            else
              if(s070607.eq.1) then
                leafg = leafv(1)*12.0d0/ftsla(ft)/18.0d0
              else
                leafg = leafv(1)/(ftsla(ft)/0.480d0)
              endif 
 
              do xi=1,ftlls(ft)
                if(s070607.eq.1) then
                 leafv_sum = leafv_sum+leafv(xi)*12.0d0/ftsla(ft)/18.0d0
                else
                 leafv_sum = leafv_sum+leafv(xi)/(ftsla(ft)/0.480d0)
                endif 
              enddo
              do xi=1,ftsls(ft)
                stemv_sum = stemv_sum+stemv(xi) 
              enddo
              do xi=1,ftrls(ft)
                rootv_sum = rootv_sum+rootv(xi) 
              enddo
            endif


            !these daily_out arrays must have the same numbering as the <otags> and 
            !associated vectors in <sdgvm1.f> 
            daily_out(1,ft,mnth,day) = lai(ft)                        ! LAI
            daily_out(2,ft,mnth,day) = roff                           ! runoff/drainage
            daily_out(3,ft,mnth,day) = evap+tran                      ! evapotranspiration
            daily_out(4,ft,mnth,day) = tran                           ! transpiration
            daily_out(5,ft,mnth,day) =(leafnpp(ft) + stemnpp(ft) +    ! NPP
     &rootnpp(ft) + nppstore(ft) - leafold - stemold - rootold - 
     &nppsold)*12.0d0
            daily_out(6,ft,mnth,day) = daygpp*12.0d0                  ! GPP
            !PRINT *,"daily_gpp",ft,mnth,day,daily_out(6,ft,mnth,day)
            daily_out(7,ft,mnth,day) = srespm/                        ! heterotrophic respiration
     &real(no_days(year,mnth,thty_dys))
            if(srespm.lt.1e-6) daily_out(7,ft,mnth,day) = 0.000       ! heterotrophic respiration
            daily_out(8,ft,mnth,day) = daily_out(5,ft,mnth,day) -     ! NEE
     &daily_out(7,ft,mnth,day)
            daily_out(9,ft,mnth,day) = tmp(mnth,day)
            daily_out(10,ft,mnth,day) = prc(mnth,day)
            daily_out(11,ft,mnth,day) = hum(mnth,day)
            daily_out(12,ft,mnth,day) = nppstore(ft)*12.0d0
            daily_out(13,ft,mnth,day) = (s1(ft) + s2(ft) + s3(ft) + 
     &s4(ft) - sw(1) - sw(2) - sw(3) - sw(4))/(sfc(1) + sfc(2) + 
     &sfc(3) + sfc(4) - sw(1) - sw(2) - sw(3) - sw(4))
*            daily_out(13,ft,mnth,day) = 
*     &(s1(ft) - sw(1))/(sfc(1) - sw(1))
            daily_out(14,ft,mnth,day) = pet
            daily_out(15,ft,mnth,day) = interc
            daily_out(16,ft,mnth,day) = evbs + sl ! this is soil evap plus sublimation from snow
            daily_out(17,ft,mnth,day) = min(1.0d0,s1(ft)/10.0d0/topsl)
            daily_out(18,ft,mnth,day) = s1(ft) + s2(ft) + s3(ft) + 
     &s4(ft)
            daily_out(19,ft,mnth,day) = 12.0d0*(resp+dayra)
            daily_out(20,ft,mnth,day) = qdirect
            daily_out(21,ft,mnth,day) = qdiff
            daily_out(22,ft,mnth,day) = nleaf
            daily_out(23,ft,mnth,day) = leaflit(ft) - lflitold
            daily_out(24,ft,mnth,day) = cld(mnth)
            !diverges from Mark's version here
            daily_out(25,ft,mnth,day) = lsn(ft)+sn(ft) ! liquid snow and snow water equivalent
            daily_out(26,ft,mnth,day) = nleaf
            daily_out(27,ft,mnth,day) = leaf_nit
            daily_out(28,ft,mnth,day) = vcmax(1,ft)
            daily_out(29,ft,mnth,day) = jmax(1,ft)
            daily_out(30,ft,mnth,day) = kg(ft)        !soil water limitation factor (unitless)
            daily_out(31,ft,mnth,day) = tassim        !topleaf assimilation (umol CO2 m-2 s-1)
            daily_out(32,ft,mnth,day) = tgs           !topleaf stomatal conductance (umol CO2 m-2 s-1)
            daily_out(33,ft,mnth,day) = tci           !topleaf internal CO2 conc (Pa)
            !transfered variables from previous FACE version 
            daily_out(34,ft,mnth,day) = fpr           ! fAPAR
            daily_out(35,ft,mnth,day) = leafresp*12.0d0   ! leaf growth respiration
            daily_out(36,ft,mnth,day) = stemresp*12.0d0   ! stem maintenance respiration
            daily_out(37,ft,mnth,day) = rootresp*12.0d0   ! root maintenance respiration
            daily_out(38,ft,mnth,day) = dayra*12.0d0      ! leaf maintentance resp in grams
            daily_out(39,ft,mnth,day) = leafv_sum         ! leaf biomass C in grams
            daily_out(40,ft,mnth,day) = stemv_sum*12.0d0  ! live sapwood biomass C in grams
            daily_out(41,ft,mnth,day) = rootv_sum*12.0d0  ! live root biomass C in grams
            daily_out(42,ft,mnth,day) = leafg             ! leaf growth in grams
            daily_out(43,ft,mnth,day) = stemv(1)*12.0d0   ! stem growth in grams
            daily_out(44,ft,mnth,day) = rootv(1)*12.0d0   ! root growth in grams
            !the following conductance outputs are in mol H2O m-2 s-1 
            daily_out(45,ft,mnth,day) = canga/(8.3144*(tmp(mnth,day)+
     &273.15)/101325.0)                                   !canopy aerodynamic conductance (H2O units)
            daily_out(46,ft,mnth,day) = ga*1.3            !leaf boundary layer conductance (H2O units)
            daily_out(47,ft,mnth,day) = can2g*1.6         !canopy stomatal conductance (H2O units)
            daily_out(73,ft,mnth,day) = qdirect + qdiff   !total PAR (H2O units)
            daily_out(74,ft,mnth,day) = ca(mnth,day)      !co2 partial pressure (Pa)
            daily_out(75,ft,mnth,day) = hrs               !daylight hours (hours)


            IF ((bb(ft).EQ.day+(mnth-1)*30).AND.(budo(ft).EQ.0))
     &budo(ft) = bb(ft)
            IF ((ss(ft).EQ.day+(mnth-1)*30).AND.(seno(ft).EQ.0))
     &seno(ft) = ss(ft)

*----------------------------------------------------------------------*
* Sumation of evaporation and transpiration.                           *
*----------------------------------------------------------------------*
            evapm(mnth,ft) = evapm(mnth,ft) + evap
            tranm(mnth,ft) = tranm(mnth,ft) + tran
            roffm(mnth,ft) = roffm(mnth,ft) + roff
            petm(mnth,ft) = petm(mnth,ft) + pet
            avmnpet(ft) = avmnpet(ft) + pet

*----------------------------------------------------------------------*
* Sumation of GPP and NPP.                                             *
*----------------------------------------------------------------------*

            photm(mnth,ft) = photm(mnth,ft) + daygpp*12
            laimnth(mnth,ft) = laimnth(mnth,ft) + 
     &           lai(ft)/real(no_days(year,mnth,thty_dys))

            DO d=1,ftlls(ft)
              leafdp(d,ft)=leafv(d)
            ENDDO
            DO d=1,ftsls(ft)
              stemdp(d,ft)=stemv(d)
            ENDDO
            DO d=1,ftrls(ft)
              rootdp(d,ft)=rootv(d)
            ENDDO

            DO d=1,30
              sm_trig(d,ft) = smtrig(d)
            ENDDO

            DO d=1,360
              sumadp(d,ft) = suma(d)
            ENDDO

          ! dolydo IF
          ELSE 
            DO i=1,douts
              daily_out(i,ft,mnth,day) = 0.0d0
            ENDDO
          ENDIF


*----------------------------------------------------------------------*
*         End of ft loop
*----------------------------------------------------------------------*
          ENDDO

*----------------------------------------------------------------------*
*       End of daily loop 
*----------------------------------------------------------------------*
        ENDDO        

*----------------------------------------------------------------------*
*      Monthly Operation                                               *
*----------------------------------------------------------------------*
        DO ft=1,nft
        IF (dolydo(ft).EQ.1) THEN

          h2o = (s1(ft) + s2(ft) + s3(ft) + s4(ft))
          IF (h2o.LT.0.0d0) h2o = 0.1d0
          DO i=1,8
            c0v(i)=c0(i,ft)
            n0v(i)=n0(i,ft)
          ENDDO
          DO i=1,3
            minnv(i)=minn(i,ft)
          ENDDO

*          print*,avmnpet(ft),avmnppt,avmnt
          CALL DOLYMONTH(ts,tc,avmnpet(ft),avmnppt,avmnt,h2o,flow1(ft),
     &flow2(ft),c0v,n0v,minnv,nfix,nci,dslc,drlc,dsln,srespm,lchm,ca,
     &site,year,yr0,yrf,speedc,soilc(ft),soiln(ft),mnth,w_scalar,
     &t_scalar,fl,cal)

          sresp(ft) = sresp(ft) + srespm
          lch(ft) = lch(ft) + lchm
               
          DO i=1,8
            c0(i,ft)=c0v(i)
            n0(i,ft)=n0v(i)
          ENDDO
          f0(:,ft) = fl(:)

          DO i=1,3
            minn(i,ft)=minnv(i)
          ENDDO

          !traceability analysis soil pool outputs 
          daily_out(71,ft,mnth,:) = w_scalar
          daily_out(72,ft,mnth,:) = t_scalar
               
        ENDIF
        ENDDO
        
*----------------------------------------------------------------------*
*                             End of month loop                         *
*----------------------------------------------------------------------*
      ENDDO                     ! end of monthly loop

      !print*, ce_ci(30,1,:)

      DO ft=1,nft

        IF (dolydo(ft).EQ.1) THEN

          IF (abs(leafnpp(ft) + stemnpp(ft) + rootnpp(ft) + 

     &nppstore(ft) - nppstoreold(ft)).GT.1.0E-6) THEN
            npl(ft) = leafnpp(ft)/(leafnpp(ft) + stemnpp(ft) + 
     &rootnpp(ft) + nppstore(ft) - nppstoreold(ft))*100.0d0
            nps(ft) = stemnpp(ft)/(leafnpp(ft) + stemnpp(ft) + 
     &rootnpp(ft) + nppstore(ft) - nppstoreold(ft))*100.0d0
            npr(ft) = rootnpp(ft)/(leafnpp(ft) + stemnpp(ft) + 
     &rootnpp(ft) + nppstore(ft) - nppstoreold(ft))*100.0d0
          ELSE
            npl(ft) = 0.0d0
            nps(ft) = 0.0d0
            npr(ft) = 0.0d0
            npp(ft) = 0.0d0
          ENDIF
          npp(ft) = (leafnpp(ft) + stemnpp(ft) + rootnpp(ft) + 
     &nppstore(ft) - nppstoreold(ft))*12.0d0
*          print'(4f8.1)',npp(ft),npl(ft),nps(ft),npr(ft)
*          print*,leafnpp(ft)*12.0d0,stemnpp(ft)*12.0d0,
*     &rootnpp(ft)*12.0d0

c     put c0,n0 and minn into vector
          DO i=1,8
            c0v(i)=c0(i,ft)
            n0v(i)=n0(i,ft)
          ENDDO
          DO i=1,3
            minnv(i)=minn(i,ft)
          ENDDO

          CALL SUMDV(tc0,tn0,tminn,ts1,ts2,ts3,ts4,tsn,tlsn,
     &c0v,n0v,minnv,s1(ft),s2(ft),s3(ft),s4(ft),
     &sn(ft),lsn(ft),ftcov(ft))

c     check water cycle closure
          prctot = 0.0d0
          yrtran = 0.0d0
          yrevap = 0.0d0
          yrroff = 0.0d0
          yrpet = 0.0d0
          gpp(ft) = 0.0d0
          DO mnth=1,12
            DO day=1,no_days(year,mnth,thty_dys)
              prctot = prctot + prc(mnth,day)
            ENDDO
            yrtran = yrtran + tranm(mnth,ft)
            yrevap = yrevap + evapm(mnth,ft)
            yrroff = yrroff + roffm(mnth,ft)
            yrpet = yrpet + petm(mnth,ft)
            gpp(ft) = gpp(ft) + photm(mnth,ft)
          ENDDO

          trn(ft) = yrtran
          evt(ft) = yrtran + yrevap
          rof(ft) = yrroff
          fpet(ft) = yrpet

          lai(ft) = laimax(ft)

*     Convert nppstore back to grams
          nppstore(ft) = nppstore(ft)*12.0d0
          nppstorx(ft) = nppstorx(ft)*12.0d0
          nppstor2(ft) = nppstor2(ft)*12.0d0

        ELSE ! not dolydo
            DO i=1,10
              DO mnth=1,12
                DO day=1,no_days(year,mnth,thty_dys)
                  daily_out(i,ft,mnth,day) = 0.0d0
                ENDDO
              ENDDO
            ENDDO
            npp(ft) = 0.0d0
            lai(ft) = 0.0d0
            laimax(ft) = 0.0d0
            evp(ft) = 0.0d0
            dof(ft) = 0.0d0
            npl(ft) = 0.0d0
            nps(ft) = 0.0d0
            npr(ft) = 0.0d0
            rof(ft) = 0.0d0
            trn(ft) = 0.0d0
            sresp(ft) = 0.0d0
            evt(ft) = 0.0d0
            fpet(ft) = 0.0d0
            gpp(ft) = 0.0d0
            lch(ft) = 0.0d0
            bb(ft) = 0
            ss(ft) = 0
            bbgs(ft) = 0
            dsbb(ft) = 0
            nppstore(ft) = 0.0d0
            nppstorx(ft) = 0.0d0
            nppstor2(ft) = 0.0d0

            DO day=1,ftlls(ft)
               leafdp(day,ft) = 0.0d0
            ENDDO
            DO day=1,ftsls(ft)
               stemdp(day,ft) = 0.0d0
            ENDDO
            DO day=1,ftrls(ft)
               rootdp(day,ft) = 0.0d0
            ENDDO
         ENDIF                  ! doly do
      ENDDO                     ! ft's
      !from here below is indented 2 columns more
      npp(1) = 0.0d0
      lai(1) = 0.0d0

      CALL GROWTH(nft,ftmor,ftwd,ftxyl,ftpd,ftgr0,ftgrf,cov,bio,
     &nppstore,npp,lai,nps,npr,evp,slc,rlc,sln,rln,stembio,rootbio,ppm,
     &hgt,leaflit)

      CALL SWAP(ic0,in0,iminn,is1,is2,is3,is4,isn,ilsn,tc0,tn0,
     &tminn,ts1,ts2,ts3,ts4,tsn,tlsn)

*----------------------------------------------------------------------*
* Average outputs by cover proportions.                                *
*----------------------------------------------------------------------*
        avlai   = 0.0d0
	avnpp   = 0.0d0
        avnppst = 0.0d0
	avdof   = 0.0d0
	avrof   = 0.0d0
	avtrn   = 0.0d0
	avsresp = 0.0d0
	avevt   = 0.0d0
	avpet   = 0.0d0
        avgpp   = 0.0d0
        avlch   = 0.0d0
        avyield = 0.0d0
        avnleaf = 0.0d0
        avleaf_nit = 0.0d0
        avvcmax = 0.0d0
        avjmax  = 0.0d0
        avsla   = 0.0d0
        respsum = 0.0d0
        kg_beta = 0.0d0
        av_hgt_all= 0.0d0
        qdirsum   = 0.0d0
        qdifsum   = 0.0d0
        yearsoilw = 0.0d0

        DO ft=1,nft
          avlai   = avlai   + ftcov(ft)*lai(ft)
          avnpp   = avnpp   + ftcov(ft)*npp(ft)
          avnppst = avnppst + ftcov(ft)*nppstore(ft)

          avdof   = avdof   + ftcov(ft)*dof(ft)
          avrof   = avrof   + ftcov(ft)*rof(ft)
          avtrn   = avtrn   + ftcov(ft)*trn(ft)
          avsresp = avsresp + ftcov(ft)*sresp(ft)
          avevt   = avevt   + ftcov(ft)*evt(ft)
          avpet   = avpet   + ftcov(ft)*fpet(ft)
          avgpp   = avgpp   + ftcov(ft)*gpp(ft)
          avlch   = avlch   + ftcov(ft)*lch(ft)
          avyield = avyield + ftcov(ft)*yield(ft)
          hi = 0
          do i = 1,ftmor(ft)
            av_hgt(ft) = av_hgt(ft) + hgt(i,ft)
            if(hgt(i,ft).gt.0.0d0) hi = hi + 1 
          enddo
          if(hi.ge.1) av_hgt(ft) = av_hgt(ft)/hi 
          av_hgt_all = av_hgt_all + ftcov(ft)*av_hgt(ft)
        ENDDO

        gi=0.000d0 ! gi tracks how many days in the year at least one PFT had LAI > 1
        wi=0.000d0 ! wi tracks how many days in the year at least one PFT had GPP > 1
        DO mnth=1,12
           DO day=1,no_days(year,mnth,thty_dys)
             active_cov  = 0.d0
             DO ft=1,nft
              ! need to use a method that doesn't just weight by cov but
              ! also consideres other PFTs and whether they are active or
              ! not. e.g. if a PFT with large fractional cover has no LAI
              ! then the small fractional cover of the PFT that does will
              ! give a low value of vcmax even if it's the only PFT with
              ! LAI. 
c             active_cov1 = 0.d0
              if(daily_out(1,ft,mnth,day).gt.1.0d0) then
                active_cov = active_cov + ftcov(ft)
              endif
c             if(daily_out(6,ft,mnth,day).gt.1.0d0) active_cov1=ftcov(ft)
             ENDDO

c             print*, active_cov

             fti=1 
             DO ft=1,nft
             ! sums the below variables when lai>1 (because weird things happen to leaf N when LAI<1)
             ! weighted by pft cover as a proportion of 'active' PFT cover, defined here as LAI > 1  
             !if(daily_out(1,ft,mnth,day).gt.1.0d0) then
             if((daily_out(1,ft,mnth,day).gt.1.0d0).and.
     &(active_cov.gt.1.d-3)) then
               avnleaf   = avnleaf    + (ftcov(ft)/active_cov) * 
     &daily_out(26,ft,mnth,day)
               avleaf_nit= avleaf_nit + (ftcov(ft)/active_cov) * 
     &daily_out(27,ft,mnth,day)
               avvcmax   = avvcmax    + (ftcov(ft)/active_cov) * 
     &daily_out(28,ft,mnth,day)
               avjmax    = avjmax     + (ftcov(ft)/active_cov) * 
     &daily_out(29,ft,mnth,day)
               kg_beta   = kg_beta    + (ftcov(ft)/active_cov) *
     &daily_out(30,ft,mnth,day)
               avsla     = avsla      + (ftcov(ft)/active_cov) * 
     &ftsla(ft)

               if(fti.eq.1) gi=gi+1.d0
               fti=fti+1
             endif
             ENDDO

             fti=1 
             DO ft=1,nft
             ! sums the soil water stress scalar when gpp>1, weighted by pft cover as a proportion of 'active' PFT cover, defined here as GPP > 1  
c             if(daily_out(6,ft,mnth,day).gt.1.0d0) then
c               kg_beta  = kg_beta  + ftcov(ft)/active_cov1 *
c     &daily_out(30,ft,mnth,day)
c               if(fti.eq.1) wi=wi+1
c               fti=fti+1
c             endif
             yearsoilw= yearsoilw+ ftcov(ft) * daily_out(18,ft,mnth,day)
             respsum  = respsum  + ftcov(ft) * daily_out(19,ft,mnth,day)
             qdirsum  = qdirsum  + ftcov(ft) * daily_out(20,ft,mnth,day)
             qdifsum  = qdifsum  + ftcov(ft) * daily_out(21,ft,mnth,day)
             ENDDO
           ENDDO 
        ENDDO
        
        oscale = 0.0d0
        DO mnth=1,12       
          oscale = oscale + real(no_days(year,mnth,thty_dys))   
        ENDDO

        if(gi.lt.1.d-1) gi = 1d0
        !if(wi.lt.1.d-1) wi = 1d0

        oscale     = 1.0/oscale
        avnleaf    = avnleaf/gi
        avleaf_nit = avleaf_nit/gi
        avvcmax    = avvcmax/gi
        avjmax     = avjmax/gi
        avsla      = avsla/gi
c       kg_beta    = kg_beta/wi
        kg_beta    = kg_beta/gi
        max_hgt    = maxval(hgt(:,:))
        yearsoilw  = yearsoilw*oscale
        !convert year soil W from mm to prop by volume
        yearsoilw  = yearsoilw / (dep*10.0d0)

*----------------------------------------------------------------------*

        tsoilc = 0.0d0
        tslc   = 0.0d0
        trlc   = 0.0d0
        tsoiln = 0.0d0
        tabglitterc = 0.0000 
        tblgc       = 0.0000
        DO i=1,8
          tsoilc = tsoilc + ic0(i)
          tsoiln = tsoiln + in0(i)
        ENDDO
        tsoiln = tsoiln + iminn(3)
        DO ft=1,nft
          !print*,tslc,trlc
          !tslc = tslc + ftcov(ft) * slc(ft)
          !trlc = trlc + ftcov(ft) * rlc(ft)
          tslc = tslc + slc(ft)
          trlc = trlc + rlc(ft)
        ENDDO
        tsoilc      = tsoilc + tslc + trlc
        !print*, tslc , ic0(1) , ic0(5)
        !print*, tslc + ic0(1) + ic0(5)
        tabglitterc = tslc + ic0(1) + ic0(5)
        tblgc       = trlc + ic0(2) + ic0(3)+ic0(4)+ic0(6)+ic0(7)+ic0(8)

        swcold = swcnew
        swcnew = is1 + is2 + is3 + is4 + isn + ilsn
        sumbio = 0.0d0
        bioind = 0
        covind = ca(1,1)
        maxbio = 0.0d0
        maxcov = 0.0d0
        leafper = 0.0d0
        stemper = 0.0d0
        rootper = 0.0d0
        soilcn = tsoilc/tsoiln

        tbioleaf = 0.0000
        DO ft=1,nft
          bioleaf(ft) = 0.0d0
          if(s070607.eq.1) then
            DO day=1,ftlls(ft)
              bioleaf(ft) = bioleaf(ft) + 
     &leafdp(day,ft)*12.0d0/ftsla(ft)/18.0d0
            ENDDO
          else
            DO day=1,ftlls(ft)
              bioleaf(ft) = bioleaf(ft) + 
     &leafdp(day,ft)/(ftsla(ft)/0.480d0)
            ENDDO
          endif
          tbioleaf = bioleaf(ft) 
        ENDDO

        avppm = 0.0
        DO ft=1,nft
          bioo(ft) = 0.0d0
          covo(ft) = 0.0d0
          DO i=1,ftmor(ft)
            bioo(ft) = bioo(ft) + (bio(i,1,ft) +
     &bio(i,2,ft) + bioleaf(ft) + nppstore(ft))*cov(i,ft)
            covo(ft) = covo(ft) + cov(i,ft)
            leafper = leafper + npl(ft)*cov(i,ft)
            stemper = stemper + nps(ft)*cov(i,ft)
            rootper = rootper + npr(ft)*cov(i,ft)
            IF (ftmor(ft).gt.50) THEN
              avppm = avppm + ppm(i,ft)*cov(i,ft)
            ENDIF
          ENDDO
          sumbio = sumbio + bioo(ft)
          IF (covo(ft).GT.maxcov) THEN
            maxcov = covo(ft)
            !covind = ft
          ENDIF
          IF (bioo(ft).GT.maxbio) THEN
            maxbio = bioo(ft)
            bioind = ft
          ENDIF
        ENDDO


      !print*,(((max(bio(i,j,ft),0.0),j=1,2),i=1,ftmor(ft)),ft=2,nft)
      !print*, 'Biomass:'
      !print*,(max(bioo(ft),0.0),ft=2,nft)
      !print*,(bioleaf(ft),ft=2,nft) 
      !print*,(ftsla(ft),ft=2,nft) 
      !print*,(nppstore(ft),ft=2,nft) 
      !print*, ''


*----------------------------------------------------------------------*
* Write outputs.                                                       *
*----------------------------------------------------------------------*
* General. Primarily (maybe all) annual.                                                             *
*----------------------------------------------------------------------*

        IF (iyear.GE.nyears-outyears+1) THEN
          WRITE(21,'('' '',f8.1,$)') avlai
          WRITE(22,'('' '',f8.1,$)') avnpp
          WRITE(23,'('' '',f8.1,$)') tsoilc
          WRITE(24,'('' '',f8.1,$)') tsoiln
          WRITE(25,'('' '',f8.1,$)') avnpp-avsresp
          WRITE(26,'('' '',f8.1,$)') min(swcnew,9999.0d0)
          WRITE(27,'('' '',f8.1,$)') sumbio
          WRITE(28,'('' '',i2,$)')   bioind
          WRITE(29,'('' '',f6.2,$)')  covind
          WRITE(30,'('' '',f8.1,$)') avdof
          WRITE(31,'('' '',f8.1,$)') avrof
          WRITE(32,'('' '',f8.2,$)') firec
          WRITE(33,'('' '',f8.1,$)') avnppst
          WRITE(34,'('' '',f8.1,$)') stembio
          WRITE(35,'('' '',f8.1,$)') rootbio
          WRITE(36,'('' '',f8.1,$)') leafper
          WRITE(37,'('' '',f8.1,$)') stemper
          WRITE(38,'('' '',f8.1,$)') rootper
          WRITE(39,'('' '',f8.1,$)') avsresp
          WRITE(40,'('' '',f8.1,$)') avevt
          WRITE(41,'('' '',f8.1,$)') avgpp
          WRITE(42,'('' '',f8.3,$)') avlch
          WRITE(43,'('' '',f8.1,$)') yearprc
          WRITE(44,'('' '',f8.1,$)') avnpp-avsresp-firec-avlch-avyield-
     &flulccc
          WRITE(45,'('' '',f8.1,$)') avtrn
          WRITE(46,'('' '',f8.5,$)') fprob
          WRITE(47,'('' '',f8.2,$)') yeartmp
          WRITE(48,'('' '',f8.2,$)') yearhum
          WRITE(49,'('' '',f8.5,$)') avppm
          WRITE(50,'('' '',f8.1,$)') avpet
          WRITE(51,'('' '',f8.1,$)') tslc 
          WRITE(52,'('' '',f8.1,$)') trlc
          WRITE(53,'('' '',f8.4,$)') tleaf_n
          WRITE(54,'('' '',f8.5,$)') tleaf_p
          WRITE(55,'('' '',f8.4,$)') avsla 
          WRITE(56,'('' '',f8.3,$)') max_hgt 
          WRITE(57,'('' '',i2,$)')   vcmax_type 
          WRITE(58,'('' '',f8.3,$)') avnleaf
          WRITE(59,'('' '',f8.3,$)') avleaf_nit                  
          WRITE(60,'('' '',f8.3,$)') avvcmax         
          WRITE(61,'('' '',f8.3,$)') avjmax
          WRITE(62,'('' '',f8.3,$)') respsum         
          WRITE(63,'('' '',f8.3,$)') yearswr         
          WRITE(64,'('' '',f8.3,$)') avgpp - avnpp - respsum
          WRITE(65,'('' '',f8.3,$)') masoilcn          
          WRITE(66,'('' '',f8.3,$)') kg_beta     
          
          IF(read_clump.eq.2) THEN
            WRITE(67,'('' '',f8.3,$)') ftcan_clump(1)
          ELSE
            WRITE(67,'('' '',f8.3,$)') -99.99
          ENDIF
          
          WRITE(68,'('' '',f8.1,$)') qdirsum          
          WRITE(69,'('' '',f8.1,$)') qdirsum + qdifsum          
          WRITE(601,'('' '',f8.6,$)') wtfc
          WRITE(602,'('' '',f8.6,$)') wtwp
          WRITE(603,'('' '',f12.2,$)') tabglitterc
          WRITE(604,'('' '',f12.2,$)') tblgc
          WRITE(605,'('' '',f10.2,$)') tbioleaf
          WRITE(606,'('' '',f10.2,$)') flulccc 
        ENDIF


*----------------------------------------------------------------------*
* Write optional cov bio bud sen.                                      *
*----------------------------------------------------------------------*
        iofn = iofngft
        IF (iyear.GE.nyears-outyears2+1) THEN
          DO ft=1,nft
            IF (out_cov) THEN
              iofn = iofn + 1
              WRITE(iofn,'('' '',f8.6,$)') covo(ft)
            ENDIF
            IF (out_bio) THEN
              iofn = iofn + 1
              WRITE(iofn,'('' '',f8.1,$)') bioo(ft)
            ENDIF
            IF (out_bud) THEN
              iofn = iofn + 1
              WRITE(iofn,'('' '',i8,$)') budo(ft)
            ENDIF
            IF (out_sen) THEN
              iofn = iofn + 1
              WRITE(iofn,'('' '',i8,$)') seno(ft)
            ENDIF
          ENDDO
        ENDIF

*----------------------------------------------------------------------*
* Write monthly PIXEL outputs.                                         *
*----------------------------------------------------------------------*
        iofn = 400
*        print*,'thty_dys'
        IF (iyear.GE.nyears-outyears1+1) THEN
          IF ((oymd.EQ.1).OR.(oymd.EQ.3)) THEN
            DO i=1,douts
              IF (otagsn(i).EQ.1) THEN 
                DO mnth=1,12
                  ans(mnth,1) = 0.0d0
                ENDDO
                DO ft=1,nft
                  DO mnth=1,12
                    IF (omav(i).eq.1) THEN
                      oscale = 1.0/real(no_days(year,mnth,thty_dys))
                    ELSE
                      oscale = 1.0
                    ENDIF
                    DO day=1,no_days(year,mnth,thty_dys)
                      ans(mnth,1) = ans(mnth,1) + 
     &daily_out(i,ft,mnth,day)*ftcov(ft)*oscale
                    ENDDO
                    IF(i.eq.6) THEN
!                     PRINT *,"MOPIX GPP",mnth,i,ft,ftcov(ft),ans(mnth,1)
                    ENDIF
                  ENDDO
                ENDDO
                iofn = iofn + 1
                sum1 = 0.0d0
                IF (omav(i).eq.1) THEN
                  oscale = 1.0d0/12.0d0
                ELSE
                  oscale = 1.0d0
                ENDIF
                DO mnth=1,12
                  sum1 = sum1 + ans(mnth,1)*oscale
                ENDDO
                WRITE(iofn,ofmt(i)) year,(ans(mnth,1),
     &mnth=1,12),sum1
              ENDIF

            ENDDO
          ENDIF

*----------------------------------------------------------------------*
* Write daily PIXEL outputs.                                           *
*----------------------------------------------------------------------*
          IF ((oymd.EQ.2).OR.(oymd.EQ.3)) THEN
            DO i=1,douts
              IF (otagsn(i).EQ.1) THEN
                iofn = iofn + 1
                WRITE(iofn,'(i4)') year
                DO mnth=1,12
                  DO day=1,no_days(year,mnth,thty_dys)
                    ans(mnth,day) = 0.0d0
                  ENDDO
                ENDDO
                DO ft=1,nft
                  DO mnth=1,12
                    DO day=1,no_days(year,mnth,thty_dys)
                      ans(mnth,day) = ans(mnth,day) + 
     &daily_out(i,ft,mnth,day)*ftcov(ft)
                    ENDDO
                  ENDDO
                ENDDO
                DO mnth=1,12
                  WRITE(iofn,ofmt(100+i)) (ans(mnth,day),
     &day=1,no_days(year,mnth,thty_dys))
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDIF

        iofn = iofnft
        IF (iyear.GE.nyears-outyears2+1) THEN
*----------------------------------------------------------------------*
* Write monthly SUBPIXEL outputs.                                      *
*----------------------------------------------------------------------*
          IF ((oymdft.EQ.1).OR.(oymdft.EQ.3)) THEN
            DO i=1,douts
              IF (otagsnft(i).EQ.1) THEN 
                DO ft=1,nft
                  DO mnth=1,12
                    ans(mnth,1) = 0.0d0
                  ENDDO
                  DO mnth=1,12
                    IF (omav(i).eq.1) THEN
                      oscale = 1.0/real(no_days(year,mnth,thty_dys))
                    ELSE
                      oscale = 1.0
                    ENDIF
                    DO day=1,no_days(year,mnth,thty_dys)
                      ans(mnth,1) = ans(mnth,1) + 
     &daily_out(i,ft,mnth,day)*oscale
                    ENDDO
                  ENDDO
                  iofn = iofn + 1
                  WRITE(iofn,ofmt(i)) year,(ans(mnth,1),
     &mnth=1,12)
                ENDDO
              ENDIF
            ENDDO
          ENDIF

*----------------------------------------------------------------------*
* Write daily SUBPIXEL outputs.                                        *
*----------------------------------------------------------------------*
          IF ((oymdft.EQ.2).OR.(oymdft.EQ.3)) THEN
            DO i=1,douts
              IF (otagsnft(i).EQ.1) THEN
                DO ft=1,nft
                  iofn = iofn + 1
                  WRITE(iofn,'(i4)') year
                  DO mnth=1,12
                    DO day=1,no_days(year,mnth,thty_dys)
                      ans(mnth,day) = 0.0d0
                    ENDDO
                  ENDDO
                  DO mnth=1,12
                    DO day=1,no_days(year,mnth,thty_dys)
                      ans(mnth,day) = ans(mnth,day) + 
     &daily_out(i,ft,mnth,day)
                    ENDDO
                  ENDDO
                  DO mnth=1,12
                    WRITE(iofn,ofmt(100+i)) (ans(mnth,day),
     &day=1,no_days(year,mnth,thty_dys))
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDIF

*----------------------------------------------------------------------*
* Write snapshots
*----------------------------------------------------------------------*
        IF (snp_no.GT.0) THEN
          IF (snp_year.LE.snp_no) THEN
            IF (year.EQ.snpshts(snp_year)) THEN
              fno = 100 + (snp_year-1)*4
              WRITE(fno+1,'(F7.3,F9.3,4500F9.1)') 
     &lat,lon,(((bio(i,j,k),j=1,2),i=1,ftmor(k)),k=1,nft)
              WRITE(fno+2,'(F7.3,F9.3,1500F12.9)') 
     &lat,lon,((cov(i,j),i=1,ftmor(j)),j=1,nft)
              WRITE(fno+3,'(F7.3,F9.3,1500F12.7)') 
     &lat,lon,((ppm(i,j),i=1,ftmor(j)),j=1,nft)
              WRITE(fno+4,'(F7.3,F9.3,1500F8.3)') 
     &lat,lon,((hgt(i,j),i=1,ftmor(j)),j=1,nft)
              snp_year = snp_year + 1
            ENDIF
          ENDIF
        ENDIF

*----------------------------------------------------------------------*
* Carbon at the end of the year.                                       *
*----------------------------------------------------------------------*

        ans1 = 0.0d0
        DO ft=1,nft
          DO i=1,ftmor(ft)
            ! adding stem_carbon=bio(i,1,ft) & root_carbon(i,2,ft) & leaf_carbon & nppstore
            ans1 = ans1 + (bio(i,1,ft) + bio(i,2,ft) + bioleaf(ft) +
     &nppstore(ft))*cov(i,ft)
          ENDDO
          ans1 = ans1 + slc(ft) + rlc(ft) ! adding stem_litter_carbon(ft) & root_litter_carbon(ft)
        ENDDO
        IF(debug .EQV. .TRUE.) THEN
          WRITE(*,'(''check0'',3f13.6)') ccheck
          WRITE(*,'(''avnpp'',3f12.6)') avnpp 
          WRITE(*,'(''ans1'',3f12.6)') ans1 
          WRITE(*,'(''tc0(1)'',3f12.6)') tc0(1) 
          WRITE(*,'(''tc0(2)'',3f12.6)') tc0(2) 
          WRITE(*,'(''tc0(3)'',3f12.6)') tc0(3) 
          WRITE(*,'(''tc0(4)'',3f12.6)') tc0(4) 
          WRITE(*,'(''tc0(5)'',3f12.6)') tc0(5) 
          WRITE(*,'(''tc0(6)'',3f12.6)') tc0(6) 
          WRITE(*,'(''tc0(7)'',3f12.6)') tc0(7) 
          WRITE(*,'(''tc0(8)'',3f12.6)') tc0(8) 
          WRITE(*,'(''avlch'',3f12.6)') avlch 
          WRITE(*,'(''avsresp'',3f12.6)') avsresp 
          WRITE(*,'(''firec'',3f12.6)') firec 
          WRITE(*,'(''avyield'',3f12.6)') avyield 
          WRITE(*,'(''flulccc'',3f12.6)') flulccc 
        ENDIF
!        ccheck = ccheck - (ans1 + tc0(1) + 
!     &tc0(2) + tc0(3) + tc0(4) + tc0(5) + tc0(6) + tc0(7) + tc0(8) - 
!     &(avnpp-avlch-avsresp-firec)) - avyield
        ccheck = ccheck + avnpp - (ans1 + tc0(1) + 
     &tc0(2) + tc0(3) + tc0(4) + tc0(5) + tc0(6) + tc0(7) + tc0(8) + 
     &avlch + avsresp + firec +  avyield + flulccc)
        IF(debug .EQV. .TRUE.) THEN
          WRITE(*,'(''check'',3f13.6)') ccheck
        ENDIF

*----------------------------------------------------------------------*
* Check carbon and water balance, write to 'DIAG' if any problems.     *
*----------------------------------------------------------------------*
        IF (abs(ccheck).GT.1.0e-3) THEN
          IF (.NOT.speedc) THEN
            WRITE(11,'(f12.6,'' C budget not closed, site '',f7.3,f9.3,
     &i6,'', year '',2i6)') ccheck,lat,lon,site,iyear,yearv(iyear)
            WRITE(*,'(''check'',3f12.6)') ccheck
          ENDIF
        ENDIF
      !if(site.gt.1629) write(*,'(360f8.2)') swr

      ENDDO
*----------------------------------------------------------------------*
*                             End of year loop                         *
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*                      End record for default output files             *
*----------------------------------------------------------------------*
        DO i=21,69
          WRITE(i,*)
        ENDDO
        DO i=601,606
          WRITE(i,*)
        ENDDO

        iofn = iofngft
        IF (outyears2.GT.0) THEN
          DO ft=1,nft
            IF (out_cov) THEN
              iofn = iofn + 1
              WRITE(iofn,*)
            ENDIF
            IF (out_bio) THEN
              iofn = iofn + 1
              WRITE(iofn,*)
            ENDIF
            IF (out_bud) THEN
              iofn = iofn + 1
              WRITE(iofn,*)
            ENDIF
            IF (out_sen) THEN
              iofn = iofn + 1
              WRITE(iofn,*)
            ENDIF
          ENDDO
        ENDIF

*----------------------------------------------------------------------*
*                             Write initialisation output              *
*----------------------------------------------------------------------*
        IF (initiseo) THEN
          WRITE(81,'(F7.3,F9.3)') lat,lon
          WRITE(82,'(F7.3,F9.3)') lat,lon
          WRITE(83,'(F7.3,F9.3)') lat,lon
          WRITE(84,'(F7.3,F9.3)') lat,lon
          WRITE(85,'(F7.3,F9.3)') lat,lon

          DO ft=2,nft
            WRITE(81,'(1000E13.6)') 
     &((max(bio(i,j,ft),0.0),j=1,2),i=1,ftmor(ft))
            WRITE(82,'(1000E13.6)') (max(cov(i,ft),0.0),i=1,ftmor(ft))
            WRITE(83,'(1000E13.6)') (max(ppm(i,ft),0.0),i=1,ftmor(ft))
            WRITE(84,'(1000E13.6)') (max(hgt(i,ft),0.0),i=1,ftmor(ft))
            WRITE(85,'(1000E13.6)') (max(wdt(i,ft),0.0),i=1,ftmor(ft))
          ENDDO

          WRITE(86,'(F7.3,F9.3)') lat,lon
          WRITE(87,'(F7.3,F9.3)') lat,lon
          WRITE(88,'(F7.3,F9.3)') lat,lon
          DO ft=2,nft
            WRITE(86,'(2000E13.6)' ) 
     &(max(leafdp(i,ft),0.0),i=1,ftlls(ft))
            WRITE(87,'(1000E13.6)') 
     &(max(stemdp(i,ft),0.0),i=1,ftsls(ft))
            WRITE(88,'(1000E13.6)') 
     &(max(rootdp(i,ft),0.0),i=1,ftrls(ft))
          ENDDO

          WRITE(89,'(F7.3,F9.3)') lat,lon
          WRITE(89,'(100E16.8)') (lai(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (npp(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (nps(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (npl(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (evp(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (dof(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (slc(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (rlc(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (sln(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (rln(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (nppstore(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (nppstorx(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (nppstor2(ft),ft=1,nft)
          WRITE(89,'(100I4)   ') (bb(ft),ft=1,nft)
          WRITE(89,'(100I4)   ') (bbgs(ft),ft=1,nft)
          WRITE(89,'(100I4)   ') (dsbb(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (maxlai(ft),ft=1,nft)
          WRITE(89,'(100E16.8)') (tmem(ft),ft=1,200)
          WRITE(89,*) maprc,maswr,masoilcn,soilp,
     &matmp,matmp_max,matmp_min,map_days,mahum,masoilc,masoilw  
          DO ft=1,nft
            WRITE(89,'(100E16.8)') (sm_trig(i,ft),i=1,30)
          ENDDO
          DO ft=3,nft
            WRITE(89,'(360E16.8)') (sumadp(i,ft),i=1,360)
          ENDDO
          WRITE(89,'(100E16.8)') (stemfr(ft),ft=1,nft)
          WRITE(89,*) ce_t 
          WRITE(89,*) ce_rh       
          DO ft=3,nft
            WRITE(89,*) (ce_ci(:,i,ft),      i=1,12)
            WRITE(89,*) (ce_ga(:,i,ft),      i=1,12)       
            WRITE(89,*) (ce_light(:,i,ft),   i=1,12)
            WRITE(89,*) (ce_maxlight(:,i,ft),i=1,12)
            WRITE(89,*) (vcmax(i,ft),i=1,12)
            WRITE(89,*) (jmax(i,ft) ,i=1,12)
            WRITE(89,*) (pnlc(i,ft) ,i=1,12)
            WRITE(89,*) (enzs(i,ft) ,i=1,12)
          ENDDO
        ENDIF

        WRITE(80,'(F7.3,F9.3)') lat,lon
        WRITE(80,'(100E16.8)') (s1(ft),ft=1,nft)
        WRITE(80,'(100E16.8)') (s2(ft),ft=1,nft)
        WRITE(80,'(100E16.8)') (s3(ft),ft=1,nft)
        WRITE(80,'(100E16.8)') (s4(ft),ft=1,nft)
        WRITE(80,'(100E16.8)') (sn(ft),ft=1,nft)
        WRITE(80,'(100E16.8)') (lsn(ft),ft=1,nft)
        WRITE(80,'(1008E16.8)') ic0
        WRITE(80,'(1008E16.8)') in0
        WRITE(80,'(1003E16.8)') iminn
        WRITE(80,'(100E16.8)') dslc
        WRITE(80,'(100E16.8)') drlc
        WRITE(80,'(100E16.8)') dsln
        WRITE(80,'(100E16.8)') drln
        WRITE(80,'(100E16.8)') avlai
        WRITE(80,'(100E16.8)') avnpp

      ELSE
C        WRITE(80,'(A4,F7.3,F9.3)') 'BAD',lat,lon        !PCM
C        WRITE(80,*) '                  clm stt ssc blk wfs dep lus'
C        WRITE(80,'(f7.3,f9.3,1x,20L4)') lat,lon,l_clim,l_stats, !PCM
C     &l_soil(1),l_soil(3),l_soil(5),l_soil(8),l_lu              !PCM

        WRITE(11,*) '                  clm stt ssc blk wfs dep lus'
        WRITE(11,'(f7.3,f9.3,1x,20L4)') lat,lon,l_clim,l_stats,
     &l_soil(1),l_soil(3),l_soil(5),l_soil(8),l_lu

*----------------------------------------------------------------------*
* End of climate data exists 'if' statement.                           *
*----------------------------------------------------------------------*
      ENDIF

*----------------------------------------------------------------------*
* End of site loop                                                     *
*----------------------------------------------------------------------*
      ENDDO

*----------------------------------------------------------------------*
* Open file to record version number, command line, input file and     *
* parameter file.                                                      *
*----------------------------------------------------------------------*
      st1 = stoutput
      OPEN(13,FILE=st1(1:blank(st1))//'/simulation.dat')
      CALL GETARG(0,buff1)
      DO j=1,80
        st1(j:j) = buff1(j:j)
      ENDDO
      st4 = stver
      CALL STRIPB(st4)
      WRITE(13,'(A)') st4(1:28)
      CLOSE(98)
      WRITE(13,*)

      WRITE(13,'(''************************************************'')')
      WRITE(13,'(''* Command line arguments                       *'')')
      WRITE(13,'(''************************************************'')')

      DO i=0,narg
        CALL GETARG(i,buff1)
        DO j=1,80
          st1(j:j) = buff1(j:j)
        ENDDO
        WRITE(13,'(A)') st1(1:blank(st1))
      ENDDO
      WRITE(13,*)

      WRITE(13,'(''************************************************'')')
      WRITE(13,'(''* Parameter file (param.dat)                   *'')')
      WRITE(13,'(''************************************************'')')

      OPEN(98,file='param.dat',STATUS='OLD')
80    READ(98,'(A)',end=90) st1
      WRITE(13,'(A)') st1(1:last_blank(st1))
      GOTO 80
90    CONTINUE 
      CLOSE(98)
      WRITE(13,*)

      WRITE(13,'(''************************************************'')')
      WRITE(13,'(''* Input file                                   *'')')
      WRITE(13,'(''************************************************'')')

      CALL GETARG(1,buff1)
*      buff1='..\..\lindsay2\lind1.dat'

      OPEN(98,file=buff1,STATUS='OLD')
65    READ(98,'(A)',end=70) st1
      WRITE(13,'(A)') st1(1:last_blank(st1))
      GOTO 65
70    CONTINUE 
      CLOSE(98)

      WRITE(13,'(''************************************************'')')
*----------------------------------------------------------------------*

      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      DO i=21,69
        CLOSE(i)
      ENDDO

      iofn = iofngft
      IF (outyears2.GT.0) THEN
        DO ft=1,nft
          IF (out_cov) THEN
            iofn = iofn + 1
            CLOSE(iofn)
          ENDIF
          IF (out_bio) THEN
            iofn = iofn + 1
            CLOSE(iofn)
          ENDIF
          IF (out_bud) THEN
            iofn = iofn + 1
            CLOSE(iofn)
          ENDIF
          IF (out_sen) THEN
            iofn = iofn + 1
            CLOSE(iofn)
          ENDIF
        ENDDO
      ENDIF

      DO i=70,89
        CLOSE(i)
      ENDDO
      DO i=91,93
        CLOSE(i)
      ENDDO
      DO i=98,99
        CLOSE(i)
      ENDDO

      IF (snp_no.GT.0) THEN
        DO i=1,snp_no
          fno = 100+(i-1)*4
          CLOSE(fno+1)
          CLOSE(fno+1)
          CLOSE(fno+3)
          CLOSE(fno+4)
        ENDDO
      ENDIF

      print*, ''
      print*, 'SDGVM End'
      STOP
      END

*----------------------------------------------------------------------*
*                          File I/O table                              *
*----------------------------------------------------------------------*
*  1:
*  2:
*  3:
*  4:
*  5:
*  6:
*  7:
*  8:
*  9:
* 10:
* 11: diag: Diagnostics file
* 12: site: Site information
* 13: simulation.dat: Record of the command, version number, input file
* and parameter file.
* 14: sites: List of land sites from the land sea mask.
* 15:
* 16:
* 17:
* 18:
* 19:
* 20: leafc: leaf Carbon (g C m-2)
* 21: lai: LAI
* 22: npp: Net Primary Productivity (g C /m^2/yr)
* 23: scn: Soil carbon (g/m^2)
* 24: snn: Soil nitrogen (g/m^2)
* 25: nep: Net Ecosystem Productivity (g C /m^2/yr)
* 26: swc: Average soil water content (mm)
* 27: biot: Total biomass (g C /m^2)
* 28: bioind: Dominant ft in terms of biomass
* 29: covind: Dominant ft in terms of cover
* 30: dof: Average number of days that veg is without leaves (days/yr)
* 31: rof: Runoff (mm/yr)
* 32: fcn: Burnt carbon (g C /m^2/yr)
* 33: nppstore: Stored npp (g C /m^2/yr)
* 34: stembio: Stem biomass (g C /m^2)
* 35: rootbio: Root biomass (g C /m^2)
* 36: leafper: Percentage of npp going to leaf production/maintenance
* 37: stemper: Percentage of npp going to stem production/maintenance
* 38: rootper: Percentage of npp going to root production/maintenance
* 39: sresp: Soil respiration (g C /m^2/yr)
* 40: evt: Evapotranspiration (mm/yr)
* 41: gpp: Gross Primary Productivity (g C /m^2/yr)
* 42: lch: Leached soil carbon (g C /m^2/yr)
* 43: prc: Yearly precipitation (mm/year)
* 44: nbp: Net Biosphere production (g C /m^2/yr)
* 45: trn: Transpiration (mm/yr)
* 46: fab: fraction of area burnt (0-1)
* 47: tmp: Average yearly temperature (degrees C)
* 48: hum: Average yearly temperature (degrees C)
* 49: full
* 50: *
* 51: *
* 52:  
* 53:
* 54:
* 55:
* 56:
* 57:
* 58:
* 59:
* 60:
* 61:
* 62:
* 63:
* 64:
* 65:
* 66:
* 67:
* 68: full
* 69: full
* 70: init: Input initialisation file
* 71: initbio: Input biomass state
* 72: initcov: Input cover state
* 73: initppm: Input plants per metre squared
* 74: inithgt: Input height
* 75: initwdt: Input width of stems
* 76: initleaf: Input leaf biomass by age in days
* 77: initstem: Input stem biomass by age in days
* 78: initroot: Input root biomass by age in days
* 79: initmisc: Input miscelaneous state parameters
* 80: init: Output initialisation file
* 81: initbio: Output biomass state
* 82: initcov: Output cover state
* 83: initppm: Output plants per metre squared
* 84: inithgt: Output height
* 85: initwdt: Output width of stems
* 86: initleaf: Output leaf biomass by age in days
* 87: initstem: Output stem biomass by age in days
* 88: initroot: Output root biomass by age in days
* 89: initmisc: Output miscelaneous state parameters
* 90:
* 91: Used to open temporary input files. 
* 92: Used to open temporary input files.
* 93: Used to open temporary input files.
* 94:
* 95:
* 96:
* 97:
* 98: Used to open temporary input files.
* 99: Used to open temporary input files.
*101-200: Snapshots of the system state.
*201-400: cov_'tag': Fractional coverage per ft (0-1)
*         bio_'tag': Biomass per ft (g C /m^2)
*         bud_'tag': Biomass per ft (g C /m^2)
*         sen_'tag': Biomass per ft (g C /m^2)
*401-: Output files for monthly/daily/pixel/subpixel variables specified
*in the input file.
*----------------------------------------------------------------------*


