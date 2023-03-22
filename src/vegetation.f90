module mod_vegetation
   !=================================================================================
   !  main subroutines  :  canopy,  respiration, plantgrowth 
   !             canopy => yrday,   xlayers,     Tsoil_simu
   !            xlayers => Radiso,  goudriaan,   agsean_day,  agsean_ngt
   !         agsean_day => photosyn
   !           photosyn => ciandA
   ! functions:
   !           sinbet, esat, VJtemp, fJQres, EnzK
   !=================================================================================
   use mod_data
   use mod_soil
   implicit none
   integer doy, hour
   
   real Gaussx(5), Gaussw(5), Gaussw_cum(5)                          ! Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
   real coszen, fbeam, Radabv(2)
   real Acan1, Acan2, Ecan1, Ecan2
   real extKb, extkd
   real flai, Qabs(3,2), emair, Rnstar(2), grdn
   real reff(3,2), kpr(3,2), scatt(2)
   real windUx, Vcmxx, eJmxx
   real Aleaf(2), Eleaf(2), Hleaf(2), Tleaf(2), co2ci(2)
   real gbleaf(2), gsleaf(2)
   real Tlk1, Tlk2
   real weighJ, gbHu, Tlk, Dleaf, co2cs, Qapar
   real gbc, gsc0, varQc, weighR
   real Aleafx, Gscx
   ! ---------------------------------------------------------------------------------
   data Gaussx/0.0469101, 0.2307534, 0.5, 0.7692465, 0.9530899/        ! 5-point
   data Gaussw/0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635/
   data Gaussw_cum/0.11846, 0.35777, 0.64222, 0.88153, 1.0/
! ------------------------------------------------------------------------------------
contains
   subroutine canopy() 
      implicit none
      real Ecanop
      doy  = iday
      hour = ihour+1

      call yrday()                                            ! calculate beam fraction in incoming solar radiation
      coszen = sinbet()                                       ! cos zenith angle of sun
      if (wind .lt. 0.01) wind = 0.01                         ! set windspeed to the minimum speed to avoid zero Gb
      if (topfws .gt. 0.5) then                               ! calculate soil albedo for NIR as a function of soil water (Garratt pp292)
         rhoS(2) = 0.18
      else
         rhoS(2) = 0.52 - 0.68*topfws
      end if
      ! assign plant biomass and leaf area index at time t
      ! assume leaf biomass = root biomass
      FLAIT     = LAI
      ! eairP     = esat(Tair) - Dair              !air water vapour pressure
      radabv(1) = 0.5*radsol                    !(1) - solar radn
      radabv(2) = 0.5*radsol                    !(2) - NIR
      ! call multilayer model of Leuning - uses Gaussian integration but radiation scheme
      call xlayers() 
      ! if (itest .gt.  847) then
      !    iitest = iitest + 1
      !    if (mod(iitest, 252) .eq. 0)then
      !       ! write(*,*)"test Acan1:", Vcmxx, Tlk, -0.0089*Vcmxx*exp(0.069*(Tlk - 293.2)), Tair, Aleaf, Aleafx
      !       write(*,*)"test Vcmxx: ", Vcmx0, scalex, extkn,flai, Vcmxx
      !    endif 
      ! endif
      ! run Tsoil simulation
      if (do_soilphy) call Tsoil_simu()                  ! *** ..int   added 'testout,Rsoilab1,Rsoilab2,QLleaf,QLair,raero,do_soilphy,G,Esoil,Hsoil'
      ! summary the results of canopy
      Acanop = Acan1 + Acan2                                
      Ecanop = Ecan1 + Ecan2
      gpp    = Acanop*3600.0*12.0                           ! every hour, g C m-2 h-1
      transp = AMAX1(Ecanop*3600.0/(1.0e6*(2.501 - 0.00236*Tair)), 0.) ! mm H2O /hour
      evap   = AMAX1(Esoil *3600.0/(1.0e6*(2.501 - 0.00236*Tair)), 0.)
      ! write(*,*)"test_new_evap: ", Esoil, Tair
      ! write(*,*)"GPP: ", Acan1, Acan2
      return
   end subroutine canopy

   ! -------- autotrophic respiration -----------------------------------------------------------------
   subroutine respiration()
      ! calculate plant and soil respiration by the following equation:
      ! RD=BM*Rd*Q10**((T-25)/10) (Sun et al. 2005. Acta Ecologica Sinica)
      implicit none
      real conv 
      conv = 3600.*12./1000000.                  ! converter from "umol C /m2/s" to "gC/m2/hour"
      ! SNRauto =exp(-(CN(1)-CN0(1))/CN0(1))
      if (LAI .gt. LAIMIN) then
         RmLeaf = Rl0*SNRauto*bmleaf*0.48*SLA*0.1*Q10**((Tair - 10.)/10.)*fnsc*conv
         RmStem = Rs0*SNRauto*StemSap*0.001*Q10**((Tair - 25.)/10.)*fnsc*conv
         RmRoot = Rr0*SNRauto*RootSap*0.001*Q10**((Tair - 25.)/10.)*fnsc*conv
      else
         RmLeaf = 0.3*GPP
         RmStem = 0.3*GPP
         RmRoot = 0.4*GPP
      end if
      Rmain = Rmleaf + Rmstem + Rmroot
      if (Rmain > 0.0015*NSC) then             ! If Total autotropic respiration greater than 0.15% of Nonstructure Carbon, rescale.
         Rmleaf = Rmleaf/Rmain*0.0015*NSC
         Rmstem = Rmstem/Rmain*0.0015*NSC
         Rmroot = Rmstem/Rmain*0.0015*NSC
         Rmain  = Rmleaf + Rmstem + Rmroot
      end if
      return
   end subroutine respiration

   ! ------------- plant growth model --------------------------------------------------------------------
   subroutine plantgrowth()
      implicit none  ! Jian: RS0 is different from Rs0 from data type
      ! real store
      real GrowthP, GrowthL, GrowthR, GrowthS
      real hmax, hl0, CNP0
      REAL LAIMAX0, la0, GPmax, acP, c1, c2
      real bmL, bmR, bmP, bmS !, StemSap, RootSap
      ! real Rgrowth, Rgroot, Rgleaf, Rgstem
      real, save :: addaccu = 0, GrowthLaccu = 0, GrowthSaccu = 0, GrowthRaccu = 0
      real St, Ss, Sn, SL_rs, SR_rs, Slai, phiN            ! scalars
      real RS, RSw, RS_0
      real gamma_W, gamma_Wmax, gamma_T, gamma_Tmax, gamma_N
      real beta_T, Tcold, Twarm, Topt
      real bW, bT, W
      real L_add, NL_fall, NL_add
      real alpha_St

      Twarm = 35.0
      Tcold = 5.0
      ! Tcold=0.0       ! For SPRUCE
      Topt = 30.
      phiN = 0.33

      bmL = bmleaf*0.48   ! Carbon
      bmR = bmRoot*0.48
      bmS = bmStem*0.48

      if (bmL .lt. NSC/0.333) bmL = NSC/0.333
      if (bmR .lt. NSC/0.333) bmR = NSC/0.333
      if (bmS .lt. NSC/0.334) bmS = NSC/0.334
      StemSap = SapS*bmS  ! Weng 12/05/2008
      RootSap = SapR*bmR
      if (StemSap .lt. 0.001) StemSap = 0.001
      if (RootSap .lt. 0.001) RootSap = 0.001

      bmP = bmL + bmR + bmS                                        ! Plant C biomass
      acP = bmL + StemSap + bmS                                        ! Plant available sapwood C
      CNp0 = bmP/(bmL/CN0(1) + bmR/CN0(3) + bmS/CN0(2))                ! Plant CN ratio

      hmax = 24.19   ! m
      hl0 = 0.00019  ! m2/kg C
      ! LAIMAX0 = 6.
      LAIMAX0=8.
      la0 = 0.2
      ht = hmax*(1.-exp(-hl0*bmP))                                ! Scaling plant C biomass to height
      LAIMAX = AMAX1(LAIMAX0*(1.-exp(-la0*ht)), LAIMIN + 0.1)  ! Scaling plant height to maximum LAI

      ! Phenology
      if ((GDD5 .gt. gddonset) .and. onset .eq. 0 .and. storage .gt. stor_use) then
         onset = 1
      end if
      if ((onset .eq. 1) .and. (storage .gt. stor_use)) then
         if (LAI .lt. LAIMAX) add = stor_use!Amin1(stor_use, 0.004*NSN*CNp0)
         ! if(LAI.lt.LAIMAX)add=stor_use/20.0
         storage = storage - add
         ! add = 0
      else
         add = 0.0
         onset = 0
      end if
      if (accumulation .lt. (NSCmax + 0.001*RootSap)) then
         store = AMAX1(0., 0.0001*NSC)                        ! 0.5% of nonstructure carbon is stored
      else
         store = 0.0
      end if
      accumulation = accumulation + store

      ! Scalars for plant growth
      ! Sps=Amin1(1.0,3.33*AMAX1(0.0,1.0 - fnsc))
      Sps = Sps*(1.-exp(-phiN*NSN))                                                        ! Sps is not assigned previous, something is wrong. -JJJJJJJJJJJJJJJJJJJJJ
      Ss = AMIN1(1.0, 2.*fnsc)
      RS_0 = 1.0
      RS = bmR/bmL
      SL_rs = RS/(RS + RS_0*(2.-W))
      SR_rs = (RS_0*(2.-W))/(RS + RS_0*(2.-W))
      Slai = amin1(1.0, 2.333*(LAIMAX - LAI)/(LAIMAX - LAIMIN))
      St = AMAX1(0.0, 1.0 - exp(-(Tair - gddonset/10.)/5.0))  !0.5 !
      ! Sw=AMAX1(0.333, 0.333+omega)
      Sw = AMIN1(0.5, AMAX1(0.333, 0.333 + omega))
      W = AMIN1(1.0, 3.333*omega)

      ! Plant growth and allocation, based on LM3V
      GPmax   = (GLmax*bmL + GSmax*StemSap + GRmax*bmR)                        !/acP
      ! Jian: limit the Growth > 0
      GrowthP = Amax1(AMIN1(GPmax*fnsc*St*(1.-exp(-NSN)),0.004*NSC,0.004*NSN*CNp0), 0.) !    ! Jian: GrowthP is used to calculate NPP, which means the minimum of group plant needs and 0.004*NSC, 0.004*NSN*CNp0
      GrowthL = MAX(0.0, GrowthP*0.5)                                          ! updated when QC leaf and wood changed due to the change of plot area for tree biomass
      GrowthR = Amax1(MIN(GrowthP*0.4, MAX(0.0, 0.75/Sw*bmL - bmR)), 0.)                  ! *c1/(1.+c1+c2)
      GrowthS = MAX(0.0, GrowthP - (GrowthL + GrowthR))                        ! *c2/(1.+c1+c2)

      NPP         = GrowthL + GrowthR + GrowthS + add       ! Modified by Jiang Jiang 2015/10/13
      addaccu     = addaccu + add
      GrowthLaccu = GrowthLaccu + GrowthL
      GrowthRaccu = GrowthRaccu + GrowthR
      GrowthSaccu = GrowthSaccu + GrowthS
      ! test_gpp = (/GrowthL,GrowthR,GrowthS,add,0.,0.,0.,0.,0./)
      test_gpp = (/GPmax*fnsc*St*(1.-exp(-NSN)),0.004*NSC,0.004*NSN*CNp0,GPmax,fnsc,nsn,CNP0,GLmax,bmL/)
      ! test_gpp = (/GPmax*fnsc*St*(1.-exp(-NSN)),0.004*NSC,0.004*NSN*CNp0,GPmax,fnsc,QC(1),bmleaf,bmleaf*0.48,bmL/)
      
      if (NPP .eq. 0.0) then
         alpha_L = 0.333
         alpha_W = 0.333
         alpha_R = 0.333
      else
         alpha_L = (GrowthL + add)/NPP
         alpha_W = GrowthS/NPP
         alpha_R = GrowthR/NPP
      end if

      ! Carbon cost for growth
      ! Rgrowth,Rgroot,Rgleaf,Rgstem, 0.5 is from IBIS and Amthor, 1984
      Rgleaf  = 0.5*GrowthL
      Rgstem  = 0.5*GrowthS
      Rgroot  = 0.5*GrowthR
      Rgrowth = Rgleaf + Rgstem + Rgroot
      ! Leaf litter
      gamma_Wmax = 0.12/24. ! maxmum leaf fall rate per hour
      gamma_Tmax = 0.12/24.
      bW = 4.0
      bT = 2.0
      if (Tair .gt. (Tcold + 10.)) then
         beta_T = 1.
      else
         if (Tair .gt. Tcold) beta_T = (Tair - Tcold)/10.
         if (Tair .LE. Tcold) beta_T = 0.0
      end if

      if (TauC(1) < 8760.) then
         gamma_W = (1.-W)**bW*gamma_Wmax
         gamma_T = (1.-beta_T)**bT*gamma_Tmax
      else
         gamma_W = 0.
         gamma_T = 0.
      end if

      gamma_N = 1.0/TauC(1)*Sw      ! Modify by Jiang Jiang 2015/10/20
      if (LAI < LAIMIN) then
         gamma_W = 0.
         gamma_T = 0.
         gamma_N = 0.
      end if
      L_fall = bmleaf*0.48*gamma_N     ! L_fall=bmleaf*0.48*AMIN1((gamma_T+gamma_N),0.99)
!       if (itest .gt. 850)then
!          iitest = iitest + 1
!          if (mod(iitest, 36) .eq. 0)then
!       !   Write(*,*)"nppFIX: ",GrowthP,GrowthL,add,NPP,GPmax*fnsc*St*(1.-exp(-NSN)),0.004*NSC,0.004*NSN*CNp0
! ! Write(*,*)"GrowthP: ",NPP,fnsc, St,(1.-exp(-NSN)),NSN,NSC,NSCmax,NSCmin,1.0-exp(-(Tair-gddonset/10.)/5.0),Tair,GDDonset
! write(*,*)"",GDD5,gddonset,onset,storage,stor_use
!         ENDIF
!         if (NSN<0.00003)STOP
!         endif
      return
   end
   
   ! --------- yrday --------------------------------------------------------------
   subroutine yrday()
      real pidiv, slatx, sindec, cosdec
      real a, b, sinbet0, solext, tmprat, tmpR, tmpK, fdiff
      pidiv = pi/180.0
      slatx = lat*pidiv
      sindec = -sin(23.4*pidiv)*cos(2.0*pi*(doy + 10.0)/365.0)
      cosdec = sqrt(1.-sindec*sindec)
      a = sin(slatx)*sindec
      b = cos(slatx)*cosdec
      sinbet0 = a + b*cos(2*pi*(hour - 12.)/24.)
      solext = 1370.0*(1.0 + 0.033*cos(2.0*pi*(doy - 10.)/365.0))*sinbet0
      tmprat = radsol/solext
      tmpR = 0.847 - 1.61*sinbet0 + 1.04*sinbet0*sinbet0
      tmpK = (1.47 - tmpR)/1.66
      if (tmprat .le. 0.22) fdiff = 1.0
      if (tmprat .gt. 0.22 .and. tmprat .le. 0.35) then
         fdiff = 1.0 - 6.4*(tmprat - 0.22)*(tmprat - 0.22)
      end if
      if (tmprat .gt. 0.35 .and. tmprat .le. tmpK) then
         fdiff = 1.47 - 1.66*tmprat
      end if
      if (tmprat .ge. tmpK) then
         fdiff = tmpR
      end if
      fbeam = 1.0 - fdiff
      if (fbeam .lt. 0.0) fbeam = 0.0
      return
   end subroutine yrday

   !-----------  used by canopy module ------------------------------------------------------
   subroutine xlayers() ! added from soil thermal ..int
      ! the multi-layered canopy model developed by
      ! Ray Leuning with the new radiative transfer scheme
      ! implemented by Y.P. Wang (from Sellers 1986)
      ! 12/Sept/96 (YPW) correction for mean surface temperature of sunlit
      ! and shaded leaves
      ! Tleaf,i=sum{Tleaf,i(n)*fslt*Gaussw(n)}/sum{fslt*Gaussw(n)}
      ! -------------------------------------------------------------------------------------
      real Rnst1 , Rnst2 , Qcan1 , Qcan2  !net rad, sunlit, vis rad
      real Rcan1 , Rcan2 , Hcan1 , Hcan2  !NIR rad !Sens heat
      real Gbwc1 , Gbwc2 , Gswc1 , Gswc2  !Boundary layer conductance; Canopy conductance
      real Tleaf1 , Tleaf2                           !Leaf Temp
      real xphi1, xphi2, funG ! more places?
      real pi180, cozen15, cozen45, cozen75, xK15, xK45, xK75
      real transd, extkn
      integer nw, ng
      real rhoc(3, 2)       !Goudriaan
      real rhoch, rhoc15, rhoc45, rhoc75
      real scalex, fslt, fshd
      real RnStL(5), QcanL(5), RcanL(5), AcanL(5), EcanL(5), HcanL(5)
      real layer1(5), layer2(5)
      real FLAIT1
      
      ! soil water conditions
      WILTPT = wsmin/100.
      FILDCP = wsmax/100.
      ! reset the vairables
      Rnst1=0.0        !net rad, sunlit
      Rnst2=0.0        !net rad, shaded
      Qcan1=0.0        !vis rad
      Qcan2=0.0
      Rcan1=0.0        !NIR rad
      Rcan2=0.0
      Acan1=0.0        !CO2
      Acan2=0.0
      Ecan1=0.0        !Evap
      Ecan2=0.0
      Hcan1=0.0        !Sens heat
      Hcan2=0.0
      Gbwc1=0.0        !Boundary layer conductance
      Gbwc2=0.0
      Gswc1=0.0        !Canopy conductance
      Gswc2=0.0
      Tleaf1=0.0       !Leaf Temp
      Tleaf2=0.0  
      ! aerodynamic resistance
      raero = 50./wind
      ! Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
      xphi1 = 0.5 - 0.633*xfang - 0.33*xfang*xfang
      xphi2 = 0.877*(1.0 - 2.0*xphi1)
      funG = xphi1 + xphi2*coszen                             !G-function: Projection of unit leaf area in direction of beam
      if (coszen .gt. 0) then                                     !check if day or night
         extKb = funG/coszen                                    !beam extinction coeff - black leaves
      else
         extKb = 100.
      end if
      ! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
      ! Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
      pi180   = 3.1416/180.
      cozen15 = cos(pi180*15)
      cozen45 = cos(pi180*45)
      cozen75 = cos(pi180*75)
      xK15    = xphi1/cozen15 + xphi2
      xK45    = xphi1/cozen45 + xphi2
      xK75    = xphi1/cozen75 + xphi2
      transd  = 0.308*exp(-xK15*FLAIT) + 0.514*exp(-xK45*FLAIT) +     &
                &       0.178*exp(-xK75*FLAIT)
      extkd   = (-1./FLAIT)*alog(transd)
      extkn   = extkd                        !N distribution coeff

      ! canopy reflection coefficients (Array indices: first;  1=VIS,  2=NIR
      !   second; 1=beam, 2=diffuse
      do nw = 1, 2  ! nw:1=VIS, 2=NIR
         scatt(nw)   = tauL(nw) + rhoL(nw)                                               ! scattering coeff
         if ((1.-scatt(nw)) < 0.0) scatt(nw) = 0.9999                                    ! Weng 10/31/2008
         kpr(nw, 1)  = extKb*sqrt(1.-scatt(nw))                                          ! modified k beam scattered (6.20)
         kpr(nw, 2)  = extkd*sqrt(1.-scatt(nw))                                          ! modified k diffuse (6.20)
         rhoch       = (1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))                   ! canopy reflection black horizontal leaves (6.19)
         rhoc15      = 2.*xK15*rhoch/(xK15 + extkd)                                      ! canopy reflection (6.21) diffuse
         rhoc45      = 2.*xK45*rhoch/(xK45 + extkd)
         rhoc75      = 2.*xK75*rhoch/(xK75 + extkd)
         rhoc(nw, 2) = 0.308*rhoc15 + 0.514*rhoc45 + 0.178*rhoc75
         rhoc(nw, 1) = 2.*extKb/(extKb + extkd)*rhoch                                    ! canopy reflection (6.21) beam
         reff(nw, 1) = rhoc(nw, 1) + (rhoS(nw) - rhoc(nw, 1))*exp(-2.*kpr(nw, 1)*FLAIT)  ! effective canopy-soil reflection coeff - beam (6.27)
         reff(nw, 2) = rhoc(nw, 2) + (rhoS(nw) - rhoc(nw, 2))*exp(-2.*kpr(nw, 2)*FLAIT)  ! effective canopy-soil reflection coeff - diffuse (6.27)
      end do
      ! isothermal net radiation & radiation conductance at canopy top - needed to calc emair
      call Radiso()           ! Jian: some parameters not initialization.
      TairK = Tair + 273.2
      Tleaf1 = 0.
      do ng = 1, 5
         flai = gaussx(ng)*FLAIT
         ! radiation absorption for visible and near infra-red
         call goudriaan()
         ! isothermal net radiation & radiation conductance at canopy top
         call Radiso()
         windUx = wind*exp(-extkU*flai)                ! windspeed at depth xi
         scalex = exp(-extkn*flai)                     ! scale Vcmx0 & Jmax0

         Vcmxx = Vcmx0*scalex                         ! Vcmx0 ---> Vcmax0
         ! write(*,*)"test Vcmxx0: ", Vcmax0, scalex, Vcmax0*scalex
         ! write(*,*)"test Vcmxx: ", Vcmx0, scalex, extkn,flai, Vcmxx
         ! test_gpp =(/Vcmx0, scalex,extkn,flai, eJmx0,SNvcmax,FLAIT,transd,0./)

         eJmxx = eJmx0*scalex
         if (radabv(1) .ge. 10.0) then                          !check solar Radiation > 10 W/m2
            ! leaf stomata-photosynthesis-transpiration model - daytime
            call agsean_day()
         else
            call agsean_ngt()
         end if
         fslt      = exp(-extKb*flai)                        !fraction of sunlit leaves
         fshd      = 1.0 - fslt                                !fraction of shaded leaves
         Rnst1     = Rnst1 + fslt*Rnstar(1)*Gaussw(ng)*FLAIT  !Isothermal net rad`
         Rnst2     = Rnst2 + fshd*Rnstar(2)*Gaussw(ng)*FLAIT
         RnstL(ng) = Rnst1 + Rnst2

         Qcan1     = Qcan1 + fslt*Qabs(1, 1)*Gaussw(ng)*FLAIT  !visible
         Qcan2     = Qcan2 + fshd*Qabs(1, 2)*Gaussw(ng)*FLAIT
         QcanL(ng) = Qcan1 + Qcan2

         Rcan1     = Rcan1 + fslt*Qabs(2, 1)*Gaussw(ng)*FLAIT  !NIR
         Rcan2     = Rcan2 + fshd*Qabs(2, 2)*Gaussw(ng)*FLAIT
         RcanL(ng) = Rcan1 + Rcan2

         if (Aleaf(1) .lt. 0.0) Aleaf(1) = 0.0      !Weng 2/16/2006
         if (Aleaf(2) .lt. 0.0) Aleaf(2) = 0.0      !Weng 2/16/2006

         Acan1      = Acan1 + fslt*Aleaf(1)*Gaussw(ng)*FLAIT*stom_n    !amphi/hypostomatous
         Acan2      = Acan2 + fshd*Aleaf(2)*Gaussw(ng)*FLAIT*stom_n
         AcanL(ng)  = Acan1 + Acan2

         layer1(ng) = Aleaf(1)
         layer2(ng) = Aleaf(2)

         Ecan1      = Ecan1 + fslt*Eleaf(1)*Gaussw(ng)*FLAIT
         Ecan2      = Ecan2 + fshd*Eleaf(2)*Gaussw(ng)*FLAIT
         EcanL(ng)  = Ecan1 + Ecan2

         Hcan1     = Hcan1 + fslt*Hleaf(1)*Gaussw(ng)*FLAIT
         Hcan2     = Hcan2 + fshd*Hleaf(2)*Gaussw(ng)*FLAIT
         HcanL(ng) = Hcan1 + Hcan2

         Gbwc1     = Gbwc1 + fslt*gbleaf(1)*Gaussw(ng)*FLAIT*stom_n
         Gbwc2     = Gbwc2 + fshd*gbleaf(2)*Gaussw(ng)*FLAIT*stom_n

         Gswc1     = Gswc1 + fslt*gsleaf(1)*Gaussw(ng)*FLAIT*stom_n
         Gswc2     = Gswc2 + fshd*gsleaf(2)*Gaussw(ng)*FLAIT*stom_n

         Tleaf1    = Tleaf1 + fslt*Tleaf(1)*Gaussw(ng)*FLAIT
         Tleaf2    = Tleaf2 + fshd*Tleaf(2)*Gaussw(ng)*FLAIT
         ! if (itest .eq.  853) then
         ! write(*,*)"test Acan1:", Vcmxx, Tlk, -0.0089*Vcmxx*exp(0.069*(Tlk - 293.2)), Tair
         ! write(*,*)"test Teaf1: ",iforcing,gpp, Acan1, Acan2, fshd, fslt, Aleaf, Gaussw(ng), flait, stom_n,Vcmxx,Tlk,extKb,flai
         ! endif
         ! test_gpp(ng,:) = (/Acan1, Acan2, fshd, fslt, Aleaf, Gaussw(ng), flait, stom_n,Vcmxx,Tlk,extKb,flai/)
      end do  ! 5 layers

      FLAIT1 = (1.0 - exp(-extKb*FLAIT))/extkb
      Tleaf1 = Tleaf1/FLAIT1
      Tleaf2 = Tleaf2/(FLAIT - FLAIT1)
      ! Soil surface energy and water fluxes
      ! Radiation absorbed by soil
      Rsoilab1 = fbeam*(1.-reff(1, 1))*exp(-kpr(1, 1)*FLAIT)        &
          &         + (1.-fbeam)*(1.-reff(1, 2))*exp(-kpr(1, 2)*FLAIT)          !visible
      Rsoilab2 = fbeam*(1.-reff(2, 1))*exp(-kpr(2, 1)*FLAIT)        &
          &         + (1.-fbeam)*(1.-reff(2, 2))*exp(-kpr(2, 2)*FLAIT)          !NIR
      Rsoilab1 = Rsoilab1*Radabv(1)
      Rsoilab2 = Rsoilab2*Radabv(2)
      Tlk1     = Tleaf1 + 273.2
      Tlk2     = Tleaf2 + 273.2
      ! temp1=-extkd*FLAIT
      QLair    = emair*sigma*(TairK**4)*exp(-extkd*FLAIT)
      QLleaf   = emleaf*sigma*(Tlk1**4)*exp(-extkb*FLAIT)           &
                 &      + emleaf*sigma*(Tlk2**4)*(1.0 - exp(-extkb*FLAIT))
      if(ISNAN(QLleaf)) then
         write(*,*)"QLleaf is NAN1111: ",QLleaf, emleaf, sigma, Tlk1, extKb, flait, Tlk2, extkd, Tleaf1, FLAIT1
         stop
      endif
      QLleaf   = QLleaf*(1.0 - exp(-extkd*FLAIT))
      QLsoil   = emsoil*sigma*(TairK**4)
      Rsoilab3 = (QLair + QLleaf)*(1.0 - rhoS(3)) - QLsoil
      if(ISNAN(QLleaf)) then
         write(*,*)"QLleaf is NAN: ", emleaf, sigma, Tlk1, extKb, flait, Tlk2, extkd
         stop
      endif
      ! Net radiation absorbed by soil
      ! the old version of net long-wave radiation absorbed by soils
      ! (with isothermal assumption)
      ! Rsoil3=(sigma*TairK**4)*(emair-emleaf)*exp(-extkd*FLAIT)         !Longwave
      ! Rsoilab3=(1-rhoS(3))*Rsoil3

      ! Total radiation absorbed by soil
      Rsoilabs = Rsoilab1 + Rsoilab2 + Rsoilab3

      ! thermodynamic parameters for air
      TairK  = Tair + 273.2
      rhocp  = cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv  = H2oLv0 - 2.365e3*Tair
      slope  = (esat(Tair + 0.1) - esat(Tair))/0.1
      psyc   = Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar = Patm/(Rconst*TairK)
      fw1    = AMIN1(AMAX1((FILDCP - wcl(1))/(FILDCP - WILTPT), 0.05), 1.0)
      Rsoil  = 30.*exp(0.2/fw1)
      rLAI   = exp(FLAIT)
      ! latent heat flux into air from soil
      ! Eleaf(ileaf)=1.0*
      ! &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
      ! &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
      Esoil  = (slope*(Rsoilabs - G) + rhocp*Dair/(raero + rLAI))/       &
               &      (slope + psyc*(rsoil/(raero + rLAI) + 1.))!*omega!*AMIN1(0.5, AMAX1(0.333, 0.333 + omega))
      ! sensible heat flux into air from soil
      Hsoil = Rsoilabs - Esoil - G
      
      ! write(*,*)"test_Esoil: ", slope, Rsoilabs, G, rhocp, Dair, raero, rLAI, psyc, Rsoil, raero
      ! write(*,*)"test_G: ", G
      ! endif
      ! if (itest .gt.  944) then
      !    iitest = iitest + 1
      !    if (mod(iitest, 252) .eq. 0)then
      !       write(*,*)"test Acan1:", Vcmxx, Tlk, -0.0089*Vcmxx*exp(0.069*(Tlk - 293.2)), Tair, Aleaf, Aleafx
      !       ! write(*,*)"test Vcmxx: ", Vcmx0, scalex, extkn,flai, Vcmxx,Vcmax0,SNvcmax

      !    endif 
      ! endif
      ! if (itest .gt.  847) then
      !    iitest = iitest + 1
      !    if (mod(iitest, 2520) .eq. 0)then
      !       ! write(*,*)"Aleafx: ",Aleafx,Gma, Bta, gsc0, X, Rd, co2Csx, gammas, Aqx
      !       ! write(*,*)"ej: ",weighJ,fJQres(eJmxT, alpha, Qapar, theta),eJmxT, alpha, Qapar, theta
      !       ! write(*,*)"eJmxT",eJmxT,VJtemp(Tlf, TminJ, TmaxJ, ToptJ, eJmxx),Tlf, TminJ, TmaxJ, ToptJ, eJmxx
      !       ! write(*,*)"eJmxx:",eJmx0,scalex,Vcmax0,SNvcmax
      !    endif 
      ! endif
      return
   end subroutine xlayers

   subroutine Radiso()
      ! Rnstar(type): type=1 for sunlit; =2 for shaded leaves (W/m2)
      ! 23 Dec 1994
      ! calculates isothermal net radiation for sunlit and shaded leaves under clear skies
      real emsky, ep8z, tau8, emcloud, Bn0, Bnxi
      TairK = Tair + 273.2
      ! thermodynamic properties of air
      rhocp   = cpair*Patm*airMa/(Rconst*TairK)   !volumetric heat capacity (J/m3/K)
      ! apparent atmospheric emissivity for clear skies (Brutsaert, 1975)
      emsky   = 0.642*(eairP/Tairk)**(1./7)       !note eair in Pa
      ! apparent emissivity from clouds (Kimball et al 1982)
      ep8z    = 0.24 + 2.98e-12*eairP*eairP*exp(3000/TairK)
      tau8    = amin1(1.0, 1.0 - ep8z*(1.4 - 0.4*ep8z))            !ensure tau8<1
      emcloud = 0.36*tau8*(1.-fbeam)*(1 - 10./TairK)**4      !10 from Tcloud = Tair-10
      ! apparent emissivity from sky plus clouds
      !      emair=emsky+emcloud
      ! 20/06/96
      emair   = emsky
      if (emair .gt. 1.0) emair = 1.0
      ! net isothermal outgoing longwave radiation per unit leaf area at canopy
      ! top & thin layer at flai (Note Rn* = Sn + Bn is used rather than Rn* = Sn - Bn in Leuning et al 1985)
      Bn0  = sigma*(TairK**4.)
      Bnxi = Bn0*extkd*(exp(-extkd*flai)*(emair - emleaf)       &
          &    + exp(-extkd*(flait - flai))*(emsoil - emleaf))
      ! isothermal net radiation per unit leaf area for thin layer of sunlit and
      ! shaded leaves
      Rnstar(1) = Qabs(1, 1) + Qabs(2, 1) + Bnxi
      Rnstar(2) = Qabs(1, 2) + Qabs(2, 2) + Bnxi
      ! radiation conductance (m/s) @ flai
      grdn = 4.*sigma*(TairK**3.)*extkd*emleaf*               &       ! corrected by Jiang Jiang 2015/9/29
          &    (exp(-extkd*flai) + exp(-extkd*(flait - flai)))       &
          &    /rhocp
      return
   end subroutine Radiso

   subroutine goudriaan()
      ! for spheric leaf angle distribution only
      ! compute within canopy radiation (PAR and near infra-red bands)
      ! using two-stream approximation (Goudriaan & vanLaar 1994)
      ! tauL: leaf transmittance
      ! rhoL: leaf reflectance
      ! rhoS: soil reflectance
      ! sfang XiL function of Ross (1975) - allows for departure from spherical LAD
      !     (-1 vertical, +1 horizontal leaves, 0 spherical)
      ! FLAI: canopy leaf area index
      ! funG: Ross' G function
      ! scatB: upscatter parameter for direct beam
      ! scatD: upscatter parameter for diffuse
      ! albedo: single scattering albedo
      ! output:
      ! Qabs(nwave,type), nwave=1 for visible; =2 for NIR,
      !                     type=1 for sunlit;   =2 for shaded (W/m2)
      real xu, xphi1, xphi2, funG, Qd0, Qb0
      integer nw
      xu = coszen                                         !cos zenith angle

      ! Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
      xphi1 = 0.5 - 0.633*xfang - 0.33*xfang*xfang
      xphi2 = 0.877*(1.0 - 2.0*xphi1)
      funG = xphi1 + xphi2*xu                             !G-function: Projection of unit leaf area in direction of beam
      if (coszen .gt. 0) then                                  !check if day or night
         extKb = funG/coszen                                   !beam extinction coeff - black leaves
      else
         extKb = 100.
      end if
      ! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
      do nw = 1, 2
         Qd0 = (1.-fbeam)*radabv(nw)                                          !diffuse incident radiation
         Qb0 = fbeam*radabv(nw)                                               !beam incident radiation
         Qabs(nw, 2) = Qd0*(kpr(nw, 2)*(1.-reff(nw, 2))*exp(-kpr(nw, 2)*FLAI)) +  & !absorbed radiation - shaded leaves, diffuse
             &      Qb0*(kpr(nw, 1)*(1.-reff(nw, 1))*exp(-kpr(nw, 1)*FLAI) -   & !beam scattered
             &      extKb*(1.-scatt(nw))*exp(-extKb*FLAI))
         Qabs(nw, 1) = Qabs(nw, 2) + extKb*Qb0*(1.-scatt(nw))                     !absorbed radiation - sunlit leaves
      end do
      return
   end subroutine goudriaan

   subroutine agsean_day()
      integer kr1, ileaf
      real Gras, gbHf, gbH, rbH, rbw, rbH_L, rrdn, Y
      real gsw, gswv, rswv
      ! thermodynamic parameters for air
      TairK = Tair + 273.2
      rhocp = cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv = H2oLv0 - 2.365e3*Tair
      slope = (esat(Tair + 0.1) - esat(Tair))/0.1
      psyc = Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar = Patm/(Rconst*TairK)
      weighJ = 1.0
      ! boundary layer conductance for heat - single sided, forced convection
      ! (Monteith 1973, P106 & notes dated 23/12/94)
      if (windUx/wleaf >= 0.0) then
         gbHu = 0.003*sqrt(windUx/wleaf)    !m/s
      else
         gbHu = 0.003 !*sqrt(-windUx/wleaf)
      end if         ! Weng 10/31/2008
      ! raero=0.0                        !aerodynamic resistance s/m
      do ileaf = 1, 2              ! loop over sunlit and shaded leaves
         ! first estimate of leaf temperature - assume air temp
         Tleaf(ileaf) = Tair
         Tlk = Tleaf(ileaf) + 273.2    !Tleaf to deg K
         ! first estimate of deficit at leaf surface - assume Da
         Dleaf = Dair                !Pa
         ! first estimate for co2cs
         co2cs = co2ca               !mol/mol
         Qapar = (4.6e-6)*Qabs(1, ileaf)
         ! ********************************************************************
         kr1 = 0                     !iteration counter for LE
         ! return point for evaporation iteration
         do ! iteration for leaf temperature
            ! single-sided boundary layer conductance - free convection (see notes 23/12/94)
            Gras   = 1.595e8*ABS(Tleaf(ileaf) - Tair)*(wleaf**3.)     !Grashof
            gbHf   = 0.5*Dheat*(Gras**0.25)/wleaf
            gbH    = gbHu + gbHf                         !m/s
            rbH    = 1./gbH                            !b/l resistance to heat transfer
            rbw    = 0.93*rbH                          !b/l resistance to water vapour
            ! Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
            rbH_L  = rbH*stom_n/2.                   !final b/l resistance for heat
            rrdn   = 1./grdn
            Y      = 1./(1.+(rbH_L + raero)/rrdn)
            ! boundary layer conductance for CO2 - single side only (mol/m2/s)
            gbc    = Cmolar*gbH/1.32            !mol/m2/s
            gsc0   = gsw0/1.57                 !convert conductance for H2O to that for CO2
            varQc  = 0.0
            weighR = 1.0
            ! -------------------------------------------
            ! write(*,*)"test ileaf: ", ileaf
            ! write(*,*)"before photosyn: ", Gbc, fwsoil, Vcmxx
            ! write(*,*)"before photosyn2: ", Cmolar, gbH, gbHu, gbHf, Dheat,Gras, wleaf
            call photosyn()  !outputs
            ! choose smaller of Ac, Aq
               
            Aleaf(ileaf) = Aleafx      !0.7 Weng 3/22/2006          !mol CO2/m2/s
            ! calculate new values for gsc, cs (Lohammer model)
            co2cs = co2ca - Aleaf(ileaf)/gbc
            co2Ci(ileaf) = co2cs - Aleaf(ileaf)/gscx
            ! scale variables
            gsw  = gscx*1.56       !gsw in mol/m2/s, oreginal:gsw=gscx*1.56,Weng20090226
            gswv = gsw/Cmolar                           !gsw in m/s
            rswv = 1./gswv
            ! calculate evap'n using combination equation with current estimate of gsw
            Eleaf(ileaf) = 1.0*(slope*Y*Rnstar(ileaf) + rhocp*Dair/(rbH_L + raero))/    &   !2* Weng 0215
                           & (slope*Y + psyc*(rswv + rbw + raero)/(rbH_L + raero))!*omega!*AMIN1(0.5, AMAX1(0.333, 0.333 + omega))!*(0.08+(1.-0.08)* Amin1(1.0, 0.3*omega))
            ! calculate sensible heat flux
            Hleaf(ileaf) = Y*(Rnstar(ileaf) - Eleaf(ileaf))
            ! calculate new leaf temperature (K)
            Tlk1 = 273.2 + Tair + Hleaf(ileaf)*(rbH/2.+raero)/rhocp
            ! calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
            Dleaf = psyc*Eleaf(ileaf)/(rhocp*gswv)
            gbleaf(ileaf) = gbc*1.32*1.075
            gsleaf(ileaf) = gsw
            ! compare current and previous leaf temperatures
            if (abs(Tlk1 - Tlk) .le. 0.1) exit ! original is 0.05 C Weng 10/31/2008
            ! update leaf temperature  ! leaf temperature calculation has many problems! Weng 10/31/2008
            Tlk = Tlk1
            Tleaf(ileaf) = Tlk1 - 273.2
            kr1 = kr1 + 1
            if (kr1 > 500) then
               Tlk = TairK
               exit
            end if
            if (Tlk < 200.) then
               Tlk = TairK
               exit
            end if                     ! Weng 10/31/2008
            ! goto 100                          !solution not found yet
         end do
         ! 10  continue
      end do
      return
   end subroutine agsean_day

   ! ****************************************************************************
   subroutine agsean_ngt()
      ! implicit real (a-z)
      integer kr1, ileaf
      real Gras, gbHf, gbH, rbH, rbw, rbH_L, rrdn, Y
      real gsw, gswv, rswv
      real gsc
      ! thermodynamic parameters for air
      TairK = Tair + 273.2
      rhocp = cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv = H2oLv0 - 2.365e3*Tair
      slope = (esat(Tair + 0.1) - esat(Tair))/0.1
      psyc = Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar = Patm/(Rconst*TairK)
      weighJ = 1.0
      ! boundary layer conductance for heat - single sided, forced convection
      ! (Monteith 1973, P106 & notes dated 23/12/94)
      gbHu = 0.003*sqrt(windUx/wleaf)    !m/s
      ! raero=0.0                        !aerodynamic resistance s/m
      do ileaf = 1, 2                  ! loop over sunlit and shaded leaves
         ! first estimate of leaf temperature - assume air temp
         Tleaf(ileaf) = Tair
         Tlk = Tleaf(ileaf) + 273.2    !Tleaf to deg K
         ! first estimate of deficit at leaf surface - assume Da
         Dleaf = Dair                !Pa
         ! first estimate for co2cs
         co2cs = co2ca               !mol/mol
         Qapar = (4.6e-6)*Qabs(1, ileaf)
         ! ********************************************************************
         kr1 = 0                     !iteration counter for LE
         do
            !100        continue !    return point for evaporation iteration
            ! single-sided boundary layer conductance - free convection (see notes 23/12/94)
            Gras = 1.595e8*abs(Tleaf(ileaf) - Tair)*(wleaf**3)     !Grashof
            gbHf = 0.5*Dheat*(Gras**0.25)/wleaf
            gbH = gbHu + gbHf                         !m/s
            rbH = 1./gbH                            !b/l resistance to heat transfer
            rbw = 0.93*rbH                          !b/l resistance to water vapour
            ! Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
            rbH_L = rbH*stom_n/2.                   !final b/l resistance for heat
            rrdn = 1./grdn
            Y = 1./(1.+(rbH_L + raero)/rrdn)
            ! boundary layer conductance for CO2 - single side only (mol/m2/s)
            gbc = Cmolar*gbH/1.32            !mol/m2/s
            gsc0 = gsw0/1.57                        !convert conductance for H2O to that for CO2
            varQc = 0.0
            weighR = 1.0
            ! respiration
            Aleafx = -0.0089*Vcmxx*exp(0.069*(Tlk - 293.2))
            gsc = gsc0
            ! choose smaller of Ac, Aq
            if (ISNAN(Aleafx)) then
               
               stop
            endif
            Aleaf(ileaf) = Aleafx                     !mol CO2/m2/s
            ! calculate new values for gsc, cs (Lohammer model)
            co2cs = co2ca - Aleaf(ileaf)/gbc
            co2Ci(ileaf) = co2cs - Aleaf(ileaf)/gsc
            ! scale variables
            gsw = gsc*1.56                              !gsw in mol/m2/s
            gswv = gsw/Cmolar                           !gsw in m/s
            rswv = 1./gswv
            ! calculate evap'n using combination equation with current estimate of gsw
            Eleaf(ileaf) = (slope*Y*Rnstar(ileaf) + rhocp*Dair/(rbH_L + raero))/   &
                &      (slope*Y + psyc*(rswv + rbw + raero)/(rbH_L + raero))!*omega!*AMIN1(0.5, AMAX1(0.333, 0.333 + omega))
            ! calculate sensible heat flux
            Hleaf(ileaf) = Y*(Rnstar(ileaf) - Eleaf(ileaf))
            ! calculate new leaf temperature (K)
            Tlk1 = 273.2 + Tair + Hleaf(ileaf)*(rbH/2.+raero)/rhocp
            ! calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
            Dleaf = psyc*Eleaf(ileaf)/(rhocp*gswv)
            gbleaf(ileaf) = gbc*1.32*1.075
            gsleaf(ileaf) = gsw

            ! compare current and previous leaf temperatures
            if (abs(Tlk1 - Tlk) .le. 0.1) exit
            if (kr1 .gt. 500) exit
            ! update leaf temperature
            Tlk = Tlk1
            Tleaf(ileaf) = Tlk1 - 273.2
            kr1 = kr1 + 1
         end do                          !solution not found yet
10       continue
      end do
      return
   end subroutine agsean_ngt


   subroutine photosyn()
      ! calculate Vcmax, Jmax at leaf temp (Eq 9, Harley et al 1992)
      ! turned on by Weng, 2012-03-13
      ! VcmxT = Vjmax(Tlkx,Trefk,Vcmx1,Eavm,Edvm,Rconst,Entrpy)
      ! eJmxT = Vjmax(Tlkx,Trefk,eJmx1,Eajm,Edjm,Rconst,Entrpy)
      ! real Sps,CO2Ca,CO2Csx,Dleaf,Tlk,Qapar,Gbc
      ! real theta,a1,Ds0,fwsoil,varQc,weighR 
      ! real gsc0,alpha                              
      ! real Vcmxx,eJmxx,weighJ,conKc0,conKo0,Ekc,Eko,o2ci 
      ! real Rconst,Trefk,Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2 
      ! real Aleafx,Gscx,gddonset

      real CO2Csx

      real TminV, TmaxV, ToptV, TminJ, TmaxJ, ToptJ, Tlf
      real VcmxT, eJmxT, eJ
      real conKcT, conKoT
      real Rd, Tdiff, gammas, gamma
      real X, Gma, Bta
      real Acx, Aqx ! from ciandA

      CO2Csx = AMAX1(CO2Cs, 0.6*CO2Ca)
      ! check if it is dark - if so calculate respiration and g0 to assign conductance
      if (Qapar .le. 0.) then                            !night, umol quanta/m2/s
         Aleafx = -0.0089*Vcmxx*exp(0.069*(Tlk - 293.2))   ! original: 0.0089 Weng 3/22/2006
         Gscx = gsc0
      end if
      ! calculate  Vcmax, Jmax at leaf temp using Reed et al (1976) function J appl Ecol 13:925
      TminV = gddonset/10.  ! original -5.        !-Jiang Jiang 2015/10/13
      TmaxV = 50.
      ToptV = 35.

      TminJ = TminV
      TmaxJ = TmaxV
      ToptJ = ToptV

      Tlf = Tlk - 273.2
      VcmxT = VJtemp(Tlf, TminV, TmaxV, ToptV, Vcmxx)
      eJmxT = VJtemp(Tlf, TminJ, TmaxJ, ToptJ, eJmxx)
      ! calculate J, the asymptote for RuBP regeneration rate at given Q
      eJ = weighJ*fJQres(eJmxT, alpha, Qapar, theta)
      ! calculate Kc, Ko, Rd gamma*  & gamma at leaf temp
      conKcT = EnzK(Tlk, Trefk, conKc0, Rconst, Ekc)
      conKoT = EnzK(Tlk, Trefk, conKo0, Rconst, Eko)
      ! following de Pury 1994, eq 7, make light respiration a fixed proportion of Vcmax
      Rd     = 0.0089*VcmxT*weighR                              !de Pury 1994, Eq7
      Tdiff  = Tlk - Trefk
      gammas = gam0*(1.+gam1*Tdiff + gam2*Tdiff*Tdiff)       !gamma*
      ! gamma = (gammas+conKcT*(1.+O2ci/conKoT)*Rd/VcmxT)/(1.-Rd/VcmxT)
      gamma = 0.0
      ! ***********************************************************************
      ! Analytical solution for ci. This is the ci which satisfies supply and demand
      ! functions simultaneously
      ! calculate X using Lohammer model, and scale for soil moisture
      a1 = 1./(1.-0.7)
      X = a1*fwsoil/((co2csx - gamma)*(1.0 + Dleaf/Ds0))
      ! calculate solution for ci when Rubisco activity limits A
      Gma = VcmxT
      Bta = conKcT*(1.0 + o2ci/conKoT)

      call ciandA(Gma, Bta, gsc0, X, Rd, co2Csx, gammas, Acx)
      ! calculate +ve root for ci when RuBP regeneration limits A
      Gma = eJ/4.
      Bta = 2.*gammas
      ! calculate coefficients for quadratic equation for ci
      call ciandA(Gma, Bta, gsc0, X, Rd, co2Csx, gammas, Aqx)
      ! choose smaller of Ac, Aq
      sps = AMAX1(0.001, sps)                  !Weng, 3/30/2006
      Aleafx = (amin1(Acx, Aqx) - Rd) !*sps     ! Weng 4/4/2006

      ! if(Aleafx.lt.0.0) Aleafx=0.0    ! by Weng 3/21/2006
      ! calculate new values for gsc, cs (Lohammer model)
      CO2csx = co2ca - Aleafx/Gbc
      Gscx = gsc0 + X*Aleafx  ! revised by Weng
      ! if (itest .gt.  944) then
      !    iitest = iitest + 1
      !    if (mod(iitest, 2520) .eq. 0)then
      !       write(*,*)"Aleafx: ",Aleafx,Gma, Bta, gsc0, X, Rd, co2Csx, gammas, Aqx
      !       ! write(*,*)"ej: ",weighJ,fJQres(eJmxT, alpha, Qapar, theta),eJmxT, alpha, Qapar, theta
      !       ! write(*,*)"eJmxT",eJmxT,VJtemp(Tlf, TminJ, TmaxJ, ToptJ, eJmxx),Tlf, TminJ, TmaxJ, ToptJ, eJmxx
      !       ! eJmx0*scalex,
      !    endif 
      ! endif
      ! test_gpp = (/Aleafx,Gma, Bta, gsc0, X, Rd, co2Csx, gammas, Aqx/)
      ! test_gpp = (/Aleafx, VcmxT, co2cs, ej,weighJ,fJQres(eJmxT, alpha, Qapar, theta),eJmxT, Qapar, theta/)
      ! test_gpp = (/Tlf, TminV, TmaxV, ToptV, Vcmxx,TminJ, TmaxJ, ToptJ, eJmxx/)
      return
   end subroutine photosyn


   subroutine ciandA(Gma,Bta,g0,X,Rd,co2Csx,gammas,Aquad)      ! Gma,Bta,g0,X,Rd,co2Cs,gammas,ciquad,Aquad
      real Gma,Bta,g0,X,Rd,co2Csx,gammas,ciquad,Aquad
      real b2, b1, b0, bx
      ! calculate coefficients for quadratic equation for ci
      b2 =  g0 + X*(Gma - Rd)
      b1 =  (1.-co2Csx*X)*(Gma - Rd) + g0*(Bta - co2Csx) - X*(Gma*gammas + Bta*Rd)
      b0 = -(1.-co2Csx*X)*(Gma*gammas + Bta*Rd) - g0*Bta*co2Csx

      bx = b1*b1 - 4.*b2*b0
      if (bx .gt. 0.0) then
         ! calculate larger root of quadratic
         ciquad = (-b1 + sqrt(bx))/(2.*b2)
      end if

      IF (ciquad .lt. 0 .or. bx .lt. 0.) THEN
         Aquad = 0.0
         ciquad = 0.7*co2Csx
      ELSE
         Aquad = Gma*(ciquad - gammas)/(ciquad + Bta)
      END IF
      return
   end


!  functions used in canopy -------------------------------
   real function sinbet()
      real rad, sinlat, coslat, sindec, cosdec, A, B
      ! sin(bet), bet = elevation angle of sun
      ! calculations according to Goudriaan & van Laar 1994 P30
      rad = pi/180.
      ! sine and cosine of latitude
      sinlat = sin(rad*lat)
      coslat = cos(rad*lat)
      ! sine of maximum declination
      sindec = -sin(23.45*rad)*cos(2.0*pi*(doy + 10.0)/365.0)
      cosdec = sqrt(1.-sindec*sindec)
      ! terms A & B in Eq 3.3
      A = sinlat*sindec
      B = coslat*cosdec
      sinbet = A + B*cos(pi*(hour - 12.)/12.)
      return
   end

   ! real function esat(T)
   !    real T
   !    ! returns saturation vapour pressure in Pa
   !    esat = 610.78*exp(17.27*T/(T + 237.3))
   !    return
   ! end

   ! ****************************************************************************
   ! Reed et al (1976, J appl Ecol 13:925) equation for temperature response
   ! used for Vcmax and Jmax
   real function VJtemp(Tlf,TminVJ,TmaxVJ,ToptVJ,VJmax0)
      real Tlf,TminVJ,TmaxVJ,ToptVJ,VJmax0
      real pwr
      if (Tlf .lt. TminVJ) Tlf = TminVJ   !constrain leaf temperatures between min and max
      if (Tlf .gt. TmaxVJ) Tlf = TmaxVJ
      pwr    = (TmaxVJ - ToptVJ)/(ToptVj - TminVj)
      VJtemp = VJmax0*((Tlf - TminVJ)/(ToptVJ - TminVJ))*     &
               &       ((TmaxVJ - Tlf)/(TmaxVJ - ToptVJ))**pwr
      return
   end

   ! ****************************************************************************
   real function fJQres(eJmx,alpha,Q,theta)
      real eJmx,alpha,Q,theta
      real AX, BX, CX
      AX = theta                                 !a term in J fn
      BX = alpha*Q + eJmx                          !b term in J fn
      CX = alpha*Q*eJmx                          !c term in J fn
      if ((BX*BX - 4.*AX*CX) >= 0.0) then
         fJQres = (BX - SQRT(BX*BX - 4.*AX*CX))/(2*AX)
      else
         fJQres = (BX)/(2*AX)                   !Weng 10/31/2008
      end if
      return
   end

   ! *************************************************************************
   real function EnzK(Tk,Trefk,EnzK0,Rconst,Eactiv)
      real Tk,Trefk,EnzK0,Rconst,Eactiv
      real temp1 
      temp1 = (Eactiv/(Rconst*Trefk))*(1.-Trefk/Tk)
      ! if (temp1<50.)then
      EnzK = EnzK0*EXP((Eactiv/(Rconst*Trefk))*(1.-Trefk/Tk))
      ! else
      ! EnzK = EnzK0*EXP(50.)                                          ! Weng 10/31/2008
      ! endif
      return
   end

end module mod_vegetation
