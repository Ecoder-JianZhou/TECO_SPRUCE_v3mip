module mod_soil
    use mod_data

    implicit none
    real(kind=8) FLDCAP, WILTPT_x  ! Jian: use WILTPT_x replace the WILTPT

    contains
    ! subroutine for soil moisture
    subroutine soilwater() 
        ! All of inputs, the unit of water is 'mm', soil moisture or soil water content is a ratio
        implicit none
        real infilt_max
        real DWCL(10), evapl(10), wupl(10), Tr_ratio(10)
        real SRDT(10), rain_new, rain_t
        integer nfr
        real infilt_dbmemo, twtadd, wtadd, omegaL(10)
        real exchangeL,supply,demand
        real Tsrdt, tr_allo
        real plantup(10), vtot
        real zmax,thetasmin,zthetasmin,az
        real zwt1,zwt2,zwt3
        real fw(10), ome(10)
        real test_a, test_b
        real omegaL_max, omegaL_min

        infilt_max = 15.
        WILTPT_x   = wsmin/100.000
        FLDCAP     = wsmax/100.000
        
        do i = 1,10
            dwcl(i)  = 0.0
            evapl(i) = 0.0
            WUPL(i)  = 0.0
            SRDT(i)  = 0.0
            DEPTH(i) = 0.0
            if (i>3) wcl(i) = FLDCAP    ! Jian: set underground 30cm to saturated.
        enddo
        ! Layer volume (cm3)
        DEPTH(1) = THKSL(1) ! Determine which layers are reached by the root system. 
        DO i=2,10
            DEPTH(i)=DEPTH(i-1)+THKSL(i)
        enddo
        do i=1,10
            IF(rdepth.GT.DEPTH(i)) nfr=i+1
        enddo
        IF (nfr.GT.10) nfr=10
 
        ! ---------- added for soil thermal    
        if (do_soilphy) then 
            rain_new = rain
            ! if (ta .lt. -4.) rain_new =0.    ! if (ta .lt. -0.4) rain_new =0.       !dbice   !tuneice     
            rain_t = melt/24+rain_new        ! here it defines how the melt water is added to water input; add melted water hourly
            infilt = infilt+rain_t
            if (ice(1) .gt. 0.0) then        ! Jian: no use in ice? 2022/11/15
                !infilt = 0.0
            endif
        else
            infilt = infilt + rain           ! ..int commented lines for soil thermal module, included in the previous loop; !mm/hour  
        endif     
        ! write(*,*)"test_infilt: ", infilt, rain, melt, melt/24
        ! if (rain>0) write(*,*)"test_infilt: ", rain, infilt
        infilt_dbmemo = infilt    
        ! water infiltration through layers; Loop over all soil layers.
        TWTADD=0
        IF(infilt.GE.0.0)THEN
            ! Add water to this layer, pass extra water to the next.
            WTADD  = AMIN1(INFILT,infilt_max,AMAX1((FLDCAP-wcl(1))*thksl(1)*10.0,0.0)) ! from cm to mm
            WCL(1) = (WCL(1)*(thksl(1)*10.0)+WTADD)/(thksl(1)*10.0)
            TWTADD = TWTADD+WTADD       !calculating total added water to soil layers (mm)
            INFILT = INFILT-WTADD       !update infilt
        ENDIF
        ! write(*,*)"TEST_WCL: ",WCL(1),thksl(1),WTADD,infilt,AMAX1((FLDCAP-wcl(1))*thksl(1)*10.0,0.0)
        if (do_soilphy) then 
            runoff= INFILT*0.005   !(infilt_rate = 0.0017 defined earlier by Yuan, changed  to 0.001 by shuang )
        else
            runoff=INFILT*0.001    ! Shuang added this elseif line! Shuang Modifed  Mar16 used to be 0.0019, the water lose too much lowest wt was >400
        endif   
        ! Modified based on Ma et al.,2022

        ! runoff = infilt*0.001
        infilt = infilt-runoff

        !----------------------------------------------------------------------------------------------------
        if (transp .gt. 0.2 .and. transp .le. 0.22) then
            infilt = infilt+transp*0.4
        else if (transp .gt. 0.22) then
            ! infilt = infilt+infilt*0.0165
            infilt = infilt+transp*0.8
            ! infilt = infilt+0.22*0.4+(transp-0.22)*0.9
        else
            infilt = infilt+transp*0.001
        endif

        if (evap .ge. 0.1 .and. evap .le. 0.15) then
            infilt = infilt+evap*0.4
        else if (evap .gt. 0.15) then
            infilt = infilt+evap*0.8
        else
            infilt = infilt+evap*0.001
        endif
        !----------------------------------------------------------------------------------------------------   
        ! water redistribution among soil layers
        do i=1,10
            wsc(i) = Amax1(0.00,(wcl(i)-WILTPT_x)*THKSL(i)*10.0)
            if (do_soilphy) then ! ..int commented lines for soil thermal 
                omegaL(i)=Amax1(0.001,(liq_water(i)*100./thksl(i)-WILTPT_x)/(FLDCAP-WILTPT_x))
            else
                omegaL(i)=Amax1(0.001,(wcl(i)-WILTPT_x)/(FLDCAP-WILTPT_x))
            endif        
        enddo

        ! write(*,*)"test_wsc: ", wsc(1),(wcl(1)-WILTPT_x)*THKSL(1)*10.0, wcl(1),WILTPT_x,THKSL(1)
        supply = 0.0
        demand = 0.0
        test_a = 0.0
        do i=1,9
            if(omegaL(i).gt.0.3)then
                supply    = wsc(i)*(omegaL(i)-0.3)   ! supply=wsc(i)*omegaL(i)
                demand    = (FLDCAP-wcl(i+1))*THKSL(i+1)*10.0*(1.0-omegaL(i+1))
                exchangeL = AMIN1(supply,demand)
                wsc(i)    = wsc(i)- exchangeL
                wsc(i+1)  = wsc(i+1)+ exchangeL
                wcl(i)    = wsc(i)/(THKSL(i)*10.0)+WILTPT_x
                wcl(i+1)  = wsc(i+1)/(THKSL(i+1)*10.0)+WILTPT_x
                if (i .eq. 1)then
                    test_a = exchangeL
                endif
            endif
            ! if (i>3) wcl(i) = 1    ! Jian: set underground 30cm to saturated.
            ! if (i<4)then ! Jian: if 30cm water too low, ground water feedback to depth 1-3?
            ! if (i .eq. 3) then
            !     if (omegaL(3)<0.6) then
            !         supply    = Amax1(wsc(4)*(omegal(4)-0.3),0.0)   ! supply=wsc(i)*omegaL(i)
            !         demand    = (FLDCAP-wcl(3))*THKSL(3)*10.0*(1-omegaL(3))
            !         exchangeL = AMIN1(supply,demand)
            !         wsc(3) = wsc(3) + exchangeL
            !         write(*,*)"exchangeL1:",exchangeL
            !     endif
            ! endif
            ! if (i .eq. 2)then
            !     if (omegaL(2)<0.3) then
            !         supply    = Amax1(wsc(3)*(omegal(3)-0.3),0.)   ! supply=wsc(i)*omegaL(i)
            !         demand    = (FLDCAP-wcl(2))*THKSL(2)*10.0*(1-omegaL(2))
            !         exchangeL = AMIN1(supply,demand)
            !         wsc(2) = wsc(2) + exchangeL
            !         write(*,*)"exchangeL2:",exchangeL
            !     endif
            ! endif
        enddo
        

        ! if (wcl(10)>FLDCAP)then 
            wsc(10) = wsc(10)-wsc(10)*0.00001     ! Shuang modifed
            runoff  = runoff+wsc(10)*0.00001     ! Shuang modifed
        ! endif
        wcl(10) = wsc(10)/(THKSL(10)*10.0)+WILTPT_x
        ! end of water redistribution among soil layers
        ! Redistribute evaporation among soil layers
        Tsrdt = 0.0
        DO i=1,10
            ! Fraction of SEVAP supplied by each soil layer
            SRDT(I) = EXP(-6.73*(DEPTH(I)-THKSL(I)/2.0)/100.0) !/1.987 ! SRDT(I)=AMAX1(0.0,SRDT(I)*(wcl(i)-wiltpt)) !*THKSL(I))
            Tsrdt   = Tsrdt+SRDT(i)  ! to normalize SRDT(i)
        enddo
        ! write(*,*)"test_wcl: ", wcl
        do i=1,10
            EVAPL(I) = Amax1(AMIN1(evap*SRDT(i)/Tsrdt,wsc(i)),0.0)  !mm
            DWCL(I)  = EVAPL(I)/(THKSL(I)*10.0) !ratio
            wcl(i)   = wcl(i)-DWCL(i)
            ! if (i>3) wcl(i) = 1    ! Jian: set underground 30cm to saturated.
            ! if (DWCL(i)>0.5)write(*,*)
        enddo

        ! write(*,*)"test_evap: ", evap
        ! write(*,*)"test_SRDT: ", SRDT
        ! write(*,*)"test_Tsrdt:", Tsrdt
        
        ! IF(iforcing .Eq. 500*24)then
            
        !     write(*,*)"zwc: ", wsc
        !     stop
        ! endif
       
        evap = 0.0       
        do i=1,10
            evap=evap+EVAPL(I)
        enddo
        ! Redistribute transpiration according to root biomass
        ! and available water in each layer
        tr_allo=0.0
        do i=1,nfr
            tr_ratio(i) = FRLEN(i)*wsc(i) !*(wcl(i)-wiltpt)) !*THKSL(I))
            tr_allo     = tr_allo+tr_ratio(i)
        enddo
        ! write(*,*)"test ...", transp,tr_ratio
        do i=1,nfr
            plantup(i) = AMIN1(transp*tr_ratio(i)/tr_allo, wsc(i)) !mm          
            wupl(i)    = plantup(i)/(thksl(i)*10.0)
            wcl(i)     = wcl(i)-wupl(i)
            ! if (i>3) wcl(i) = 1    ! Jian: set underground 30cm to saturated.
        enddo
        transp = 0.0
        do i=1,nfr
            transp=transp+plantup(i)
        enddo
        ! if(transp > 10.)then
        !     write(*,*)"test ...", iforcing/24, transp, plantup
        !     stop
        ! endif

        ! if (transp .gt. wsc(1) .and. transp .le. (wsc(1)+wsc(2))) then
        !     infilt = infilt+transp*0.4
        ! else if (transp .gt. (wsc(1)+wsc(2))) then
        !     ! infilt = infilt+infilt*0.0165
        !     infilt = infilt+transp*0.6
        !     ! infilt = infilt+0.22*0.4+(transp-0.22)*0.9
        ! else if (transp .gt. (wsc(1)+wsc(2)+wsc(3))) then
        !     infilt = infilt+transp*0.8
        ! else 
        !     infilt = infilt+transp*0.001
        ! endif

        ! if (evap .ge. wsc(1) .and. evap .le. (wsc(1)+wsc(2))) then
        !     infilt = infilt+evap*0.4
        ! else if (transp .gt. (wsc(1)+wsc(2))) then
        !     infilt = infilt+evap*0.6
        ! else if (transp .gt. (wsc(1)+wsc(2)+wsc(3))) then
        !     infilt = infilt+evap*0.8
        ! else
        !     infilt = infilt+evap*0.001
        ! endif

        ! Jian: readd water according to evaporation, transpiration and soil moisture.
        infilt = infilt + amax1(0.0,(1-omega)*transp)
        infilt = infilt + amax1(0.0,(1-omega)*evap)

        ! do i=4,10
        !     wcl(i) = FLDCAP
        ! enddo
        ! IF(infilt.GE.0.0)THEN
        !     ! Add water to this layer, pass extra water to the next.
        !     WTADD  = AMIN1(INFILT,infilt_max,AMAX1((FLDCAP-wcl(1))*thksl(1)*10.0,0.0)) ! from cm to mm
        !     WCL(1) = (WCL(1)*(thksl(1)*10.0)+WTADD)/(thksl(1)*10.0)
        !     TWTADD = TWTADD+WTADD       !calculating total added water to soil layers (mm)
        !     INFILT = INFILT-WTADD       !update infilt
        ! ENDIF
        ! supply = 0.0
        ! demand = 0.0
        ! test_a = 0.0
        ! do i=1,9
        !     if(omegaL(i).gt.0.3)then
        !         supply    = wsc(i)*(omegaL(i)-0.3)   ! supply=wsc(i)*omegaL(i)
        !         demand    = (FLDCAP-wcl(i+1))*THKSL(i+1)*10.0*(1.0-omegaL(i+1))
        !         exchangeL = AMIN1(supply,demand)
        !         wsc(i)    = wsc(i)- exchangeL
        !         wsc(i+1)  = wsc(i+1)+ exchangeL
        !         wcl(i)    = wsc(i)/(THKSL(i)*10.0)+WILTPT_x
        !         wcl(i+1)  = wsc(i+1)/(THKSL(i+1)*10.0)+WILTPT_x
        !         if (i .eq. 1)then
        !             test_a = exchangeL
        !         endif
        !     endif
        ! enddo
        
        do i=1,10
            wsc(i) = Amax1(0.00,(wcl(i)-WILTPT_x)*THKSL(i)*10.0)
        enddo


        ! ---------------------------------------------------------------------------    
        ! water table module starts here
        ! vtot = MAX(145.,wsc(1)+wsc(2)+wsc(3)+infilt)!+wsc(4)+wsc(5)   !total amount of water in top 500mm of soil  mm3/mm2 infilt here is standing water   infilt has distributed to wsc?
        if (do_soilphy) then
            ! vtot = wsc(1)+wsc(2)+wsc(3)+infilt+ice(1)*1000.*(10./9.)+ice(2)*1000.*(10./9.)+ice(3)*1000.*(10./9.)
            ! vtot = wsc(1)+wsc(2)+wsc(3)+ice(1)*1000.*(10./9.)+ice(2)*1000.*(10./9.)+ice(3)*1000.*(10./9.)
            ! vtot = wsc(1)+wsc(2)+wsc(3)+infilt
            ! vtot = (liq_water(1)+liq_water(2)+liq_water(3))*1000+(ice(1)+ice(2)+ice(3))*1000+infilt
            vtot = wsc(1)+wsc(2)+wsc(3)+infilt ! modified based on Ma et al., 2022
            ! vtot = wsc(1)+wsc(2)+wsc(3)+infilt+ice(1)*1000.*(9./10.)+ice(2)*1000.*(9./10.)+ice(3)*1000.*(9./10.)
            ! write(*,*) ice(1)*1000.,ice(2)*1000.,ice(3)*1000.,wsc(1),liq_water(1)*1000.
        else 
            vtot = wsc(1)+wsc(2)+wsc(3)+infilt!+wsc(4)+wsc(5) 
        endif
        ! infilt means standing water according to jiangjiang
        ! vtot = MAX(145.,vtot+145.+rain-evap-transp-runoff)         ! vtot should not be smaller than 145, which is the water content when wt is at -300mm
        ! phi        = 0.56                           ! soil porosity   mm3/mm3   the same unit with theta
        phi        = 0.85
        zmax       = 300                            ! maximum water table depth   mm
        thetasmin  = 0.25                           ! minimum volumetric water content at the soil surface   cm3/cm3
        zthetasmin = 100                            ! maximum depth where evaporation influences soil moisture   mm
        az         = (phi-thetasmin)/zthetasmin     ! gradient in soil moisture resulting from evaporation at the soil surface    mm-1

        zwt1 = -sqrt(3.0*(phi*zmax-vtot)/(2.0*az))
        zwt2 = -(3.0*(phi*zmax-vtot)/(2.0*(phi-thetasmin)))
        zwt3 = vtot-phi*zmax           
        ! write(*,*)"test: ", vtot, zwt1, zwt2, zwt3                        
        if ((zwt1 .ge. -100) .and. (zwt1 .le. 0))   zwt = zwt1  !the non-linear part of the water table changing line
        if (zwt2 .lt. -100)                         zwt = zwt2  !the linear part of the water table changing line
        ! if ((zwt2 .lt. -100) .and. (zwt2 .ge. -300))zwt = zwt2 !the linear part of the water table changing line valid when Vtot>145mm
        ! if (zwt2 .le. -300)                         zwt = -300
        if (phi*zmax .lt. vtot)                     zwt = zwt3  !the linear part when the water table is above the soil surface   
   
        ! water table module ends here
        ! ---------------------------------------------------------------------------------------------------------------
        ! Output fwsoil, omega, and topfws
        ! ..int commented lines below for soil thermal module
        ! do i=1,nfr       
        !     ome(i)=(wcl(i)-WILTPT)/(FLDCAP-WILTPT)
        !     ome(i)=AMIN1(1.0,AMAX1(0.0,ome(i)))
        !     fw(i)=amin1(1.0,3.333*ome(i))
        ! enddo
        ! topfws=amin1(1.0,(wcl(1)-WILTPT)/((FLDCAP-WILTPT)))    
        ! ..int new lines added for soil thermal module 
        do i=1,nfr       
            if (do_soilphy) then 
                ome(i)=(liq_water(i)*100./thksl(i)-WILTPT_x)/(FLDCAP-WILTPT_x)
                ome(i)=AMIN1(1.0,AMAX1(0.0,ome(i)))
            else 
                ome(i)=(wcl(i)-WILTPT_x)/(FLDCAP-WILTPT_x)
                ome(i)=AMIN1(1.0,AMAX1(0.0,ome(i)))
            endif 
            fw(i)=amin1(1.0,3.333*ome(i))
        enddo
        ! write(*,*)"liq_water: ", liq_water
        ! write(*,*)"wcl: ", wcl
        ! write(*,*)"thksl: ", thksl
        ! write(*,*)"FLDCAP, WILTPT: ", FLDCAP, WILTPT_x
        if (do_soilphy) then 
            topfws=amax1(0.0,topfws)
        else 
            topfws=amin1(1.0,(wcl(1)-WILTPT_x)/((FLDCAP-WILTPT_x)))
        endif     

        fwsoil = 0.0
        omega  = 0.0
        do i=1,nfr
            fwsoil = fwsoil+fw(i)*frlen(i)
            omega  = omega+ome(i)*frlen(i)
        enddo
        ! if (omega>1 .or. omega<0) write(*,*)"omega_ome: ", omega
        ! write(*,*)"omega_frlen: ", FRLEN
        ! write(*,*)"liq_water: ", liq_water
        ! write(*,*)"thksl: ", thksl
        ! write(*,*)"FLDCAP, WILTPT: ", FLDCAP, WILTPT
    ! write(*,*) "testa: ", iforcing, wcl(1), DWCL(1),WTADD,test_a,wupl(1), wtadd-DWCL(1)-test_a-wupl(1),rain,melt,infilt
    
        return
    end subroutine soilwater

    ! subroutine snow_d(rain_d,lat,days,ta,snow_dsim,fa,fsub,rho_snow,melt,dcount,decay_m)
    subroutine snow_d()
        ! real lat,tr,daylength,dec,melt,fa,sublim,dsnow,snow_in,decay_m,fsub
        real tr,daylength,dec,dsnow, in_snow
        ! real rain_d,snow_dsim,rho_snow,dcount,ta
        ! integer days
        real snow_dsim_pre

        tr        = 0.0174532925
        dec       = sin(((real(iday)-70.)/365.)*360.*tr)*23.44  ! dec=sin(((real(days)-70.)/365.)*360.*tr)*23.44
        daylength = acos(-tan(lat*tr)*tan(dec*tr))/7.5 
        daylength = daylength/tr/24.
                
        if (snow_dsim .ge. 0.) then
            dcount = dcount +1.
        else 
            dcount =0.
        endif
        sublim=0.
        if (ta .gt. 0. .and. snow_dsim .gt. 0.) sublim=fsub*715.5*daylength*esat(ta)/(ta+273.2)*0.001   ! 0.001 from Pa to kPa
        melt=0.
        ! if (ta .gt. 0. .and. snow_dsim .gt. 0.) melt=fa*(2.63+2.55*ta+0.0912*ta*rain_d)       !yy version
        if (ta .gt. 1.0e-10 .and. snow_dsim .gt. 0.) melt=fa*(2.63+2.55*ta+0.0912*ta*rain_d)   !dbmemo updated version
        ! write(*,*) 'melt=fa*(2.63+2.55*ta+0.0912*ta*rain_d)','fa',fa,'ta',ta,'rain_d',rain_d
        ! if (ta .gt. 0. .and. snow_dsim .gt. 0.) melt=fa*(0.55*ta)     
        if (dcount .gt.0. .and. ta .lt.5.) then
            ! write(*,*)'melt_befor',melt         !dbmemo dbice
            melt=melt*EXP(-decay_m*dcount/365.)  !dbmemo dbice
            ! write(*,*)'melt_after',melt
        endif
        ! write(*,*),EXP(-3.*dcount/365.)
        ! melt=AMIN1(melt, snow_dsim*rho_snow-sublim)
        ! if (melt .lt. 2.) melt=0.
        if (ta .le. 0.) then         ! dbmemo second bug in dbmemo
            in_snow =rain_d
        else
            in_snow = 0.
        endif
        dsnow         = in_snow-sublim-melt 
        snow_dsim_pre = snow_dsim
        snow_dsim     = snow_dsim + dsnow/rho_snow 
        if (snow_dsim .le. 0.0) then 
            snow_dsim=0.0 
            melt = snow_dsim_pre*rho_snow +in_snow-sublim    !! for water part
        endif 
        melt=AMAX1(melt, 0.)
        return
    end subroutine snow_d


    subroutine Tsoil_simu()          
        implicit none 
        real difsv2,difsv1
        real delta
        !real thksl(10)
        real ufw(10),frac_ice1,frac_ice2
        ! real,dimension(10):: Tsoill,liq_water
        ! real,dimension(11)::testout
        real temph1,temph2     
        real thkns1,thkns2
        real flux_snow
        real condu_water,shcap_water,shcap_ice              !,shcap_snow,condu_snow,depth_ex
        real albedo_water,ice_incr,heat_excess,heat_adjust  !,ice_tw,water_tw
        real inter_var,latent_heat_fusion                   !,QLsoil,Rsoilab3
        real resdh,dnr,dsh,dgh,dle,drsdh                    ! rhocp,slope,psyc,Cmolar,fw1,resoil,rLAI,resht,
        ! real f_om,theta_sat_om,b_om,b_min,phi_om,phi_min,theta_sat,b_tot,phi_sat,gravi
        real water_table_depth,temph_water,temph_snow
        real condu_air,shcap_air,condu(10), shcap(10), condu_ice,tsoill_pre, thd_t
        real ice_density,condu_soil,shcap_soil 
        real resht_lai,snow_depth_t
        real condu_s,tsoill_0,diff_air,d_cor                !,condu_b,dcount,dcount_soil
        real sftmp_pre
        integer n_layers
        real, allocatable ::depth_z(:) 
        n_layers=10
        allocate(depth_z(n_layers))      
        ! write(*,*),thd_t
        ! soil thermal conductivity W m-2 K-1
        ice_density = 916.!916.
        thkns1      = thksl(1)/4.             ! thkns1=thksl(1)/2.
        shcap_ice   = 2117.27*ice_density
        condu_ice   = 2.29
        condu_water = 0.56!0.56
        shcap_water = 4188000.
        condu_soil  = 0.25
        shcap_soil  = 2600000.
        condu_s     = 0.25
        ! thd_t=0.0
        thd_t       = -1.0        
        diff_snow   = 3600.*condu_snow/shcap_snow*10000.
        diff_s      = 3600.*condu_b/shcap_soil*10000.
        latent_heat_fusion = 333700.   ! j kg-1
        condu_air    = 0.023
        shcap_air    = 1255.8
        diff_air     = 3600.*condu_air/shcap_air*10000.      
        water_tw     = zwt*0.001-ice_tw ! might means total water that is liquid, add up all layers
        water_table_depth=zwt*0.1
        snow_depth_t = snow_depth - 0.46*0.0     ! warming in Tair impact on snow_depth
                                                   ! in unit cm 0.46 based on snow_depth vs. tair regression     
        if (snow_depth_t .lt. thd_snow_depth) snow_depth_t =0.0
        if (snow_depth_t .gt. 0.) then
            dcount_soil = dcount_soil +1./24.
        else 
            dcount_soil =0.
        endif
    
        if (water_table_depth .lt. 4. .and. water_table_depth .gt. 0.0) water_table_depth =0.    ! avoid numerical issues when 
        ! if (water_table_depth .lt. -99.) water_table_depth =-30.    ! temporary for NaN    
        ! if (water_table_depth .lt. -299.) water_table_depth =-30.    ! -299 is a more reasonable value Shuang Ma         
        albedo_water = 0.1      
        ! soil water conditions
        WILTPT       = wsmin/100.
        FILDCP       = wsmax/100.
        TairK        = Tair+273.2            
        flux_snow    = 0.0   
        depth_z      = (/0., 0., 0., 0., 0., 0., 0.,0.,0.,0./) 
        ! ..int add unfrozen water ratio
        ! ufw=(/0.0042,0.0063,0.0063,0.0063,0.0063,0.0063,0.0063,0.0063,0.0063,0.0063/)
        ufw=(/0.0163,0.0263,0.0563,0.0563,0.0563,0.1162,0.1162,0.1162,0.1162,0.1162/)
        ! ufw=(/0.0042,0.009,0.009,0.0563,0.0563,0.1162,0.1162,0.1162,0.1162,0.1162/)
        frac_ice1 = 0.01!0.015
        frac_ice2 = 0.001!0.01
        ! if (snow_depth_t .gt. 0.0) then 
        !     emsoil =0.98
        ! elseif (water_table_depth .gt. 0.0) then
        !     emsoil =0.99
        ! endif
        QLsoil   = emsoil*sigma*((sftmp+273.2)**4)
        Rsoilab3 = (QLair+QLleaf)*(1.0-rhoS(3))-QLsoil 
        if(ISNAN(Rsoilab3)) then
            write(*,*)"Rsoilab3 is NaN: ", Rsoilab3, QLair, QLleaf, rhoS, QLsoil
            stop
        endif      
        ! Total radiation absorbed by soil
        if (snow_depth_t .gt. 0.0) then 
            Rsoilabs = (Rsoilab1+Rsoilab2)*(1-albedo_snow)/(1-0.1)+Rsoilab3  
        elseif (water_table_depth .gt. 0.0) then 
            Rsoilabs = (Rsoilab1+Rsoilab2)*(1-albedo_water)/(1-0.1)+Rsoilab3  
        else
            Rsoilabs = Rsoilab1+Rsoilab2+Rsoilab3
        endif
          
        ! thermodynamic parameters for air
        rhocp  = cpair*Patm*AirMa/(Rconst*TairK)      
        H2OLv  = H2oLv0-2.365e3*Tair
        slope  = (esat(Tair+0.01)-esat(Tair))/0.01   
    
        psyc   = Patm*cpair*AirMa/(H2OLv*H2OMw)
        Cmolar = Patm/(Rconst*TairK)
        fw1    = AMIN1(AMAX1((FILDCP-wcl(1))/(FILDCP-WILTPT),0.3),1.0)    

        if (water_table_depth .gt. 0.0) then 
            Rsoil = 0. 
        else 
            Rsoil=30.*exp(0.2/fw1)
        endif 
        ! Rsoil=40.
        ! Rsoil=5.
        rLAI=exp(FLAIT)     
        ! latent heat flux into air from soil
        !       Eleaf(ileaf)=1.0*
        ! &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
        ! &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
    
        Esoil=(slope*(Rsoilabs-G)+rhocp*Dair/(raero+rLAI))/       &
               &      (slope+psyc*(Rsoil/(raero+rLAI)+1.))
        resht_lai=resht*FLAIT
        ! resht_lai= resht*exp(FLAIT)/15. ! need improvement, should be a function of LAI 
        ! if (water_table_depth .gt. 0.0) resht_lai=resht/FLAIT*0.2

        ! resht_lai=200.      
        Hsoil=rhocp*(sftmp-Tair)/resht_lai   
        ! Hsoil=1010.*1.17*(sftmp-Tair)/resht_lai
        i=1;
        condu(i) = (FILDCP-wcl(i))*condu_air+liq_water(i)/(thksl(i)*0.01)*condu_water+ &
                    &  ice(i)/(thksl(i)*0.01)*condu_ice +(1-FILDCP)*condu_soil
        shcap(i) = (FILDCP-wcl(i))*shcap_air+liq_water(i)/(thksl(i)*0.01)*shcap_water+ &
                    &  ice(i)/(thksl(i)*0.01)*shcap_ice +(1-FILDCP)*shcap_soil
        difsv1   = 3600.*condu(i)/shcap(i)*10000.    
        G        = condu(1)*(sftmp-tsoill(1))/(thksl(1)/2.*0.01)
        if (snow_depth_t .gt. 0.0) then 
            ! write(*,*)"test here ..."
            G = condu_snow*(sftmp-Tsnow)/(snow_depth_t/2.*0.01)
        endif
        ! write(*,*)"test_G: ", condu(1), sftmp, tsoill(1), thksl(1), condu_snow, Tsnow, snow_depth_t
        ! thksl(1)
        ! G=0.   
        ! Residual heat energy.
        RESDH = Rsoilabs-Hsoil-Esoil-G
        ! First derivative of net radiation; sensible heat; ground heat;
        DNR   = 4.*emsoil*sigma*(sftmp+273.2)**3
        DSH   = rhocp/resht_lai 
        DGH   = condu_s/(thksl(1)/2.*0.01)
        DLE   = (DNR+DGH)*slope/(slope+psyc*(Rsoil/(raero+rLAI)+1.))      
        drsdh = -DNR-DSH-DGH-DLE
        ! Calculate increment DELTA.
        DELTA     = resdh/drsdh
        sftmp_pre = sftmp
        sftmp     = sftmp-DELTA
        if (ABS(sftmp_pre -sftmp) .gt. 20. ) sftmp=sftmp_pre  
        tsoill_0  = sftmp
        if (isnan(Rsoilabs)) then
            write(*,*) "Rsoilabs is NAN: ", albedo_snow, Rsoilab1, Rsoilab2, Rsoilab3
        endif

        ! if (isnan(Esoil)) then
        !     write(*,*) "Esoil is NAN: ", albedo_snow, Rsoilab1, Rsoilab2, Rsoilab3
        ! endif
        ! write(*,*)"test sftmp: ", sftmp, delta, sftmp_pre, resdh, drsdh, Rsoilabs, Hsoil, Esoil, G
        if(isnan(sftmp)) then
            write(*,*)"test sftmp: ", sftmp, delta, sftmp_pre, resdh, drsdh, Rsoilabs, Hsoil, Esoil, G
            stop
        endif
        do i=1,10
            Tsoill_pre = tsoill(i)  
            ! if (i .eq. 1) then 
            !     depth_z(1) = thksl(1)
            ! else 
            !     depth_z(i) = depth_z(i-1)+thksl(i)
            ! endif  

            if (water_table_depth .lt. 0.0 .and. -water_table_depth .lt. depth_z(i)) then
                liq_water(i) = FILDCP*thksl(i)*0.01-ice(i)
            else
                liq_water(i) = wcl(i)*thksl(i)*0.01-ice(i)
            endif          
            if (i .eq. 1) then 
                depth_z(1) = thksl(1)
            else 
                depth_z(i) = depth_z(i-1)+thksl(i)
            endif
            
            ! Jian: some error in here, thksl is 10 layers
            ! thkns2 = (thksl(i)+thksl(i+1))/2.
            ! change to below ...
            if (i<10) then
                thkns2 = (thksl(i)+thksl(i+1))/2.
            else
                thkns2 = (thksl(i-1)+thksl(i))/2.
            endif

            if (i .eq. 10) then
                difsv2 = 3600.*condu(i)/shcap(i)*10000. 
            else
                condu(i+1) = (FILDCP-wcl(i+1))*condu_air+liq_water(i+1)/(thksl(i+1)*0.01)*condu_water+ &
                                &  ice(i+1)/(thksl(i+1)*0.01)*condu_ice +(1-FILDCP)*condu_soil
                shcap(i+1) = (FILDCP-wcl(i+1))*shcap_air+liq_water(i+1)/(thksl(i+1)*0.01)*shcap_water+ &
                                &  ice(i+1)/(thksl(i+1)*0.01)*shcap_ice +(1-FILDCP)*shcap_soil 
                difsv2     = 3600.*condu(i+1)/shcap(i+1)*10000.
            endif    
            
            ! Jian: also error here, layer is 10...
            ! temph2 = (difsv1+difsv2)*(Tsoill(i)-Tsoill(i+1))/thkns2 
            if (i<10)then
                temph2 = (difsv1+difsv2)*(Tsoill(i)-Tsoill(i+1))/thkns2 
            else
                temph2 = (difsv1+difsv2)*(Tsoill(i-1)-Tsoill(i))/thkns2
            endif 
            !!!!!!!!!!!!!!!!!!!! start first layer !!!!!!!!!!!!!!!!!!!!!!
            !!!!!! adjust if there are snow or water layer above !!!!!!!!!!!!!!!!!!!!
            if(i.eq.1) then
                if (snow_depth_t .gt. 0.) then   
                    temph_snow = Amin1(diff_snow,difsv1)*(Tsnow-Tsoill(1))/((snow_depth_t+thksl(1))/2.)
                    Tsnow      = Tsnow+(exp(-depth_ex*snow_depth_t)*diff_snow*(sftmp-Tsnow)/(snow_depth_t/2.) &
                                    &    -temph_snow)/(snow_depth_t/2.+(snow_depth_t+thksl(1))/2.) 
                    Tsoill(1)  = Tsoill(1)+(temph_snow &
                                    &    -temph2)/((snow_depth_t+thksl(1))/2.+thkns2) 
                    if (Tsnow .gt.0.0) then 
                        Tsnow     = 0.0   
                        Tsoill(1) = 0.
                    endif
                    drsdh    = 0.0    ! temporarily set drsdh =0 for heat adjustment of soil when  
                    tsoill_0 = (Tsoill(1)+Tsnow)/2.
                elseif (water_table_depth .gt. 0.) then  
                    temph_water = (3600.*condu_water/shcap_water*10000.+difsv1)*(Twater-Tsoill(1))/((water_table_depth+thksl(1))/2.)! there is snow layer 
                    Twater      = Twater+(2.*3600.*condu_water/shcap_water*10000.*(sftmp-Twater)/(water_table_depth/2.) &
                                    &        -temph_water)/(water_table_depth/2.+(water_table_depth+thksl(1))/2.) 
                    !!!!!!!!!!!!!!!!!!  Phase change surface water !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    if (Twater .lt. 0.0 .and. water_tw .gt. 0.0) then  ! freeze 
                        heat_excess = -(shcap_water/360000.*water_tw*100.-drsdh)*Twater
                        ice_incr    = heat_excess*3600./latent_heat_fusion/ice_density
                        !!-----------int add mechanism of unfrozen water in frozen soil layers, typically happens in high latitude region
                        !!              according to obs soil water content, winter water never goes below 0.063 at -20cm and 0.042 at surface layer
                        !!       !tuneice
                        !                if (ice_incr .lt. 0.) then   
                        !                   if (i .eq. 1.) then
                        !                       if (liq_water(i) .le. ufw(i)) then
                        !!                           ice_incr = 0. !ice_incr*0.1
                        !                          ice_incr = ice_incr*frac_ice1
                        !                       endif
                        !    !               elseif (i .eq. 2) then
                        !    !                   if (liq_water(i) .le. 0.063) then
                        !    !                       ice_incr = 0.
                        !    !                   endif
                        !                   elseif ( i .gt. 1.) then
                        !                       if (liq_water(i) .le. ufw(i)) then
                        !!                           ice_incr = ice_incr*0. !0.9
                        !                           ice_incr = ice_incr*frac_ice2
                        !                       endif
                        !                   endif
                        !                endif
                        !! $$$$$$$$$$$   
                       
                        ! write(*,*)'water_tw',water_tw
                        if (ice_incr .lt. water_tw) then
                            ice_tw   = ice_tw +ice_incr
                            water_tw = water_tw-ice_incr
                            Twater   = 0.0
                            Tice     = 0.0
                        else
                            ice_tw   = ice_tw +water_tw
                            water_tw = 0.0
                            Tice     = Tice - latent_heat_fusion*(ice_incr-water_tw)*ice_density/(shcap_ice*ice_tw)
                        endif     
                    elseif (Twater .gt. 0.0 .and. ice_tw .gt. 0.0) then    ! thraw              
                        heat_excess  = (shcap_water/360000.*ice_tw*100.-drsdh)*Twater
                        ice_incr     = heat_excess*3600./latent_heat_fusion/ice_density
                        !! $$$$$$$$$$$
                        !! $$$$$$$$$$$      !tuneice
                        !                if (ice_incr .lt. 0.) then   
                        !                   if (i .eq. 1.) then
                        !                       if (liq_water(i) .le. ufw(i)) then
                        !!                           ice_incr = 0. !ice_incr*0.1
                        !                          ice_incr = ice_incr*frac_ice1
                        !                       endif
                        !    !               elseif (i .eq. 2) then
                        !    !                   if (liq_water(i) .le. 0.063) then
                        !    !                       ice_incr = 0.
                        !    !                   endif
                        !                   elseif ( i .gt. 1.) then
                        !                       if (liq_water(i) .le. ufw(i)) then
                        !!                           ice_incr = ice_incr*0. !0.9
                        !                           ice_incr = ice_incr*frac_ice2
                        !                       endif
                        !                   endif
                        !                endif
                        !! $$$$$$$$$$$                   
                        !! $$$$$$$$$$$                   
             
                        if (ice_incr .lt. ice_tw) then
                            ice_tw   = ice_tw -ice_incr
                            water_tw = water_tw+ice_incr
                            Twater   = 0.0
                            Tice     = 0.0
                        else
                            water_tw = water_tw +ice_tw
                            ice_tw   = 0.0
                            Twater   = Twater + latent_heat_fusion*(ice_incr-ice_tw)*ice_density/(shcap_water*water_tw)
                        endif
                    endif                       
                    !!!!!!!!!!!!!!!!!!!!!!!!! end of phase change for surface layer !!!!!!!!!!!!!!!!!!!  
                    temph2=(difsv1+3600.*condu_water/shcap_water*10000.)*(Tsoill(i)-Tsoill(i+1))/thkns2 
                    if (water_tw .eq. 0.0 .and. ice_tw .gt. 0.0) then 
                        Tsoill(1) = Tsoill(1)+(2.*3600.*condu_ice/shcap_ice*10000.*(Tice-Tsoill(1))/thkns1 &
                                    &     -temph2)/(thkns1+thkns2) 
                    else 
                        Tsoill(1) = Tsoill(1)+(2.*3600.*condu_water/shcap_water*10000.*(Twater-Tsoill(1))/thkns1 &
                                    &     -temph2)/(thkns1+thkns2) 
                    endif
                    drsdh = 0.0    ! temporarily set drsdh =0 for heat adjustment of soil       
                else   
                    Tsoill(1) = Tsoill(1)+(diff_s*(sftmp-Tsoill(1))/thkns1 &
                                &     -temph2)/(thkns1+thkns2)
                endif
                !!!!!  phase change in top soil       
                heat_excess = drsdh*(thd_t-Tsoill(i))+shcap(i)*thksl(i)*(Tsoill(i)-thd_t)/360000.         
                ice_incr    = heat_excess*3600./latent_heat_fusion/ice_density        
                !
                !! $$$$$$$$$$$
                !! $$$$$$$$$$$      !tuneice
                !                if (ice_incr .lt. 0.) then   
                !                   if (i .eq. 1.) then
                !                       if (liq_water(i) .le. ufw(i)) then
                !!                           ice_incr = 0. !ice_incr*0.1
                !                          ice_incr = ice_incr*frac_ice1
                !                       endif
                !    !               elseif (i .eq. 2) then
                !    !                   if (liq_water(i) .le. 0.063) then
                !    !                       ice_incr = 0.
                !    !                   endif
                !                   elseif ( i .gt. 1.) then
                !                       if (liq_water(i) .le. ufw(i)) then
                !!                           ice_incr = ice_incr*0. !0.9
                !                           ice_incr = ice_incr*frac_ice2
                !                       endif
                !                   endif
                !                endif
                !! $$$$$$$$$$$   
                !! $$$$$$$$$$$                           
                !!          
                inter_var = ice(i)   
                if (ice_incr .lt. 0.) then     ! freeze             
                    ice(i) = Amin1(liq_water(i)+inter_var,ice(i)-ice_incr)            
                else 
                    ice(i) = Amax1(ice(i)-ice_incr,0.0)              
                endif
                !! readjust energy and temp 
                heat_adjust = heat_excess-latent_heat_fusion*(inter_var-ice(i))*ice_density/3600.
                Tsoill(i)   = thd_t+heat_adjust/(shcap(i)*thksl(i)/360000.-drsdh)      
            else
                ! if ( i .gt. 9) then 
                !     temph2=0
                !     thkns2=500  ! boundary conditions, rethink
                ! endif
                if ( i .gt. 9) then 
                    temph2 = 0.00003
                    thkns2 = 500  ! boundary conditions, rethink
                endif            
                Tsoill(i)   = Tsoill(i)+(temph1-temph2)/(thkns1+thkns2)    
                heat_excess = shcap(i)*thksl(i)*(Tsoill(i)-thd_t)/360000.        
                ice_incr    = heat_excess*3600./latent_heat_fusion/ice_density         
                !! $$$$$$$$$$$
                !! $$$$$$$$$$$      !tuneice
                !                if (ice_incr .lt. 0.) then   
                !                   if (i .eq. 1.) then
                !                       if (liq_water(i) .le. ufw(i)) then
                !!                           ice_incr = 0. !ice_incr*0.1
                !                          ice_incr = ice_incr*frac_ice1
                !                       endif
                !    !               elseif (i .eq. 2) then
                !    !                   if (liq_water(i) .le. 0.063) then
                !    !                       ice_incr = 0.
                !    !                   endif
                !                   elseif ( i .gt. 1.) then
                !                       if (liq_water(i) .le. ufw(i)) then
                !!                           ice_incr = ice_incr*0. !0.9
                !                           ice_incr = ice_incr*frac_ice2
                !                       endif
                !                   endif
                !                endif
                !! $$$$$$$$$$$   
                !! $$$$$$$$$$$                   
                !            
                
                inter_var = ice(i) 
                if (ice_incr .lt. 0.) then     ! freeze             
                    ice(i) = Amin1(liq_water(i)+inter_var,ice(i)-ice_incr)             
                else 
                    ice(i) = Amax1(ice(i)-ice_incr,0.0)              
                endif         
                !! readjust energy and temp 
                heat_adjust = heat_excess-latent_heat_fusion*(inter_var-ice(i))*ice_density/3600.
                Tsoill(i)   = thd_t+heat_adjust/(shcap(i)/360000.*thksl(i))
            endif
           
            if (ABS(tsoill_pre -tsoill(i)) .gt. 5. ) Tsoill(i)=tsoill_pre
            TEMPH1 = TEMPH2
            THKNS1 = THKNS2
            DIFSV1 = DIFSV2        
        enddo
        ! write(*,*)"new_FILDCP: ", FILDCP
        ! write(*,*)"new_thksl: " , thksl
        ! write(*,*)"new_ice: "   , ice
        ! write(*,*)"ice_incr:", ice_incr,heat_excess,latent_heat_fusion,ice_density
        ! write(*,*)"test_new1:", shcap
        ! write(*,*)"test_new2:", Tsoil
        ! write(*,*)"test_new: ", drsdh, thd_t
        ! write(*,*)"test_shcap_wcl: ", wcl
        ! write(*,*)"test_shcap_air: ", shcap_air
        ! write(*,*)"test_shcap_water: ", shcap_water, shcap_ice, shcap_soil
        ! write(*,*)"trace_G1: ", condu
        ! write(*,*)"trace_wcl: ", condu_air, liq_water,condu_water, ice, condu_ice, FILDCP, condu_soil 
        ! write(*,*)"trace_G2: ", condu_snow, sftmp, tsoil
        testout(1)    = tsoill_0       ! testout(1)=tsnow
        testout(2:11) = tsoill(1:10) 
        ! write(*,*)"test_testout: ", tsoill_0, sftmp, Tsoill(1), Tsnow, (Tsoill(1)+Tsnow)/2.  !Tsoill
        deallocate(depth_z)
        return 
    end subroutine Tsoil_simu
    
    subroutine methane()  !update single value of Rh_pools,Tsoil,zwt,wsc 
        !                                                           in a hourly loop when MEMCMC=1
        !******************************************************************************************************************************
        !****!introduce variables and constants used in this subroutine
        !****************************************************************************************************************************** 
        !     set soil layers
        !****************************************************************************************************************************** 
        implicit none
        ! integer i,MEMCMC
        ! integer i
        ! integer,parameter :: nlayers=10       !use this statement to set the parameter value
        ! real zwt    
        real consum
        !****************************************************************************************************************************** 
        !     set values for MEMCMC
        !******************************************************************************************************************************       
        ! integer,parameter :: miterms=17
        integer,parameter :: miterms=29         ! modified based on Ma et al., 2022
        integer,parameter :: ilines=9000      
        !******************************************************************************************************************************
        !     CH4 Production      
        !******************************************************************************************************************************      
        !      real Rhetero
        real Rh_resp(nlayers)!,Rh_pools(5),Rh_h,ProCH4(nlayers),Pro_sum
        ! real r_me         !release ratio of CH4 to CO2
        ! real Q10pro
        real fSTP(nlayers)         !CH4 production factor of soil temperature
        real vt,xt
        real Tmax_me!,Tpro_me
        real fpH          !CH4 production factor of soil pH
        real fEhP         !CH4 production factor of soil redox potential
        ! real FRLEN(nlayers)        !fraction of root in each layer
        real FRLEN_PMT(nlayers)     ! modified based on Ma et al., 2022
        ! real Tsoil
        !****************************************************************************************************************************** 
        !     CH4 Oxidation
        !******************************************************************************************************************************      
        !      real CH4(nlayers),CH4_V(nlayers+1)          !both are CH4 concentration: CH4(nlayers)unit gC/m2, CH4_V(nlayers) unit g C/ m3
        ! real CH4(nlayers),CH4_V(nlayers)          !both are CH4 concentration: CH4(nlayers)unit gC/m2, CH4_V(nlayers) unit g C/ m3
        ! real wsc(nlayers)      
        ! real OxiCH4(nlayers)!,Oxi_sum       !CH4 oxidation
        real Omax_layers(nlayers)!,Omax       !maximum oxidation rate
        real kCH4_layers(nlayers)!,kCH4       !system specific half saturation constant
        real Q10oxi
        real fCH4(nlayers)         !CH4 oxidation factor of CH4 concentration
        real fSTO(nlayers)         !CH4 oxidation factor of soil temperature
        real fEhO         !CH4 oxidation factor of soil redox potential
        ! real Toxi
        !****************************************************************************************************************************** 
        !     CH4 Diffusion
        !******************************************************************************************************************************      
        real Deff(nlayers)     !CH4 effective diffusivity !!used in mineral soil  v1.1 
        real D_CH4_a           !CH4 diffusion coefficient in air  !unit cm2 s-1   diffusivity of CH4 in air
        real D_CH4_w           !CH4 diffusion coefficient in water  !unit cm2 s-1   diffusivity of CH4 in water
        ! real phi          !soil porosity  also used in water table part
        real fwater(nlayers),fair(nlayers)
        real D_CH4_soil(nlayers),D_CH4_soil_a(nlayers),D_CH4_soil_b(nlayers)      !!used in organic peat soil  v1.2
        real fcoarse      !relative volume of coarse pores depending on soil texture  Zhuang 2004
        real ftort        !tortuosity coefficient with a value of 0.66    Walter and Heimann 2000
        !suggesting that the distance covered by diffusion is about two thirds of the length of the real average path
        real SAND         !relative contents of sand (%) in the soil
        real PVSAND       !relative volume of coarse pores in sandy soils     set to 0.45     value from Walter 2001
        real SILT         !relative contents of silt (%) in the soil
        real PVSILT       !relative volume of coarse pores in silty soils     set to 0.20     value from Walter 2001
        real CLAY         !relative contents of clay (%) in the soil
        real PVCLAY       !relative volume of coarse pores in clayish soils     set to 0.14   value from Walter 2001
        ! real DEPTH(10)        !depth in soil  will define it inside this subroutine again      resolution 100mm 200mm
        ! real THKSL(10)        !will define it inside this subroutine again  
        ! real Fdifu(nlayers+1)
        ! real Fdifu(nlayers)
        real CH4_atm      !concentration of CH4 in atmosphere     seen as 0 cause the value is too low someone use 0.076
        ! real simuCH4      !simulated CH4 emission
        ! ***********  Boundary condition parameters    *************      
        real ScCH4                                 !Schmidt numbers for methane Wania
        real pistonv                               !Piston velocity
        real Ceq                                   !equilibrium concentration of gas in the atmosphere
        real kHinv                                 !Henry's coefficient dependent variable on left side of equation, T is the independent variable
        real kH_CH4         !Henry's constant at standard temperature (CH4) Unit L atm mol-1
        real CHinv          !Coefficient in Henry's Law Unit K      
        real Tsta           !standard temperature Unit K
        real Ppartial       !CH4 partial pressure in air Unit atm
        real Cold           !last time step of ch4 concentration on first soil layer
        integer jwt    !index of the soil layer right above the water table (-)
        !****************************************************************************************************************************** 
        !     Ebullition 
        !******************************************************************************************************************************      
        real CH4_thre_ly(nlayers),EbuCH4(nlayers),Kebu !CH4_thre,
        ! real Ebu_sum_unsat,Ebu_sum_sat!,Ebu_sum          !sum value one dimension is enough 
        integer wtlevelindex
        real rouwater,V   ! local paras for ebullition concentration threshold calculation
        real mebu_out2(nlayers)
        !****************************************************************************************************************************** 
        !     Plant transport
        !******************************************************************************************************************************      
        real PlaCH4(nlayers)!,Pla_sum
        ! real LAIMIN,LAIMAX
        real Tgr,Tmat,fgrow,Pox,Kpla !Tveg,
        real rh_layert
        !****************************************************************************************************************************** 
        !******************************************************************************************************************************
        ! Yuan added for soil temp  
        ! logical do_soilphy
        ! real  tsoil_layer(11)
        !****************************************************************************************************************************** 
        !******************************************************************************************************************************
        !      MEMCMC=0   ! note here, any changes here result unexpected bug 
        ! write(*,*)"here0: ", CH4(1), OxiCH4(1), Fdifu(1), ProCH4(1)  
        ! write(*,*)"check0: ", CH4(1)
        Rh_h=Rh_pools(1)+Rh_pools(2)+Rh_pools(3)+Rh_pools(4)+Rh_pools(5)  !hourly Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
        tsoil_layer = testout
        ! FRLEN = (/0.75,0.2,0.02,0.015,0.005,0.0,0.0,0.0,0.0,0.0/)             
        ! FRLEN = (/0.1,0.25,0.25,0.2,0.1,0.05,0.025,0.015,0.005,0.005/)
        ! FRLEN = (/0.05,0.1,0.1,0.1,0.15,0.25,0.25,0.0,0.00,0.00/)
        FRLEN = (/0.75,0.2,0.02,0.02,0.01,0.0,0.0,0.0,0.0,0.0/)
        FRLEN_PMT = (/0.75,0.2,0.02,0.02,0.01,0.0,0.0,0.0,0.0,0.0/)

        thksl = (/10.,10.,10.,10.,10.,20.,20.,20.,20.,20./)
        simuCH4 = 0.0                 ! v1.2 
        rh_layert=0.0
        do i = 1, nlayers
            ! !!!!!!!put it out of the subroutine
            ! !****************************************************
            ! !* Rh weighed according to the distribution of root *
            ! !****************************************************
            ! if (i .LE. 3) then                                 ! the empirical method used here is from CLM4.5
            !     Rh_resp(i)= 0.5*Rh_h*FRLEN(i)+((0.5*Rh_h)/0.3)*0.1   
            !     ! Rh(h,i)Rh produced by each layer per hour  unit should be g C m-2 h-1 
            ! else                                               ! i*10: depth of ith soil layers
            !     Rh_resp(i)= 0.5*Rh_h*FRLEN(i)
            ! endif
            ! ! Rh(i) = Rh(i) + OxiCH4(i)*(11/4)
            ! --------------------------------------------------------
            ! Modified based on Ma et al., 2022
            !        *************   **********    **************    *****************    *********************
            ! the empirical method used here is from CLM4.5
            Rh_resp(i)= 0.5*Rh_h*FRLEN(i) + 0.5*Rh_h*(thksl(i)/depth(10))
            rh_layert=rh_layert+Rh_resp(i)
        enddo   
           
        !****************************************************************************************************************************** 
        !******************************************************************************************************************************            
    
        !****************************************************          
        !A. methane production     hourly  gC m-2 hour-1
        !Methane production is modeled as an anaerobic process that occurs in the saturated zone of the soil profile ZHUANG
        !****************************************************
        !Rh_h=Rh_pools(1)+Rh_pools(2)+Rh_pools(3)+Rh_pools(4)+Rh_pools(5)  !hourly Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
        !r assignment
        ! r_me=0.3      !find in parafile
        Tmax_me=45.0
        ! Tpro_me=10.0
        ! Q10pro=3.0    !find in parafile
        do i = 1,nlayers          
            if (do_soilphy) then
                if (tsoil_layer(i+1) .lt. 0.0) then
                    fSTP(i) = 0.0
                else if (tsoil_layer(i+1) .gt. Tmax_me) then
                    fSTP(i) = 0.0
                else if (tsoil_layer(i+1) .ge. 0.0 .and. tsoil_layer(i) .le. Tmax_me) then
                    fSTP(i) = Q10pro**((tsoil_layer(i+1)-Tpro_me)/10)        !Tsoil is the only variable
                endif
            else 
                if (Tsoil .lt. 0.0) then
                    fSTP(i) = 0.0
                else if (Tsoil .gt. Tmax_me) then
                    fSTP(i) = 0.0
                else if (Tsoil .ge. 0.0 .and. Tsoil .le. Tmax_me) then
                    fSTP(i) = Q10pro**((Tsoil-Tpro_me)/10)        !Tsoil is the only variable
                endif
            endif
        enddo
        ! fpH assignment
        fpH=1.0
        ! fEhP assignment
        fEhP=1.0
        
        depth(1)=10.0                                  !calculate soil depth unit cm
        do i=2,nlayers
            depth(i)=depth(i-1)+THKSL(i)
        enddo
          
        Pro_sum=0.0
        do i = 1,nlayers
            ! (depth(i)*10)                   !convert unit from cm to mm
            ! (THKSL(i)*10)                   !convert unit from cm to mm convert the unit in each of the equations
            if ((depth(i)*10) .le. -zwt) then
                ProCH4(i)=0.0
            else
                if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
                    ProCH4(i)=Rh_resp(i)*r_me*fSTP(i)*fpH*fEhP*(((depth(i)*10.0)-(-zwt))/(THKSL(i)*10.0))     ! *percent
                elseif (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
                    ProCH4(i)=Rh_resp(i)*r_me*fSTP(i)*fpH*fEhP
                endif
            endif
            Pro_sum=Pro_sum+ProCH4(i)
        enddo
    
        !**************************************************
        !Add CH4 production to CH4 pool    (gC layer -1)=(gC m-2)
        !**************************************************
        do i=1,nlayers
            CH4(i) = CH4(i) + ProCH4(i)
        enddo  
        ! write(*,*)"check3: ","CH4(1)", CH4(1),ProCH4(1),"ch4(2)",ch4(2)
        ! END OF METHANE PRODUCTION
        ! ********************************************************************************************************************
        ! B. methane oxidation      hourly  unit gC m-2 h-1     !!!!!!!!!!!method of CLM and Zhuang!!!!!!!!!!!!
        ! Methane oxidation is modeled as an aerobic process that occurs in the unsaturated zone of the soil profile ZHUANG
        ! ********************************************************************************************************************
        ! fSTO assignment
        ! ***************          
        Q10oxi=2.0      !Zhu 2014 results from previous studies  unit 1  also used by zhang
        ! Toxi=10.0       !Zhuang 2004 table1 Boreal Forest Wetland
        do i=1,nlayers
            if (do_soilphy) then
                fSTO(i)=Q10oxi**((tsoil_layer(i+1)-Toxi)/10.0)
            else
                fSTO(i)=Q10oxi**((Tsoil-Toxi)/10.0)
            endif
        enddo
        ! fEhO assignment
        fEhO=1.0        !Walter 2000  did not consider it, equal to value of 1
    
        ! Omax assignment
        ! ***************
        Oxi_sum=0.0
        do i = 1,nlayers
            ! Omax=1.5
            ! Omax=15.0  !!find in parafile Zhuang 2004 table1 Boreal Forest Wetland μmol L-1 h-1 system specific maximum oxidation coefficient
            ! convert the unit of Omax from μmol L-1 h-1 to gC m-2 h-1
            ! /1000,000 to get mol
            ! *12 cmass to get gC
            ! *1000 to get from dm-3(L) to m-3
            ! *(wsc*0.001) to get unit of omax_layers from m-3 to m-2     !caution that wsc unit is mm
            ! ** w  /   (w/t)           CLM used 
            Omax_layers(i)=(Omax/(1000000))*12*1000*(wsc(i)*0.001)     !convert the unit of Omax from μmol L-1 h-1 to gC m-2 h-1
            ! Omax_layers(i)=(Omax/(1000000))*12*1000*(THKSL(i)*10.0)*0.001     !modified on 11/27/2016 no sig change in oxidation and emission
            ! in unsaturated part of oxidation in CLM, they used the Omax/10 but did not expained why   they also /water volume
            ! fCH4 assignment
            ! ***************          
            ! kCH4=5.0     !!find in parafile Zhuang 2004 range between 1 and 66.2μmol L-1 system specific half saturation constant  1.0e
            ! convert the unit of kCH4 from μmol L-1 to gC m-2
            ! Jian: add limit
            if (wsc(i).le.0) then
                kCH4_layers = 0.
            else
                kCH4_layers(i)=(kCH4/(1000000))*12*1000*(wsc(i)*0.001)    !convert the unit of kCH4 from μmol L-1 to gC m-2
            endif
                ! then calculate fCH4 with CH4(i) and kCH4_layers(i) 
            if ((kCH4_layers(i)+CH4(i)) .eq. 0)then
                fCH4(i)=0.
            else
                fCH4(i)=CH4(i)/(kCH4_layers(i)+CH4(i))   !  CH4 concentration factor
            endif

            if ((depth(i)*10.0) .le. -zwt) then                !unit of Omax: gC m-2 h-1
                OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO(i)*fEhO!*0.1      !wrong:*(THKSL(i)/1000)!mm to m account for the thickness
                ! OxiCH4(i)=CH4(i)*0.001
            else
                if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
                    if (i .eq. 1) then
                        OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO(i)*fEhO*((-zwt)/(THKSL(i)*10.0))
                    else
                        OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO(i)*fEhO*(((-zwt)-(depth(i-1)*10.0))/(THKSL(i)*10.0))      !  *percent
                    endif
                    ! OxiCH4(i)=CH4(i)*0.001*(((-zwt)-(depth(i-1)*10.0))/(THKSL(i)*10.0))
                else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
                    OxiCH4(i)= 0.0
                endif
            endif            
            if (OxiCH4(i) .gt. CH4(i)) then
                OxiCH4(i)=CH4(i)
            endif  
            Oxi_sum=Oxi_sum+OxiCH4(i)   
        enddo 
        ! WRITE(*,*)"test oxiCH4: ",Omax_layers(1),fCH4(1),fSTO(1),fEhO,zwt,thksl, depth(1),depth(2)
        ! WRITE(*,*)"test oxiCH4: ",iforcing,fCH4(1),CH4(1),kCH4_layers(1),kCH4,wsc(1),rain
    
        !*******************************************************************
        !minus CH4 oxidation from CH4 pool     
        !*******************************************************************
        do i=1,nlayers
            CH4(i) = CH4(i) - OxiCH4(i)               !minus CH4 oxidation from CH4 pool
            ! if (wsc(i) .lt. (THKSL(i)*10.0)*WILTPT_x) then
            if (wsc(i) .le. 0) then
                CH4_V(i) = 0.
            else
                CH4_V(i) = CH4(i)/(wsc(i)*0.001)          !convert concentration from gC/m2 to gC/m3
            endif                                        !CH4_V(i) can be used for DA with observation data in soil layers
        enddo
        ! write(*,*)"here2: ", CH4(1), OxiCH4(1), Fdifu(1), ProCH4(1)  
        ! END OF METHANE OXIDATION
          
        ! ****************************************************
        ! C. methane diffusion
        ! ****************************************************
        ! Parameters assignment 
        D_CH4_a=0.2                             !unit cm2 s-1   D_CH4_a is the molecular diffusion coefficient of methane in air
        D_CH4_a=(D_CH4_a/10000.0)*3600.0        !unit m2 h-1
        D_CH4_w=0.00002                         !unit cm2 s-1   D_CH4_a is the molecular diffusion coefficient of methane in water
        D_CH4_w=(D_CH4_w/10000.0)*3600.0        !unit m2 h-1          
        ftort=0.66                              !tortuosity coefficient with a value of 0.66    Walter and Heimann 2000
        ! parameters for fcoarse algorithm      
        SAND=0.4          !   %   SPRUCE site value    0.4
        SILT=0.4          !   %   SPRUCE site value   0.4
        CLAY=0.2          !   %   SPRUCE site value   0.2
        PVSAND=0.45       !relative volume of coarse pores in sandy soils       set to 0.45     value from Walter 2001 zhuang
        PVSILT=0.20       !relative volume of coarse pores in silty soils       set to 0.20     value from Walter 2001 zhuang
        PVCLAY=0.14       !relative volume of coarse pores in clayish soils     set to 0.14     value from Walter 2001 zhuang  
        fcoarse=SAND*PVSAND+SILT*PVSILT+CLAY*PVCLAY
        CH4_atm=0.076       !unit umol L-1
        ! CH4_atm=0.0       !unit umol L-1      
        ! ******************************************************************************************************
        ! * Peat soil solution for diffusion coefficient: Equations for D_CH4_soil *         v1.2    Millington and Quirk Model
        ! ******************************************************************************************************
        do i=1,nlayers
            fwater(i) = wsc(i)/(THKSL(i)*10)      
            fair(i) = phi-fwater(i)
                        
            D_CH4_soil_a(i) = (((fair(i))**(10/3))/((phi)**2))*D_CH4_a
            D_CH4_soil_b(i) = D_CH4_W
            if (fair(i) .ge. 0.05) then
                D_CH4_soil(i) = D_CH4_soil_a(i)
            else
                D_CH4_soil(i) = D_CH4_soil_b(i)
            endif
            ! D_CH4_soil(i) = ge(fair,0.05)*D_CH4_soil_a(i) + lt(fair,0.05)*D_CH4_soil_b(i)        
                    
            ! Here I divided into saturated layer and unsaturated layer conditions because in most cases fair is > 0.05 and that there might be too much diffusion v1.2
            ! or maybe I can adjust the value of threshold 0.05 to around 0.08 as in most cases fwater=0.88 fair=0.07
            ! if (zwt .ge. 0.0) then                                  !when water table is above the soil surface
            !     Deff(i) = D_CH4_W
            ! elseif (zwt .lt. 0.0) then                                  !when water table is below the soil surface
            !     if ((depth(i)*10.0) .le. -zwt) then               !acrotelm layers
            !         Deff(i) = D_CH4_soil(i)
            !     elseif (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then       !partly acrotelm layer
            !         Deff(i) = D_CH4_soil(i)
            !     elseif (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then   !catotelm layers
            !         Deff(i) = D_CH4_W
            !     endif
            ! endif
            ! in this case diffusion should be more
            Deff(i) = D_CH4_soil(i)
        enddo 
    
            
        ! ******************************************************************************************************
        ! * Mineral soil solution for diffusion coefficient: Equations for D_CH4_soil *         v1.1   Three-porosity-model
        ! ******************************************************************************************************
        ! do i = 1,nlayers
        !     fwater = wsc(i)/(THKSL(i)*10)
        !     fair = phi-fwater
        !     fwater = 0.68         ! switch on when testing the effect of fwater on diffusion 0.6 crash 0.7fine  02172017
        !     Deff(i) = D_CH4_a*fcoarse*ftort*phi*(phi-fwater)+D_CH4_w*fwater           !
        ! enddo
    
          
        !convert the unit of CH4_atm from μmol L-1 to gC m-3
        !/1000,000 to get mol
        !*12 cmass to get gC
        !*1000 to get from dm-3(L) to m-3
        CH4_atm = (CH4_atm/1000000)*12*1000
              
        ! Fdifu(1) = Deff(1)*(CH4_V(1)-CH4_atm)/(THKSL(1)*0.01)         !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        ! New improvement 2017: Boundary condition     
        kH_CH4 = 714.29
        CHinv = 1600.0
        Tsta = 298.15

        ! ----------------------------------------------------------------------------------------
        ! ! Ppartial = 1.7E-6
        ! Ppartial = 1.7E-20 
          
        ! ScCH4 = 1898 - 110.1*Tsoil + 2.834*Tsoil**2 - 0.02791*Tsoil**3
        ! pistonv = 2.07 * (ScCH4/600)**(-1/2)
        ! kHinv = kH_CH4 /((exp(CHinv*(1/(Tsoil+273.15)-1/Tsta))))
        ! ! write (*,*) kHinv,pistonv
        ! Ceq = Ppartial / kHinv    ! Ceq: mol L-1   p_partial: atm  kHinv：L atm mol-1
        ! ! Fdifu(1) = Deff(1)*(CH4_V(1)-CH4_atm)/(THKSL(1)*0.01)         !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        ! Fdifu(1) =  pistonv * (CH4_V(1) - Ceq)
        ! ! Fdifu(1) =  -pistonv * (CH4_V(1) - Ceq)
        ! ! if (zwt .ge. -100.0) then
        ! ! Fdifu(1) = Deff(1)*(CH4_V(1)-CH4_atm)/(THKSL(1)*0.01)         !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        ! ! Fdifu(1) = - pistonv * (CH4_V(1) - Ceq)                         !switch on/off
        ! ! else
        ! ! Fdifu(1) = Deff(1)*(CH4_V(1)-CH4_atm)/(THKSL(1)*0.01)         !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        ! ! endif
             
        ! ! if (zwt .ge. -200 .and. zwt .le. 100.0) then
        ! ! Fdifu(2) = - pistonv * (CH4_V(2) - Ceq) 
        ! ! else
        ! ! Fdifu(2) = Deff(2)*(CH4_V(2)-CH4_V(1))/(THKSL(2)*0.01)         !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        ! ! endif   
    
        ! do i = 2,nlayers                                  !refer to flux from layer ii to ii-1 
        !     Fdifu(i)= Deff(i)*(CH4_V(i)-CH4_V(i-1))/(THKSL(i)*0.01)      !the unit of Fdifu is gC/m-2/h
        ! enddo
        ! --------------------------------------------------------------------------------------
        ! Modified based on Ma et al., 2022
        !
        Ppartial = 1.7E-6     !unit atm  partial pressure of methane in atmosphere
        ! *****  before Tsoil layers were added
        if (.not. do_soilphy) then
            ScCH4 = 1898 - 110.1*Tsoil + 2.834*Tsoil**2 - 0.02791*Tsoil**3
            pistonv = 2.07 * (ScCH4/600)**(-1/2)  !n=-1/2   pistonv unit=m/s in Wania's paper but cm/h in the code
            pistonv = pistonv*0.01   ! convert from cm/h to m/h
            kHinv = kH_CH4 /(exp(CHinv*(1.0/(Tsoil+273.15)-1.0/Tsta)))
            Ceq = Ppartial / kHinv    ! Ceq: mol L-1   p_partial: atm  kHinv：L atm mol-1
            Ceq = Ceq * 12 * 1000    ! Ceq mol/L to g/m3
        endif
        ! ********************************
        if (zwt .ge. -100) then   !index j, water table is right below the jth layer
            jwt=0.
        elseif (zwt .ge. -500.0) then  !layer 1-5
            ! jwt=int(-zwt/100)-1
            jwt=int(-zwt/100)
        else
            ! jwt=int((-zwt-500)/200+5)-1
            jwt=int((-zwt-500)/200+5)
        endif
        ! write (*,*) 'jwt',jwt,'zwt',zwt
        ! *******2*************************
        ! write(*,*)"here4: ", CH4(1), CH4(2), Fdifu(1), Fdifu(2)
        do i=1,nlayers
        
        ! *****  after Tsoil layers were added
            if (do_soilphy) then
                ScCH4 = 1898 - 110.1*tsoil_layer(i+1) + 2.834*(tsoil_layer(i+1))**2 - 0.02791*(tsoil_layer(i+1))**3
                pistonv = 2.07 * (ScCH4/600)**(-1/2)  !n=-1/2   pistonv unit=m/s in Wania's paper but cm/h in the code
                pistonv = pistonv*0.01   ! convert from cm/h to m/h
                kHinv   = kH_CH4 /((exp(CHinv*(1.0/(tsoil_layer(i+1)+273.15)-1.0/Tsta))))
                Ceq     = Ppartial / kHinv    ! Ceq: mol L-1   p_partial: atm  kHinv：L atm mol-1
                Ceq     = Ceq * 12 * 1000    ! Ceq mol/L to g/m3
            endif
        ! ********************************
            if (i .eq. 1 .and. jwt .ge. 1) then
                ! if (wsc(i) .lt. (THKSL(i)*10.0)*WILTPT_x) then
                if (wsc(i) .le. 0) then
                    Fdifu(i) = 0.
                else
                    Fdifu(i) = Deff(i)*(CH4_V(i)-CH4_atm)/(wsc(i)*0.01)
                endif
                ! write (*,*) 'jwt=1',CH4_V(1), wsc(1), wsc(1)*0.001,CH4(1), CH4(1)/(wsc(1)*0.001)
            elseif (i .eq. 1) then !.and. jwt .eq. 0
                Cold = CH4_V(i)
                ! if (wsc(i) .lt. (THKSL(i)*10.0)*WILTPT_x) then
                if (wsc(i) .le. 0) then
                    CH4_V(i) = 0.
                    Fdifu(i) = 0.
                else
                    CH4_V(i) = Ceq + (Cold-Ceq)*exp(-pistonv/(wsc(i)*0.001)) !pistonv/wsc m/m unit=1
                    Fdifu(i) = (Cold-CH4_V(i))*(wsc(i)*0.001)
                endif
                ! write (*,*) 'jwt=0',CH4_V(1)
            elseif (i .le. nlayers .and. i .ne. jwt+1) then
                ! if (wsc(i) .lt. (THKSL(i)*10.0)*WILTPT_x) then
                if (wsc(i) .le. 0) then
                    Fdifu(i)=0.
                else
                    Fdifu(i)= Deff(i)*(CH4_V(i)-CH4_V(i-1))/(wsc(i)*0.01)
                endif
            elseif (i .le. nlayers .and. i .eq. jwt+1) then
                Cold = CH4_V(i)
                ! if (wsc(i) .lt.(THKSL(i)*10.0)*WILTPT_x) then
                if (wsc(i) .le. 0) then
                    CH4_V(i) = 0.
                    Fdifu(i) = 0.
                else
                    CH4_V(i) = Ceq + (Cold-Ceq)*exp(-pistonv/(wsc(i)*0.001)) !pistonv/wsc m/m unit=1
                    Fdifu(i) = (Cold-CH4_V(i))*(wsc(i)*0.001)
                endif
            endif
        enddo
        ! write(*,*)"here3: ", CH4(1), CH4(2), Fdifu(1), Fdifu(2), Deff(1)*(CH4_V(1)-CH4_atm)/(wsc(1)*0.01), wsc(1)
        ! write(*,*)"here3: ", Deff(1)*(CH4_V(1)-CH4_atm)/(wsc(1)*0.01), wsc(1),Deff(1),CH4_V(1),CH4_atm
        ! write (*,*) 'Fdifu(1)',Fdifu(1),'CH4_V(1)',CH4_V(1),'Deff(1)',Deff(1),'Deff(2)',Deff(2),'wsc(1)',wsc(1),'wsc(2))',wsc(2)
        ! ********************************
        !! *******1*************************
        !       do i=1,nlayers
        !           if (i .eq. 1 .and. jwt .ge. 1) then
        !               Fdifu(i) = Deff(i)*(CH4_V(i)-CH4_atm)/(wsc(i)*0.01)
        !           elseif (i .eq. 1) then !.and. jwt .eq. 0
        !               Fdifu(i) =  pistonv * (CH4_V(i) - Ceq)
        !           elseif (i .le. nlayers .and. i .ne. jwt+1) then
        !               Fdifu(i)= Deff(i)*(CH4_V(i)-CH4_V(i-1))/(wsc(i)*0.01)
        !           elseif (i .le. nlayers .and. i .eq. jwt+1) then
        !               Fdifu(i) =  pistonv * (CH4_V(i) - Ceq)
        !           endif
        !       enddo
        !!!!! ********3************************
        !      Fdifu(1) =  pistonv * (CH4_V(1) - Ceq)
        !      do i = 2,nlayers                                  !refer to flux from layer ii to ii-1
        !          Fdifu(i)= Deff(i)*(CH4_V(i)-CH4_V(i-1))/(wsc(i)*0.01)      !the unit of Fdifu is gC/m-2/h
        !      enddo
        !!! ********************************
        !below I try to keep the CH4 flux no larger than the amount of CH4 that exist at the moment   V1.1 V1.2



        ! CH4_V(11) = CH4_V(10)
        ! Fdifu(11) = Deff(10)*(CH4_V(11)-CH4_V(10))/(THKSL(10)*0.01)       !MODIFIED ON 2017 inserted  switch depend on hypothesis: the bottom boundary is a no-flux boundary or the 11th layer concentration is 0
        ! below I try to keep the CH4 flux no larger than the amount of CH4 that exist at the moment   V1.1 V1.2
        ! do i=1,nlayers+1 ! Jian: CH4 just has nlayers, but it reachs nlayers+1
        do i=1,nlayers
            if (Fdifu(i) .gt. 0.0 .and. (Fdifu(i)) .gt. CH4(i)) then
                Fdifu(i)=0.999*CH4(i)
            endif
            ! if (Fdifu(i) .lt. 0.0 .and. (abs(Fdifu(i))) .gt. CH4(i-1)) then
            if (Fdifu(i) .lt. 0.0) then
                if (i>1) then
                    if ((abs(Fdifu(i))) .gt. CH4(i-1)) Fdifu(i)=-0.999*CH4(i-1)
                else
                    Fdifu(i)=0
                endif
            endif
            IF (i>1)THEN
                if (CH4(i-1) .EQ.0) Fdifu=0.
            else
                if (CH4(1) .Eq. 0) Fdifu(i)=0
            endif
        enddo
        ! write(*,*)"here11: ", CH4(1), Fdifu(1), Fdifu(2), CH4(2)
        ! CH4(1) = CH4(1) + (0.0+Fdifu(1))/(THKSL(1)*0.01)
        do i = 1,nlayers-1                                  !loop of time
            CH4(i) = CH4(i) + (Fdifu(i+1)-Fdifu(i))*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
            ! CH4(i) = CH4(i) - 0.1*CH4(i)
            if (CH4(i) .lt. 0.0) then                     ! this part need to be improved until deleted   V1.2
                CH4(i) = 0.0
            endif
        enddo  
        ! write(*,*)"here1: ", CH4(1), Fdifu(1), (Fdifu(2)-Fdifu(1))*1
        CH4(10) = CH4(10) - Fdifu(10)                                   !MODIFIED ON 07/25/2016
        if (CH4(10) .lt. 0.0) then                                    !defined the Fdifu(11) to be 0.0
            CH4(10)= 0.0                                              ! switch on/off
        endif
        ! CH4(10) = CH4(10) +(Fdifu(11)-Fdifu(10))*1                      !MODIFIED ON 05/04/2018
        ! if (CH4(10) .lt. 0.0) then                                    
        !     CH4(10)= 0.0                                              
        ! endif        
        simuCH4 = simuCH4 + (Fdifu(1)-0.0) 
        ! write(*,*)"check2: ","CH4(1)", CH4(1),(Fdifu(2)-Fdifu(1))*1,"Fdifu(2)", Fdifu(2),"Fdifu(1)",Fdifu(1)
        ! ********************************************************************************************************************      
        ! D. methane ebullition     !assume bubbles can reach the water table within 1 h&
                                    !& the bubbles is added to the methane concentration in the soil layer just above the wt
                                    !& and then diffused through layers   ??not correct
        ! this subroutine is modified on 02132017 by deleting the unsat from bubble and add unsat to concentration so as to increase diffusion
        ! just by searching "switch" you can switch from old to new mode by adding or deleting "!"
        ! modified threshold value to 100 for testing
        ! ********************************************************************************************************************
        !& and then diffused through layers
        !mechanisms selectable: ECT (concentration threshold) and EBG (bubble growth)
        !EBG is added as a subroutine on Nov 26th 2018, so that minimum changes to the original code
        ! Modified based on Ma et al., 2022
        Ebu_sum_unsat=0.0
        Ebu_sum_sat=0.0                                      !initial value
        Ebu_sum = 0.0

        rouwater = 1000.   !kg/m3
        g = 9.81; ! m/s2
        Rgas = 8.3145; ! m3 Pa/K/mol
        ! write(*,*)"Test CH4: ", CH4(1), ProCH4(1), oxiCH4(1)
        ! if ((iforcing>1000) .and. (ProCH4(1).eq.0)) stop
        if (do_EBG) then
            call ebullition_EBG(rouwater, mebu_out2)
            ! write(*,*)"test: ", iforcing, CH4
        else
            !use ECT mechanisms for ebullition
            Kebu=1.0                    !unit  h-1   rate constant               
        
        
            ! ---------------------------------------------------------------------------------------
            ! d     o i=1,nlayers
            !     CH4_thre=1000.0  !!find in parafile  !unit  umol L-1 according to Walter's 500-1000
            !     CH4_thre_ly(i)=(CH4_thre*1.0e-6)*12*1000*(wsc(i)*0.001)    !convert the unit of CH4_thre from µmol L-1 to gC m-2
            ! enddo
            ! ---------------------------------------------------------------------------------------
            ! Modified based on Ma et al., 2022
            do i=1,nlayers
                V = 0.05708 - 0.001545*MAX(0.0, tsoil_layer(i+1)) + 0.00002069*(MAX(0.0, tsoil_layer(i+1)))**2 ![m3 CH4 m-3 H2O in soil]
                ! V = 0.05708 - 0.001545*tsoil_layer(i+1) + 0.00002069*tsoil_layer(i+1)**2  ![m3 CH4 m-3 H2O in soil]
                ! **************************************************************************************************************************
                !         pwater = rouwater*g*(wsc(i)*0.001)    !kg/m3    convert wsc from % to m
                ! ******   water height not correct, change to this:
                if (depth(i) .le. (-zwt)*0.1) then
                    pwater(i) = rouwater*g*(depth(i)*0.01-(-zwt)*0.001)   ![kg m s-2 or N], pwater/m2 = N/m2 = Pa
                else
                    pwater(i) = 0.
                endif
                ! **************************************************************************************************************************
                CH4_thre = ((dpatm + pwater(i))*V/(Rgas*(tsoil_layer(i+1)+273.15)))*1000!-200  ! mol CH4 /m3 to mmol CH4 /m3
                ! CH4_thre=700.0  !!find in parafile  !unit  umol L-1 according to Walter's 500-1000
                CH4_thre_ly(i)=(CH4_thre*1.0e-6)*12*1000*(wsc(i)*0.001)    !convert the unit of CH4_thre from µmol L-1 to gC m-2
                if (CH4_thre .lt. 1300.) then
                    ! write (*,*) CH4_thre
                endif
            enddo
        
            if (zwt .ge. 0.0) then                                  !when water table is above the soil surface
                do i=1,nlayers
                    if (CH4(i) .gt. CH4_thre_ly(i)) then                  
                        EbuCH4(i)=Kebu*(CH4(i)-CH4_thre_ly(i))     !only if the concentration is larger than threshold
                    else !if (CH4(i) .le. CH4_thre_ly(i)) then
                        EbuCH4(i)=0.0
                    endif
                    Ebu_sum_sat=Ebu_sum_sat+EbuCH4(i)               !& the bubbles are directly added into CH4 efflux into atmosphere
                    CH4(i)=CH4(i)- EbuCH4(i)                        !& update the concentration at the end of this hour in each layers
                enddo
            endif
            ! write (*,*) CH4(1),CH4_thre_ly(1),EbuCH4(1),Ebu_sum_sat 
            if (zwt .lt. 0.0) then                                  !when water table is below the soil surface
                do i=1,nlayers
                    if ((depth(i)*10.0) .le. -zwt) then               !acrotelm layers
                        EbuCH4(i)=0.0
                        Ebu_sum_unsat=Ebu_sum_unsat+EbuCH4(i)         
                        CH4(i)=CH4(i)- EbuCH4(i) 
                    else
                        if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then       !partly acrotelm layer
                            wtlevelindex = i
                            if (CH4(i) .gt. CH4_thre_ly(i)) then                  
                                EbuCH4(i)=Kebu*(CH4(i)-CH4_thre_ly(i))!*(((depth(i)*10.0)-(-zwt))/(THKSL(i)*10.0))        ! * percent
                            else !if (CH4(i) .le. CH4_thre_ly(i)) then                     ??????????,??????????????????
                                EbuCH4(i)=0.0
                            endif 
                            CH4(i)=CH4(i)- EbuCH4(i)
                            Ebu_sum_unsat=Ebu_sum_unsat+EbuCH4(i)                ! !  modified by Mary on 02132017
                            ! Jian: not sure why wtlevelindex-1? 20230123
                            if (wtlevelindex>1)then ! Jian: maybe > 1????
                                CH4(wtlevelindex-1)=CH4(wtlevelindex-1)+EbuCH4(i)    !!!!!-1-!!!! !switch on in new mode should be added add burst bubbles below surface to diffusion modified by Mary on 02132017
                            end if
                            ! ：the problem is the resolution of soil layer is 10cm and EbuCH4(i) is directly added to the upper layer of boundary layer   02152017                    
                        else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then   !catotelm layers
                            if (CH4(i) .gt. CH4_thre_ly(i)) then                  
                                EbuCH4(i)=Kebu*(CH4(i)-CH4_thre_ly(i))
                            else !if (CH4(i) .le. CH4_thre_ly(i)) then
                                EbuCH4(i)=0.0
                            endif 
                            CH4(i)=CH4(i)- EbuCH4(i)    
                            ! Jian: also wtlevelindex -1 error. may be >1
                            ! if (wtlevelindex>1) then    
                            !     CH4(wtlevelindex-1)=CH4(wtlevelindex-1)+EbuCH4(i)     !!!!!-2-!!!! !switch on in new mode should be added     modified by Mary on 02152017
                            ! end if
                            ! Ebu_sum_unsat=Ebu_sum_unsat+EbuCH4(i)                  ! modified by Mary on 02132017
                            !                  wtlevelindex equal to the last layer (wt layer)) that give a value to it, in all the layers below the wt layer
                            if (wtlevelindex>1) then
                            CH4(wtlevelindex-1)=CH4(wtlevelindex-1)+EbuCH4(i)     !!!!!-2-!!!! !switch on in new mode should be added
                            endif
                            !modified by Mary on 02152017
                            Ebu_sum_unsat=Ebu_sum_unsat+EbuCH4(i)                  ! modified by Mary on 02132017
                        endif
                    endif
                enddo
            endif
          
            Ebu_sum= Ebu_sum_sat
            simuCH4=simuCH4+Ebu_sum_sat                         !& the bubbles are directly added into CH4 efflux into atmosphere
        endif    
        ! write(*,*)"check1: ", CH4(1), EbuCH4(1)
        ! WRITE(*,*)"Here1.5", ch4(1)
        ! write (*,*) Ebu_sum        
        ! ******************************************************************************************************
        ! E. plant mediated methane transportation      totoally used Walter's model also used by Zhuang et. al
        ! ******************************************************************************************************
        Kpla=0.01         !unit h-1
        ! Kpla=0.01         !unit h-1
        ! Tveg=0.3 ! a factor describing the quality of plant-mediated transport depend on the density of plant stands and plant types 
        ! 0 for boreal forest and 0.5 for tundra
        ! find in parafile !
        ! the Tsoil used here would be better if refer to the 20cm soil temperature after &
        ! & the accomplishment of soil heat dynamics module. according to Zhuang. however Walter used 50cm soil temp.
        Tgr=2.0               !unit degree Celsius if annual mean temp is below 5 (otherwise 7)
        Tmat=Tgr+10.0         !unit degree Celsius
        Pox=0.5               !50% of mediated methane are oxidised 
        ! define fgrow
        if (.not. do_soilphy) then
            if (Tsoil .lt. Tgr) then
                fgrow=LAIMIN
            else if (Tsoil .ge. Tgr .and. Tsoil .le. Tmat) then
                fgrow=LAIMIN+LAIMAX*(1-((Tmat-Tsoil)/(Tmat-Tgr))**2)
            else if (Tsoil .gt. Tmat) then
                fgrow=LAIMAX
            endif
        else
            if (tsoil_layer(3) .lt. Tgr) then
                fgrow=LAIMIN
            else if (tsoil_layer(3) .ge. Tgr .and. tsoil_layer(3) .le. Tmat) then
                ! fgrow=LAIMIN+(LAIMAX+2)*(1-((Tmat-Tsoil)/(Tmat-Tgr))**2)
                fgrow=LAIMIN+(LAIMAX)*(1-((Tmat-tsoil_layer(3))/(Tmat-Tgr))**2)
                ! WRITE (*,*) 'fgrow',fgrow,'LAIMIN',LAIMIN
            else if (tsoil_layer(3) .gt. Tmat) then
                ! fgrow=LAIMAX+2
                fgrow=LAIMAX
            endif
        endif  
        Pla_sum=0.0
        do i=1,nlayers
            PlaCH4(i)=Kpla*Tveg*FRLEN(i)*fgrow*CH4(i)*(1-Pox)         !not sensitive at all to this change, but better
            ! PlaCH4(i)=Kpla*Tveg*FRLEN(i)*fgrow*CH4(i)
            Pla_sum=Pla_sum+PlaCH4(i)
            ! CH4(i)=CH4(i)-PlaCH4(i)
            CH4(i)=CH4(i)-Kpla*Tveg*FRLEN_PMT(i)*fgrow*CH4(i)!PlaCH4(i)/(1-Pox)
            ! if (wsc(i) .lt. (THKSL(i)*10.0)*WILTPT_x) then
            if (wsc(i) .le. 0) then
                CH4_V(i) = 0.
            else
                CH4_V(i) = CH4(i)/(wsc(i)*0.001)          !convert concentration from gC/m2 to gC/m3
            endif
        enddo  
        simuCH4=simuCH4+Pla_sum
        consum=simuCH4+OxiCH4(1)+OxiCH4(2)+OxiCH4(3)+OxiCH4(4)+OxiCH4(5)+OxiCH4(6)+OxiCH4(7)+OxiCH4(8)+OxiCH4(9)+OxiCH4(10)
        ! write(*,*)"test CH4: ", iforcing, Kpla*Tveg*FRLEN_PMT(1)*fgrow*CH4(1), Kpla,Tveg,FRLEN_PMT(1),fgrow,CH4(1)
        ! write(*,*)"test : ",  CH4(1), OxiCH4(1), Fdifu(1), ProCH4(1),EbuCH4(1)
        ! IF(CH4(1) > 1000.) then
        !     write(*,*)"test CH4-1: ", iforcing, Kpla*Tveg*FRLEN_PMT(1)*fgrow*CH4(1), Kpla,Tveg,FRLEN_PMT(1),fgrow,CH4(1)
        !     stop
        ! endif 
        if (isnan(CH4(1))) then
            write(*,*)"test CH4: ", iforcing, Kpla*Tveg*FRLEN_PMT(1)*fgrow*CH4(1), Kpla,Tveg,FRLEN_PMT(1),fgrow,CH4(1)
            stop
        endif  
        return
    end subroutine methane

    ! *************************************************************************************
    ! subroutine EBG used by methane submodel
    subroutine ebullition_EBG(rouwater,mebu_out2)
        ! INPUT:
        !		CH4 = porewater methane in the whole layer, unit=gC/m2, size=nlayers,CH4_V(i) = CH4(i)/(wsc(i)*0.001)
        !		CH4_V(nlayers) = porewater methane concentration, unit=gC/m3, size=nlayers,CH4_V(i) = CH4(i)/(wsc(i)*0.001)
        !		simuCH4 = methane efflux, added up with diffusion, ebullition, and PMT
        !		nlayers = number of layers, 10 by default
        !		zwt = water table level, unit = mm, size=1, below surface when zwt<0
        !		dpatm = atm pressure, input from the atmospheric forcing file, unit=Pa, estimated value from CNCEPT, used by CLM
        !		tsoil_layer = soil temperature of different layers, unit=Celsius, size = nlayers
        !       wsc = volumetric content of soil water in a layer, unit=mm, size = nlayers
        !		THKSL = thickness of the (i)th soil layer in the model, unit=cm, size = nlayers
        !     	depth = depth of the (i)th soil layer, unit=cm, size = nlayers
        !		rouwater = density of water, 1000. kg/m3
        !		g = 9.81 ! m/s2
        !		Rgas = 8.3145 ! m3 Pa/K/mol
        !OUTPUT:
        !		Ebu_sum = sum of ebullited methane finally get into the atm, size = 1,unit=gC/m2
        !		EbuCH4(nlayers) = ebullited methane each time step from one layer, size = nlayers,unit=gC/m2
        !		Ebu_sum_sat = ebullited methane each time step from one layer when water table is higher than soil surface, added
            !up to calculate Ebu_sum
        !		Ebu_sum_unsat = ebullited methane each time step from one layer when water table is below the soil surface, add to
            !CH4 and diffused up

        implicit none
        integer mwt,ind
        ! INPUT parameters:
        real rouwater!,g,Rgas
        ! integer, parameter :: nlayers=10
        ! real CH4(nlayers),CH4_V(nlayers),simuCH4,zwt,dpatm,tsoil_layer(11),wsc(nlayers),THKSL(nlayers),depth(nlayers)
        ! OUTPUT parameters:
        ! real Ebu_sum,Ebu_sum_sat,Ebu_sum_unsat,EbuCH4(nlayers)
        ! INTERNAL parameters:
        ! xP: amount of x in each layer [mol], really [mol] not [mol m-3]
        ! real Vp(nlayers),
        real Vtmp(nlayers)	!total volume of the bubbles in a layer, size(nlayers), unit m3/layer
        ! real methanebP(nlayers) !amount of CH4 moles in bubbles, size(nlayers), unit mol/layer, not concentration, not m-3
        ! real methaneP(nlayers)  !amount of CH4 moles in porewater, size (nlayers), unit mol/layer  CH4(i) gC/layer
        real mrateP(nlayers)	!The rate at which gas transfer with bubbles modifies gas at each layer, [mol s-1 layer-1]
        real met_db(nlayers)	!calculated in the exceeded bubble module, pore water CH4 concentration change rate due to
        !interaction with bubbles, unit = mol hr-1, size=nz
        real mCwP				!CH4 molar density in the water (mol m-3), size=1

        real dc
        real Dw,Hcp,mebu_out,r,randbub
        real, parameter :: peat_coeff_w = 0.9
        real, parameter :: Dw_298 = 1.5e-9 ! [m2 s-1]
        real mnv! CH4 molar density in the bubble (mol m-3)           ! Shuang mnv=cb in their paper

        ! real bubble_methane_tot	!amount of CH4 stored in bubbles
        real mebu_rate_out      !pore water CH4 concentration change rate due to interaction with bubbles, unit = mol s-1, size=nz
        real mebu_out2(nlayers) ! ebullition from the layers are temporarily saved here
        real Vout		! in the equation 8), the exceeded volume of methane in bubbles
        real nout		!nout = Febu in equation 8) EBG, unit mol. = bubble conc * volume
        !real gases_out(2)  !1 for ebu into air, 2 for ebu into the lowest air layer
        ! real presP(nlayers) ! total pressure of air and water in layers
        ! real pwater(nlayers)! water pressure in layers
        real tempP(nlayers)  ! Kalvin soil temperature
        real met_alphaP(nlayers) !methane solubility
        real a,b,c  !constants used to calculate methane solubility
        ! half-life of supersaturated dissolved methane
        !integer, parameter :: dp = selected_real_kind(15, 300)
        real, parameter :: ebu_hl = 1800. ! 30 minutes, 1800 s
        real, parameter :: k = log(2.0)/ebu_hl  ! turnover rate s-1
        ! real, parameter :: pi = 3.141592653589793
        !          !  1 = methane ebullition into air
        !          !  2 = methane ebullition total (when WTD is below peat surface, ebullition is
        !          !      released in the lowest air layer inside peat)
        !        integer, parameter :: mebuair = 1, mebutot = 2
        !*******************
        ! input parameters
        !*******************

        ! real f !threshold fraction of pore space filled with gas bubbles needed for ebullition or
                        ! CH4 mixing ratio in bubbles (mol mol-1)
        ! real Nbub ! Amount of bubbles in one model layer. This directly impacts the CH4 exchange
                            ! rate between pore water and bubbles
        ! real Vmaxfraction ! Maximum fraction of volume occupied by bubbles (10 %)
                                                    ! If this is exceeded, then the extra volume is released
        ! real bubprob ! probability that a bubble will get stuck at one layer
        !********************************************************************************


        !******************************************************************************
        !******************************************************************************
        ! inside the time loop starts here, add to methane module
        !******************************************************************************

        !****  #2. add mwt as index for the wt layer  *************
        if (zwt .ge. -100) then   !index mwt, water table is right at the mwt th layer
            mwt=1.
        elseif (zwt .ge. -500.0) then  !layer 1-5
            mwt=int(-zwt/100)+1
        else
            mwt=int((-zwt-500)/200+5)+1
        endif
        !       if (zwt .le. -100.) then
        !           write (*,*) 'mwt',mwt,'zwt',zwt
        !       endif
        !****  #2. end of mwt as index for the wt layer  *************


        !****  #3. update value, might not be necessary
        do i=1,nlayers
            methaneP(i) = CH4(i)/12		!##need to make sure CH4 is updated at the end of EBG
            !				  gC/layer  /12   unit molC/layer
            tempP(i) = tsoil_layer(i+1)+273.15 !unit Kelvin

        if (i .ge. mwt) then
            pwater(i) = rouwater*g*(depth(i)*0.01-(-zwt)*0.001)
        else
            pwater(i)  = 0.
        endif
        ! write(*,*)"before: ", CH4(1),dpatm
        presP(i) = dpatm + pwater(i)  ! unit Pa
        ! write (*,*) 'i',i,'presP(i)', presP(i),'pwater(i))',pwater(i)
        Vp(i) = methanebP(i) * Rgas * tempP(i)/(presP(i) * f)  !Vp is updated since methanebP is updated in the middle of the EBG
        !******************************************************************************
        !        write (*,*) 'Vp update1', Vp(i)

        !*******#4. bubble_methane_tot, bubble_methane_tot is accumulated total of all the 10 layers,
            !I assume air peat layer Vp=0, yes that's right ****
        !amount of CH4 stored in bubbles   !Shuang this is equation (10), update number of CH4 in bubbles
        !bubble_methane_tot = sum(f * presP(i)/(Rgas * (tsoil_layer(i+1)+273.15)) * Vp)
        bubble_methane_tot = bubble_methane_tot+f * presP(i)/(Rgas * tempP(i)) * Vp(i)
        !
        !******************************************************************************


        !*******#5. methane_solubility (alpha) ****************************************
        !calculate methane_solubility (alpha) using tempP,a,b,c,Rgas and Hcp
        !alpha: unit [mol(CH4,water) m(water)-3 mol(CH4,air)-1 m(air)3]
        !Hcp: unit [mol(CH4) m(H2O)-3 Pa-1]
        !Rgas = 8.3144621_dp ! [J mol-1 K-1]
        a = 1.3e-3
        b = 1700.
        c = 298.0
        ! Tang et al. 2010, eq. A4
        ! see also
        Hcp = a * exp(b * (1./tempP(i) - 1./c)) * 9.86923266716e-3
        met_alphaP(i) = Rgas * tempP(i) * Hcp
        !            write (*,*) 'i',i,'methaneP(i)',methaneP(i),'CH4(i)',CH4(i),'methanebP(i)',methanebP(i)
        enddo
        !******************************************************************************
        !    gases_out = 0
        mebu_rate_out = 0
        !*******#6. *******************************************************************
        ! methane_to_bubble - Transferring CH4 between bubbles and the surrounding pore water
        ! input
        ! INPUT:
        !       presP = pressure vector, unit=Pa, size=nz
        !       tempP = pore water temperature vector, unit=K, size=nz
        !       methaneP = CH4 pore water concentration vector, unit=mol/layer, size=nz
        !       met_alphaP = CH4 solubility (dimensionless Henry solubility) vector, unit=-, size=nz
        !       nz = amount of model layers, size=1
        !       geom = model geometry
        !       por = peat porosity, unitless, size=1
        !       Vp = bubble volume vector, unit=m3, size=nz
        ! OUTPUT:
        !       mrateP = mebu_rateP = methane_d = met_d = pore water CH4 concentration change rate due to interaction with bubbles,
        !unit = mol s-1, size=nz
        !   grows (and shrinks) bubbles in each layer based on Henry's law
        !  equation 6
        !    write (*,*) "1 stopping point"
        mrateP = 0
        !! if layers with water peat mixture exist
        do i = mwt,nlayers
        ! CH4 molar density in the water (mol m-3)
            if (wsc(i) .eq. 0) then
                mCwP = 0.
            else
                mCwP = methaneP(i) / ((wsc(i)*0.001)*1);
            endif
        ! concentration difference between pore water and bubble surface
        dc = (mCwP - met_alphaP(i) * f * presP(i)/(Rgas * tempP(i))) !  Shuang this is part 2 of equation (5)
        !        mol m-3
        !        write (*,*) 'dc i',i,'mCwP',mCwP,'methaneP(i))',methaneP(i)

        if (Vp(i) == 0) then										!Shuang equation (6)
        ! creating a bubble from the excess CH4
        if (dc > 0) then
        ! making sure that CH4 exchange FROM bubble TO water does not exist when there is no bubble
        mrateP(i) = -k * dc * ((wsc(i)*0.001)*1)
        !            write (*,*) 'Vp=0, dc>0 i',i, 'mrateP(i)',mrateP(i)
        !(geom % dzP(i) * por) thickness of the layer * porosity (*1 m2) = water volume
        !			mol s-1		s-1	 mol m-3	m3
        !when mrateP is negative it means transfer from water to bubble
        end if
        !"Nbub" no need to consider Nbub because mrateP is already the total exchange amount

        else
        ! growing the bubble only if it already exists

        ! radius of one bubble (m)
        r = (3.0/4.0 * Vp(i)/Nbub/pi)**(1.0/3.0)		! r updated with Vp, calculated from the last time step in EBG

        ! CH4 diffusion coefficient in water & peat mixture
        !          Dw = methane_D_water(tempP(i)) !this function outputs Dw
        !

        Dw = peat_coeff_w * Dw_298 * (tempP(i)/298.)
        ! if (temperature < 273.15) Dw = 0

        ! change in pore water CH4 due to mass transfer with one bubble   !Shuang part of equation (5)
        mrateP(i) = -4.0 * pi * r * Dw * dc
        !          write (*,*) 'mrateP(i)',mrateP(i),'r',r,'dc',dc
        !            write (*,*) 'Vp/=0 i',i, 'mrateP(i)',mrateP(i)
        ! Nbub bubbles
        mrateP(i) = mrateP(i) * Nbub			!mol s-1
        end if
        ! ***** end of #6. methane_to_bubble output is mrateP = mebu_rateP = methane_d=met_d
        !******************************************************************************

        !        write (*,*) '6 i',i,'methaneP(i)',methaneP(i),'CH4(i)',CH4(i),'methanebP(i)',methanebP(i),'mrateP(i)',mrateP(i)
        !******************************************************************************
        ! ***** #7. now we are calculating methanebP, amount of CH4 moles in bubbles, size(nlayers), unit mol
        ! removing/adding the same amount of CH4 to the bubbles that was added/removed
        ! to the pore water due to interaction with the bubbles
        !        write (*,*) '7 i',i,'Vp update2', Vp(i),'methanebP(i)',methanebP(i),'mrateP(i)',mrateP(i),'r',r,'dc',dc
        methanebP(i) = methanebP(i) - 3600 * mrateP(i) ! %% the unit of mrateP is mol s-1  %%%%%%%%%%%
        !								s in a hr   mol/s       because the model time step is a hour
        !										        bubble to water when mrateP is positive
        !        write (*,*) 'i',i,'methaneP(i)',methaneP(i),'methanebP(i)',methanebP(i),'mrateP(i)',mrateP(i)
        ! making sure that the concentrations do not go negative
        if (methanebP(i) < 0) then
        methaneP(i) = methaneP(i) + methanebP(i)
        methanebP(i) = 0
        end if
        !        write (*,*) '+bub i',i,'methaneP(i)',methaneP(i),'mrateP(i)',mrateP(i)
        ! *****************************************************************************
        ! updating bubble volumes
        Vp(i) = methanebP(i) * Rgas * tempP(i)/(presP(i) * f)
        ! *****************************************************************************
        !        write (*,*) 'i',i,'Vp update2', Vp(i),'methanebP(i)',methanebP(i)

        ! ***** #8. *********************************************************************
        ! pore water CH4 concentration change rate due to interaction with bubbles    negative mrateP: lost methane from water to bubbles
        methaneP(i) = methaneP(i) + 3600 * mrateP(i)		! %%%%%%%%%%%%%%%%%%%%%%
        !        write (*,*) '+mrateP i',i,'methaneP(i)',methaneP(i),'mrateP(i)',mrateP(i)
        mebu_rate_out = mebu_rate_out + 3600 * mrateP(i)	! total of all the layers
        !mebu_rate_out=pore water CH4 concentration change rate due to interaction with bubbles, unit = mol s-1, size=nz

        ! change rate of pore water CH4 due to exchange with bubble
        mebu_rate_out = mebu_rate_out/1  !1 hr for big_dt... it is only an output
        !written as mebu_rate_out in the bigstep module, changed to mebu_rateP, but never used, just an output
        ! *****************************************************************************
        enddo
        !    write (*,*) "2 stopping point"
        ! *****************************************************************************
        ! releasing part of the gas volume if it exceeds certain size
        ! INPUT:
        !       Vp = bubble volume vector, unit=m3, size=nz
        !       presP = pressure vector, unit=Pa, size=nz
        !       tempP = pore water temperature vector, unit=K, size=nz
        !       geom = model geometry
        !       por = peat porosity, unitless, size=1
        !       dt = model time step, unit=s,size=1
        !       nz = amount of model layers, size=1
        !       methanebP = CH4 bubble concentration vector, unit=mol, size=nz
        ! OUTPUT:
        !       met_db = pore water CH4 concentration change rate due to interaction with bubbles, unit = mol hr-1, size=nz
        !       mebu_out = CH4 released in bubbles to deepest air layer, unit = mol hr-1, size=1

        ! If bubble in a certain layer grows larger than Vmax then it is transported to the deepest air layer
        ! Currently considers only CH4
        do i = 1,nlayers
        met_db(i) = 0	!met_db: due to bubble
        mebu_out = 0	!mebu_out: one dimension
        mebu_out2(i) = 0;    ! mebu_out2 = ebullition from the layers are temporarily saved here
        Vtmp(i) = Vp(i)
        enddo

        ! layers with water peat mixture exist
        ! looping from bottom to top layer of water-peat mix
        do i = 10, mwt, -1
        ! CH4 molar density in the bubble (mol m-3)           ! Shuang mnv=cb in their paper
        mnV = f * presP(i)/(Rgas * tempP(i))
                !! 		 CH4 molar density in the water (mol m-3)
                !        mCwP = methaneP(i) / ((wsc(i)*0.001)*1);
                !	     methanebP(i) = f * presP(i) * Vp/(Rgas * tempP(i))  !unit mol/layer

        ! releasing the bubble if it exceeds certain size   ! threshold Vmax
        !        write (*,*) 'i',i,'Vtmp(i)',Vtmp(i),'Vmaxfraction', Vmaxfraction,'wsc(i)',wsc(i)
        if ((Vtmp(i) - Vmaxfraction * (wsc(i)*0.001)*1) > 1e-11) then

        ! bubble was released

        !Equation 8)
        Vout = Vtmp(i) - Vmaxfraction * ((wsc(i)*0.001)*1)
        Vtmp(i) = Vmaxfraction * ((wsc(i)*0.001)*1)

        ! the size of the bubble increases as it ascends, but the amount of moles in the bubble stay the same
        nout = mnV * Vout                             ! Shuang nout=Febu in the paper, unit mol, mnV=cb, Vout=the rest
        methanebP(i) = methanebP(i) - nout
        !          write (*,*) 'mnV * Vout  i',i,'mnV',mnV,'Vout',Vout,'nout',nout
        !***** this is 'if bubble got stuck loop'
        ind = i
        do while (Vout > 0 .AND. ind >= mwt + 1)
        ind = ind - 1

        call RANDOM_NUMBER(randbub)
        !            write (*,*) 'randbub',randbub
        if (randbub <= bubprob) then
            ! bubble got stuck, and bubble is added back to the upper i-1 layer
        methanebP(ind) = methanebP(ind) + nout
        Vtmp(ind) = methanebP(ind) * Rgas * tempP(ind)/(f * presP(ind))

        Vout = 0
        nout = 0
        end if

        end do

        ! bubble did not get stuck => it is released to the lowest air layer
        ! This 'if' loop and the 'do while' loop is either or, Vout will be 0 if the do while loop worked
        if (Vout > 0) then
        mebu_out2(i) = nout/1	!dt  here I assign 1 hour

        end if
        !     write (*,*) 'Vout > 0  i',i,'Vout',Vout,'mebu_out2(i)',mebu_out2(i)
        end if
        !        write (*,*) 'i',i,'Vout',Vout,'mnV',mnV,'nout',nout
        end do   !end of do i = 10, 2, -1
        !     write (*,*) "3 stopping point"


        ! updating the bubble volume profile
        Vp = Vtmp
        ! a1: First box of peat-air. If 0, there is no peat-air, a2 has no meaning.
        ! a2: Last box of peat-air.
        ! w1: First box of peat-water. If 0, no peat-water, w2 has no meaning.
        ! w2: Last box of peat-water.
        if (mwt .gt. 1) then! if air-peat layer is present
        met_db(mwt-1) = sum(mebu_out2) ! bubbles released to deepest air layer, sum of all the water layers
        ! ebu_out is already 0, no need to set again
        else			      ! if all layers are flooded, porewater CH4 was already updated
        mebu_out = sum(mebu_out2);
        !        write (*,*) 'mwt mebu_out',mebu_out,mebu_out2
        end if

        ! end of releasing part of the gas volume if it exceeds certain size
        ! *****************************************************************************


        ! *****************************************************************************
        do i = 1,nlayers
        ! bubbles are released to the lowest air layer if wtd is below surface
        !   if wtd is above surfacr all the met_db(i)=0
        methaneP(i) = methaneP(i) + 1 * met_db(i) !big_dt replaced with 1hr  %%%%%%%%%%%% need to edit: add to with layer
        !        write (*,*) '+db i',i,'methaneP(i)',methaneP(i)
        !     methaneP = methaneP + 1 * met_db

        !  integer, parameter :: mebuair = 1, mebutot = 2
        !  1 = methane ebullition into air
        !  2 = methane ebullition total (when WTD is below peat surface, ebullition is
        !      released in the lowest air layer inside peat)

        ! for output
        CH4(i)=12*methaneP(i)
        if (wsc(i) .eq. 0) then
            CH4_V(i) = 0.
        else
            CH4_V(i)=CH4(i)/(wsc(i)*0.001)
        endif
        ! write(*,*)"CH44: ", i, methaneP(i),  met_db(i)
        enddo
        ! write(*,*)"CH4 components:",methaneP(1), mrateP(1), methanebP(1)

        ! mebu_out and met_db are larger than 0. the amount of CH4 gets out
        ! amount of CH4 released to the atmosphere  !%%%% I changed gases_out_bigstep into gases_out=Ebu_sum_unsat  sat
        Ebu_sum_sat = Ebu_sum_sat + (1 * mebu_out)*12		!big_dt replaced with 1hr  %%%% mol/layer to gC/layer %%%%%%%%
        ! amount of CH4 released to the lowest air layer if WTD is below the surface
        Ebu_sum_unsat = Ebu_sum_unsat + 1 * sum(met_db)*12	!big_dt replaced with 1hr  %%%% mol/layer to gC/layer %%%%%%%%

        if (Ebu_sum_sat .ne. 0.) then
        !        write (*,*) 'Ebu_sum_sat', Ebu_sum_sat
        endif
        if (Ebu_sum_unsat .ne. 0.) then
        !        write (*,*) 'Ebu_sum_unsat', Ebu_sum_unsat
        endif

        simuCH4=simuCH4+Ebu_sum_sat
        !    write (*,*) 'mebu_out',mebu_out,'Ebu_sum_sat',Ebu_sum_sat
        ! 	gases_out=gases_out vector containing CH4 flux caused by ebullition
        ! *****************************************************************************
        !    gases_out_rates = gases_out / dt    !Shuang dt=24*3600

        ! amount of CH4 stored in bubbles
        !	bubble_methane_tot = bubble_methane_tot+f * presP(i)/(Rgas * tempP(i)) * Vp(i)
        bubble_methane_tot = sum(f * presP/(Rgas * tempP) * Vp)   !this is equation (10)
        !    methanebub2 = bubble_methane_tot


        !    write (69,609) CH4(1),CH4(2),CH4(3),CH4(4),CH4(5),CH4(6),CH4(7),CH4(8),CH4(9),CH4(10), &	! concentrations
        !    & CH4_V(1),CH4_V(2),CH4_V(3),CH4_V(4),CH4_V(5),CH4_V(6),CH4_V(7),CH4_V(8),CH4_V(9),CH4_V(10), &
        !    & Vp(1),Vp(2),Vp(3),Vp(4),Vp(5),Vp(6),Vp(7),Vp(8),Vp(9),Vp(10), &		! bubble volumes
        !    & mrateP(1),mrateP(2),mrateP(3),mrateP(4),mrateP(5), &	! flux of methane in and out from the bubble
        !    & mrateP(6),mrateP(7),mrateP(8),mrateP(9),mrateP(10), &	! methane fluxes
        !    & Ebu_sum_sat,Ebu_sum_unsat,simuCH4,bubble_methane_tot
        !609	format(43(f11.4,","),f11.4)
        !    write (*,*) "4 stopping point"
        return
    end subroutine ebullition_EBG
    !   *************************************************************************************

    real function esat(T)   ! returns saturation vapour pressure in Pa
        real T
        esat = 610.78*exp(17.27*T/(T+237.3))
        return
    end

end module mod_soil