module driver
    ! simulation based on the forcing data
    use mod_data
    use mod_vegetation
    use mod_soil
    use mod_transfer
    use mod_upAndSum
    implicit none
    ! some paramaters for cycle and saving results
    integer daysOfyear, hoursOfmonth, hoursOfYear, daysOfmonth(12) 
    integer iTotHourly, iTotDaily, iTotMonthly  ! used to cycle the record of different frequences of the results

    contains
    subroutine teco_simu()
        ! Jian: change the cycle according to the forcing data. year, doy, hour
        implicit none
        integer year0, first_year                               ! year0: record the current year to judge whether a new year
        real    Difference                                      ! GPP-Rauto-NPP. Jian: no sure whether for balance? 
        real    RaLeaf,RaStem,RaRoot                            ! for summary the automatic respiration of leaf, stem, root
        integer dlayer                                          ! to run cycle of each layer.
        real    Q_soil                                          ! total soil carbon
        real    RECOh                                           ! ecosytem respiration
        real    ETh, Th, Eh                                     ! record hourly ET, transp, evap. Jian: why not use original variable?
        real    INTh,ROh,DRAINh,LEh,SHh                         ! 
        real    VPDh, LWH
        real    esat1 

        ! Jian: start the cycle of the forcing data
        first_year  = forcing%year(1)
        do iforcing = 1, nforcing 
            ! write(*,*)"iforcing: ", iforcing, N_deposit
            if (iforcing .eq. 1) then
                year0 = first_year             ! Jian: record whether it is a new year.
                iTotHourly  = 1
                iTotDaily   = 1
                iTotMonthly = 1
                call init_year()                                ! initilization the variables of a year
            endif
            iyear = forcing%year(iforcing)                      ! force%year
            iday  = forcing%doy(iforcing)                    
            ihour = forcing%hour(iforcing)
            ! if it is a new year
            if (iyear > year0) then                             ! a new year loop    
                call init_year()                                ! initilization the variables of a year
            endif
            ! leap year
            call isLeap_update_daysOfyear()

            ! for update the results of monthly and yearly
            call update_summary_init_monthly()

            ! initialize the daily variables to run hourly simulaiton.
            if (ihour .eq. 0) then
                ! a new day simulation.
                if (do_snow) then 
                    if (iyear .eq. first_year .and. iday .eq. 1.) then
                        ta     = -12.85                         ! since changed the ta criteria (0. to 1.e-10)) in calculating melt
                        rain_d = 0.                             ! dbmemo
                    endif
                    call snow_d()                               ! Jian: update snow_dsim snow_d(rain_d,lat,days,ta,snow_dsim,fa,fsub,rho_snow,melt,dcount,decay_m)                            
                    snow_depth_e = snow_dsim
                endif
                StemSap = AMIN1(Stemmax,SapS*bmStem)            ! Stemmax and SapS were input from parameter file, what are they? Unit? Maximum stem biomass? -JJJJJJJJJJJJJJJJJJJJJJ 
                RootSap = AMIN1(Rootmax,SapR*bmRoot)
                NSCmax  = 0.05*(StemSap+RootSap+QC(1))          ! Jian: update the NSCmax each step? and fixed NSCmin  = 5.? 
                if(Ta.gt.5.0) GDD5 = GDD5+Ta
                call init_day()                                 ! Jian: initilize the daily data.
            endif

            ! forcing data --------------------------------------------------------------------------------
            Tair  = forcing%Tair(iforcing)                      ! Tair
            Tsoil = forcing%Tsoil(iforcing)                     ! SLT
            co2ca = forcing%CO2(iforcing)*1.0E-6                ! CO2 concentration,ppm-->1.0E-6
            if (co2ca .lt. 0) co2ca = 380.0*1.0E-6              ! Jian: if no CO2 (-9999), then use the default value 
            Tair  = Tair  + Ttreat                              ! Jian: whether it has the treatment
            Tsoil = Tsoil + Ttreat
            if (CO2treat .ne. 0.) co2ca = CO2treat*1.0E-6 
            ! ----------------------------------------------------------                
            RH     = forcing%RH(iforcing)
            Dair   = forcing%VPD(iforcing)                      ! air water vapour defficit? Unit Pa
            rain   = forcing%Rain(iforcing)                     ! rain fal per hour
            wind   = ABS(forcing%WS(iforcing))                  ! wind speed m s-1
            PAR    = forcing%PAR(iforcing)                      ! Unit ? umol/s/m-2
            radsol = forcing%PAR(iforcing)                      ! unit ? PAR actually  Jian: or incoming shortwave/longwave radiation?
            dpatm  = forcing%PBOT(iforcing)
            if (do_ndep) N_deposit = forcing%Ndep(iforcing)*3600
            ! Ajust some unreasonable values 
            RH     = AMAX1(0.01,AMIN1(99.99,RH))                ! relative humidity
            esat1  = 610.78*exp(17.27*Tair/(Tair + 237.3))      ! intermediate parameter
            eairP  = esat1*RH/100.                              ! Added for SPRUCE, due to lack of VPD data. Jian: ? SPRUCE has the data? !air water vapour pressure
            Dair   = esat1-eairP                                ! Jian: confused that SPRUCE has the VPD data, why calculate it again?
            radsol = AMAX1(radsol,0.01)

            ! intially added for soil thermal/ soil water
            if (do_snow) then
                snow_depth = snow_depth_e
            else
                snow_depth = snow_in(iforcing)                  ! read from input file
            endif
            if (snow_depth .lt. 0.0) snow_depth = 0.0   
            snow_depth = snow_depth*100.                        ! change from m to cm  
            
            ! Jian: G and Esoil?
            if (do_soilphy) then 
                GOTO 160
            endif
            if(radsol.gt.10.0) then
                G = -25.0
            else
                G = 20.5
            endif
            if (isnan(G)) stop
            Esoil = 0.05*radsol
            if(radsol.LE.10.0) Esoil = 0.5*G
160 continue        
            ! for daily mean conditions 
            ta     = ta + tair/24.0                             ! sum of a day, for calculating daily mean temperature, snow_d and soilwater
            rain_d = rain_d+rain                                
            ! calculating scaling factor of NSC
            if(NSC.le.NSCmin)fnsc=0.0
            if(NSC.ge.NSCmax)fnsc=1.0
            if((NSC.lt.NSCmax).and.(NSC.gt.NSCmin))then 
                fnsc=(NSC-NSCmin)/(NSCmax-NSCmin)
            endif
            ! update vcmx0 and eJmx0 according to C/N of leaves
            Vcmx0 = Vcmax0*SNvcmax*1.0e-6
            eJmx0 = 1.67*Vcmx0 ! Weng 02/21/2011 Medlyn et al. 2002     
            ! run canopy module
            call canopy()
            ! run soil water processes
            call soilwater()                      
            ET = evap+transp
            
            ! ! Jian: to update module
            call respiration()
            ! THE Third Part: update LAI
            call plantgrowth()

            ! THE Fourth PART: simulating C influx allocation in pools
            call TCS_CN()  
            ! if (do_matrix) call matrix_struct() 
            call methane()       !update single value of Rh_pools,Tsoil,zwt,wsc 
           
            ! update NSC
            Rauto      = Rmain + Rgrowth + Rnitrogen
            NSC        = NSC + GPP - Rauto - (NPP-add)-store
            Difference = GPP - Rauto - NPP
            if(NSC<0)then
                bmstem = bmstem + NSC/0.48
                NPP    = NPP    + NSC
                NSN    = NSN    - NSC/CN(2)
                NSC    = 0.
            endif
            GL_d    = GL_d + NPP*alpha_L
            GW_d    = GW_d + NPP*alpha_W
            GR_d    = GR_d + NPP*alpha_R
            LFALL_d = LFALL_d + L_fall
            ! update
            RaLeaf  = RgLeaf + RmLeaf
            RaStem  = RgStem + RmStem
            RaRoot  = RgRoot + RmRoot + Rnitrogen
            WFALL_d = WFALL_d+ OutC(2) !_wood
            RFALL_d = RFALL_d+ OutC(3) !_root
            N_LG_d  = N_LG_d + N_leaf
            N_WG_d  = N_WG_d + N_wood
            N_RG_d  = N_RG_d + N_root
            N_LF_d  = N_LF_d + N_LF
            N_WF_d  = N_WF_d + N_WF
            N_RF_d  = N_RF_d + N_RF

            N_up_d    = N_up_d  + N_uptake
            N_fix_d   = N_fix_d + N_fixation
            N_dep_d   = N_dep_d + N_deposit
            N_leach_d = N_leach_d + N_leach
            N_vol_d   = N_vol_d + N_vol

            N_up_yr    = N_up_yr+N_uptake
            N_fix_yr   = N_fix_yr+N_fixation
            N_dep_yr   = N_dep_yr+N_deposit
            N_leach_yr = N_leach_yr+N_leach
            N_vol_yr   = N_vol_yr+N_vol

            R_Ntr_yr   = R_Ntr_yr + Rnitrogen

            ! *** ..int 
            do dlayer=1,10
                ice_d_simu(dlayer) = ice_d_simu(dlayer)+ice(dlayer) 
            enddo       
            do dlayer=1,11
                soilt_d_simu(dlayer) = soilt_d_simu(dlayer)+testout(dlayer)  
                ! first = surface soil temperature 2:11=1:10 layer soil temperatures 
            enddo                    
            do dlayer=1,10
                CH4V_d(dlayer) = CH4V_d(dlayer) + CH4_V(dlayer) 
            enddo                  
            zwt_d=zwt_d+zwt    ! ..int I doubt it... mean for zwt?     check later  Shuang 
            !   *** 
            ! ==================== test variables
            topfws_yr = topfws_yr + topfws/8760.
            omega_yr  = omega_yr  + omega/8760.
            fwsoil_yr = fwsoil_yr + fwsoil/8760.
            omega_d   = omega_d   + omega/24.

            ! Rhetero=Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
            Rhetero = Rh_pools(1) + Rh_pools(2) + Rh_pools(3) &
                &   + Rh_pools(4) + Rh_pools(5)
            Rsoil   = Rhetero+RmRoot+RgRoot+Rnitrogen
            NEE     = Rauto+Rhetero - GPP
            Q_soil  = QC(6) + QC(7) + QC(8)
            bmleaf  = QC(1)/0.48
            bmstem  = QC(2)/0.48
            bmroot  = QC(3)/0.48
            bmplant = bmleaf+bmroot+bmstem
            LAI     = bmleaf*SLA
            NMIN_d  = NMIN_d+N_miner
            ! output hourly
            Recoh   = Rhetero+Rauto
            ETh     = ET !*1000.
            Th      = transp !*1000.
            Eh      = evap !*1000.
            INTh    = -9999
            VPDh    = Dair/1000.
            ROh     = runoff !*1000.
            DRAINh  = -9999
            LEh     = ETh*((2.501-0.00236*Tair)*1000.0)/3600.
            SHh     = -9999
            LWh     = -9999
            NEP     = -NEE

            call updateHourly()
            call summaryHourly(iTotHourly)
            !  Write(*,*)"TEST_ALL_:", wetlandCH4_h(iTotHourly-1,:)
            ! sums of a day
            diff_d   = diff_d   + difference
            gpp_d_old    = gpp_d_old    + GPP
            gpp_ra   = gpp_ra   + Rauto
            npp_d_old    = npp_d_old    + NPP
            NEP_d    = NEP_d    + NEP
            NEE_d    = NEE_d    + NEE
            RECO_d   = RECO_d   + Recoh
            rh_d_old     = rh_d_old     + Rhetero
            ra_d_old     = Reco_d   - rh_d_old
            mat_Rh_d = mat_Rh_d + mat_Rh                        ! Jian: add for matrix form
            RLEAV_d  = RLEAV_d  + RmLeaf + RgLeaf
            RWOOD_d  = RWOOD_d  + RmStem + RgStem
            RROOT_d  = RROOT_d  + RmRoot + RgRoot + Rnitrogen
            Rsoil_d  = rh_d_old     + RROOT_d
            NUP_d    = NUP_d    + N_uptake
            NVOL_d   = NVOL_d   + N_vol
            NLEACH_d = NLEACH_d + N_leach
            transp_d = transp_d + transp*(24./dtimes)
            evap_d   = evap_d   + evap*(24./dtimes)
            ET_d     = transp_d + evap_d
            LE_d     = LE_d     + LEh/24.
            Hcanop_d = Hcanop_d + Hcanop/(24./dtimes)
            runoff_d = runoff_d + runoff
            ! added for MEMCMC also for generation of daily methane emission                  
            simuCH4_d = simuCH4_d + simuCH4
            Pro_sum_d = Pro_sum_d + Pro_sum
            Oxi_sum_d = Oxi_sum_d + Oxi_sum
            Fdifu1_d  = Fdifu1_d  + Fdifu(1)
            Ebu_sum_d = Ebu_sum_d + Ebu_sum
            Pla_sum_d = Pla_sum_d + Pla_sum

            call updateDaily()
           
            call updateMonthly(hoursOfmonth)

            ! sum of the whole year
            diff_yr = diff_yr + difference
            gpp_yr  = gpp_yr  + gpp
            NPP_yr  = NPP_yr  + NPP
            Rh_yr   = Rh_yr   + Rhetero
            Ra_yr   = Ra_yr   + Rauto
            Rh4_yr  = Rh4_yr  + Rh_pools(1)
            Rh5_yr  = Rh5_yr  + Rh_pools(2)
            Rh6_yr  = Rh6_yr  + Rh_pools(3)
            Rh7_yr  = Rh7_yr  + Rh_pools(4)
            Rh8_yr  = Rh8_yr  + Rh_pools(5)
            Pool1   = Pool1   + QC(1)/8760.
            Pool2   = Pool2   + QC(2)/8760.
            Pool3   = Pool3   + QC(3)/8760.
            Pool4   = Pool4   + QC(4)/8760.
            Pool5   = Pool5   + QC(5)/8760.
            Pool6   = Pool6   + QC(6)/8760.
            Pool7   = Pool7   + QC(7)/8760.
            Pool8   = Pool8   + QC(8)/8760.
            out1_yr = out1_yr + OutC(1)
            out2_yr = out2_yr + OutC(2)
            out3_yr = out3_yr + OutC(3)
            out4_yr = out4_yr + OutC(4)
            out5_yr = out5_yr + OutC(5)
            out6_yr = out6_yr + OutC(6)
            out7_yr = out7_yr + OutC(7)
            out8_yr = out8_yr + OutC(8)
            NEE_yr  = NEE_yr  + NEE
            GL_yr   = GL_yr   + NPP*alpha_L
            GW_yr   = GW_yr   + NPP*alpha_W
            GR_yr   = GR_yr   + NPP*alpha_R

            ! call updateYearly(hoursOfYear)

            ! ! added for soil thermal      unknown function check later   Shuang
            ! if((yr+first_year-1).eq.obs_soilwater(1,k1) .and.    &
            !     &     days .eq. obs_soilwater(2,k1) .and.      &
            !     &     (i-1).eq. obs_soilwater(3,k1))then
            !     Simu_soilwater(1:10,k1)=wcl(1:10)
            !     Simu_soiltemp(1:11,k1)=testout
            !     Simu_watertable(1,k1)=zwt
            !     k1=k1+1
            ! endif
            ! enddo              ! end of dtimes ! Jian: mark for end of dtimes in original model.
            if(ihour .eq. 23) then
                if((GDD5.gt.gddonset) .and. phenoset.eq.0) then
                    pheno    = iday    ! pheno=days
                    phenoset = 1
                endif
            endif

            if (iforcing .eq. 1) then
                i_record = 1
                nday4out = 0
            endif
            if (ihour .eq. 23) then
                call summaryDaily(iTotDaily)
                nday4out = nday4out+1
                Simu_dailyflux14(1,i_record)  = gpp_d_old 
                Simu_dailyflux14(2,i_record)  = NEE_d                   
                Simu_dailyflux14(3,i_record)  = Reco_d	!Rh             
                Simu_dailyflux14(4,i_record)  = npp_d_old!*1.5
                Simu_dailyflux14(5,i_record)  = ra_d_old!*1.5
                ! Simu_dailyflux14(6,i_record)  = QC(1)!mat_x(1)!QC(1)
                ! Simu_dailyflux14(7,i_record)  = QC(2)!mat_x(2)!QC(2)
                ! Simu_dailyflux14(8,i_record)  = QC(3)!mat_x(3)!QC(3)
                ! Simu_dailyflux14(9,i_record)  = QC(4)!mat_x(4)!QC(4)
                ! Simu_dailyflux14(10,i_record) = QC(5)!mat_x(5)!QC(5)
                ! Simu_dailyflux14(11,i_record) = QC(6)!mat_x(6)!QC(6)
                ! Simu_dailyflux14(12,i_record) = QC(7)!mat_x(7)!QC(7)
                ! Simu_dailyflux14(13,i_record) = QC(8)!mat_x(8)!QC(8)!*0.48
                Simu_dailyflux14(6,i_record)  = mat_x(1,1)!QC(1)
                Simu_dailyflux14(7,i_record)  = mat_x(2,1)!QC(2)
                Simu_dailyflux14(8,i_record)  = mat_x(3,1)!QC(3)
                Simu_dailyflux14(9,i_record)  = mat_x(4,1)!QC(4)
                Simu_dailyflux14(10,i_record) = mat_x(5,1)!QC(5)
                Simu_dailyflux14(11,i_record) = mat_x(6,1)!QC(6)
                Simu_dailyflux14(12,i_record) = mat_x(7,1)!QC(7)
                Simu_dailyflux14(13,i_record) = mat_x(8,1)!QC(8)!*0.48
                Simu_dailyflux14(14,i_record) = rh_d_old

                Simu_dailyflux14_2023(1,i_record)  = gpp_d_old 
                Simu_dailyflux14_2023(2,i_record)  = NEE_d                   
                Simu_dailyflux14_2023(3,i_record)  = Reco_d	!Rh             
                Simu_dailyflux14_2023(4,i_record)  = npp_d_old!*1.5
                Simu_dailyflux14_2023(5,i_record)  = ra_d_old!*1.5
                Simu_dailyflux14_2023(6,i_record)  = QC(1)!mat_x(1)!QC(1)
                Simu_dailyflux14_2023(7,i_record)  = QC(2)!mat_x(2)!QC(2)
                Simu_dailyflux14_2023(8,i_record)  = QC(3)!mat_x(3)!QC(3)
                Simu_dailyflux14_2023(9,i_record)  = QC(4)!mat_x(4)!QC(4)
                Simu_dailyflux14_2023(10,i_record) = QC(5)!mat_x(5)!QC(5)
                Simu_dailyflux14_2023(11,i_record) = QC(6)!mat_x(6)!QC(6)
                Simu_dailyflux14_2023(12,i_record) = QC(7)!mat_x(7)!QC(7)
                Simu_dailyflux14_2023(13,i_record) = QC(8)!mat_x(8)!QC(8)!*0.48
                ! Simu_dailyflux14_2023(6,i_record)  = mat_x(1) ! QC(1)
                ! Simu_dailyflux14_2023(7,i_record)  = mat_x(2) ! QC(2)
                ! Simu_dailyflux14_2023(8,i_record)  = mat_x(3) ! QC(3)
                ! Simu_dailyflux14_2023(9,i_record)  = mat_x(4) ! QC(4)
                ! Simu_dailyflux14_2023(10,i_record) = mat_x(5) ! QC(5)
                ! Simu_dailyflux14_2023(11,i_record) = mat_x(6) ! QC(6)
                ! Simu_dailyflux14_2023(12,i_record) = mat_x(7) ! QC(7)
                ! Simu_dailyflux14_2023(13,i_record) = mat_x(8) ! QC(8)!*0.48
                Simu_dailyflux14_2023(14,i_record) = rh_d_old
                Simu_dailyflux14_2023(15,i_record) = ta       ! for environmental scalar
                Simu_dailyflux14_2023(16,i_record) = omega_d  !  
                Simu_dailyflux14_2023(17,i_record) = mat_B(1,1)
                Simu_dailyflux14_2023(18,i_record) = mat_B(2,1)       ! for environmental scalar
                Simu_dailyflux14_2023(19,i_record) = mat_B(3,1)  !  
                Simu_dailyflux14_2023(20,i_record) = mat_x(1,1) ! QC(1)
                Simu_dailyflux14_2023(21,i_record) = mat_x(2,1) ! QC(2)
                Simu_dailyflux14_2023(22,i_record) = mat_x(3,1) ! QC(3)
                Simu_dailyflux14_2023(23,i_record) = mat_x(4,1) ! QC(4)
                Simu_dailyflux14_2023(24,i_record) = mat_x(5,1) ! QC(5)
                Simu_dailyflux14_2023(25,i_record) = mat_x(6,1) ! QC(6)
                Simu_dailyflux14_2023(26,i_record) = mat_x(7,1) ! QC(7)
                Simu_dailyflux14_2023(27,i_record) = mat_x(8,1) ! QC(8)!*0.48
                Simu_dailyflux14_2023(28,i_record) = mat_Rh_d   
                ! record_yr(i_record) = iyear
                i_record = i_record+1
            end if
            ! if ((iyear > year0) .or. ((year0 .eq. first_year) .and. (iday .eq.1))) then      
            if (iforcing < nforcing)then
                if (forcing%year(iforcing+1)>iyear) then            
                    year0        = iyear                                     ! update the record of year (year0)
                    storage      = accumulation
                    stor_use     = Storage/720.
                    accumulation = 0.0
                    onset=0
                endif
            else
                year0        = iyear                                     ! update the record of year (year0)
                storage      = accumulation
                stor_use     = Storage/720.
                accumulation = 0.0
                onset=0
            endif
        enddo
    end subroutine teco_simu

    subroutine isLeap_update_daysOfyear()    
        implicit none
        if (MOD(iyear, 4) .eq. 0)then
            if (MOD(iyear, 100) .eq. 100)then
                if (MOD(iyear, 400) .eq. 400)then
                    daysOfyear = 366
                else
                    daysOfyear = 365
                endif
            else
                daysOfyear = 366
            endif
        else
            daysOfyear = 365
        endif
    end subroutine isLeap_update_daysOfyear

    subroutine update_summary_init_monthly()
        implicit none
        if (daysOfyear .eq. 365) then ! common year
            hoursOfYear = 365*24
            daysOfmonth = (/31,59,90,120,151,181,212,243,273,304,334,365/)
            if((iday .eq. 365) .and. (ihour .eq. 23))  call summaryMonthly(iTotMonthly)
        else
            hoursOfYear = 366*24
            daysOfmonth = (/31,60,91,121,152,182,213,244,274,305,335,366/)
            if((iday .eq. 366) .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        endif

        ! select case (iday)
        if (iday .eq. 1) then
            ! case (1)  ! first day of Jan.
                hoursOfmonth = (daysOfmonth(1)-0)*24
                if (ihour .eq. 0) call init_monthly()
        endif
        if (iday .eq. daysOfmonth(1)+1)then
            ! case (daysOfmonth(1)+1)
                hoursOfmonth = (daysOfmonth(2)-daysOfmonth(1))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        endif
        if (iday .eq. daysOfmonth(2)+1)then
            ! case (daysOfmonth(2)+1)
                hoursOfmonth = (daysOfmonth(3)-daysOfmonth(2))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        endif
        if (iday .eq. daysOfmonth(3)+1)then
            ! case (daysOfmonth(3)+1)
                hoursOfmonth = (daysOfmonth(4)-daysOfmonth(3))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        ENDIF
        if (iday .eq. daysOfmonth(4)+1)then
            ! case (daysOfmonth(4)+1)
                hoursOfmonth = (daysOfmonth(5)-daysOfmonth(4))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        endif
        if (iday .eq. daysOfmonth(5)+1)then
            ! case (daysOfmonth(5)+1)
                hoursOfmonth = (daysOfmonth(6)-daysOfmonth(5))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        endif
        if(iday .eq. daysOfmonth(6)+1)then
            ! case (daysOfmonth(6)+1)
                hoursOfmonth = (daysOfmonth(7)-daysOfmonth(6))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        endif
        if(iday .eq. daysOfmonth(7)+1)then
            ! case (daysOfmonth(7)+1)
                hoursOfmonth = (daysOfmonth(8)-daysOfmonth(7))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        endif
        if(iday .eq. daysOfmonth(8)+1)then
            ! case (daysOfmonth(8)+1)
                hoursOfmonth = (daysOfmonth(9)-daysOfmonth(8))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        endif
        if(iday .eq. daysOfmonth(9)+1)then
            ! case (daysOfmonth(9)+1)
                hoursOfmonth = (daysOfmonth(10)-daysOfmonth(9))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        endif
        if(iday .eq. daysOfmonth(10)+1)then
            ! case (daysOfmonth(10)+1)
                hoursOfmonth = (daysOfmonth(11)-daysOfmonth(10))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        endif
        if(iday .eq. daysOfmonth(11)+1)then
            ! case (daysOfmonth(11)+1)
                hoursOfmonth = (daysOfmonth(12)-daysOfmonth(11))*24
                if (ihour .eq. 0) then
                    call summaryMonthly(iTotMonthly)
                    call init_monthly()
                endif
        endif
    end subroutine update_summary_init_monthly
end module driver