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
        integer aaatest

        aaatest = 0

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

            if ((iday .eq. 1) .and. (ihour .eq. 0)) call init_update_year()
            if (do_simu .and. (iday .eq. 1) .and. (ihour .eq. 0)) write(*,*)iyear
            if ((iyear .eq. 1974) .and. (iday .eq. 1) .and. (ihour .eq. 0))then
                ! 1974 remove 99% of tree biomass
                QC(1)    = 0.1 * QC(1)
                QC(2)    = 0.1 * QC(2)
                QC(3)    = 0.1 * QC(3)
                QN(1)    = 0.1 * QN(1)
                QN(2)    = 0.1 * QN(2)
                QN(3)    = 0.1 * QN(3)
                bmleaf   = 0.1 * bmleaf
                bmstem   = 0.1 * bmstem
                bmroot   = 0.1 * bmroot
                nsc      = 0.1 * nsc
                nsn      = 0.1 * nsn
                storage  = 0.1 * storage
                lai      = LAIMIN!0.1 * lai
                stor_use = 0.1 * stor_use
            endif

            ! leap year
            if (do_leap) then
                if(iday .eq. 1) call isLeap_update_daysOfyear()
            else
                daysOfyear = 365
            endif

            ! for update the results of monthly and yearly
            call update_hoursOfYear_daysOfmonth_initMonthly()
            ! call update_summary_init_monthly()

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
                ! if (GDD5 .ne. 0.) write(*,*)"GDD: ",iforcing, GDD5
                ! if (onset .eq. 1) write(*,*)"GDD: ",iforcing, GDD5
                ! write(*,*) aaatest
                if(Ta.gt.5.0) GDD5 = GDD5+Ta
                call init_day()                                 ! Jian: initilize the daily data.
            endif
            ! if ((itest .eq.1) .and. (iforcing .eq. 180)) stop

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
                NPP    = Amax1(NPP + NSC, 0.) 
                NSN    = NSN    - NSC/CN(2)
                NSC    = 0.
            endif

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
            call updateDaily()
            call updateMonthly(hoursOfmonth)

            call updateYearly(hoursOfYear)

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
            call update_summary_monthly()
                 
            if (iforcing < nforcing)then
                if (forcing%year(iforcing+1)>iyear) then            
                    year0 = iyear                                   ! update the record of year (year0)
                    storage      = accumulation
                    stor_use     = Storage/times_storage_use
                    accumulation = 0.0
                    onset        = 0
                endif
            else
                year0        = iyear                                     ! update the record of year (year0)
                storage      = accumulation
                stor_use     = Storage/times_storage_use
                accumulation = 0.0
                onset        = 0
            endif
        enddo
    end subroutine teco_simu

    subroutine isLeap_update_daysOfyear()    
        implicit none
        if (MOD(iyear, 4) .eq. 0)then
            if (MOD(iyear, 100) .eq. 0)then
                if (MOD(iyear, 400) .eq. 0)then
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

    subroutine update_hoursOfYear_daysOfmonth_initMonthly()
        implicit none
        if (daysOfyear .eq. 365) then ! common year
            hoursOfYear = 365*24
            daysOfmonth = (/31,59,90,120,151,181,212,243,273,304,334,365/)
        else
            hoursOfYear = 366*24
            daysOfmonth = (/31,60,91,121,152,182,213,244,274,305,335,366/)
        endif
        ! hours of month
        ! January:
        if (iday .eq. 1)then 
            hoursOfmonth = (daysOfmonth(1)-0)*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! Feburay:
        if (iday .eq. daysOfmonth(1)+1)then
            hoursOfmonth = (daysOfmonth(2)-daysOfmonth(1))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! March
        if (iday .eq. daysOfmonth(2)+1)then
            hoursOfmonth = (daysOfmonth(3)-daysOfmonth(2))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! April
        if (iday .eq. daysOfmonth(3)+1)then
            hoursOfmonth = (daysOfmonth(4)-daysOfmonth(3))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! May
        if (iday .eq. daysOfmonth(4)+1)then
            hoursOfmonth = (daysOfmonth(5)-daysOfmonth(4))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! June
        if (iday .eq. daysOfmonth(5)+1)then
            hoursOfmonth = (daysOfmonth(6)-daysOfmonth(5))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! July
        if(iday .eq. daysOfmonth(6)+1)then
            hoursOfmonth = (daysOfmonth(7)-daysOfmonth(6))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! Auguest
        if(iday .eq. daysOfmonth(7)+1)then
            hoursOfmonth = (daysOfmonth(8)-daysOfmonth(7))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! Septemble
        if(iday .eq. daysOfmonth(8)+1)then
            hoursOfmonth = (daysOfmonth(9)-daysOfmonth(8))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! October
        if(iday .eq. daysOfmonth(9)+1)then
            hoursOfmonth = (daysOfmonth(10)-daysOfmonth(9))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! November
        if(iday .eq. daysOfmonth(10)+1)then
            hoursOfmonth = (daysOfmonth(11)-daysOfmonth(10))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! December
        if(iday .eq. daysOfmonth(11)+1)then
            hoursOfmonth = (daysOfmonth(12)-daysOfmonth(11))*24
            if (ihour .eq. 0) call init_monthly()
        endif
    end subroutine update_hoursOfYear_daysOfmonth_initMonthly

    subroutine update_summary_monthly()
        implicit none
        ! January
        if ((iday .eq. daysOfmonth(1))  .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(2))  .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(3))  .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(4))  .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(5))  .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(6))  .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(7))  .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(8))  .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(9))  .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(10)) .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(11)) .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
        if ((iday .eq. daysOfmonth(12)) .and. (ihour .eq. 23)) call summaryMonthly(iTotMonthly)
    end subroutine update_summary_monthly
end module driver