! This module is used to summary the values of hourly, daily, monthly and yearly
module mod_upAndSum
    use mod_data
    implicit NONE
    real convert_g2kg, convert_h2s

    contains
    subroutine updateHourly()
        implicit none
        ! integer iTotHourly
        convert_g2kg = 0.001
        convert_h2s  = 1/3600.
        ! carbon fluxes (KgC m-2 s-1) Jian: TECO unit is gC m-2 h-1
        gpp_h             = gpp*convert_g2kg*convert_h2s
        npp_h             = npp*convert_g2kg*convert_h2s
        nppLeaf_h         = NPP_L*convert_g2kg*convert_h2s
        nppWood_h         = NPP_W*convert_g2kg*convert_h2s  
        nppStem_h         = NPP_W*convert_g2kg*convert_h2s                   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues. Jian: TECO has no above ground woody tissues, set to equit wood
        nppRoot_h         = NPP_R*convert_g2kg*convert_h2s
        nppOther_h        = 0*convert_g2kg*convert_h2s                       ! Jian: no other storage, NSC seems different from other NPP.
        ra_h              = Rauto*convert_g2kg*convert_h2s
        raLeaf_h          = Rmleaf*convert_g2kg*convert_h2s
        raStem_h          = Rmstem*convert_g2kg*convert_h2s
        raRoot_h          = Rmroot*convert_g2kg*convert_h2s
        raOther_h         = Rnitrogen *convert_g2kg*convert_h2s               ! Total C cost for nitrogen
        rMaint_h          = Rmain *convert_g2kg*convert_h2s                   ! maintenance respiration
        rGrowth_h         = Rgrowth *convert_g2kg*convert_h2s                 ! growth respiration
        rh_h              = Rhetero *convert_g2kg*convert_h2s                 ! heterotrophic respiration
        nbp_h             = (gpp - Rhetero - Rauto) *convert_g2kg*convert_h2s   ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        wetlandCH4_h      = simuCH4 *convert_g2kg*convert_h2s                 ! wetland net fluxes of CH4
        wetlandCH4prod_h  = Pro_sum *convert_g2kg*convert_h2s                 ! wetland net fluxes of CH4 production
        wetlandCH4cons_h  = Oxi_sum *convert_g2kg*convert_h2s                 ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        cLeaf_h           = QC(1)*convert_g2kg
        cStem_h           = QC(2)*convert_g2kg
        cRoot_h           = QC(3)*convert_g2kg
        cOther_h          = NSC*convert_g2kg                        ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        cLitter_h         = QC(4)*convert_g2kg                      ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        cLitterCwd_h      = QC(5)*convert_g2kg                      ! cLitterCwd: carbon in coarse woody debris
        cSoil_h           = (QC(6) + QC(7) + QC(8))*convert_g2kg    ! cSoil: soil organic carbon (Jian: total soil carbon);
        cSoilLevels_h     = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)                                       ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        cSoilFast_h       = QC(6)*convert_g2kg                      ! cSoilPools (different pools without depth)
        cSoilSlow_h       = QC(7)*convert_g2kg 
        cSoilPassive_h    = QC(8)*convert_g2kg 
        cCH4_h            = CH4*convert_g2kg                        ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        fBNF_h            = N_fixation*convert_g2kg*convert_h2s                 ! fBNF: biological nitrogen fixation;
        ! fN2O_h            = N_miner*convert_g2kg*convert_h2s                    ! fN2O: loss of nitrogen through emission of N2O;
        ! fNloss_h          = N_loss*convert_g2kg*convert_h2s                     ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        ! fNnetmin_h        = fNnetmin*convert_g2kg*convert_h2s                   ! net mineralizaiton
        ! fNdep_h           = N_deposit*convert_g2kg*convert_h2s                  ! deposition of N
        fN2O_h            = (N_transfer+N_uptake+N_fixation)*convert_g2kg*convert_h2s                    ! fN2O: loss of nitrogen through emission of N2O;
        fNloss_h          = (N_leaf+N_wood+N_root)*convert_g2kg*convert_h2s                     ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        fNnetmin_h        = ((N_transfer+N_uptake+N_fixation)-(N_leaf+N_wood+N_root))*convert_g2kg*convert_h2s                   ! net mineralizaiton
        fNdep_h           = N_wood*convert_g2kg*convert_h2s                  ! deposition of N
        ! Nitrogen pools (kgN m-2)
        nLeaf_h           = QN(1)*convert_g2kg
        nStem_h           = QN(2)*convert_g2kg
        nRoot_h           = QN(3)*convert_g2kg
        nOther_h          = NSN*convert_g2kg         ! other N pool
        nLitter_h         = QN(4)*convert_g2kg
        nLitterCwd_h      = QN(5)*convert_g2kg
        nSoil_h           = (QN(6)+QN(7)+QN(8))*convert_g2kg
        nMineral_h        = QNminer*convert_g2kg                                 ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        hfls_h            = Hsoil ! Sensible heat flux;
        hfss_h            = Esoil ! Latent heat flux;
        SWnet_h           = 0 ! Net shortwave radiation;
        LWnet_h           = 0 ! Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        ec_h              = 0!evap*convert_g2kg*convert_h2s        ! Canopy evaporation;
        tran_h            = transp*convert_g2kg*convert_h2s      ! Canopy transpiration;
        es_h              = evap*convert_g2kg*convert_h2s ! Soil evaporation
        hfsbl_h           = sublim*convert_g2kg*convert_h2s ! Snow sublimation
        mrro_h            = runoff*convert_g2kg*convert_h2s
        mrros_h           = forcing%Rain(iforcing)    
        mrrob_h           = 0 ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        mrso_h            = liq_water*1000                                   ! Kg m-2, soil moisture in each soil layer
        tsl_h             = tsoil_layer(1:10)+273.15                            ! K, soil temperature in each soil layer Jian: not sure the tsoil_layer is correct or not
        tsland_h          = Tair+273.15                                   ! K, surface temperature
        wtd_h             = zwt/1000                                       ! m, Water table depth
        snd_h             = snow_depth/100                                ! m, Total snow depth, Jian: change from m to cm in code, and now change from cm to m
        lai_h             = LAI                                           ! m2 m-2, Leaf area index
        ! not used in SPRUCE-MIP
        
    end subroutine updateHourly

    subroutine updateDaily()
        implicit none
        ! carbon fluxes
        gpp_d             = gpp_d            + gpp_h/24
        npp_d             = npp_d            + npp_h/24
        nppLeaf_d         = nppLeaf_d        + nppLeaf_h/24
        nppWood_d         = nppWood_d        + nppWood_h/24
        nppStem_d         = nppStem_d        + nppStem_h/24                      ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        nppRoot_d         = nppRoot_d        + nppRoot_h/24
        nppOther_d        = nppOther_d       + nppOther_h/24
        ra_d              = ra_d             + ra_h/24
        raLeaf_d          = raLeaf_d         + raLeaf_h/24
        raStem_d          = raStem_d         + raStem_h/24
        raRoot_d          = raRoot_d         + raRoot_h/24
        raOther_d         = raOther_d        + raOther_h/24
        rMaint_d          = rMaint_d         + rMaint_h/24                           ! maintenance respiration
        rGrowth_d         = Rgrowth_d        + rGrowth_h/24                         ! growth respiration
        rh_d              = rh_d             + rh_h/24                                   ! heterotrophic respiration
        nbp_d             = nbp_d            + nbp_h/24                                 ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        wetlandCH4_d      = wetlandCH4_d     + wetlandCH4_h/24                   ! wetland net fluxes of CH4
        wetlandCH4prod_d  = wetlandCH4prod_d + wetlandCH4prod_h/24           ! wetland net fluxes of CH4 production
        wetlandCH4cons_d  = wetlandCH4cons_d + wetlandCH4cons_h/24           ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        cLeaf_d           = cLeaf_d          + cLeaf_h/24
        cStem_d           = cStem_d          + cStem_h/24
        cRoot_d           = cRoot_d          + cRoot_h/24
        cOther_d          = cOther_d         + cOther_h/24                                    ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        cLitter_d         = cLitter_d        + cLitter_h/24                                   ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        cLitterCwd_d      = cLitterCwd_d     + cLitterCwd_h/24                                ! cLitterCwd: carbon in coarse woody debris
        cSoil_d           = cSoil_d          + cSoil_h/24                                     ! cSoil: soil organic carbon (Jian: total soil carbon);
        cSoilLevels_d     = cSoilLevels_d    + cSoilLevels_h/24                               ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        cSoilFast_d       = cSoilFast_d      + cSoilFast_h/24                                 ! cSoilPools (different pools without depth)
        cSoilSlow_d       = cSoilSlow_d      + cSoilSlow_h/24 
        cSoilPassive_d    = cSoilPassive_d   + cSoilPassive_h/24 
        cCH4_d            = cCH4_d           + cCH4_h/24                                      ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        fBNF_d            = fBNF_d           + fBNF_h/24                               ! fBNF: biological nitrogen fixation;
        fN2O_d            = fN2O_d           + fN2O_h/24                               ! fN2O: loss of nitrogen through emission of N2O;
        fNloss_d          = fNloss_d         + fNloss_h/24                           ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        fNnetmin_d        = fNnetmin_d       + fNnetmin_h/24                       ! net mineralizaiton
        fNdep_d           = fNdep_d          + fNdep_h/24                             ! deposition of N
        ! Nitrogen pools (kgN m-2)
        nLeaf_d           = nLeaf_d          + nLeaf_h/24
        nStem_d           = nStem_d          + nStem_h/24
        nRoot_d           = nRoot_d          + nRoot_h/24
        nOther_d          = nOther_d         + nOther_h/24
        nLitter_d         = nLitter_d        + nLitter_h/24
        nLitterCwd_d      = nLitterCwd_d     + nLitter_h/24
        nSoil_d           = nSoil_d          + nSoil_h/24
        nMineral_d        = nMineral_d       + nMineral_h/24                                  ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        hfls_d            = hfls_d           + hfls_h/24                               ! Sensible heat flux;
        hfss_d            = hfss_d           + hfss_h/24                               ! Latent heat flux;
        SWnet_d           = SWnet_d          + SWnet_h/24                             ! Net shortwave radiation;
        LWnet_d           = LWnet_d          + LWnet_h/24                             ! Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        ec_d              = ec_d             + ec_h/24
        tran_d            = tran_d           + tran_h/24
        es_d              = es_d             + es_h/24                                   ! Canopy evaporation; Canopy transpiration; Soil evaporation
        hfsbl_d           = hfsbl_d          + hfsbl_h/24                             ! Snow sublimation
        mrro_d            = mrro_d           + mrro_h/24
        mrros_d           = mrros_d          + mrros_h/24
        mrrob_d           = mrrob_d          + mrrob_h/24                              ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        mrso_d            = mrso_d           + mrro_h/24                  ! Kg m-2, soil moisture in each soil layer
        tsl_d             = tsl_d            + tsl_h/24                   ! K, soil temperature in each soil layer
        tsland_d          = tsland_d         + tsland_h/24                ! K, surface temperature
        wtd_d             = wtd_d            + wtd_h/24                   ! m, Water table depth
        snd_d             = snd_d            + snd_h/24                   ! m, Total snow depth
        lai_d             = lai_d            + lai_h/24                   ! m2 m-2, Leaf area index
        ! not used in SPRUCE-MIP    
    end subroutine updateDaily

    subroutine updateMonthly(hoursOfmonth)
        implicit none
        integer hoursOfmonth
        ! carbon fluxes
        gpp_m             = gpp_m            + gpp_h            /hoursOfmonth
        npp_m             = npp_m            + npp_h            /hoursOfmonth
        nppLeaf_m         = nppLeaf_m        + nppLeaf_h        /hoursOfmonth
        nppWood_m         = nppWood_m        + nppWood_h        /hoursOfmonth
        nppStem_m         = nppStem_m        + nppStem_h        /hoursOfmonth ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        nppRoot_m         = nppRoot_m        + nppRoot_h        /hoursOfmonth
        nppOther_m        = nppOther_m       + nppOther_h       /hoursOfmonth
        ra_m              = ra_m             + ra_h             /hoursOfmonth
        raLeaf_m          = raLeaf_m         + raLeaf_h         /hoursOfmonth
        raStem_m          = raStem_m         + raStem_h         /hoursOfmonth
        raRoot_m          = raRoot_m         + raRoot_h         /hoursOfmonth
        raOther_m         = raOther_m        + raOther_h        /hoursOfmonth
        rMaint_m          = rMaint_m         + rMaint_h         /hoursOfmonth ! maintenance respiration
        rGrowth_m         = rGrowth_m        + rGrowth_h        /hoursOfmonth ! growth respiration
        rh_m              = rh_m             + rh_h             /hoursOfmonth ! heterotrophic respiration
        nbp_m             = nbp_m            + nbp_h            /hoursOfmonth ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        wetlandCH4_m      = wetlandCH4_m     + wetlandCH4_h     /hoursOfmonth ! wetland net fluxes of CH4
        wetlandCH4prod_m  = wetlandCH4prod_m + wetlandCH4prod_h /hoursOfmonth ! wetland net fluxes of CH4 production
        wetlandCH4cons_m  = wetlandCH4cons_m + wetlandCH4cons_h /hoursOfmonth ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        cLeaf_m           = cLeaf_m          + cLeaf_h          /hoursOfmonth     
        cStem_m           = cStem_m          + cStem_h          /hoursOfmonth
        cRoot_m           = cRoot_m          + cRoot_h          /hoursOfmonth      
        cOther_m          = cOther_m         + cOther_h         /hoursOfmonth ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        cLitter_m         = cLitter_m        + cLitter_h        /hoursOfmonth ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        cLitterCwd_m      = cLitterCwd_m     + cLitterCwd_h     /hoursOfmonth ! cLitterCwd: carbon in coarse woody debris
        cSoil_m           = cSoil_m          + cSoil_h          /hoursOfmonth ! cSoil: soil organic carbon (Jian: total soil carbon);
        cSoilLevels_m     = cSoilLevels_m    + cSoilLevels_h    /hoursOfmonth  ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        cSoilFast_m       = cSoilFast_m      + cSoilFast_h      /hoursOfmonth ! cSoilPools (different pools without depth)
        cSoilSlow_m       = cSoilSlow_m      + cSoilSlow_h      /hoursOfmonth
        cSoilPassive_m    = cSoilPassive_m   + cSoilPassive_h   /hoursOfmonth
        cCH4_m            = cCH4_m           + cCH4_h           /hoursOfmonth ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        fBNF_m            = fBNF_m           + fBNF_h           /hoursOfmonth ! fBNF: biological nitrogen fixation;
        fN2O_m            = fN2O_m           + fN2O_h           /hoursOfmonth ! fN2O: loss of nitrogen through emission of N2O;
        fNloss_m          = fNloss_m         + fNloss_h         /hoursOfmonth ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        fNnetmin_m        = fNnetmin_m       + fNnetmin_h       /hoursOfmonth ! net mineralizaiton
        fNdep_m           = fNdep_m          + fNdep_h          /hoursOfmonth ! deposition of N
        ! Nitrogen pools (kgN m-2)
        nLeaf_m           = nLeaf_m          + nLeaf_h          /hoursOfmonth
        nStem_m           = nStem_m          + nStem_h          /hoursOfmonth
        nRoot_m           = nRoot_m          + nRoot_h          /hoursOfmonth
        nOther_m          = nOther_m         + nOther_h         /hoursOfmonth
        nLitter_m         = nLitter_m        + nLitter_h        /hoursOfmonth
        nLitterCwd_m      = nLitterCwd_m     + nLitterCwd_h     /hoursOfmonth
        nSoil_m           = nSoil_m          + nSoil_h          /hoursOfmonth
        nMineral_m        = nMineral_m       + nMineral_h       /hoursOfmonth ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        hfls_m            = hfls_m           + hfls_h           /hoursOfmonth ! Sensible heat flux;
        hfss_m            = hfss_m           + hfss_h           /hoursOfmonth ! Latent heat flux;
        SWnet_m           = SWnet_m          + SWnet_h          /hoursOfmonth ! Net shortwave radiation;
        LWnet_m           = LWnet_m          + LWnet_h          /hoursOfmonth !    Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        ec_m              = ec_m             + ec_h             /hoursOfmonth
        tran_m            = tran_m           + tran_h           /hoursOfmonth
        es_m              = es_m             + es_h             /hoursOfmonth ! Canopy evaporation; Canopy transpiration; Soil evaporation
        hfsbl_m           = hfsbl_m          + hfsbl_h          /hoursOfmonth ! Snow sublimation
        mrro_m            = mrro_m           + mrro_h           /hoursOfmonth
        mrros_m           = mrros_m          + mrros_h          /hoursOfmonth
        mrrob_m           = mrrob_m          + mrrob_h          /hoursOfmonth ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        mrso_m            = mrso_m           + mrro_h           /hoursOfmonth                  ! Kg m-2, soil moisture in each soil layer
        tsl_m             = tsl_m            + tsl_h            /hoursOfmonth                   ! K, soil temperature in each soil layer
        tsland_m          = tsland_m         + tsland_h         /hoursOfmonth                ! K, surface temperature
        wtd_m             = wtd_m            + wtd_h            /hoursOfmonth                   ! m, Water table depth
        snd_m             = snd_m            + snd_h            /hoursOfmonth                   ! m, Total snow depth
        lai_m             = lai_m            + lai_h            /hoursOfmonth                   ! m2 m-2, Leaf area index
        ! not used in SPRUCE-MIP

        
    end subroutine updateMonthly

    subroutine updateYearly(hoursOfYear)
        implicit none
        integer hoursOfYear
        ! carbon fluxes
        gpp_y             = gpp_y            + gpp_h            /hoursOfYear
        npp_y             = npp_y            + npp_h            /hoursOfYear
        nppLeaf_y         = nppLeaf_y        + nppLeaf_h        /hoursOfYear
        nppWood_y         = nppWood_y        + nppWood_h        /hoursOfYear
        nppStem_y         = nppStem_y        + nppStem_h        /hoursOfYear        ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        nppRoot_y         = nppRoot_y        + nppRoot_h        /hoursOfYear
        nppOther_y        = nppOther_y       + nppOther_h       /hoursOfYear
        ra_y              = ra_y             + ra_h             /hoursOfYear
        raLeaf_y          = raLeaf_y         + raLeaf_h         /hoursOfYear
        raStem_y          = raStem_y         + raStem_h         /hoursOfYear
        raRoot_y          = raRoot_y         + raRoot_h         /hoursOfYear
        raOther_y         = raOther_y        + raOther_h        /hoursOfYear
        rMaint_y          = rMaint_y         + rMaint_h         /hoursOfYear        ! maintenance respiration
        rGrowth_y         = rGrowth_y        + rGrowth_h        /hoursOfYear        ! growth respiration
        rh_y              = rh_y             + rh_h             /hoursOfYear        ! heterotrophic respiration
        nbp_y             = nbp_y            + nbp_h            /hoursOfYear        ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        wetlandCH4_y      = wetlandCH4_y     + wetlandCH4_h     /hoursOfYear        ! wetland net fluxes of CH4
        wetlandCH4prod_y  = wetlandCH4prod_y + wetlandCH4prod_h /hoursOfYear        ! wetland net fluxes of CH4 production
        wetlandCH4cons_y  = wetlandCH4cons_y + wetlandCH4cons_h /hoursOfYear        ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        cLeaf_y           = cLeaf_y          + cLeaf_h          /hoursOfYear
        cStem_y           = cStem_y          + cStem_h          /hoursOfYear
        cRoot_y           = cRoot_y          + cRoot_h          /hoursOfYear
        cOther_y          = cOther_y         + cOther_h         /hoursOfYear        ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        cLitter_y         = cLitter_y        + cLitter_h        /hoursOfYear        ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        cLitterCwd_y      = cLitterCwd_y     + cLitterCwd_h     /hoursOfYear        ! cLitterCwd: carbon in coarse woody debris
        cSoil_y           = cSoil_y          + cSoil_h          /hoursOfYear        ! cSoil: soil organic carbon (Jian: total soil carbon);
        cSoilLevels_y     = cSoilLevels_y    + cSoilLevels_h    /hoursOfYear        ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        cSoilFast_y       = cSoilFast_y      + cSoilFast_h      /hoursOfYear        ! cSoilPools (different pools without depth)
        cSoilSlow_y       = cSoilSlow_y      + cSoilSlow_h      /hoursOfYear
        cSoilPassive_y    = cSoilPassive_y   + cSoilPassive_h   /hoursOfYear
        cCH4_y            = cCH4_y           + cCH4_h           /hoursOfYear        ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        fBNF_y            = fBNF_y           + fBNF_h           /hoursOfYear        ! fBNF: biological nitrogen fixation;
        fN2O_y            = fN2O_y           + fN2O_h           /hoursOfYear        ! fN2O: loss of nitrogen through emission of N2O;
        fNloss_y          = fNloss_y         + fNloss_h         /hoursOfYear        ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        fNnetmin_y        = fNnetmin_y       + fNnetmin_h       /hoursOfYear        ! net mineralizaiton
        fNdep_y           = fNdep_y          + fNdep_h          /hoursOfYear        ! deposition of N
        ! Nitrogen pools (kgN m-2)
        nLeaf_y           = nLeaf_y          + nLeaf_h          /hoursOfYear
        nStem_y           = nStem_y          + nStem_h          /hoursOfYear
        nRoot_y           = nRoot_y          + nRoot_h          /hoursOfYear
        nOther_y          = nOther_y         + nOther_h         /hoursOfYear
        nLitter_y         = nLitter_y        + nLitter_h        /hoursOfYear
        nLitterCwd_y      = nLitterCwd_y     + nLitterCwd_h     /hoursOfYear
        nSoil_y           = nSoil_y          + nSoil_h          /hoursOfYear
        nMineral_y        = nMineral_y       + nMineral_h       /hoursOfYear        ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        hfls_y            = hfls_y           + hfls_h           /hoursOfYear        ! Sensible heat flux;
        hfss_y            = hfss_y           + hfss_h           /hoursOfYear        ! Latent heat flux;
        SWnet_y           = SWnet_y          + SWnet_h          /hoursOfYear        ! Net shortwave radiation;
        LWnet_y           = LWnet_y          + LWnet_h          /hoursOfYear        ! Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        ec_y              = ec_y             + ec_h             /hoursOfYear
        tran_y            = tran_y           + tran_h           /hoursOfYear
        es_y              = es_y             + es_h             /hoursOfYear        ! Canopy evaporation; Canopy transpiration; Soil evaporation
        hfsbl_y           = hfsbl_y          + hfsbl_h          /hoursOfYear        ! Snow sublimation
        mrro_y            = mrro_y           + mrro_h           /hoursOfYear
        mrros_y           = mrros_y          + mrros_h          /hoursOfYear
        mrrob_y           = mrrob_y          + mrrob_h          /hoursOfYear        ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        mrso_y            = mrso_y           + mrro_h           /hoursOfYear                  ! Kg m-2, soil moisture in each soil layer
        tsl_y             = tsl_y            + tsl_h            /hoursOfYear                   ! K, soil temperature in each soil layer
        tsland_y          = tsland_y         + tsland_h         /hoursOfYear                ! K, surface temperature
        wtd_y             = wtd_y            + wtd_h            /hoursOfYear                   ! m, Water table depth
        snd_y             = snd_y            + snd_h            /hoursOfYear                   ! m, Total snow depth
        lai_y             = lai_y            + lai_h            /hoursOfYear
        test_gpp_y        = test_gpp_y       + test_gpp         /hoursOfYear

        ! not used in SPRUCE-MIP
        rain_yr   = rain_yr   + rain
        transp_yr = transp_yr + transp
        evap_yr   = evap_yr   + evap
        runoff_yr = runoff_yr + runoff
    end subroutine updateYearly

    subroutine summaryHourly(iTotHourly)
        implicit NONE
        integer iTotHourly
        ! summary in the total results
        all_gpp_h(iTotHourly)            = gpp_h         
        all_npp_h(iTotHourly)            = npp_h
        all_nppLeaf_h(iTotHourly)        = nppLeaf_h
        all_nppWood_h(iTotHourly)        = nppWood_h
        all_nppStem_h(iTotHourly)        = nppStem_h
        all_nppRoot_h(iTotHourly)        = nppRoot_h
        all_nppOther_h(iTotHourly)       = nppOther_h                  ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        all_ra_h(iTotHourly)             = ra_h
        all_raLeaf_h(iTotHourly)         = raLeaf_h
        all_raStem_h(iTotHourly)         = raStem_h
        all_raRoot_h(iTotHourly)         = raRoot_h
        all_raOther_h(iTotHourly)        = raOther_h
        all_rMaint_h(iTotHourly)         = rMaint_h
        all_rGrowth_h(iTotHourly)        = rGrowth_h                                             ! maintenance respiration and growth respiration
        all_rh_h(iTotHourly)             = rh_h
        all_nbp_h(iTotHourly)            = nbp_h                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        all_wetlandCH4_h(iTotHourly)     = wetlandCH4_h 
        all_wetlandCH4prod_h(iTotHourly) = wetlandCH4prod_h 
        all_wetlandCH4cons_h(iTotHourly) = wetlandCH4cons_h                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        all_cLeaf_h(iTotHourly)          = cLeaf_h
        all_cStem_h(iTotHourly)          = cStem_h
        all_cRoot_h(iTotHourly)          = cRoot_h
        all_cOther_h(iTotHourly)         = cOther_h                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        all_cLitter_h(iTotHourly)        = cLitter_h
        all_cLitterCwd_h(iTotHourly)     = cLitterCwd_h                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        all_cSoil_h(iTotHourly)          = cSoil_h
        all_cSoilLevels_h(iTotHourly,:)  = cSoilLevels_h
        all_cSoilFast_h(iTotHourly)      = cSoilFast_h
        all_cSoilSlow_h(iTotHourly)      = cSoilSlow_h
        all_cSoilPassive_h(iTotHourly)   = cSoilPassive_h                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        all_cCH4_h(iTotHourly,:)         = cCH4_h                                                        ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        all_fBNF_h(iTotHourly)           = fBNF_h
        all_fN2O_h(iTotHourly)           = fN2O_h
        all_fNloss_h(iTotHourly)         = fNloss_h
        all_fNnetmin_h(iTotHourly)       = fNnetmin_h
        all_fNdep_h(iTotHourly)          = fNdep_h                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        all_nLeaf_h(iTotHourly)          = nLeaf_h
        all_nStem_h(iTotHourly)          = nStem_h
        all_nRoot_h(iTotHourly)          = nRoot_h
        all_nOther_h(iTotHourly)         = nOther_h
        all_nLitter_h(iTotHourly)        = nLitter_h
        all_nLitterCwd_h(iTotHourly)     = nLitterCwd_h
        all_nSoil_h(iTotHourly)          = nSoil_h
        all_nMineral_h(iTotHourly)       = nMineral_h                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        all_hfls_h(iTotHourly)           = hfls_h
        all_hfss_h(iTotHourly)           = hfss_h
        all_SWnet_h(iTotHourly)          = SWnet_h
        all_LWnet_h(iTotHourly)          = LWnet_h                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        all_ec_h(iTotHourly)             = ec_h
        all_tran_h(iTotHourly)           = tran_h
        all_es_h(iTotHourly)             = es_h                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        all_hfsbl_h(iTotHourly)          = hfsbl_h                                                         ! Snow sublimation
        all_mrro_h(iTotHourly)           = mrro_h
        all_mrros_h(iTotHourly)          = mrros_h
        all_mrrob_h(iTotHourly)          = mrrob_h                                        ! Total runoff; Surface runoff; Subsurface runoff
        ! Other
        all_mrso_h(iTotHourly,:)         = mrso_h                                                   ! Kg m-2, soil moisture in each soil layer
        all_tsl_h(iTotHourly,:)          = tsl_h                                             ! K, soil temperature in each soil layer
        all_tsland_h(iTotHourly)         = tsland_h                                                 ! K, surface temperature
        all_wtd_h(iTotHourly)            = wtd_h                                                 ! m, Water table depth
        all_snd_h(iTotHourly)            = snd_h                                                 ! m, Total snow depth
        all_lai_h(iTotHourly)            = lai_h                     ! m2 m-2, Leaf area index
        all_gdd5_h(iTotHourly)           = GDD5
        all_onset_h(iTotHourly)           = onset
        all_storage_h(iTotHourly)        = storage
        all_add_h(iTotHourly)            = add
        all_accumulation_h(iTotHourly)   = accumulation
        ! all_test_h(5*(iTotHourly-1)+1:5*(iTotHourly),:) = test_gpp
        all_test_h(iTotHourly,:) = test_gpp

        iTotHourly = iTotHourly+1
    end subroutine summaryHourly
       
    subroutine summaryDaily(iTotDaily)
        implicit none
        integer iTotDaily
        ! Summary 
        ! daily: 
        ! ---------------------------------------------------------------------
        ! carbon fluxes (Kg C m-2 s-1)
        all_gpp_d(iTotDaily)            = gpp_d
        all_npp_d(iTotDaily)            = npp_d
        all_nppLeaf_d(iTotDaily)        = nppLeaf_d
        all_nppWood_d(iTotDaily)        = nppWood_d
        all_nppStem_d(iTotDaily)        = nppStem_d
        all_nppRoot_d(iTotDaily)        = nppRoot_d
        all_nppOther_d(iTotDaily)       = nppOther_d   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        all_ra_d(iTotDaily)             = ra_d
        all_raLeaf_d(iTotDaily)         = raLeaf_d
        all_raStem_d(iTotDaily)         = raStem_d
        all_raRoot_d(iTotDaily)         = raRoot_d
        all_raOther_d(iTotDaily)        = raOther_d
        all_rMaint_d(iTotDaily)         = rMaint_d
        all_rGrowth_d(iTotDaily)        = Rgrowth_d                                            ! maintenance respiration and growth respiration
        all_rh_d(iTotDaily)             = rh_d
        all_nbp_d(iTotDaily)            = nbp_d                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        all_wetlandCH4_d(iTotDaily)     = wetlandCH4_d
        all_wetlandCH4prod_d(iTotDaily) = wetlandCH4prod_d
        all_wetlandCH4cons_d(iTotDaily) = wetlandCH4cons_d                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        all_cLeaf_d(iTotDaily)          = cLeaf_d
        all_cStem_d(iTotDaily)          = cStem_d
        all_cRoot_d(iTotDaily)          = cRoot_d
        all_cOther_d(iTotDaily)         = cOther_d                           ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        all_cLitter_d(iTotDaily)        = cLitter_d
        all_cLitterCwd_d(iTotDaily)     = cLitterCwd_d                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        all_cSoil_d(iTotDaily)          = cSoil_d
        all_cSoilLevels_d(iTotDaily,:)    = cSoilLevels_d
        all_cSoilFast_d(iTotDaily)      = cSoilFast_d
        all_cSoilSlow_d(iTotDaily)      = cSoilSlow_d
        all_cSoilPassive_d(iTotDaily)   = cSoilPassive_d                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        all_cCH4_d(iTotDaily,:)         = cCH4_d                                                         ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        all_fBNF_d(iTotDaily)           = fBNF_d
        all_fN2O_d(iTotDaily)           = fN2O_d
        all_fNloss_d(iTotDaily)         = fNloss_d
        all_fNnetmin_d(iTotDaily)       = fNnetmin_d
        all_fNdep_d(iTotDaily)          = fNdep_d                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        all_nLeaf_d(iTotDaily)          = nLeaf_d
        all_nStem_d(iTotDaily)          = nStem_d
        all_nRoot_d(iTotDaily)          = nRoot_d
        all_nOther_d(iTotDaily)         = nOther_d
        all_nLitter_d(iTotDaily)        = nLitter_d
        all_nLitterCwd_d(iTotDaily)     = nLitterCwd_d
        all_nSoil_d(iTotDaily)          = nSoil_d
        all_nMineral_d(iTotDaily)       = nMineral_d                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        all_hfls_d(iTotDaily)           = hfls_d
        all_hfss_d(iTotDaily)           = hfss_d
        all_SWnet_d(iTotDaily)          = SWnet_d
        all_LWnet_d(iTotDaily)          = LWnet_d                              ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        all_ec_d(iTotDaily)             = ec_d
        all_tran_d(iTotDaily)           = tran_d
        all_es_d(iTotDaily)             = es_d                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        all_hfsbl_d(iTotDaily)          = hfsbl_d                                                         ! Snow sublimation
        all_mrro_d(iTotDaily)           = mrro_d
        all_mrros_d(iTotDaily)          = mrros_d
        all_mrrob_d(iTotDaily)          = mrrob_d                                        ! Total runoff; Surface runoff; Subsurface runoff
        ! Other
        all_mrso_d(iTotDaily,:)         = mrso_d                                          ! Kg m-2, soil moisture in each soil layer
        all_tsl_d(iTotDaily,:)          = tsl_d                                          ! K, soil temperature in each soil layer
        all_tsland_d(iTotDaily)         = tsland_d                                                 ! K, surface temperature
        all_wtd_d(iTotDaily)            = wtd_d                                                 ! m, Water table depth
        all_snd_d(iTotDaily)            = snd_d                                                 ! m, Total snow depth
        all_lai_d(iTotDaily)            = lai_d                                                 ! m2 m-2, Leaf area index

        iTotDaily = iTotDaily+1
    end subroutine summaryDaily

    subroutine summaryMonthly(iTotMonthly)
        implicit none
        integer iTotMonthly
        ! monthly
        ! carbon fluxes (Kg C m-2 s-1)
        all_gpp_m(iTotMonthly)            = gpp_m
        all_npp_m(iTotMonthly)            = npp_m
        all_nppLeaf_m(iTotMonthly)        = nppLeaf_m
        all_nppWood_m(iTotMonthly)        = nppWood_m
        all_nppStem_m(iTotMonthly)        = nppStem_m
        all_nppRoot_m(iTotMonthly)        = nppRoot_m
        all_nppOther_m(iTotMonthly)       = nppOther_m   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        all_ra_m(iTotMonthly)             = ra_m
        all_raLeaf_m(iTotMonthly)         = raLeaf_m
        all_raStem_m(iTotMonthly)         = raStem_m
        all_raRoot_m(iTotMonthly)         = raRoot_m
        all_raOther_m(iTotMonthly)        = raOther_m
        all_rMaint_m(iTotMonthly)         = rMaint_m
        all_rGrowth_m(iTotMonthly)        = rGrowth_m                                             ! maintenance respiration and growth respiration
        all_rh_m(iTotMonthly)             = rh_m
        all_nbp_m(iTotMonthly)            = nbp_m                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        all_wetlandCH4_m(iTotMonthly)     = wetlandCH4_m
        all_wetlandCH4prod_m(iTotMonthly) = wetlandCH4prod_m
        all_wetlandCH4cons_m(iTotMonthly) = wetlandCH4cons_m                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        all_cLeaf_m(iTotMonthly)          = cLeaf_m
        all_cStem_m(iTotMonthly)          = cStem_m
        all_cRoot_m(iTotMonthly)          = cRoot_m
        all_cOther_m(iTotMonthly)         = cOther_m                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        all_cLitter_m(iTotMonthly)        = cLitter_m
        all_cLitterCwd_m(iTotMonthly)     = cLitterCwd_m                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        all_cSoil_m(iTotMonthly)          = cSoil_m
        all_cSoilLevels_m(iTotMonthly,:)    = cSoilLevels_m
        all_cSoilFast_m(iTotMonthly)      = cSoilFast_m
        all_cSoilSlow_m(iTotMonthly)      = cSoilSlow_m
        all_cSoilPassive_m(iTotMonthly)   = cSoilPassive_m                           ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        all_cCH4_m(iTotMonthly,:)         = cCH4_m                                                      ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        all_fBNF_m(iTotMonthly)           = fBNF_m
        all_fN2O_m(iTotMonthly)           = fN2O_m
        all_fNloss_m(iTotMonthly)         = fNloss_m
        all_fNnetmin_m(iTotMonthly)       = fNnetmin_m
        all_fNdep_m(iTotMonthly)          = fNdep_m                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        all_nLeaf_m(iTotMonthly)          = nLeaf_m
        all_nStem_m(iTotMonthly)          = nStem_m
        all_nRoot_m(iTotMonthly)          = nRoot_m
        all_nOther_m(iTotMonthly)         = nOther_m
        all_nLitter_m(iTotMonthly)        = nLitter_m
        all_nLitterCwd_m(iTotMonthly)     = nLitterCwd_m
        all_nSoil_m(iTotMonthly)          = nSoil_m
        all_nMineral_m(iTotMonthly)       = nMineral_m                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        all_hfls_m(iTotMonthly)           = hfls_m
        all_hfss_m(iTotMonthly)           = hfss_m
        all_SWnet_m(iTotMonthly)          = SWnet_m
        all_LWnet_m(iTotMonthly)          = LWnet_m                                ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        all_ec_m(iTotMonthly)             = ec_m
        all_tran_m(iTotMonthly)           = tran_m
        all_es_m(iTotMonthly)             =  es_m                                             ! Canopy evaporation; Canopy transpiration; Soil evaporation
        all_hfsbl_m(iTotMonthly)          = hfsbl_m                                                         ! Snow sublimation
        all_mrro_m(iTotMonthly)           = mrro_m
        all_mrros_m(iTotMonthly)          = mrros_m
        all_mrrob_m(iTotMonthly)          = mrrob_m 
        ! Other
        all_mrso_m(iTotMonthly,:)         = mrso_m                                         ! Kg m-2, soil moisture in each soil layer
        all_tsl_m(iTotMonthly,:)          = tsl_m                                         ! K, soil temperature in each soil layer
        all_tsland_m(iTotMonthly)         = tsland_m                                                 ! K, surface temperature
        all_wtd_m(iTotMonthly)            = wtd_m                                                 ! m, Water table depth
        all_snd_m(iTotMonthly)            = snd_m                                                 ! m, Total snow depth
        all_lai_m(iTotMonthly)            = lai_m                                                 ! m2 m-2, Leaf area index

        iTotMonthly = iTotMonthly + 1
    end subroutine summaryMonthly
end module mod_upAndSum