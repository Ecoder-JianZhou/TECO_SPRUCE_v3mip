module mod_ncd_io
    use netcdf
    use mod_data
    implicit none
    CHARACTER(len=4) :: str_startyr, str_endyr
    
    contains
    subroutine spruce_mip_cmip6Format()
        use netcdf
        implicit none
        ! Daily and monthly
        ! carbon flux (KgC m-2 s-1): gpp, npp, nppLeaf, nppWood, nppRoot, nppOther,
        !              ra, raLeaf, raStem, raRoot, raOther, rMaint, rGrowth, rh
        !              nbp (=gpp - Rh - Ra - other losses)
        !              wetlandCH4, wetlandCH4prod, wetlandCH4cons
        ! carbon pools (KgC m-2): cLeaf, cStem, cRoot, cOther, cLitter (excluding coarse wood debris), cLitterCwd
        !              cSoil, cSoilLevels, cSoilPools (soil organic carbon for each pool), cCH4 (Methane concentration)
        ! Nitrogen flux (KgN m-2 s-1) : fBNF(biological nitrogen fixation), fN2O, fNloss, fNnetmin, fNdep
        ! Nitrogen pools (KgN m-2): nleaf, nStem, nRoot, nOther, nLitter, nLitterCwd, nSoil, nMineral
        ! Energy Fluxes (W m-2): hfls(sensible heat flux), hfss(Latent heat flux), SWnet (Net Shortwave radiation), LWnet(Net Longwave radiation)
        ! Water Fluxes  (Kg m-2 s-1): ec(canopy evaporation), tran(canopy transpiration), es(soil evaporation), hfsbl (snow sublimation), mrro(total runoff),
        !                 mrros (surface runoff), mrrob(subsurface runoff)
        ! other         : mrso (soil moisture in each soil layer, Kg m-2), tsl(soil temperature in each soil layer, K), tsland(surface temperature, K),
        !                 wtd (Water table depth, m), snd (total snow depth, m), lai(m2 m-2) 
        ! ===================================================================================================================================================
        ! carbon fluxes variables
        ! ----------:-----------:----------:-----------------------
        write(str_startyr,"(I4)")forcing%year(1)
        write(str_endyr,"(I4)")forcing%year(nforcing)
        

        ! hourly outputs
        ! :----------:-----------:----------:-----------:----------:-----------:----------:-----------
        ! hourly GPP
        call write_nc(outDir_h,nHours,all_gpp_h,"gpp","kgC m-2 s-1", "gross primary productivity","hourly",1)
        ! hourly NPP
        call write_nc(outDir_h,nHours,all_npp_h,"npp","kgC m-2 s-1", "Total net primary productivity","hourly",1)
        ! hourly leaf NPP
        call write_nc(outDir_h,nHours,all_nppLeaf_h,"nppLeaf","kgC m-2 s-1", "NPP allocated to leaf tissues","hourly",1)
        ! Hourly wood NPP
        call write_nc(outDir_h,nHours,all_nppWood_h,"nppWood","kgC m-2 s-1", &
            & "NPP allocated to above ground woody tissues","hourly",1)
        ! Hourly stem NPP
        call write_nc(outDir_h,nHours,all_nppStem_h,"nppStem","kgC m-2 s-1", "NPP allocated to stem tissues","hourly",1)
        ! Hourly root NPP
        call write_nc(outDir_h,nHours,all_nppRoot_h,"nppRoot","kgC m-2 s-1", "NPP allocated to root tissues","hourly",1)
        ! Hourly other NPP
        call write_nc(outDir_h,nHours,all_nppOther_h,"nppOther","kgC m-2 s-1", &
            & "NPP allocated to other plant organs (reserves, fruits, exudates)","hourly",1)
        ! Hourly ra 
        call write_nc(outDir_h,nHours,all_ra_h,"ra","kgC m-2 s-1", "Plant Autotrophic Respiration","hourly",1)
        ! Hourly leaf ra
        call write_nc(outDir_h,nHours,all_raLeaf_h,"raLeaf","kgC m-2 s-1", "Ra from leaves","hourly",1)
        ! Hourly stem ra
        call write_nc(outDir_h,nHours,all_raStem_h,"raStem","kgC m-2 s-1", "Ra from above ground woody tissues","hourly",1)
        ! Hourly raRoot_h
        call write_nc(outDir_h,nHours,all_raRoot_h,"raRoot","kgC m-2 s-1", "Ra from fine roots","hourly",1)
        ! Hourly raOther_h
        call write_nc(outDir_h,nHours,all_raOther_h,"raOther","kgC m-2 s-1", &
            & "Ra from other plant organs (reserves, fruits, exudates)","hourly",1)
        ! Hourly rMaint_h
        call write_nc(outDir_h,nHours,all_rMaint_h,"rMaint","kgC m-2 s-1", "Maintenance respiration","hourly",1)
        ! Hourly rGrowth_h                                             ! maintenance respiration and growth respiration
        call write_nc(outDir_h,nHours,all_rGrowth_h,"rGrowth","kgC m-2 s-1", "Growth respiration","hourly",1)
        ! Hourly rh_h
        call write_nc(outDir_h,nHours,all_rh_h,"rh","kgC m-2 s-1", "Heterotrophic respiration rate","hourly",1)
        ! Hourly nbp_h                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        call write_nc(outDir_h,nHours,all_nbp_h,"nbp","kgC m-2 s-1", &
            &"Net Biome productivity (NBP = GPP - Rh - Ra - other losses)","hourly",1)
        ! Hourly wetlandCH4_h
        call write_nc(outDir_h,nHours,all_wetlandCH4_h,"wetlandCH4","kgC m-2 s-1", "Net fluxes of CH4","hourly",1)
        ! Hourly wetlandCH4prod_h
        call write_nc(outDir_h,nHours,all_wetlandCH4prod_h,"wetlandCH4prod","kgC m-2 s-1", "CH4 production","hourly",1)
        ! Hourly wetlandCH4cons_h                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        call write_nc(outDir_h,nHours,all_wetlandCH4cons_h,"wetlandCH4cons","kgC m-2 s-1", "CH4 consumption","hourly",1)

        ! Carbon Pools  (KgC m-2)
        ! Hourly cLeaf_h
        call write_nc(outDir_h,nHours,all_cLeaf_h,"cLeaf","kgC m-2", "Carbon biomass in leaves","hourly",1)
        ! Hourly cStem_h
        call write_nc(outDir_h,nHours,all_cStem_h,"cStem","kgC m-2", "Carbon above ground woody biomass","hourly",1)
        ! Hourly cRoot_h
        call write_nc(outDir_h,nHours,all_cRoot_h,"cRoot","kgC m-2", "Carbon biomass in roots","hourly",1)
        ! Hourly cOther_h                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        call write_nc(outDir_h,nHours,all_cOther_h,"cOther","kgC m-2", &
            & "Carbon biomass in other plant organs (reserves, fruits)","hourly",1)
        ! Hourly cLitter_h
        call write_nc(outDir_h,nHours,all_cLitter_h,"cLitter","kgC m-2", &
            & "Carbon in litter (excluding coarse woody debris)","hourly",1)
        ! Hourly cLitterCwd_h                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        call write_nc(outDir_h,nHours,all_cLitterCwd_h,"cLitterCwd","kgC m-2", "Carbon in coarse woody debris","hourly",1)
        ! Hourly cSoil_h
        call write_nc(outDir_h,nHours,all_cSoil_h,"cSoil","kgC m-2", "Total soil organic carbon","hourly",1)

        ! Hourly cSoilLevels_h
        call write_nc(outDir_h,nHours,all_cSoilLevels_h,"cSoilLevels","kgC m-2", &
            & "Depth-specific soil organic carbon","hourly",nlayers)
        
        ! Hourly cSoilFast_h
        call write_nc(outDir_h,nHours,all_cSoilFast_h,"cSoilFast","kgC m-2", "Fast soil organic carbon","hourly",1)
        ! Hourly cSoilSlow_h
        call write_nc(outDir_h,nHours,all_cSoilSlow_h,"cSoilSlow","kgC m-2 s-1", "Slow soil organic carbon","hourly",1)
        ! Hourly cSoilPassive_h                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        call write_nc(outDir_h,nHours,all_cSoilPassive_h,"cSoilPassive","kgC m-2 s-1", "Passive soil organic carbon","hourly",1)
        ! Hourly cCH4_h                                                        ! methane concentration
        call write_nc(outDir_h,nHours,all_cCH4_h,"cCH4","kgC m-2 s-1", "Methane concentration","hourly",nlayers)
        
        ! Nitrogen fluxes (kgN m-2 s-1)
        ! Hourly fBNF_h
        call write_nc(outDir_h,nHours,all_fBNF_h,"fBNF","kgN m-2 s-1", "biological nitrogen fixation","hourly",1)
        ! Hourly fN2O_h
        call write_nc(outDir_h,nHours,all_fN2O_h,"fN2O","kgN m-2 s-1", "loss of nitrogen through emission of N2O","hourly",1)
        ! Hourly fNloss_h
        call write_nc(outDir_h,nHours,all_fNloss_h,"fNloss","kgN m-2 s-1", &
            & "Total loss of nitrogen to the atmosphere and from leaching","hourly",1)
        ! Hourly fNnetmin_h
        call write_nc(outDir_h,nHours,all_fNnetmin_h,"fNnetmin","kgN m-2 s-1", "net mineralization of N","hourly",1)
        ! Hourly fNdep_h                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        call write_nc(outDir_h,nHours,all_fNdep_h,"fNdep","kgN m-2 s-1", "Nitrogen deposition","hourly",1)
        
        ! Hourly  Nitrogen pools (kgN m-2)
        ! Hourly nLeaf_h
        call write_nc(outDir_h,nHours,all_nLeaf_h,"nLeaf","kgN m-2", "Nitrogen in leaves","hourly",1)
        ! Hourly nStem_h
        call write_nc(outDir_h,nHours,all_nStem_h,"nStem","kgN m-2", "Nitrogen in stems","hourly",1)
        ! Hourly nRoot_h
        call write_nc(outDir_h,nHours,all_nRoot_h,"nRoot","kgN m-2", "Nirogen in roots","hourly",1)
        ! Hourly nOther_h
        call write_nc(outDir_h,nHours,all_nOther_h,"nOther","kgN m-2", &
            & "nitrogen in other plant organs (reserves, fruits)","hourly",1)
        ! Hourly nLitter_h
        call write_nc(outDir_h,nHours,all_nLitter_h,"nLitter","kgN m-2", &
            & "Nitrogen in litter (excluding coarse woody debris)","hourly",1)
        ! Hourly nLitterCwd_h
        call write_nc(outDir_h,nHours,all_nLitterCwd_h,"nLitterCwd","kgN m-2", &
            & "Nitrogen in coarse woody debris","hourly",1)
        ! Hourly nSoil_h
        call write_nc(outDir_h,nHours,all_nSoil_h,"nSoil","kgN m-2", "Nitrogen in soil organic matter","hourly",1)
        ! Hourly nMineral_h                    ! nMineral: Mineral nitrogen pool
        call write_nc(outDir_h,nHours,all_nMineral_h,"nMineral","kgN m-2", "Mineral nitrogen pool","hourly",1)

        ! Hourly ! energy fluxes (W m-2)
        ! Hourly hfls_h
        call write_nc(outDir_h,nHours,all_hfls_h,"hfls","W m-2", "Sensible heat flux","hourly",1)
        ! Hourly hfss_h
        call write_nc(outDir_h,nHours,all_hfss_h,"hfss","W m-2", "Latent heat flux","hourly",1)
        ! Hourly SWnet_h
        call write_nc(outDir_h,nHours,all_SWnet_h,"SWnet","W m-2", "Net shortwave radiation","hourly",1)
        ! Hourly LWnet_h                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        call write_nc(outDir_h,nHours,all_LWnet_h,"LWnet","W m-2", "Net longwave radiation","hourly",1)

        ! Hourly ! water fluxes (kg m-2 s-1)
        ! Hourly ec_h
        call write_nc(outDir_h,nHours,all_ec_h,"ec","kg m-2 s-1", "Canopy evaporation","hourly",1)
        ! Hourly tran_h
        call write_nc(outDir_h,nHours,all_tran_h,"tran","kg m-2 s-1", "Canopy transpiration","hourly",1)
        ! Hourly es_h                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        call write_nc(outDir_h,nHours,all_es_h,"es","kg m-2 s-1", "Soil evaporation","hourly",1)
        ! Hourly hfsbl_h                                                         ! Snow sublimation
        call write_nc(outDir_h,nHours,all_hfsbl_h,"hfsbl","kg m-2 s-1", "Snow sublimation","hourly",1)
        ! Hourly mrro_h
        call write_nc(outDir_h,nHours,all_mrro_h,"mrro","kg m-2 s-1", "Total runoff","hourly",1)
        ! Hourly mrros_h
        call write_nc(outDir_h,nHours,all_mrros_h,"mrros","kg m-2 s-1", "Surface runoff","hourly",1)
        ! Hourly mrrob_h                                        ! Total runoff; Surface runoff; Subsurface runoff
        call write_nc(outDir_h,nHours,all_mrrob_h,"mrrob","kg m-2 s-1", "Subsurface runoff","hourly",1)
        ! Hourly Other
        ! Hourly mrso_h
        call write_nc(outDir_h,nHours,all_mrso_h,"mrso","kg m-2", "soil moisture in each soil layer","hourly",nlayers)      ! Kg m-2, soil moisture in each soil layer
        ! Hourly tsl_h 
        call write_nc(outDir_h,nHours,all_tsl_h,"tsl","K", "soil temperature in each soil layer","hourly",nlayers)
        ! Hourly tsland_h
        call write_nc(outDir_h,nHours,all_tsland_h,"tsland","K", "surface temperature","hourly",1)
        ! Hourly wtd_h
        call write_nc(outDir_h,nHours,all_wtd_h,"wtd","m", "Water table depth","hourly",1)
        ! Hourly snd_h
        call write_nc(outDir_h,nHours,all_snd_h,"snd","m", "Total snow depth","hourly",1)
        ! Hourly lai_h
        call write_nc(outDir_h,nHours,all_lai_h,"lai","m2 m-2", "Leaf area index","hourly",1)
        call write_nc(outDir_h,nHours,all_gdd5_h,"GDD5","m2 m-2", "GDD5","hourly",1)
        call write_nc(outDir_h,nHours,all_onset_h,"onset","m2 m-2", "onset","hourly",1)
        call write_nc(outDir_h,nHours,all_storage_h,"storage","m2 m-2", "onset","hourly",1)
        call write_nc(outDir_h,nHours,all_add_h,"add","m2 m-2", "onset","hourly",1)
        call write_nc(outDir_h,nHours,all_accumulation_h,"accumulation","m2 m-2", "accumulation","hourly",1)
        call write_nc(outDir_h,nHours,all_test_h,"test_gpp","m2 m-2", "test_gpp","hourly",9)



        ! daily: 
        ! ---------------------------------------------------------------------
        ! carbon fluxes (KgC m-2 s-1)
        ! Daily gpp_d
        call write_nc(outDir_d,nDays,all_gpp_d,"gpp","kgC m-2 s-1", "gross primary productivity","daily",1)
        ! daily NPP
        call write_nc(outDir_d,nDays,all_npp_d,"npp","kgC m-2 s-1", "Total net primary productivity","daily",1)
        ! daily leaf NPP
        call write_nc(outDir_d,nDays,all_nppLeaf_d,"nppLeaf","kgC m-2 s-1", "NPP allocated to leaf tissues","daily",1)
        ! daily wood NPP
        call write_nc(outDir_d,nDays,all_nppWood_d,"nppWood","kgC m-2 s-1", &
            & "NPP allocated to above ground woody tissues","daily",1)
        ! daily stem NPP
        call write_nc(outDir_d,nDays,all_nppStem_d,"nppStem","kgC m-2 s-1", "NPP allocated to stem tissues","daily",1)
        ! daily root NPP
        call write_nc(outDir_d,nDays,all_nppRoot_d,"nppRoot","kgC m-2 s-1", "NPP allocated to root tissues","daily",1)
        ! daily other NPP
        call write_nc(outDir_d,nDays,all_nppOther_d,"nppOther","kgC m-2 s-1", &
            & "NPP allocated to other plant organs (reserves, fruits, exudates)","daily",1)
        ! daily ra 
        call write_nc(outDir_d,nDays,all_ra_d,"ra","kgC m-2 s-1", "Plant Autotrophic Respiration","daily",1)
        ! daily leaf ra
        call write_nc(outDir_d,nDays,all_raLeaf_d,"raLeaf","kgC m-2 s-1", "Ra from leaves","daily",1)
        ! daily stem ra
        call write_nc(outDir_d,nDays,all_raStem_d,"raStem","kgC m-2 s-1", "Ra from above ground woody tissues","daily",1)
        ! daily raRoot_d
        call write_nc(outDir_d,nDays,all_raRoot_d,"raRoot","kgC m-2 s-1", "Ra from fine roots","daily",1)
        ! daily raOther_d
        call write_nc(outDir_d,nDays,all_raOther_d,"raOther","kgC m-2 s-1", &
            & "Ra from other plant organs (reserves, fruits, exudates)","daily",1)
        ! daily rMaint_d
        call write_nc(outDir_d,nDays,all_rMaint_d,"rMaint","kgC m-2 s-1", "Maintenance respiration","daily",1)
        ! daily rGrowth_d                                             ! maintenance respiration and growth respiration
        call write_nc(outDir_d,nDays,all_rGrowth_d,"rGrowth","kgC m-2 s-1", "Growth respiration","daily",1)
        ! daily rh_d
        call write_nc(outDir_d,nDays,all_rh_d,"rh","kgC m-2 s-1", "Heterotrophic respiration rate","daily",1)
        ! daily nbp_d                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        call write_nc(outDir_d,nDays,all_nbp_d,"nbp","kgC m-2 s-1", &
            & "Net Biome productivity (NBP = GPP - Rh - Ra - other losses)","daily",1)
        ! daily wetlandCH4_d
        call write_nc(outDir_d,nDays,all_wetlandCH4_d,"wetlandCH4","kgC m-2 s-1", "Net fluxes of CH4","daily",1)
        ! daily wetlandCH4prod_d
        call write_nc(outDir_d,nDays,all_wetlandCH4prod_d,"wetlandCH4prod","kgC m-2 s-1", "CH4 production","daily",1)
        ! daily wetlandCH4cons_d                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        call write_nc(outDir_d,nDays,all_wetlandCH4cons_d,"wetlandCH4cons","kgC m-2 s-1", "CH4 consumption","daily",1)

        ! Carbon Pools  (KgC m-2)
        ! daily cLeaf_d
        call write_nc(outDir_d,nDays,all_cLeaf_d,"cLeaf","kgC m-2", "Carbon biomass in leaves","daily",1)
        ! daily cStem_d
        call write_nc(outDir_d,nDays,all_cStem_d,"cStem","kgC m-2", "Carbon above ground woody biomass","daily",1)
        ! daily cRoot_d
        call write_nc(outDir_d,nDays,all_cRoot_d,"cRoot","kgC m-2", "Carbon biomass in roots","daily",1)
        ! daily cOther_d                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        call write_nc(outDir_d,nDays,all_cOther_d,"cOther","kgC m-2", &
            & "Carbon biomass in other plant organs (reserves, fruits)","daily",1)
        ! daily cLitter_d
        call write_nc(outDir_d,nDays,all_cLitter_d,"cLitter","kgC m-2", &
            & "Carbon in litter (excluding coarse woody debris)","daily",1)
        ! daily cLitterCwd_d                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        call write_nc(outDir_d,nDays,all_cLitterCwd_d,"cLitterCwd","kgC m-2", "Carbon in coarse woody debris","daily",1)
        ! daily cSoil_d
        call write_nc(outDir_d,nDays,all_cSoil_d,"cSoil","kgC m-2", "Total soil organic carbon","daily",1)

        ! daily cSoilLevels_d
        call write_nc(outDir_d,nDays,all_cSoilLevels_d,"cSoilLevels","kgC m-2", &
            & "Depth-specific soil organic carbon","daily",nlayers)
        
        ! daily cSoilFast_d
        call write_nc(outDir_d,nDays,all_cSoilFast_d,"cSoilFast","kgC m-2", "Fast soil organic carbon","daily",1)
        ! daily cSoilSlow_d
        call write_nc(outDir_d,nDays,all_cSoilSlow_d,"cSoilSlow","kgC m-2 s-1", "Slow soil organic carbon","daily",1)
        ! daily cSoilPassive_d                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        call write_nc(outDir_d,nDays,all_cSoilPassive_d,"cSoilPassive","kgC m-2 s-1", "Passive soil organic carbon","daily",1)
        ! daily cCH4_d                                                        ! methane concentration
        call write_nc(outDir_d,nDays,all_cCH4_d,"cCH4","kgC m-2 s-1", "Methane concentration","daily",nlayers)
        
        ! Nitrogen fluxes (kgN m-2 s-1)
        ! daily fBNF_d
        call write_nc(outDir_d,nDays,all_fBNF_d,"fBNF","kgN m-2 s-1", "biological nitrogen fixation","daily",1)
        ! daily fN2O_d
        call write_nc(outDir_d,nDays,all_fN2O_d,"fN2O","kgN m-2 s-1", "loss of nitrogen through emission of N2O","daily",1)
        ! daily fNloss_d
        call write_nc(outDir_d,nDays,all_fNloss_d,"fNloss","kgN m-2 s-1", &
            & "Total loss of nitrogen to the atmosphere and from leaching","daily",1)
        ! daily fNnetmin_d
        call write_nc(outDir_d,nDays,all_fNnetmin_d,"fNnetmin","kgN m-2 s-1", "net mineralization of N","daily",1)
        ! daily fNdep_d                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        call write_nc(outDir_d,nDays,all_fNdep_d,"fNdep","kgN m-2 s-1", "Nitrogen deposition","daily",1)
        
        ! daily  Nitrogen pools (kgN m-2)
        ! daily nLeaf_d
        call write_nc(outDir_d,nDays,all_nLeaf_d,"nLeaf","kgN m-2", "Nitrogen in leaves","daily",1)
        ! daily nStem_d
        call write_nc(outDir_d,nDays,all_nStem_d,"nStem","kgN m-2", "Nitrogen in stems","daily",1)
        ! daily nRoot_d
        call write_nc(outDir_d,nDays,all_nRoot_d,"nRoot","kgN m-2", "Nirogen in roots","daily",1)
        ! daily nOther_d
        call write_nc(outDir_d,nDays,all_nOther_d,"nOther","kgN m-2", &
            &"nitrogen in other plant organs (reserves, fruits)","daily",1)
        ! daily nLitter_d
        call write_nc(outDir_d,nDays,all_nLitter_d,"nLitter","kgN m-2",&
            & "Nitrogen in litter (excluding coarse woody debris)","daily",1)
        ! daily nLitterCwd_d
        call write_nc(outDir_d,nDays,all_nLitterCwd_d,"nLitterCwd","kgN m-2", "Nitrogen in coarse woody debris","daily",1)
        ! daily nSoil_d
        call write_nc(outDir_d,nDays,all_nSoil_d,"nSoil","kgN m-2", "Nitrogen in soil organic matter","daily",1)
        ! daily nMineral_d                    ! nMineral: Mineral nitrogen pool
        call write_nc(outDir_d,nDays,all_nMineral_d,"nMineral","kgN m-2", "Mineral nitrogen pool","daily",1)

        ! daily ! energy fluxes (W m-2)
        ! daily hfls_d
        call write_nc(outDir_d,nDays,all_hfls_d,"hfls","W m-2", "Sensible heat flux","daily",1)
        ! daily hfss_d
        call write_nc(outDir_d,nDays,all_hfss_d,"hfss","W m-2", "Latent heat flux","daily",1)
        ! daily SWnet_d
        call write_nc(outDir_d,nDays,all_SWnet_d,"SWnet","W m-2", "Net shortwave radiation","daily",1)
        ! daily LWnet_d                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        call write_nc(outDir_d,nDays,all_LWnet_d,"LWnet","W m-2", "Net longwave radiation","daily",1)

        ! daily ! water fluxes (kg m-2 s-1)
        ! daily ec_d
        call write_nc(outDir_d,nDays,all_ec_d,"ec","kg m-2 s-1", "Canopy evaporation","daily",1)
        ! daily tran_d
        call write_nc(outDir_d,nDays,all_tran_d,"tran","kg m-2 s-1", "Canopy transpiration","daily",1)
        ! daily es_d                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        call write_nc(outDir_d,nDays,all_es_d,"es","kg m-2 s-1", "Soil evaporation","daily",1)
        ! daily hfsbl_d                                                         ! Snow sublimation
        call write_nc(outDir_d,nDays,all_hfsbl_d,"hfsbl","kg m-2 s-1", "Snow sublimation","daily",1)
        ! daily mrro_d
        call write_nc(outDir_d,nDays,all_mrro_d,"mrro","kg m-2 s-1", "Total runoff","daily",1)
        ! daily mrros_d
        call write_nc(outDir_d,nDays,all_mrros_d,"mrros","kg m-2 s-1", "Surface runoff","daily",1)
        ! daily mrrob_d                                        ! Total runoff; Surface runoff; Subsurface runoff
        call write_nc(outDir_d,nDays,all_mrrob_d,"mrrob","kg m-2 s-1", "Subsurface runoff","daily",1)
        ! daily Other
        ! daily mrso_d
        call write_nc(outDir_d,nDays,all_mrso_d,"mrso","kg m-2", "soil moisture in each soil layer","daily",nlayers)      ! Kg m-2, soil moisture in each soil layer
        ! daily tsl_d 
        call write_nc(outDir_d,nDays,all_tsl_d,"tsl","K", "soil temperature in each soil layer","daily",nlayers)
        ! daily tsland_d
        call write_nc(outDir_d,nDays,all_tsland_d,"tsland","K", "surface temperature","daily",1)
        ! daily wtd_d
        call write_nc(outDir_d,nDays,all_wtd_d,"wtd","m", "Water table depth","daily",1)
        ! daily snd_d
        call write_nc(outDir_d,nDays,all_snd_d,"snd","m", "Total snow depth","daily",1)
        ! daily lai_d
        call write_nc(outDir_d,nDays,all_lai_d,"lai","m2 m-2", "Leaf area index","daily",1)

        ! monthly
        ! ---------------------------------------------------------------------
        ! carbon fluxes (KgC m-2 s-1)
        ! monthly gpp_m
        call write_nc(outDir_m,nMonths,all_gpp_m,"gpp","kgC m-2 s-1", "gross primary productivity","monthly",1)
        ! monthly NPP
        call write_nc(outDir_m,nMonths,all_npp_m,"npp","kgC m-2 s-1", "Total net primary productivity","monthly",1)
        ! monthly leaf NPP
        call write_nc(outDir_m,nMonths,all_nppLeaf_m,"nppLeaf","kgC m-2 s-1", "NPP allocated to leaf tissues","monthly",1)
        ! monthly wood NPP
        call write_nc(outDir_m,nMonths,all_nppWood_m,"nppWood","kgC m-2 s-1", &
            & "NPP allocated to above ground woody tissues","monthly",1)
        ! monthly stem NPP
        call write_nc(outDir_m,nMonths,all_nppStem_m,"nppStem","kgC m-2 s-1", "NPP allocated to stem tissues","monthly",1)
        ! monthly root NPP
        call write_nc(outDir_m,nMonths,all_nppRoot_m,"nppRoot","kgC m-2 s-1", "NPP allocated to root tissues","monthly",1)
        ! monthly other NPP
        call write_nc(outDir_m,nMonths,all_nppOther_m,"nppOther","kgC m-2 s-1", &
            & "NPP allocated to other plant organs (reserves, fruits, exudates)","monthly",1)
        ! monthly ra 
        call write_nc(outDir_m,nMonths,all_ra_m,"ra","kgC m-2 s-1", "Plant Autotrophic Respiration","monthly",1)
        ! monthly leaf ra
        call write_nc(outDir_m,nMonths,all_raLeaf_m,"raLeaf","kgC m-2 s-1", "Ra from leaves","monthly",1)
        ! monthly stem ra
        call write_nc(outDir_m,nMonths,all_raStem_m,"raStem","kgC m-2 s-1", "Ra from above ground woody tissues","monthly",1)
        ! monthly raRoot_m
        call write_nc(outDir_m,nMonths,all_raRoot_m,"raRoot","kgC m-2 s-1", "Ra from fine roots","monthly",1)
        ! monthly raOther_m
        call write_nc(outDir_m,nMonths,all_raOther_m,"raOther","kgC m-2 s-1", &
            & "Ra from other plant organs (reserves, fruits, exudates)","monthly",1)
        ! monthly rMaint_m
        call write_nc(outDir_m,nMonths,all_rMaint_m,"rMaint","kgC m-2 s-1", "Maintenance respiration","monthly",1)
        ! monthly rGrowth_m                                             ! maintenance respiration and growth respiration
        call write_nc(outDir_m,nMonths,all_rGrowth_m,"rGrowth","kgC m-2 s-1", "Growth respiration","monthly",1)
        ! monthly rh_m
        call write_nc(outDir_m,nMonths,all_rh_m,"rh","kgC m-2 s-1", "Heterotrophic respiration rate","monthly",1)
        ! monthly nbp_m                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        call write_nc(outDir_m,nMonths,all_nbp_m,"nbp","kgC m-2 s-1",&
            & "Net Biome productivity (NBP = GPP - Rh - Ra - other losses)","monthly",1)
        ! monthly wetlandCH4_m
        call write_nc(outDir_m,nMonths,all_wetlandCH4_m,"wetlandCH4","kgC m-2 s-1", "Net fluxes of CH4","monthly",1)
        ! monthly wetlandCH4prod_m
        call write_nc(outDir_m,nMonths,all_wetlandCH4prod_m,"wetlandCH4prod","kgC m-2 s-1", "CH4 production","monthly",1)
        ! monthly wetlandCH4cons_m                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        call write_nc(outDir_m,nMonths,all_wetlandCH4cons_m,"wetlandCH4cons","kgC m-2 s-1", "CH4 consumption","monthly",1)

        ! Carbon Pools  (KgC m-2)
        ! monthly cLeaf_m
        call write_nc(outDir_m,nMonths,all_cLeaf_m,"cLeaf","kgC m-2", "Carbon biomass in leaves","monthly",1)
        ! monthly cStem_m
        call write_nc(outDir_m,nMonths,all_cStem_m,"cStem","kgC m-2", "Carbon above ground woody biomass","monthly",1)
        ! monthly cRoot_m
        call write_nc(outDir_m,nMonths,all_cRoot_m,"cRoot","kgC m-2", "Carbon biomass in roots","monthly",1)
        ! monthly cOther_m                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        call write_nc(outDir_m,nMonths,all_cOther_m,"cOther","kgC m-2", &
            & "Carbon biomass in other plant organs (reserves, fruits)","monthly",1)
        ! monthly cLitter_m
        call write_nc(outDir_m,nMonths,all_cLitter_m,"cLitter","kgC m-2", &
            & "Carbon in litter (excluding coarse woody debris)","monthly",1)
        ! monthly cLitterCwd_m                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        call write_nc(outDir_m,nMonths,all_cLitterCwd_m,"cLitterCwd","kgC m-2", "Carbon in coarse woody debris","monthly",1)
        ! monthly cSoil_m
        call write_nc(outDir_m,nMonths,all_cSoil_m,"cSoil","kgC m-2", "Total soil organic carbon","monthly",1)

        ! monthly cSoilLevels_m
        call write_nc(outDir_m,nMonths,all_cSoilLevels_m,"cSoilLevels","kgC m-2", &
            & "Depth-specific soil organic carbon","monthly",nlayers)
        
        ! monthly cSoilFast_m
        call write_nc(outDir_m,nMonths,all_cSoilFast_m,"cSoilFast","kgC m-2", "Fast soil organic carbon","monthly",1)
        ! monthly cSoilSlow_m
        call write_nc(outDir_m,nMonths,all_cSoilSlow_m,"cSoilSlow","kgC m-2", "Slow soil organic carbon","monthly",1)
        ! monthly cSoilPassive_m                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        call write_nc(outDir_m,nMonths,all_cSoilPassive_m,"cSoilPassive","kgC m-2", "Passive soil organic carbon","monthly",1)
        ! monthly cCH4_m                                                        ! methane concentration
        call write_nc(outDir_m,nMonths,all_cCH4_m,"cCH4","kgC m-2", "Methane concentration","monthly",nlayers)
        
        ! Nitrogen fluxes (kgN m-2 s-1)
        ! monthly fBNF_m
        call write_nc(outDir_m,nMonths,all_fBNF_m,"fBNF","kgN m-2 s-1", "biological nitrogen fixation","monthly",1)
        ! monthly fN2O_m
        call write_nc(outDir_m,nMonths,all_fN2O_m,"fN2O","kgN m-2 s-1", "loss of nitrogen through emission of N2O","monthly",1)
        ! monthly fNloss_m
        call write_nc(outDir_m,nMonths,all_fNloss_m,"fNloss","kgN m-2 s-1", &
            & "Total loss of nitrogen to the atmosphere and from leaching","monthly",1)
        ! monthly fNnetmin_m
        call write_nc(outDir_m,nMonths,all_fNnetmin_m,"fNnetmin","kgN m-2 s-1", "net mineralization of N","monthly",1)
        ! monthly fNdep_m                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        call write_nc(outDir_m,nMonths,all_fNdep_m,"fNdep","kgN m-2 s-1", "Nitrogen deposition","monthly",1)
        
        ! monthly  Nitrogen pools (kgN m-2)
        ! monthly nLeaf_m
        call write_nc(outDir_m,nMonths,all_nLeaf_m,"nLeaf","kgN m-2", "Nitrogen in leaves","monthly",1)
        ! monthly nStem_m
        call write_nc(outDir_m,nMonths,all_nStem_m,"nStem","kgN m-2", "Nitrogen in stems","monthly",1)
        ! monthly nRoot_m
        call write_nc(outDir_m,nMonths,all_nRoot_m,"nRoot","kgN m-2", "Nirogen in roots","monthly",1)
        ! monthly nOther_m
        call write_nc(outDir_m,nMonths,all_nOther_m,"nOther","kgN m-2", &
            &"nitrogen in other plant organs (reserves, fruits)","monthly",1)
        ! monthly nLitter_m
        call write_nc(outDir_m,nMonths,all_nLitter_m,"nLitter","kgN m-2",&
            & "Nitrogen in litter (excluding coarse woody debris)","monthly",1)
        ! monthly nLitterCwd_m
        call write_nc(outDir_m,nMonths,all_nLitterCwd_m,"nLitterCwd","kgN m-2", "Nitrogen in coarse woody debris","monthly",1)
        ! monthly nSoil_m
        call write_nc(outDir_m,nMonths,all_nSoil_m,"nSoil","kgN m-2", "Nitrogen in soil organic matter","monthly",1)
        ! monthly nMineral_m                    ! nMineral: Mineral nitrogen pool
        call write_nc(outDir_m,nMonths,all_nMineral_m,"nMineral","kgN m-2", "Mineral nitrogen pool","monthly",1)

        ! monthly ! energy fluxes (W m-2)
        ! monthly hfls_m
        call write_nc(outDir_m,nMonths,all_hfls_m,"hfls","W m-2", "Sensible heat flux","monthly",1)
        ! monthly hfss_m
        call write_nc(outDir_m,nMonths,all_hfss_m,"hfss","W m-2", "Latent heat flux","monthly",1)
        ! monthly SWnet_m
        call write_nc(outDir_m,nMonths,all_SWnet_m,"SWnet","W m-2", "Net shortwave radiation","monthly",1)
        ! monthly LWnet_m                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        call write_nc(outDir_m,nMonths,all_LWnet_m,"LWnet","W m-2", "Net longwave radiation","monthly",1)

        ! monthly ! water fluxes (kg m-2 s-1)
        ! monthly ec_m
        call write_nc(outDir_m,nMonths,all_ec_m,"ec","kg m-2 s-1", "Canopy evaporation","monthly",1)
        ! monthly tran_m
        call write_nc(outDir_m,nMonths,all_tran_m,"tran","kg m-2 s-1", "Canopy transpiration","monthly",1)
        ! monthly es_m                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        call write_nc(outDir_m,nMonths,all_es_m,"es","kg m-2 s-1", "Soil evaporation","monthly",1)
        ! monthly hfsbl_m                                                         ! Snow sublimation
        call write_nc(outDir_m,nMonths,all_hfsbl_m,"hfsbl","kg m-2 s-1", "Snow sublimation","monthly",1)
        ! monthly mrro_m
        call write_nc(outDir_m,nMonths,all_mrro_m,"mrro","kg m-2 s-1", "Total runoff","monthly",1)
        ! monthly mrros_m
        call write_nc(outDir_m,nMonths,all_mrros_m,"mrros","kg m-2 s-1", "Surface runoff","monthly",1)
        ! monthly mrrob_m                                        ! Total runoff; Surface runoff; Subsurface runoff
        call write_nc(outDir_m,nMonths,all_mrrob_m,"mrrob","kg m-2 s-1", "Subsurface runoff","monthly",1)
        ! monthly Other
        ! monthly mrso_m
        call write_nc(outDir_m,nMonths,all_mrso_m,"mrso","kg m-2", "soil moisture in each soil layer","monthly",nlayers)      ! Kg m-2, soil moisture in each soil layer
        ! monthly tsl_m 
        call write_nc(outDir_m,nMonths,all_tsl_m,"tsl","K", "soil temperature in each soil layer","monthly",nlayers)
        ! monthly tsland_m
        call write_nc(outDir_m,nMonths,all_tsland_m,"tsland","K", "surface temperature","monthly",1)
        ! monthly wtd_m
        call write_nc(outDir_m,nMonths,all_wtd_m,"wtd","m", "Water table depth","monthly",1)
        ! monthly snd_m
        call write_nc(outDir_m,nMonths,all_snd_m,"snd","m", "Total snow depth","monthly",1)
        ! monthly lai_m
        call write_nc(outDir_m,nMonths,all_lai_m,"lai","m2 m-2", "Leaf area index","monthly",1)
        
    end subroutine spruce_mip_cmip6Format
    ! ----------------------------------------------

    subroutine write_restart()
        implicit none
        integer(KIND=4) :: ncid, dimid_1, dimid_n, dimid_8
        integer(kind=4) :: id_eJmx0, id_QC, id_Storage
        integer(kind=4) :: id_nsc, id_accumulation, id_SNvcmax, id_LAI
        integer(kind=4) :: id_alphaN, id_NSN, id_QNminer, id_N_deficit
        integer(kind=4) :: id_liq_water, id_fwsoil, id_topfws, id_omega
        integer(kind=4) :: id_wcl, id_zwt, id_infilt, id_sftmp, id_Tsnow, id_Twater, id_Tice
        integer(kind=4) :: id_G, id_snow_dsim, id_dcount, id_dcount_soil, id_ice_tw, id_Tsoill, id_ice
        integer(kind=4) :: id_shcap_snow
        integer(kind=4) :: id_condu_snow
        integer(kind=4) :: id_condu_b
        integer(kind=4) :: id_depth_ex
        integer(kind=4) :: id_diff_s
        integer(kind=4) :: id_diff_snow
        integer(kind=4) :: id_albedo_snow
        integer(kind=4) :: id_thd_snow_depth
        integer(kind=4) :: id_b_bound
        integer(kind=4) :: id_infilt_rate
        integer(kind=4) :: id_fa
        integer(kind=4) :: id_fsub
        integer(kind=4) :: id_rho_snow
        integer(kind=4) :: id_CH4_V
        integer(kind=4) :: id_CH4
        integer(kind=4) :: id_Vp
        integer(kind=4) :: id_bubble_methane_tot
        integer(kind=4) :: id_Nbub

        integer(kind=4) :: id_stor_use, id_water_tw, id_CN, id_QN
        integer(kind=4) :: id_Esoil,id_pwater,id_presP,id_methanebP,id_methaneP
  
        !Create the netCDF file.
        CALL check(nf90_create(outFile_restart, NF90_CLOBBER, ncid))
        !Define the dimensions.
        CALL check(nf90_def_dim(ncid, "value",          1, dimid_1))
        call check(nf90_def_dim(ncid, "n_layers", nlayers, dimid_n))
        call check(nf90_def_dim(ncid, "n_pools",        8, dimid_8))

        CALL check(nf90_def_var(ncid, "eJmx0",              NF90_FLOAT, dimid_1, id_eJmx0))
        CALL check(nf90_def_var(ncid, "QC",                 NF90_FLOAT, dimid_8, id_QC))
        CALL check(nf90_def_var(ncid, "Storage",            NF90_FLOAT, dimid_1, id_Storage))
        CALL check(nf90_def_var(ncid, "stor_use",           NF90_FLOAT, dimid_1, id_stor_use))
        CALL check(nf90_def_var(ncid, "nsc",                NF90_FLOAT, dimid_1, id_nsc))
        CALL check(nf90_def_var(ncid, "accumulation",       NF90_FLOAT, dimid_1, id_accumulation))
        CALL check(nf90_def_var(ncid, "SNvcmax",            NF90_FLOAT, dimid_1, id_SNvcmax))
        CALL check(nf90_def_var(ncid, "LAI",                NF90_FLOAT, dimid_1, id_LAI))
        CALL check(nf90_def_var(ncid, "alphaN",             NF90_FLOAT, dimid_1, id_alphaN))
        CALL check(nf90_def_var(ncid, "NSN",                NF90_FLOAT, dimid_1, id_NSN))
        CALL check(nf90_def_var(ncid, "QNminer",            NF90_FLOAT, dimid_1, id_QNminer))
        CALL check(nf90_def_var(ncid, "N_deficit",          NF90_FLOAT, dimid_1, id_N_deficit))
        CALL check(nf90_def_var(ncid, "CN",                 NF90_FLOAT, dimid_8, id_CN))
        CALL check(nf90_def_var(ncid, "QN",                 NF90_FLOAT, dimid_8, id_QN))

        CALL check(nf90_def_var(ncid, "liq_water",          NF90_FLOAT, dimid_n, id_liq_water))
        CALL check(nf90_def_var(ncid, "fwsoil",             NF90_FLOAT, dimid_1, id_fwsoil))
        CALL check(nf90_def_var(ncid, "topfws",             NF90_FLOAT, dimid_1, id_topfws))
        CALL check(nf90_def_var(ncid, "omega",              NF90_FLOAT, dimid_1, id_omega))
        CALL check(nf90_def_var(ncid, "wcl",                NF90_FLOAT, dimid_n, id_wcl))
        CALL check(nf90_def_var(ncid, "zwt",                NF90_FLOAT, dimid_1, id_zwt))
        CALL check(nf90_def_var(ncid, "water_tw",           NF90_FLOAT, dimid_1, id_water_tw))
        
        CALL check(nf90_def_var(ncid, "infilt",             NF90_FLOAT, dimid_1, id_infilt))
        CALL check(nf90_def_var(ncid, "sftmp",              NF90_FLOAT, dimid_1, id_sftmp))
        CALL check(nf90_def_var(ncid, "Tsnow",              NF90_FLOAT, dimid_1, id_Tsnow))
        CALL check(nf90_def_var(ncid, "Twater",             NF90_FLOAT, dimid_1, id_Twater))
        CALL check(nf90_def_var(ncid, "Tice",               NF90_FLOAT, dimid_1, id_Tice))
        CALL check(nf90_def_var(ncid, "G",                  NF90_FLOAT, dimid_1, id_G))
        CALL check(nf90_def_var(ncid, "snow_dsim",          NF90_FLOAT, dimid_1, id_snow_dsim))
        CALL check(nf90_def_var(ncid, "dcount",             NF90_FLOAT, dimid_1, id_dcount))
        CALL check(nf90_def_var(ncid, "dcount_soil",        NF90_FLOAT, dimid_1, id_dcount_soil))
        CALL check(nf90_def_var(ncid, "ice_tw",             NF90_FLOAT, dimid_1, id_ice_tw))
        CALL check(nf90_def_var(ncid, "Tsoill",             NF90_FLOAT, dimid_n, id_Tsoill))
        CALL check(nf90_def_var(ncid, "ice",                NF90_FLOAT, dimid_n, id_ice))
        CALL check(nf90_def_var(ncid, "shcap_snow",         NF90_FLOAT, dimid_1, id_shcap_snow))
        CALL check(nf90_def_var(ncid, "condu_snow",         NF90_FLOAT, dimid_1, id_condu_snow))
        CALL check(nf90_def_var(ncid, "condu_b",            NF90_FLOAT, dimid_1, id_condu_b))
        CALL check(nf90_def_var(ncid, "depth_ex",           NF90_FLOAT, dimid_1, id_depth_ex))
        CALL check(nf90_def_var(ncid, "diff_s",             NF90_FLOAT, dimid_1, id_diff_s))
        CALL check(nf90_def_var(ncid, "diff_snow",          NF90_FLOAT, dimid_1, id_diff_snow))
        CALL check(nf90_def_var(ncid, "albedo_snow",        NF90_FLOAT, dimid_1, id_albedo_snow))
        CALL check(nf90_def_var(ncid, "thd_snow_depth",     NF90_FLOAT, dimid_1, id_thd_snow_depth))
        CALL check(nf90_def_var(ncid, "b_bound",            NF90_FLOAT, dimid_1, id_b_bound))
        CALL check(nf90_def_var(ncid, "infilt_rate",        NF90_FLOAT, dimid_1, id_infilt_rate))
        CALL check(nf90_def_var(ncid, "fa",                 NF90_FLOAT, dimid_1, id_fa))
        CALL check(nf90_def_var(ncid, "fsub",               NF90_FLOAT, dimid_1, id_fsub))
        CALL check(nf90_def_var(ncid, "rho_snow",           NF90_FLOAT, dimid_1, id_rho_snow))
        CALL check(nf90_def_var(ncid, "CH4_V",              NF90_FLOAT, dimid_n, id_CH4_V))
        CALL check(nf90_def_var(ncid, "CH4",                NF90_FLOAT, dimid_n, id_CH4))
        CALL check(nf90_def_var(ncid, "Vp",                 NF90_FLOAT, dimid_n, id_Vp))
        CALL check(nf90_def_var(ncid, "bubble_methane_tot", NF90_FLOAT, dimid_1, id_bubble_methane_tot))
        CALL check(nf90_def_var(ncid, "Nbub",               NF90_FLOAT, dimid_1, id_Nbub))

        ! id_Esoil,id_pwater,id_presP,id_methanebP,id_methaneP
        CALL check(nf90_def_var(ncid, "Esoil",              NF90_FLOAT, dimid_1, id_Esoil))
        CALL check(nf90_def_var(ncid, "pwater",             NF90_FLOAT, dimid_n, id_pwater))
        CALL check(nf90_def_var(ncid, "presP",              NF90_FLOAT, dimid_n, id_presP))
        CALL check(nf90_def_var(ncid, "methanebP",          NF90_FLOAT, dimid_n, id_methanebP))
        CALL check(nf90_def_var(ncid, "methaneP",           NF90_FLOAT, dimid_n, id_methaneP))

        CALL check(nf90_enddef(ncid)) !End Definitions

        CALL check(nf90_put_var(ncid, id_eJmx0,              eJmx0))
        CALL check(nf90_put_var(ncid, id_QC,                 QC))
        CALL check(nf90_put_var(ncid, id_Storage,            Storage))
        CALL check(nf90_put_var(ncid, id_stor_use,           stor_use))
        CALL check(nf90_put_var(ncid, id_nsc,                nsc))
        CALL check(nf90_put_var(ncid, id_accumulation,       accumulation))
        CALL check(nf90_put_var(ncid, id_SNvcmax,            SNvcmax))
        CALL check(nf90_put_var(ncid, id_LAI,                LAI))
        CALL check(nf90_put_var(ncid, id_alphaN,             alphaN))
        CALL check(nf90_put_var(ncid, id_NSN,                NSN))
        CALL check(nf90_put_var(ncid, id_QNminer,            QNminer))
        CALL check(nf90_put_var(ncid, id_N_deficit,          N_deficit))
        CALL check(nf90_put_var(ncid, id_CN,                 CN))
        CALL check(nf90_put_var(ncid, id_QN,                 QN))

        CALL check(nf90_put_var(ncid, id_liq_water,          liq_water))
        CALL check(nf90_put_var(ncid, id_fwsoil,             fwsoil))
        CALL check(nf90_put_var(ncid, id_topfws,             topfws))
        CALL check(nf90_put_var(ncid, id_omega,              omega))
        CALL check(nf90_put_var(ncid, id_wcl,                wcl))
        CALL check(nf90_put_var(ncid, id_zwt,                zwt))
        CALL check(nf90_put_var(ncid, id_water_tw,           water_tw))
        
        CALL check(nf90_put_var(ncid, id_infilt,             infilt))
        CALL check(nf90_put_var(ncid, id_sftmp,              sftmp))
        CALL check(nf90_put_var(ncid, id_Tsnow,              Tsnow))
        CALL check(nf90_put_var(ncid, id_Twater,             Twater))
        CALL check(nf90_put_var(ncid, id_Tice,               Tice))
        CALL check(nf90_put_var(ncid, id_G,                  G))
        CALL check(nf90_put_var(ncid, id_snow_dsim,          snow_dsim))
        CALL check(nf90_put_var(ncid, id_dcount,             dcount))
        CALL check(nf90_put_var(ncid, id_dcount_soil,        dcount_soil))
        CALL check(nf90_put_var(ncid, id_ice_tw,             ice_tw))
        CALL check(nf90_put_var(ncid, id_Tsoill,             Tsoill))
        CALL check(nf90_put_var(ncid, id_ice,                ice))
        CALL check(nf90_put_var(ncid, id_shcap_snow,         shcap_snow))
        CALL check(nf90_put_var(ncid, id_condu_snow,         condu_snow))
        CALL check(nf90_put_var(ncid, id_condu_b,            condu_b))
        CALL check(nf90_put_var(ncid, id_depth_ex,           depth_ex))
        CALL check(nf90_put_var(ncid, id_diff_s,             diff_s))
        CALL check(nf90_put_var(ncid, id_diff_snow,          diff_snow))
        CALL check(nf90_put_var(ncid, id_albedo_snow,        albedo_snow))
        CALL check(nf90_put_var(ncid, id_thd_snow_depth,     thd_snow_depth))
        CALL check(nf90_put_var(ncid, id_b_bound,            b_bound))
        CALL check(nf90_put_var(ncid, id_infilt_rate,        infilt_rate))
        CALL check(nf90_put_var(ncid, id_fa,                 fa))
        CALL check(nf90_put_var(ncid, id_fsub,               fsub))
        CALL check(nf90_put_var(ncid, id_rho_snow,           rho_snow))
        CALL check(nf90_put_var(ncid, id_CH4_V,              CH4_V))
        CALL check(nf90_put_var(ncid, id_CH4,                CH4))
        CALL check(nf90_put_var(ncid, id_Vp,                 Vp))
        CALL check(nf90_put_var(ncid, id_bubble_methane_tot, bubble_methane_tot))
        CALL check(nf90_put_var(ncid, id_Nbub,               Nbub))

        CALL check(nf90_put_var(ncid, id_Esoil,              Esoil))
        CALL check(nf90_put_var(ncid, id_pwater,             pwater))
        CALL check(nf90_put_var(ncid, id_presP,              presP))
        CALL check(nf90_put_var(ncid, id_methanebP,          methanebP))
        CALL check(nf90_put_var(ncid, id_methaneP,           methaneP))
        ! endif
        CALL check(nf90_close(ncid))
    end subroutine write_restart

    !------------------------------------------------
    subroutine read_restart(infile)
        implicit none
        CHARACTER(LEN=*), INTENT(IN) :: infile
        integer(KIND=4) :: ncid
        INTEGER(KIND=4) :: varid

        !Open netCDF file
        !:-------:-------:-------:-------:-------:-------:-------:-------:
        CALL check(nf90_open(infile, nf90_nowrite, ncid))

        CALL check(nf90_inq_varid(ncid, "eJmx0", varid))
        call check(nf90_get_var(ncid,     varid, eJmx0))

        CALL check(nf90_inq_varid(ncid,  "QC", varid))
        call check(nf90_get_var(ncid,   varid,    QC))

        CALL check(nf90_inq_varid(ncid, "Storage",varid))
        call check(nf90_get_var(ncid, varid,            Storage))

        CALL check(nf90_inq_varid(ncid, "nsc",varid))
        call check(nf90_get_var(ncid, varid,                nsc))

        CALL check(nf90_inq_varid(ncid, "accumulation",varid))
        call check(nf90_get_var(ncid, varid,       accumulation))

        CALL check(nf90_inq_varid(ncid, "SNvcmax",varid))
        call check(nf90_get_var(ncid, varid,            SNvcmax))

        CALL check(nf90_inq_varid(ncid, "LAI",varid))
        call check(nf90_get_var(ncid, varid,                LAI))

        CALL check(nf90_inq_varid(ncid, "alphaN",varid))
        call check(nf90_get_var(ncid, varid,             alphaN))

        CALL check(nf90_inq_varid(ncid, "NSN",varid))
        call check(nf90_get_var(ncid, varid,                NSN))

        CALL check(nf90_inq_varid(ncid, "QNminer",varid))
        call check(nf90_get_var(ncid, varid,            QNminer))

        CALL check(nf90_inq_varid(ncid, "N_deficit",varid))
        call check(nf90_get_var(ncid, varid,          N_deficit))

        CALL check(nf90_inq_varid(ncid, "liq_water",varid))
        call check(nf90_get_var(ncid, varid,          liq_water))

        CALL check(nf90_inq_varid(ncid, "fwsoil",varid))
        call check(nf90_get_var(ncid, varid,             fwsoil))

        CALL check(nf90_inq_varid(ncid, "topfws",varid))
        call check(nf90_get_var(ncid, varid,             topfws))

        CALL check(nf90_inq_varid(ncid, "omega",varid))
        call check(nf90_get_var(ncid, varid,              omega))

        CALL check(nf90_inq_varid(ncid, "wcl",varid))
        call check(nf90_get_var(ncid, varid,                wcl))

        CALL check(nf90_inq_varid(ncid, "zwt",varid))
        call check(nf90_get_var(ncid, varid,                zwt))

        CALL check(nf90_inq_varid(ncid, "infilt",varid))
        call check(nf90_get_var(ncid, varid,             infilt))

        CALL check(nf90_inq_varid(ncid, "sftmp",varid))
        call check(nf90_get_var(ncid, varid,              sftmp))

        CALL check(nf90_inq_varid(ncid, "Tsnow",varid))
        call check(nf90_get_var(ncid, varid,              Tsnow))

        CALL check(nf90_inq_varid(ncid, "Twater",varid))
        call check(nf90_get_var(ncid, varid,             Twater))

        CALL check(nf90_inq_varid(ncid, "Tice",varid))
        call check(nf90_get_var(ncid, varid,               Tice))

        CALL check(nf90_inq_varid(ncid, "G",varid))
        call check(nf90_get_var(ncid, varid,                  G))

        CALL check(nf90_inq_varid(ncid, "snow_dsim",varid))
        call check(nf90_get_var(ncid, varid,          snow_dsim))

        CALL check(nf90_inq_varid(ncid, "dcount",varid))
        call check(nf90_get_var(ncid, varid,             dcount))

        CALL check(nf90_inq_varid(ncid, "dcount_soil",varid))
        call check(nf90_get_var(ncid, varid,        dcount_soil))

        CALL check(nf90_inq_varid(ncid, "ice_tw",varid))
        call check(nf90_get_var(ncid, varid,             ice_tw))

        CALL check(nf90_inq_varid(ncid, "Tsoill",varid))
        call check(nf90_get_var(ncid, varid,             Tsoill))

        CALL check(nf90_inq_varid(ncid, "ice",varid))
        call check(nf90_get_var(ncid, varid,                ice))

        CALL check(nf90_inq_varid(ncid, "shcap_snow",varid))
        call check(nf90_get_var(ncid, varid,         shcap_snow))

        CALL check(nf90_inq_varid(ncid, "condu_snow",varid))
        call check(nf90_get_var(ncid, varid,         condu_snow))

        CALL check(nf90_inq_varid(ncid, "condu_b",varid))
        call check(nf90_get_var(ncid, varid,            condu_b))

        CALL check(nf90_inq_varid(ncid, "depth_ex",varid))
        call check(nf90_get_var(ncid, varid,           depth_ex))

        CALL check(nf90_inq_varid(ncid, "diff_s",varid))
        call check(nf90_get_var(ncid, varid,             diff_s))

        CALL check(nf90_inq_varid(ncid, "diff_snow",varid))
        call check(nf90_get_var(ncid, varid,          diff_snow))

        CALL check(nf90_inq_varid(ncid, "albedo_snow",varid))
        call check(nf90_get_var(ncid, varid,        albedo_snow))

        CALL check(nf90_inq_varid(ncid, "thd_snow_depth",varid))
        call check(nf90_get_var(ncid, varid,     thd_snow_depth))

        CALL check(nf90_inq_varid(ncid, "b_bound",varid))
        call check(nf90_get_var(ncid, varid,            b_bound))

        CALL check(nf90_inq_varid(ncid, "infilt_rate",varid))
        call check(nf90_get_var(ncid, varid,        infilt_rate))

        CALL check(nf90_inq_varid(ncid, "fa",varid))
        call check(nf90_get_var(ncid, varid,                 fa))

        CALL check(nf90_inq_varid(ncid, "fsub",varid))
        call check(nf90_get_var(ncid, varid,               fsub))

        CALL check(nf90_inq_varid(ncid, "rho_snow",varid))
        call check(nf90_get_var(ncid, varid,           rho_snow))

        CALL check(nf90_inq_varid(ncid, "CH4_V",varid))
        call check(nf90_get_var(ncid, varid,              CH4_V))

        CALL check(nf90_inq_varid(ncid, "CH4",varid))
        call check(nf90_get_var(ncid, varid,                CH4))

        CALL check(nf90_inq_varid(ncid, "Vp",varid))
        call check(nf90_get_var(ncid, varid,                 Vp))

        CALL check(nf90_inq_varid(ncid, "bubble_methane_tot",varid))
        call check(nf90_get_var(ncid, varid, bubble_methane_tot))

        CALL check(nf90_inq_varid(ncid, "Nbub",varid))
        call check(nf90_get_var(ncid, varid,               Nbub))

        !----------------------------------------------------------------
        CALL check(nf90_inq_varid(ncid, "stor_use", varid))
        call check(nf90_get_var(ncid,   varid,      stor_use))

        CALL check(nf90_inq_varid(ncid, "water_tw", varid))
        call check(nf90_get_var(ncid, varid,        water_tw))

        CALL check(nf90_inq_varid(ncid, "CN",varid))
        call check(nf90_get_var(ncid, varid,  CN))
        
        CALL check(nf90_inq_varid(ncid, "QN",varid))
        call check(nf90_get_var(ncid, varid,  QN))

        CALL check(nf90_inq_varid(ncid, "Esoil",varid))
        call check(nf90_get_var(ncid, varid,    Esoil))
        
        CALL check(nf90_inq_varid(ncid, "pwater",varid))
        call check(nf90_get_var(ncid, varid,    pwater))
        
        CALL check(nf90_inq_varid(ncid, "presP",varid))
        call check(nf90_get_var(ncid, varid,    presP))
        
        CALL check(nf90_inq_varid(ncid, "methanebP",varid))
        call check(nf90_get_var(ncid, varid,    methanebP))
        
        CALL check(nf90_inq_varid(ncid, "methaneP",varid))
        call check(nf90_get_var(ncid, varid,    methaneP))


        ! CALL check(nf90_inq_varid(ncid, "eJmx0",              varid))
        ! call check( nf90_get_var(ncid, varid, eJmx0) )
        ! CALL check(nf90_inq_varid(ncid, "QC",                 varid))
        ! call check( nf90_get_var(ncid, varid, eJmx0) )
        ! CALL check(nf90_inq_varid(ncid, "Storage",            varid))
        ! CALL check(nf90_inq_varid(ncid, "nsc",                varid))
        ! CALL check(nf90_inq_varid(ncid, "accumulation",       varid))
        ! CALL check(nf90_inq_varid(ncid, "SNvcmax",            varid))
        ! CALL check(nf90_inq_varid(ncid, "LAI",                varid))
        ! CALL check(nf90_inq_varid(ncid, "alphaN",             varid))
        ! CALL check(nf90_inq_varid(ncid, "NSN",                varid))
        ! CALL check(nf90_inq_varid(ncid, "QNminer",            varid))
        ! CALL check(nf90_inq_varid(ncid, "N_deficit",          varid))
        ! CALL check(nf90_inq_varid(ncid, "liq_water",          varid))
        ! CALL check(nf90_inq_varid(ncid, "fwsoil",             varid))
        ! CALL check(nf90_inq_varid(ncid, "topfws",             varid))
        ! CALL check(nf90_inq_varid(ncid, "omega",              varid))
        ! CALL check(nf90_inq_varid(ncid, "wcl",                varid))
        ! CALL check(nf90_inq_varid(ncid, "zwt",                varid))
        ! CALL check(nf90_inq_varid(ncid, "infilt",             varid))
        ! CALL check(nf90_inq_varid(ncid, "sftmp",              varid))
        ! CALL check(nf90_inq_varid(ncid, "Tsnow",              varid))
        ! CALL check(nf90_inq_varid(ncid, "Twater",             varid))
        ! CALL check(nf90_inq_varid(ncid, "Tice",               varid))
        ! CALL check(nf90_inq_varid(ncid, "G",                  varid))
        ! CALL check(nf90_inq_varid(ncid, "snow_dsim",          varid))
        ! CALL check(nf90_inq_varid(ncid, "dcount",             varid))
        ! CALL check(nf90_inq_varid(ncid, "dcount_soil",        varid))
        ! CALL check(nf90_inq_varid(ncid, "ice_tw",             varid))
        ! CALL check(nf90_inq_varid(ncid, "Tsoill",             varid))
        ! CALL check(nf90_inq_varid(ncid, "ice",                varid))
        ! CALL check(nf90_inq_varid(ncid, "shcap_snow",         varid))
        ! CALL check(nf90_inq_varid(ncid, "condu_snow",         varid))
        ! CALL check(nf90_inq_varid(ncid, "condu_b",            varid))
        ! CALL check(nf90_inq_varid(ncid, "depth_ex",           varid))
        ! CALL check(nf90_inq_varid(ncid, "diff_s",             varid))
        ! CALL check(nf90_inq_varid(ncid, "diff_snow",          varid))
        ! CALL check(nf90_inq_varid(ncid, "albedo_snow",        varid))
        ! CALL check(nf90_inq_varid(ncid, "thd_snow_depth",     varid))
        ! CALL check(nf90_inq_varid(ncid, "b_bound",            varid))
        ! CALL check(nf90_inq_varid(ncid, "infilt_rate",        varid))
        ! CALL check(nf90_inq_varid(ncid, "fa",                 varid))
        ! CALL check(nf90_inq_varid(ncid, "fsub",               varid))
        ! CALL check(nf90_inq_varid(ncid, "rho_snow",           varid))
        ! CALL check(nf90_inq_varid(ncid, "CH4_V",              varid))
        ! CALL check(nf90_inq_varid(ncid, "CH4",                varid))
        ! CALL check(nf90_inq_varid(ncid, "Vp",                 varid))
        ! CALL check(nf90_inq_varid(ncid, "bubble_methane_tot", varid))
        ! CALL check(nf90_inq_varid(ncid, "Nbub",               varid))

        !Close netCDF file
        !:-------:-------:-------:-------:-------:-------:-------:-------:
        CALL check(nf90_close(ncid))
        !:=========================================================================
    end subroutine read_restart

    ! ----------------------------------------------
    subroutine write_spinup_res()
        implicit none
        integer(KIND=4) :: ncid, dimid_nloops, dimid_nlayer,dimid_ntest
        integer(kind=4) :: id_sp_gpp
        integer(kind=4) :: id_sp_npp
        integer(kind=4) :: id_sp_ra
        integer(kind=4) :: id_sp_rh
        integer(kind=4) :: id_sp_wetlandCH4
        integer(kind=4) :: id_sp_wetlandCH4prod
        integer(kind=4) :: id_sp_wetlandCH4cons
        ! Carbon Pools  (KgC m-2)
        integer(kind=4) :: id_sp_cLeaf
        integer(kind=4) :: id_sp_cStem
        integer(kind=4) :: id_sp_cRoot
        integer(kind=4) :: id_sp_cOther        
        integer(kind=4) :: id_sp_cLitter
        integer(kind=4) :: id_sp_cLitterCwd      
        integer(kind=4) :: id_sp_cSoil
        integer(kind=4) :: id_sp_cSoilFast
        integer(kind=4) :: id_sp_cSoilSlow
        integer(kind=4) :: id_sp_cSoilPassive
        integer(kind=4) :: id_sp_cCH4
        ! Nitrogen fluxes (kgN m-2 s-1)
        integer(kind=4) :: id_sp_fBNF
        integer(kind=4) :: id_sp_fN2O
        integer(kind=4) :: id_sp_fNloss
        integer(kind=4) :: id_sp_fNnetmin
        integer(kind=4) :: id_sp_fNdep 
        ! Nitrogen pools (kgN m-2)
        integer(kind=4) :: id_sp_nLeaf
        integer(kind=4) :: id_sp_nStem
        integer(kind=4) :: id_sp_nRoot
        integer(kind=4) :: id_sp_nOther
        integer(kind=4) :: id_sp_nLitter
        integer(kind=4) :: id_sp_nLitterCwd
        integer(kind=4) :: id_sp_nSoil
        integer(kind=4) :: id_sp_nMineral
        ! energy fluxes (W m-2)
        integer(kind=4) :: id_sp_hfls 
        integer(kind=4) :: id_sp_hfss
        ! water fluxes (kg m-2 s-1)
        integer(kind=4) :: id_sp_ec
        integer(kind=4) :: id_sp_tran
        integer(kind=4) :: id_sp_es
        integer(kind=4) :: id_sp_hfsbl
        integer(kind=4) :: id_sp_mrro
        integer(kind=4) :: id_sp_mrros
        integer(kind=4) :: id_sp_mrrob
        integer(kind=4) :: id_sp_lai
        integer(kind=4) :: id_test

        !Create the netCDF file.
        CALL check(nf90_create(outFile_sp, NF90_CLOBBER, ncid))
        !Define the dimensions.
        CALL check(nf90_def_dim(ncid, "time",     nloops, dimid_nloops))
        call check(nf90_def_dim(ncid, "n_layers", nlayers, dimid_nlayer))
        call check(nf90_def_dim(ncid, "n_Test",   9, dimid_ntest))
        ! Define the variables

        CALL check(nf90_def_var(ncid, "gpp",            NF90_FLOAT, dimid_nloops, id_sp_gpp))
        CALL check(nf90_def_var(ncid, "npp",            NF90_FLOAT, dimid_nloops, id_sp_npp))
        CALL check(nf90_def_var(ncid, "ra",             NF90_FLOAT, dimid_nloops, id_sp_ra))
        CALL check(nf90_def_var(ncid, "rh",             NF90_FLOAT, dimid_nloops, id_sp_rh))
        CALL check(nf90_def_var(ncid, "wetlandCH4",     NF90_FLOAT, dimid_nloops, id_sp_wetlandCH4))
        CALL check(nf90_def_var(ncid, "wetlandCH4prod", NF90_FLOAT, dimid_nloops, id_sp_wetlandCH4prod))
        CALL check(nf90_def_var(ncid, "wetlandCH4cons", NF90_FLOAT, dimid_nloops, id_sp_wetlandCH4cons))
        ! Carbon Pools  (KgC m-2)
        CALL check(nf90_def_var(ncid, "cLeaf",          NF90_FLOAT, dimid_nloops, id_sp_cLeaf))
        CALL check(nf90_def_var(ncid, "cStem",          NF90_FLOAT, dimid_nloops, id_sp_cStem))
        CALL check(nf90_def_var(ncid, "cRoot",          NF90_FLOAT, dimid_nloops, id_sp_cRoot))
        CALL check(nf90_def_var(ncid, "cOther",         NF90_FLOAT, dimid_nloops, id_sp_cOther))
        CALL check(nf90_def_var(ncid, "cLitter",        NF90_FLOAT, dimid_nloops, id_sp_cLitter))
        CALL check(nf90_def_var(ncid, "cLitterCwd",     NF90_FLOAT, dimid_nloops, id_sp_cLitterCwd))
        CALL check(nf90_def_var(ncid, "cSoil",          NF90_FLOAT, dimid_nloops, id_sp_cSoil))
        CALL check(nf90_def_var(ncid, "cSoilFast",      NF90_FLOAT, dimid_nloops, id_sp_cSoilFast))
        CALL check(nf90_def_var(ncid, "cSoilSlow",      NF90_FLOAT, dimid_nloops, id_sp_cSoilSlow))
        CALL check(nf90_def_var(ncid, "cSoilPassive",   NF90_FLOAT, dimid_nloops, id_sp_cSoilPassive))
        CALL check(nf90_def_var(ncid, "cCH4",           NF90_FLOAT, (/dimid_nloops,dimid_nlayer/), id_sp_cCH4))
        ! Nitrogen fluxes (kgN m-2 s-1)
        CALL check(nf90_def_var(ncid, "fBNF",           NF90_FLOAT, dimid_nloops, id_sp_fBNF))
        CALL check(nf90_def_var(ncid, "fN2O",           NF90_FLOAT, dimid_nloops, id_sp_fN2O))
        CALL check(nf90_def_var(ncid, "fNloss",         NF90_FLOAT, dimid_nloops, id_sp_fNloss))
        CALL check(nf90_def_var(ncid, "fNnetmin",       NF90_FLOAT, dimid_nloops, id_sp_fNnetmin))
        CALL check(nf90_def_var(ncid, "fNdep",          NF90_FLOAT, dimid_nloops, id_sp_fNdep))
        ! Nitrogen pools (kgN m-2)
        CALL check(nf90_def_var(ncid, "nLeaf",          NF90_FLOAT, dimid_nloops, id_sp_nLeaf))
        CALL check(nf90_def_var(ncid, "nStem",          NF90_FLOAT, dimid_nloops, id_sp_nStem))
        CALL check(nf90_def_var(ncid, "nRoot",          NF90_FLOAT, dimid_nloops, id_sp_nRoot))
        CALL check(nf90_def_var(ncid, "nOther",         NF90_FLOAT, dimid_nloops, id_sp_nOther))
        CALL check(nf90_def_var(ncid, "nLitter",        NF90_FLOAT, dimid_nloops, id_sp_nLitter))
        CALL check(nf90_def_var(ncid, "nLitterCwd",     NF90_FLOAT, dimid_nloops, id_sp_nLitterCwd))
        CALL check(nf90_def_var(ncid, "nSoil",          NF90_FLOAT, dimid_nloops, id_sp_nSoil))
        CALL check(nf90_def_var(ncid, "nMineral",       NF90_FLOAT, dimid_nloops, id_sp_nMineral))
        ! energy fluxes (W m-2)
        CALL check(nf90_def_var(ncid, "hfls",           NF90_FLOAT, dimid_nloops, id_sp_hfls))
        CALL check(nf90_def_var(ncid, "hfss",           NF90_FLOAT, dimid_nloops, id_sp_hfss))
        ! water fluxes (kg m-2 s-1)
        CALL check(nf90_def_var(ncid, "ec",             NF90_FLOAT, dimid_nloops, id_sp_ec))
        CALL check(nf90_def_var(ncid, "tran",           NF90_FLOAT, dimid_nloops, id_sp_tran))
        CALL check(nf90_def_var(ncid, "es",             NF90_FLOAT, dimid_nloops, id_sp_es))
        CALL check(nf90_def_var(ncid, "hfsbl",          NF90_FLOAT, dimid_nloops, id_sp_hfsbl))
        CALL check(nf90_def_var(ncid, "mrro",           NF90_FLOAT, dimid_nloops, id_sp_mrro))
        CALL check(nf90_def_var(ncid, "mrros",          NF90_FLOAT, dimid_nloops, id_sp_mrros))
        CALL check(nf90_def_var(ncid, "mrrob",          NF90_FLOAT, dimid_nloops, id_sp_mrrob))
        CALL check(nf90_def_var(ncid, "lai",            NF90_FLOAT, dimid_nloops, id_sp_lai))
        CALL check(nf90_def_var(ncid, "test",           NF90_FLOAT, (/dimid_nloops,dimid_ntest/), id_test))

        CALL check(nf90_enddef(ncid)) !End Definitions

        CALL check(nf90_put_var(ncid, id_sp_gpp,            sp_gpp_y))
        CALL check(nf90_put_var(ncid, id_sp_npp,            sp_npp_y))
        CALL check(nf90_put_var(ncid, id_sp_ra,             sp_ra_y))
        CALL check(nf90_put_var(ncid, id_sp_rh,             sp_rh_y))
        CALL check(nf90_put_var(ncid, id_sp_wetlandCH4,     sp_wetlandCH4_y))
        CALL check(nf90_put_var(ncid, id_sp_wetlandCH4prod, sp_wetlandCH4prod_y))
        CALL check(nf90_put_var(ncid, id_sp_wetlandCH4cons, sp_wetlandCH4cons_y))
        ! Carbon Pools  (KgC m-2)
        CALL check(nf90_put_var(ncid, id_sp_cLeaf,          sp_cLeaf_y))
        CALL check(nf90_put_var(ncid, id_sp_cStem,          sp_cStem_y))
        CALL check(nf90_put_var(ncid, id_sp_cRoot,          sp_cRoot_y))
        CALL check(nf90_put_var(ncid, id_sp_cOther,         sp_cOther_y))
        CALL check(nf90_put_var(ncid, id_sp_cLitter,        sp_cLitter_y))
        CALL check(nf90_put_var(ncid, id_sp_cLitterCwd,     sp_cLitterCwd_y))
        CALL check(nf90_put_var(ncid, id_sp_cSoil,          sp_cSoil_y))
        CALL check(nf90_put_var(ncid, id_sp_cSoilFast,      sp_cSoilFast_y))
        CALL check(nf90_put_var(ncid, id_sp_cSoilSlow,      sp_cSoilSlow_y))
        CALL check(nf90_put_var(ncid, id_sp_cSoilPassive,   sp_cSoilPassive_y))
        CALL check(nf90_put_var(ncid, id_sp_cCH4,           sp_cCH4_y))
        ! Nitrogen fluxes (kgN m-2 s-1)
        CALL check(nf90_put_var(ncid, id_sp_fBNF,           sp_fBNF_y))
        CALL check(nf90_put_var(ncid, id_sp_fN2O,           sp_fN2O_y))
        CALL check(nf90_put_var(ncid, id_sp_fNloss,         sp_fNloss_y))
        CALL check(nf90_put_var(ncid, id_sp_fNnetmin,       sp_fNnetmin_y))
        CALL check(nf90_put_var(ncid, id_sp_fNdep,          sp_fNdep_y))
        ! Nitrogen pools (kgN m-2)
        CALL check(nf90_put_var(ncid, id_sp_nLeaf,          sp_nLeaf_y))
        CALL check(nf90_put_var(ncid, id_sp_nStem,          sp_nStem_y))
        CALL check(nf90_put_var(ncid, id_sp_nRoot,          sp_nRoot_y))
        CALL check(nf90_put_var(ncid, id_sp_nOther,         sp_nOther_y))
        CALL check(nf90_put_var(ncid, id_sp_nLitter,        sp_nLitter_y))
        CALL check(nf90_put_var(ncid, id_sp_nLitterCwd,     sp_nLitterCwd_y))
        CALL check(nf90_put_var(ncid, id_sp_nSoil,          sp_nSoil_y))
        CALL check(nf90_put_var(ncid, id_sp_nMineral,       sp_nMineral_y))
        ! energy fluxes (W m-2)
        CALL check(nf90_put_var(ncid, id_sp_hfls,           sp_hfls_y))
        CALL check(nf90_put_var(ncid, id_sp_hfss,           sp_hfss_y))
        ! water fluxes (kg m-2 s-1)
        CALL check(nf90_put_var(ncid, id_sp_ec,             sp_ec_y))
        CALL check(nf90_put_var(ncid, id_sp_tran,           sp_tran_y))
        CALL check(nf90_put_var(ncid, id_sp_es,             sp_es_y))
        CALL check(nf90_put_var(ncid, id_sp_hfsbl,          sp_hfsbl_y))
        CALL check(nf90_put_var(ncid, id_sp_mrro,           sp_mrro_y))
        CALL check(nf90_put_var(ncid, id_sp_mrros,          sp_mrros_y))
        CALL check(nf90_put_var(ncid, id_sp_mrrob,          sp_mrrob_y))
        CALL check(nf90_put_var(ncid, id_sp_lai,            sp_lai_y))
        CALL check(nf90_put_var(ncid, id_test,              sp_test_y))
        ! endif
        CALL check(nf90_close(ncid))
    end subroutine write_spinup_res

    subroutine write_nc(outfile, lenTime, data, varName, unit, description, freq, nSoilLayer)
        IMPLICIT NONE
        real(kind=4), Dimension(lenTime,nSoilLayer), intent(in) :: data
        integer(kind=4) :: nSoilLayer
        integer(KIND=4) :: ncid, timid, dp_dimid, timvarid
        integer(kind=4) :: varid
        integer(kind=4), intent(in) :: lenTime
        CHARACTER(LEN=*), INTENT(IN) :: outfile, freq
        CHARACTER(len=*), intent(in) :: varName, unit, description
        character(len=:), allocatable :: nc_fileName
        character(len=100) :: timeUnit
        integer itime
        real, dimension(lenTime) :: time_values 
        integer :: start(1), count(1)
        
        allocate(character(len=200+len(outfile)) :: nc_fileName)
        nc_fileName = adjustl(trim(outfile))//"/"//adjustl(trim(varName))//"_"//freq//"_TECO-SPRUCE_"//&
            & adjustl(trim(experiment))//"_"//adjustl(trim(str_startyr))//"-"//adjustl(trim(str_endyr))//".nc"   
        
        !Create the netCDF file.
        CALL check(nf90_create(nc_fileName, NF90_CLOBBER, ncid))

        !Define the dimensions.
        CALL check(nf90_def_dim(ncid, "time", lenTime, timid))
    
        if (nSoilLayer>1)then
            call check(nf90_def_dim(ncid, "depth", nSoilLayer, dp_dimid))
            CALL check(nf90_def_var(ncid = ncid, name = varName,  xtype = NF90_FLOAT, &
                & dimids = (/timid, dp_dimid/),  varID =varid))
        else
            CALL check(nf90_def_var(ncid, varName, NF90_FLOAT, timid, varid))
        endif
        call check(nf90_def_var(ncid, "time", NF90_DOUBLE, timid, timvarid))
        !Define data variable
        
        !Add attributes
        if (freq .eq. "hourly") then
            timeUnit = "hours since "//adjustl(trim(str_startyr))//"-01-01 00:00:00"
        else if (freq .eq. "daily") then
            timeUnit = "days since "//adjustl(trim(str_startyr))//"-01-01 00:00:00"
        else if (freq .eq. "monthly") then
            timeUnit = "months since "//adjustl(trim(str_startyr))//"-01-01 00:00:00"
        end if
        
        call check(nf90_put_att(ncid,timvarid,"units",adjustl(trim(timeUnit))))
        CALL check(nf90_put_att(ncid,varid,"units",unit))
        CALL check(nf90_put_att(ncid,varid,"description",description))
        CALL check(nf90_enddef(ncid)) 
        !End Definitions

        !Write Data
        ! if (nSoilLayer>1)then
        !     do i = 1, nSoilLayer
        !         CALL check(nf90_put_var(ncid, varid, data, start=[1,i], count=[lenTime,1]))
        !     enddo
        ! else

        do itime = 1, lenTime
            time_values(itime) = itime-1
        enddo
        start = 1
        count = lenTime
        CALL check(nf90_put_var(ncid, timvarid, time_values,start,count))
        CALL check(nf90_put_var(ncid, varid, data))
        
        CALL check(nf90_close(ncid))
    end subroutine write_nc


    
    ! check (ever so slightly modified from www.unidata.ucar.edu)
    subroutine check(istatus)
        ! use netcdf
        implicit none
        integer, intent(in) :: istatus
        if(istatus /= nf90_noerr) then
            write(*,*) trim(adjustl(nf90_strerror(istatus)))
        end if
    end subroutine check
end module mod_ncd_io