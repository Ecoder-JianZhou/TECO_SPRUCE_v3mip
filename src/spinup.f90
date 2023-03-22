module mod_spinup
    use mod_data
    use driver

    !-------------------------------------------------------------------------------------------------------
    integer iloop

    contains
    subroutine run_spinup()
        implicit none
        write(*,*)"This is spinup", nloops
        do iloop = 1, nloops
            write(*,*) "iloop: ", iloop
            call teco_simu()
            call update_spinup_values()
        enddo
        
    end subroutine run_spinup
    
    !-----------------------------------------------------------------------
    subroutine update_spinup_values()
        implicit none 
        sp_gpp_y(iloop)            = gpp_y
        sp_npp_y(iloop)            = npp_y
        sp_ra_y(iloop)             = ra_y
        sp_rh_y(iloop)             = rh_y
        sp_wetlandCH4_y(iloop)     = wetlandCH4_y
        sp_wetlandCH4prod_y(iloop) = wetlandCH4prod_y
        sp_wetlandCH4cons_y(iloop) = wetlandCH4cons_y
        ! Carbon Pools  (KgC m-2)
        sp_cLeaf_y(iloop)          = cLeaf_y
        sp_cStem_y(iloop)          = cStem_y
        sp_cRoot_y(iloop)          = cRoot_y
        sp_cOther_y(iloop)         = cOther_y             
        sp_cLitter_y(iloop)        = cLitter_y
        sp_cLitterCwd_y(iloop)     = cLitterCwd_y                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        sp_cSoil_y(iloop)          = cSoil_y
        sp_cSoilFast_y(iloop)      = cSoilFast_y
        sp_cSoilSlow_y(iloop)      = cSoilSlow_y
        sp_cSoilPassive_y(iloop)   = cSoilPassive_y                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        sp_cCH4_y(iloop,:)         = cCH4_y                                                          ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        sp_fBNF_y(iloop)           = fBNF_y
        sp_fN2O_y(iloop)           = fN2O_y
        sp_fNloss_y(iloop)         = fNloss_y
        sp_fNnetmin_y(iloop)       = fNnetmin_y
        sp_fNdep_y(iloop)          = fNdep_y                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        sp_nLeaf_y(iloop)          = nLeaf_y
        sp_nStem_y(iloop)          = nStem_y
        sp_nRoot_y(iloop)          = nRoot_y
        sp_nOther_y(iloop)         = nOther_y
        sp_nLitter_y(iloop)        = nLitter_y
        sp_nLitterCwd_y(iloop)     = nLitterCwd_y
        sp_nSoil_y(iloop)          = nSoil_y
        sp_nMineral_y(iloop)       = nMineral_y                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        sp_hfls_y(iloop)           = hfls_y
        sp_hfss_y(iloop)           = hfss_y                                                       ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        sp_ec_y(iloop)             = ec_y
        sp_tran_y(iloop)           = tran_y
        sp_es_y(iloop)             = es_y                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        sp_hfsbl_y(iloop)          = hfsbl_y                                                         ! Snow sublimation
        sp_mrro_y(iloop)           = mrro_y
        sp_mrros_y(iloop)          = mrros_y
        sp_mrrob_y(iloop)          = mrrob_y
        sp_lai_y(iloop)            = lai_y
        sp_test_y(iloop,:)         = test_gpp_y
    end subroutine update_spinup_values

    !-----------------------------------------------------------------------
    subroutine init_spinup_variables()
        implicit none
        allocate(sp_gpp_y(nloops))
        allocate(sp_npp_y(nloops))
        allocate(sp_ra_y(nloops))
        allocate(sp_rh_y(nloops)) 
        allocate(sp_wetlandCH4_y(nloops))
        allocate(sp_wetlandCH4prod_y(nloops))
        allocate(sp_wetlandCH4cons_y(nloops))
        ! Carbon Pools  (KgC m-2)
        allocate(sp_cLeaf_y(nloops))
        allocate(sp_cStem_y(nloops))
        allocate(sp_cRoot_y(nloops))
        allocate(sp_cOther_y(nloops))                             
        allocate(sp_cLitter_y(nloops))
        allocate(sp_cLitterCwd_y(nloops))                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        allocate(sp_cSoil_y(nloops))
        allocate(sp_cSoilFast_y(nloops))
        allocate(sp_cSoilSlow_y(nloops))
        allocate(sp_cSoilPassive_y(nloops))                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        allocate(sp_cCH4_y(nloops,nlayers))                                                          ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        allocate(sp_fBNF_y(nloops))
        allocate(sp_fN2O_y(nloops))
        allocate(sp_fNloss_y(nloops))
        allocate(sp_fNnetmin_y(nloops))
        allocate(sp_fNdep_y(nloops))                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        allocate(sp_nLeaf_y(nloops))
        allocate(sp_nStem_y(nloops))
        allocate(sp_nRoot_y(nloops))
        allocate(sp_nOther_y(nloops))
        allocate(sp_nLitter_y(nloops))
        allocate(sp_nLitterCwd_y(nloops))
        allocate(sp_nSoil_y(nloops))
        allocate(sp_nMineral_y(nloops))                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        allocate(sp_hfls_y(nloops))
        allocate(sp_hfss_y(nloops))                                                       ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        allocate(sp_ec_y(nloops))
        allocate(sp_tran_y(nloops))
        allocate(sp_es_y(nloops))                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        allocate(sp_hfsbl_y(nloops))                                                         ! Snow sublimation
        allocate(sp_mrro_y(nloops))
        allocate(sp_mrros_y(nloops))
        allocate(sp_mrrob_y(nloops))
        allocate(sp_lai_y(nloops))
        allocate(sp_test_y(nloops,9))
    end subroutine init_spinup_variables
    
    !-----------------------------------------------------------------------
    subroutine deallo_spinup_variables()
        implicit none
        deallocate(sp_gpp_y)
        deallocate(sp_npp_y)
        deallocate(sp_ra_y)
        deallocate(sp_rh_y) 
        deallocate(sp_wetlandCH4_y)
        deallocate(sp_wetlandCH4prod_y)
        deallocate(sp_wetlandCH4cons_y)
        ! Carbon Pools  (KgC m-2)
        deallocate(sp_cLeaf_y)
        deallocate(sp_cStem_y)
        deallocate(sp_cRoot_y)
        deallocate(sp_cOther_y)                             
        deallocate(sp_cLitter_y)
        deallocate(sp_cLitterCwd_y)                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        deallocate(sp_cSoil_y)
        deallocate(sp_cSoilFast_y)
        deallocate(sp_cSoilSlow_y)
        deallocate(sp_cSoilPassive_y)                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        deallocate(sp_cCH4_y)                                                          ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        deallocate(sp_fBNF_y)
        deallocate(sp_fN2O_y)
        deallocate(sp_fNloss_y)
        deallocate(sp_fNnetmin_y)
        deallocate(sp_fNdep_y)                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        deallocate(sp_nLeaf_y)
        deallocate(sp_nStem_y)
        deallocate(sp_nRoot_y)
        deallocate(sp_nOther_y)
        deallocate(sp_nLitter_y)
        deallocate(sp_nLitterCwd_y)
        deallocate(sp_nSoil_y)
        deallocate(sp_nMineral_y)                    ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        deallocate(sp_hfls_y)
        deallocate(sp_hfss_y)                                                       ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        deallocate(sp_ec_y)
        deallocate(sp_tran_y)
        deallocate(sp_es_y)                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        deallocate(sp_hfsbl_y)                                                         ! Snow sublimation
        deallocate(sp_mrro_y)
        deallocate(sp_mrros_y)
        deallocate(sp_mrrob_y)
        deallocate(sp_lai_y)
        deallocate(sp_test_y)
    end subroutine deallo_spinup_variables

    !-----------------------------------------------------------------------
end module mod_spinup