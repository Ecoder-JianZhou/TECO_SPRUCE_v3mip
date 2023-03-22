program TECO
    use mod_data
    use mod_spinup
    use mod_mcmc
    use driver
    use mod_ncd_io
    ! to run TECO simulation, spin-up and data simulation
    implicit none
    ! character(len=32) :: cmdChar_folder, out_csv, out_nc, out_nc_daily, out_nc_hourly, out_nc_monthly
    character(len=1000) :: outDir_nc, outfile 
    
    call CreateFolder(adjustl(trim(outdir)))
    outDir_nc     = adjustl(trim(outdir))//"/outputs_nc"
    outDir_csv    = adjustl(trim(outdir))//"/outputs_csv"
    outDir_h      = adjustl(trim(outDir_nc))//"/Hourly"
    outDir_d      = adjustl(trim(outDir_nc))//"/Daily"
    outDir_m      = adjustl(trim(outDir_nc))//"/Monthly"
    outFile_restart = adjustl(trim(outdir))//"/restart.nc"

    call CreateFolder(adjustl(trim(outDir_nc)))
    call CreateFolder(adjustl(trim(outDir_h)))
    call CreateFolder(adjustl(trim(outDir_d)))
    call CreateFolder(adjustl(trim(outDir_m)))
    call CreateFolder(adjustl(trim(outDir_csv)))

    if (do_spinup)then
        outDir_sp  = adjustl(trim(outdir))//"/res_spinup_nc"
        outFile_sp = adjustl(trim(outDir_sp))//"/results_spinup.nc"
        call CreateFolder(adjustl(trim(outDir_sp)))
        call init_spinup_variables()
    endif

    ! write(outfile,*) adjustl(trim(outDir_csv)),"/Simu_dailyflux_23.csv"
    ! outfile       = trim(outfile)
    ! outfile       = adjustl(outfile)
    ! open(662,file=outfile)
    ! ! write(662,*)'seqday,year,doy, GPP_d, NEE_d, Reco_d, NPP_d, Ra_d, QC1, QC2, QC3, QC4, QC5, QC6, QC7, QC8, Rh_d' 
    ! write(662,*)'seqday,year,doy, GPP_d, NEE_d, Reco_d, NPP_d, Ra_d, QC1, QC2, QC3, QC4, QC5, QC6, QC7, QC8, Rh_d, ta_d, omega_d, &
    !             & b1,b2,b3,mat_QC1,mat_QC2,mat_QC3,mat_QC4,mat_QC5,mat_QC6,mat_QC7,mat_QC8,mat_Rh'

    call get_params()                           ! read parameters values
    call get_forcingdata()                      ! read forcing data
    nHours  = nforcing
    nDays   = nHours/24.
    nYears  = int(nforcing/(365*24))
    nMonths = nYears*12 
    Write(*,*)nHours, nDays, nMonths, nYears
    call assign_all_results(nHours, nDays, nMonths, nYears)
    if (.not. do_snow) call get_snowdepth()
    
    call initialize()                           ! initializations
    if (do_restart)then
        call read_restart(in_restart)
        call initialize_with_restart()
    endif
    itest = 0
    if (do_spinup) then
        call run_spinup()
        call write_spinup_res()
        call deallo_spinup_variables()
    endif 
    if (do_mcmc) then
        call read_obs()
        call run_mcmc()
    endif
    
    itest = 1
    call teco_simu()
    ! writing output file...
!     i_record = 1
!     do idayOfnyear=1,nday4out
!         ! write(662,6602)i,record_yr(i_record),idayOfnyear,(Simu_dailyflux14(j,i_record),j=1,14)
!         write(662,6602)i_record,record_yr(i_record),idayOfnyear,(Simu_dailyflux14_2023(j,i_record),j=1,28)
!         i_record=i_record+1
!     enddo
    
! 6602     format(3(i7,","),27(f15.4,","),(f15.4))
!     close(662)
    call spruce_mip_cmip6Format()
    call write_restart()
    ! end of the simulation, then deallocate the forcing_data
    deallocate(forcing%year)
    deallocate(forcing%doy)
    deallocate(forcing%hour)
    deallocate(forcing%Tair)
    deallocate(forcing%Tsoil)
    deallocate(forcing%RH)
    deallocate(forcing%VPD)
    deallocate(forcing%Rain)
    deallocate(forcing%WS)
    deallocate(forcing%PAR)
    deallocate(forcing%CO2)
    deallocate(forcing%PBOT)
    deallocate(forcing%Ndep)
    ! deallocate the snow_in
    if (.not. do_snow) deallocate(snow_in)
    call deallocate_all_results()
end program TECO

subroutine CreateFolder(path_new)
    implicit none
    character(len=*), INTENT(in) :: path_new
    character (len=:), allocatable :: cmdChar
    logical :: dirExists
    ! ----------------------------------------------------
    allocate(character(len=6+len(path_new)) :: cmdChar)
    cmdChar = "mkdir "//path_new
    inquire( file=trim(path_new)//'/.', exist=dirExists )  ! Works with gfortran, but not ifort
    ! inquire( directory=newDirPath, exist=dirExists )         ! Works with ifort, but not gfortran
    if (.not. dirExists) call system(cmdChar)
    deallocate(cmdChar)
end subroutine CreateFolder