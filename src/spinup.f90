module mod_spinup
    use mod_data
    use driver
    contains
    subroutine run_spinup()
        write(*,*)"This is spinup", nloops
        do iloop = 1, nloops
            write(*,*) "iloop: ", iloop
            itest = iloop
            call teco_simu()
            ! write(*,*) sum(all_gpp_d), all_cOther_d(size(all_cOther_d))
            if (sum(all_gpp_d) .eq. 0.) stop
        enddo
    end subroutine run_spinup
end module mod_spinup