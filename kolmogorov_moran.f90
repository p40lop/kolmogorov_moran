module lb_consts
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    !!! Costanti del modello

    real(dp), parameter :: pi = 3.1415926535897932_dp
    real(dp), parameter :: cs2 = 1.0_dp/3.0_dp, &
                           cs4 = cs2*cs2, &
                           invcs2 = 3.0_dp, &
                           invcs4 = invcs2*invcs2, &
                           h_invcs2 = 0.5_dp*invcs2, &
                           h_invcs4 = 0.5_dp*invcs4
    real(dp), parameter, dimension(0:1):: gravity = (/0.0_dp, 9.8067_dp/)
    ! real(dp),parameter, dimension(0:1):: gravity = (/0.0_dp, 0.0_dp/)
    real(dp), parameter, dimension(0:8) :: weights = (/4.0_dp/9.0_dp, 1.0_dp/9.0_dp, &
                                                       1.0_dp/9.0_dp, 1.0_dp/9.0_dp, &
                                                       1.0_dp/9.0_dp, 1.0_dp/36.0_dp, &
                                                       1.0_dp/36.0_dp, 1.0_dp/36.0_dp, &
                                                       1.0_dp/36.0_dp/)
    integer, parameter, dimension(0:8) :: bb = (/0, 3, 4, 1, 2, 7, 8, 5, 6/)
    real(dp), parameter, dimension(0:8, 0:1) :: &
        ! guardare in colonna : prima colonna è c0
        velocities = reshape((/ &
                             !c0      c1      c2      c3      c4      c5      c6      c7      c8
                             +0.0_dp, +1.0_dp, +0.0_dp, -1.0_dp, +0.0_dp, +1.0_dp, -1.0_dp, -1.0_dp, +1.0_dp, &
                             +0.0_dp, +0.0_dp, +1.0_dp, +0.0_dp, -1.0_dp, +1.0_dp, +1.0_dp, -1.0_dp, -1.0_dp/), &
                             (/9, 2/))

contains

end module lb_consts

program kolmogorov_moran
    use lb_consts
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, parameter        :: xsize = 4, &
                                 ysize = 16
    integer                   :: t = 0, &
                                 save_interval = 2e0, &
                                 tmax, ix, iy, id
    real(dp)                  :: rho_0 = 1.0_dp, &
                                 tau, &
                                 dampening = 1e-7_dp, &
                                 invtau, &
                                 h_invtau
    !  umax = 1e-2, &
    real(dp)                  :: r, &
                                 nu, &
                                 energy = 0.0_dp, &
                                 energy_old = -1.0_dp
    real(dp), allocatable     :: u(:, :, :), forcing(:, :, :)
    real(dp), allocatable     :: f(:, :, :), feq(:, :, :), source_term(:, :, :)
    real(dp), allocatable     :: rho(:, :), moran(:, :)
    character(:), allocatable :: dirpath, filename, command

    if (get_os_type() == 3) then
        command = "powershell New-Item -ItemType Directory -Force -Path '"
    else
        command = "mkdir -p '"
    end if
    command = trim(adjustl(command))

    tau = 0.50009_dp

    tau_loop: do while (tau > 0.5_dp)

        allocate (u(0:xsize - 1, 0:ysize - 1, 0:1))
        allocate (forcing(0:xsize - 1, 0:ysize - 1, 0:1))
        allocate (f(0:xsize - 1, 0:ysize - 1, 0:8))
        allocate (feq(0:xsize - 1, 0:ysize - 1, 0:8))
        allocate (source_term(0:xsize - 1, 0:ysize - 1, 0:8))
        allocate (rho(0:xsize - 1, 0:ysize - 1))
        allocate (moran(0:xsize - 1, 0:ysize - 1))

        t = 0
        energy_old = -1.0_dp
        energy = 0.0_dp
        invtau = 1.0_dp/tau
        h_invtau = 0.5_dp/tau
        nu = cs2*(tau - 0.5_dp)
        ! dampening=nu*umax*(2*pi/ysize)**2
        ! tmax = nint(ysize**2/nu)*xsize*ysize
        save_interval = 2e0
        tmax = 1e5 - 1
        ! call random_number(rho)
        rho = rho_0 !+ rho
        moran = 0.0_dp
        u = 0.0_dp
        ! call random_number(u)
        ! u = u * 1e-4
        feq = equilibrium(u, rho)
        f = feq

        forcing = 0.0_dp
        r = 2*pi/ysize
        do iy = 0, ysize - 1
            forcing(:, iy, 0) = dampening*sin(iy*r)
        end do
        ! open (unit=30, file="forcing.dat")
        ! do ix = 0, xsize - 1
        !     do iy = 0, ysize - 1
        !         write (30, "(2(I0.3,X),SP,*(es15.7e3,x))") ix, iy, forcing(ix, iy, 0)
        !     end do
        ! end do
        ! close (30)
        ! stop
        energy = sum(rho*(u(:, :, 0)**2 + u(:, :, 1)**2))

        allocate (character(1000)::dirpath)
        allocate (character(6)::filename)
        write (dirpath, "(a,f7.5,a)") "data/tau_", tau, "/"
        dirpath = trim(adjustl(dirpath))

        !!! WINDOWS
        call execute_command_line(command//dirpath//"'")
        call execute_command_line(command//dirpath//"/time/'")
        !!! UNIX
        ! call execute_command_line("mkdir -p '"//dirpath//"/time/'")

        open (unit=20, file=dirpath//'macro.dat', access="append")
        moran = local_moran(u(:, :, 0))
        print *, sum(moran)
        call write_data(t)

        print *, sum(rho)
        time_loop: do while (t < tmax)

            energy_old = energy
            t = t + 1

            feq = equilibrium(u, rho)
            ! print"(SP,*(es15.7e3,x))", source(u,forcing)
            ! stop
            f = f + collision(f, feq) + source(u, forcing)
            f = streaming(f)
            rho = sum(f, 3)
            do iy = 0, ysize - 1
                do ix = 0, xsize - 1
                    u(ix, iy, 0) = (f(ix, iy, 1) + f(ix, iy, 5) + f(ix, iy, 8) - (f(ix, iy, 3) + f(ix, iy, 6) + f(ix, iy, 7)) + (forcing(ix, iy, 0)*0.5_dp))/rho(ix, iy)
                    u(ix, iy, 1) = (f(ix, iy, 2) + f(ix, iy, 5) + f(ix, iy, 6) - (f(ix, iy, 4) + f(ix, iy, 7) + f(ix, iy, 8)) + (forcing(ix, iy, 1)*0.5_dp))/rho(ix, iy)
                end do
            end do
            if (mod(t, save_interval) == 0) then
                moran = local_moran(u(:, :, 0))
                write (20, "(I0.6,x,SP,*(es15.7e3,x))") t, sum(moran)
                call write_data(t)
            end if
            if (t == save_interval*5) save_interval = save_interval*10
        end do time_loop
        print *, sum(rho)

        moran = local_moran(u(:, :, 0))
        print *, sum(moran)
        write (20, "(I0.6,x,SP,*(es15.7e3,x))") t, sum(moran)
        call write_data(t)

        deallocate (filename)
        deallocate (dirpath)
        deallocate (u)
        deallocate (forcing)
        deallocate (f)
        deallocate (feq)
        deallocate (source_term)
        deallocate (rho)
        deallocate (moran)

        ! if (tau > 0.55_dp) then
        !     tau = tau - 0.05_dp
        ! else
        !     tau = tau - 0.01_dp
        ! end if
        tau = tau - 0.00001_dp

        close (20)
    end do tau_loop

contains

    subroutine write_data(time)
        integer, intent(in) :: time
        integer :: iix, iiy

        write (filename, "(I0.6)") time
        filename = trim(adjustl(filename))

        open (unit=10, file=dirpath//"/time/"//filename//".dat")
        write (10, "(a,x,a,x,a,x,a,x,a,x,a)") "#1:ix", "2:iy", "3:rho", "4:ux", "5:uy", "6:moran"

        do iix = 0, xsize - 1
            do iiy = 0, ysize - 1
                write (10, "(2(I0.3,X),SP,*(es15.7e3,x))") iix, iiy, rho(iix, iiy), u(iix, iiy, 0), u(iix, iiy, 1), moran(iix, iiy), forcing(iix, iiy, 0), forcing(iix, iiy, 1)
            end do
        end do

        close (10)

    end subroutine write_data

    function equilibrium(global_velocity, density) result(equilibrium_distribution)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        use lb_consts
        implicit none
        real(dp), intent(in), allocatable ::global_velocity(:, :, :)
        real(dp), intent(in), allocatable :: density(:, :)
        real(dp), allocatable :: scalar(:, :), global_squared(:, :)
        real(dp), allocatable :: equilibrium_distribution(:, :, :)
        integer :: iid, iix, iiy

        allocate (scalar(0:xsize - 1, 0:ysize - 1))
        allocate (global_squared(0:xsize - 1, 0:ysize - 1))
        allocate (equilibrium_distribution(0:xsize - 1, 0:ysize - 1, 0:8))

        scalar = 0.0_dp
        global_squared = 0.0_dp
        equilibrium_distribution = 0._dp
    !!$OMP PARALLEl DO
        do iiy = 0, ysize - 1
            do iix = 0, xsize - 1
                global_squared(iix, iiy) = dot_product(global_velocity(iix, iiy, :), global_velocity(iix, iiy, :))
            end do
        end do

        equilibrium_distribution(:, :, 0) = weights(0)*density*(1.0_dp - global_squared*h_invcs2)

        do iid = 1, 8
            scalar = velocities(iid, 0)*global_velocity(:, :, 0) + velocities(iid, 1)*global_velocity(:, :, 1)
            equilibrium_distribution(:, :, iid) = weights(iid)*density*( &
                                                  1.0_dp + (scalar*invcs2) + (scalar*scalar*(h_invcs4)) - ( &
                                                  global_squared*(h_invcs2)))
        end do

    end function equilibrium

    function collision(distribution, equilibrium_distribution) result(collided_distribution)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        use lb_consts
        implicit none
        real(dp), allocatable :: distribution(:, :, :), equilibrium_distribution(:, :, :), collided_distribution(:, :, :)

        allocate (collided_distribution(0:xsize - 1, 0:ysize - 1, 0:8))
        collided_distribution = -invtau*(distribution - equilibrium_distribution)

    end function collision

    function streaming(distribution) result(streamed_distribution)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        use lb_consts
        implicit none
        real(dp), intent(in), allocatable :: distribution(:, :, :)
        real(dp), allocatable :: streamed_distribution(:, :, :)
        integer :: iid, iix, iiy, xnew, ynew
        allocate (streamed_distribution(0:xsize - 1, 0:ysize - 1, 0:8))
        streamed_distribution = 0.0_dp
        streamed_distribution(:, :, 0) = distribution(:, :, 0)

        do iid = 1, 8
            do iiy = 0, ysize - 1
                do iix = 0, xsize - 1
                    xnew = modulo(iix + nint(velocities(iid, 0)), xsize)
                    ynew = modulo(iiy + nint(velocities(iid, 1)), ysize)
                    streamed_distribution(xnew, ynew, iid) = distribution(iix, iiy, iid)
                end do
            end do
        end do

    end function streaming

    function source(global_velocity, forcing_term) result(source_res)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        use lb_consts
        implicit none
        real(dp), intent(in), allocatable :: global_velocity(:, :, :)
        real(dp), intent(in), allocatable :: forcing_term(:, :, :)
        real(dp), allocatable :: source_res(:, :, :)
        integer :: iid
        allocate (source_res(0:xsize - 1, 0:ysize - 1, 0:8))
        source_res = 0._dp
        do iid = 0, 8
            source_res(:, :, iid) = source_res(:, :, iid) + forcing_term(:, :, 0)*(invcs2*velocities(iid, 0) + invcs4*(velocities(iid, 0)*(velocities(iid, 0)*global_velocity(:, :, 0) + velocities(iid, 1)*global_velocity(:, :, 1)) - cs2*global_velocity(:, :, 0)))

            source_res(:, :, iid) = source_res(:, :, iid) + forcing_term(:, :, 1)*(invcs2*velocities(iid, 1) + invcs4*(velocities(iid, 1)*(velocities(iid, 0)*global_velocity(:, :, 0) + velocities(iid, 1)*global_velocity(:, :, 1)) - cs2*global_velocity(:, :, 1)))
            source_res(:, :, iid) = source_res(:, :, iid)*weights(iid)
        end do
        source_res = source_res*(1.0_dp - h_invtau)
    end function source

    function local_moran(scalar) result(local_moran_index)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        use lb_consts
        implicit none
        real(dp), intent(in) :: scalar(0:xsize - 1, 0:ysize - 1)
        real(dp) :: local_moran_index(0:xsize - 1, 0:ysize - 1)
        integer :: iix, iiy, iid, xnew, ynew
        real(dp) :: mean, variance

        mean = sum(scalar)/size(scalar)
        variance = sum((scalar - mean)**2)/(size(scalar) - 1)
        local_moran_index = 0.0_dp

        do iiy = 0, ysize - 1
            do iix = 0, xsize - 1
                do iid = 1, 8
                    xnew = modulo(iix + nint(velocities(id, 0)), xsize)
                    ynew = modulo(iiy + nint(velocities(id, 1)), ysize)
                    local_moran_index(iix, iiy) = local_moran_index(iix, iiy) + weights(id)*(scalar(xnew, ynew) - mean)
                end do

                local_moran_index(iix, iiy) = local_moran_index(iix, iiy)*(scalar(iix, iiy) - mean)/variance

            end do
        end do

    end function local_moran

    integer function get_os_type() result(result)
        !!
        !! Returns one of OS_UNKNOWN, OS_LINUX, OS_MACOS, OS_WINDOWS, OS_CYGWIN,
        !! OS_SOLARIS, OS_FREEBSD, OS_OPENBSD.
        !!
        !! At first, the environment variable `OS` is checked, which is usually
        !! found on Windows. Then, `OSTYPE` is read in and compared with common
        !! names. If this fails too, check the existence of files that can be
        !! found on specific system types only.
        !!
        !! Returns OS_UNKNOWN if the operating system cannot be determined.
        character(len=255) :: val
        integer            :: length, rc
        logical            :: file_exists
        logical, save      :: first_run = .true.
        integer, save      :: ret = 0
        !$omp threadprivate(ret, first_run)

        if (.not. first_run) then
            result = ret
            return
        end if

        first_run = .false.
        result = 0

        ! Check environment variable `OSTYPE`.
        call get_environment_variable('OSTYPE', val, length, rc)

        if (rc == 0 .and. length > 0) then
            ! Linux
            if (index(val, 'linux') > 0) then
                result = 1
                ret = result
                return
            end if

            ! macOS
            if (index(val, 'darwin') > 0) then
                result = 2
                ret = result
                return
            end if

            ! Windows, MSYS, MinGW, Git Bash
            if (index(val, 'win') > 0 .or. index(val, 'msys') > 0) then
                result = 3
                ret = result
                return
            end if

            ! Cygwin
            if (index(val, 'cygwin') > 0) then
                result = 4
                ret = result
                return
            end if

            ! Solaris, OpenIndiana, ...
            if (index(val, 'SunOS') > 0 .or. index(val, 'solaris') > 0) then
                result = 5
                ret = result
                return
            end if

            ! FreeBSD
            if (index(val, 'FreeBSD') > 0 .or. index(val, 'freebsd') > 0) then
                result = 6
                ret = result
                return
            end if

            ! OpenBSD
            if (index(val, 'OpenBSD') > 0 .or. index(val, 'openbsd') > 0) then
                result = 7
                ret = result
                return
            end if
        end if

        ! Check environment variable `OS`.
        call get_environment_variable('OS', val, length, rc)

        if (rc == 0 .and. length > 0 .and. index(val, 'Windows_NT') > 0) then
            result = 3
            ret = result
            return
        end if

        ! Linux
        inquire (file='/etc/os-release', exist=file_exists)

        if (file_exists) then
            result = 1
            ret = result
            return
        end if

        ! macOS
        inquire (file='/usr/bin/sw_vers', exist=file_exists)

        if (file_exists) then
            result = 2
            ret = result
            return
        end if

        ! FreeBSD
        inquire (file='/bin/freebsd-version', exist=file_exists)

        if (file_exists) then
            result = 6
            ret = result
            return
        end if
    end function get_os_type

end program kolmogorov_moran
