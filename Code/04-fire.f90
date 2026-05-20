module mo_fire
    !
    !  FIRE method of local minimization, see more at paper https://doi.org/10.1103/PhysRevLett.97.170201
    !
    use mo_syst
    use mo_config
    use mo_list
    use mo_network
    use mo_force
    implicit none

    !-- fire var
    real(8), private, parameter :: fmax    = 1.d-14
    integer, private, parameter :: stepmax = 1e7
    integer, private :: step
    !--
    real(8), private, parameter :: dt0   = 3.d-2
    real(8), private, parameter :: dtmax = 3.d-1
    real(8), private, parameter :: beta0 = 1.d-1
    real(8), private, parameter :: finc  = 1.1d0
    real(8), private, parameter :: fdec  = 0.5d0
    real(8), private, parameter :: fbeta = 0.99d0
    integer, private, parameter :: nmin  = 5
    !--
    real(8), private :: dt, dt2, dt22, beta, power, fn, vn
    integer, private :: count


    type(con_t)  :: confire, confire1, confire2
    type(list_t) :: nbfire

    procedure(abstract_force), pointer :: calc_force_h => null()

contains

    subroutine init_confire( tconfire, tcon, tnet )
        implicit none

        ! para list
        type(con_t),     intent(inout)        :: tconfire
        type(con_t),     intent(in)           :: tcon
        type(network_t), intent(in), optional :: tnet

        ! copy con to confire
        tconfire = tcon

        ! if net don't exist, allocate list
        if ( .not. present( tnet ) ) then
            call init_list( nbfire, tconfire )
        end if

        calc_force_h => calc_force
    end subroutine

    subroutine check_system_force( tcon, tnet )
        implicit none

        ! para list
        type(con_t),     intent(inout)           :: tcon
        type(network_t), intent(inout), optional :: tnet

        if ( .not. present( tnet ) ) then
            call make_list( nbfire, tcon )
            call calc_force_h( tcon, nbfire )
        else
           !call make_network( tnet, tcon )
            call calc_force_spring( tcon, tnet )
        end if
    end subroutine

    ! main 2
    subroutine mini_fire_cp( tcon, tnet, opboxp_set )
        implicit none

        ! para list
        type(con_t),     intent(inout)        :: tcon
        type(network_t), intent(in), optional :: tnet
        real(8),         intent(in), optional :: opboxp_set    ! target press

        ! local
        real(8) :: dstrain
        logical :: nonnetwork_flag
        real(8) :: boxp_set
        logical :: cp_flag
        logical :: boxp_flag

        real(8) :: lainv(free)
        real(8) :: dt, beta, temp
        real(8) :: onembeta, betavndfn
        integer :: cumn
        integer :: i, ipin

        ! 1. network or not
        nonnetwork_flag = .true.
        if ( present( tnet ) ) nonnetwork_flag = .false.


        boxp_flag = .false.
        if ( present( opboxp_set ) ) then
            boxp_set  = opboxp_set
            boxp_flag = .true.
        end if

        ! cp or not
        cp_flag = .false.
        if ( boxp_flag ) cp_flag = .true.


        associate(                     &
            natom    => tcon%natom,    &
            ra       => tcon%ra,       &
            va       => tcon%va,       &
            fa       => tcon%fa,       &
            press    => tcon%press,    &
            strain   => tcon%strain,   &
            strainv  => tcon%strainv,  &
            strainf  => tcon%strainf,  &
            stress   => tcon%stress,   &
            la       => tcon%la,       &
            lav      => tcon%lav,      &
            laf      => tcon%laf       &
            )

            lainv = 1.d0 / la

            ! initial sets
            fa      = 0.d0 ; va      = 0.d0
            lav     = 0.d0 ; laf     = 0.d0
            strainv = 0.d0 ; strainf = 0.d0
            dt      = dt0  ; beta    = beta0
            cumn = 0

            ! pre
            if ( nonnetwork_flag ) then
                call make_list( nbfire, tcon )
                call calc_force_h( tcon, nbfire )
            else
                call calc_force_spring( tcon, tnet )
            end if


            if ( cp_flag ) then
                if ( boxp_flag ) laf = press - boxp_set
            end if

            do step=1, stepmax

                ! step length
                dt2  = dt    * 0.5d0
                dt22 = dt**2 * 0.5d0

                ! velocity verlet method / move 1 / config
                ra = ra + va * dt + fa * dt22
                va = va + fa * dt2

                ! velocity verlet method / move 1 / box
                if ( cp_flag ) then
                    !v affine deformation
                    if ( boxp_flag ) then
                        la(1)  = la(1)  + lav(1) * dt + laf(1) * dt22
                        lav(1) = lav(1) + laf(1) * dt2
                        !
                        la(2:free) = la(2:free) * ( la(1)*lainv(1) )
                        ra = ra * ( la(1)*lainv(1) )
                    end if
                    !^
                    lainv = 1.d0 / la
                end if

                ! check list
                if ( nonnetwork_flag .and. check_list( nbfire, tcon ) ) then
                    call make_list( nbfire, tcon )
                end if

                ! velocity verlet method / force
                if ( nonnetwork_flag ) then
                    call calc_force_h( tcon, nbfire )
                else
                    call calc_force_spring( tcon, tnet )
                end if

                if ( cp_flag ) then
                    if ( boxp_flag ) laf = press - boxp_set
                end if

                ! velocity verlet method / move 2
                va = va + fa * dt2
                if ( cp_flag ) lav = lav + laf * dt2

                ! fire
                cumn  = cumn + 1
                power = sum( fa * va ) + sum( laf*lav ) + strainv * strainf

                fn = sqrt( sum( fa**2 ) + sum(laf**2) + strainf**2 )
                vn = sqrt( sum( va**2 ) + sum(lav**2) + strainv**2 )

                onembeta  = 1.d0 - beta
                betavndfn = beta*vn/fn

                va      = onembeta * va      + betavndfn * fa
                lav     = onembeta * lav     + betavndfn * laf
                strainv = onembeta * strainv + betavndfn * strainf

                if ( power > 0.d0 .and. cumn > nmin ) then
                    dt   = min( dt*finc, dtmax )
                    beta = beta * fbeta
                end if

                if ( power < 0.d0 ) then
                    cumn    = 0
                    dt      = dt * fdec
                    beta    = beta0
                    va      = 0.d0
                    lav     = 0.d0
                    strainv = 0.d0
                end if

                temp = maxval( abs(fa) )
                if ( cp_flag ) then
                    if ( boxp_flag ) temp = max( temp, abs( press - boxp_set  ) )
                end if

                if ( temp < fmax ) then
                    !write(*,'(5es16.6)') 1.0, tcon%press, tcon%pressx, tcon%pressy, tcon%stress
                    exit
                end if

                if ( step == stepmax ) then
                    write(*,*) "subroutine reached step maximum and existed with no force balance"
                    stop
                end if

            end do


        end associate
    end subroutine

end module
