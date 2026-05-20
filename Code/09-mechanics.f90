module mo_mechanics
    use mo_syst
    use mo_config
    use mo_list
    use mo_fire
    use mo_mode
    use mo_math
    use mo_network
    implicit none


contains
    subroutine calc_ParticleNet_Bond_Wili_stress( tnetwork, tcon, opcase, tdelta )
        implicit none

        ! para list
        type(con_t),     intent(in)    :: tcon
        type(network_t), intent(inout) :: tnetwork
        integer,         optional      :: opcase
        real(8),         optional      :: tdelta

        ! local
        integer :: ii, i, j, k, casenu
        real(8) :: lnow
        real(8) :: ks, l0, delta
        real(8), dimension(free) :: rai, raj, dra
        real(8), dimension(tnetwork%nsps) :: temp
        real(8) :: ri, rj, rij2, rij, dij, a 
        real(8) :: lainv(free)
        integer :: iround(free), cory

        associate(                    &
            nsps   => tnetwork%nsps,  &
            sps    => tnetwork%sps,   &
            natom    => tcon%natom,   &
            radius   => tcon%r,       &
            ra       => tcon%ra,      &
            la       => tcon%la,      &
            strain   => tcon%strain   &
            )

            if ( present(opcase) ) then
                casenu = opcase
            else
                casenu = 0
            end if

            if ( .not. present( tdelta ) ) then
                delta = 1.d0
            else
                delta = tdelta
            end if
            
            lainv = 1.d0 / la
            
            do ii=1, nsps

                i = sps(ii)%i
                j = sps(ii)%j
                l0 = sps(ii)%l0

                rai = ra(:,i)
                ri  = radius(i)

                raj = ra(:,j)
                rj  = radius(j)

                dra = raj - rai

                cory = nint( dra(free) * lainv(free) )
                dra(1) = dra(1) - strain * la(free) * cory

                do k=1, free-1
                    iround(k) = nint( dra(k) * lainv(k) )
                end do
                iround(free) = cory

                do k=1, free
                    dra(k) = dra(k) - iround(k) * la(k)
                end do

                rij2 = sum( dra**2 )
                dij = ri + rj

                rij = sqrt( rij2 )

                sps(ii)%Es = ( 1.d0 - rij/dij )**alpha/alpha
                sps(ii)%wili = delta * (1.d0 - rij/dij)**(alpha-1) * rij / dij
                sps(ii)%press = sps(ii)%wili / product(la) / free
                sps(ii)%fij  = sps(ii)%wili / rij
                sps(ii)%stress = - dra(1) * dra(2) * sps(ii)%wili / rij2 / product(la)
                
                temp(ii) = (rij-l0)/l0
           
            end do
                
            select case( casenu )
            ! 0. normal case
            case(0)
                a = 0.d0 
            ! 1. for compression
            case(1)
                sps(:)%cuij = temp
            ! 2. for shear 
            case(2)
                sps(:)%suij = temp
            end select
        end associate
    end subroutine
 
    function calc_rattler_number( tcon )  result( count_rattler )
        implicit none
        !--- para list
        type(con_t) :: tcon
        integer :: count_rattler

        !--- local para list
        type(matrix_t) :: dymode
        integer, allocatable, dimension(:) :: rattler_flag
        integer :: i

        !--- initial config rattler
        rattler_flag = calc_rattler( tcon )
        count_rattler = 0
        do i = 1, tcon%natom
            if( rattler_flag(i) == 1 ) then
                count_rattler = count_rattler + 1
            end if
        end do
    end function calc_rattler_number 


end module
