module mo_network
    use mo_syst
    use mo_config
    implicit none

    type spring_t
        integer :: i, j, sort
        integer :: cory, iround(free)
        real(8) :: l0, lvec(free)
        real(8) :: vdl
        real(8) :: ks
        real(8) :: Bi, Gi, Gixy, Gis
        real(8) :: Es, wili, press, stress, fij, Bij, Gij, cuij, suij
        logical :: existflag
    end type

    type network_t
        integer :: nsps, max_of_springs
        integer :: natom, nnode, nrattler, ndangle, nridge
        real(8) :: mb, mg, mgs, mgxy
        real(8) :: bcorr, gcorr, kscorr, ksmeancorr
        real(8), allocatable, dimension(:,:) :: kvec
        type(spring_t), allocatable, dimension(:) :: sps
    contains
        procedure :: dra    => calc_spring_dra
        procedure :: len    => calc_spring_len
    end type

    type(network_t) :: net

contains

    subroutine init_network( tnetwork, tcon )
        implicit none

        ! para list
        type(con_t),     intent(in)    :: tcon
        type(network_t), intent(inout) :: tnetwork

        associate(                                    &
            natom          => tnetwork%natom,         &
            max_of_springs => tnetwork%max_of_springs &
            )

            natom = tcon%natom
            max_of_springs = natom * 15
            allocate( tnetwork%sps(max_of_springs) )
            allocate( tnetwork%kvec(free,natom) )

        end associate
    end subroutine

    subroutine make_particle_network( tnetwork, tcon )
        use mo_math, only: dot
        implicit none

        ! para list
        type(network_t), intent(inout) :: tnetwork
        type(con_t),     intent(in)    :: tcon

        ! local
        integer :: i, j, k
        real(8) :: lainv(free), dra(free), rij2, dij
        integer :: cory, iround(free)

        associate(                         &
            natom    => tcon%natom,        &
            ra       => tcon%ra,           &
            r        => tcon%r,            &
            la       => tcon%la,           &
            strain   => tcon%strain,       &
            nsps     => tnetwork%nsps,     &
            sps      => tnetwork%sps      &
            )

            lainv = 1.d0 / la

            nsps = 0
            do i=1, natom
                do j=i+1, natom

                    dra = ra(:,j) - ra(:,i)

                    cory = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - cory * strain * la(free)

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = r(i) + r(j)

                    if ( rij2 > dij**2 ) cycle
                    nsps             = nsps + 1
                    sps(nsps)%i      = i
                    sps(nsps)%j      = j
                    sps(nsps)%cory   = cory
                    sps(nsps)%iround = iround
                    sps(nsps)%lvec   = dra
                    sps(nsps)%l0     = sqrt(rij2)
                    sps(nsps)%ks     = ( alpha-1) * ( 1.d0 - sqrt(rij2)/dij )**(alpha-2) / dij**2
                    sps(nsps)%Bi     = 0.d0
                    sps(nsps)%Gi     = 0.d0
                    sps(nsps)%Gis    = 0.d0
                    sps(nsps)%Gixy   = 0.d0
                    sps(nsps)%existflag = .true.

                end do
            end do

        end associate
      
        ! reallocate array of sps
        tnetwork%sps = tnetwork%sps(1:tnetwork%nsps)
    end subroutine
    

    pure function calc_spring_dra( tnetwork, tcon, tibond ) result(tdra)
        implicit none

        ! para list
        class(network_t), intent(in) :: tnetwork
        type(con_t), intent(in)      :: tcon
        integer, intent(in)          :: tibond

        ! result
        real(8), dimension(free)     :: tdra

        ! local
        integer :: i, j, k
        integer :: cory, iround(free)

        associate(                   &
            natom  => tcon%natom,    &
            ra     => tcon%ra,       &
            r      => tcon%r,        &
            la     => tcon%la,       &
            strain => tcon%strain,   &
            nsps   => tnetwork%nsps, &
            sps    => tnetwork%sps   &
            )

            i = sps(tibond)%i
            j = sps(tibond)%j

            tdra = ra(:,j) - ra(:,i)

            cory = nint( tdra(free) / la(free) )
            tdra(1) = tdra(1) - cory * strain * la(free)

            do k=1, free-1
                iround(k) = nint( tdra(k) / la(k) )
            end do
            iround(free) = cory

            do k=1, free
                tdra(k) = tdra(k) - iround(k) * la(k)
            end do

        end associate
    end function

    pure function calc_spring_len( tnetwork, tcon, tibond ) result(tl)
        implicit none

        ! para list
        class(network_t), intent(in) :: tnetwork
        type(con_t), intent(in)      :: tcon
        integer, intent(in)          :: tibond

        ! result
        real(8) :: dra(free), tl

        dra = calc_spring_dra( tnetwork, tcon, tibond )
        tl  = norm2(dra)
    end function



    function countb(this) result(n)
        implicit none

        !--- para list
        class(network_t) :: this

        !--- result
        integer :: n

        n = count(this%sps(:)%existflag)
    end function

end module
