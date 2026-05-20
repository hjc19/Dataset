module mo_list
    !
    !  verlet list, see more at page 147 of <computer simulation of liquids>
    !
    use mo_syst
    use mo_config
    implicit none

    ! global constants.
    !! skin
    real(8), private, parameter :: set_nlcut = 0.70d0
    !! for sake of memory saving, we consider listmax neighbor of one particle at most
    !! enlarge this if you study 3D system or high volume fraction system
    integer, private, parameter :: listmax = 200

    ! neighbor list of one particle
    type listone_t
        ! neighbor number
        integer    :: nbsum
        ! neighbor index
        integer    :: nblist(listmax)
        ! precalculated relative relation, used for distance calculation in perodical cell
        integer(1) :: iround(free,listmax)
        integer(1) :: cory(listmax)
        !
        real(4)    :: con0(free)
    end type

    ! neighbor list struct of system
    type list_t
        integer                                    :: natom
        type(listone_t), allocatable, dimension(:) :: list
        ! contact number
        integer, allocatable, dimension(:)         :: nbi
        ! tag particle is rattler or not
        integer, allocatable, dimension(:)         :: rattlerflag
        real(8)    :: nlcut
    end type

    type(list_t) :: nb


    ! voronoi cell
    type voro_one_t
        integer :: nbsum
        integer :: nblist(listmax)
        integer :: vid1(listmax), vid2(listmax)
    end type

contains

    subroutine init_list(tnb, tcon)
        !
        !  allocate memory
        !
        implicit none

        ! para list
        type(list_t), intent(inout) :: tnb
        type(con_t),  intent(in)    :: tcon

        ! local
        integer :: tnatom

        tnatom    = tcon%natom
        tnb%natom = tnatom
        tnb%nlcut = set_nlcut

        if ( allocated(tnb%list) ) then
            if ( size(tnb%list) /= tnatom ) then
                deallocate( tnb%list )
                allocate( tnb%list(tnatom) )
            end if
        else
            allocate( tnb%list(tnatom) )
        end if
    end subroutine

    subroutine make_list(tnb, tcon)
        !
        !  make list of system
        !
        implicit none

        ! para list
        type(con_t),  intent(in)    :: tcon
        type(list_t), intent(inout) :: tnb

        ! local
        real(8) :: lainv(free), dra(free), rai(free), raj(free), ri, rj, rij2, dij
        integer :: cory, iround(free)
        integer :: i, j, k, itemp

        associate(                 &
            natom  => tcon%natom,  &
            ra     => tcon%ra,     &
            r      => tcon%r,      &
            la     => tcon%la,     &
            strain => tcon%strain, &
            list   => tnb%list,    &
            nlcut  => tnb%nlcut    &
            )

            nlcut = set_nlcut
            lainv = 1.d0 / la

            ! set nbsum to zero
            list(:)%nbsum = 0

            do i=1, natom

                list(i)%con0 = ra(:,i)
                rai          = ra(:,i)
                ri           = r(i)

                do j=i+1, natom

                    raj = ra(:,j)
                    rj  = r(j)

                    dra = raj - rai

                    cory   = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij  = ri + rj

                    if ( rij2 > ( dij+nlcut )**2 ) cycle

                    if ( list(i)%nbsum < listmax ) then
                        itemp                   = list(i)%nbsum
                        itemp                   = itemp + 1
                        list(i)%nbsum           = itemp
                        list(i)%nblist(itemp)   = j
                        list(i)%iround(:,itemp) = iround
                        list(i)%cory(itemp)     = cory
                    end if

                end do
            end do

        end associate
    end subroutine

    subroutine calc_z( tnb, tcon )
        !
        !  calculate coordination number
        !
        implicit none

        ! para list
        type(con_t),  intent(inout)    :: tcon
        type(list_t), intent(inout) :: tnb

        ! local
        real(8) :: dra(free), rai(free), raj(free), ri, rj, rij2, dij
        integer :: cory, iround(free)
        integer :: i, j, k, jj

        if ( allocated( tnb%nbi ) .and. size(tnb%nbi) /= tcon%natom ) then
            deallocate( tnb%nbi )
        end if

        if ( .not. allocated( tnb%nbi ) ) then
            allocate( tnb%nbi(tcon%natom) )
        end if

        associate(                 &
            natom  => tcon%natom,  &
            ra     => tcon%ra,     &
            r      => tcon%r,      &
            la     => tcon%la,     &
            strain => tcon%strain, &
            list   => tnb%list,    &
            nbi    => tnb%nbi      &
            )

            ! set nbsum to zero
            nbi = 0

            do i=1, natom

                rai          = ra(:,i)
                ri           = r(i)

                do jj=1, list(i)%nbsum

                    j = list(i)%nblist(jj)

                    raj = ra(:,j)
                    rj  = r(j)

                    dra = raj - rai

                    cory = list(i)%cory(jj)
                    iround = list(i)%iround(:,jj)

                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij  = ri + rj

                    if ( rij2 > ( dij )**2 ) cycle

                    nbi(i) = nbi(i) + 1
                    nbi(j) = nbi(j) + 1

                end do
            end do

            tcon%z = tnb%nbi
        end associate
    end subroutine

    function check_list( tnb, tcon ) result(flag)
        !
        !  determine remake list or not
        !
        implicit none

        ! para list
        type(list_t), intent(in) :: tnb
        type(con_t),  intent(in) :: tcon

        ! result
        logical :: flag

        ! local
        real(8) :: maxdis, dra(free), dr2
        integer :: i

        associate(               &
            natom => tcon%natom, &
            ra    => tcon%ra,    &
            nlcut => tnb%nlcut   &
            )

            maxdis = 0.d0
            do i=1, tcon%natom
                dra = tcon%ra(:,i) - tnb%list(i)%con0
                dr2 = sum( dra**2 )
                if ( maxdis < dr2 ) maxdis = dr2
            end do

        flag = .false.
        if ( maxdis > 0.25 * nlcut**2 ) flag = .true.

        end associate
    end function

    ! ToDo
    ! subroutine check_rattler
    function calc_rattler( tcon, tnblist ) result(flag)
        implicit none

        ! para list
        type(con_t),  intent(inout)           :: tcon
        type(list_t), intent(in), optional :: tnblist

        ! result
        integer, dimension(tcon%natom) :: flag

        ! local
        type(list_t) :: lclist
        integer      :: i

        if ( present(tnblist) ) then
            lclist = tnblist
        else
            call init_list( lclist, tcon )
            call make_list( lclist, tcon )
        end if

        call calc_z( lclist, tcon )

        flag = 0
        do i=1, lclist%natom
            if ( lclist%nbi(i) == 0 ) then
                flag(i) = 1
            elseif ( lclist%nbi(i) <= free ) then
                write(*,*) "There exist unstable particle(s)"
                stop
            end if
        end do
    end function

end module
