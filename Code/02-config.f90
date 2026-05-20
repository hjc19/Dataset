module mo_config
    !
    !  base struct and method of configuration generating, storing...
    !
    use mo_syst
    use mo_math, only: init_rand
    implicit none

    type con_t
        ! Atom
        !! number of particles
        integer :: natom
        !! configuration, velocity, force
        real(8), allocatable, dimension(:,:) :: ra, va, fa
        !! radius
        real(8), allocatable, dimension(:)   :: r
        !! contact number
        real(8), allocatable, dimension(:)   :: z 

        ! box
        !! volume fraction
        real(8) :: phi
        !! lengthes of box, la = (lx, ly, [lz])
        real(8) :: la(free)
        !! v for vary lengthes of box
        real(8) :: lav(free)
        real(8) :: laf(free)
        !! strain, for shear
        real(8) :: strain
        real(8) :: strainv, strainf

        ! sets
        real(8) :: T
        ! properties
        !! energy
        real(8) :: Ea, Ek, Ev
        !! stress tensor
        real(8) :: stress, press 
        !!! pressures of x, y and z
        real(8) :: pressxyz(free)
    contains
        ! con%dra(i,j) => (xij, yij, zij)
        procedure :: dra         => calc_dra
        ! con%vec(i, [xk,yk,zk]) => round( xj-xk, yj-yk, zj-zk )
        procedure :: vec         => calc_vec
        ! con%calc_phi() => phi
        procedure :: calc_phi    => calc_phi
        ! compress and shear
    end type

    type(con_t), target :: con, con0, contemp, contemp1, contemp2

contains

    subroutine init_system( tcon, tnatom, tphi )
        !
        ! initiate system, allocate memory of system
        !
        implicit none

        ! para list
        type(con_t),       intent(inout) :: tcon
        integer,           intent(in)    :: tnatom
        real(8), optional, intent(in)    :: tphi

        tcon%natom = tnatom
        if ( present( tphi ) ) tcon%phi = tphi

        allocate( tcon%ra(free,tnatom), tcon%r(tnatom), &
            & tcon%va(free,tnatom), tcon%fa(free,tnatom) )
        
        if (free==3) then
            print*, "Warning !!! stress just calculated xz-direction !"
        end if
    end subroutine

    subroutine gen_rand_config( tcon, tseed, tphi )
        !
        !  generate random configuration with given seed and phi
        !
        implicit none

        ! para list
        type(con_t), intent(inout)        :: tcon
        integer,     intent(inout)        :: tseed
        real(8),     intent(in), optional :: tphi

        ! local
        integer :: i, j

        ! initialized rand
        !call init_rand(tseed)
        !tseed = 0

        if ( present( tphi ) ) tcon%phi = tphi

        associate(                &
            natom  => tcon%natom, &
            ra     => tcon%ra,    &
            fa     => tcon%fa,    &
            va     => tcon%va,    &
            r      => tcon%r,     &
            la     => tcon%la,    &
            strain => tcon%strain &
            )

            ! radius
            r(1:natom/2)       = 0.5d0
            r(natom/2+1:natom) = 0.5d0 * ratio

            ! box length
            la     = calc_box_length(tcon)
            strain = 0.d0

            ! config
            call srand(tseed)
            if ( free == 2 ) then
                do i = 1, natom
                    ra(1, i) = rand(0)
                    ra(2, i) = rand(0)
                end do
            else if ( free == 3 ) then
                do i = 1, natom
                    ra(1, i) = rand(0)
                    ra(2, i) = rand(0)
                    ra(3, i) = rand(0)
                end do
            end if

            do i=1, natom
                do j=1, free
                    ra(j,i) = ( ra(j,i) - 0.5d0 ) * la(j)
                end do
            end do

            ! f v
            va = 0.d0
            fa = 0.d0

        end associate
    end subroutine

    subroutine read_config_simple_version( tcon, tfilename, tnatom, tphi )
        !
        !  read configuration from text file
        !
        implicit none

        ! para list
        type(con_t),  intent(inout)        :: tcon
        character(*), intent(in)           :: tfilename
        integer,      intent(in)           :: tnatom
        real(8),      intent(in), optional :: tphi

        ! local
        integer :: i
        real(8) :: a, b

        ! allocate array of tcon
        if ( present( tphi ) ) then
            call init_system( tcon, tnatom, tphi )
        else
            call init_system( tcon, tnatom )
        end if

        ! read config
        open(901,file=tfilename)
            read(901, *) tcon%phi, tcon%press, tcon%stress
            read(901, *) tcon%la, tcon%strain
            do i=1, tnatom
                read(901,*) tcon%ra(:,i), tcon%r(i) 
            end do
        close(901)

        !!!!
        !tcon%r = tcon%r * 0.5

        ! phi
        tcon%phi = pi * sum(tcon%r**2) / product(tcon%la)
    end subroutine

    pure function calc_box_length(tcon) result(l)
        implicit none

        ! para list
        type(con_t), intent(in) :: tcon

        ! result
        real(8) :: l

        ! local
        real(8) :: phi
        real(8) :: sdisk, volume

        phi = tcon%phi

        ! V_n(r) = pi^(n/2) / Gamma( n/2 + 1 ) * r^n
        sdisk = sqrt(pi**free) / gamma(dble(free)/2.d0+1) * sum(tcon%r**free)

        ! box length
        volume = sdisk / phi
        l      = volume ** ( 1.d0 / dble(free) )
    end function


    pure function calc_dra( this, ti, tj ) result(dra)
        !
        !  dra = rj - ri = rij = [xij, yij, zij]
        !
        implicit none

        ! para list
        class(con_t), intent(in) :: this
        integer,      intent(in) :: ti
        integer,      intent(in) :: tj

        ! result
        real(8), dimension(free) :: dra

        ! local
        real(8) :: rai(free), raj(free)
        integer :: k, cory, iround(free)

        associate(                &
            ra     => this%ra,    &
            la     => this%la,    &
            strain => this%strain &
            )

            rai = ra(:,ti)
            raj = ra(:,tj)

            dra = raj - rai

            cory = nint( dra(free) / la(free) )
            dra(1) = dra(1) - strain * la(free) * cory

            do k=1, free-1
                iround(k) = nint( dra(k) / la(k) )
            end do
            iround(free) = cory

            do k=1, free
                dra(k) = dra(k) - iround(k) * la(k)
            end do

        end associate
    end function

    pure function calc_vec( this, ti, traj ) result(dra)
        !
        !  calculate vector distance of particle i and point[x,y,z]
        !
        implicit none

        ! para list
        class(con_t), intent(in) :: this
        integer,      intent(in) :: ti
        real(8), dimension(free), intent(in) :: traj

        ! result
        real(8), dimension(free) :: dra

        ! local
        real(8) :: rai(free)
        integer :: k, cory, iround(free)

        associate(                &
            ra     => this%ra,    &
            la     => this%la,    &
            strain => this%strain &
            )

            rai = ra(:,ti)

            dra = rai - traj

            cory = nint( dra(free) / la(free) )
            dra(1) = dra(1) - strain * la(free) * cory

            do k=1, free-1
                iround(k) = nint( dra(k) / la(k) )
            end do
            iround(free) = cory

            do k=1, free
                dra(k) = dra(k) - iround(k) * la(k)
            end do

        end associate
    end function

    subroutine trim_config( tcon, opsumxyz )
        !
        !  make sure all the partiles locate in the periodical cell
        !
        implicit none

        ! para list
        type(con_t), intent(inout)        :: tcon
        logical,     intent(in), optional :: opsumxyz  ! set center of mass to zero

        ! local
        real(8) :: lainv(free), temp
        integer :: iround(free), cory
        integer :: i, k

        associate(                &
            natom  => tcon%natom, &
            ra     => tcon%ra,    &
            la     => tcon%la,    &
            strain => tcon%strain &
            )

            lainv = 1.d0 / la

            do i=1, natom

                cory = nint( ra(free,i) * lainv(free) )
                ra(1,i) = ra(1,i) - strain * la(free) * cory

                do k=1, free-1
                    iround(k) = nint( ra(k,i) * lainv(k) )
                end do
                iround(free) = cory

                do k=1, free
                    ra(k,i) = ra(k,i) - iround(k) * la(k)
                end do

            end do

            if ( present( opsumxyz ) ) then
                if ( opsumxyz .eqv. .true. ) then
                    do i=1, free
                        temp = sum(ra(i,:)) / natom
                        ra(i,:) = ra(i,:) - temp
                    end do
                end if
            end if

        end associate
    end subroutine

    pure function calc_phi(this) result(re)
        !
        !  calculate volume fraction
        !  phi = volume of particles / volume of box
        !
        implicit none

        ! para list
        class(con_t), intent(in) :: this

        ! result
        real(8) :: re

        ! local
        real(8) :: sdisk, volume

        volume = product( this%la(1:free) )

        !if ( free == 2 ) then
        !    sdisk = pi * sum( this%r(1:this%natom)**2 )
        !elseif ( free == 3 ) then
        !    sdisk = 4.d0/3.d0 * pi * sum( this%r(1:this%natom)**3 )
        !end if
        sdisk = sqrt(pi**free) / gamma(dble(free)/2.d0+1) * sum(this%r**free)

        re = sdisk / volume
    end function

end module

subroutine save_config_simple_version( tcon, tfilename )
    use mo_config
    implicit none

    ! para list
    type(con_t),  intent(in) :: tcon
    character(*), intent(in) :: tfilename

    ! local
    integer :: i

    associate(                  &
        natom   => tcon%natom,  &
        ra      => tcon%ra,     &
        r       => tcon%r,      &
        la      => tcon%la,     &
        strain  => tcon%strain  &
        )

        open(901,file=tfilename)
            write(901,'(3es26.16)') tcon%phi, tcon%press, tcon%stress
            write(901,'(3es26.16)') la, strain
            do i=1, natom
                write(901,*) ra(:,i), r(i)
            end do
        close(901)
    end associate
end subroutine


