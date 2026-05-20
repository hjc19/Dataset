module mo_force
    use mo_syst
    use mo_config
    use mo_list
    use mo_network
    implicit none

    abstract interface
        subroutine abstract_force( tcon, tnb )
            import :: con_t, list_t
            type(con_t),  intent(inout)  :: tcon
            type(list_t), intent(in), optional :: tnb
        end subroutine
        logical function abstract_fun_force( tcon, tnb )
            import :: con_t, list_t
            type(con_t),  intent(inout)  :: tcon
            type(list_t), intent(in), optional :: tnb
        end function
    end interface

contains

    subroutine calc_force( tcon, tnb )
        implicit none

        ! para list
        type(con_t),  intent(inout) :: tcon
        type(list_t), intent(in), optional :: tnb

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilixyz(free)
        real(8) :: lainv(free)
        integer :: iround(free), cory
        integer :: i, j, k, jj

        associate(                    &
            natom    => tcon%natom,   &
            radius   => tcon%r,       &
            ra       => tcon%ra,      &
            fa       => tcon%fa,      &
            Ea       => tcon%Ea,      &
            la       => tcon%la,      &
            strain   => tcon%strain,  &
            stress   => tcon%stress,  &
            press    => tcon%press,   &
            pressxyz => tcon%pressxyz &
            )

        Ea     = 0.d0
        fa     = 0.d0
        stress = 0.d0
        wili   = 0.d0; wilixyz = 0.d0

        if ( present(tnb) ) then

            associate(list => tnb%list)

            do i=1, natom

                rai = ra(:,i)
                ri  = radius(i)

                do jj=1, list(i)%nbsum

                    j = list(i)%nblist(jj)
                    iround = list(i)%iround(:,jj)
                    cory = list(i)%cory(jj)

                    raj = ra(:,j)
                    rj  = radius(j)

                    dra = raj - rai
                    dra(1) = dra(1) - cory * strain * la(free)

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = ri + rj

                    if ( rij2 > dij**2 ) cycle

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    wilixyz = wilixyz + fr * dra(:)**2
                    
                    ! 3d - just xz shear stress                    
                    stress = stress - dra(1) * dra(free) * fr

                end do

            end do

            end associate

        else

            lainv = 1.d0 / la

            do i=1, natom

                rai = ra(:,i)
                ri  = radius(i)

                do j=i+1, natom

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

                    if ( rij2 > dij**2 ) cycle

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    wilixyz = wilixyz + fr * dra(:)**2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    ! 3d - just xz shear stress                    
                    stress = stress - dra(1) * dra(free) * fr

                end do

            end do

        end if

        stress   = stress  / product(la)
        press    = wili    / product(la) / free
        pressxyz = wilixyz / product(la)

        end associate
    end subroutine

    subroutine calc_force_spring( tcon, tnet )
        implicit none

        ! para list
        type(con_t),     intent(inout)          :: tcon
        type(network_t), intent(in), optional   :: tnet

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: rij2, rij, l0, fr, wij, wili, wilixyz(free), ks
        integer :: iround(free), cory
        integer :: i, j, k, ii

        associate(                     &
            natom    => tcon%natom,    &
            radius   => tcon%r,        &
            ra       => tcon%ra,       &
            fa       => tcon%fa,       &
            Ea       => tcon%Ea,       &
            la       => tcon%la,       &
            strain   => tcon%strain,   &
            stress   => tcon%stress,   &
            press    => tcon%press,    &
            pressxyz => tcon%pressxyz, &
            list     => tnet%sps,      &
            nlist    => tnet%nsps      &
            )

            Ea     = 0.d0
            fa     = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilixyz = 0.d0

            do ii=1, nlist

                if ( list(ii)%existflag .eqv. .false. ) cycle
                       
                i      = list(ii)%i
                j      = list(ii)%j
                cory   = list(ii)%cory
                iround = list(ii)%iround
                l0     = list(ii)%l0
                ks     = list(ii)%ks

                rai = ra(:,i)
                raj = ra(:,j)

                dra = raj - rai
                dra(1) = dra(1) - cory * strain * la(free)

                do k=1, free
                    dra(k) = dra(k) - iround(k) * la(k)
                end do

                rij2 = sum( dra**2 )
                rij  = sqrt( rij2 )

                Ea = Ea + 0.5d0 * ks * ( l0 - rij )**2

                wij  = ks * ( l0 - rij ) * rij
                wili = wili + wij

                fr = wij / rij2

                wilixyz = wilixyz + fr * dra(:)**2

                fa(:,j) = fa(:,j) + fr * dra
                fa(:,i) = fa(:,i) - fr * dra

                stress = stress - dra(1) * dra(free) * fr  ! 3d ? just xz shear stress 

            end do

            stress   = stress  / product(la)
            press    = wili    / product(la) / free
            pressxyz = wilixyz / product(la)
        end associate
    end subroutine

end module
