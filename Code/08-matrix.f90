module mo_mode
    use mo_syst
    use mo_config
    use mo_network
    use mo_force
    implicit none

    type matrix_t
        ! dynamic matrix; inverse matrix
        real(8), allocatable, dimension(:,:)   :: dymatrix, dymatrix0, invmatrix
        ! eigenvalues of matrix; participation ratio
        real(8), allocatable, dimension(:)     :: egdymatrix, pw
        ! vars related to works of Maloney
        real(8), allocatable, dimension(:)     :: varXi_x, varXi_y, varXi_s
        ! Tong Hua's \psi
        real(8), allocatable, dimension(:)     :: psi_th
        ! third order matrix
        real(8), allocatable, dimension(:,:,:) :: trimatrix
        ! ndim = free * m
        ! mdim = ndim + freedom of box
        integer :: natom, ndim, mdim
        ! flags of box's freedom
        logical :: boxflag   = .false.      ! can box length change
        logical :: xyflag    = .false.      ! can box length change seperately
        logical :: shearflag = .false.      ! can box be changed by shear
        logical :: allflag   = .false.      ! can box be changed by shear and compression

    contains
        ! calculate eigen-problem
        procedure :: solve   => solve_mode
        ! calculate inverse matrix
        procedure :: inv     => calc_inverse_matrix
    end type

    type(matrix_t) :: mode, mode0, mode1, mode2

contains
    subroutine init_mode( tmode, tcon, opboxflag, opxyflag, opshearflag )
        implicit none

        ! para list
        type(matrix_t), intent(inout)        :: tmode
        type(con_t),    intent(in)           :: tcon
        logical,        intent(in), optional :: opboxflag, opxyflag, opshearflag

        ! local

        associate(                       &
            natom     => tmode%natom,    &
            mdim      => tmode%mdim,     &
            ndim      => tmode%ndim,     &
            boxflag   => tmode%boxflag,  &
            xyflag    => tmode%xyflag,   &
            shearflag => tmode%shearflag &
            )

            boxflag   = .false.
            xyflag    = .false.
            shearflag = .false.
            if ( present(opboxflag) .and. present(opxyflag) ) then
                stop "boxflag and xyflag can not exist simultaneously"
            end if
            if ( present(opboxflag) )   boxflag   = opboxflag
            if ( present(opxyflag) )    xyflag    = opxyflag
            if ( present(opshearflag) ) shearflag = opshearflag

            natom = tcon%natom
            ndim  = free * tcon%natom
            mdim  = ndim

            if ( boxflag   ) mdim = mdim + 1
            if ( xyflag    ) mdim = mdim + free
            if ( shearflag ) mdim = mdim + 1

            if ( .not. allocated(tmode%dymatrix) ) then
                allocate( tmode%dymatrix(mdim,mdim), tmode%egdymatrix(mdim) )
               !allocate( tmode%dymatrix0(mdim,mdim) )
               !allocate( tmode%trimatrix(mdim,mdim,mdim) )
            end if

        end associate
    end subroutine

    subroutine kernel_matrix_compress_box( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
        implicit none

        ! para list
        integer, intent(in)    :: mdim
        real(8), intent(inout) :: dymatrix(mdim,mdim)
        integer, intent(in)    :: i, j
        integer, intent(in)    :: natom
        real(8), intent(in)    :: vr, vrr
        real(8), intent(in)    :: xij, yij, rij

        ! local
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: rex, rey
        real(8) :: rexx, rexy, reyx, reyy
        real(8) :: rexex, rexey, reyey

        real(8) :: mij(free,free)
        real(8) :: miex(free), miey(free)
        !real(8) :: mee(2,2)
        real(8) :: mee

        ! \partial r_ij / \partial x_ij
        rx = xij / rij
        ry = yij / rij
        ! \partial r_ij / \partial \epsilon_x = rx * x_ij
        rex = rx * xij
        rey = ry * yij

        ! \partial r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]
        rexx = rxx * xij
        rexy = rxy * xij
        reyx = rxy * yij
        reyy = ryy * yij
        !
        rexex = rxx * xij**2
        reyey = ryy * yij**2
        rexey = rxy * xij*yij

        mij(1,1) = - ( vrr*rx**2 + vr*rxx )
        mij(2,2) = - ( vrr*ry**2 + vr*ryy )
        mij(1,2) = - ( vrr*rx*ry + vr*rxy )
        mij(2,1) = mij(1,2)


        miex(1) = - ( vrr*rx*rex + vr*rexx )
        miex(2) = - ( vrr*ry*rex + vr*rexy )
        miey(1) = - ( vrr*rx*rey + vr*reyx )
        miey(2) = - ( vrr*ry*rey + vr*reyy )

        mee = 0.d0
        mee = mee + vrr*rex*rex + vr*rexex
        mee = mee + vrr*rey*rey + vr*reyey
        mee = mee + 2 * ( vrr*rex*rey + vr*rexey )

        ! con ij
          dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
        & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
          dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
        & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij

          dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
          dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

          ! ex ey
          dymatrix( free*(i-1)+1:free*i, free*natom+1 ) = &
        & dymatrix( free*(i-1)+1:free*i, free*natom+1 ) + miex + miey
          dymatrix( free*natom+1, free*(i-1)+1:free*i ) = &
        & dymatrix( free*natom+1, free*(i-1)+1:free*i ) + miex + miey

          dymatrix( free*(j-1)+1:free*j, free*natom+1 ) = &
        & dymatrix( free*(j-1)+1:free*j, free*natom+1 ) - miex - miey
          dymatrix( free*natom+1, free*(j-1)+1:free*j ) = &
        & dymatrix( free*natom+1, free*(j-1)+1:free*j ) - miex - miey
          dymatrix(free*natom+1,free*natom+1) = &
        & dymatrix(free*natom+1,free*natom+1) + mee
    end subroutine

    subroutine keybond_modify_matrix_compress_box( tmode, tcon, i, j, tdelta )
        implicit none

        ! para list
        type(matrix_t),  intent(inout) :: tmode
        type(con_t),     intent(in)    :: tcon
        integer,         intent(in)    :: i, j
        real(8),              optional :: tdelta

        ! local config
        real(8) :: vr, vrr
        real(8) :: xij, yij, rij, wili_ij, press_ij
        real(8) :: dra(free), rij2, dij, delta

        ! local matrix
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: rex, rey
        real(8) :: rexx, rexy, reyx, reyy
        real(8) :: rexex, rexey, reyey

        real(8) :: mij(free,free)
        real(8) :: miex(free), miey(free)
        real(8) :: mee

        associate(                      &
            dymatrix => tmode%dymatrix, &
            radius   => tcon%r,         &
            natom    => tmode%natom,    &
            mdim     => tmode%mdim,     &
            ndim     => tmode%ndim      &
            )
            
            if ( .not. present(tdelta) ) then
                delta = 1.d0
            else
                delta = tdelta
            end if

            dra  = calc_dra( tcon, i, j )
            rij2 = sum(dra**2)

            dij = radius(i) + radius(j)

            if ( rij2 > dij**2 ) then
                print*, "something error in key bond, Please check !!!"
                stop
            end if

            rij  = sqrt(rij2)

            ! \partial V_ij / \partial r_ij
            vr = - delta * ( 1.d0 - rij/dij )**(alpha-1) / dij
            ! \partial^2 V_ij / \partial r_ij^2
            vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2

            xij = dra(1)
            yij = dra(2)
            
            ! \partial r_ij / \partial x_ij
            rx = xij / rij
            ry = yij / rij
            ! \partial r_ij / \partial \epsilon_x = rx * x_ij
            rex = rx * xij
            rey = ry * yij

            ! \partial r_ij / [ \partial x_j \partial x_j ]
            rxx = 1.d0/rij - xij**2 /rij**3
            ryy = 1.d0/rij - yij**2 /rij**3
            rxy =          - xij*yij/rij**3
            ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
            rxx = 1.d0/rij - xij**2 /rij**3
            ryy = 1.d0/rij - yij**2 /rij**3
            rxy =          - xij*yij/rij**3
            ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]
            rexx = rxx * xij
            rexy = rxy * xij
            reyx = rxy * yij
            reyy = ryy * yij
            !
            rexex = rxx * xij**2
            reyey = ryy * yij**2
            rexey = rxy * xij*yij

            mij(1,1) = - ( vrr*rx**2 + vr*rxx )
            mij(2,2) = - ( vrr*ry**2 + vr*ryy )
            mij(1,2) = - ( vrr*rx*ry + vr*rxy )
            mij(2,1) = mij(1,2)

            miex(1) = - ( vrr*rx*rex + vr*rexx )
            miex(2) = - ( vrr*ry*rex + vr*rexy )
            miey(1) = - ( vrr*rx*rey + vr*reyx )
            miey(2) = - ( vrr*ry*rey + vr*reyy )

            mee = 0.d0
            mee = mee + vrr*rex*rex + vr*rexex
            mee = mee + vrr*rey*rey + vr*reyey
            mee = mee + 2 * ( vrr*rex*rey + vr*rexey )

            ! con ij
              dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
            & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) + mij
              dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
            & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) + mij

              dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = & 
            & dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) - mij
              dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = & 
            & dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) - mij
            
            ! ex ey
              dymatrix( free*(i-1)+1:free*i, free*natom+1 ) = &
            & dymatrix( free*(i-1)+1:free*i, free*natom+1 ) - miex - miey
              dymatrix( free*natom+1, free*(i-1)+1:free*i ) = &
            & dymatrix( free*natom+1, free*(i-1)+1:free*i ) - miex - miey

              dymatrix( free*(j-1)+1:free*j, free*natom+1 ) = &
            & dymatrix( free*(j-1)+1:free*j, free*natom+1 ) + miex + miey
              dymatrix( free*natom+1, free*(j-1)+1:free*j ) = &
            & dymatrix( free*natom+1, free*(j-1)+1:free*j ) + miex + miey
              dymatrix(free*natom+1,free*natom+1) = &
            & dymatrix(free*natom+1,free*natom+1) - mee
        end associate
    end subroutine

    subroutine stability_make_dymatrix( tmode, tcon, tnet, opflag, tdelta )
        implicit none

        ! para list
        type(matrix_t),  intent(inout)        :: tmode
        type(con_t),     intent(in)           :: tcon
        type(network_t), intent(in), optional :: tnet
        real(8),                     optional :: tdelta
        integer,         intent(in)           :: opflag

        ! local
        integer :: i, j, ii
        real(8) :: dra(free), rij, rij2, dij
        real(8) :: vr, vrr
        real(8) :: xij, yij
        real(8) :: ks, l0, delta

        associate(                      &
            dymatrix => tmode%dymatrix, &
            radius   => tcon%r,         &
            natom    => tmode%natom,    &
            mdim     => tmode%mdim,     &
            ndim     => tmode%ndim      &
            )

            dymatrix = 0.d0

            if ( .not. present( tdelta ) ) then
                delta = 1.d0
            else
                delta = tdelta
            end if

            if ( .not. present( tnet ) ) then

                do i=1, natom-1
                    do j=i+1, natom

                        dra  = calc_dra( tcon, i, j )
                        rij2 = sum(dra**2)

                        dij = radius(i) + radius(j)

                        if ( rij2 > dij**2 ) cycle

                        rij  = sqrt(rij2)

                        ! \partial V_ij / \partial r_ij
                        vr = - delta * ( 1.d0 - rij/dij )**(alpha-1) / dij
                        ! \partial^2 V_ij / \partial r_ij^2
                        vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2

                        xij = dra(1)
                        yij = dra(2)

                        if ( opflag == 11 ) then
                            call kernel_matrix_compress_box( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        end if

                    end do
                end do

    
            else
                associate(             &
                    list  => tnet%sps, &
                    nlist => tnet%nsps &
                    )

                    do ii=1, nlist
                        
                        !if ( .not. list(ii)%existflag) cycle

                        i  = list(ii)%i
                        j  = list(ii)%j
                        l0 = list(ii)%l0
                        ks = list(ii)%ks

                        dra = calc_dra( tcon, i, j )
                        rij2 = sum(dra**2)

                        rij  = sqrt(rij2)

                        ! \partial V_ij / \partial r_ij
                        vr = - delta * ( 1.d0 - rij/dij )**(alpha-1) / dij
                        !vr = ks * ( rij - l0 )
                        ! \partial^2 V_ij / \partial r_ij^2
                        vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2
                        !vrr = ks

                        xij = dra(1)
                        yij = dra(2)

                        if ( opflag == 11 ) then
                            call kernel_matrix_compress_box( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
                        end if

                    end do
      
                end associate

            end if

        end associate
    end subroutine

    subroutine solve_mode( this, oprange )
        implicit none

        ! para list
        class(matrix_t) :: this
        integer, optional :: oprange

        ! local
        integer :: rangevar

        rangevar = 0
        if ( present( oprange ) ) rangevar = oprange

        this%dymatrix0 = this%dymatrix

        associate(                        &
            mdim       => this%mdim,      &
            dymatrix   => this%dymatrix,  &
            egdymatrix => this%egdymatrix &
            )

            call solve_matrix( dymatrix, mdim, egdymatrix, rangevar )

        end associate
    end subroutine

    subroutine calc_inverse_matrix( this, nu_ratter )
        implicit none

        ! para list
        class(matrix_t), intent(inout) :: this
        integer, intent(in)            :: nu_ratter

        ! local
        integer :: i

        this%invmatrix = this%dymatrix

        ! m0 = Matrix
        ! m1 = eig(m0)
        ! m_inv = \sum_{i_real} eigenvalue(i) * e_{\omega_i} * e_{\omega_i}^T
        do i=1, this%mdim
            if ( i<=free*(1+nu_ratter) ) then
                this%invmatrix(:,i) = 0.d0
            else
                this%invmatrix(:,i) = ( 1.d0 / this%egdymatrix(i) ) * this%invmatrix(:,i)
            end if
        end do

        this%invmatrix = matmul(this%invmatrix, transpose(this%dymatrix))
    end subroutine

    subroutine solve_matrix(a,order,b, rangevar)
        implicit none

        ! para list
        integer :: order
        real(8),dimension(1:order,1:order) :: a
        real(8),dimension(1:order) :: b
        integer :: rangevar

        ! local
        character :: jobz = 'V'
        character :: range = 'A'
        character :: uplo = 'U'

        integer :: n
        !real, dimension  :: a
        integer :: lda
        integer :: vl, vu, il, iu   ! will not referenced when range=A
        real(8) :: abstol           ! important variable
        integer :: m
        !real, dimension :: w        ! use b above ! output eigenvalues
        real(8), allocatable, dimension(:,:) :: z
        integer :: ldz
        integer, allocatable, dimension(:) :: isuppz
        real(8), allocatable, dimension(:) :: work
        integer :: lwork
        integer, allocatable, dimension(:) :: iwork
        integer :: liwork
        integer :: info

        n = order
        lda = order
        ldz= order

        if( rangevar == 0 ) then
            range = 'A'
        elseif( rangevar > 0 .and. rangevar <= order ) then
            range = 'I'
            il = 1
            iu = rangevar
        else
            write(*,*) "error rangevar set"
            stop
        end if

        !--
        abstol = -1
        !abstol = dlamch('S')        ! high precision

        !allocate(a(order,order)); allocate(b(order))
        allocate(z(order,order))
        allocate(isuppz(2*order))


        !- query the optimal workspace
        allocate(work(1000000)); allocate(iwork(1000000))
        lwork  = -1
        liwork = -1

        call dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,b,z,       &
                   & ldz, isuppz, work, lwork, iwork, liwork, info)

        lwork  = int(work(1))
        liwork = int(iwork(1))

        deallocate(work); deallocate(iwork)

        allocate(work(lwork)); allocate(iwork(liwork))

        !vvv- main

        call dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,b,z,       &
                   & ldz, isuppz, work, lwork, iwork, liwork, info)

        if(info .ne. 0) then
            write(*,*) '** Find error in subroutine diago'
            stop
        end if
        !^^^- done

        deallocate(work); deallocate(iwork)

        !- output results via a
        a = z

        !deallocate(a); deallocate(b)
        deallocate(z); deallocate(isuppz)
    end subroutine

end module
