module mo_syst
    implicit none

    type set_t
        integer :: natom, num
        real(8) :: phi
        integer :: seed
        real(8) :: press
        !--- force constant for matrix
        real(8) :: beta
    end type
    type(set_t) :: sets

    real(8), parameter :: pi    = 3.1415926535897932d0
    real(8), parameter :: ratio = 1.4d0
    real(8), parameter :: alpha = 2d0
    integer, parameter :: free  = 2

end module mo_syst


module mo_var
    implicit none

    integer :: i, j, k, ii, jj, kk
    integer :: step, step0, step1, step2, nstep
    real(8) :: temp, temp0, temp1, temp2, temp3, temp4, temp5
    integer :: itemp, itemp0, itemp1, itemp2, itemp3, itemp4
    real(8), allocatable, dimension(:)   :: tempvec
    real(8), allocatable, dimension(:,:) :: temparray
    complex(16), allocatable, dimension(:) :: tempcpxvec

    ! main
    real(8) :: testp
    real(8) :: Mk_x, Mk_y, mu_xy, mu_yx, MB, MG_s, MG_xy
    real(8) :: epl_xx, epl_yy

    logical :: exist_flag, eof_flag

    character(250) :: filename, chtemp
    character(250) :: faverage, fresults, fdump, ftemp
end module
