program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_fire
    use mo_mode
    use mo_mechanics
    use mo_math
    
    implicit none

    type(con_t)   :: con_trim, kbcon
    type(matrix_t) :: exmode, kbmode, mode_tmp
    type(list_t) :: list
   
    character(80) ::  f0, f3, f4, f5, t1, t2, t3
    logical :: file_exist
   
    integer        :: count_rattler, nonzero_mode_index, i_index, j_index, m, n, neg
    integer, allocatable, dimension(:) :: neg_index, pos_index, neg_index_i, pos_index_i, neg_index_j, pos_index_j
    integer, allocatable, dimension(:) :: rattler_flag

    real(8), allocatable, dimension(:) :: negfre, posfre
    real(8) :: press, deltap

    !------ vars
    call testvar
    
    allocate( rattler_flag(sets%natom) )
    write(t1,*), sets%natom
    write(t2,*), sets%seed
    write(t3,'(es26.1)'), sets%press
!-------------------- generate initial random config -----------------!   
    f0 = "con_"//trim(adjustl(t1))//"_"//trim(adjustl(t2))//"_"//trim(adjustl(t3))//".txt"
    inquire( file=f0, exist=file_exist)
    if ( file_exist ) then
        call read_config_simple_version( con, f0, sets%natom )
   
        call init_list(list, con)
        call make_list(list, con)
        call calc_z(list, con)

        call calc_force( con, list )
        
        print*,maxval(abs(con%fa))
      
        if ( maxval(abs(con%fa)) > 1e-12 ) then
            print*, "The force is unstable !"
            stop
        end if
    else
        print*, "there no config exist"
        !stop
        call init_system( con, sets%natom, sets%phi )
        call gen_rand_config( con, sets%seed )
    
        call init_confire( confire, con )
        call mini_fire_cp( confire, opboxp_set=sets%press )
        con = confire

        call init_list(list, con)
        call make_list(list, con)
        call calc_z(list, con)

        call calc_force( con, list )
        
    end if
!-------------------------------------------------------------------------!
    
    count_rattler = calc_rattler_number(con)
    nonzero_mode_index = free * ( 1 + count_rattler ) + 1
    
    con_trim = con
    call save_config_simple_version( con_trim, f0 )

    !--- 2D extended Hessian Matrix
    call init_mode( exmode, con, opboxflag=.true. )
    call stability_make_dymatrix( exmode, con, opflag=11, tdelta=sets%beta )

    mode_tmp = exmode
!----------------------------------------------------------------------------!

!--------------------------- Network Analysis -------------------------------!
    call init_network( net, con )
    call make_particle_network( net, con )

    !--- for check 
    if ( net%nsps*2 /= sum(list%nbi(:)) ) then
        print*, "There are some unknown error during finding k-bond !"
    end if

    m=0
    n=0

    allocate(negfre(net%nsps),posfre(net%nsps))
    allocate(neg_index(net%nsps),pos_index(net%nsps))
    allocate(neg_index_i(net%nsps),pos_index_i(net%nsps),neg_index_j(net%nsps),pos_index_j(net%nsps))
    
    do k = 1, net%nsps
        kbcon = con
        kbmode = mode_tmp
       
        i_index = net%sps(k)%i
        j_index = net%sps(k)%j
        
        call keybond_modify_matrix_compress_box(kbmode, kbcon, i_index, j_index, tdelta=sets%beta)
        
        call solve_mode( kbmode, 1 )
        if ( kbmode%egdymatrix(1)< - 1.d-12 ) then
            m=m+1
        
            neg_index(m) = k
            neg_index_i(m) = net%sps(k)%i
            neg_index_j(m) = net%sps(k)%j
            negfre(m) = kbmode%egdymatrix(1)
            
        else
            n=n+1
            
            pos_index(n) = k
            pos_index_i(n) = net%sps(k)%i
            pos_index_j(n) = net%sps(k)%j
            posfre(n) = kbmode%egdymatrix(nonzero_mode_index)
        end if
        
    end do

    if ( m==0 ) goto 888
    f3 = "packing_kb_"//trim(adjustl(t1))//"_"//trim(adjustl(t2))//"_"//trim(adjustl(t3))//".txt"
    open(unit=91,file=f3,Access='Append')
        do i = 1, m
            write(91,*) neg_index(i), neg_index_i(i), neg_index_j(i), negfre(i)
        end do
    close(91)
 
888 if ( n==0 ) goto 999
    f4 = "packing_nonkb_"//trim(adjustl(t1))//"_"//trim(adjustl(t2))//"_"//trim(adjustl(t3))//".txt"
    open(unit=92, file=f4, Access='Append')
        do i = 1, n
            write(92,*) pos_index(i), pos_index_i(i), pos_index_j(i), posfre(i)
        end do
    close(92)

999 f5 = "packing_kb_ratio_"//trim(adjustl(t1))//"_"//trim(adjustl(t2))//"_"//trim(adjustl(t3))//".txt"
    open(unit=93,file=f5,Access='Append')
        write(93,*) m, n, net%nsps, dble(m)/dble(net%nsps)
    close(93)
            
print*, sets%seed, dble(m)/dble(net%nsps)
!===========================================================================================!


!===========================================================================================!
contains

    subroutine testvar
        implicit none

        !sets%natom = 256 
        !read(*,*) sets%num
        !sets%num = 1
        sets%phi = 0.86d0
        !sets%press = 1.d-2
        !sets%seed = 10
        sets%beta = 1.d0
        read(*,*) sets%seed
        read(*,*) sets%natom
        read(*,*) sets%press
        !read(*,*) sets%beta
    end subroutine testvar

end program main
