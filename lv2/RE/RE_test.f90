program main

  use scatter
  use mtmod  
  
  implicit none
!  include 'mpif.h'

 ! integer :: myrank, numprocs
 integer :: ierr
! integer :: wcomm = MPI_COMM_WORLD

 integer :: fileNum,iter,iter_max
 character(50) :: fileName

! integer(8) :: exec_t0, exec_t1, t_rate, count_max, t_diff  ! for exec time

    real(8) :: z_total,x_total,y_total
    integer :: t,t_max, i,j,number, m
    real(8) :: mu_k, q_in_post, k0_e_post_dummy, cross_section_norm, r
    real(8), dimension(4) :: k_e_post, k_f_post
    real(8), dimension(4) :: k_e_pre, k_f_pre
    real(8) :: k0_e_pre
    real(8), dimension(3) :: stokes_pre
    real(8), dimension(3) :: stokes_post
    real(8),dimension(100) :: num,mu_total
    real(8),dimension(200) :: q,u,n,n_total
    real(8),dimension(100) :: q_total,u_total,w_total
    real(8),dimension(100) :: q_norm, u_norm
    real(8),dimension(100) :: mu_ave
    real(8) :: delta_tau, delta_section,max_a,min_a,true_a,tau_slab
    real(8) :: x_p,y_p,z_p,s_p
    real(8) :: z_delta,tent_section
    real(8) :: alpha,alpha_0,alpha_a
    real(8) :: x_delta,y_delta ! shuoud be fixed
    real(8) :: pm,p
    real(8) :: i_,depth

    integer :: rank_number
    integer :: access

 !  call system_clock(exec_t0, t_rate, count_max)
  

    !bug? deback option check !FFLAGS = -O2 -g  -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow

    write(fileName, '("bremsstrahlung_start.csv")')
    ! write(fileName, '("lv1_x_layer.csv")')                                                                     
    fileName = trim(adjustl(fileName))

    fileNum  = 1

    open(fileNum, file=fileName, status="replace")

    write(fileNum, '("z, N")')
   ! t_max = 2 ! test to be fixed 
    iter_max = 100000000 !should be fixed
    ! iter_max = 100

    depth = 10.0d0

!    tau_slab = sqrt(pi/2)*alpha_0*depth*erf(1.0/sqrt(2.0))

    do i = 1,200
       n(i) = 0.0
    end do
    
       
    do iter = 1, iter_max
       
       x_total = 0.0d0
       y_total = 0.0d0

       do
          
          z_total = (2.0d0*grnd()-1.0d0)*depth

          if(exp(-(z_total/depth)**2.0d0) .ge. grnd())exit
       end do
       



       do i=1,200
          i_=real(i)

          if(z_total.ge.(i_-101.0d0)*depth/100.0d0 .and. z_total.lt.(i_-100.0d0)*depth/100)then
             n(i) = n(i) + 1.0d0
          end if
       end do
    end do
    
          
    do i = 1,200
       i_ =real(i)
    write(fileNum, '(f5.2,",",f20.10)')&
         (i_-100.5d0)*depth/100.0d0 ,n(i)
 end do
 
        
 close(fileNum)

    
end program main
       

          

    

    
