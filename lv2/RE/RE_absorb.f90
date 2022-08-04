program main

  use scatter
  use mtmod  
  
  implicit none
  include 'mpif.h'

  integer :: myrank, numprocs
 integer :: ierr
 integer :: wcomm = MPI_COMM_WORLD

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
    real(8),dimension(100) :: q,u,n,n_total
    real(8),dimension(100) :: q_total,u_total,w_total
    real(8),dimension(100) :: q_norm, u_norm
    real(8),dimension(100) :: mu_ave
    real(8) :: chi_pre, chi_post, cos_chi_post
    real(8) :: theta_pre, theta_post, phi_pre, phi_post
    real(8) :: q_ein_pre, q_eex_pre, u_eex_pre
    real(8) :: x,y,z,tau,mu,mu_
    real(8) :: pi,a,b,c,d
    real(8) :: p,o,q_i,u_i
    real(8) :: theta_e, phi_e
    real(8) :: depth,i_,j_
    real(8) :: theta_, phi_
    real(8) :: theta_s, s, mu_dummy
    real(8) :: count,number_
    real(8) :: section,random
    real(8) :: delta_tau, delta_section,max_a,min_a,true_a,tau_slab
    real(8) :: x_p,y_p,z_p,s_p
    real(8) :: z_delta,tent_section
    real(8) :: alpha,alpha_0,alpha_a
    real(8) :: x_delta,y_delta ! shuoud be fixed
    real(8) :: radius,radius_xy
    real(8) :: w,tau_a,section_total,tau_ratio,all_n,all_w

    integer :: rank_number
    integer :: access

 !  call system_clock(exec_t0, t_rate, count_max)
  

    !bug? deback option check !FFLAGS = -O2 -g  -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow

    write(fileName, '("RE_absorb_tau=20.csv")')
    ! write(fileName, '("lv1_x_layer.csv")')                                                                     
    fileName = trim(adjustl(fileName))

    fileNum  = 1

    open(fileNum, file=fileName, status="replace")

    write(fileNum, '("mu, q, u, number, w_total")')

    t_max = 100000000
    pi = 4.0d0*atan(1.0d0)
    alpha_0 = 4.33*10.0d0**(-6.0d0)!scatter conficient !outer
    alpha_a = 7.89d0*10.0d0**(-11.0d0)!absorb conficient(T=10^9)
   ! alpha_0 = 1.0
   ! depth = 34464290.0d0 !outer
     !depth = 10.0d0**7.0d0
    ! tent_section = depth/100
    ! tau_slab = sqrt(pi/2)*alpha_0*depth*erf(1.0/sqrt(2.0))
     tau_slab = 20.0d0
     depth = tau_slab/(sqrt(pi/2)*alpha_0*erf(1.0/sqrt(2.0)))
     tent_section = depth/100 
     tau_ratio = alpha_a/alpha_0 ! tau_a/tau_s
    ! tau_ratio = 0.0d00

     
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(wcomm, myrank, ierr)
    call MPI_COMM_SIZE(wcomm, numprocs, ierr)
    
    rank_number = numprocs
   ! t_max = 2 ! test to be fixed 
    iter_max = 100000000/rank_number !should be fixed
   ! iter_max = 100

    tau_slab = sqrt(pi/2)*alpha_0*depth*erf(1.0/sqrt(2.0))
       

    if(myrank==0)then
       print'("tau="f11.5",depth="f20.5",N_per_core="i10",tau_a/tau_s="f11.5)',tau_slab,depth,iter_max,tau_ratio
    end if
 
    call sgrnd(myrank+1)
    
    do i = 1,100
          q(i) = 0.0d0
          u(i) = 0.0d0
          q_total(i) = 0.0d0
          u_total(i) = 0.0d0
          mu_total(i) = 0.0d0
          num(i) = 0.0d0
          n_total(i) = 0.0d0
       end do


       do iter = 1, iter_max

          w = 1.0d0
     
       x_total = 0.0d0
       y_total = 0.0d0
       z_total = (grnd()-0.5d0)*depth/0.5d0
   
       s = grnd()
       theta_s = 2.0d0*pi*s

       stokes_pre(1) = cos(theta_s)
       stokes_pre(2) = sin(theta_s) 

      ! stokes_pre(1) = 0.0d0
      ! stokes_pre(2) = 0.0d0 !nonpolarization
       stokes_pre(3) = 0.0d0 ! Initial condition
        mu_k = 1.0d0


       do t = 0, t_max
                 
          if (t .eq. 0) then

             c = grnd()
             d = grnd()
         
             mu_= c
             !mu_ = 2.0d0*c-1.0d0
             theta_ = acos(mu_)
             phi_= 2.0d0*pi*d

             k_f_pre(1) = 0.001d0
             k_f_pre(2) = k_f_pre(1)*sin(theta_)*cos(phi_)
             k_f_pre(3) = k_f_pre(1)*sin(theta_)*sin(phi_)
             k_f_pre(4) = k_f_pre(1)*cos(theta_) ! Initial conditions

            ! k_f_pre(2) = 0.01d0*k_f_pre(1)
            ! k_f_pre(3) = 0.0d0
            ! k_f_pre(4) = sqrt(k_f_pre(1)**2.0d0-k_f_pre(2)**2.0d0)     
          
             k_f_post(1) = k_f_pre(1)
             k_f_post(2) = k_f_pre(2)
             k_f_post(3) = k_f_pre(3)
             k_f_post(4) = k_f_pre(4)

             stokes_post(1) = stokes_pre(1)
             stokes_post(2) = stokes_pre(2)         
             
          else

             q_eex_pre = stokes_pre(1)
             u_eex_pre = stokes_pre(2)

       
             k_e_pre(2) = k_f_pre(2)
             k_e_pre(3) = k_f_pre(3)
             k_e_pre(4) = k_f_pre(4)
             k_e_pre(1) = sqrt(k_e_pre(2)*k_e_pre(2)+k_e_pre(3)*k_e_pre(3)+k_e_pre(4)*k_e_pre(4)) !cleaning
             k0_e_pre = k_e_pre(1)
             

             
             do

                call sample_scattered_photons(k0_e_pre, k_e_pre, k_e_post, mu_k,theta_pre, theta_post, phi_pre, phi_post)

                if(phi_pre/=phi_pre)then
                   print '("chi=NaN")'
                end if
                

      
                
                call calc_chi(k_e_pre, k_e_post, mu_k, &
                     chi_pre, chi_post, theta_pre, theta_post, phi_pre, phi_post, cos_chi_post)

                call rotate_q_chi(q_eex_pre, u_eex_pre, q_ein_pre, chi_pre, chi_post)
                                           

               ! if (k0_e_pre .lt. 1.0d-4) then  ! Thomson Scattering Limit
                   call calc_thomson_diff(k0_e_pre, q_ein_pre, mu_k, &
                        k0_e_post_dummy, cross_section_norm) !3/11 test3
               ! else  ! Compton Scattering
                !   call calc_kn_diff(k0_e_pre, q_ein_pre, mu_k, &
                 !       k0_e_post_dummy, cross_section_norm)
               ! end if

                r = grnd()


                if (cross_section_norm .gt. r) exit !rejection


             end do
             

             call update_stokes(stokes_pre, k_e_pre, k_e_post, mu_k, chi_pre, chi_post, stokes_post)

          !   call lorentz_inverse(k_f_post, k_e_post, beta, gamma, n) ! lorentz method


             k_f_post(1) = k_e_post(1)
             k_f_post(2) = k_e_post(2)
             k_f_post(3) = k_e_post(3)
             k_f_post(4) = k_e_post(4)
             
          end if
          
          section = 0.0d0

          do 
              random = grnd()

              alpha = alpha_0*exp(-abs(z_total)**2.0d0/(2.0d0*depth**2.0d0))

             tau = tent_section*alpha

              if(exp(-tau) .lt. random)exit
             
              section = section + tent_section

              tau_a = tent_section*alpha*tau_ratio! to be fixed

              w = w*exp(-tau_a)

                z = k_f_post(4)*tent_section/k_f_post(1) ! shuoud be fixed                                                                                        
               z_total = z_total + z ! shuoud be fixed 

               if(abs(z_total) .gt. depth)exit
          end do

        ! delta_tau = log(1/(1-random))
         delta_tau = log(1/random)
          
         delta_section = delta_tau/alpha

         s_p = section + delta_section

         tau_a = delta_tau*tau_ratio ! to be fixed

         w = w*exp(-tau_a)

 !        print'(f100.90,i5)',w,myrank

         if (w .lt. 10.0d0**(-3.0d0))exit
          
          x = k_f_post(2)*s_p/k_f_post(1)! shuoud be fixed
          y = k_f_post(3)*s_p/k_f_post(1)! shuoud be fixed 
          z_delta = k_f_post(4)*delta_section/k_f_post(1)! shuoud be fixed
          

       z_total = z_total + z_delta ! shuoud be fixed
       y_total = y_total + y
       x_total = x_total + x ! shuoud be fixed

       mu = k_f_post(4)/sqrt(k_f_post(2)**2.0d0 + k_f_post(3)**2.0d0 + k_f_post(4)**2.0d0)

       stokes_pre(1) = stokes_post(1)/sqrt(stokes_post(1)*stokes_post(1)+stokes_post(2)*stokes_post(2))
       stokes_pre(2) = stokes_post(2)/sqrt(stokes_post(1)*stokes_post(1)+stokes_post(2)*stokes_post(2)) ! cleaning of stokes parameters

       if(stokes_post(1) .eq. 0.0d0 .and. stokes_post(2) .eq. 0.0d0)then
          stokes_pre(1) = 0.0d0
          stokes_pre(2) = 0.0d0
       end if
       

     !  stokes_pre(1) = stokes_post(1)
       !  stokes_pre(2) = sqrt(1-stokes_pre(1)*stokes_pre(1))  !another cleaning

       k_f_pre(1) = k_f_post(1)
       k_f_pre(2) = k_f_post(2)
       k_f_pre(3) = k_f_post(3)
       k_f_pre(4) = k_f_post(4)

       if(abs(z_total) .gt. depth) exit !Escape Conditions
     
    end do
    
    
    if (abs(z_total).ge.depth)then
       
       do i = 1,100
          i_ = real(i)
          if(abs(mu) .ge. (i_- 1.0d0)*0.01d0 .and. abs(mu) .lt. i_*0.01d0)then
          q(i) = q(i) + stokes_post(1)*w
          u(i) = u(i) + stokes_post(2)*w
          !    mu_total(i) = mu_total(i) + abs(mu)
          n(i) = n(i) + 1.0d0
          num(i) = num(i) + w
       end if
       
    end do
 end if
 
   
   ! count = real(iter)
    if(myrank==0)then
    do number = 1,10
     !  number_ = real(number)
       if (iter .eq. number*(iter_max/10))then
          print '(i9,i5)',iter,myrank
       end if
    end do
 end if
 
    
 end do

! print'(f100.90,i5)',w,myrank

 call MPI_REDUCE(q, q_total, 100, MPI_DOUBLE_PRECISION,MPI_SUM,0, wcomm, ierr)
 call MPI_REDUCE(u, u_total, 100, MPI_DOUBLE_PRECISION,MPI_SUM,0, wcomm, ierr)
 call MPI_REDUCE(num, w_total, 100, MPI_DOUBLE_PRECISION,MPI_SUM,0, wcomm, ierr)
 call MPI_REDUCE(n, n_total, 100, MPI_DOUBLE_PRECISION,MPI_SUM,0, wcomm, ierr)
 
 if (myrank ==0) then
    all_n = 0.0d0
    all_w = 0.0d0
 do j = 1,100
    q_norm(j) = q_total(j)/n_total(j)
    u_norm(j) = u_total(j)/n_total(j)
    i_ = real(j)
    write(fileNum, '(f5.2,",",f20.10,",",f20.10,",",f15.0,",",f20.10)')&
         i_*0.01d0 ,q_norm(j),u_norm(j),n_total(j),w_total(j)

    all_n = all_n + n_total(j)
    all_w = all_w + w_total(j)
 end do
end if
 
 call MPI_FINALIZE(ierr)
 
        
 close(fileNum)

! call system_clock(exec_t1)
 if(myrank==0)then
    print'("tau="f11.5",N="f11.0",w="f15.5)',tau_slab,all_n,all_w
 end if
 

    
end program main
       

          

    

    
