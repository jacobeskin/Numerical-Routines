module tridiag

  ! This module contains subroutines for calculating the bound states of 
  ! quarkonia withthe Cornell potential. Three numerical schemes were written.
  ! The successive over relaxation does not work for this problem since 
  ! all the convergence criterion (especially strickt diagonal dominance) are
  ! not met. The Thomas algorithm doesn't seem to be able to differentiate 
  ! properly between the convergent and divergent solutions. The only one that 
  ! does properly its job is the Numerov method with shooting to a fitting
  ! point from both sides of the domain. I have nevertheless kept this module 
  ! intact for possible future use of the SOR and Thomas algorithm in better
  ! suited problems. 

  ! The name "tridiag" comes from this problem being tridiagonal and so these
  ! routines are written that in mind.

  ! Since this module was designated as a source for possible future use for 
  ! for different problems, the subroutines can be optimized quite a bit for 
  ! speed. 

contains

  !----------------------------------------------------------------------------
  !--------------------------- SOR routine ------------------------------------
  !----------------------------------------------------------------------------

  subroutine qSOR(m_c, a, b, l, E_0, n, u)
    implicit none

    integer, intent(in) :: n, l                    ! # of elements and l
    double precision, intent(in) :: m_c, a, b, E_0 ! parameters for model
    double precision, intent(inout) :: u(n)  ! arrays

    double precision, allocatable :: u_0(:)
    double precision :: uu, m2_h, j, ai, x, ai_x, eps, b_i, w, d_x, pi
    integer :: i, k, m

    ! Allocate arrays
    allocate(u_0(n))

    ! Define parameter values
    m2_h = 51.3637898*m_c ! 1 divided by hbar/2m
    eps = 1d-5      ! SOR accuracy limit
    w = 0.9  ! omega for SOR
    d_x = 0.0005*0.0005  ! step size to the second power
    pi = 4*atan(1.0)

    ! Build first guess vector from the asymptotic behavior of wavefunction
    ! Near zero u~x**(l+1) and near infinity u~Ai, and for this Airy function 
    ! will be used its asymptotic aproximation. Search for the closest point
    ! these functions will intersect and build our first trial function u_0.
    do i=2,n-1
       j = (i-1)*0.0005
       x = j**(l+1)
       ai = (1.0/(2*sqrt(pi)))*(j**(-1.0/4))*exp((-2.0/3)*(j**(3.0/2)))
       if (i==2) then
          ai_x = abs(ai-x)
       else if (i>2) then
          if (abs(ai-x)<ai_x) then
             ai_x = abs(ai-x)
             k = i                   ! k is the point where they intersect
          end if
       end if
    end do

    print *,
    print *, "k is:", k

    ! Now build the actual trial function
    u_0(1) = 0.0 ! Boundary condition 1
    !u_0(n) = 0.0 ! Boundary condition 2
    do i=2,n
       j = (i-1)*0.0005
       if (i<=k) then
          u_0(i) = j**(l+1)
       else
          u_0(i)=(1.0/(2*sqrt(pi)))*(j**(-1.0/4))*exp((-2.0/3)*(j**(3.0/2)))
       end if
    end do

    ! For testing
    open(unit=2,file='startti',status='unknown')
    do i=1,n
       write(2,*) u_0(i)
    end do
    close(2)
    

    ! Now use SOR to find best solution
    m = 0
    do while (m<1000000) ! Upper limit for iterations if no convergense
       u(1) = 0.0
       u(n) = u_0(n)
       uu = 0.0
       do i = 2,n-1
          j = (i-1)*0.0005
          ! Vector b element
          b_i = d_x*(((l*(l+1))*1.0/(j*j))+b*m2_h*j-a*m2_h/j-E_0*m2_h)*u_0(i)
          u(i) = (1-w)*u_0(i)-(w/2.0)*(b_i-u(i-1)-u_0(i+1)) ! Cosntruct next u
          ! Calculate (u-u_0)**2 for evaluation
          uu = uu+(u(i)-u_0(i))**2
       end do
       m = m+1

       if (mod(m,10000)==0) then
          print *,
          print *, "# m:", m
          print *, "||u-u_0||:", sqrt(uu)
       end if
       

       ! Evaluate ||u-u_0||
       if (sqrt(uu)<eps) then
          exit
       else
          u_0 = u
       end if
    end do

    deallocate(u_0)

    print *,
    print *, "Number of iterations:", m
    print *, "||u-u_0||:", sqrt(uu)
    print *,
  
  end subroutine qSOR

  !----------------------------------------------------------------------------
  !--------------------------- Thomas algorithm -------------------------------
  !----------------------------------------------------------------------------

  subroutine qTHOMAS(m_c, a, b, l, E_0, n, u)
    implicit none

    integer, intent(in) :: n, l                    ! # of elements and l
    double precision, intent(in) :: m_c, a, b, E_0 ! parameters for model
    double precision, intent(inout) :: u(n)  ! arrays

    double precision, allocatable :: a_i(:)
    double precision :: m2_h, j,  d_x, pi, r, aii
    integer :: i

    ! Allocate arrays
    allocate(a_i(n-2))

    ! Define parameter values
    m2_h = 51.3637898*m_c ! 1 divided by hbar/2m
    d_x = 0.0005**2  ! step size to the second power
    pi = 4*atan(1.0)
    r = 0.0005

    ! Thomas algorithm 
    u(1) = 0.0
    u(n) = (1.0/(2*sqrt(pi)))*((n*r)**(-1.0/4))*exp((-2.0/3)*((n*r)**(3.0/2)))
    print *, u(n)
       
    ! Forward sweep
    open(unit=99, file='a_i', status='unknown')
    a_i(1) = -2-d_x*(l*(l+1)*1.0/(r**2)-a*m2_h/r+b*m2_h*r-m2_h*E_0)
    write(99,*) a_i(1)
    do i=2,n-2
       j = i*r
       aii = -2-d_x*(l*(l+1)*1.0/(j**2)-a*m2_h/j+b*m2_h*j-m2_h*E_0)
       a_i(i) = aii-(1.0/a_i(i-1))
       write(99,*) a_i(i)
    end do
    close(99)
    
    ! Back substitution 
    u(n-1) = -u(n)/a_i(n-2)
    do i=n-2, 2, -1
       u(i) = -u(i+1)/a_i(i-1)
    end do

    deallocate(a_i)

  end subroutine qTHOMAS

  !----------------------------------------------------------------------------
  !--------------------------- Numerov Shooting -------------------------------
  !----------------------------------------------------------------------------

  subroutine qshooting(m_c, a, b, l, E_0, n, u, E)
    implicit none
    
    integer, intent(in) :: n, l                    ! # of elements 
    double precision, intent(in) :: m_c, a, b, E_0 ! parameters for model
    double precision, intent(inout) :: u(n)
    double precision, intent(out) :: E
    double precision, allocatable :: u_l(:), u_r(:)
    
    integer :: i, k, m, n_l, n_r
    double precision ::  m2_h, j, r, rr, d_E, pi
    double precision :: g_1, g_2, g_3, f, f_0, ff, gam_l, gam_r
    double precision :: c1, c2, c3, c4, V_i, V_0
    
    ! Define parameter values
    m2_h = 51.3637898*m_c ! 1 divided by hbar/2m
    r = 0.0005
    rr = r**2
    pi = 4*atan(1.0)
    
    allocate(u_l(n),u_r(n))
    
    ! Shooting & bisection method
    E = E_0
    d_E = 1.d-3
    u_l(1) = 0.0
    u_l(2) = r**(1.0*l+1.0)
    u_r(n) = 0.0
    u_r(n-1) = (1.0/(2*sqrt(pi)))*(((n-1)*r)**(-1.0/4))*exp((-2.0/3)*(((n-1)*r)**(3.0/2)))
    
    c1 = 1.0*l*(l+1.0)
    c2 = m2_h*a
    c3 = m2_h*b
    
    
    ! Determine the "turning point"
    do i=1,n
       j = i*r
       if (i==1) then
          V_0 = (c1/(j**2))-(c2/j)+(c3*j)
       else
          V_i = (c1/(j**2))-(c2/j)+(c3*j)
          if ((V_i>V_0).and.(V_i<=0)) then
             m = i
             V_0 = V_i
          end if
       end if
    end do
    
    k = 0
    
    do while (k<100000)
       c4 = m2_h*E
       ! Construct the wave function from both ends
       
       do i=3,n
          
          g_1=c1/((i*r)**2)-c2/(i*r)+c3*(i*r)-c4
          g_2=c1/(((i-1)*r)**2)-c2/((i-1)*r)+c3*(i-1)*r-c4
          g_3=c1/(((i-2)*r)**2)-c2/((i-2)*r)+c3*(i-2)*r-c4
          u_l(i) = (2.0*(1+5*rr*g_2/12)*u_l(i-1)-&
               (1-rr*g_3/12)*u_l(i-2))/(1-rr*g_1/12)
          
       end do
       
       do i=n-2,1,-1
          
          g_1=c1/((i*r)**2)-c2/(i*r)+c3*(i*r)-c4
          g_2=c1/(((i+1)*r)**2)-c2/((i+1)*r)+c3*(i+1)*r-c4
          g_3=c1/(((i+2)*r)**2)-c2/((i+2)*r)+c3*(i+2)*r-c4
          
          u_r(i) = (2.0*(1+5*rr*g_2/12)*u_r(i+1)-&
               (1-rr*g_3/12)*u_r(i+2))/(1-rr*g_1)
          
       end do
       
       ! Calculate the slope/u on the turning point
       gam_l = (u_l(m+1)-u_l(m))/(r*u_l(m))
       gam_r = (u_r(m)-u_r(m-1))/(r*u_r(m))
       f_0 = gam_l-gam_r
       
       ! Bisection
       if (k==0) then
          E = E+d_E
          f = f_0
          k = k+1
       else if (k>0) then
          if(f*f_0>0) then
             E = E+d_E
             f = f_0
             k = k+1
          else if (f*f_0<=0) then
             if (abs(f_0)>=r) then
                d_E = d_E*(-0.5)
                E = E+d_E
                f = f_0
                k = k+1
             else if (abs(f_0)<r) then
                k = k+1
                exit
             end if
          end if
       end if
       
       if (mod(k,100)==0) then
          print *, k, f_0, E, d_E
       end if
       
    end do
    
    print *,
    print *, k, f_0, E, d_E
    print *, u_l(m-1), u_l(m), u_l(m+1)
    print *, u_r(m-1)*(u_l(m)/u_r(m)), u_r(m)*(u_l(m)/u_r(m)), &
         u_r(m+1)*(u_l(m)/u_r(m))
    
    
    
    open(unit=1, file='u_l', status='unknown')
    open(unit=2, file='u_r', status='unknown')
    open(unit=3, file='u_tot', status='unknown')
    do i=1,n
       j=i*r
       write(1,*) u_l(i)
       write(2,*) u_r(i)*(u_l(m)/u_r(m))
       if (i<=m) then
          write(3,*) u_l(i)
       else
          write(3,*) u_r(i)*(u_l(m)/u_r(m))
       end if
    end do
    close(1)
    close(2)
    close(3)
    
    
    
    
  end subroutine qshooting
    
  
end module tridiag
       
       
       
    
