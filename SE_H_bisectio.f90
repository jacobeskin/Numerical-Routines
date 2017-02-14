! Written by Jacob Eskin 26.5.2016
! Compilation gfortran -Wall -g -o SE_H_bisectio.exe SE_H_bisectio.f90

program SE_H_bisectio
implicit none

! This program calculates the spectrum of the Hydrogen atom
! using the simple bisection method, i.e. it calculates the 
! values of E for which the radial wavefunction goes to zero 
! far from the origin. Second derivative in the Schrödinger 
! equation is evaluated with the finite difference method leading
! to a recursion relation which is used when calculating the value
! the radial equation far from origin. The iteration of the radial 
! solution is separated into its own subroutine.


double precision :: u_2, u_e              ! Variables for SE iteration
double precision :: E_1, d_E, E_p1, E_p2, s, e, E_0 ! Energy variables & stuff 
integer :: n, l, m                        ! Quantum numbers & stuff

call cpu_time(s)

! Initial values for some of the loop variables
E_p1 = -14d+0
E_1 = E_p1
n = 1

! Begin the search for the energy levels
open(unit=1, file='Energylevels', status='unknown')
open(unit=2, file='Energylevelsminus', status='unknown')
do while (E_1<0)! Increase principal quantum number while energy is negatve
   write(1, '(a2,x,i2)'), "n=", n
   write(2, '(a2,x,i2)'), "n=", n
   print *, "n= ", n
   E_1 = E_p1
   m = n-1
   d_E = 0.001
   l = 0
   u_e = 0
   E_0 = -13.60569179989
   do while ((l<n).and.(E_1<0)) ! Loop over all azimuthal quantum numbers
      !print *, E_1, n, l
      call iteration(E_1, l, u_2) ! Iterate the SE solution
      
      ! If radial wf does not cross zero
      if ((u_2*u_e)>=0) then
         u_e = u_2
         E_1 = E_1+d_E
      end if
      
      ! If radial wf crosses zero
      if ((u_2*u_e)<0) then
         
         ! If desired precision is not reached, then adjust the energy level
         if (abs(d_E)>1d-15) then
            d_E = -d_E/2
            E_1 = E_1+d_E
            u_e = u_2
         end if
         
         ! If desired precision is reached for the enrgy level, write it 
         ! down and reset
         if (abs(d_E)<=1d-15) then
            print *, "E=", E_1-E_0, "l=", l
            write(1,'(a2,i1,a1,2x,f15.11, x, a2)'), 'l=',l,':',E_1-E_0, "eV"
            write(2,'(a2,i1,a1,2x,f15.11, x, a2)'), 'l=',l,':',E_1, "eV"
            d_E = 0.001
            u_e = 0

            if (n==1) then
               E_0 = E_1
               E_p1 = E_1+10
            else if (n>=2) then
               if (l==0) then
                  E_p2 = E_1
                  E_1 = E_p1
               else if ((l/=0).and.(l/=(n-1))) then
                  E_1 = E_p1
               else if ((l/=0).and.(l==(n-1))) then
                  E_p1 = E_p2+0.01
               end if
            endif
            l = l+1
         end if
      end if
      
   end do
   n = n+1
end do

close(1, status='keep')
close(2, status='keep')
call cpu_time(e)
print *, "Elapsed time: ", (e-s)/60

contains 
  
  subroutine iteration(E_1, l, u_2)
    implicit none

    double precision, intent(in) :: E_1
    integer, intent(in) :: l
    double precision, intent(out) :: u_2
    double precision :: u_0, u_1, e, V, Y, X, j
    integer :: i

    ! Values for the natural constants in units of eV and fm 
    e = 1.439964471d6
    X = 2.62468426d-11
    
    ! Iterating the finite difference solution up to 20000 fm for u_2
    u_0 = 0.0
    u_1 = 1.0
    u_2 = 0.0
    
    do i=1,10000000
       
       ! Calculate potential
       V = -e/i
       
       ! Calculate l-dependence
       if (l==0) then
          Y = 0.0
       else if (l/=0) then
          j = 1.0*i
          Y = (l**2+l)/(j**2)
       end if
       ! Calculate u_2
       u_2 = (2+Y+X*(V-E_1))*u_1-u_0
       ! Update u_0 and u_1
       u_0 = u_1
       u_1 = u_2
       
    end do
    
  end subroutine iteration

end program SE_H_bisectio



     

         

     
