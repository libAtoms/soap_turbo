! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   soap_turbo
! HND X
! HND X   soap_turbo is copyright (c) 2019-2023, Miguel A. Caro
! HND X
! HND X   soap_turbo is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   soap_turbo is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


module soap_turbo_radial_op

  contains

!**************************************************************************
!
! This function returns the normalization coefficients for the polynomial
! basis functions used to construct the orthonormal basis.
!
  function N_a(rcut, a)
    implicit none

    integer :: a, b
    real*8 :: rcut, N_a

    b = 2*a + 5

!   **************** New basis ******************
!    N_a = dsqrt( rcut**b / dfloat(b) )
    N_a = dsqrt( rcut / dfloat(b) )
!   *********************************************

  return
  end function
!**************************************************************************


!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set and a polynomial piecewise representation for 
! the atomic density
!
  subroutine get_radial_exp_coeff_operator_poly3(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, &
                                                 rcut_hard_in, atom_sigma_in, atom_sigma_scaling, &
                                                 amplitude_scaling, W, scaling_mode, mask, &
                                                 radial_enhancement, do_derivatives, do_central, &
                                                 central_weight, exp_coeff)

    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real*8, intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, atom_sigma_scaling
    real*8, intent(in) :: amplitude_scaling, central_weight
    real*8 :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, amplitude
    logical, intent(in) :: mask(:), do_derivatives, do_central
    character(*), intent(in) :: scaling_mode
!
    integer :: i, j, k
    real*8 :: pi, rj, s2, atom_width
    real*8 :: lim_soft(1:3), lim_buffer(1:3), B(1:7)
    real*8 :: M_left(1:7, 1:2), M_right(1:7, 1:2)
    real*8, allocatable :: A(:,:), I0(:,:), I_left(:), I_right(:)
    real*8 :: W(:,:)
!   Results will be stored in exp_coeff, which is an array of dimension (alpha_max, n_atom_pairs)
    real*8 :: exp_coeff(:,:)
    real*8, allocatable :: exp_coeff_soft(:), exp_coeff_buffer(:)
    logical, save :: print_basis = .false.
    real*8 :: amplitude_der

!   NOTE: the derivatives ARE MISSING !!!!!!!!
    allocate( exp_coeff_soft(1:alpha_max) )
    allocate( exp_coeff_buffer(1:alpha_max) )
    allocate( A(1:alpha_max, 1:7) )
    allocate( I0(1:alpha_max + 4, 1:3) )
    allocate( I_left(1:alpha_max) )
    allocate( I_right(1:alpha_max) )

!   This is for debugging. It prints the basis set to plot it with Gnuplot (gfortran only)
    if( print_basis )then
      print_basis = .false.
      write(*,*) "p(x,n,rcut) = (1.-x/rcut)**(n+2) / sqrt( rcut / (2.*n+5.) ) "
      do j=1, alpha_max
        write(*,"(A,I0,A)",advance="no") "p", j, "(x) = "
        do i = 1, alpha_max
          write(*,"(A,I2,A,F16.10,A,E16.8,A)",advance="no") "p(x,", i, ",", rcut_hard_in, ") *", W(j,i), "+"
        end do
        write(*,*) "0."
      end do
    end if

!   *************** New basis *******************
!   Normalized with rcut_hard to avoid instabilities
    rcut_soft = rcut_soft_in/rcut_hard_in
    rcut_hard = 1.d0
    atom_sigma = atom_sigma_in/rcut_hard_in
!   *********************************************
    pi = dacos(-1.d0)
    exp_coeff = 0.d0


    call get_constant_poly_coeff(alpha_max, rcut_hard, A)

    k = 0
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k = k + 1
!       Check if we need to do the central atom
        if( j == 1 .and. .not. do_central )then
          cycle
        end if
        if( rjs_in(k) < rcut_hard_in .and. mask(k) )then
!   **************** New basis ******************
          rj = rjs_in(k)/rcut_hard_in
!   *********************************************
          atom_sigma_scaled = atom_sigma + atom_sigma_scaling*rj
          s2 = atom_sigma_scaled**2
          atom_width = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaled
!   ~~~~~~~~~~~~~~~ Amplitude ~~~~~~~~~~~~~~~~~~~~~~~
          if( scaling_mode == "polynomial" )then
!           WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!           be taken into account in the corresponding part of the code.
!       
!           WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!           from the expressions
            if( amplitude_scaling == 0.d0 )then
              amplitude = 1.d0 / atom_sigma_scaled
              amplitude_der = - atom_sigma_scaling / s2 
            else if( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 <= 1.d-10 )then
              amplitude = 0.d0
              amplitude_der = 0.d0
            else
              if( amplitude_scaling == 1.d0 )then
                amplitude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )
                amplitude_der = 6.d0 / atom_sigma_scaled * (rj**2 - rj) &
                                - atom_sigma_scaling / atom_sigma_scaled * amplitude
              else
                amplitude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**amplitude_scaling
                amplitude_der = 6.d0*amplitude_scaling / atom_sigma_scaled * (rj**2 - rj) &
                                * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaled * amplitude
              end if
            end if
          end if
!         The central atom needs to be scaled by central_weight
          if( j == 1 )then
            amplitude = central_weight * amplitude
            amplitude_der = central_weight * amplitude_der
          end if
!         The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!         of atom j.
          if( radial_enhancement == 1 )then
            amplitude_der = amplitude * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                            amplitude_der * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
            amplitude = amplitude * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
          else if( radial_enhancement == 2 )then
            amplitude_der = amplitude*( 2.d0*rj + 2.d0*atom_sigma_scaled*atom_sigma_scaling + &
                                        dsqrt(8.d0/pi)*atom_sigma_scaled + dsqrt(8.d0/pi)*rj*atom_sigma_scaling ) + &
                            amplitude_der*( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
            amplitude = amplitude * ( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
          end if     
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          exp_coeff_soft = 0.d0
          exp_coeff_buffer = 0.d0
!         contribution inside rcut_soft from left and right part of the piecewise polynomial
          lim_soft = 0.d0
          if ( rj - atom_width < rcut_soft )then
!           integration limits inside r_soft
            lim_soft(1) = max( 0.d0, rj - atom_width )         ! lower limit left
            lim_soft(2) = min( max(0.d0, rj), rcut_soft )             ! upper limit left / lower limit right
            lim_soft(3) = min( rcut_soft, rj + atom_width ) ! upper limit right
!            write(*,'(A)') '- Inner zone:'
!            print *, lim_soft
!           contribution to the expansion coeff.
            I0 = transpose( M_radial_poly(lim_soft, alpha_max + 4, rcut_hard) ) ! +4 because we include poly alpha_max
            M_left(1:4, 1) = I0(1:4, 1) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(1), rj, atom_width, "left" )
            M_left(1:4, 2) = I0(1:4, 2) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(2), rj, atom_width, "left" )
            I_left = matmul( A(1:alpha_max, 1:4), M_left(1:4, 2) ) * I0(5:alpha_max + 4, 2) - &
                     matmul( A(1:alpha_max, 1:4), M_left(1:4, 1) ) * I0(5:alpha_max + 4, 1)
            M_right(1:4, 1) = I0(1:4, 2) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(2), rj, atom_width, "right" )
            M_right(1:4, 2) = I0(1:4, 3) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(3), rj, atom_width, "right" )
            I_right = matmul( A(1:alpha_max, 1:4), M_right(1:4, 2) ) * I0(5:alpha_max + 4, 3) - &
                      matmul( A(1:alpha_max, 1:4), M_right(1:4, 1) ) * I0(5:alpha_max + 4, 2)
            exp_coeff_soft = I_left + I_right
          end if
!         contribution in the buffer zone from left and right part of the piecewise polynomial
          lim_buffer = 0.d0
          if ( rcut_soft < rcut_hard .and. rj + atom_width > rcut_soft )then
!           integration limits inside r_soft
            lim_buffer(1) = max( rcut_soft, rj - atom_width ) ! lower limit left
            lim_buffer(2) = max( rj, rcut_soft )                     ! upper limit left / lower limit right
            lim_buffer(3) = min( rcut_hard, rj + atom_width ) ! upper limit right
!            write(*,'(A)') '* Buffer zone:'
!           contribution to the expansion coeff.
            I0 = transpose( M_radial_poly(lim_buffer, max(7, alpha_max + 4), rcut_hard) )
            call get_constant_poly_filter_coeff(rj, atom_width, rcut_soft, rcut_hard, 'left', B)
            M_left(1:7, 1) = matmul( -B(1:7), M_radial_monomial(lim_buffer(1), 6) ) * I0(1:7, 1)
            M_left(1:7, 2) = matmul( -B(1:7), M_radial_monomial(lim_buffer(2), 6) ) * I0(1:7, 2)
            I_left = matmul( A(1:alpha_max, 1:7), M_left(1:7, 2) ) * I0(5:alpha_max + 4, 2) - &
                     matmul( A(1:alpha_max, 1:7), M_left(1:7, 1) ) * I0(5:alpha_max + 4, 1)
            call get_constant_poly_filter_coeff(rj, atom_width, rcut_soft, rcut_hard, 'right', B)
            M_right(1:7, 1) = matmul( -B(1:7), M_radial_monomial(lim_buffer(2), 6) ) * I0(1:7, 2)
            M_right(1:7, 2) = matmul( -B(1:7), M_radial_monomial(lim_buffer(3), 6) ) * I0(1:7, 3)
            I_right = matmul( A(1:alpha_max, 1:7), M_right(1:7, 2) ) * I0(5:alpha_max + 4, 3) - &
                      matmul( A(1:alpha_max, 1:7), M_right(1:7, 1) ) * I0(5:alpha_max + 4, 2)
            exp_coeff_buffer = I_left + I_right
          end if
!         Transform from g_alpha to g_n (the orthonormal basis)
          exp_coeff(1:alpha_max, k) = amplitude * matmul( W, exp_coeff_soft(1:alpha_max) + exp_coeff_buffer(1:alpha_max) )
        end if
      end do
    end do

!   **************** New basis ******************
!   This results from the change of variable in the
!   overlap integrals. We only need this if we want to
!   know the actual value of the expansion coefficients.
!   Since this is a global factor, once we have normalized
!   the SOAP vectors it does not have an effect anymore.
    exp_coeff = exp_coeff * dsqrt(rcut_hard_in)
!   *********************************************

!   This is for debugging
    if( .false. )then
      open(10, file="radial_expansion_coefficients.dat", status="unknown")
      do i=1, size(exp_coeff, 2)
        write(10,*) exp_coeff(1:alpha_max, i)
      end do
      close(10)
    end if

    deallocate( exp_coeff_soft, exp_coeff_buffer, A, I0, I_left, I_right)

  return
  end subroutine
!**************************************************************************


!**************************************************************************
!
! Auxiliary function that computes a piecewise polynomial function and
! its derivatives.
!
  function g_aux(r, r0, width, piece) result(poly)
    implicit none
 
    real*8, intent(in) :: r, r0, width
    character(*), intent(in) :: piece
    real*8 :: x
    real*8 :: poly(1:4)

    x = (r - r0)/width

    if ( piece == "left" )then
      poly = [ 1.d0 - 3.d0*x**2 - 2.d0*x**3, -6.d0*(x**2 + x)/width, -3.d0*(2.d0*x + 1)/width**2, -2.d0/width**3 ]
    else if ( piece == "right" )then
      poly = [ 1.d0 - 3.d0*x**2 + 2.d0*x**3, 6.d0*(x**2 - x)/width, 3.d0*(2.d0*x - 1)/width**2, 2.d0/width**3 ]
    end if

  return
  end function
!**************************************************************************

!**************************************************************************
!
! This subroutine returns the matrix containing the atom-independent
! coefficients for the polynomial radial basis.
! (in the notes: A, A*)
!
  subroutine get_constant_poly_coeff(alpha_max, rcut, A)
    implicit none

    integer, intent(in) :: alpha_max
    real*8, intent(in) :: rcut
    integer :: i, j, l(1:2)
    integer, allocatable :: factor(:)
    real*8, intent(inout) :: A(:,:)

    allocate( factor(1:alpha_max) )
    
    do i = 1, alpha_max
      A(i, 1) = rcut / N_a( rcut, i ) / dfloat(i + 3)
      factor(i) = (i + 3)
    end do

    l = shape(A)
    do j = 2, l(2)
        A(:, j) = A(:, j - 1) * rcut / dfloat(factor + j - 1) 
    end do

  end subroutine
!**************************************************************************


!**************************************************************************
!
! This subroutine returns the constant coefficients for the 
! contribution from the buffer zone, coming from the dot products
! of poly. radial basis and smoothing function (filter)
! (in the notes: B_l, B_r)
!
  subroutine get_constant_poly_filter_coeff(rj, sigma_j, rcut_soft, rcut_hard, piece, B)
    implicit none

    real*8, intent(in) :: rj, sigma_j, rcut_soft, rcut_hard
    character(*), intent(in) :: piece
    integer :: i
    real*8 :: C_filter(1:4), col_poly(1:4), C_poly(1:7,1:4)
    real*8, intent(inout) :: B(:)

!   coeff. from the filter
    C_filter = g_aux(0.d0, rcut_soft, 2.d0*sqrt(2.d0*log(2.d0))*(rcut_hard - rcut_soft), "right")

!   build Toeplitz matrix a.k.a. diagonal-constant matrix
    C_poly = 0.d0
    col_poly(1:4) = g_aux(0.d0, rj, sigma_j, piece)

    do i = 1, 4
      C_poly(i:i+3, i) = col_poly
    end do

    B = matmul( C_poly, C_filter )

  end subroutine
!**************************************************************************


!**************************************************************************
!
! This subroutine returns the radial terms coming from the
! polynomial radial basis 
! (in the notes: I0, R_hard)
!
  function M_radial_poly(r, alpha_max, rcut) result(radial_terms)
    implicit none

    integer, intent(in) :: alpha_max
    real*8, intent(in) :: rcut, r(:)
    integer :: i
    real*8, dimension(1:size(r), 1:alpha_max) :: radial_terms

    radial_terms = 1.d0
    do i = 2, alpha_max
      radial_terms(:, i) = radial_terms(:, i - 1) * (1.d0 - r/rcut)
    end do

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns a vector of monomial terms (1, x, x**2, ...)
! for a given polynomial degree
! (in the notes: R, R_ext) 
!
  function radial_monomial(r, degree) result(radial_terms)
    implicit none

    integer, intent(in) :: degree
    real*8, intent(in) :: r
    integer :: p
    real*8, dimension(1:degree + 1) :: radial_terms 

    radial_terms = 1.d0
    do p = 2, degree + 1
      radial_terms(p) = r * radial_terms(p - 1) 
    end do

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns the matrix of monomial terms and their
! derivatives, for a given polynomial degree:
!        ( 1  x  x**2    ...  x**degree              )
!        ( 0  1  2 * x   ...  degree * x**(degree-1) )
!        ( .   .                                     )
!        ( .       .                                 )
!        ( .             .    degree!                )
!
! (in the notes: R*) 
!
  function M_radial_monomial(r, degree) result(M_radial)
    implicit none

    integer, intent(in) :: degree
    real*8, intent(in) :: r
    real*8, dimension(1:degree + 1) :: radial_terms
    real*8, dimension(1:degree + 1, 1:degree + 1) :: M_radial

    radial_terms = radial_monomial(r, degree)

    M_radial = 0.d0
    M_radial(1:degree + 1, 1) = radial_terms(1:7)
    M_radial(2:degree + 1, 2) = radial_terms(1:6) * [1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0]
    M_radial(3:degree + 1, 3) = radial_terms(1:5) * [2.d0, 6.d0, 12.d0, 20.d0, 30.d0]
    M_radial(4:degree + 1, 4) = radial_terms(1:4) * [6.d0, 24.d0, 60.d0, 120.d0]
    M_radial(5:degree + 1, 5) = radial_terms(1:3) * [24.d0, 120.d0, 360.d0]
    M_radial(6:degree + 1, 6) = radial_terms(1:2) * [120.d0, 720.d0]
    M_radial(7:degree + 1, 7) = 720.d0

  end function
!**************************************************************************


end module soap_turbo_radial_op
