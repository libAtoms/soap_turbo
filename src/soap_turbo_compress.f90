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

module soap_turbo_compress_module

  contains

!**************************************************************************
  subroutine get_compress_indices( compress_mode, alpha_max, l_max, dim, P_nonzero, P_i, P_j, P_el, what_to_do )

  implicit none

! Input variables
  integer, intent(in) :: alpha_max(:), l_max
  character(*), intent(in) :: compress_mode, what_to_do

! Input-Output variables
  real*8, intent(inout) :: P_el(:)
  integer, intent(inout) :: dim, P_nonzero, P_i(:), P_j(:)

! Internal variables
  integer, allocatable :: pivot(:)
  integer :: n_species, n_max, i, counter, n, m, l, k
  logical :: set_indices
  real*8 :: f
  integer, allocatable :: ns(:,:), dense_to_sparse(:,:)
  integer :: nu_R, nu_S, N1, N2, S1, S2, n_compressed, n_uncompressed, n_1, n_2, alpha, beta, z1, z2, &
             a1, i_c, j, j_c, a2, ind, comp_ind, i_temp, alpha_2, beta_2, z1_2, z2_2, n_1_2, n_2_2, &
             a1_2, a2_2, i_c2, j_c2, comp_ind2, k2
  logical :: sym
  character*3 :: cm

  if( what_to_do == "get_dim" )then
    set_indices = .false.
  else if( what_to_do == "set_indices" )then
    set_indices = .true.
  else
    write(*,*) "ERROR: incorrect what_to_do in get_compress_indices() subroutine!"
    stop
  end if

  n_species = size(alpha_max)
  n_max = sum(alpha_max)

  if( set_indices )then
    P_i = 0
    P_j = 0
    P_el = 0.d0
  end if

  if( compress_mode == "trivial" .or. compress_mode == "Trivial" )then
!   Our traditional trivial compression recipes which has been used with all
!   GAPs with soap_turbo support released until at least early 2022
    allocate( pivot(1:n_species) )
    pivot = 0
    pivot(1) = 1

    do i = 1, n_species-1
      pivot(i+1) = pivot(i) + alpha_max(i)
    end do

    counter = 0
    k = 1

    do n = 1, n_max    
      do m = n, n_max    
        do l = 0, l_max    
          if( any(n == pivot) .or. any(m == pivot) )then
            counter = counter + 1
            if( set_indices )then
              P_i(counter) = counter
              P_j(counter) = k
              P_el(counter) = 1.d0
            end if
          end if
          k = k + 1
        end do
      end do
    end do
!   For trivial compression the SOAP vector dimension and the number of non-zero elements
!   in the transformation matrix are the same
    dim = counter
    P_nonzero = counter
    deallocate( pivot )
  else if( compress_mode == "0_0" .or. &
           compress_mode == "0_1" .or. &
           compress_mode == "0_2" .or. &
           compress_mode == "1_0" .or. &
           compress_mode == "1_1" .or. &
           compress_mode == "1_2" .or. &
           compress_mode == "2_0" .or. &
           compress_mode == "2_1" .or. &
           compress_mode == "2_2" )then
!   The compression recipes proposed by James Darby using nu_R = [0,1,2] (radially-sensitive
!   correlation index) and nu_S = [0,1,2] (species-sensitive correlation index)
    cm = adjustl(compress_mode)
    read(cm(1:1),*) nu_R
    read(cm(3:3),*) nu_S
!   At the moment this only works if all species have the same alpha_max
    do i = 1, n_species
      if( alpha_max(1) /= alpha_max(i) )then
        write(*,'(A,A,A)') "ERROR: compress_mode=", trim(adjustl(compress_mode)), &
                           " requires all species to have the same number of radial basis functions!"
        stop
      end if
    end do
!   The original code was written by James Darby in Python and ported to Fortran by Miguel Caro
!   Combined index for species (S) and radial basis functions (N)
    allocate( ns(1:2, 0:n_species*alpha_max(1)-1) )
    k = 0
    do i = 0, n_species-1
      do n = 0, alpha_max(1)-1
        ns(1, k) = i
        ns(2, k) = n
        k = k + 1
      end do
    end do
!   Limits for compressed SOAP
    sym = mod(nu_R,2) == 0 .and. mod(nu_S,2) == 0
    if( nu_R > 0 )then
      N1 = alpha_max(1)
    else
      N1 = 1
    end if
    if( nu_R == 2 )then
      N2 = alpha_max(1)
    else
      N2 = 1
    end if
    if( nu_S > 0 )then
      S1 = n_species
    else
      S1 = 1
    end if
    if( nu_S == 2 )then
      S2 = n_species
    else
      S2 = 1
    end if
    n_uncompressed = alpha_max(1)*n_species * (alpha_max(1)*n_species+1) / 2 * (l_max+1)
    if( .not. sym )then
      n_compressed = N1*S1 * N2*S2 * (l_max+1)
    else
      n_compressed = N1*S1 * (N1*S1+1) / 2 * (l_max+1)
    end if

    allocate( dense_to_sparse(0:n_compressed-1, 0:n_uncompressed-1) )
    dense_to_sparse = 0

    if( sym )then
      ind = 0
      counter = 0
      do i = 0, size(ns,2)-1
        alpha = ns(1, i)
        z1 = mod(alpha, S1)
        n_1 = ns(2, i)
        a1 = mod(n_1, N1)
        do j = 0, size(ns,2)-1
          beta = ns(1, j)
          z2 = mod(beta, S2)
          n_2 = ns(2, j)
          a2 = mod(n_2, N2)
!         The original power spectrum is symmetric
          if( j >= i )then
!           Compressions are sums across full power spectrum without symmetry
!           Non-trivial pre-factors to handle symmetry correctly
            do l = 0, l_max
              i_c = z1*N1 + a1
              j_c = z2*N2 + a2
!             Populate correct half of symmetric compressed power spectrum
              if( j_c < i_c )then
                i_temp = i_c
                i_c = j_c
                j_c = i_temp
              end if
              comp_ind = (i_c * (2*N1*S1+1-i_c) * (l_max+1))/2 + (j_c-i_c) * (l_max+1)+l
              if( dense_to_sparse(comp_ind, ind) == 0 )then
                counter = counter + 1
                dense_to_sparse(comp_ind, ind) = counter
              end if
              k = dense_to_sparse(comp_ind, ind)
              if( i /= j .and. i_c == j_c .and. set_indices )then
                P_el(k) = P_el(k) + dsqrt(2.d0)
                P_i(k) = comp_ind + 1
                P_j(k) = ind + 1
              else if( set_indices )then
                P_el(k) = P_el(k) + 1.d0
                P_i(k) = comp_ind + 1
                P_j(k) = ind + 1
              end if
              ind = ind + 1
            end do
          end if         
        end do
      end do
    else if( .not. sym )then
      ind = 0
      counter = 0
      do i = 0, size(ns,2)-1
        alpha = ns(1, i)
        z1 = mod(alpha, S1)
        n_1 = ns(2, i)
        a1 = mod(n_1, N1)
        i_c = z1*N1 + a1
        do j = 0, size(ns,2)-1
          beta = ns(1, j)
          z2 = mod(beta, S2)
          n_2 = ns(2, j)
          a2 = mod(n_2, N2)
          j_c = z2*N2 + a2

!         Compute i_c2 and j_c2, compressed indices for the symmetric term in the original power spectrum
!         that gets "skipped"
          alpha_2 = ns(1, j)
          z1_2 = mod(alpha_2, S1)
          n_1_2 = ns(2, j)
          a1_2 = mod(n_1_2, N1)
          beta_2 = ns(1, i)
          z2_2 = mod(beta_2, S2)
          n_2_2 = ns(2, i)
          a2_2 = mod(n_2_2, N2)
          i_c2 = z1_2*N1 + a1_2
          j_c2 = z2_2*N2 + a2_2
!         The  original power spectrum is symmetric
          if( j >= i )then
            do l = 0, l_max
              comp_ind = (i_c*S2*N2 + j_c) * (l_max+1) + l
              comp_ind2 = (i_c2*S2*N2 + j_c2) * (l_max+1) + l
              if( dense_to_sparse(comp_ind, ind) == 0 )then
                counter = counter + 1
                dense_to_sparse(comp_ind, ind) = counter
              end if
              k = dense_to_sparse(comp_ind, ind)
              if( dense_to_sparse(comp_ind2, ind) == 0 .and. i /= j .and. comp_ind /= comp_ind2 )then
                counter = counter + 1
                dense_to_sparse(comp_ind2, ind) = counter
              end if
              k2 = dense_to_sparse(comp_ind2, ind)
!             This will only contribute to a single element of the compressed SOAP vector
              if( i == j .and. set_indices )then
                P_el(k) = P_el(k) + 1.d0
                P_i(k) = comp_ind + 1
                P_j(k) = ind + 1
!             The diagonal image may contribute to 1 or 2 elements in the compressed SOAP vector
              else if( set_indices )then
                if( comp_ind == comp_ind2 )then
                  P_el(k) = P_el(k) + dsqrt(2.d0)
                  P_i(k) = comp_ind + 1
                  P_j(k) = ind + 1
                else
                  P_el(k) = P_el(k) + 1.d0/dsqrt(2.d0)
                  P_i(k) = comp_ind + 1
                  P_j(k) = ind + 1
                  P_el(k2) = P_el(k2) + 1.d0/dsqrt(2.d0)
                  P_i(k2) = comp_ind2 + 1
                  P_j(k2) = ind + 1
                end if
              end if
              ind = ind + 1
            end do
          end if
        end do
      end do
    end if

    P_nonzero = counter
    dim = n_compressed
    deallocate( ns, dense_to_sparse )
  else
    write(*,*) "ERROR: I don't understand compress_mode =", compress_mode
    stop
  end if

  end subroutine
!**************************************************************************

end module soap_turbo_compress_module
