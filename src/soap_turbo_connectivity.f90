! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   soap_turbo
! HND X
! HND X   soap_turbo is copyright (c) 2019-2024, Miguel A. Caro
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


module connectivity_module

  implicit none

  contains

    recursive subroutine find_neighbors(i, j, n_atoms, bonded, atom_visited, atom_belongs_to_cluster, cluster)

      implicit none

!     Input variables
      integer, intent(in) :: i, j, n_atoms, cluster
      logical, intent(in) :: bonded(:,:)
!     Inout variables
      integer, intent(inout) :: atom_belongs_to_cluster(:)
      logical, intent(inout) :: atom_visited(:)
!     Internal variables
      integer :: k

      if( .not. atom_visited(j) .and. i /= j .and. bonded(i,j) )then
        atom_visited(j) = .true.
        atom_belongs_to_cluster(j) = cluster
        do k = 1, n_atoms
          call find_neighbors(j, k, n_atoms, bonded, atom_visited, atom_belongs_to_cluster, cluster)
        end do
      end if
    end subroutine




    subroutine cluster_atoms(positions, cutoff_hard, cutoff_soft, connectivity_weights)

      implicit none

!     Input variables
      real*8, intent(in) :: positions(:,:), cutoff_hard(:), cutoff_soft(:)
!     Output variables
      real*8, intent(out), dimension(1:size(positions,2)) :: connectivity_weights
!     Internal variables
      integer, dimension(1:size(positions,2)) :: atom_belongs_to_cluster_hard
      integer, dimension(1:size(positions,2)) :: atom_belongs_to_cluster_soft
      logical, dimension(1:size(positions,2),1:size(positions,2)) :: bonded_hard, bonded_soft
      logical, dimension(1:size(positions,2)) :: atom_visited
      real*8, dimension(1:size(positions,2)) :: cluster_weights
      real*8, dimension(1:size(positions,2),1:size(positions,2)) :: pairwise_connectivity_weights
      real*8 :: d, cw, cs, ch, x, y, z
      integer :: n_atoms, cluster, i, j

      n_atoms = size(positions,2)

      bonded_hard = .false.
      bonded_soft = .false.
      pairwise_connectivity_weights = 0.d0

      !$omp parallel do private(i,j,d)
      do i = 1, n_atoms
        do j = i+1, n_atoms
          x = positions(1, j) - positions(1, i)
          y = positions(2, j) - positions(2, i)
          z = positions(3, j) - positions(3, i)
          d = sqrt( x**2 + y**2 + z**2)
          cs = (cutoff_soft(i)+cutoff_soft(j))/2.d0
          ch = (cutoff_hard(i)+cutoff_hard(j))/2.d0
          if( d < ch )then
!           Partly connected atoms
            bonded_hard(i, j) = .true.
            bonded_hard(j, i) = .true.
            if( d < cs )then
              cw = 1.d0
!             Fully connected atoms
              bonded_soft(i, j) = .true.
              bonded_soft(j, i) = .true.
            else
              cw = 1.d0-3.d0*((d-cs)/(ch-cs))**2+2.d0*((d-cs)/(ch-cs))**3
            end if
            pairwise_connectivity_weights(i, j) = cw
            pairwise_connectivity_weights(j, i) = cw
          end if
        end do
      end do

!     Construct clusters according to hard cutoff
      atom_belongs_to_cluster_hard = 0
      atom_visited = .false.
      cluster = 0
      do i = 1, n_atoms
!     Find all the atoms that can be connected to i and put them in the same cluster
        if( .not. atom_visited(i) )then
          cluster = cluster + 1
          atom_belongs_to_cluster_hard(i) = cluster
          atom_visited(i) = .true.
          do j = 1, n_atoms
            call find_neighbors(i, j, n_atoms, bonded_hard, atom_visited, atom_belongs_to_cluster_hard, cluster)
          end do
        end if
      end do

!     Construct clusters according to soft cutoff
      atom_belongs_to_cluster_soft = 0
      atom_visited = .false.
      cluster = 0
      do i = 1, n_atoms
!     Find all the atoms that can be connected to i and put them in the same cluster
        if( .not. atom_visited(i) )then
          cluster = cluster + 1
          atom_belongs_to_cluster_soft(i) = cluster
          atom_visited(i) = .true.
          do j = 1, n_atoms
            call find_neighbors(i, j, n_atoms, bonded_soft, atom_visited, atom_belongs_to_cluster_soft, cluster)
          end do
        end if
      end do

!     Check which atoms are:
!         1) fully connected to central cluster (within soft cutoff)
!         2) completely unconnected to central cluster (beyond hard cutoff)
!         3) partly connected to central cluster
      cluster_weights(1) = 1.d0
      cluster_weights(2:n_atoms) = 0.d0
      do i = 1, n_atoms
        if( atom_belongs_to_cluster_soft(i) == 1 )then
          continue
        else if( atom_belongs_to_cluster_hard(i) /= 1 )then
          continue
        else
          cluster = atom_belongs_to_cluster_soft(i)
!         This can be made more efficient by looping only over atoms that are (soft) in cluster 1
          do j = 1, n_atoms
            if( atom_belongs_to_cluster_soft(j) == 1 .and. bonded_hard(i,j) )then
              if( pairwise_connectivity_weights(i,j) > cluster_weights(cluster) )then
                cluster_weights(cluster) = pairwise_connectivity_weights(i,j)
              end if
            end if
          end do
        end if
      end do
      do i = 1, n_atoms
        if( atom_belongs_to_cluster_soft(i) == 1 )then
          connectivity_weights(i) = 1.d0
        else if( atom_belongs_to_cluster_hard(i) /= 1 )then
          connectivity_weights(i) = 0.d0
        else
          cluster = atom_belongs_to_cluster_soft(i)
          connectivity_weights(i) = cluster_weights(cluster)
        end if
      end do

    end subroutine


end module
