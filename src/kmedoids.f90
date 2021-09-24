! Copyright (c) 2021 Miguel A. Caro
! Based on C. Bauckhage's Python implementation of kmedoids

module kmedoids_module

contains

subroutine fast_kmedoids(D, tmax, M_init, M, C)

  implicit none

! Input variables
  real*8, intent(in) :: D(:, :)
  integer, intent(in) :: tmax
  integer, intent(in) :: M_init(:)
! Output variables
  integer, intent(out) :: C(1:size(D,1)), M(1:size(M_init))
! Internal variables
  real*8, allocatable :: D_intracluster(:,:), average_D(:)
  real*8 :: d_min
  integer, allocatable :: M_new(:), cluster_size(:), cluster_map(:)
  integer :: t, npts, nclusters, i, j, k, i2, j2

  if( size(D,1) /= size(D,2) )then
    write(*,*) "ERROR: distance matrix D must be square!"
    return
  else if( size(M) > size(D,1) )then
    write(*,*) "ERROR: number of medoids can't be larger than number of data points!"
    return
  else if( any(M < 0) .or. any(M >= size(D,1)) )then
    write(*,*) "ERROR: bad medoids initialization (elements in M out of range)!"
    return
  else if( tmax < 0 )then
    write(*,*) "WARNING: tmax < 0, nothing to do!"
    return
  end if

  M = M_init

  npts = size(D, 1)
  nclusters = size(M)

  allocate( cluster_size(1:nclusters) )
  allocate( M_new(1:nclusters) )
  C = 0
! Fix Python convention internally (we change it back at the end)
  M = M + 1

  if( tmax == 0 )then
    t = 0
  else
    t = 1
  end if

  do while(t <= tmax)
!   Find the closest medoid to each point
    cluster_size = 0
    do i = 1, npts
      d_min = 1.d100
      do k = 1, nclusters
        j = M(k)
        if( D(j,i) < d_min )then
          d_min = D(j,i)
          C(i) = k
        end if
      end do
      cluster_size(C(i)) = cluster_size(C(i)) + 1
    end do

    if( t == 0 )return

!   Update clusters
    do k = 1, nclusters
!     Update cluster map
      allocate( cluster_map(1:cluster_size(k)) )
      j = 0
      do i = 1, npts
        if( C(i) == k )then
          j = j + 1
          cluster_map(j) = i
        end if
      end do
!     Compute intracluster distances
      allocate( D_intracluster(1:cluster_size(k), 1:cluster_size(k)) )
      allocate( average_D(1:cluster_size(k)) )
      do i = 1, cluster_size(k)
        i2 = cluster_map(i)
        do j = 1, cluster_size(k)
          j2 = cluster_map(j)
          D_intracluster(j, i) = D(j2, i2)
        end do
        average_D(i) = sum( D_intracluster(1:cluster_size(k),i) )
      end do
!     Find the new medoid of each cluster
      i = minloc( average_D, 1 )
      M_new(k) = cluster_map(i)
      deallocate( D_intracluster, cluster_map, average_D )
    end do
    if( all(M_new == M) )then
      M = M_new
      M = M - 1
      C = C - 1
      return
    end if

    M = M_new
    t = t + 1
  end do

  M = M - 1
  C = C - 1

end subroutine

end module
