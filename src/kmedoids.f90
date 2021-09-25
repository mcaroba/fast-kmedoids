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
  else if( size(M_init) > size(D,1) )then
    write(*,*) "ERROR: number of medoids can't be larger than number of data points!"
    return
  else if( any(M_init < 0) .or. any(M_init >= size(D,1)) )then
    write(*,*) "ERROR: bad medoids initialization (elements in M_init out of range)!"
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






subroutine find_isolated(D, nclusters, niso, M)

  implicit none

! Input variables
  real*8, intent(in) :: D(:,:)
  integer, intent(in) :: nclusters, niso
! Output variables
  integer, intent(out) :: M(1:nclusters)
! Internal variables
  real*8, allocatable :: summed_D(:)
  real*8 :: av_D, rand
  integer :: npts, i, j, k, l, tries
  logical :: is_j_in_M, is_j_too_close

  if( niso < 1 )then
    write(*, *) "WARNING: niso < 1, so nothing to do!"
    return
  end if

  npts = size(D,2)
  allocate( summed_D(1:npts) )

  av_D = sum(D) / float(npts)**2

! Find the first one
  do i = 1, npts
    summed_D(i) = sum( D(1:npts, i) )
  end do
  i = maxloc(summed_D, 1)
  M(1) = i

! Find the rest
  do i = 2, niso
    summed_D = 0.d0
    do j = 1, i-1
      k = M(j)
      do l = 1, npts
        if( l /= k )then
          summed_D(l) = summed_D(l) + D(l, k)
        end if
      end do
    end do
    j = maxloc(summed_D, 1)
    M(i) = j
  end do

! Randomize the rest
  do i = niso+1, nclusters
    tries = 0
    do
      call random_number(rand)
      j = floor(rand * float(npts)) + 1
!     Check if j is already in M
      is_j_in_M = .false.
      do k = 1, i-1
        if( M(k) == j )then
          is_j_in_M = .true.
          exit
        end if
      end do
!     Check if j is too close
      if( .not. is_j_in_M )then
        is_j_too_close = .false.
        do k = 1, i-1
          if( D(j,M(k)) < av_D )then
            is_j_too_close = .true.
            tries = tries + 1
            exit
          end if
        end do      
      end if
!     Relax search criteria
      if( tries == 100 )then
        tries = 0
        av_D = 0.9d0 * av_D
      end if
!     In case of success, we break the loop
      if( .not. is_j_in_M .and. .not. is_j_too_close )then
        M(i) = j
        exit
      end if
    end do
  end do

  M = M - 1

  deallocate( summed_D )

end subroutine






subroutine find_random(npts, nclusters, M_input, is_empty, M)

  implicit none

! Input variables
  integer, intent(in) :: npts, nclusters
  integer, intent(in) :: M_input(:)
  logical, intent(in) :: is_empty
! Output variables
  integer, intent(out) :: M(1:nclusters)
! Internal variables
  real*8 :: rand
  integer :: nstart, i, j, k
  logical :: is_j_in_M

  if( .not. is_empty )then
    nstart = size(M_input)
    if( size(M_input) > nclusters )then
      write(*,*) "ERROR: size(M_input) can't be larger than nclusters"
      return
    end if
  else
    nstart = 0
  end if

! Fix Python convention; we fix it back at the end
  if( nstart > 0 )then
    M(1:nstart) = M_input(1:nstart) + 1
  end if

! Randomize the rest
  do i = nstart+1, nclusters
    do
      call random_number(rand)
      j = floor(rand * float(npts)) + 1
!     Check if j is already in M
      is_j_in_M = .false.
      do k = 1, i-1
        if( M(k) == j )then
          is_j_in_M = .true.
          exit
        end if
      end do
!     In case of success, we break the loop
      if( .not. is_j_in_M )then
        M(i) = j
        exit
      end if
    end do
  end do

  M = M - 1


end subroutine





subroutine fast_kmedoids_iterative(D, nclusters, tmax, n_inits, n_iso, incoherence, M_in, M_init_mode, verbosity, M, C)

  implicit none

! Input variables
  real*8, intent(in) :: D(:, :)
  integer, intent(in) :: tmax, n_inits, n_iso, nclusters
  integer, intent(in) :: M_in(:), verbosity
  character(len=*), intent(in) :: incoherence
  character(len=*), intent(in) :: M_init_mode
! Output variables
  integer, intent(out) :: C(1:size(D,1)), M(1:nclusters)
! Internal variables
  real*8 :: R, R_prev
  integer, allocatable :: M_init(:)
  integer :: C_i(1:size(D,1)), M_i(1:nclusters), M_ii(1:nclusters)
  integer :: nstart, npts, i
  logical :: is_empty

  npts = size(D,1)

  if( M_init_mode == "random" )then
    is_empty = .true.
  else if( M_init_mode == "random+" )then
    nstart = size(M_in)
    allocate( M_init(1:nstart) )
    M_init(1:nstart) = M_in(1:nstart)
    is_empty = .false.
  end if

  if( M_init_mode(1:6) /= "random" .and. M_init_mode(1:8) /= "isolated" )then
    write(*,*) "ERROR: I don't understand M_init_mode, it must be random|isolated"
    return
  end if

  R_prev = 1.d100
  if( verbosity >= 1 )then
    write(*,"(A,A,A)") "# Iteration; Residual ('", trim(adjustl(incoherence)), "' incoherence criterion)"
  end if
!$OMP parallel do schedule(static,1) private(M_ii,M_i,C_i,R)
  do i = 1, n_inits
    if( M_init_mode(1:6) == "random" )then
!     Set up initial M
      call find_random(npts, nclusters, M_init, is_empty, M_ii)
    else if( M_init_mode == "isolated" )then
!     Set up initial M
      call find_isolated(D, nclusters, n_iso, M_ii)
    end if
    call fast_kmedoids(D, tmax, M_ii, M_i, C_i)
    call get_incoherence(D, M_i, C_i, incoherence, R)
!$OMP critical
    if( R < R_prev )then
      M = M_i
      C = C_i
      R_prev = R
      if( verbosity >= 1 )then
        write(*,*) i, R
      end if
    end if
!$OMP end critical
  end do
!$OMP end parallel do

  if( allocated(M_init) )deallocate( M_init )

end subroutine






subroutine get_incoherence(D, M, C, incoherence, R)

  implicit none

! Input variables
  real*8, intent(in) :: D(:,:)
  integer, intent(in) :: M(:), C(:)
  character(len=*) :: incoherence
! Output variables
  real*8, intent(out) :: R
! Internal variables
  integer, allocatable :: cluster_size(:)
  integer :: npts, nclusters, i, j, k, k2

  npts = size(D, 1)
  nclusters = size(M)

  if( incoherence == "rel" )then
    allocate( cluster_size(1:nclusters) )

    cluster_size = 0
    do i = 1, npts
      cluster_size(C(i)+1) = cluster_size(C(i)+1) + 1
    end do

    R = 0.d0
    do i = 1, nclusters
      k = M(i)+1
      do j = 1, npts
        k2 = C(j) + 1
        R = R + D(k2, k)/float(cluster_size(i)) 
      end do
    end do
    deallocate( cluster_size )
  else if( incoherence == "tot" )then
    R = 0.d0
    do i = 1, nclusters
      k = M(i)+1
      do j = 1, npts
        k2 = C(j) + 1
        R = R + D(k2, k) 
      end do
    end do
  end if

end subroutine

end module
