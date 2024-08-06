! CMA-ES Module
! The Covariance Matrix Adaptation Evolution Strategy (CMA-ES) with
! - Separable CMA for large scale problems [Ros 2008]
! - (Improved) Box Constraint Handling [Hansen 2009]
! - (Improved) Active-CMA [Jastrebski 2006]
!
! The current version can solve box constrained optimization problem of the form
!   Minimize f(x), x = (x_1,...,x_n)
!   Subject to lb <= x_i <= ub for i = 1,...,n
!
! Written and maintained by Youhei Akimoto @ Shinshu University (yohe.aki@gmail.com)
!
! BUG fixed (2018/01/19): array initializations (to prevent nan * 0)
! BUG fixed (2018/01/16): indexsort (local array memory allocation)
! BUG fixed (2017/04/19): tolf
! Interface Change (2017/04/19): Input File
! BUG fixed (2016/12/26): effective_median_iqr.
! BUG fixed (2016/11/14): inout -> out. Thanks to Prof. Fujii and Mr. Takahashi.

module cmaes
  use random, only : normal, uniform, uniformint, randperm
  use eigen, only : myddot, mydger, mydsyr, mydgemv, mydgemm, mydsyev 
  implicit none
  
  type CmaParam
     ! Static Parameters
     integer :: dim, lambda, mu, flag_diagonal, flag_active, fhistory_size
     double precision :: flg_min
     double precision :: mueff, cm, cc, cone, cmu, cs, ds
     double precision :: dble_dim, chidim
     double precision, allocatable :: weights(:), neg_weights(:) ! mu-D
     ! Dynamic Parameters
     integer :: iter, neval, eig_iter
     double precision :: sigma
     double precision, allocatable :: xmean(:), vecD(:), ps(:), pc(:), diagSqrtC(:)
     double precision, allocatable :: matC(:, :), matB(:, :), matInvSqrtC(:, :)
     ! Box Constraint Parameters
     integer :: idx_fhistory, flg_fhistory, flg_winit, cnt_winit
     double precision, allocatable :: ub(:), lb(:), penalty(:)  ! dim-D
     double precision, allocatable :: fhistory(:)
     ! Termination Condition
     integer :: flg_dsyev_err
     ! Work Space
     integer :: lwork
     double precision, allocatable :: work(:)
     double precision, allocatable :: workmat(:, :)
     ! Mixed-integer optimization Parameters
     integer :: flg_mixed_integer
     double precision :: sigmaCorrectionFactor
     double precision, allocatable :: vecStairWidth(:), bestsolutionX(:)
  end type CmaParam
    
contains

  subroutine writecp( cp, fname )
    implicit none
    type(CmaParam), intent(in) :: cp
    character(*), intent(in) :: fname    
    integer, parameter :: fid = 10

    open( fid, file=fname, action='write', status='replace', form='unformatted' )    
    ! Static paramters
    write( fid ) cp%dim
    write( fid ) cp%lambda
    write( fid ) cp%mu
    write( fid ) cp%flag_diagonal
    write( fid ) cp%flag_active
    write( fid ) cp%fhistory_size
    write( fid ) cp%flg_min
    write( fid ) cp%mueff
    write( fid ) cp%cm
    write( fid ) cp%cc
    write( fid ) cp%cone
    write( fid ) cp%cmu
    write( fid ) cp%cs
    write( fid ) cp%ds
    write( fid ) cp%dble_dim
    write( fid ) cp%chidim
    write( fid ) cp%weights
    write( fid ) cp%neg_weights
    ! Dynamic parameters
    write( fid ) cp%iter
    write( fid ) cp%neval
    write( fid ) cp%eig_iter
    write( fid ) cp%sigma
    write( fid ) cp%xmean
    write( fid ) cp%vecD
    write( fid ) cp%ps
    write( fid ) cp%pc
    write( fid ) cp%diagSqrtC
    if ( cp%flag_diagonal == 0 ) then
       write( fid ) cp%matC
       write( fid ) cp%matB
       write( fid ) cp%matInvSqrtC
    end if
    ! Box Constraint
    write( fid ) cp%idx_fhistory
    write( fid ) cp%flg_fhistory    
    write( fid ) cp%flg_winit
    write( fid ) cp%cnt_winit
    write( fid ) cp%ub
    write( fid ) cp%lb
    write( fid ) cp%penalty
    write( fid ) cp%fhistory
    ! Mixed-Integer optimization
    write( fid ) cp%flg_mixed_integer
    write( fid ) cp%sigmaCorrectionFactor
    write( fid ) cp%vecStairWidth
    write( fid ) cp%bestsolutionX
    ! Others
    write( fid ) cp%flg_dsyev_err
    if ( cp%flag_diagonal == 0 ) then
       write( fid ) cp%lwork
       write( fid ) cp%work
       write( fid ) cp%workmat
    end if
    close( fid )
  end subroutine writecp

  subroutine readcp( cp, fname )
    implicit none
    type(CmaParam), intent(out) :: cp
    character(*), intent(in) :: fname    
    integer, parameter :: fid = 10
    
    open( fid, file=fname, action='read', form='unformatted' )
    ! Static paramters
    read( fid ) cp%dim
    read( fid ) cp%lambda
    read( fid ) cp%mu
    read( fid ) cp%flag_diagonal
    read( fid ) cp%flag_active
    read( fid ) cp%fhistory_size
    allocate(cp%weights(cp%mu))
    if ( cp%flag_active == 1 ) then
       allocate(cp%neg_weights(cp%mu))
    endif
    allocate(cp%xmean(cp%dim), cp%ps(cp%dim), cp%pc(cp%dim), cp%vecD(cp%dim), cp%diagSqrtC(cp%dim))
    allocate(cp%lb(cp%dim), cp%ub(cp%dim), cp%penalty(cp%dim))
    allocate(cp%fhistory(cp%fhistory_size))
    read( fid ) cp%flg_min
    read( fid ) cp%mueff
    read( fid ) cp%cm
    read( fid ) cp%cc
    read( fid ) cp%cone
    read( fid ) cp%cmu
    read( fid ) cp%cs
    read( fid ) cp%ds
    read( fid ) cp%dble_dim
    read( fid ) cp%chidim
    read( fid ) cp%weights
    read( fid ) cp%neg_weights
    ! Dynamic parameters
    read( fid ) cp%iter
    read( fid ) cp%neval
    read( fid ) cp%eig_iter
    read( fid ) cp%sigma
    read( fid ) cp%xmean
    read( fid ) cp%vecD
    read( fid ) cp%ps
    read( fid ) cp%pc
    read( fid ) cp%diagSqrtC
    if ( cp%flag_diagonal == 0 ) then
       allocate(cp%matB(cp%dim, cp%dim), cp%matC(cp%dim, cp%dim), &
            cp%matInvSqrtC(cp%dim, cp%dim))
       read( fid ) cp%matC
       read( fid ) cp%matB
       read( fid ) cp%matInvSqrtC
    end if
    ! Box Constraint
    read( fid ) cp%idx_fhistory
    read( fid ) cp%flg_fhistory    
    read( fid ) cp%flg_winit
    read( fid ) cp%cnt_winit
    read( fid ) cp%ub
    read( fid ) cp%lb
    read( fid ) cp%penalty
    read( fid ) cp%fhistory
    ! Mixed-Integer optimization
    allocate(cp%vecStairWidth(cp%dim), cp%bestsolutionX(cp%dim))
    read( fid ) cp%flg_mixed_integer
    read( fid ) cp%sigmaCorrectionFactor
    read( fid ) cp%vecStairWidth
    read( fid ) cp%bestsolutionX
    ! Others
    read( fid ) cp%flg_dsyev_err
    if ( cp%flag_diagonal == 0 ) then
       read( fid ) cp%lwork
       allocate(cp%work(cp%lwork))
       allocate(cp%workmat(cp%dim, max(cp%dim, cp%lambda)))       
       read( fid ) cp%work
       read( fid ) cp%workmat
    end if
    close( fid )
  end subroutine readcp
  
  subroutine printcp( cp )
    implicit none
    type(CmaParam), intent(inout) :: cp
    
    ! Static paramters
    write( *, * )  'dim = ', cp%dim
    write( *, * )  'lambda = ', cp%lambda
    write( *, * )  'mu = ', cp%mu
    write( *, * )  'flag_diagonal = ', cp%flag_diagonal
    write( *, * )  'flag_active = ', cp%flag_active
    write( *, * )  'fhistory_size', cp%fhistory_size    

    write( *, * )  'flg_min = ', cp%flg_min
    write( *, * )  'mueff = ', cp%mueff
    write( *, * )  'cm = ', cp%cm
    write( *, * )  'cc = ', cp%cc
    write( *, * )  'cone = ', cp%cone
    write( *, * )  'cmu = ', cp%cmu
    write( *, * )  'cs = ', cp%cs
    write( *, * )  'ds = ', cp%ds
    write( *, * )  'dble_dim = ', cp%dble_dim
    write( *, * )  'chidim = ', cp%chidim        
    write( *, * )  'weights = ', cp%weights
    write( *, * )  'neg_weights = ', cp%neg_weights
    ! Dynamic parameters
    write( *, * )  'iter = ', cp%iter
    write( *, * )  'neval = ', cp%neval
    write( *, * )  'eig_iter = ', cp%eig_iter
    write( *, * )  'sigm = ', cp%sigma
    write( *, * )  'xmean = ', cp%xmean
    write( *, * )  'vecD = ', cp%vecD
    write( *, * )  'ps = ', cp%ps
    write( *, * )  'pc = ', cp%pc
    write( *, * )  'diagSqrtC = ', cp%diagSqrtC
    if ( cp%flag_diagonal == 0 ) then
       write( *, * )  'matC = ', cp%matC
       write( *, * )  'matB = ', cp%matB
       write( *, * )  'matInvSqrtC = ', cp%matInvSqrtC
    end if
    ! Box Constraint
    write( *, * )  'idx_fhistory = ', cp%idx_fhistory
    write( *, * )  'flg_fhistory = ', cp%flg_fhistory
    write( *, * )  'flg_winit = ', cp%flg_winit
    write( *, * )  'cnt_winit = ', cp%cnt_winit
    write( *, * )  'ub = ', cp%ub
    write( *, * )  'lb = ', cp%lb
    write( *, * )  'penalty = ', cp%penalty
    write( *, * )  'fhistory = ', cp%fhistory
    ! Mixed-Integer Optimization
    write( *, * )  'flg_mixed_integer = ', cp%flg_mixed_integer
    write( *, * )  'sigmaCorrectionFactor = ', cp%sigmaCorrectionFactor
    write( *, * )  'vecStairWidth = ', cp%vecStairWidth
    write( *, * )  'bestsolutionX = ', cp%bestsolutionX
    ! Others
    write( *, * )  'flg_dsyev_err = ', cp%flg_dsyev_err
    if ( cp%flag_diagonal == 0 ) then
       write( *, * )  'lwork = ', cp%lwork
       write( *, * )  'work = ', cp%work
       write( *, * )  'workmat = ', cp%workmat
    end if

  end subroutine printcp

  Subroutine clear( cp )
    implicit none
    type(CmaParam), intent(inout) :: cp
    deallocate(cp%weights)
    deallocate(cp%xmean, cp%vecD, cp%ps, cp%pc, cp%diagSqrtC)
    deallocate(cp%ub, cp%lb, cp%penalty)
    deallocate(cp%fhistory)
    if (cp%flag_diagonal == 0) then    
       deallocate(cp%matC, cp%matB, cp%matInvSqrtC)
       deallocate(cp%work, cp%workmat)       
    end if
    if (cp%flag_active == 1) then
       deallocate(cp%neg_weights)
    endif
    deallocate(cp%vecStairWidth, cp%bestsolutionX)
  end Subroutine clear

  Subroutine initialize_cmaes( cp, num_param, lb, ub, init_mean, init_std, &
       stairwidth, flag_diagonal, popsize, flag_active, sigma_correction_factor )

    implicit none
    
    type(CmaParam), intent(inout) :: cp    
    integer, intent(in) :: num_param
    double precision, intent(in) :: lb(:)
    double precision, intent(in) :: ub(:)    
    double precision, intent(in), optional :: init_mean(:) ! default = (ub + lb) / 2
    double precision, intent(in), optional :: init_std(:)  ! default = (ub - lb) / 4
    double precision, intent(in), optional :: stairwidth(:) ! default = 0
    integer, intent(in), optional :: flag_diagonal ! 0: CMA (default), 1: Separable CMA
    integer, intent(in), optional :: popsize       ! default = 4 + |3 * ln(N)|
    integer, intent(in), optional :: flag_active   ! 0: CMA (default), 1: Active CMA
    double precision, intent(in), optional :: sigma_correction_factor ! default = 1.0d0

    integer :: i
    double precision :: neg_mueff, neg_amu

    ! Compute the default static parameters
    cp%dim = num_param
    cp%dble_dim = dble(cp%dim)
    if (present(flag_diagonal)) then
       cp%flag_diagonal = flag_diagonal
    else
       cp%flag_diagonal = 0 
    endif
    cp%lambda = 4 + int(3.0d0 * log(cp%dble_dim))    
    if (present(popsize)) then
       if  ( popsize < cp%lambda ) then
          print *, 'Warning:'
          print *, 'Initial Population Size "popsize" = ', popsize, '.'
          print *, 'It is smaller than the default value 4 + int(3 * log(dble(N))) = ', cp%lambda, '.'
          print *,  'popsize >= 4 + int(3 * log(dble(N))) recommented.' 
       end if
       cp%lambda = popsize
    endif    
    cp%mu = cp%lambda / 2
    if (present(flag_active)) then
       cp%flag_active = flag_active
    else
       cp%flag_active = 0
    endif
    cp%fhistory_size = 20+(3*cp%dim)/cp%lambda

    ! Memory Allocation
    allocate(cp%weights(cp%mu))
    if ( cp%flag_active == 1 ) then
       allocate(cp%neg_weights(cp%mu))
    endif
    allocate(cp%vecStairWidth(cp%dim), cp%bestsolutionX(cp%dim))
    allocate(cp%xmean(cp%dim), cp%ps(cp%dim), cp%pc(cp%dim), cp%vecD(cp%dim), cp%diagSqrtC(cp%dim))
    allocate(cp%lb(cp%dim), cp%ub(cp%dim), cp%penalty(cp%dim))
    allocate(cp%fhistory(cp%fhistory_size))
    if (cp%flag_diagonal == 0) then
       allocate(cp%matB(cp%dim, cp%dim), cp%matC(cp%dim, cp%dim), &
            cp%matInvSqrtC(cp%dim, cp%dim))
    endif

    cp%lb(:) = lb(:)
    cp%ub(:) = ub(:)
    cp%penalty(:) = 0.0d0
    if (present(init_mean)) then
       cp%xmean(:) = init_mean(:)
    else
       cp%xmean(:) = (ub(:) + lb(:)) / 2.0d0
    end if
    if (present(init_std)) then
       cp%vecD(:) = init_std(:)       
    else
       cp%vecD(:) = (ub(:) - lb(:)) / 4.0d0
    end if
    cp%diagSqrtC(:) = cp%vecD
    if (cp%flag_diagonal == 0) then
       cp%matB(:, :) = 0.0d0
       cp%matC(:, :) = 0.0d0       
       cp%matInvSqrtC(:, :) = 0.0d0
       do i = 1, cp%dim
          cp%matB(i, i) = 1.0d0
          cp%matC(i, i) = cp%vecD(i) * cp%vecD(i)
          cp%matInvSqrtC(i, i) = 1.0d0 / cp%vecD(i)
       end do
    endif
    cp%sigma = 1.0d0
    cp%pc(:) = 0.0d0
    cp%ps(:) = 0.0d0

    ! Compute the positive weights
    do i = 1, cp%mu
       cp%weights(i) = log(dble(cp%lambda + 1) / 2.0d0) - log(dble(i))
    end do
    cp%weights(:) = cp%weights(:) / sum(cp%weights)
    cp%mueff = sum(cp%weights)**2 / sum(cp%weights * cp%weights)

    ! Learning Rates
    cp%cm = 1.0d0
    cp%cc = (4.0d0 + cp%mueff / cp%dble_dim) / (cp%dble_dim + 4.0d0 + 2.0d0 * cp%mueff / cp%dble_dim)
    cp%cs = (cp%mueff + 2.0d0) / (cp%mueff + cp%dble_dim + 5.0d0)
    cp%ds = 1.0d0 + cp%cs + &
         2.0d0 * max(0.0d0, sqrt((cp%mueff - 1.0d0) / (cp%dble_dim + 1.0d0)) - 1.0d0)
    if (cp%flag_diagonal == 0) then
       cp%cone = 2.0d0 / ((cp%dble_dim + 1.3d0)**2 + cp%mueff)
       cp%cmu = min(1.0d0 - cp%cone, &
            2.0d0 * (cp%mueff - 2.0d0 + 1.0d0/cp%mueff) / ((cp%dble_dim + 2.0d0)**2 + cp%mueff))
    else if (cp%flag_diagonal == 1) then
       cp%cone = 2.0d0 / ((cp%dble_dim + 1.3d0)**2 + cp%mueff) * (cp%dble_dim + 2.0d0) / 3.0d0
       cp%cmu = min(1.0d0 - cp%cone, &
            2.0d0 * (cp%mueff - 2.0d0 + 1.0d0/cp%mueff) / ((cp%dble_dim + 2.0d0)**2 + cp%mueff) * &
            (cp%dble_dim + 2.0d0) / 3.0d0)
    end if

    ! Compute the negative weights
    if ( cp%flag_active == 1 ) then
       do i = 1, cp%mu
          cp%neg_weights(i) = log(dble(cp%lambda + 1) / 2.0d0) - log(dble(cp%lambda + 1 - i))
       end do
       cp%neg_weights(:) = cp%neg_weights(:) / sum(abs(cp%neg_weights))
       neg_mueff = sum(cp%neg_weights)**2 / sum(cp%neg_weights * cp%neg_weights)
       neg_amu = 1.0d0 + min(cp%cone / cp%cmu, 2.0d0 * neg_mueff / (neg_mueff + 2.0d0))
       cp%neg_weights(:) = cp%neg_weights(:) * min(neg_amu, (1.0d0 - cp%cone - cp%cmu) / (cp%dble_dim * cp%cmu))
    endif
    
    ! Set parameters of mixed-integer optimization    
    if ( present(sigma_correction_factor) ) then
       cp%sigmaCorrectionFactor = sigma_correction_factor
    else
       cp%sigmaCorrectionFactor = 0.0d0
    end if
    cp%vecStairWidth(:) = stairwidth(:)
    cp%bestsolutionX(:) = cp%xmean(:)
    cp%flg_mixed_integer = 0
    do i = 1, cp%dim
       if ( cp%vecStairWidth(i) .ne. 0.0d0 ) then
          cp%flg_mixed_integer = 1
          exit
       end if
    end do
    
    ! Others
    cp%iter = 1
    cp%neval = 0
    cp%flg_winit = 0
    cp%cnt_winit = 3
    cp%idx_fhistory = 1
    cp%flg_fhistory = 0
    cp%eig_iter = 0
    cp%chidim = sqrt(cp%dble_dim) * &
         (1.0d0 - 0.25d0 / cp%dble_dim + 1.0d0 / (21.0d0 * cp%dble_dim * cp%dble_dim))
    cp%flg_dsyev_err = 0

    ! Work Space
    if (cp%flag_diagonal == 0) then    
       cp%lwork = 3 * cp%dim - 1
       allocate(cp%work(cp%lwork))
       allocate(cp%workmat(cp%dim, max(cp%dim, cp%lambda)))       
    end if

  end subroutine initialize_cmaes
  
  subroutine randn( outmat )
    implicit none
    double precision, intent(out) :: outmat(:, :)

    integer :: i, j, m, n

    m = size(outmat, 1)
    n = size(outmat, 2)    
    do j = 1, n
       do i = 1, m
          outmat(i, j) = normal(0.0d0, 1.0d0)
       end do
    end do
  end subroutine randn

  subroutine integer_mutation( cp, inmat, outmat )
    implicit none
    
    type(CmaParam), intent(inout) :: cp
    double precision, intent(in) :: inmat(:, :)
    double precision, intent(out) :: outmat(:, :)    
    
    integer :: i, j, m, n

    ! Mixed-Integer optimization parameters
    integer :: IR(cp%dim), randdim(cp%dim), idxdim(cp%dim)
    integer, allocatable :: Iplusminus(:, :)
    integer, allocatable :: Rint(:, :), Rp(:, :), Rpp(:, :)
    integer, allocatable :: rndidx(:)
    integer :: absIR, lamint
    
    m = size(inmat, 1)
    n = size(inmat, 2)

    ! for mixed-integer optimization
    ! following steps prepare for Rint.
    allocate(Rint(m, n))

    ! compute IR, it is masking vector for the components, where an integer mutation will be conducted
    absIR = 0
    do i = 1, m
       if ( 2.0d0 * coordinate_length(cp, i) < cp%vecStairWidth(i) ) then
          IR(i) = 1
          absIR = absIR + 1
          idxdim(absIR) = i
       else
          IR(i) = 0
       end if
    end do

    ! compute lambda_int, it is number of solutions that be applied an integer mutation
    if ( absIR == 0 ) then
       lamint = 0
    else if ( absIR == m ) then
       lamint = n / 2
    else
       lamint = int(min( dble(n)/1.0d1 + dble(absIR) + 1.0d0, dble(n/2 - 1)))
    end if

    if ( lamint == 0 ) then
       Rint(:, :) = 0
    else
       allocate(Rp(m, lamint), Rpp(m, lamint), Iplusminus(m, lamint))

       ! compute Rp
       randdim = randperm(absIR)          
       Rp(1:m, 1:lamint) = 0
       do i = 1, lamint
          Rp(idxdim(randdim(mod(i - 1, absIR) + 1)), i) = 1
       end do

       ! compute Rpp
       Rpp(:, :) = 0
       do j = 1, lamint
          do i = 1, m
             if ( IR(i) == 1 ) then
                do while ( uniform(0.0d0, 1.0d0) >= 0.7d0 ** (1.0d0 / dble(absIR)) )
                   Rpp(i, j) = Rpp(i, j) + 1
                end do
             end if
          end do
       end do

       ! compute Iplusminus, it is sign-switching (dim, lamint)-array for Rint, with each element being +-1
       do j = 1, lamint
          do i = 1, m
             if (uniformint(0, 1) == 0) then
                Iplusminus(i, j) = -1
             else
                Iplusminus(i, j) = 1
             end if
          end do
       end do

       ! compute Rint, it is the integer mutation
       ! type: int
       Rint(1:m, 1:lamint) = Iplusminus(1:m, 1:lamint) * (Rp(1:m, 1:lamint) + Rpp(1:m, 1:lamint))
       do i = 1, cp%dim
          if ( cp%vecStairWidth(i) == 0.0d0 ) then
             Rint(i, lamint+1) = 0
          else
             Rint(i, lamint+1) = floor(cp%bestsolutionX(i) / cp%vecStairWidth(i)) - &
                  floor(cp%xmean(i) / cp%vecStairWidth(i))
          end if
       end do
       Rint(1:m, lamint+2:n) = 0
    end if

    ! apply integer mutation
    allocate(rndidx(n))
    rndidx = randperm(n)
    do i = 1, n
       outmat(1:m, rndidx(i)) = inmat(1:m, rndidx(i)) + (cp%vecStairWidth(1:m) * dble(Rint(1:m, i)))
    end do
  end subroutine integer_mutation
  
  subroutine generate_candidate_solutions( cp, outmat )
    implicit none
    
    type(CmaParam), intent(inout) :: cp
    double precision, intent(out) :: outmat(:, :)

    integer :: i, m, n

    m = size(outmat, 1)
    n = size(outmat, 2)
    
    call randn(outmat)
    if (cp%flag_diagonal == 0) then
       
       do i = 1, n
          cp%workmat(1:m, i) = outmat(1:m, i) * cp%vecD(1:m)
          outmat(1:m, i) = cp%xmean(1:m)
       end do
       call mydgemm( 'N', 'N', m, n, m, cp%sigma, cp%matB(1:m, 1:m), m, &
            cp%workmat(1:m, 1:n), m, 1.0d0, outmat(1:m, 1:n), m)
    else if (cp%flag_diagonal == 1) then
       do i = 1, n
          outmat(1:m, i) = cp%xmean(1:m) + cp%sigma * cp%vecD(1:m) * outmat(1:m, i)
       end do
    end if

  end subroutine generate_candidate_solutions

  subroutine update_parameters( cp, arx, arxmod, arf )
    implicit none

    type(CmaParam), intent(inout) :: cp
    double precision, intent(in) :: arx(:, :)
    double precision, intent(in) :: arxmod(:, :)    
    double precision, intent(in) :: arf(:)
    
    integer :: i, lam, dim
    integer :: idx(cp%lambda)
    double precision :: dx(cp%dim)
    double precision :: hsig, normps, sum_weights, scaled_nw
    double precision :: ymean(cp%dim)
    double precision :: sary(cp%dim, cp%lambda)

    ! Mixed-Integer optimization parameters
    integer :: n_int_val
    double precision :: sig_fac
    double precision :: chidim_mixed_integer, count_Isig
    
    dim = cp%dim
    lam = cp%lambda
    call indexsort(lam, arf, idx)
    do i = 1, lam
       sary(1:dim, i) = (arx(1:dim, idx(i)) - cp%xmean) / cp%sigma
    end do
    ! Logging best solution for Mixed-integer optimization
    cp%bestsolutionX(1:dim) = arxmod(1:dim, idx(1))
    ymean(:) = 0.0d0  ! fix: 2017/10/03, thanks to Mr. Banno@OsakaUniv.
    call mydgemv('N', dim, cp%mu,1.0d0, sary(1:dim, 1:lam), dim, &
         cp%weights(:), 1, 0.0d0, ymean(1:dim), 1)
    if (cp%flag_diagonal == 0) then
       call mydgemv('N', dim, dim, sqrt(cp%cs * (2.0d0 - cp%cs) * cp%mueff), &
            cp%matInvSqrtC(1:dim, 1:dim), dim, ymean(1:dim), 1, (1.0d0 - cp%cs), &
            cp%ps(1:dim), 1)
       ! cp%ps = (1.0d0 - cp%cs) * cp%ps + &
       !      sqrt(cp%cs * (2.0d0 - cp%cs) * cp%mueff) * matmul(cp%matInvSqrtC, ymean)
    else if (cp%flag_diagonal == 1) then
       cp%ps(1:dim) = (1.0d0 - cp%cs) * cp%ps(1:dim) + &
            sqrt(cp%cs * (2.0d0 - cp%cs) * cp%mueff) * (ymean(1:dim) / cp%vecD(1:dim))
    end if

    normps = sqrt(myddot(cp%ps(1:dim), cp%ps(1:dim)))
    if (normps < (1.5d0 + 1.0d0 / (cp%dble_dim - 0.5d0)) * cp%chidim) then
       hsig = 1.0d0  !TODO:
    else
       hsig = 0.0d0
    end if
    cp%pc(1:dim) = (1.0d0 - cp%cc) * cp%pc(1:dim) + &
         hsig * sqrt(cp%cc * (2.0d0 - cp%cc) * cp%mueff) * ymean(1:dim)

    ! Modified by YA on Sep 13, 2017
    dx(1:dim) = 0.0d0
    do i = 1, cp%mu
       dx(1:dim) = dx(1:dim) + cp%weights(i) * (arxmod(1:dim, idx(i)) - cp%xmean)
    end do
    cp%xmean = cp%xmean + cp%cm * dx

    if ( cp%flg_mixed_integer == 1 ) then
       normps = 0.0d0
       count_Isig = 0.0d0
       do i = 1, dim
          if ( coordinate_length(cp, i) / sqrt(cp%cs) >= 0.2d0 * cp%vecStairWidth(i) ) then
             normps = normps + cp%ps(i) ** 2
             count_Isig = count_Isig + 1.0d0
          end if
       end do
       normps = sqrt(normps)
       
       if ( int(count_Isig) > 0 ) then
          chidim_mixed_integer = sqrt(count_Isig) * &
               (1.0d0 - 0.25d0 / count_Isig + 1.0d0 / (21.0d0 * count_Isig * count_Isig))
          cp%sigma = cp%sigma * exp((cp%cs / cp%ds) * (normps / chidim_mixed_integer - 1.0d0))
       end if
       ! sigma correction for mixed integer
       n_int_val = 0
       sig_fac = 1.0d0
       do i = 1, cp%dim
          if ( cp%vecStairWidth(i) > 0.0d0 ) then
             n_int_val = n_int_val + 1
          end if
       end do
       do i = 1, cp%dim
          sig_fac = max( sig_fac, cp%sigmaCorrectionFactor * &
               cp%vecStairWidth(i) / coordinate_length(cp, i) / dble(2 * n_int_val + 1))
       end do
       cp%sigma = cp%sigma * sig_fac
    else
       cp%sigma = cp%sigma * exp((cp%cs / cp%ds) * (normps / cp%chidim - 1.0d0))
    end if

    ! CMA
    if ( cp%flag_active == 1 ) then
       sum_weights = sum(cp%weights(:)) + sum(cp%neg_weights(:))
    else if ( cp%flag_active == 0 ) then
       sum_weights = sum(cp%weights(:))
    endif

    ! Full CMA
    if (cp%flag_diagonal == 0) then
       cp%matC(1:dim, 1:dim) = &
            (1.0d0 - hsig * cp%cone - cp%cmu * sum_weights) * cp%matC(1:dim, 1:dim)
       call mydsyr('U', dim, (hsig * cp%cone), cp%pc(1:dim), 1, &
            cp%matC(1:dim, 1:dim), dim)
       do i = 1, cp%mu

          ! Negative Update          
          if ( cp%flag_active == 1 ) then
             ! The following initialization is necessary to prevent 0 * nan
             cp%work(1:dim) = 0.0d0  ! BUG fixed on Jan 18, 2018
             call mydgemv('N', dim, dim, 1.0d0, cp%matInvSqrtC(1:dim, 1:dim), dim, &
                  sary(1:dim, lam + 1 - i), 1,0.0d0, cp%work(1:dim), 1)
             scaled_nw = cp%neg_weights(i) * cp%dble_dim &
                  / myddot(cp%work(1:dim), cp%work(1:dim))
             !scaled_nw = cp%neg_weights(i) * cp%dble_dim / sum(matmul(cp%matInvSqrtC, sary(:, cp%lambda + 1 - i)) ** 2)
             call mydsyr('U', dim, (cp%cmu * scaled_nw), sary(:, lam + 1 - i), 1, &
                  cp%matC(1:dim, 1:dim), dim)                    
          endif

          ! Positive Update
          call mydsyr('U', dim, (cp%cmu * cp%weights(i)), sary(1:dim, i), 1, &
               cp%matC(1:dim, 1:dim), dim)
       end do

       ! Perform Eigendecomposition every after 1 / (cone + cmu) / (10 * N) iterations
       cp%eig_iter = cp%eig_iter + 1
       if (dble(cp%eig_iter * cp%dim * 10) * (cp%cone + cp%cmu) > 1.0d0) then
          cp%eig_iter = 0  ! initialize the counter

          ! Eigendecomposition of matC.
          cp%matB(1:dim, 1:dim) = cp%matC(1:dim, 1:dim)          
          ! Note that matC is replaced with its eigenvector matrix.
          call mydsyev( 'V', 'U', dim, cp%matB(1:dim, 1:dim), dim, cp%vecD(1:dim), &
               cp%work(1:cp%lwork), cp%lwork, cp%flg_dsyev_err)
          ! (Re)compute C = B * D * B'
          do i = 1, dim
             cp%workmat(1:dim, i) = cp%matB(1:dim, i) * cp%vecD(i)
          end do
          call mydgemm('N', 'T', dim, dim, dim, 1.0d0, cp%matB(1:dim, 1:dim), dim, &
               cp%workmat(1:dim, 1:dim), dim, 0.0d0, cp%matC(1:dim, 1:dim), dim)
          ! D <- sqrt of eigenvalues
          cp%vecD(1:dim) = sqrt(cp%vecD(1:dim))

          ! Compute C^{-1/2}
          do i = 1, dim
             cp%workmat(1:dim, i) = cp%matB(1:dim, i) / cp%vecD(i)
          end do
          call mydgemm('N', 'T', dim, dim, dim, 1.0d0, cp%matB(1:dim, 1:dim), dim, &
               cp%workmat(1:dim, 1:dim), dim, 0.0d0, cp%matInvSqrtC(1:dim, 1:dim), dim)
          ! diag(C)^{1/2} for stopping criteria and visualization
          do i = 1, dim
             cp%diagSqrtC(i) = sqrt(cp%matC(i, i))
          end do
       end if

    ! Separable CMA
    else if (cp%flag_diagonal == 1) then
       cp%vecD(:) = (1.0d0 - hsig * cp%cone - cp%cmu * sum_weights) * cp%vecD(:) * cp%vecD(:)
       cp%vecD(:) = cp%vecD(:) + (hsig * cp%cone) * cp%pc(:) * cp%pc(:)
       do i = 1, cp%mu
          ! Negative Update
          if ( cp%flag_active == 1 ) then
             scaled_nw = cp%neg_weights(i) * cp%dble_dim / sum((sary(:, cp%lambda + 1 - i) / cp%diagSqrtC(:)) ** 2) 
             cp%vecD(:) = cp%vecD(:) + cp%cmu * scaled_nw * sary(:, cp%lambda + 1 - i) * sary(:, cp%lambda + 1 - i)
          endif
          ! Positive Update
          cp%vecD(:) = cp%vecD(:) + cp%cmu * cp%weights(i) * sary(:, i) * sary(:, i)
       end do

       ! D <- sqrt(D)
       cp%vecD(:) = sqrt(cp%vecD)
       cp%diagSqrtC(:) = cp%vecD(:)
    end if

  end subroutine update_parameters

  subroutine handle_box_boundary(cp, arx, arxfeas)
    implicit none
    !! argument declaration 
    type(CmaParam), intent(in) :: cp    
    double precision, intent(in) :: arx(:, :)
    double precision, intent(out) :: arxfeas(:, :)
    !! variable declaration
    integer :: i, j 

    arxfeas(:, :) = arx(:, :)
    do j = 1, cp%lambda
       do i = 1, cp%dim
          if (arx(i, j) < cp%lb(i)) then
             arxfeas(i, j) = cp%lb(i)
          else if (arx(i, j) > cp%ub(i)) then
             arxfeas(i, j) = cp%ub(i)
          end if
       end do
    end do
  end subroutine handle_box_boundary
  
  ! subroutine compute_penalty(cp, arxfeas, arffeas, arx, arf)
  !   !! Box Constraint Handling with Adaptive Penalty Function
  !   !! Author Version of "A Method for Handling Uncertainty ..." by Niko (2009)
  !   implicit none
  !   !! argument declaration 
  !   type(CmaParam), intent(inout) :: cp
  !   double precision, intent(in) :: arxfeas(:, :)    
  !   double precision, intent(in) :: arffeas(:)    
  !   double precision, intent(in) :: arx(:, :)
  !   double precision, intent(out) :: arf(:)    
  !   !! variable declaration
  !   integer :: i, j, k, flg
  !   double precision :: delta, tmp, delth, dgam, delm
  !   integer :: fhidx(cp%fhistory_size), fidx(cp%lambda)
  !   double precision :: xi(cp%dim)

  !   ! Modified Penalty Box Constraint (based on the author's version of Hansen et al. 2009)
  !   cp%cnt_winit = max(cp%cnt_winit - 1, 0)
    
  !   ! Compute the trace of C
  !   tmp = trace( cp )
        
  !   ! Update fhistory
  !   call indexsort(cp%lambda, arffeas, fidx)    
  !   !! CHECK POINT I: Normalized IQR
  !   cp%fhistory(cp%idx_fhistory) = &
  !        (arffeas(fidx(cp%lambda*3/4)) - arffeas(fidx(cp%lambda/4))) / &
  !        (cp%sigma * cp%sigma * tmp / cp%dble_dim)
    
  !   cp%idx_fhistory = cp%idx_fhistory + 1
  !   if (cp%idx_fhistory > cp%fhistory_size) then
  !      cp%flg_fhistory = 1
  !      cp%idx_fhistory = 1 !! BUG: fixed on Dec. 7, 2015, (thanks to Mr. Takahashi)
  !   end if

  !   ! Check the feasibility of xmean
  !   flg = 0
  !   do i = 1, cp%dim
  !      if (cp%xmean(i) < cp%lb(i)) then
  !         flg = 1
  !      else if (cp%xmean(i) > cp%ub(i)) then
  !         flg = 1
  !      end if
  !   end do

  !   ! Set the penalty coefficient            
  !   if (flg == 1) then
  !      if (cp%flg_fhistory == 1) then
  !         k = size(cp%fhistory)
  !      else
  !         k = cp%idx_fhistory - 1
  !      end if

  !      !! CHECK POINT II: delta computed only once or twice       
  !      if ((cp%flg_winit == 0) .or. (cp%cnt_winit > 0)) then
  !         call indexsort(k, cp%fhistory(1:k), fhidx(1:k))
  !         delta = cp%fhistory(fhidx(k/2))

  !         !! CHECK POINT I: Normalized IQR
  !         cp%penalty(:) = 2.0d0 * delta          
  !         cp%flg_winit = 1
  !      end if
  !   end if

  !   !! CHECK POINT III: different update
  !   ! Increase and decrease the penalty coefficients
  !   delth = 3.0d0 * max(1.0d0, sqrt(cp%dble_dim) / cp%mueff)
  !   dgam = min(1.0d0, cp%mueff / (1.0d1 * cp%dble_dim))
  !   do i = 1, cp%dim
  !      delm = max(cp%lb(i) - cp%xmean(i), 0.0d0) + max(cp%xmean(i) - cp%ub(i), 0.0d0)
  !      delm = delm / coordinate_length(cp, i) ! BUG: fixed on Dec. 9, 2016
  !      cp%penalty(i) = cp%penalty(i) * exp((dgam / 2.0d0) * tanh(max(0.0d0, delm - delth) / 3.0d0))
  !      if (cp%penalty(i) > 5.0d0 * delta) then
  !         cp%penalty(i) = cp%penalty(i) * exp(- dgam / 3.0d0)
  !      end if
  !   end do

  !   ! Compute the penalized fitness
  !   ! if (cp%flag_diagonal == 0) then 
  !   !    do i = 1, cp%dim
  !   !       xi(i) = log(cp%matC(i, i))
  !   !    end do
  !   ! else if (cp%flag_diagonal == 1) then
  !   !    do i = 1, cp%dim
  !   !       xi(i) = 2.0d0 * log(cp%vecD(i))
  !   !    end do
  !   ! end if
  !   do i = 1, cp%dim
  !      xi(i) = 2.0d0 * log(coordinate_length(cp, i))
  !   end do
  !   xi(:) = exp(0.9d0 * (xi(:) - sum(xi(:))/cp%dble_dim))
       
  !   arf(:) = 0.0d0
  !   do j = 1, cp%lambda
  !      do i = 1, cp%dim
  !         arf(j) = arf(j) + (arxfeas(i, j) - arx(i, j))**2 * cp%penalty(i) / xi(i)
  !      end do
  !   end do
  !   arf(:) = arffeas(:) + arf(:) / cp%dble_dim
  ! end subroutine compute_penalty


  ! subroutine compute_penalty_v1(cp, arxfeas, arffeas, arx, arf)
  !   !! Box Constraint Handling with Adaptive Penalty Function
  !   !! Based on "A Method for Handling Uncertainty ..." by Niko (2009),
  !   !! but the penalty coefficients are updated every iteration using IQR history.
  !   implicit none

  !   type(CmaParam), intent(inout) :: cp
  !   double precision, intent(in) :: arxfeas(:, :)    
  !   double precision, intent(in) :: arffeas(:)    
  !   double precision, intent(in) :: arx(:, :)
  !   double precision, intent(out) :: arf(:)    

  !   integer :: i, j, k, flg
  !   double precision :: delta, tmp
  !   integer :: fhidx(size(cp%fhistory)), fidx(cp%lambda)
  !   double precision :: xi(cp%dim)

  !   ! Update fhistory
  !   call indexsort(cp%lambda, arffeas, fidx)

  !   !! CHECK POINT I
  !   cp%fhistory(cp%idx_fhistory) = arffeas(fidx(cp%lambda*3/4)) - arffeas(fidx(cp%lambda/4))

  !   cp%idx_fhistory = cp%idx_fhistory + 1
  !   if (cp%idx_fhistory > size(cp%fhistory)) then
  !      cp%flg_fhistory = 1
  !      !cp%idx_fhistory = 0 !! BUG: fixed on Dec. 7, 2015, (thanks to Mr. Takahasi)
  !      cp%idx_fhistory = 1 
  !   end if

  !   flg = 0
  !   ! Set the penalty coefficient    
  !   do i = 1, cp%dim
  !      if (cp%xmean(i) < cp%lb(i)) then
  !         flg = 1
  !      else if (cp%xmean(i) > cp%ub(i)) then
  !         flg = 1
  !      end if
  !   end do
    
  !   if (flg == 1) then
  !      if (cp%flg_fhistory == 1) then
  !         k = size(cp%fhistory)
  !      else
  !         k = cp%idx_fhistory - 1
  !      end if

  !      !! CHECK POINT II
  !      call indexsort(k, cp%fhistory(1:k), fhidx(1:k))
  !      delta = cp%fhistory(fhidx(k/2))

  !      tmp = trace( cp )
  !      ! tmp = 0.0d0
  !      ! if (cp%flag_diagonal == 0) then
  !      !    do i = 1, cp%dim
  !      !       tmp = tmp + cp%matC(i, i)
  !      !    end do
  !      ! else if (cp%flag_diagonal == 1) then
  !      !    do i = 1, cp%dim
  !      !       tmp = tmp + cp%vecD(i) * cp%vecD(i)
  !      !    end do
  !      ! end if

  !      !! CHECK POINT I
  !      cp%penalty(:) = 2.0d0 * delta / (cp%sigma * cp%sigma * tmp / cp%dble_dim)
  !   end if
    
  !   !! CHECK POINT III: different update
  !   ! Increase the penalty coefficient
  !   do i = 1, cp%dim
  !      if ((cp%xmean(i) < cp%lb(i)) .or. (cp%xmean(i) > cp%ub(i))) then
  !         if (max(cp%lb(i) - cp%xmean(i), cp%xmean(i) - cp%ub(i)) > &
  !              3.0d0 * coordinate_length(cp, i) * &
  !              max(1.0d0, sqrt(cp%dble_dim)/cp%mueff)) then
  !            cp%penalty(i) = cp%penalty(i) * 1.1d0**min(1.0d0, cp%mueff/dble(10*cp%dim))
  !         end if
  !      end if
  !   end do    
  !   ! do i = 1, cp%dim
  !   !    if ((cp%xmean(i) < cp%lb(i)) .or. (cp%xmean(i) > cp%ub(i))) then
  !   !       if (cp%flag_diagonal == 0) then
  !   !          if (max(cp%lb(i) - cp%xmean(i), cp%xmean(i) - cp%ub(i)) > &
  !   !               3.0d0 * cp%sigma * sqrt(cp%matC(i, i)) * &
  !   !               max(1.0d0, sqrt(cp%dble_dim)/cp%mueff)) then
  !   !             cp%penalty(i) = cp%penalty(i) * 1.1d0**max(1.0d0, cp%mueff/dble(10*cp%dim))
  !   !          end if
  !   !       else if (cp%flag_diagonal == 1) then
  !   !          if (max(cp%lb(i) - cp%xmean(i), cp%xmean(i) - cp%ub(i)) > &
  !   !               3.0d0 * cp%sigma * cp%vecD(i) * &
  !   !               max(1.0d0, sqrt(cp%dble_dim)/cp%mueff)) then
  !   !             cp%penalty(i) = cp%penalty(i) * 1.1d0**max(1.0d0, cp%mueff/dble(10*cp%dim))
  !   !          end if
  !   !       end if
  !   !    end if
  !   ! end do

  !   ! Compute the penalized fitness
  !   ! if (cp%flag_diagonal == 0) then 
  !   !    do i = 1, cp%dim
  !   !       xi(i) = log(cp%matC(i, i))
  !   !    end do
  !   ! else if (cp%flag_diagonal == 1) then
  !   !    do i = 1, cp%dim
  !   !       xi(i) = 2.0d0 * log(cp%vecD(i))
  !   !    end do
  !   ! end if
  !   do i = 1, cp%dim
  !      xi(i) = 2.0d0 * log(coordinate_length(cp, i))
  !   end do    
  !   xi(:) = exp(0.9d0 * (xi(:) - sum(xi(:))/cp%dble_dim))
       
  !   arf(:) = 0.0d0
  !   do j = 1, cp%lambda
  !      do i = 1, cp%dim
  !         arf(j) = arf(j) + (arxfeas(i, j) - arx(i, j))**2 * cp%penalty(i) / xi(i)
  !      end do
  !   end do
  !   arf(:) = arffeas(:) + arf(:) / cp%dble_dim
  ! end subroutine compute_penalty_v1

  subroutine compute_penalty(cp, arxfeas, arffeas, arx, arf)
    !! Box Constraint Handling with Adaptive Penalty Function [Sakamoto GECCO 17]
    !! Based on the author version of "A Method for Handling Uncertainty ..." by Niko (2009)
    !! Two main changes are:
    !! 1. median computation of the IQR history is replaced with the effective median
    !! 2. penalty coefficient is updated if their mean value is greater than three times
    !!    the effective median of the IQR.
    !! These modification improves the performance when the objective function is far from
    !! quadratic function such as an exponential function.
    implicit none
    !! argument declaration 
    type(CmaParam), intent(inout) :: cp
    double precision, intent(in) :: arxfeas(:, :)    
    double precision, intent(in) :: arffeas(:)    
    double precision, intent(in) :: arx(:, :)
    double precision, intent(out) :: arf(:)    
    !! variable declaration
    integer :: i, j, k, flg
    double precision :: delta, tmp, delth, dgam, delm, iqr_arffeas
    integer :: fidx(cp%lambda) !, fhidx(cp%fhistory_size)

    cp%cnt_winit = max(cp%cnt_winit - 1, 0)
    
    ! Compute the trace of C
    tmp = trace( cp )
        
    ! Update fhistory
    call indexsort(cp%lambda, arffeas, fidx)
    iqr_arffeas = arffeas(fidx(cp%lambda*3/4)) - arffeas(fidx(cp%lambda/4))
    !print *, 'check A'
    call append_iqr_hist( cp, iqr_arffeas / (cp%sigma * cp%sigma * tmp / cp%dble_dim) )
    ! cp%fhistory(cp%idx_fhistory) = iqr_arffeas / (cp%sigma * cp%sigma * tmp / cp%dble_dim)
    ! cp%idx_fhistory = cp%idx_fhistory + 1
    ! if (cp%idx_fhistory > size(cp%fhistory)) then
    !    cp%flg_fhistory = 1
    !    cp%idx_fhistory = 1 !! BUG: fixed on Dec. 7, 2015, (thanks to Mr. Takahashi)
    ! end if

    !print *, 'check B'
    ! Compute the effective median of IQR history
    delta = effective_median_iqr( cp )

    !print *, 'check C'    
    ! Check the feasibility of xmean
    flg = 0
    do i = 1, cp%dim
       if (cp%xmean(i) < cp%lb(i)) then
          flg = 1
       else if (cp%xmean(i) > cp%ub(i)) then
          flg = 1
       end if
    end do

    !print *, 'check D'
    ! Set the penalty coefficient            
    if (flg == 1) then
       k = get_iqr_hist_size( cp )

       if ((cp%flg_winit == 0) .or. (cp%cnt_winit > 0)) then
          ! fhidx(1:k) = indexsort(cp%fhistory(1:k))
          ! delta = cp%fhistory(fhidx(k/2))
          cp%penalty(:) = 2.0d0 * delta          
          cp%flg_winit = 1
       end if
    end if

    !print *, 'check E'
    ! Increase the penalty coefficients
    delth = 3.0d0 * max(1.0d0, sqrt(cp%dble_dim) / cp%mueff)
    dgam = min(1.0d0, cp%mueff / (1.0d1 * cp%dble_dim))
    do i = 1, cp%dim
       delm = max(cp%lb(i) - cp%xmean(i), 0.0d0) + max(cp%xmean(i) - cp%ub(i), 0.0d0)
       delm = delm / coordinate_length(cp, i) ! BUG: fixed on Dec. 9, 2016
       cp%penalty(i) = cp%penalty(i) * exp((dgam / 2.0d0) * tanh(max(0.0d0, delm - delth) / 3.0d0))
    end do

    !print *, 'check F'
    ! Decrease the penalty coefficients
    tmp = sum( cp%penalty(:) ) / cp%dble_dim
    if ( tmp > 3.0d0 * delta ) then
       cp%penalty(:) = cp%penalty(:) * (3.0d0 * delta / tmp)
    end if

    ! Compute the penalty factor
    arf(:) = 0.0d0
    do j = 1, cp%lambda
       do i = 1, cp%dim
          arf(j) = arf(j) + (arxfeas(i, j) - arx(i, j))**2 * cp%penalty(i)
       end do
    end do
    arf(:) = arffeas(:) + arf(:) / cp%dble_dim
  end subroutine compute_penalty

  function effective_median_iqr( cp ) result(emi)
    implicit none
    type(CmaParam), intent(inout) :: cp
    double precision :: emi
    ! variable declaration
    integer :: i, t, k
    double precision :: med, tmp1, tmp2
    double precision :: iqrhist(cp%fhistory_size)
    integer :: idx(cp%fhistory_size)

    !print *, 'check I'
    ! step 0.
    k = get_iqr_hist_size( cp )    
    if ( k == 1 ) then
       emi = get_iqr_hist( cp, 1 )
       return
    else if ( k == 2 ) then
       tmp1 = get_iqr_hist( cp, 1 )
       tmp2 = get_iqr_hist( cp, 2 )          
       emi = (tmp1 + tmp2) / 2.0d0
       return
    end if

    !print *, 'check II'    
    ! step 1. compute the median of latest three iterations
    iqrhist(1) = get_iqr_hist( cp, 1 )
    iqrhist(2) = get_iqr_hist( cp, 2 )
    iqrhist(3) = get_iqr_hist( cp, 3 )    
    if ( iqrhist(1) < iqrhist(2) ) then
       if ( iqrhist(3) < iqrhist(1) ) then
          med = iqrhist(1)
       else if ( iqrhist(2) < iqrhist(3) ) then
          med = iqrhist(2)
       else
          med = iqrhist(3)          
       end if
    else if ( iqrhist(2) < iqrhist(1) ) then
       if ( iqrhist(1) < iqrhist(3) ) then
          med = iqrhist(1)
       else if ( iqrhist(3) < iqrhist(2) ) then
          med = iqrhist(2)
       else
          med = iqrhist(3)                    
       end if
    end if

    !print *, 'check III'
    ! step 2. compute effective iterations
    if ( k == 3 ) then
       emi = med
    else 
       tmp1 = log(med)
       t = k !! BUG: fixed on Dec. 26, 2016, (thanks to Prof. Fujii)
       do i = 4, k
          iqrhist(i) = get_iqr_hist( cp, i )
          tmp2 = log(iqrhist(i))
          if ( abs(tmp2 - tmp1) > log(5.0d0) ) then
             t = i - 1
             exit
          end if
       end do

       ! step 3. comupte
       call indexsort( t, iqrhist(1:t), idx(1:t) )
       !print *, 'check V'
       !print *, size(idx), t
       if (mod(t, 2) == 1) then
          !print *, size(idx), t          
          emi = iqrhist( idx((t + 1) / 2) )
       else
          !print *, iqrhist( idx(1:t) ), t
          emi = ( iqrhist( idx(t / 2) ) + iqrhist( idx(t / 2 + 1) ) ) / 2.0d0
          ! Added "/ 2.0d0" to compute mean value by Sakamoto 20170911
       end if
    end if
    !print *, 'check V'
  end function effective_median_iqr

  subroutine append_iqr_hist( cp, iqr )
    implicit none
    type(CmaParam), intent(inout) :: cp
    double precision, intent(in) :: iqr

    cp%fhistory(cp%idx_fhistory) = iqr
    ! update the next index
    cp%idx_fhistory = cp%idx_fhistory + 1
    if (cp%idx_fhistory > cp%fhistory_size) then
       cp%flg_fhistory = 1
       cp%idx_fhistory = 1 !! BUG: fixed on Dec. 7, 2015, (thanks to Mr. Takahashi)
    end if
  end subroutine append_iqr_hist

  function get_iqr_hist_size( cp ) result(histsize)
    implicit none
    type(CmaParam), intent(in) :: cp
    integer :: histsize
    
    if (cp%flg_fhistory == 1) then
       histsize = cp%fhistory_size
    else
       histsize = cp%idx_fhistory - 1
    end if
  end function get_iqr_hist_size
    
  function get_iqr_hist( cp, idx ) result(iqr)
    ! idx >= 1: return the IQR of the (t + 1 - idx)th iteration
    implicit none
    type(CmaParam), intent(in) :: cp
    integer, intent(in) :: idx
    double precision :: iqr
    integer :: i, k

    k = get_iqr_hist_size( cp )
    i = cp%idx_fhistory - idx
    
    if ( idx < 1 .or. idx > k ) then
       ! should not occur
       ! TODO: error
       print *, 'error in get_iqr_hist', idx, k
       iqr = -1.0d0
    else if ( i >= 1 ) then
       iqr = cp%fhistory( i )
    else
       iqr = cp%fhistory( k + i )
    end if

    return
  end function get_iqr_hist
    
  
  function check_convergence( cp, tolx, tolf ) result(flg)
    implicit none

    type(CmaParam), intent(inout) :: cp
    double precision, intent(in) :: tolx, tolf
    logical :: flg

    integer :: i, k
    double precision :: fdiff
    double precision :: dx(cp%dim)

    flg = .false.
    if (cp%flg_dsyev_err .ne. 0) then
       flg = .true.
       print *, 'dsyev_err_', cp%flg_dsyev_err
    end if
    ! check tolf
    if (cp%flg_fhistory == 1) then
       fdiff = maxval(cp%fhistory(:)) * cp%sigma * cp%sigma * trace(cp) / cp%dble_dim
       if (fdiff <= tolf) then
          flg = .true.
          print *, 'tolf'
       end if
    end if
    ! if (cp%flg_fhistory == 1) then
    !    k = size(cp%fhistory)
    ! else
    !    k = cp%idx_fhistory - 1
    ! end if
    ! fdiff = maxval(cp%fhistory(1:k)) * cp%sigma * cp%sigma * trace(cp) / cp%dble_dim
    ! if (fdiff <= tolf) then
    !    flg = .true.
    !    print *, 'tolf'
    ! end if
    ! check tolx
    do i = 1, cp%dim
       dx(i) = coordinate_length(cp, i)
    enddo
    if (maxval(dx(:)) <= tolx) then
       flg = .true.
       print *, 'tolx'       
    end if
    if (minval(cp%vecD) <= 1.0d-7) then
       flg = .true.
       print *, 'min eigen val'       
    end if
    if (maxval(cp%vecD)/minval(cp%vecD) >= 1.0d8) then
       flg = .true.
       print *, 'max condition number'       
    end if
  end function check_convergence

  subroutine iswap( a, i, j )
    implicit none
    integer, intent(inout) :: a(:)
    integer, intent(in) :: i, j
    integer :: tmp
    tmp = a(i)
    a(i) = a(j)
    a(j) = tmp
  end subroutine iswap

  function nanleq( a, b ) result(flg)
    implicit none
    double precision, intent(in) :: a, b
    logical flg
    if ( (a > b) .or. (b == b) .and. (a /= a) ) then
       flg = .false.
    else
       flg = .true.
    end if
  end function nanleq    
  
  recursive subroutine mergesort( farray, idx, l, u, work )
    implicit none
    double precision, intent(in) :: farray(:)    
    integer, intent(inout) :: idx(:)
    integer, intent(in) :: l, u
    integer, intent(out) :: work(:)
    integer p, i, j, k
    if ( u == l ) then
       return
    else if ( u == 1 + l ) then
       if ( .not. nanleq( farray(idx(l)), farray(idx(u)) ) )  then
          call iswap( idx, l, u )
          return
       end if
    else
       p = l + (u - l) / 2
       call mergesort( farray, idx, l, p, work )
       call mergesort( farray, idx, p+1, u, work )
       i = l       
       j = p + 1
       do k = l, u
          if ( .not. nanleq( farray(idx(i)), farray(idx(j)) ) ) then
             work(k) = idx(j)
             j = j + 1
             if ( j > u ) then
                work(k+1:u) = idx(i:p)
                exit
             end if
          else
             work(k) = idx(i)
             i = i + 1
             if ( i > p ) then
                work(k+1:u) = idx(j:u)
                exit
             end if
          end if
       end do
       idx(l:u) = work(l:u)
    end if
  end subroutine mergesort

  !----- indexsort: Return Sorted Indeces -----
  subroutine indexsort( n, farray, out )
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: farray(n)
    integer, intent(out) :: out(n)
    integer :: work(n)
    integer :: i
    do i = 1, n
       out(i) = i
    end do
    call mergesort( farray, out, 1, n, work )
  end subroutine indexsort
  
  function trace( cp ) result(tr)
    implicit none

    type(CmaParam), intent(in) :: cp
    double precision :: tr
    integer :: k

    tr = 0.0d0
    if (cp%flag_diagonal == 0) then            
       do k = 1, cp%dim
          tr = tr + cp%matC(k, k)
       end do
    else if (cp%flag_diagonal == 1) then
       do k = 1, cp%dim
          tr = tr + cp%vecD(k) * cp%vecD(k)
       end do
    end if
    return
  end function trace
      
  function coordinate_length( cp, i ) result(cl)
    implicit none
    
    type(CmaParam), intent(in) :: cp
    integer, intent(in) :: i
    double precision :: cl

    if (cp%flag_diagonal == 0) then            
       cl = cp%sigma * sqrt(cp%matC(i, i))
    else if (cp%flag_diagonal == 1) then
       cl = cp%sigma * cp%vecD(i)
    end if
  end function coordinate_length

end module cmaes
