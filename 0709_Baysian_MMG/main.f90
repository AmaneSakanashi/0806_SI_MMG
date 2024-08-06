!----------------------------------------------------------
! Dependency Tree
! - cmaes
!   - eigen (module for eigendecomposition)
!   - random (module for random number generator, MT is used)
! - random (only to set the random seed)

! - benchmark (only for benchmark test)
!----------------------------------------------------------
program main
  use cmaes, only : CmaParam, initialize_cmaes, generate_candidate_solutions, &
       update_parameters, handle_box_boundary, check_convergence, &
       trace, coordinate_length, clear, printcp, randn, &
       compute_penalty, integer_mutation, writecp, readcp
  use eigen, only : myddot
  use random, only : uniform, setseed, rand_load, rand_save
  use problem_SI, only: initialize_problem, fobj2d
!   use problem, only : SystemIdentification,initialize_problem, fobj2d
  implicit none

  type CmaOption
     !! Input parameters
     logical :: flg_min, flg_restart, flg_resume, flg_save, flg_visualization, flg_integer_mutation
     integer ::  dim, initpopsize, maxpopsize, print_span, flag_diagonal, flag_active, maxeval, maxiter, seed
     double precision ::  ftarget, tolx, tolf, sigma_correction_factor
     double precision, allocatable :: lb(:), ub(:), init_mean(:), init_std(:)
     double precision, allocatable :: stairwidth(:) ! for Mixed-Integer optimization
  end type CmaOption

  !!===============================================================================================  
  !! Constant and Declaration
  logical :: flg_debug = .false. ! DEBUG MODE: compile with -Wall -g option
  integer, parameter :: FID = 50 ! Output file ID
  integer, parameter :: VFID = 60 ! Visualization file ID
  integer :: i, j, k, popsize, itr, neval
  double precision :: df, dx
  double precision, allocatable :: arx(:, :), arxfeas(:, :), arxmod(:, :)
  double precision, allocatable :: arfraw(:), arf(:), arffeas(:)
  double precision, allocatable :: init_mean(:), init_std(:)
  type(CmaParam) :: cp
  type(CmaOption) :: co
  character(99) :: input_filename, input_filename2, input_filename3
  character(*), parameter :: OUTPUT = 'OutputCmaes.txt'
  character(*), parameter :: RESUME = 'ResumeCmaes.dat'
  character(*), parameter :: SAVEMT = 'ResumeMT.dat'
  character(*), parameter :: OUTVIS = 'VisualizeCmaes.txt'
  character(*), parameter :: OPTION = 'DefaultInitCmaes.txt'
  character(*), parameter :: DoubleFormat = '(ES25.15E3)'
  character(*), parameter :: IntFormat = '(I25)'
  !!===============================================================================================  
  if ( iargc() < 3 ) then
     write(*, *) 'input files required.'
     stop
  else
     call getarg( 1_4, input_filename ) !InitCmaesXX.txt
     call getarg( 2_4, input_filename2 ) !InitProblemXX.txt
     call getarg( 3_4, input_filename3 ) ! XX_train_file_list.txt
     write(*,*)iargc(), input_filename, input_filename2, input_filename3
     call load_input_file( co, input_filename )
     call initialize_problem( input_filename2, input_filename3 )
     write(*,*)"###FILE LOADED###"
  end if
  allocate( init_mean(co%dim), init_std(co%dim) )
  init_mean(:) = co%init_mean(:)
  init_std(:) = co%init_std(:)  
  
  !!===============================================================================================  
  !! Initialization
  if (co%flg_resume) then
     !! Resume
     call readcp( cp, RESUME )
     call rand_load( SAVEMT )
     open( FID, file=OUTPUT, access='append', status='old' )
     if (co%flg_visualization) then
        open( VFID, file=OUTVIS, access='append', status='old' )
     endif
     if ( (co%maxiter .le. cp%iter) .or. (co%maxeval .le. cp%neval) ) then
        co%maxiter = co%maxiter + cp%iter
        co%maxeval = co%maxeval + cp%neval
     end if
     cp%iter = cp%iter + 1
     cp%neval = cp%neval + 1
  else
     !! New Run
     call setseed( co%seed )
     popsize = co%initpopsize
     call initialize_cmaes( cp, co%dim, co%lb, co%ub, init_mean, init_std, co%stairwidth, &
          co%flag_diagonal, popsize, co%flag_active, co%sigma_correction_factor )
     !! Open log filesj
     open( FID, file=OUTPUT, status='replace' )
     !write( FID, * ) '#Dim:', co%dim, 'lambda:', popsize, 'ub:', co%ub, 'lb:', co%lb
     write( FID, '(A100)' ) '#iter, minval(arffeas), minval(arf), df, dx, cp%sigma, cp%vecD, cp%xmean, cp%diagSqrtC'
     if (co%flg_visualization) then
        open( VFID, file=OUTVIS, status='replace' )
        write( VFID, '(A100)' ) '#arf, arffeas, arx, arxfeas, sigma, xmean, matC'        
     endif
  end if
  !! Memory Allocation
  allocate(arx(cp%dim, cp%lambda), arxfeas(cp%dim, cp%lambda), arxmod(cp%dim, cp%lambda))
  allocate(arfraw(cp%lambda), arffeas(cp%lambda), arf(cp%lambda))  
   ! write(*,*)"###INITIALIZATION COMPLETED###"
  !!===============================================================================================
  !! Main Loop
  ! $omp parallel
  ! $omp do
  do i = cp%iter, co%maxiter

     call debug_print(flg_debug, i, 'new_iter')
      ! write(*,*)"###DEGUB NO.",i,"###"
     ! Sampling Step
     ! arx : candidate points, possibly infeasible, being used to update param. 
     ! arxfeas : feasible candidate points, to be evaluated
     call generate_candidate_solutions(cp, arx)
     if ( co%flg_integer_mutation ) then
        call integer_mutation(cp, arx, arxmod)
     else
        arxmod(:, :) = arx(:, :)
     end if
     call handle_box_boundary(cp, arxmod, arxfeas)
     call debug_print(flg_debug, i, 'sampling')
     
     ! Evaluation Step
     ! arfraw : objective function value of 'arxfeas'
   !   write(*,*)"################"
     
     call fobj2d( arxfeas(:, :), arfraw(:) )
     
     cp%neval = cp%neval + cp%lambda
     call debug_print(flg_debug, i, 'evaluation')     
     
     ! Force Minimization
     ! arffeas = + or - arfraw (depending on whether it is min or max)
     if (co%flg_min) then
        arffeas(:) = arfraw(:)
     else
        arffeas(:) = - arfraw(:)
     end if
     
     ! Update Step
     ! arffeas : f(arxfeas), function value to be minimized
     ! arf : f(arxfeas) + penalty, used as a virtual f(arx) to update param.
     call compute_penalty(cp, arxfeas, arffeas, arxmod, arf)
     call debug_print(flg_debug, i, 'up_penalty')     
     call update_parameters(cp, arx, arxmod, arf)
     call debug_print(flg_debug, i, 'up_param')
     
     ! Print and Write Internal State
     if (mod(i, co%print_span) == 0) then
        if (cp%flg_fhistory == 1) then
           k = size(cp%fhistory)
        else
           k = cp%idx_fhistory - 1
        end if
        df = maxval(cp%fhistory(1:k)) * cp%sigma * cp%sigma * trace(cp) / cp%dble_dim
        dx = 0.0d0
        do k = 1, cp%dim
           dx = max(dx, coordinate_length(cp, k))
        end do
        
        print *, cp%iter, minval(arffeas), minval(arf), df, dx, cp%sigma, maxval(cp%vecD)/minval(cp%vecD)    
        write( FID, IntFormat, advance='no' ) cp%iter
        write( FID, DoubleFormat, advance='no' ) minval(arffeas)
        write( FID, DoubleFormat, advance='no' ) minval(arf)
        write( FID, DoubleFormat, advance='no' ) df
        write( FID, DoubleFormat, advance='no' ) dx
        write( FID, DoubleFormat, advance='no' ) cp%sigma
        do j = 1, cp%dim
           write( FID, DoubleFormat, advance='no' ) cp%vecD(j)
        end do
        do j = 1, cp%dim
           write( FID, DoubleFormat, advance='no' ) cp%xmean(j)
        end do
        do j = 1, cp%dim
           write( FID, DoubleFormat, advance='no' ) cp%diagSqrtC(j)
        end do
        write( FID, * )
     end if
     call debug_print(flg_debug, i, 'write_out')     
     
     ! Save Internal State to Resume
     !if (co%flg_save) then
        !call writecp( cp, 'ResumeCmaes.dat' )
        !call rand_save( 'ResumeMT.dat' )
     !end if
     call debug_print(flg_debug, i, 'save')          

     ! Write Internal State for Visualization
     if (co%flg_visualization) then
        do j = 1, cp%lambda        
           write( VFID, DoubleFormat, advance='no' ) arf(j)
        end do
        do j = 1, cp%lambda                
           write( VFID, DoubleFormat, advance='no' ) arffeas(j)
        end do
        do k = 1, cp%lambda
           do j = 1, cp%dim
              write( VFID, DoubleFormat, advance='no' ) arx(j, k)
           enddo
        enddo
        do k = 1, cp%lambda
           do j = 1, cp%dim                   
              write( VFID, DoubleFormat, advance='no' ) arxfeas(j, k)
           enddo
        enddo
        write( VFID, DoubleFormat, advance='no' ) cp%sigma
        do j = 1, cp%dim
           write( VFID, DoubleFormat, advance='no' ) cp%xmean(j)
        end do
        if ( cp%flag_diagonal == 0 ) then
           do k = 1, cp%dim
              do j = 1, cp%dim           
                 write( VFID, DoubleFormat, advance='no' ) cp%matC(j, k)
              enddo
           enddo
        else if ( cp%flag_diagonal == 1 ) then
           do k = 1, cp%dim
              write( VFID, DoubleFormat, advance='no' ) cp%vecD(k)
           enddo
        end if
        write( VFID, * )
     endif
     call debug_print(flg_debug, i, 'visualize_out')               
          
     ! Check Stopping Criteria
     if ( minval(arffeas) <= co%ftarget ) then
        print *, 'ftarget'
        exit
     else if ( cp%neval >= co%maxeval ) then
        print *, 'maxeval'
        exit
     else if ( cp%iter >= co%maxiter ) then
        print *, 'maxiter'
        exit
     end if
     call debug_print(flg_debug, i, 'check_stop')                    

     ! Restart if converged
     if (co%flg_restart .and. check_convergence(cp, co%tolx, co%tolf)) then
        ! Restart with Increasing Population Size
        popsize = min(2 * cp%lambda, co%maxpopsize)
        print *, 'Restart with popsize = ', popsize

        ! Bug Fix: 20180222
        itr = cp%iter
        neval = cp%neval

        ! Re-initialize the distribution parameters
        do j = 1, cp%dim
           if ( co%lb(j) > - huge(co%lb(j)) .and. co%ub(j) < huge(co%lb(j)) ) then
              ! random initialization 
              init_mean(j) = uniform(co%lb(j) + 0.1d0 * (co%ub(j) - co%lb(j)), co%ub(j) - 0.1d0 * (co%ub(j) - co%lb(j)))
              init_std(j) = (co%ub(j) - co%lb(j)) / 6.0d0
           else
              ! starting from the initial given parameters
              init_mean(j) = co%init_mean(j)
              init_std(j) = co%init_std(j)
           end if
        end do

        ! Deallocate
        call clear( cp )
        deallocate(arx, arxfeas, arxmod) ! Bug fixed on Jan. 16, 2018
        deallocate(arfraw, arffeas, arf)
        ! Re-initialize strategy parameters
        call initialize_cmaes( cp, co%dim, co%lb, co%ub, init_mean, init_std, co%stairwidth, &
             co%flag_diagonal, popsize, co%flag_active, co%sigma_correction_factor )
        ! Re-allocate ! Bug fixed on Jan. 16, 2018
        allocate(arx(cp%dim, cp%lambda), arxfeas(cp%dim, cp%lambda), arxmod(cp%dim, cp%lambda))
        allocate(arfraw(cp%lambda), arffeas(cp%lambda), arf(cp%lambda))
        ! Reset itr, neval
        cp%iter = itr
        cp%neval = neval        
        
     end if
     cp%iter = cp%iter + 1
  end do
  ! $omp enddo
  ! $omp end parallel
  
contains

  subroutine load_input_file( co, fname )
    implicit none
    type(CmaOption), intent(inout) :: co
    character(*), intent(in) :: fname
    character(99) varname
    open( 10, file=fname, action='read', status='old')
    read( 10, * ) varname, co%flg_min
    read( 10, * ) varname, co%flg_restart
    read( 10, * ) varname, co%flg_resume
    read( 10, * ) varname, co%flg_save
    read( 10, * ) varname, co%flg_visualization
    read( 10, * ) varname, co%flg_integer_mutation
    read( 10, * ) varname, co%dim
    read( 10, * ) varname, co%initpopsize
    read( 10, * ) varname, co%maxpopsize
    read( 10, * ) varname, co%print_span
    read( 10, * ) varname, co%flag_diagonal
    read( 10, * ) varname, co%flag_active
    read( 10, * ) varname, co%maxeval
    read( 10, * ) varname, co%maxiter
    read( 10, * ) varname, co%seed
    read( 10, * ) varname, co%ftarget
    read( 10, * ) varname, co%tolx
    read( 10, * ) varname, co%tolf
    read( 10, * ) varname, co%sigma_correction_factor    
    allocate(co%lb(co%dim), co%ub(co%dim), co%init_mean(co%dim), co%init_std(co%dim), co%stairwidth(co%dim))
    read( 10, * ) varname, co%lb
    read( 10, * ) varname, co%ub
    read( 10, * ) varname, co%init_mean
    read( 10, * ) varname, co%init_std
    read( 10, * ) varname, co%stairwidth
    close( 10 )
  end subroutine load_input_file

  subroutine debug_print( flg_debug, iter, keyword )
    implicit none
    logical, intent(in) :: flg_debug
    integer, intent(in) :: iter
    character(*), intent(in) :: keyword
    if ( flg_debug ) then
       write(0, *) iter, keyword
    end if

  end subroutine debug_print
  
  
  
end program main

    
!!----------------------------------------------------------
!! Dependency Tree
!! - cmaes
!!   - eigen (module for eigendecomposition)
!!   - random (module for random number generator, MT is used)
!! - random (only to set the random seed)
!! - benchmark (only for benchmark test)
!!----------------------------------------------------------
!program main
!  use cmaes, only : CmaParam, initialize_cmaes, generate_candidate_solutions, &
!       update_parameters, handle_box_boundary, check_convergence, &
!       trace, coordinate_length, clear, printcp, randn, &
!       compute_penalty, integer_mutation, writecp, readcp
!  use eigen, only : myddot
!  use random, only : uniform, setseed, rand_load, rand_save
!  use problem, only : initialize_problem, fobj2d
!  implicit none
!
!  type CmaOption
!     !! Input parameters
!     logical :: flg_min, flg_restart, flg_resume, flg_save, flg_visualization, flg_integer_mutation
!     integer ::  dim, initpopsize, maxpopsize, print_span, flag_diagonal, flag_active, maxeval, maxiter, seed
!     double precision ::  ftarget, tolx, tolf, sigma_correction_factor
!     double precision, allocatable :: lb(:), ub(:), init_mean(:), init_std(:)
!     double precision, allocatable :: stairwidth(:) ! for Mixed-Integer optimization
!  end type CmaOption
!
!  !!===============================================================================================  
!  !! Constant and Declaration
!  logical :: flg_debug = .false. ! DEBUG MODE: compile with -Wall -g option
!  integer, parameter :: FID = 50 ! Output file ID
!  integer, parameter :: VFID = 60 ! Visualization file ID
!  integer :: i, j, k, popsize, itr, neval
!  double precision :: df, dx
!  double precision, allocatable :: arx(:, :), arxfeas(:, :), arxmod(:, :)
!  double precision, allocatable :: arfraw(:), arf(:), arffeas(:)
!  double precision, allocatable :: init_mean(:), init_std(:)
!  type(CmaParam) :: cp
!  type(CmaOption) :: co
!  character(99) :: input_filename, input_filename2
!  character(*), parameter :: OUTPUT = 'OutputCmaes.txt'
!  character(*), parameter :: RESUME = 'ResumeCmaes.dat'
!  character(*), parameter :: SAVEMT = 'ResumeMT.dat'
!  character(*), parameter :: OUTVIS = 'VisualizeCmaes.txt'
!  character(*), parameter :: OPTION = 'DefaultInitCmaes.txt'
!  character(*), parameter :: DoubleFormat = '(ES25.15E3)'
!  character(*), parameter :: IntFormat = '(I25)'
!  !!===============================================================================================  
!  if ( iargc() < 2 ) then
!     write(*, *) 'input files required.'
!     stop
!  else
!     write(*,*)iargc()
!     call getarg( 1_4, input_filename )
!     call getarg( 2_4, input_filename2 )     
!     call load_input_file( co, input_filename )
!     call initialize_problem( input_filename2 )
!  end if
!  allocate( init_mean(co%dim), init_std(co%dim) )
!  init_mean(:) = co%init_mean(:)
!  init_std(:) = co%init_std(:)  
!  
!  !!===============================================================================================  
!  !! Initialization
!  if (co%flg_resume) then
!     !! Resume
!     call readcp( cp, RESUME )
!     call rand_load( SAVEMT )
!     open( FID, file=OUTPUT, access='append', status='old' )
!     if (co%flg_visualization) then
!        open( VFID, file=OUTVIS, access='append', status='old' )
!     endif
!     if ( (co%maxiter .le. cp%iter) .or. (co%maxeval .le. cp%neval) ) then
!        co%maxiter = co%maxiter + cp%iter
!        co%maxeval = co%maxeval + cp%neval
!     end if
!     cp%iter = cp%iter + 1
!     cp%neval = cp%neval + 1
!  else
!     !! New Run
!     call setseed( co%seed )
!     popsize = co%initpopsize
!     call initialize_cmaes( cp, co%dim, co%lb, co%ub, init_mean, init_std, co%stairwidth, &
!          co%flag_diagonal, popsize, co%flag_active, co%sigma_correction_factor )
!     !! Open log files
!     open( FID, file=OUTPUT, status='replace' )
!     !write( FID, * ) '#Dim:', co%dim, 'lambda:', popsize, 'ub:', co%ub, 'lb:', co%lb
!     write( FID, '(A100)' ) '#iter, minval(arffeas), minval(arf), df, dx, cp%sigma, cp%vecD, cp%xmean, cp%diagSqrtC'
!     if (co%flg_visualization) then
!        open( VFID, file=OUTVIS, status='replace' )
!        write( VFID, '(A100)' ) '#arf, arffeas, arx, arxfeas, sigma, xmean, matC'        
!     endif
!  end if
!  !! Memory Allocation
!  allocate(arx(cp%dim, cp%lambda), arxfeas(cp%dim, cp%lambda), arxmod(cp%dim, cp%lambda))
!  allocate(arfraw(cp%lambda), arffeas(cp%lambda), arf(cp%lambda))  
!
!  !!===============================================================================================
!  !! Main Loop
!  do i = cp%iter, co%maxiter
!
!     call debug_print(flg_debug, i, 'new_iter')
!
!     ! Sampling Step
!     ! arx : candidate points, possibly infeasible, being used to update param. 
!     ! arxfeas : feasible candidate points, to be evaluated
!     call generate_candidate_solutions(cp, arx)
!     if ( co%flg_integer_mutation ) then
!        call integer_mutation(cp, arx, arxmod)
!     else
!        arxmod(:, :) = arx(:, :)
!     end if
!     call handle_box_boundary(cp, arxmod, arxfeas)
!     call debug_print(flg_debug, i, 'sampling')
!     
!     ! Evaluation Step
!     ! arfraw : objective function value of 'arxfeas'
!     call fobj2d( arxfeas(:, :), arfraw(:) )
!     cp%neval = cp%neval + cp%lambda

!     call debug_print(flg_debug, i, 'evaluation')     
!     
!     ! Force Minimization
!     ! arffeas = + or - arfraw (depending on whether it is min or max)
!     if (co%flg_min) then
!        arffeas(:) = arfraw(:)
!     else
!        arffeas(:) = - arfraw(:)
!     end if
!     
!     ! Update Step
!     ! arffeas : f(arxfeas), function value to be minimized
!     ! arf : f(arxfeas) + penalty, used as a virtual f(arx) to update param.
!     call compute_penalty(cp, arxfeas, arffeas, arxmod, arf)
!     call debug_print(flg_debug, i, 'up_penalty')     
!     call update_parameters(cp, arx, arxmod, arf)
!     call debug_print(flg_debug, i, 'up_param')
!     
!     ! Print and Write Internal State
!     if (mod(i, co%print_span) == 0) then
!        if (cp%flg_fhistory == 1) then
!           k = size(cp%fhistory)
!        else
!           k = cp%idx_fhistory - 1
!        end if
!        df = maxval(cp%fhistory(1:k)) * cp%sigma * cp%sigma * trace(cp) / cp%dble_dim
!        dx = 0.0d0
!        do k = 1, cp%dim
!           dx = max(dx, coordinate_length(cp, k))
!        end do
!        
!        print *, cp%iter, minval(arffeas), minval(arf), df, dx, cp%sigma, maxval(cp%vecD)/minval(cp%vecD)    
!        write( FID, IntFormat, advance='no' ) cp%iter
!        write( FID, DoubleFormat, advance='no' ) minval(arffeas)
!        write( FID, DoubleFormat, advance='no' ) minval(arf)
!        write( FID, DoubleFormat, advance='no' ) df
!        write( FID, DoubleFormat, advance='no' ) dx
!        write( FID, DoubleFormat, advance='no' ) cp%sigma
!        do j = 1, cp%dim
!           write( FID, DoubleFormat, advance='no' ) cp%vecD(j)
!        end do
!        do j = 1, cp%dim
!           write( FID, DoubleFormat, advance='no' ) cp%xmean(j)
!        end do
!        do j = 1, cp%dim
!           write( FID, DoubleFormat, advance='no' ) cp%diagSqrtC(j)
!        end do
!        write( FID, * )
!     end if
!     call debug_print(flg_debug, i, 'write_out')     
!     
!     ! Save Internal State to Resume
!     if (co%flg_save) then
!        call writecp( cp, 'ResumeCmaes.dat' )
!        call rand_save( 'ResumeMT.dat' )
!     end if
!     call debug_print(flg_debug, i, 'save')          
!
!     ! Write Internal State for Visualization
!     if (co%flg_visualization) then
!        do j = 1, cp%lambda        
!           write( VFID, DoubleFormat, advance='no' ) arf(j)
!        end do
!        do j = 1, cp%lambda                
!           write( VFID, DoubleFormat, advance='no' ) arffeas(j)
!        end do
!        do k = 1, cp%lambda
!           do j = 1, cp%dim
!              write( VFID, DoubleFormat, advance='no' ) arx(j, k)
!           enddo
!        enddo
!        do k = 1, cp%lambda
!           do j = 1, cp%dim                   
!              write( VFID, DoubleFormat, advance='no' ) arxfeas(j, k)
!           enddo
!        enddo
!        write( VFID, DoubleFormat, advance='no' ) cp%sigma
!        do j = 1, cp%dim
!           write( VFID, DoubleFormat, advance='no' ) cp%xmean(j)
!        end do
!        if ( cp%flag_diagonal == 0 ) then
!           do k = 1, cp%dim
!              do j = 1, cp%dim           
!                 write( VFID, DoubleFormat, advance='no' ) cp%matC(j, k)
!              enddo
!           enddo
!        else if ( cp%flag_diagonal == 1 ) then
!           do k = 1, cp%dim
!              write( VFID, DoubleFormat, advance='no' ) cp%vecD(k)
!           enddo
!        end if
!        write( VFID, * )
!     endif
!     call debug_print(flg_debug, i, 'visualize_out')               
!          
!     ! Check Stopping Criteria
!     if ( minval(arffeas) <= co%ftarget ) then
!        print *, 'ftarget'
!        exit
!     else if ( cp%neval >= co%maxeval ) then
!        print *, 'maxeval'
!        exit
!     else if ( cp%iter >= co%maxiter ) then
!        print *, 'maxiter'
!        exit
!     end if
!     call debug_print(flg_debug, i, 'check_stop')                    
!
!     ! Restart if converged
!     if (co%flg_restart .and. check_convergence(cp, co%tolx, co%tolf)) then
!        ! Restart with Increasing Population Size
!        popsize = min(2 * cp%lambda, co%maxpopsize)
!        print *, 'Restart with popsize = ', popsize
!
!        ! Bug Fix: 20180222
!        itr = cp%iter
!        neval = cp%neval
!
!        ! Re-initialize the distribution parameters
!        do j = 1, cp%dim
!           if ( co%lb(j) > - huge(co%lb(j)) .and. co%ub(j) < huge(co%lb(j)) ) then
!              ! random initialization 
!              init_mean(j) = uniform(co%lb(j) + 0.1d0 * (co%ub(j) - co%lb(j)), co%ub(j) - 0.1d0 * (co%ub(j) - co%lb(j)))
!              init_std(j) = (co%ub(j) - co%lb(j)) / 6.0d0
!           else
!              ! starting from the initial given parameters
!              init_mean(j) = co%init_mean(j)
!              init_std(j) = co%init_std(j)
!           end if
!        end do
!
!        ! Deallocate
!        call clear( cp )
!        deallocate(arx, arxfeas, arxmod) ! Bug fixed on Jan. 16, 2018
!        deallocate(arfraw, arffeas, arf)
!        ! Re-initialize strategy parameters
!        call initialize_cmaes( cp, co%dim, co%lb, co%ub, init_mean, init_std, co%stairwidth, &
!             co%flag_diagonal, popsize, co%flag_active, co%sigma_correction_factor )
!        ! Re-allocate ! Bug fixed on Jan. 16, 2018
!        allocate(arx(cp%dim, cp%lambda), arxfeas(cp%dim, cp%lambda), arxmod(cp%dim, cp%lambda))
!        allocate(arfraw(cp%lambda), arffeas(cp%lambda), arf(cp%lambda))
!        ! Reset itr, neval
!        cp%iter = itr
!        cp%neval = neval        
!        
!     end if
!     cp%iter = cp%iter + 1
!  end do
!  
!contains
!
!  subroutine load_input_file( co, fname )
!    implicit none
!    type(CmaOption), intent(inout) :: co
!    character(*), intent(in) :: fname
!    character(99) varname
!    open( 10, file=fname, action='read', status='old')
!    read( 10, * ) varname, co%flg_min
!    read( 10, * ) varname, co%flg_restart
!    read( 10, * ) varname, co%flg_resume
!    read( 10, * ) varname, co%flg_save
!    read( 10, * ) varname, co%flg_visualization
!    read( 10, * ) varname, co%flg_integer_mutation
!    read( 10, * ) varname, co%dim
!    read( 10, * ) varname, co%initpopsize
!    read( 10, * ) varname, co%maxpopsize
!    read( 10, * ) varname, co%print_span
!    read( 10, * ) varname, co%flag_diagonal
!    read( 10, * ) varname, co%flag_active
!    read( 10, * ) varname, co%maxeval
!    read( 10, * ) varname, co%maxiter
!    read( 10, * ) varname, co%seed
!    read( 10, * ) varname, co%ftarget
!    read( 10, * ) varname, co%tolx
!    read( 10, * ) varname, co%tolf
!    read( 10, * ) varname, co%sigma_correction_factor    
!    allocate(co%lb(co%dim), co%ub(co%dim), co%init_mean(co%dim), co%init_std(co%dim), co%stairwidth(co%dim))
!    read( 10, * ) varname, co%lb
!    read( 10, * ) varname, co%ub
!    read( 10, * ) varname, co%init_mean
!    read( 10, * ) varname, co%init_std
!    read( 10, * ) varname, co%stairwidth
!    close( 10 )
!  end subroutine load_input_file
!
!  subroutine debug_print( flg_debug, iter, keyword )
!    implicit none
!    logical, intent(in) :: flg_debug
!    integer, intent(in) :: iter
!    character(*), intent(in) :: keyword
!    if ( flg_debug ) then
!       write(0, *) iter, keyword
!    end if
!  end subroutine debug_print
!  
!  

!  
!end program main
