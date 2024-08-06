module problem_SI
    use MMG_models
    use random 

    IMPLICIT NONE
    double precision , allocatable :: hydro_derivative(:), upperlim(:),lowerlim(:)

    contains
    
  !------------------------------------------------------------------------
  ! 2018/05/25 (YA) Problem Interface: Initialization
    subroutine initialize_problem( filename, filename2 )
        implicit none
        character(*), intent(in) :: filename, filename2
        ! double precision , allocatable :: hydro_derivative(:), upperlim(:),lowerlim(:)
        !----------------------------------------------------------------------
        character(20) varname
        integer :: i, j, no_row, N
        character(*), parameter :: settings = '../settings.txt';
        character(*), parameter :: paramslim = '../MMG_params.csv';
        integer, parameter :: FIDST = 20 ! setting file ID
        integer, parameter :: FIDPR = 30 ! params file ID
        integer, parameter :: FIDTRFL = 40 ! filename file ID
        ! parameters
        integer, parameter :: datacolumn = 22 ! number of column of train data

        ! define dimension
        open( 10, file=filename, action='read', status='old')
        ! read( 10, * ) varname, N
        read( 10, * )  N
        ! print *, "N : ", N
        close( 10 )
        
        ! Condition of manuevering simulation
        open( FIDST, file = settings, action='read', status='old');
        read(FIDST,*)GOMI
        read(FIDST,'(I11,6X,A17)')switch_wind;
        read(FIDST,'(I11,6X,A17)')type_si; ! type of SI  :  1: all variable, 2: Forward variable only, 3: backward variable only
        read(FIDST,'(I11,6X,A17)')number_files;
        read(FIDST,'(F11.7,6X,A17)')integration_period
        read(FIDST,'(I11,6X,A17)')type_obj_func; ! type of objective func. 1: velocity, 2: location and velocity, 3: location
        close( FIDST );

        !!! read limit of coeffifients
        open(FIDPR, file=paramslim, action='read', status='old')
        ! read(FIDPR,*)GOMI
        ! read number of hydro derivative
        no_row = -1
        ! read(FIDPR,*)GOMI
        do 
            read (FIDPR, *, end=100)GOMI, GOMI, GOMI, GOMI
            no_row = no_row + 1
        end do
        100 continue
        ! write(*,*)no_row
        
        ALLOCATE(valuename(no_row), hydro_derivative(no_row), upperlim(no_row), lowerlim(no_row))
        rewind (FIDPR)  
        
        ! === read values ===
        read(FIDPR,*)GOMI
        do i = 1, no_row
            read (FIDPR, *)valuename(i), hydro_derivative(i),  lowerlim(i),  upperlim(i)
        end do
        close (FIDPR)
        
       

        ! read train data list
        ALLOCATE(trainfilename(number_files))

        open( FIDTRFL, file = filename2, action='read', status='old');
        do i =1,number_files
            read (FIDTRFL, "(a)")trainfilename(i)
        end do
        ! print *, trainfilename
        ! read train file
        allocate(step_max(number_files))
        
        do i =1,number_files
            open( 99, file = trainfilename(i), action='read', status='old');
            step_max(i) = 0
            read(99,*)GOMI
            do
                read (99, *, end=200)GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,&
                &                      GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,GOMI,GOMI
                step_max(i) = step_max(i) + 1
            end do
            200 continue
            close(99)
        end do
        ! print *, "step_max",step_max

        ALLOCATE(traindata(maxval(step_max),datacolumn))
        ALLOCATE(time_input(maxval(step_max), number_files))
        ALLOCATE(x_input(maxval(step_max), number_files), u_input(maxval(step_max), number_files))
        ALLOCATE(y_input(maxval(step_max), number_files), vm_input(maxval(step_max), number_files))
        ALLOCATE(psi_input(maxval(step_max), number_files), r_input(maxval(step_max), number_files))
        ALLOCATE(n_input(maxval(step_max), number_files), delta_Input(maxval(step_max), number_files))
        ALLOCATE(Wind_Velocity(maxval(step_max), number_files), Wind_Direction(maxval(step_max), number_files))
        
        ! initialize
        traindata(:,:) =0.0d0
        x_input(:,:)=0.0d0; 
        u_input(:,:)=0.0d0; 
        y_input(:,:)=0.0d0; 
        vm_input(:,:)=0.0d0
        psi_input(:,:)=0.0d0; 
        r_input(:,:)=0.0d0;
        Wind_Velocity(:,:) = 0.0d0; 
        Wind_Direction(:,:) = 0.0d0
        n_input(:,:) = 0.0d0; 
        delta_Input(:,:) = 0.0d0

        ! read train data values
        do i =1, number_files

            open( 99, file = trainfilename(i), action='read', status='old');
            read(99,*)GOMI
            do j = 1, step_max(i)
                read (99, *) time_input(j,i), x_input(j,i), u_input(j,i), y_input(j,i), &
                &           vm_input(j,i), psi_input(j,i), r_input(j,i),&
                            n_input(j,i),delta_Input(j,i), &
                            traindata(j,10),traindata(j,11), Wind_Velocity(j,i), Wind_Direction(j,i), &
                            traindata(j,14),traindata(j,15),traindata(j,16),traindata(j,17),traindata(j,18),&
                            traindata(j,19),traindata(j,20),traindata(j,21),traindata(j,22)
            end do
            close(99)
        end do 
        allocate(XX(N), XXmean(N)) 
        !----------------------------------------------------------------------
    end subroutine initialize_problem

    ! 2018/05/25 (YA) Problem Interface: Evaluation Main
    subroutine fobj1d(x, f)
        implicit none
        double precision, intent(in) :: x(:)
        double precision, intent(out) :: f
        integer :: N
        !----------------------------------------------------------------------
        N = 57
        call SystemIdentification( f , N,  x )
        
        !
    end subroutine fobj1d
    
    ! 2018/05/25 (YA) Problem Interface: Evaluation Called in CMA-ES
    subroutine fobj2d(x, f)
        implicit none
        double precision, intent(in) :: x(:, :)
        double precision, intent(out) :: f(:)
        !----------------------------------------------------------------------
        ! Compute the objective function value for each x(:, i) and store the
        ! objective value in f(i).
        integer i, n
        n = size(f)
        do i = 1, n
            call fobj1d(x(:, i), f(i))
            ! stop
        end do
        !----------------------------------------------------------------------
    end subroutine fobj2d
    !------------------------------------------------------------------------  
    !------------------------------------------------------------------------

    subroutine SystemIdentification(Obj, dim, XX)
            IMPLICIT NONE
            type(hydroderivative) :: mmgparams
            integer, intent(in) :: dim
            ! double precision, intent(in) :: XX(dim)
            double precision, intent(in) :: XX(dim) 
            double precision, intent(out) :: Obj  
            ! internal
            integer :: i,j,k, reset_freq
            double precision, parameter :: pi=4.0d0*datan(1.0d0);
            integer :: update_index(dim)
            ! integer :: update_index_mean(dim/2), update_index_var(dim/2)
            ! integer, allocatable :: update_index(:)
            ! double precision , allocatable ::  hydro_derivative(:),upperlim(:),lowerlim(:)
            double precision :: update_params(dim)
            ! double precision :: update_params_mean(dim/2), update_params_var(dim/2)
            double precision :: time_step_size
            double precision :: rand_n
            double precision :: k1(6),k2(6),k3(6),k4(6),right(6), state_instant(6)
            double precision :: average_input(6), average_sim(6),stdev_input(6), stdev_sim(6)
            double precision :: average_input_sin, average_sim_sin, average_input_cos, average_sim_cos
            double precision :: stdev_input_sin, stdev_sim_sin, stdev_input_cos, stdev_sim_cos 
            double precision :: ship_length, acc_limit, rdot_limit, velo_limit, angvelo_limit
            double precision , allocatable :: state_input(:,:), state_sim(:,:)
            double precision , allocatable :: standardrized_input(:,:), standardrized_sim(:,:)
            double precision , allocatable :: standardrized_input_sin(:), standardrized_sim_sin(:),&
                                             standardrized_input_cos(:), standardrized_sim_cos(:) 
            
            ! initialize mmg params (hydro derivatives)
            mmgparams%massx_nd=hydro_derivative(1)
            mmgparams%massy_nd= hydro_derivative(2)
            mmgparams%IzzJzz_nd= hydro_derivative(3)
            mmgparams%xuu_nd = hydro_derivative(4)
            mmgparams%Xvr_nd = hydro_derivative(5)
            mmgparams%Yv_nd = hydro_derivative(6)
            mmgparams%Yr_nd = hydro_derivative(7)
            mmgparams%Nv_nd = hydro_derivative(8)
            mmgparams%Nr_nd= hydro_derivative(9)
            mmgparams%CD = hydro_derivative(10)
            mmgparams%C_rY = hydro_derivative(11)
            mmgparams%C_rN= hydro_derivative(12)
            mmgparams%X_0F_nd = hydro_derivative(13)
            mmgparams%X_0A_nd= hydro_derivative(14)
            mmgparams%kt_coeff0 = hydro_derivative(15)
            mmgparams%kt_coeff1 = hydro_derivative(16)
            mmgparams%kt_coeff2 = hydro_derivative(17)
            mmgparams%t_prop = hydro_derivative(18)
            mmgparams%wP0 = hydro_derivative(19)
            mmgparams%tau_prop = hydro_derivative(20)
            mmgparams%CP_nd = hydro_derivative(21)
            mmgparams%xP_nd = hydro_derivative(22)
            mmgparams%A1 = hydro_derivative(23)
            mmgparams%A2 = hydro_derivative(24)
            mmgparams%A3 = hydro_derivative(25)
            mmgparams%A4 = hydro_derivative(26)
            mmgparams%A5 = hydro_derivative(27)
            mmgparams%A6 = hydro_derivative(28)
            mmgparams%A7 = hydro_derivative(29)
            mmgparams%A8= hydro_derivative(30)
            mmgparams%B1 = hydro_derivative(31)
            mmgparams%B2 = hydro_derivative(32)
            mmgparams%B3 = hydro_derivative(33)
            mmgparams%B4 = hydro_derivative(34)
            mmgparams%B5 = hydro_derivative(35)
            mmgparams%B6 = hydro_derivative(36)
            mmgparams%B7 = hydro_derivative(37)
            mmgparams%B8= hydro_derivative(38)
            mmgparams%C3 = hydro_derivative(39)
            mmgparams%C6 = hydro_derivative(40)
            mmgparams%C7 = hydro_derivative(41)
            mmgparams%C10= hydro_derivative(42)
            mmgparams%Jmin = hydro_derivative(43)
            mmgparams%alpha_prop= hydro_derivative(44)
            mmgparams%tR = hydro_derivative(45)
            mmgparams%aH_rudder = hydro_derivative(46)
            mmgparams%xh_rudder_nd = hydro_derivative(47)
            mmgparams%lR_nd = hydro_derivative(48)
            mmgparams%kx_rudder = hydro_derivative(49)
            mmgparams%kx_rudder_reverse = hydro_derivative(50)
            mmgparams%epsilon_rudder = hydro_derivative(51)
            mmgparams%cpr_rudder = hydro_derivative(52)
            mmgparams%gammaN = hydro_derivative(53)
            mmgparams%gammaP= hydro_derivative(54)
            mmgparams%XX0 = hydro_derivative(55)
            mmgparams%XX1 = hydro_derivative(56)
            mmgparams%XX3 = hydro_derivative(57)
            mmgparams%XX5 = hydro_derivative(58)
            mmgparams%YY1 = hydro_derivative(59)
            mmgparams%YY3 = hydro_derivative(60)
            mmgparams%YY5 = hydro_derivative(61)
            mmgparams%NN1 = hydro_derivative(62)
            mmgparams%NN2 = hydro_derivative(63)
            mmgparams%NN3 = hydro_derivative(64)

            
            ! print *, "XX : ", XX
                    
                    ! update all valiable including wind
                    ! choose values to be updated
                    !! Ndim == 57
                    update_index =(/1, 2, 3, 5, 6, 7, 8, 9, 10, 11,&
                                    12,14,18, 19, 20, 21, 22, 23, 24, 25,&
                                    26, 27, 28, 29, 30, 31, 32, 33, 34, 35,&
                                    36, 37, 38, 39, 40, 41, 42, 45, 46, 47,&
                                    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, &
                                    58, 59, 60, 61, 62, 63, 64/)   
                    ! ! Ndim == 31
                    ! update_index =(/5, 6, 7, 8, 9, 10, 11, 12, &
                    !                 14, &
                    !                 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, &
                    !                 50, 52/)    

                ! update mmg params(hydro derivatives) by CMA-ES
                do i = 1, dim
                    update_params(i) = (lowerlim(update_index(i))) +(upperlim(update_index(i)) - lowerlim(update_index(i))) & 
                                         & * ( XX(i) + 1.0d0 ) * 0.5d0
                end do

                !!! Ndim == 57
                mmgparams%massx_nd = update_params(1)
                mmgparams%massy_nd= update_params(2)
                mmgparams%IzzJzz_nd= update_params(3)
                mmgparams%Xvr_nd = update_params(4)
                mmgparams%Yv_nd = update_params(5)
                mmgparams%Yr_nd = update_params(6)
                mmgparams%Nv_nd = update_params(7)
                mmgparams%Nr_nd= update_params(8)
                mmgparams%CD = update_params(9)
                mmgparams%C_rY = update_params(10)
                mmgparams%C_rN= update_params(11)
                mmgparams%X_0A_nd= update_params(12)
                mmgparams%t_prop = update_params(13)
                mmgparams%wP0 = update_params(14)
                mmgparams%tau_prop = update_params(15)
                mmgparams%CP_nd = update_params(16)
                mmgparams%xP_nd = update_params(17)
                mmgparams%A1 = update_params(18)
                mmgparams%A2 = update_params(19)
                mmgparams%A3 = update_params(20)
                mmgparams%A4 = update_params(21)
                mmgparams%A5 = update_params(22)
                mmgparams%A6 = update_params(23)
                mmgparams%A7 = update_params(24)
                mmgparams%A8= update_params(25)
                mmgparams%B1 = update_params(26)
                mmgparams%B2 = update_params(27)
                mmgparams%B3 = update_params(28)
                mmgparams%B4 = update_params(29)
                mmgparams%B5 = update_params(30)
                mmgparams%B6 = update_params(31)
                mmgparams%B7 = update_params(32)
                mmgparams%B8= update_params(33)
                mmgparams%C3 = update_params(34)
                mmgparams%C6 = update_params(35)
                mmgparams%C7 = update_params(36)
                mmgparams%C10= update_params(37)
                mmgparams%tR = update_params(38)
                mmgparams%aH_rudder = update_params(39)
                mmgparams%xh_rudder_nd = update_params(40)
                mmgparams%lR_nd = update_params(41)
                mmgparams%kx_rudder = update_params(42)
                mmgparams%kx_rudder_reverse = update_params(43)
                mmgparams%epsilon_rudder = update_params(44)
                mmgparams%cpr_rudder = update_params(45)
                mmgparams%gammaN = update_params(46)
                mmgparams%gammaP= update_params(47)
                mmgparams%XX0 = update_params(48)
                mmgparams%XX1 = update_params(49)
                mmgparams%XX3 = update_params(50)
                mmgparams%XX5 = update_params(51)
                mmgparams%YY1 = update_params(52)
                mmgparams%YY3 = update_params(53)
                mmgparams%YY5 = update_params(54)
                mmgparams%NN1 = update_params(55)
                mmgparams%NN2 = update_params(56)
                mmgparams%NN3 = update_params(57)

                ! Ndim == 31
                ! mmgparams%Xvr_nd = update_params(1)
                ! mmgparams%Yv_nd = update_params(2)
                ! mmgparams%Yr_nd = update_params(3)
                ! mmgparams%Nv_nd = update_params(4)
                ! mmgparams%Nr_nd= update_params(5)
                ! mmgparams%CD = update_params(6)
                ! mmgparams%C_rY = update_params(7)
                ! mmgparams%C_rN= update_params(8)
                ! mmgparams%X_0A_nd= update_params(9)
                ! mmgparams%A1 = update_params(10)
                ! mmgparams%A2 = update_params(11)
                ! mmgparams%A3 = update_params(12)
                ! mmgparams%A4 = update_params(13)
                ! mmgparams%A5 = update_params(14)
                ! mmgparams%A6 = update_params(15)
                ! mmgparams%A7 = update_params(16)
                ! mmgparams%A8= update_params(17)
                ! mmgparams%B1 = update_params(18)
                ! mmgparams%B2 = update_params(19)
                ! mmgparams%B3 = update_params(20)
                ! mmgparams%B4 = update_params(21)
                ! mmgparams%B5 = update_params(22)
                ! mmgparams%B6 = update_params(23)
                ! mmgparams%B7 = update_params(24)
                ! mmgparams%B8= update_params(25)
                ! mmgparams%C3 = update_params(26)
                ! mmgparams%C6 = update_params(27)
                ! mmgparams%C7 = update_params(28)
                ! mmgparams%C10= update_params(29)
                ! mmgparams%kx_rudder_reverse = update_params(30)
                ! mmgparams%cpr_rudder = update_params(31)
            
            ! allocate variables
            allocate(state_input(maxval(step_max),6), state_sim(maxval(step_max),6))
            allocate( standardrized_input(maxval(step_max),6), standardrized_sim(maxval(step_max),6))
            allocate( standardrized_input_sin(maxval(step_max)), standardrized_sim_sin(maxval(step_max)))
            allocate( standardrized_input_cos(maxval(step_max)), standardrized_sim_cos(maxval(step_max)))
            ! initialize Object func.
            
            Obj = 0.0d0
            ship_length = 3.0d0

            do i=1,number_files
                
                !!! initalize
                average_input(:) = 0.0d0
                average_sim(:) = 0.0d0
                stdev_input(:) = 0.0d0
                stdev_sim(:) = 0.0d0
                standardrized_input(:,:) = 0.0d0
                standardrized_sim(:,:) = 0.0d0

                
                stdev_input_sin = 0.0d0
                stdev_sim_sin = 0.0d0
                stdev_input_cos = 0.0d0
                stdev_sim_cos = 0.0d0
                
                standardrized_input(:,:) = 0.0d0
                standardrized_sim(:,:) = 0.0d0
                standardrized_input_sin(:) = 0.0d0
                standardrized_sim_sin(:) = 0.0d0
                standardrized_input_cos(:) = 0.0d0
                standardrized_sim_cos(:) = 0.0d0
                
                ! mk state variables of train data
                state_input(:,1) = x_input(:,i)
                state_input(:,2) = u_input(:,i)
                state_input(:,3) = y_input(:,i)
                state_input(:,4) = vm_input(:,i)
                state_input(:,5) = psi_input(:,i)
                state_input(:,6) = r_input(:,i)
                
                ! define time step size
                ! time_step_size = nint((time_input(3,i)-time_input(2,i))*10000.0d0)/10000.0d0
                time_step_size = 10.0d0
                ! define reset frequency
                reset_freq = nint(integration_period/time_step_size) 
                
                do j= 1, step_max(i)-1
                    if(mod(j, reset_freq)==1)then
                        state_sim(j,:) = state_input(j,:)
                    endif
                    
                    state_instant(1:6) = state_sim(j,1:6) 

                    !   Runge-Kutta method
                    call MMG_LowSpeed_model(time_input(j,i), k1, state_instant,                               &
                                            delta_input(j,i),n_input(j,i), mmgparams, switch_wind, &
                                            Wind_Direction(j,i), Wind_Velocity(j,i))
                    call MMG_LowSpeed_model(time_input(j,i), k2, state_instant + 0.5d0 * k1 * time_step_size,&
                                            delta_input(j,i), n_input(j,i), mmgparams, switch_wind,&
                                             Wind_Direction(j,i), Wind_Velocity(j,i))
                    call MMG_LowSpeed_model(time_input(j,i), k3, state_instant + 0.5d0 * k2 * time_step_size, &
                                            delta_input(j,i), n_input(j,i), mmgparams, switch_wind,&
                                             Wind_Direction(j,i), Wind_Velocity(j,i))
                    call MMG_LowSpeed_model(time_input(j,i), k4, state_instant +         k3 * time_step_size, &
                                            delta_input(j,i), n_input(j,i), mmgparams, switch_wind, &
                                            Wind_Direction(j,i), Wind_Velocity(j,i))
                    
                                    
                    

                    right= (1.0d0 * k1 + 2.0d0 * k2 + 2.0d0 * k3 + 1.0d0 * k4) / 6.0d0

                    !!! limit max velocity and acceleration 
                    acc_limit = grav*1.0d0
                    rdot_limit = acc_limit/(ship_length/2.0d0)
                    velo_limit= 1.0d0
                    angvelo_limit= velo_limit/(ship_length/2.0d0)
                    
                    do k =1,6
                        if (k==1 .or. k==3)then
                            if(abs(right(k))>velo_limit)then
                                right(k) = sign(velo_limit *(2.0d0-(dble(mod(j, reset_freq))/dble(reset_freq))), right(k))
                            endif
                        elseif (k==2 .or. k==4) then
                            if (abs(right(k)) >acc_limit) then
                                right(k) = sign(acc_limit, right(k))
                            end if
                    elseif(k==5)then
                            if(abs(right(k))>angvelo_limit)then
                                right(k) = sign(angvelo_limit * (2.0d0-(dble(mod(j, reset_freq))/dble(reset_freq))),right(k))
                            endif
                    elseif(k == 6) then
                            if (abs(right(k)) > rdot_limit)then
                                right(k) = sign(rdot_limit,right(k))
                            endif
                        end if
                    end do

                    do k = 1,6 
                        state_sim(j+1,k) = state_instant(k) + right(k) * time_step_size
                    end do

                    

               
                !!! compute objecti func !!!      
                do k =1,6
                    ! average of state val
                    if(k==5)then
                    
                        average_input_sin = average_input_sin/step_max(i)
                        average_sim_sin = average_sim_sin/step_max(i)
                        average_input_cos = average_input_cos/step_max(i)
                        average_sim_cos = average_sim_cos/step_max(i)
                        ! standard deviation of val

                        do j = 1,step_max(i)
                            stdev_input_sin = stdev_input_sin + (sin(state_input(j,k))-average_input_sin)**2.0d0
                            stdev_sim_sin   = stdev_sim_sin   + (sin(state_sim(j,k))  -average_sim_sin)**2.0d0
                            stdev_input_cos = stdev_input_cos + (cos(state_input(j,k))-average_input_cos)**2.0d0
                            stdev_sim_cos   = stdev_sim_cos   + (cos(state_sim(j,k))  -average_sim_cos)**2.0d0
                        
                        end do
                        
                        
                        stdev_input_sin = sqrt(abs(stdev_input_sin/step_max(i)))
                        stdev_sim_sin = sqrt(abs(stdev_sim_sin/step_max(i)))
                        stdev_input_cos = sqrt(abs(stdev_input_cos/step_max(i)))
                        stdev_sim_cos = sqrt(abs(stdev_sim_cos/step_max(i)))
                       
                        
                        !!! standarization
                        do j = 1,step_max(i)
                            standardrized_input_sin(j) = (sin(state_input(j,k))-average_input_sin)/stdev_input_sin
                            standardrized_sim_sin(j) = (sin(state_sim(j,k))-average_sim_sin)/stdev_sim_sin
                            standardrized_input_cos(j) = (cos(state_input(j,k))-average_input_cos)/stdev_input_cos
                            standardrized_sim_cos(j) = (cos(state_sim(j,k))-average_sim_cos)/stdev_sim_cos
                        end do
                    else                
                        average_input(k)= average_input(k)/step_max(i)
                        average_sim(k) = average_sim(k)/step_max(i)
                        ! standard deviation of val
                        do j = 1,step_max(i)
                            stdev_input(k) = stdev_input(k) + (state_input(j,k)-average_input(k))**2.0d0
                            stdev_sim(k) = stdev_sim(k) + (state_sim(j,k)-average_sim(k))**2.0d0
                        end do
                        stdev_input(k) = sqrt(abs(stdev_input(k)/step_max(i)))
                        stdev_sim(k) = sqrt(abs(stdev_sim(k)/step_max(i)))
                        
                        !!! standarization
                        do j = 1,step_max(i)
                            standardrized_input(j,k) = (state_input(j,k)-average_input(k))/stdev_input(k)
                            standardrized_sim(j,k) = (state_sim(j,k)-average_sim(k))/stdev_sim(k)
                        end do
                    endif
                end do
                
                
                ! if(type_obj_func==1)then
                    ! compute obj func form standardrized u, vm and r     
                    do j = 1,step_max(i)
                        Obj = Obj + ( (standardrized_input(j,2) - standardrized_sim(j,2)) **2.0d0 + &
                                    (standardrized_input(j,4) - standardrized_sim(j,4)) **2.0d0 + &
                                    (standardrized_input(j,6) - standardrized_sim(j,6)) **2.0d0 ) * time_step_size
                    end do
                    
                ! endif
            enddo


        end subroutine SystemIdentification

 end module problem_SI


