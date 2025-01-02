
! This is the main code file for Hanif Kawousi's Msc-project.


! This is the main code file for Hanif Kawousi's Msc-project.
! To run it, make sure to also have the Makefile and main.f90 in the same directory.
! Then run "make" in the terminal for running in evolutionary mode, or
! run "make gen_(your evolutionary generation number) to run in ecological mode.
! For the ecological mode, make sure to have the csv file produced by the evolutionary mode
! for your wanted generation in the same directory.
! ex: For the thesis, I have run 300 generations, producing 300 csv files.
! For ecological experimentations, in the terminal, I ran the following commands: 
!  "make gen_300" - this compiles the code with the generation number 300.
! If changes are made in the code and you do not want to recompile with all the
! outputs in the terminal, you can run "make model". 
! If you want to clear all produced data from previous compilation, run 
! "make clean". 

!      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!        NOTE: This is the 2Gene Model WITH BACKGROUND MORTALITY AND RISK SIGNALLING
!      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



!------------------------------------
!--------PARAMETER MODULE------------
!------------------------------------
module params
    implicit none

! NON-OBJECT: change these to access or avoid data that will affect compilation.

    logical, public, parameter :: IS_DEBUG_SCREEN = .FALSE.
    logical, public, parameter :: IS_SAVE_POPULATION = .TRUE.

    !Both must be the same: .TRUE. or .FALSE., otherwise error-message
    logical, public, parameter :: IS_ZIP_DATA =  .FALSE. ! compress large data before it takes a lot of space
    logical, public, parameter :: IS_DEBUG_DATA= .FALSE. !Run all recording of a value at each timestep: 
                                                        !weight, emotion_state, food_availability, food_eaten

    character(*), public, parameter :: ZIP_PROG = 'gzip' !program can be changed in the future.

    ! Numerical code for missing value
    real, parameter, public :: MISSING = -9999.0
    integer, parameter, public :: UNDEFINED = -9999

! ============= BIRD ====================
    integer, parameter, public :: HUNGER_MIN = 0, HUNGER_MAX = 10 
    integer, parameter, public :: FEAR_MIN = 0, FEAR_MAX = 10  

    real, parameter, public :: BIRD_INITIAL_WEIGHT = 20.0 !20 grams 
    real, parameter, public :: BIRD_MAXIMUM_WEIGHT_MULTIPLICATOR = 1.25 !This defines birds maximum weight.
    real, parameter, public :: BIRD_MINIMUM_WEIGHT_MULTIPLICATOR = 0.75 !This defines birds minimum weight. 
    real, parameter, public :: WEIGHT_REDUCTION_CONSTANT = 0.01 
    !Only for when the bird encounters predator, but escapes. 
    !The added stress of escape will manifest in excess weight reduction.

    real, parameter, public :: METABOLIC_COST_CONSTANT = 0.0275 
    ! The cost of life to be multiplied with birds weight.


    real, parameter, public :: HUNGER_MAX_DECREMENT = 0.5
    real, parameter, public :: HUNGER_DECREMENT_MULTIPLICATOR = 0.075 !decreased hunger of the bird in feeding situations.
    real, parameter, public :: FEAR_INCREMENT_AFTER_ATTACK = 0.35  !increased fear of the bird in attack situations.
    real, parameter, public :: HUNGER_STEP_FRACTION =  0.2 
    ! This fraction is used to increment the bird's fear-hunger state when it chooses to feed, 
    !reflecting the bird's strategy to balance its fear-hunger state by prioritizing fear over hunger. 
    !This mechanism ensures that the bird's behavior is dynamic and responsive to its current emotional state, 
    !potentially influencing its survival in the environment.

    real, parameter, public :: HUNGER_PLASTICITY = 0.20 ! Param for variation in hunger between birds. 
    real, parameter, public :: FEAR_PLASTICITY = 0.20

    real, parameter, public :: EXPLORATION_PROBABILITY = 0.001 !Prob. of exlploring new patches

    !Variables for fraction in fuction bird_eat_fraction_when_fear
    real, parameter, public :: FRAC_S1 = 0.5
    real, parameter, public :: FRAC_S2 = 1.0
! ========================================


! =========== PREDATOR ====================
   real, parameter, public :: FREQUENCY_OF_PREDATOR_MIN = 0.001
   real, parameter, public :: FREQUENCY_OF_PREDATOR_MAX = 0.015

   ! For signal generations only:
   ! if 1 then we are in evolutionary gens. (1-6)
   real, parameter, public :: PRED_FREQUENCY_MULTIPLIER_SIGNAL_GENS_ONLY = 1

!parameter for predation_risk, risk of the bird to meet predator in a particular environment.
! if the bird sees the predator, the predator sees the bird.

  real, parameter, public :: ATTACK_RATE = 0.15 
  !risk of bird being eaten by predator if attacked.
  !is 0.0 when ecological experiments are being conducted.
  
!  ========================================

! =========== ENVIRONMENT ====================
    integer, parameter, public :: ENVIRONMENT_SIZE = 100 !size of envornment, nr of cells = 100

    real, parameter, public :: FOOD_AVAILABILITY_MEAN = 0.40, FOOD_AVAILABILITY_VARIANCE = 0.25 !parameter for food, measured in weight grams added to the birds mass
    real, parameter, public :: FOOD_AVAILABILITY_MIN = 0.15, &
        FOOD_AVAILABILITY_MAX = FOOD_AVAILABILITY_MEAN + FOOD_AVAILABILITY_VARIANCE *1.5
! ===========================================


!! =========== OUTPUT ====================
    ! File name that keeps all output data from all generations of the model
    character(len=*), parameter, public :: MODEL_OUTPUT_FILE = "my_model_output_2G.csv"

    character(len=*), parameter, public :: FOOD_PRED_DISTRIBUTION_FILE = "food_pred_distribution_2G.csv"
! =========================================

!!=== CREATING DIRECTORIES FOR OUTPUT FILES =====

    character(*), parameter, public :: GENERATION_SUBFOLDER = "generations/"

    character(*), parameter, public :: WEIGHT_HISTORY_SUBFOLDER = "weight_history_2G/"

    character(*), parameter, public :: EMOTION_STATE_SUBFOLDER = "emotion_history_2G/"

!================================================

!=================== GA ===================
    !Population size
    integer, parameter, public :: GENERATIONS = 300 ! 30 for testing, 300 for sim. 
    integer, parameter, public :: SIGNAL_ONLY_GENERATIONS = 10 
    !Generations for when attack_rate is turned off and only signal remains
    
    integer, parameter, public :: EASY_GENERATIONS = 20 !10 for testing, 20 for sim

    integer, parameter, public :: POP_SIZE = 100000

    !Proportion of the best reproducing birds of selection
    real, parameter :: GA_REPRODUCE_PR = 0.25
    !exact number of birds reproducing
    integer, parameter, public :: GA_REPRODUCE_N = int(POP_SIZE * GA_REPRODUCE_PR)

    real, parameter, public :: GA_PROB_REPR_MAX = 0.9
    real, parameter, public :: GA_PROB_REPR_MIN = 0.05

    integer, parameter, public :: NUM_EGGS_RANGE_MAX = 9
    integer, parameter, public :: NUM_EGGS_RANGE_MIN = 1

    integer, parameter, public :: TOTAL_TIME_STEP = 365

    real, parameter, public :: MUTATION_PROBABILITY = 0.15
    real, parameter, public :: BACKGROUND_MORTALITY = 0.001 !if Is_Evolutionary_Generations = True
    !Background mortality is a common word birds killed by events as disease, hostile weather-conditions, parasitism, etc...
!==========================================

!==================== Global Variables ===========================
    integer, public :: Current_Generation 
    logical, public :: Is_Evolutionary_Generations 
    
    character(len=:), allocatable, public :: Genome_File_Name

    real, public :: Full_Pop_Factor 
    real, public :: Last_Full_Pop_Factor
    ! the real value one must multiply
    ! with survivors in order to return to POP_SIZE
    ! or size of the population at it's full. 



! ================================================================



    contains

  !-----------------------------------------------------------------------------
  !> Force a value within the range set by the vmin and vmax dummy parameter
  !! values. If the value is within the range, it does not change, if it
  !! falls outside, the output force value is obtained as
  !! min( max( value, FORCE_MIN ), FORCE_MAX )
  !! @param[in] value_in Input value for forcing transformation.
  !! @param[in] vmin minimum value of the force-to range (lower limit), if
  !!            not present, a lower limit of 0.0 is used.
  !! @param[in] vmax maximum value of the force-to range (upper limit)
  !! @returns   an input value forced to the range.
  !! @note      Note that this is the **real** precision version of the
  !!            generic `within` function.
  elemental function within(value_in, vmin, vmax) result (value_out)
    real, intent(in) :: value_in
    real, optional, intent(in) :: vmin
    real, intent(in) :: vmax
    real :: value_out

    ! Local copies of optionals.
    real :: vmin_here

    ! Check optional minimum value, if absent, set to a default value 0.0.
    if (present(vmin)) then
      vmin_here =  vmin
    else
      vmin_here = 0.0
    end if

    value_out = min( max( value_in, vmin_here ), vmax )

  end function within

  !-----------------------------------------------------------------------------
  !> @brief    Rescale a real variable with the range A:B to have the new
  !!           range A1:B1.
  !! @details  Linear transformation of the input value `value_in` such
  !!           `k * value_in + beta`, where the `k` and `beta` coefficients
  !!           are found by solving a simple linear system:
  !!           @f$  \left\{\begin{matrix}
  !!                A_{1}= k \cdot A + \beta;   \\
  !!                B_{1}= k \cdot B + \beta
  !!                \end{matrix}\right. @f$. It has this solution:
  !!           @f$ k=\frac{A_{1}-B_{1}}{A-B},
  !!               \beta=-\frac{A_{1} \cdot B-A \cdot B_{1}}{A-B} @f$
  !! @warning  The function does not check if `value_in` lies within [A:B].
  !! @note     Code for wxMaxima equation solve:
  !! @code
  !!           solve(  [A1=A*k+b, B1=B*k+b] ,[k,b] );
  !! @endcode
  elemental function rescale(value_in, A, B, A1, B1) result(rescaled)
    real, intent(in) :: value_in
    real, intent(in) :: A, B, A1, B1
    real :: rescaled

    ! Local variables
    real :: ck, cb

    !> ### Implementation details ###

    !> First, find the linear coefficients `ck` and `cb`
    !! from the simple linear system.
    ck = (A1-B1) / (A-B)
    cb = -1.0 * ((A1*B - A*B1) / (A-B))

    !> Second, do the actual linear rescale of the input value.
    rescaled = value_in*ck + cb

  end function rescale


  !-----------------------------------------------------------------------------
  !> Calculate an average value of a real array, excluding MISSING values.
  !! @param vector_in The input data vector
  !! @param missing_code Optional parameter setting the missing data code,
  !!        to be excluded from the calculation of the mean.
  !! @param undef_ret_null Optional parameter, if TRUE, the function returns
  !!        zero rather than undefined if the sample size is zero.
  !! @returns The mean value of the vector.
  !! @note This is a real array version.
  pure function average (array_in, missing_code, undef_ret_null)            &
                                                              result (mean_val)

    ! @param vector_in The input data vector
    real, dimension(:), intent(in) :: array_in

    ! @param missing_code Optional parameter setting the missing data code,
    !        to be excluded from the calculation of the mean.
    real, optional, intent(in) :: missing_code

    ! @param undef_ret_null Optional parameter, if TRUE, the function returns
    ! zero rather than undefined if the sample size is zero.
    logical, optional, intent(in) :: undef_ret_null

    ! @returns The mean value of the vector.
    real :: mean_val

    ! Local missing code.
    real :: missing_code_here

    ! Local sample size, N of valid values.
    integer :: count_valid

    ! Define high precision kind for very big value
    integer, parameter :: HRP = selected_real_kind(33, 4931)

    ! Big arrays can result in huge sum values, to accommodate them,
    ! use commondata::hrp and commondata::long types
    real(HRP) :: bigsum, bigmean

    !> ### Implementation details ###

    !> Check if missing data code is provided from dummy input.
    !! If not, use global parameter.
    if (present(missing_code)) then
      missing_code_here = missing_code
    else
      missing_code_here = MISSING
    end if

    !> Fist, count how many valid values are there in the array.
    count_valid = count(array_in /= missing_code_here)

    !> If there are no valid values in the array, mean is undefined.
    if (count_valid==0) then
      if (present(undef_ret_null)) then
        if (undef_ret_null) then
          mean_val = 0.0    !> still return zero if undef_ret_null is TRUE.
        else
          mean_val = MISSING
        end if
      else
        mean_val = MISSING
      end if
      return
    end if

    bigsum = sum( real(array_in, HRP), array_in /= missing_code_here )
    bigmean = bigsum / count_valid

    mean_val =  real(bigmean)

  end function average


!   pure function update_pop_regain_size_factor(pop_survived, pop_full) result(Full_Pop_Factor)
!   integer, dimension(:), intent(in) :: pop_survived
!   real, dimension(:), intent(out) :: pop_full ! Assuming pop_full is meant to be returned as well
!   real :: Full_Pop_Factor


  !> Calculates the factor to regain the full population size from the survived population.
  !!
  !! This function takes the number of individuals that survived and the full population size,
  !! and calculates the factor that can be used to scale up the survived population to the
  !! full population size.
  !!
  !! @param pop_survived The number of individuals that survived.
  !! @param pop_full The full population size.
  !! @return The factor to scale up the survived population to the full population size.
  pure function update_pop_regain_size_factor(pop_survived, pop_full) result(factor)
  integer, intent(in) :: pop_survived
  integer, intent(in) :: pop_full ! Assuming pop_full is meant to be returned as well
  real :: factor

    factor = real(pop_full) / real(pop_survived)


  end function update_pop_regain_size_factor


  !For finding probability of reproduction
  !> Calculates the probability of reproduction for an individual based on its current 
  !! mass (m) and the mass thresholds for reproduction (m_0 and m_max).
  !!
  !! This function uses a linear interpolation between the minimum and maximum 
  !! probabilities of reproduction (GA_PROB_REPR_MIN and GA_PROB_REPR_MAX) to determine the probability of reproduction for the given mass.
  !!
  !! @param m The current mass of the individual.
  !! @param m_0 The minimum mass threshold for reproduction.
  !! @param m_max The maximum mass threshold for reproduction.
  !! @return The probability of reproduction for the individual.

  function prob_repr(m, m_0, m_max) result(prob)

  real, intent(in) :: m, m_0, m_max
  real  :: prob
  real ::  k, b
  real, parameter :: P  = GA_PROB_REPR_MAX
  real, parameter :: P_0 = GA_PROB_REPR_MIN

  call solve_linear( k, b, m_0, P_0, m_max, P )

  prob = k * m + b

  prob = within(prob, P_0, P)

  end function prob_repr

  
  ! General solver for linear equation
  !
  ! Given the x_min, x_max and y_min y_max,
  !   determine the coefficients for the linear
  !   equation k and b
  ! y_max +       *
  !       |     * .    y = k x + b
  !       |   *   .        ?     ?
  !       | *     .
  ! y_min +-------+
  !     x_min   x_max
  
  !> Solves for the coefficients k and b of a linear equation y = kx + b, 
  !! given the minimum and maximum values of x and y.
  !!
  !! This subroutine takes the minimum and maximum values of x and y, 
  !! and calculates the coefficients k and b for the linear equation y = kx + b 
  !! that passes through those points.
  !!
  !! @param k The slope of the linear equation.
  !! @param b The y-intercept of the linear equation.
  !! @param x_min The minimum value of x.
  !! @param y_min The minimum value of y.
  !! @param x_max The maximum value of x.
  !! @param y_max The maximum value of y.
  pure subroutine solve_linear(k, b, x_min, y_min, x_max, y_max)
    real, intent(out) :: k,b
    real, intent(in)  :: x_min, y_min, x_max, y_max

    k = (y_min - y_max) / (x_min-x_max)
    b = -1 * (x_max*y_min - x_min*y_max) / ( x_min - x_max )

  end subroutine solve_linear



!  elemental function prob_repr_to_egg(prob) result(num_eggs)
!   real, intent(in) :: prob
!   integer :: num_eggs
!   real :: k, b
!   integer, parameter :: max_eggs = NUM_EGGS_RANGE_MAX 
!   integer, parameter :: min_eggs = NUM_EGGS_RANGE_MIN 
!   real, parameter :: prob_max = GA_PROB_REPR_MAX
!   real, parameter :: prob_min = GA_PROB_REPR_MIN

!   !y = kx + b
!   !y = num_eggs
!   !k = slope of line
!   !x = input val of repr_prob for current bird = prob
!   !b = y-intercept

!   !Setting this conditional due to wanting egg_prod to be shown
!   !even though probability is low. In this way, we show a more
!   !gradual increase in minimum weight for reproducing mothers.
   
!   !call solve_linear(k, b, prob_min, real(min_eggs), prob_max, real(max_eggs))

	 	
!   !Setting this conditional due to wanting egg_prod to be shown 	 	
!   !even though probability is low. In this way, we show a more 	 	
!   !gradual increase in minimum weight for reproducing mothers. 
  
!   !linear interpolation: To find b, we use one known point on the line. 	 	
!   ! This known point can be found by using minimum values known: min_eggs and 	 	
!   ! prob_min. 	 	

!    k = real(max_eggs - min_eggs) / (prob_max - prob_min)
!    b = real(min_eggs) - k * prob_min


!    num_eggs = (k * prob + b)
!    num_eggs =  nint(real(num_eggs) * Full_Pop_Factor)
 


!  end function prob_repr_to_egg


  elemental function prob_repr_to_egg(prob) result(num_eggs)
   real, intent(in) :: prob
   integer :: num_eggs
   real :: k, b
   integer, parameter :: max_eggs = NUM_EGGS_RANGE_MAX
   integer, parameter :: min_eggs = NUM_EGGS_RANGE_MIN
   real, parameter :: prob_max = GA_PROB_REPR_MAX
   real, parameter :: prob_min = GA_PROB_REPR_MIN

   !y = kx + b
   !y = num_eggs
   !k = slope of line
   !x = input val of repr_prob for current bird = prob
   !b = y-intercept

   if (prob == GA_PROB_REPR_MIN) then
     num_eggs = 0
   else
    k = (max_eggs - min_eggs) / (prob_max - prob_min)
  !linear interpolation:
    b = min_eggs - k * prob_min

    num_eggs = nint(k*prob + b) 
   end if

  end function prob_repr_to_egg


  !function and subsequent subroutine for calculating median value. 
  ! https://rosettacode.org/wiki/Averages/Median#Fortran
  function median (array_in) result(median_value)
     implicit none
     real, dimension(:), intent(in) :: array_in
  
     real, dimension(:), allocatable :: sorted_array
     real :: median_value
     integer :: i, n, mid 

     integer :: n_valid, index_valid


     ! make a smaller array from the input array that excludes missing values: 
     ! Make size of sorted array to 
 
     !Copy array to avoid modifying original data
     n_valid = 0

     do i=1, size(array_in)
       if (array_in(i) /= MISSING) n_valid = n_valid + 1
     end do

     allocate(sorted_array(n_valid))
     

     
     index_valid = 0 
     do i=1, size(array_in)
        if (array_in(i) /= MISSING) then 
          index_valid = index_valid + 1 
          sorted_array(index_valid) = array_in(i)
        end if
     end do
     
     n = size(sorted_array)

     call sort(sorted_array)

     !Calculate the median 
     
     if (mod(n, 2) == 0) then 
       mid = n / 2 
       median_value = (sorted_array(mid) + sorted_array(mid + 1)) / 2.0
     else
       mid = (n + 1) / 2
       median_value = sorted_array(mid)
     end if 

  end function median
  

  subroutine sort(array)
    implicit none
    real, intent(inout) :: array(:)
    integer :: n
    integer :: i, j
    real :: temp

    n = size(array)

  ! Making a bubble sort algorithm: 
  ! Bubble sort algorithm, also known as sinking sort,
  !  is the simplest sorting algorithm that runs through the list repeatedly,
  !   compares adjacent elements, and swaps them if they are out of order.
  ! Code inspired by https://rosettacode.org/wiki/Sorting_algorithms/Bubble_sort#Fortran 
  !
   do i = 1, n - 1
     do j = i + 1, n
       if (array(i) > array(j)) then
         temp = array(i)
         array(i) = array(j)
         array(j) = temp
       end if
     end do
   end do

  end subroutine sort


  
  !-----------------------------------------------------------------------------
  !> Calculate standard deviation using trivial formula:
  !! @f[ \sigma=\sqrt{\frac{\sum (x-\overline{x})^{2}}{N-1}} . @f]
  !! @note This is a real array version.
  function std_dev(array_in, missing_code, undef_ret_null) result (stddev)
    !> @param vector_in The input data vector
    real, dimension(:), intent(in) :: array_in
    !> @param missing_code Optional parameter setting the missing data code,
    !!        to be excluded from the calculation of the mean.
    real, optional, intent(in) :: missing_code
    !> @param undef_ret_null Optional parameter, if TRUE, the function returns
    !! zero rather than undefined if the sample size is zero.
    logical, optional, intent(in) :: undef_ret_null
    !> @returns The standard deviation of the data vector.
    real :: stddev

    ! Local missing code.
    real :: missing_code_here

    ! Local sample size, N of valid values.
    integer :: count_valid

    ! Minimum sample size resulting in real calculation, everythin less
    ! than this returns invalid.
    integer, parameter :: MIN_N = 4

    ! An array of squared deviations
    real, dimension(size(array_in)) :: dev2

    ! Mean value
    real :: sample_mean

    !> Check if missing data code is provided from dummy input.
    !! If not, use global parameter.
    if (present(missing_code)) then
      missing_code_here = missing_code
    else
      missing_code_here = MISSING
    end if

    count_valid = count(array_in /= missing_code_here)

    !> If there are no valid values in the array, std. dev. is undefined.
    if (count_valid <  MIN_N) then
      if (present(undef_ret_null)) then
        if (undef_ret_null) then
          stddev = 0.0    !> still return zero if undef_ret_null is TRUE.
        else
          stddev = MISSING
        end if
      else
        stddev = MISSING
      end if
      return
    end if

    sample_mean = average ( array_in, missing_code_here )

    where ( array_in /= missing_code_here )
      dev2 = ( array_in - sample_mean ) ** 2
    elsewhere
      dev2 = missing_code_here
    end where

    stddev = sqrt( sum( dev2, dev2 /= missing_code_here ) / (count_valid - 1) )

  end function std_dev

  !-----------------------------------------------------------------------------
  !> Converts logical to standard (kind SRP) real, .FALSE. => 0, .TRUE. => 1
  !! @note Note that this function is required to place logical data
  !!       like the survival status (alive) into the reshape array that is
  !!       saved as the CSV data file.
  !!       **Example: **
  !!       @code
  !!          call CSV_MATRIX_WRITE ( reshape(                              &
  !!                     [ habitat_safe%food%food%x,                        &
  !!                       habitat_safe%food%food%y,                        &
  !!                       habitat_safe%food%food%depth,                    &
  !!                       conv_l2r(habitat_safe%food%food%eaten),          &
  !!                       habitat_safe%food%food%size],                    &
  !!                     [habitat_safe%food%number_food_items, 5]),         &
  !!                     "zzz_food_s" // MODEL_NAME // "_" // MMDD //       &
  !!                       "_gen_" // TOSTR(generat, GENERATIONS) // csv,   &
  !!                     ["X   ","Y   ", "D   ", "EATN", "SIZE"]            &
  !!                     ) @endcode
  elemental function l2r(flag, code_false, code_true) result (num_out)
    logical, intent(in) :: flag
    real, optional, intent(in) :: code_false
    real, optional, intent(in) :: code_true
    real :: num_out
    ! Local copies of optionals.
    real :: code_false_loc, code_true_loc
    ! Local real parameters for the default FALSE and TRUE.
    real, parameter :: FALSE_DEF=0.0, TRUE_DEF=1.0

    !> First, check optional parameters.
    if (present(code_false)) then
      code_false_loc = code_false
    else
      code_false_loc = FALSE_DEF
    end if

    if (present(code_true)) then
      code_true_loc = code_true
    else
      code_true_loc = TRUE_DEF
    end if

    ! Second, do the actual conversion.
    if (flag .eqv. .TRUE.) then
      num_out = code_true_loc
    else
      num_out = code_false_loc
    end if
  end function l2r



! Function for calculating linear increase in predation to ensure soft pressure from the environment in the initial 20 generations.

  function approx_easy(generation, min_y, max_y) result (y_approx)
    real :: y_approx
    integer, intent(in) :: generation
    real, intent(in) :: min_y, max_y
    !when M is at it's highest value, the attack_rate is also at it's highest,
    !when M is lower than it's highest value, the attack rate is somewhere along
    !the linear slope that goes from minimum attack rate (min_y) to maximum attack rate
    !(max_y).

    real :: k, b ! parameters of the linear equation

    integer, parameter :: M = EASY_GENERATIONS !Parameter for the generations where the environment is less harsh and allows for an individual's mistakes.
    

     k = -1 * (min_y - max_y) / (M - 1)
     b = (M*min_y - max_y) / (M - 1)
  
     if (generation < M) then !if the generation nr is less than M(=20)
       y_approx = k * real(generation) + b !y_approx, the environments predation is k * (real)generation.
     else                                  !"Real" due to temporal continuity within and between generation (slope of the linear eq)
       y_approx = max_y !if generation is 20 or above, continue with universal parameters.
     end if


  end function approx_easy




  subroutine get_global_runtime_params()
! Purpose
!   The subroutine get_global_runtime_params is designed to retrieve and 
! process command-line arguments passed to the program. It sets global 
! parameters based on the presence and values of these arguments.
! Description:
!   This subroutine fetches the number of command-line arguments 
! and stores them in an allocatable array. It then processes the 
! first argument, if available, to set the values of two global 
! parameters: Genome_File_Name and Is_Evolutionary_Generations.
! Usage:
!   This subroutine does not take any input arguments nor return any output. 
! It relies on command-line arguments passed when the Fortran program is executed.
! Example, in terminal, if we run gfortran -o example example.f90:
    !./example
    ! expected output should be:      Is_Evolutionary_Generations:  T
    !if we run ./example my_genome_file.txt
    ! expected output:                Is_Evolutionary_Generations:  F

  integer :: i, n_args
  character(100), dimension(:), allocatable :: cmd_args

  ! Get the number of command-line arguments
  n_args = command_argument_count()
  
  ! Allocate array to hold command-line arguments
  allocate(cmd_args(n_args))
  
  ! Retrieve each command-line argument and store it in cmd_args
  do i = 1, n_args
    call get_command_argument(number = i, value = cmd_args(i))
  end do

  ! Process the command-line arguments to set global parameters
  if (n_args > 0) then
    Genome_File_Name = trim(cmd_args(1))
    Is_Evolutionary_Generations = .FALSE.
  else
    Genome_File_Name = " "
    Is_Evolutionary_Generations = .TRUE.
  end if

end subroutine get_global_runtime_params

  

end module params

!------------------------------------
!------ENVIRONMENT MODULE------------
!------------------------------------

module environment

use params
implicit none

! Spatial location at one cell, in the initial model it is
! just a single cell, integer number
! Note: defining an elementary spatial unit makes it easy to
!       extend the model to 2d or 3d
type, public :: location
    integer :: x
    contains
    procedure, public :: put_random => location_place_random
    procedure, public :: place => location_place_object_location !for bird and predator
    procedure, public :: walk => location_random_walk !only for bird


end type location

! Each spatial cell has additional data characteristics:
! - food availability
! - predation risk
type, extends(location) :: env_cell
  real :: food_availability
  real :: predator_frequency
  contains
  procedure, public :: init => cell_init
end type env_cell

! The whole environment is a collection of cells
! Note: note that the size of the environment is an allocatable (dynamic)
!       array.
type whole_environ
  type(env_cell), allocatable, dimension(:) :: point
  contains
  procedure, public :: init => env_init
  procedure, public :: save => environemnt_save_csv
  procedure, public :: load_env_csv => load_environment_from_csv

end type whole_environ

contains

! This subroutine places the basic spatial object (location class) to a
! random place.
subroutine location_place_random(this, min_pos, max_pos)
    use BASE_RANDOM
    class(location), intent(inout) :: this
    ! Optional parameters defining a range of positions to place the spatial
    ! object, but by default the position is from 1 to ENVIRONMENT_SIZE
    integer, optional, intent(in) :: min_pos, max_pos

    integer :: min_loc, max_loc

    if (present(min_pos)) then
        min_loc = min_pos
    else
        min_loc = 1
    end if

    if (present(max_pos)) then
        max_loc = max_pos
    else
        max_loc = ENVIRONMENT_SIZE
    end if

    this%x = RAND(min_loc, max_loc)

end subroutine location_place_random

! Place a spatial object to a specific location (cell) within the environment
subroutine location_place_object_location(this, where)
    class(location), intent(inout) :: this
    integer, intent(in) :: where

    this%x = where

end subroutine location_place_object_location



subroutine location_random_walk(this, min_pos, max_pos)
  use BASE_RANDOM
  class(location), intent(inout) :: this

  integer, optional, intent(in) :: min_pos, max_pos

  !deside if left or right-movement
  logical :: is_going_left

  integer :: min_loc, max_loc

  ! Process optional parameters using the min_pos and max_pos from
  ! the subroutine call or default values
  if (present(min_pos)) then
        min_loc = min_pos
  else
        min_loc = 1
  end if

  if (present(max_pos)) then
        max_loc = max_pos
  else
        max_loc = ENVIRONMENT_SIZE
  end if

  if ( RAND() > 0.5 ) then
    is_going_left = .TRUE.
  else
    is_going_left = .FALSE.
  end if

  if ( is_going_left .and. this%x == min_loc ) then
    is_going_left = .FALSE.
  elseif ( .NOT. is_going_left .and. this%x == max_loc) then
    is_going_left = .TRUE.
  end if

  if (is_going_left) then
    this%x = this%x - 1
  else
    this%x = this%x + 1
  end if

end subroutine location_random_walk

! Initialize a single cell within the environment
subroutine cell_init(this, x, env)
  use BASE_RANDOM
  class(env_cell), intent(out) :: this
  integer, intent(in) :: x
  type(whole_environ), optional, intent(inout) :: env

  this%x = x ! set location
  this%food_availability = within(RNORM(FOOD_AVAILABILITY_MEAN, FOOD_AVAILABILITY_VARIANCE), &
                                  FOOD_AVAILABILITY_MIN, FOOD_AVAILABILITY_MAX)
  if (Is_Evolutionary_Generations) then
      this%predator_frequency = RAND(FREQUENCY_OF_PREDATOR_MIN, FREQUENCY_OF_PREDATOR_MAX)
  else
      ! if (present(env)) then
      !     call env%load_env_csv()
      !     this%predator_frequency = env%point(x)%predator_frequency
      !  else
    this%predator_frequency = RAND(FREQUENCY_OF_PREDATOR_MIN * PRED_FREQUENCY_MULTIPLIER_SIGNAL_GENS_ONLY,     &
        FREQUENCY_OF_PREDATOR_MAX * PRED_FREQUENCY_MULTIPLIER_SIGNAL_GENS_ONLY)
      ! end if
      ! this%predator_frequency = this%predator_frequency * PRED_FREQUENCY_MULTIPLIER_SIGNAL_GENS_ONLY
      !  end if
  end if
end subroutine cell_init




! Initialize the whole environment
subroutine env_init(this, max_size) !max size of the environ, but optional. Made for testing with shorter arrays.
  class(whole_environ), intent(inout) :: this
  integer, intent(in), optional :: max_size

  integer :: i, max_size_loc

  ! process optional parameter defining the maximum size of the environment
  if (present(max_size)) then
    max_size_loc = max_size
  else
    max_size_loc = ENVIRONMENT_SIZE
  end if

  ! First, allocate the dynamic environment array with its size
  ! Note that we test if the array has already been allocated to guard against
  !      possible error, e.g. if init is called more than once
  if (.not. allocated(this%point) ) then
    allocate (this%point(max_size_loc))
  else
    ! Repeated initialization of the environment
    ! would be strange so we report this
    write(*,*) "WARNING: repeated initialization of the environment detected"
    deallocate (this%point)
    allocate (this%point(max_size_loc))
  end if

  do i=1, size(this%point)
    call this%point(i)%init(i, this) ! we initialize the food and predation data
  end do

end subroutine env_init



!Subroutine for saving info on food-distribution across the one-dimensional matrix that is our birds array of habitats.
!Subroutine logic:
!Using only CSV_ARRAY_WRITE from CSV_IO, an array is created. This array makes an output array(output_array)
!which is basically a one dimensional array of cells that correspond the habitats our bird walks in. This array
!consist of however many cells we want ("this"). When we loop over it from i = 1 until total array size, for each loop
!the food that is within each cell is recorded. This only needs to be done once, since the food in cells are a non-changing
!value. This is then saved to output_array, that gets sent to FOOD_PRED_DISTRIBUTION_FILE
subroutine environemnt_save_csv(this)
  use CSV_IO, only : CSV_MATRIX_WRITE
  class(whole_environ), intent(in) :: this
  integer  :: i
  real, allocatable, dimension(:,:) :: output_array


  allocate(output_array(size(this%point), 2))

  do i=1, size(this%point)
   output_array(i,1) = this%point(i)%food_availability
   output_array(i,2) = this%point(i)%predator_frequency
  end do

  call CSV_MATRIX_WRITE(output_array, FOOD_PRED_DISTRIBUTION_FILE, colnames = ["FOOD_AVAILABILITY", "RISK_PREDATION   "])

end subroutine environemnt_save_csv



!> Loads the environment data from a CSV file and initializes the environment cells.
!>
!> If the CSV file is successfully opened, the function reads the food availability and predator frequency
!> for each cell in the environment and stores them in the corresponding `this%point` elements. If the file
!> does not contain data for all cells, the remaining cells are initialized using the `cell_init` subroutine.
!>
!> If the CSV file cannot be opened, the function generates a new environment by initializing all cells
!> using the `cell_init` subroutine.
!>
!> @param this The `whole_environ` object representing the environment.
subroutine load_environment_from_csv(this)
  class(whole_environ), intent(inout) :: this
  integer :: i, io_status
  character(len=100) :: line
  
  open(unit=10, file=FOOD_PRED_DISTRIBUTION_FILE, status='old', action='read', iostat=io_status)
  
  if (io_status == 0) then
      read(10, *, iostat=io_status) ! Skip header line
      
      do i = 1, size(this%point)
          read(10, *, iostat=io_status) this%point(i)%food_availability !, this%point(i)%predator_frequency
          if (io_status /= 0) exit
          this%point(i)%x = i
      end do
      
      close(10)
      
      if (i <= size(this%point)) then
          print*, "Warning: Not all environment data loaded. Generating remaining cells."
          do i = i, size(this%point)
              call cell_init(this%point(i), i)
          end do
      end if
  else
      print*, "Error opening CSV file. Generating new environment."
      do i = 1, size(this%point)
          call cell_init(this%point(i), i)
      end do
  end if
end subroutine load_environment_from_csv

end module environment



!--------------------------------------
!-------GENOME MODULE------------------
!--------------------------------------
module organism

use params
use environment
use BASE_RANDOM
implicit none

!describe how we define genes
!make an integer range that describes fear vs hunger

!--------------------------------
!-----------LOCATION/GENE--------
!--------------------------------
type, extends(location), public :: GENOME
    integer :: gene_fear !the higher the gene_fear the lower the gene, opposite proportional
    integer :: gene_hunger !DOCUMENT ME! 
    contains
    procedure, public :: init_genome => genome_init_all
    procedure, public :: mutate => gene_mutate


end type GENOME


!--------------------------------
!-------------BIRD/GENOME----------
!--------------------------------
type, extends(GENOME), public :: BIRD
    real :: weight  !In this model set by the params BIRD_INITAL_WEIGHT and BIRD_MAXIMUM_WEIGHT_MULTIPLICATOR/BIRD_MINIMUM_WEIGHT_MULTIPLICATOR
    real, dimension(TOTAL_TIME_STEP) :: weight_history !recording weight in matrix by timestep and generation
    real, dimension(TOTAL_TIME_STEP) :: fear_history ! recording emotion in matrix by timestep and generation - previously emotion_history
    real, dimension(TOTAL_TIME_STEP) :: hunger_history
    real :: state_fear !values between 0 and 10
    real :: state_hunger !values between 0 and 10
    logical :: is_alive 
    logical :: is_killed_by_predator
    integer :: bird_meets_predator_counter
    


    !integer :: death_count !this relates to subroutine count_the_dead_birds.
    contains
    procedure, public :: init => bird_init  !updated
    procedure, public :: killed_from_starvation => bird_dead_from_starvation !updated
    procedure, public :: is_starved => bird_is_starved  !updated
    procedure, public :: fly => bird_do_fly ! updated
    procedure, public :: do_feed => bird_feeds !updated
    procedure, public :: do_explore => bird_do_explore

    procedure, public :: is_killed_by_background => bird_dead_from_background_mortality !NEW! SUBROUTINE
    procedure, public :: environment_kills => background_mortality_probability !NEW! FUNCTION

    procedure, public :: genetic_lower_hunger => bird_hunger_min_genetic_limit !NEW!
    procedure, public :: genetic_upper_hunger => bird_hunger_max_genetic_limit !NEW!

    procedure, public :: genetic_lower_fear => bird_fear_min_genetic_limit   !NEW!
    procedure, public :: genetic_upper_fear => bird_fear_max_genetic_limit   !NEW!

    procedure, public :: hunger_genetic_limit => bird_force_hunger_within_genetic_limits !updated!
    procedure, public :: fear_genetic_limit => bird_force_fear_within_genetic_limits   !NEW!

    procedure, public :: add_to_history => bird_weight_add_to_history !UPDATED!
    procedure, public :: emo_state_to_history => bird_emotion_state_add_to_history !UPDATED!
    procedure, public :: pay_for_life => bird_subtract_metabolic_cost  !Keep as is!

    procedure, public :: is_within => bird_is_within_environment_check !checking if bird is within the environment with correct number of cells 

end type BIRD

!--------------------------------
!-----------PREDATOR-------------
!--------------------------------
type, public :: PREDATOR
    real :: risk !If PREDATOR sees the bird according to the risk_of_the, the risk of attack is this


  contains

    procedure, public :: init => predator_init_new !keep as is!

    procedure, public :: predator_attack_bird

end type PREDATOR


!--------------------------------
!-----------POPULATION-----------
!--------------------------------

type, public :: POPULATION
!population size:
  type(BIRD), dimension(:), allocatable :: birds
  integer :: num_dead !global counter of dead birds
  logical :: is_alive

!any other properties of the population?
!-----
  contains
  !subroutines and functions that apply for the population
  procedure, public :: init_pop => population_init_new
  procedure, public :: dead_count => population_get_dead_count
  procedure, public :: alive_count => population_get_alive_count
  procedure, public :: time_steps => population_time_steps
  procedure, public :: sort => sort_by_fitness
  procedure, public :: save => population_save_csv
  procedure, public :: load => population_load_genome_csv
  procedure, public :: smallest_weight => population_minimum_weight_alive
  procedure, public :: biggest_weight => population_maximum_weight_alive





end type POPULATION



!--------------------------------
!---------SUBROUTINES-----------
!--------------------------------
contains
!subroutines for POPULATION


!----------------------------------------

! The "population_init_new" subroutine initializes a population of birds within a POPULATION class object. 
! It allows for the optional specification of an alternative population size (alt_size) 
! and an alternative gene initialization (is_gene_init) for testing and debugging purposes.

! "this": A POPULATION class object passed by reference (intent(out)). This parameter is the target for the population initialization.
! "alt_size": An optional integer parameter (intent(in)) that specifies an alternative population size. 
! If not provided, the default population size (POP_SIZE) defined in the parameters module is used.

! "is_gene_init": An optional logical parameter (intent(in)) that indicates whether to use an alternative gene initialization. 
! If not provided, the default value is .FALSE..!


subroutine population_init_new(this, alt_size, is_gene_init)
  class(POPULATION), intent(out) :: this! intent(out) because the initialization will be only an output.
  integer, optional, intent(in) :: alt_size! if we want some other value than POP_SIZE
  logical, optional, intent(in) :: is_gene_init !alternative gene value used for testing and debugging

  integer :: size
  integer :: i

  logical :: is_gene_init_here

  if (present(is_gene_init)) then
    is_gene_init_here = is_gene_init
  else
    is_gene_init_here = .FALSE.
  end if

  !process for optional parameter where we can define alternative pop_size:
  if (present(alt_size)) then
    size = alt_size
  else
    size = POP_SIZE !POP_SIZE defined in parameters mod.
  end if

  allocate (this%birds(size)) !allocate the array to a specific size

  do i = 1, size
    call this%birds(i)%init(is_gene_init_here) !initialize array of birds

  end do

end subroutine population_init_new

! Creating a csv-file with data under the column names seen below. 
subroutine population_save_csv(this, file_name) 
  use CSV_IO
  use params
  


  class(POPULATION), intent(in) :: this
  character(len=*), intent(in) :: file_name
  

  integer, parameter :: GENPOP_COL_NUMBER = 10
  character(len=*), dimension(GENPOP_COL_NUMBER), parameter :: GENPOP_COLUMNS = &
    ["BIRD                 ",                      &  ! 1
     "GENE_FEAR            ",                      &  ! 2
     "GENE_HUNGER          ",                      &  ! 3
     "IS_ALIVE             ",                      &  ! 4 !0 = FALSE, 1 = TRUE
     "WEIGHT               ",                      &  ! 5
     "BIRD_MEETS_PRED_COUNT",                      &  ! 6
     "STATE_FEAR           ",                      &  ! 7
     "STATE_HUNGER         ",                      &  ! 8
     "REPRODUCTION_PROB.   ",                      &  ! 9
     "NUMBER OF EGGS       "]                        ! 10


  real, dimension(POP_SIZE, GENPOP_COL_NUMBER) :: out_data

  integer :: i

   do i = 1 ,POP_SIZE
     out_data(i,1) = i
     out_data(i,2) = real(this%birds(i)%gene_fear)
     out_data(i,3) = real(this%birds(i)%gene_hunger)
     out_data(i,4) = l2r(this%birds(i)%is_alive) !0 = FALSE, 1 = TRUE
     out_data(i,5) = this%birds(i)%weight
     out_data(i,6) = real(this%birds(i)%bird_meets_predator_counter)
     out_data(i,7) = real(this%birds(i)%state_fear)
     out_data(i,8) = real(this%birds(i)%state_hunger)
     out_data(i,9) = (prob_repr(this%birds(i)%weight, this%smallest_weight(), this%biggest_weight()))
     out_data(i,10)= prob_repr_to_egg(out_data(i,9))
   end do



  call CSV_MATRIX_WRITE(out_data, file_name, colnames=GENPOP_COLUMNS)

end subroutine population_save_csv




subroutine population_load_genome_csv(this, file_name)
  use CSV_IO
  class(POPULATION), intent(inout) :: this
  character(len=*), intent(in) :: file_name


  logical, parameter :: IS_KILL_BIRDS_DEAD_IN_FILE = .FALSE. 
  !Start with false, TRUE if time for further experimentation
  ! if true, all the birds from the data are alive. 
   

  integer, parameter :: GENPOP_COL_NUMBER = 10
  integer, parameter :: FEAR_GENE_COL = 2, HUNGER_GENE_COL = 3, ALIVE_COL = 4
  real, parameter :: WEIGHT_COL = 5

  logical :: is_success

  integer :: inputdata_rows, inputdata_cols, i
  
  real, allocatable, dimension(:,:) :: in_data 

  in_data = CSV_MATRIX_READ(file_name, is_success, .FALSE.) !see documentation at AHAmodel https://ahamodel.uib.no/
  !errorflag is by default .TRUE. => everything works as it should. 

  inputdata_rows = size(in_data, 1)
  inputdata_cols = size(in_data, 2)

  
  if (.NOT. is_success) then
    print*, "ERROR: file error for ", file_name
    stop
  end if



  
  if(size(in_data,1) /= POP_SIZE) then 
   print*, "WARNING: input data wrong N of rows: n", size(in_data,1),      &
   ";Population size = ", POP_SIZE
  end if

  if (inputdata_cols /= GENPOP_COL_NUMBER) then 
   print*, "WARNING: input data wrong N of columns: n", inputdata_cols,      &
   ", must be ", GENPOP_COL_NUMBER
   stop
  end if

  this%birds(1:inputdata_rows)%gene_fear = int(in_data(:,FEAR_GENE_COL))
  this%birds(1:inputdata_rows)%gene_hunger = int(in_data(:,HUNGER_GENE_COL))
  this%birds(1:inputdata_rows)%weight = real(in_data(:,WEIGHT_COL))
  
  !if the data shows that the bird is dead by value 0, then we kill the bird again
  ! due to the bird being initialized as alive here. Thereby making sure that dead
  ! birds from original simulation stay dead in reading csv file again. 


  if (IS_KILL_BIRDS_DEAD_IN_FILE) then
   do i=1, inputdata_rows
    if (in_data(i, ALIVE_COL) == 1.0) this%birds(i)%is_alive = .FALSE.
   end do
  end if

  this%birds(inputdata_rows+1 : POP_SIZE)%is_alive = .FALSE.  
  !make code to exclude dead birds 
  !make warning if incorrect data is running 

end subroutine population_load_genome_csv





! Functions for counting dead birds, updating it. Iterating over the array.
! 
function population_get_dead_count(this) result(get_dead_count)! this = population, num_dead is the number of dead birds.
!this subroutine initialize a population of birds within a POPULATION class object. It allows
!for the optional specification of an alternative population size (alt_size), and an alternative
!gene initialization (is_gene_init) for testing and debugging.

  class(POPULATION), intent(in) :: this
  integer :: get_dead_count

  get_dead_count = count(.not.this%birds%is_alive) !birdS = population of birds.

end function population_get_dead_count


function population_get_alive_count(this) result(get_alive_count)! this = population, num_dead is the number of dead birds.
!this subroutine initialize a population of birds within a POPULATION class object. It allows
!for the optional specification of an alternative population size (alt_size), and an alternative
!gene initialization (is_gene_init) for testing and debugging.

  class(POPULATION), intent(in) :: this
  integer :: get_alive_count

  get_alive_count = count(this%birds%is_alive) !birdS = population of birds.

end function population_get_alive_count




function population_minimum_weight_alive(this) result(min_size)
  class(POPULATION), intent(in) :: this
  real :: min_size
  integer :: i 
  
  !IF time try to make more generic 


  i = this%alive_count()
  min_size = this%birds( i )%weight
  
end function population_minimum_weight_alive




function population_maximum_weight_alive(this) result(max_size)
  class(POPULATION), intent(in) :: this
  real :: max_size

  max_size = this%birds( 1 )%weight
  
end function population_maximum_weight_alive


! Subroutines for PREDATOR
subroutine predator_init_new(this, generation)

! This subroutine sets a condition for the attack rate: If the time past is less than 20
! generations, the risk is lower - somewhere on a linear slope - towards it's max which is in
! turn defined by ATTACK_RATE parameter. The conditions given are defined by the function aprox_y.

  class(PREDATOR), intent(inout) :: this
  integer, intent(in) :: generation
  real :: attack_rate_signal
    
    if (Is_Evolutionary_Generations) then 
      this%risk = approx_easy(generation, 0.0, ATTACK_RATE)
    else 
      attack_rate_signal = PRED_FREQUENCY_MULTIPLIER_SIGNAL_GENS_ONLY
      this%risk = attack_rate_signal
      
    end if
   
end subroutine predator_init_new



function bird_is_within_environment_check(this, environment_in) result(is_within)

  class(BIRD), intent(in) :: this
  class(whole_environ), intent(in) :: environment_in
  logical :: is_within
  
  integer :: env_size

  env_size = size(environment_in%point)

  if (this%location%x < env_size) then
    is_within = .TRUE. 
  else
   is_within = .FALSE. 
  end if
  

end function bird_is_within_environment_check




subroutine predator_attack_bird(this, bird_prey, environment_in,              &
                                predator_is_present, prey_is_killed)

 class(PREDATOR), intent(in) :: this
 class(BIRD), intent (inout) :: bird_prey
 class(whole_environ), intent(in) :: environment_in

 logical, optional, intent(out) :: predator_is_present
 logical, optional, intent(out) :: prey_is_killed


 logical :: p_is_present, p_prey_dies


 p_is_present = .FALSE.
 p_prey_dies = .FALSE.


 if (bird_prey%is_alive) then
   
   !Commenting out below lines due to bird being avoidant of predator, and therefore not entering "environment_in". 
   !Warning is mistakingly triggered due to this. 

     if (IS_DEBUG_DATA .AND. .NOT. bird_prey%is_within(environment_in)) then
       print*, "WARNING: Bird not in correct environment" 
     end if

    ! Check if the predator is present in the cell, given the risk of predation.
    p_is_present = (RAND() < environment_in%point( bird_prey%location%x )%predator_frequency)
    
    
    if (p_is_present) then
      ! We check if the predator will attack given it is in the same cell as the prey.
      if (Is_Evolutionary_Generations) p_prey_dies = ( RAND() < this%risk )

      if (p_prey_dies) then
        bird_prey%is_alive = .FALSE.
        bird_prey%bird_meets_predator_counter = 0
        bird_prey%is_killed_by_predator = .TRUE.
      else
        bird_prey%is_alive = .TRUE.
        bird_prey%bird_meets_predator_counter = bird_prey%bird_meets_predator_counter + 1
        bird_prey%weight = bird_prey%weight - bird_prey%weight *                                &
        (WEIGHT_REDUCTION_CONSTANT**(bird_prey%bird_meets_predator_counter))
        bird_prey%state_fear = within(                                                             &
          bird_prey%state_fear + bird_prey%state_fear * FEAR_INCREMENT_AFTER_ATTACK,    &
                                 bird_prey%genetic_lower_fear(), bird_prey%genetic_upper_fear() )
      end if
    end if
 end if

 if (present(predator_is_present)) predator_is_present = p_is_present
 if (present(prey_is_killed)) prey_is_killed = p_prey_dies


end subroutine predator_attack_bird



!****************************
!subroutine for initial gene
!****************************

!subroutine for initializing gene.
!The initial gene is randomly generated from the parameters HUNGER_MIN, HUNGER_MAX, set to -10 and +10.
!The gene can be exemplified as:
! -3: A gene coding for small to moderate degree of fear (or more prone to engage fear-like response in threat-detection ).
! 7 : A gene coding for moderate to high degree of boldness, but also same degree of hunger. in effect, a hungry bird not prone to fleeing easily.
subroutine genome_init_all(this, gene_num)
    use BASE_RANDOM !importing this.
    class(GENOME), intent(out) :: this
    integer, optional, intent(in) :: gene_num !gene number made for gene seperately initialized in gene_mutate(), so that genes are
                                              ! mutated seperately, independently of each other. 
    integer :: gene_num_loc

    if (present(gene_num)) then
      gene_num_loc = gene_num
    else
      gene_num_loc =  0
    end if

    if (gene_num_loc == 1) then
      this%gene_fear = RAND(FEAR_MIN, FEAR_MAX)       ! 1
    else if (gene_num_loc == 2) then
      this%gene_hunger = RAND(HUNGER_MIN, HUNGER_MAX)     ! 2
    else
    !classname%attribute_name
    this%gene_fear = RAND(FEAR_MIN, FEAR_MAX)       ! 1
    this%gene_hunger = RAND(HUNGER_MIN, HUNGER_MAX)     ! 2
    end if
end subroutine genome_init_all



! The gene given to each individual bird can be mutated. Since we are working with a single gene, the overall chance of mutations are high (set by parameter "MUTATION_PROBABILITY").
! altered randomly (by using function RAND()).
subroutine gene_mutate(this, is_mutate)
  class(GENOME), intent(inout) :: this
  logical, optional, intent(out) :: is_mutate !this local parameter is useful for debugging and testing purposes. 

  logical :: is_mutate_loc

  is_mutate_loc = .FALSE. ! The standard logical value for mutaten is "false", if mutation occurs by the logic explained below then it will be "true". 

  if (RAND() < MUTATION_PROBABILITY) then ! if probability for mutation is higher than the randomly generated value (RAND()), then mutation of the gene occurs. 
    call this%init_genome(1)
    is_mutate_loc = .TRUE.
  end if

  if (RAND() < MUTATION_PROBABILITY) then ! if probability for mutation is higher than the randomly generated value (RAND()), then mutation of the gene occurs. 
    call this%init_genome(2)
    is_mutate_loc = .TRUE.
  end if

  if (present(is_mutate)) is_mutate = is_mutate_loc

end subroutine gene_mutate



subroutine bird_init(this, is_gene_init)

!The bird_init subroutine initializes a BIRD class object with default attributes and optionally initializes its genetic attributes.
!It sets various properties of the bird, such as its weight, fear/hunger state, survival status, and other relevant attributes.

!this: A BIRD class object passed by reference (intent(out)). This parameter represents the bird to be initialized.
!is_gene_init: An optional logical parameter (intent(in)) that indicates whether to initialize the bird's genetic attributes. 
!  If not provided, the default value is .FALSE..


  class(BIRD), intent(out) :: this !by calling class(BIRD), we inherit all qualities of BIRD and supply it with new atributes through new subrout.
  logical, optional, intent(in) :: is_gene_init

  logical :: is_gene_init_here

  if (present(is_gene_init)) then
    is_gene_init_here = is_gene_init
  else
    is_gene_init_here = .FALSE.
  end if

  if (is_gene_init_here) then
    call this%init_genome()
   ! print*, "GENE INIT"
  end if



  call this%put_random()
  this%weight = BIRD_INITIAL_WEIGHT
  this%weight_history = MISSING !When initial weight is = 20 grams, the whole history of the array = MISSING, until the value is filled by the time_step, 
                                !revealing the development of the bird by time_step and generations.
  this%state_fear = real(this%gene_fear)
  this%state_hunger = real(this%gene_hunger)
  this%is_alive = .TRUE.
  this%is_killed_by_predator = .FALSE. 
  this%bird_meets_predator_counter = 0
  this%fear_history = MISSING
  this%hunger_history = MISSING


end subroutine bird_init


function bird_is_starved(this, threshold) result (is_dead) 

!The bird_is_starved function checks if a BIRD class object has died of starvation by comparing its weight to a specified threshold. 
!If the bird's weight drops below the threshold, the function returns .TRUE., indicating that the bird has died of starvation.


!this: A BIRD class object passed by reference (intent(in)). This parameter represents the bird whose starvation status is to be determined.
!threshold: An optional real parameter (intent(in)) that specifies the weight threshold below which the bird is considered to have died of starvation. 
!If not provided, the default threshold is calculated as BIRD_INITIAL_WEIGHT * BIRD_MINIMUM_WEIGHT_MULTIPLICATOR.

  class(BIRD), intent(in) :: this
  real, optional, intent(in) :: threshold
  logical :: is_dead

  real :: local_threshold

  if (present(threshold)) then
    local_threshold = threshold
  else
    local_threshold = BIRD_INITIAL_WEIGHT * BIRD_MINIMUM_WEIGHT_MULTIPLICATOR
  end if

  if (this%weight < local_threshold) then
    is_dead = .TRUE.
  else
    is_dead = .FALSE.
  end if

end function bird_is_starved





function background_mortality_probability(this) result(is_killed_by_background_mortality)
 class(BIRD), intent(inout) :: this
 logical :: is_killed_by_background_mortality

 is_killed_by_background_mortality = .FALSE.

 if (RAND() < BACKGROUND_MORTALITY) then
   is_killed_by_background_mortality = .TRUE. 
 end if
 
end function background_mortality_probability




subroutine bird_feeds(this, in_environment)
!The bird_feeds subroutine simulates a bird feeding by adjusting its weight and fear/hunger state based on the food availability 
!in its current cell and its emotional state. The bird's feeding behavior is influenced by its fear and hunger levels, 
!with the subroutine accounting for how these emotions affect the bird's feeding decisions.


!this: A BIRD class object passed by reference (intent(inout)). This parameter represents the bird that is feeding.
!in_environment: A whole_environ class object passed by reference (intent(in)). This parameter represents the environment 
!  in which the bird is located, including the food availability in the bird's current cell.

  class(BIRD), intent(inout) :: this
  class(whole_environ), intent(in) :: in_environment
  real :: current_cell_food

 if (this%is_alive) then 
  current_cell_food = in_environment%point(this%location%x)%food_availability
  this%weight = this%weight + current_cell_food   
  this%state_hunger = this%state_hunger - this%state_hunger * hunger_decrement(current_cell_food)
  call this%hunger_genetic_limit()
 end if
  
end subroutine bird_feeds




function hunger_decrement(food_eaten) result (decrement_out)
 real, intent(in) :: food_eaten
 real :: decrement_out
 
 !Function describes how much food is eaten according to hunger value in bird. 
 ! If bird has hunger = 10, it will eat all the food in the environment and therefore 
 ! reduce it by decrement_out

 ! if HUNGER_MAX_DECREMENT = 0.5 then: 
!  0.5 +      *
!      !    * |
!      !  *   |
!      !*     |
!      +------+-------------
!             FOO_AVAILABILITY_MAX

 decrement_out = HUNGER_MAX_DECREMENT/FOOD_AVAILABILITY_MAX * food_eaten

 end function hunger_decrement




function bird_eat_fraction_when_fear(emotion, min, max) result(fract)

! Calculate the proportion of food availability eaten at specific state of
! fear (state_hunger<0)
! The proportion of food eaten is 0.0 at state_fear = HUNGER_MIN
!                             and 1.0 at state_fear = HUNHER_MAX
!
!  S = k * emotion + b
!
!   | S1 = k * H_MIN + b
!   | S2 = k * M_MAX + b
!
!   solution:   k = (S2 - S1) / (H_MIN - H_MAX)
!               b = (S2 * H_MIN - S1 * H_MAX) / ( H_MIN-H_MAX )
!

!The bird_eat_fraction_when_fear function calculates the fraction of food a bird eats when it is experiencing fear, 
!based on its emotional state and optional minimum and maximum values for the emotional state. 
!The function uses a linear transformation to map the emotional state to a fraction of food eaten, 
!with the transformation parameters (k and b) determined by the minimum and maximum emotional states.

!emotion: A real parameter (intent(in)) representing the bird's emotional state, which influences how much food the bird eats.
!min: An optional real parameter (intent(in)) representing the minimum emotional state. If not provided, the default value is HUNGER_MIN.
!max: An optional real parameter (intent(in)) representing the maximum emotional state. If not provided, the default value is HUNGER_MAX.

!fract: A real value representing the fraction of food the bird eats when it is experiencing fear. This value is constrained to be within the 
!   range [0.0, 1.0].

!Parameter Handling: The function checks if the min and max parameters are present. If so, it assigns their values to min_loc and max_loc; 
!    otherwise, it uses the default values HUNGER_MIN and HUNGER_MAX.

!Transformation Parameters Calculation: The function calculates the transformation parameters k and b based on the minimum and maximum 
!   emotional states (min_loc and max_loc), and the predefined constants S1 and S2. S1 and S2 being determined by parameters FRAC_S1 and 
!   FRAC_S2 respectively (see param_info.txt). 


!Fraction Calculation: The function calculates the fraction of food eaten (fract) using the linear transformation defined by k and b, 
!   and ensures that the result is within the range [0.0, 1.0] using the within function.


  real, intent(in) :: emotion
  real, intent(in), optional :: min, max

  real :: fract
  real :: k, b
  real, parameter :: S1 = FRAC_S1, S2 = FRAC_S2

  real :: min_loc, max_loc

  if (present(min)) then
    min_loc = min
  else
    min_loc = HUNGER_MIN
  end if

  if (present(max)) then
    max_loc = max
  else
    max_loc = HUNGER_MAX
  end if

  k = (S2 - S1) / (max_loc - min_loc)
  b = (S2 * min_loc - S1 * max_loc) / ( min_loc-max_loc )

  fract = within(k * emotion + b, 0.0, 1.0)

end function bird_eat_fraction_when_fear




subroutine bird_dead_from_starvation(this)

!The bird_dead_from_starvation subroutine checks if a BIRD class object has died from starvation by calling the is_starved method on 
!the BIRD object. If the bird is starved, the subroutine sets the bird's is_alive attribute to .FALSE., indicating that the bird has died.

  class(BIRD), intent(inout) :: this

  if (this%is_starved()) this%is_alive = .FALSE.

end subroutine bird_dead_from_starvation


subroutine bird_dead_from_background_mortality(this)
 
 class(BIRD), intent(inout) :: this
 integer :: environment_kill_counter
 
  environment_kill_counter = 0
  if (this%environment_kills()) then
    this%is_alive = .FALSE. 
    environment_kill_counter = environment_kill_counter + 1

  end if

end subroutine bird_dead_from_background_mortality



subroutine bird_subtract_metabolic_cost(this)

!The bird_subtract_metabolic_cost subroutine adjusts the weight of a BIRD class object by subtracting a fraction of its current weight,
!which represents the metabolic cost of the bird's activities in relation to it's current weight. This subroutine simulates the energy 
!expenditure of the bird, which can affect its survival and overall health. The cost is applied at every timestep.

  class(BIRD), intent(inout) :: this

  this%weight = this%weight - this%weight * METABOLIC_COST_CONSTANT

end subroutine bird_subtract_metabolic_cost





subroutine bird_do_fly(this, in_environment)
!This subroutine is designed to simulate the behavior of a bird (this) in response to its environment (in_environment). 
!The bird's behavior is influenced by its hunger state, specifically its fear and hunger levels.


!The subroutine first checks if the bird is alive by evaluating this%is_alive. If the bird is not alive, the subroutine
!   ends without performing any actions.

!If the bird is alive, it calculates min_loc and max_loc based on the size of the point array in the in_environment object.

!The subroutine then checks the bird's fear and hunger state by evaluating this%state_fear. 
! If the birds state of fear has a value higher than fear_fly_threshold (range 0-10), the bird will fly and not eat. 
! This is simulated by calling the walk method of the BIRD class with min_loc and max_loc as arguments.

!If the bird's fear and hunger state is not within the specified range, the subroutine does not perform any action, 
!implying that the bird will eat instead of moving. 

  class(BIRD), intent(inout) :: this
  class(whole_environ), intent(in) :: in_environment

  real :: fear_fly_threshold

  integer :: min_loc, max_loc

  fear_fly_threshold = this%genetic_upper_fear() * 0.95

  
  if (this%is_alive) then

    min_loc = 1 
    max_loc = size(in_environment%point) 

    if (this%state_fear > fear_fly_threshold ) then
      call this%walk(min_loc, max_loc)
    end if

  end if

end subroutine bird_do_fly

subroutine bird_do_explore (this, in_environment)
  class(BIRD), intent (inout) :: this
  class(whole_environ), intent(in) :: in_environment

  integer :: min_loc, max_loc

  min_loc = 1
  max_loc = size(in_environment%point)

  if (this%is_alive) then
    if(RAND() < EXPLORATION_PROBABILITY) then
      call this%walk(min_loc, max_loc)
    end if
  end if


end subroutine bird_do_explore



!The two functions below calculate the minimum and maximum genetic limits for a bird's hunger, 
!which is used to determine the bird's hunger state. 
!The function is designed to be called on a BIRD class object and returns a 
!real number representing the lower limit of the bird's genetic hunger range.

!this: A class object of type BIRD. This parameter is passed with the intent(in) attribute, 
!indicating that the function can only read from this object but cannot modify it.

!lower_limit: A real number representing the minimum genetic limit for the bird's hunger.

!range: A real number representing the difference between the maximum and minimum hunger values (HUNGER_MAX and HUNGER_MIN).

!The function calculates the range by subtracting HUNGER_MIN from HUNGER_MAX.
!It then calculates the lower_limit by subtracting the product of range and HUNGER_PLASTICITY from the 
!bird's genetic fear of hunger (this%gene_fear).
!The lower_limit is then adjusted to ensure it falls within the defined range of hunger values 
!(HUNGER_MIN and HUNGER_MAX) using a hypothetical within function. This function is not defined 
!within the provided code snippet, so its specific implementation and behavior are not detailed here.
!Finally, the lower_limit is returned as the result of the function.

function bird_hunger_min_genetic_limit(this) result(lower_limit)
  class(BIRD), intent(in) :: this
  real :: lower_limit

  real :: range

  range = (HUNGER_MAX - HUNGER_MIN) 

  lower_limit = this%gene_hunger - range * HUNGER_PLASTICITY

  lower_limit = within(lower_limit, real(HUNGER_MIN), real(HUNGER_MAX))


end function bird_hunger_min_genetic_limit


function bird_hunger_max_genetic_limit(this) result(upper_limit)
  class(BIRD), intent(in) :: this
  real :: upper_limit

  real :: range

  range = (HUNGER_MAX - HUNGER_MIN) 

  upper_limit = this%gene_hunger + range * HUNGER_PLASTICITY

  upper_limit = within(upper_limit, real(HUNGER_MIN), real(HUNGER_MAX))

end function bird_hunger_max_genetic_limit


!! FEAR
function bird_fear_min_genetic_limit(this) result(lower_limit)
  class(BIRD), intent(in) :: this
  real :: lower_limit

  real :: range

  range = (FEAR_MAX - FEAR_MIN) 

  lower_limit = this%gene_fear - range * FEAR_PLASTICITY

  lower_limit = within(lower_limit, real(FEAR_MIN), real(FEAR_MAX))


end function bird_fear_min_genetic_limit



function bird_fear_max_genetic_limit(this) result(upper_limit)
  class(BIRD), intent(in) :: this
  real :: upper_limit

  real :: range

  range = (FEAR_MAX - FEAR_MIN) 

  upper_limit = this%gene_fear + range * FEAR_PLASTICITY

  upper_limit = within(upper_limit, real(FEAR_MIN), real(FEAR_MAX))


end function bird_fear_max_genetic_limit






subroutine bird_force_fear_within_genetic_limits(this)
! setting a limited variance for genetic dispositions in birds.
! This way there will be heterogeneity in the genetic compositions.
 class(BIRD), intent(inout) :: this

  this%state_fear = within(this%state_fear,                         &
                                  this%genetic_lower_fear(),                    &
                                  this%genetic_upper_fear() )

end subroutine bird_force_fear_within_genetic_limits



subroutine bird_force_hunger_within_genetic_limits(this)
! setting a limited variance for genetic dispositions in birds.
! This way there will be heterogeneity in the genetic compositions.
 class(BIRD), intent(inout) :: this

  this%state_hunger = within(this%state_hunger,                         &
                                  this%genetic_lower_hunger(),                    &
                                  this%genetic_upper_hunger() )

end subroutine bird_force_hunger_within_genetic_limits


!Subroutine for history array

subroutine bird_weight_add_to_history(this, time_step)
  class(BIRD), intent(inout) :: this
  integer, intent(in) :: time_step

  associate ( A => this%weight_history, W => this%weight, i => time_step)
   A(i) = W
  end associate

end subroutine bird_weight_add_to_history





subroutine bird_emotion_state_add_to_history(this, time_step)
  class(BIRD), intent(inout) :: this
  integer, intent(in) :: time_step

  this%fear_history(time_step) = this%state_fear
  this%hunger_history(time_step) = this%state_hunger

end subroutine bird_emotion_state_add_to_history


subroutine population_time_steps(this, environment_in, predator_in) !bird is in this environment = environment_in
  class(POPULATION), intent(inout) :: this
  class(whole_environ), intent(in) :: environment_in
  class(PREDATOR), intent(inout) :: predator_in

  integer :: current_step
  integer :: bird_current !this current bird
  integer :: environment_kill_counter


  logical :: debug_check , debug_ckeck2
   
  environment_kill_counter = 0


  do bird_current = 1, size(this%birds)

    do current_step=1, TOTAL_TIME_STEP

      if (this%birds(bird_current)%is_alive) then


        DEBUG_GET_INFO: if (IS_DEBUG_SCREEN) then
                        print*, " "
                        print *, "Bird", bird_current, "at step", current_step
                        if (this%birds(bird_current)%is_starved()) then
                          print*, "STARVED"
                        end if
                        print *, "Location:", this%birds(bird_current)%x, ",   food",                   &
                                  environment_in%point( this%birds(bird_current)%x )%food_availability

                        print *, "Location:", this%birds(bird_current)%x, ",   risk",           &
                                  environment_in%point( this%birds(bird_current)%x )%predator_frequency


                        print *, "state fear   ", this%birds(bird_current)%state_fear, &
                                          "[", this%birds(bird_current)%genetic_lower_fear(),     &
                                              this%birds(bird_current)%genetic_upper_fear(), "]"


                        print *, "state hunger ", this%birds(bird_current)%state_hunger, &
                                          "[", this%birds(bird_current)%genetic_lower_hunger(),     &
                                              this%birds(bird_current)%genetic_upper_hunger(), "]"

                        print *, " mass", this%birds(bird_current)%weight
                        
                      
                        if (this%birds(bird_current)%environment_kills()) then
                          print*, "*************************************************************************************"
                          print*, "                  BIRD IS KILLED BY BACKGROUND MORTALITY                             "
                          print*, "*************************************************************************************"
                          environment_kill_counter = environment_kill_counter + 1
                        end if
                        print*, "TOTAL BACKGROUND MORTALITY KILL: ", environment_kill_counter

                        

                    end if DEBUG_GET_INFO

        call this%birds(bird_current)%killed_from_starvation() !is the bird dead or alive from starvation?
        !what birds do if is alive
       
        call this%birds(bird_current)%fly(environment_in) !bird is flying

        call predator_in%predator_attack_bird(this%birds(bird_current), environment_in,                   &
              predator_is_present=debug_ckeck2, prey_is_killed=debug_check)! predator attacks bird

        DEBUG_PRINT_KILL: if (IS_DEBUG_SCREEN) then
                                  print *, "Predator ", debug_ckeck2
                                  if (debug_check) print *, "KILLED by predator"
                                end if DEBUG_PRINT_KILL
        call this%birds(bird_current)%do_feed(environment_in)!birds is feeding, and weight is incrementing
        call this%birds(bird_current)%pay_for_life() !bird looses some weight per time step due to metabolism.
        call this%birds(bird_current)%do_explore(environment_in) !bird may explore by moving to another habitat. 
        call this%birds(bird_current)%is_killed_by_background()!background mortality
        


        DEBUG_GET_HISTORY: if (IS_DEBUG_DATA) then
                                    call this%birds(bird_current)%add_to_history(current_step)
                                    call this%birds(bird_current)%emo_state_to_history(current_step) !both fear and hunger within emotion
                            end if DEBUG_GET_HISTORY

        else
                         if (IS_DEBUG_SCREEN) print *, bird_current , "DEAD"
      end if


    end do
  end do

end subroutine population_time_steps

subroutine sort_by_fitness(this)!sorting by bodymass
  class(POPULATION), intent(inout) :: this

  !-----------------------------------------------------------------------------
  ! Sort parents by fitness:

  ! This subroutine sorts the birds in a POPULATION object by their fitness, 
  ! which is represented by their body mass. The sorting is performed in descending
  ! order, so that the fittest birds (those with the highest body mass) are placed 
  ! first in the array.

   ! Initial Sorting: The subroutine begins by calling QsortC on the birds array of the POPULATION object. 
   ! This initiates the QuickSort algorithm, which sorts the birds based on their weight.

   ! Reversing the Order: After the initial sorting, the subroutine reverses the order of the birds array. 
   ! This is done to ensure that the fittest birds (those with the highest body mass) are placed first in the array.

   ! QuickSort Algorithm: The QsortC subroutine is a recursive implementation of the QuickSort algorithm, 
   ! adapted for sorting BIRD objects based on their weight. It uses a helper subroutine, Partition, to partition 
   ! the array around a pivot point, which is the weight of the first bird in the array.

   ! Partition Subroutine: The Partition subroutine partitions the array of BIRD objects around a pivot point, 
   ! which is the weight of the first bird in the array. It uses two indices, i and j, to track the boundaries 
   ! of the partitions. The subroutine iterates through the array, swapping elements as necessary to ensure that 
   ! all elements less than the pivot are to its left, and all elements greater than the pivot are to its right.
  !-----------------------------------------------------------------------------

  call QsortC(this%birds)

  !Reverse the sorting so that the fittest birds(highest mass) go first (10 -> 1).
  this%birds = this%birds(size(this%birds):1:-1)

  contains
    ! The two subroutines below are a variant of the recursive QuickSort
    ! algorithm adapted for BIRD fitness (real)

    recursive subroutine QsortC(A)

      type(BIRD), intent(in out), dimension(:) :: A
      integer :: iq

      if(size(A) > 1) then
        call Partition(A, iq)
        call QsortC(A(:iq-1))
        call QsortC(A(iq:))
      endif

    end subroutine QsortC

    subroutine Partition(A, marker)

      type(BIRD), intent(in out), dimension(:) :: A
      integer, intent(out) :: marker
      integer :: i, j
      type(BIRD) :: temp
      real :: x      ! pivot point

      x = A(1)%weight
      i= 0
      j= size(A) + 1

      do
        j = j-1
        do
            if (A(j)%weight <= x) exit
            j = j-1
        end do
        i = i+1
        do
            if (A(i)%weight >= x) exit
            i = i+1
        end do
        if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
        elseif (i == j) then
            marker = i+1
            return
        else
            marker = i
            return
        endif
      end do

    end subroutine Partition

end subroutine sort_by_fitness

end module organism


!Module for managing simulations
module GA !genetic algorythm
use params
use environment
use BASE_RANDOM
use CSV_IO
use organism
use, intrinsic :: ISO_FORTRAN_ENV, only : OUTPUT_UNIT, ERROR_UNIT
implicit none


type(whole_environ) :: habitat

type(POPULATION),target :: generation_one
type(POPULATION),target :: generation_two

type(POPULATION), pointer :: parent_generation
type(POPULATION), pointer :: offspring_generation

type(PREDATOR) :: predator_in_habitat


!---------------------------------------------
! Making the output data:
!---------------------------------------------
! collumns
! 1: generation nr (int)
! 2: Num alive (int)
! 3: average mass all (real)
! 4: average mass alive (real)
! 5: best mass(real)
! 6: average mass top 25% (real)
! 7: Average gene alive (real)
! 8: Average emotion alive (real)
! 8: standard deviation of mass for alive (real))
! 9: Frequency of predator (real)



integer, parameter :: CSV_COL_NUMBER = 25
character(len=*), dimension(CSV_COL_NUMBER), parameter :: CSV_COLUMNS =    &
  ["GENERATION                   ",                 &  !1
   "NUM-ALIVE                    ",                 &  !2
   "AVERAGE-MASS-ALL             ",                 &  !3
   "AVERAGE-MASS-ALIVE           ",                 &  !4
   "BEST-MASS                    ",                 &  !5
   "AVERAGE-MASS-DOMINANTS       ",                 &  !6
   "MEDIAN MASS ALIVE            ",                 &  !7
   "AVERAGE-GENE-FEAR-ALIVE      ",                 &  !8
   "AVERAGE-GENE-HUNGER-ALIVE    ",                 &  !9
   "AVERAGE-FEAR-STATE-ALIVE     ",                 &  !10
   "AVERAGE-HUNGER-STATE-ALIVE   ",                 &  !11
   "STD-DEV-MASS-ALIVE           ",                 &  !12
   "BIRD-MEETS-PREDATOR-AVERAGE  ",                 &  !13
   "NUM-KILLED-BY-PREDATOR       ",                 &  !14
   "NUM-KILLED-BY-PREDATOR-AVG   ",                 &  !15
   "MUTATION PROBABILITY         ",                 &  !16
   "AVERAGE-GENE-FEAR-DEAD       ",                 &  !17
   "AVERAGE-GENE-HUNGER-DEAD     ",                 &  !18
   "AVERAGE-STATE-FEAR-DEAD      ",                 &  !19
   "AVERAGE-STATE-HUNGER-DEAD    ",                 &  !20
   "STD-DEV-GENE_FEAR-ALIVE      ",                 &  !21
   "STD-DEV-HUNGER-ALIVE         ",                 &  !22
   "AVEREAGE-GENE-FEAR-ALL       ",                 &  !23
   "AVERAGE-PROB-REPRODUCTION    ",                 &  !24
   "Regain Pop when multiply with"]



real, dimension(GENERATIONS, CSV_COL_NUMBER) :: csv_generation_output

private :: CSV_COL_NUMBER, CSV_COLUMNS, csv_generation_output




contains

subroutine genetic_algorithm()
!This subroutine simulates the evolution of a population of birds over multiple generations using a genetic algorithm. 
!The algorithm involves initializing the population, running through time steps, sorting by fitness, selecting the best 
!individuals for reproduction, and swapping generations.

!Initialization: The subroutine begins by initializing the random seed, the habitat, and the predator presence within 
!  the habitat. It also initializes two generations of birds, parent_generation and offspring_generation, and creates 
!  directories for saving weight history and emotion state data.


!Main Loop: The subroutine enters a loop that runs for a specified number of generations. In each generation, it performs the following steps:
  !Initialization: If it's not the first generation, the parent generation is re-initialized.
  !Time Steps: The birds in the parent generation go through time steps, interacting with the habitat and predators.
  !Predator Initialization: The predator presence is updated for the current generation.
  !Sorting by Fitness: The birds are sorted by their fitness, which is determined by their body mass.
  !Data Saving: If specified, the current generation's data is saved to a CSV file.
  !Output Data Building: Characteristics of the parent generation are recorded for output.
  !Check for Surviving Birds: The subroutine checks if there are enough surviving birds. If not, it writes a warning to the error unit and stops the simulation.
  !Selection and Reproduction: The best individuals are selected for reproduction, and the next generation is produced.
  !Generation Swap: The roles of the parent and offspring generations are swapped for the next generation.

!Finalization: After all generations have been processed, the subroutine writes the final generation's data to a CSV file 
!   optionally saves the habitat's nutritional value data if debugging is enabled.


  use params

  use BASE_UTILS, only : TOSTR

  real, dimension(POP_SIZE) :: fitness

  integer :: n_alive_get, n_of_generations

  csv_generation_output = MISSING

  call RANDOM_SEED_INIT()

  call habitat%init()

  Full_Pop_Factor = MISSING
  Last_Full_Pop_Factor = MISSING

  !!!!!call predator_in_habitat%init(Current_Generation)
  
  Is_Evolutionary_Generations = .TRUE. 

  call get_global_runtime_params()

  if (Is_Evolutionary_Generations) then
    n_of_generations = GENERATIONS
  else
    n_of_generations = SIGNAL_ONLY_GENERATIONS
  end if

  parent_generation => generation_one
  offspring_generation => generation_two


  call parent_generation%init_pop(is_gene_init = .TRUE.)
  call offspring_generation%init_pop(is_gene_init = .TRUE.)


!First generation of experimental generations
  if (.NOT. Is_Evolutionary_Generations) then
    call parent_generation%load(Genome_File_Name)
    !call offspring_generation%load(Genome_File_Name
  end if


  call FS_MKDIR(WEIGHT_HISTORY_SUBFOLDER)
  call FS_MKDIR(EMOTION_STATE_SUBFOLDER)


  do Current_Generation = 1, n_of_generations
    
   ! if (.NOT. Is_Evolutionary_Generations) call parent_generation%load(Genome_File_Name)

    write(OUTPUT_UNIT, *) "Generation: ", Current_Generation

    call predator_in_habitat%init(Current_Generation)


    if (Current_Generation > 1 .AND. Is_Evolutionary_Generations)      &
       call parent_generation%init_pop(is_gene_init = .TRUE.)


    call parent_generation%time_steps(habitat, predator_in_habitat)
      
    !sorting by fitness
    call parent_generation%sort()

    SAVE: block

    character(len=500) :: generation_data_file
    character(*), parameter :: GENERATION_SUBFOLDER = "generations/"
    call FS_MKDIR(GENERATION_SUBFOLDER)
    generation_data_file = GENERATION_SUBFOLDER // "gen_" // TOSTR(Current_Generation) // ".csv"
    if (IS_SAVE_POPULATION) call parent_generation%save(generation_data_file)

    end block SAVE


    !record parent characteristics to output file
    call build_output_data(Current_Generation)

    n_alive_get = parent_generation%alive_count()

    write(OUTPUT_UNIT, *)                                                                   &
        "|Num alive: ", n_alive_get,                                                        &
        "|  avg. mass alive:         ", csv_generation_output(Current_Generation, 4),       &
        "|  avg. gene fear alive:    ", csv_generation_output(Current_Generation, 8),       &
        "|  avg. gene hunger alive:  ", csv_generation_output(Current_Generation, 9),      &                                                     
        "|  avg. fear-state alive:   ", csv_generation_output(Current_Generation, 10),       &
        "|  avg. hunger-state alive: ", csv_generation_output(Current_Generation, 11)!,    &
    !    "     Mutation probability: ", MUTATION_PROBABILITY,                                &
    !    "     attack rate: ", predator_in_habitat%risk
 

      call save_history_array(parent_generation, Current_Generation)


    if (n_alive_get < 10) then
      write(ERROR_UNIT, *) "WARNING: Too few birds remain in the generation"

      call CSV_MATRIX_WRITE(csv_generation_output,                      &
                            MODEL_OUTPUT_FILE,colnames=CSV_COLUMNS)
      stop
    end if
    


    call select_and_reproduce()

    if (Is_Evolutionary_Generations) then
      call generation_swap()
    else
      ! If experimental generations 
      parent_generation = reset_generation()
      call minimum_weight_exp_gens
    end if

  

    

  end do

call CSV_MATRIX_WRITE(csv_generation_output,                      &
                      MODEL_OUTPUT_FILE,colnames=CSV_COLUMNS)


! calling function to save nutritional value of each habitat upon environment initiation.
                      ! currently restarting with both predation and nutritional cloning, 
                      ! effectively canceling the signal-increase in experimental generations.
                      ! Therefore habitat cloning is not used during experimental generations per august 23rd 2024.
if (IS_DEBUG_DATA) then
  call habitat%save()
end if


print*, "SIMULATION COMPLETE"

end subroutine genetic_algorithm



! Renamed from save_history_array_weight, since the former is less general.
! We did this since we wanted a general subroutine for saving "history" of
! birds walk and all its updated parameters for each step.
subroutine save_history_array(population_birds, gen_n)
!This subroutine is designet to save the weight and emotion numerical data depicting the birds mass and emotional state
!at each time step. It takes two arguments: 
!  population_birds, which is of the POPULATION class, and gen_n which is an integer value representing the generation. 

! 1. The subroutine begins by using the BASE_UTILS module for string conversion (TOSTR) and declaring several variables 
! and arrays. These include: 
!   var_name: An array of character strings to hold the names of the time steps
!   file_name_weight_history and file_name_emo_history: Allocatable character strings for the file names where the weight 
!      and emotion history will be saved. 
!   weigh_history_population and emotion_history_population: 2D arrays to store the weight and emotion history of each bird 
!      in the population. 
!
! 2. Debug Data Check: The subroutine check if IS_DEBUG_DATA is true. If it is, the subroutine proceeds to generate file 
! names for the weight and emotion history files, appending the generation number to the file names. 
!
! 3. Generating Variable Arrays: For each time step, a variable name is generated in the format STEP_ followed by the 
! time step number. This is done using the TOSTR function, which is presumably defined in the BASE_UTILS module.
!
! 4. Population History Arrays: The subroutine then iterates over each bird in the population and each time step, populating
! the weight_history_population and emotion_history_population arrays with the corresponding data from the population_birds 
! object. 
!
! 5. Writing to CSV Files: the CSV_MATRIX_WRITE subroutine (not shown in this subroutine) is called to write the weight and 
! emotion history arrays to CSV files. The file name and column names (generated variable names) are passed as arguments. 
!
! 6. File Compression: If IS_ZIP_DATA is true, the subroutine compresses the generated CSV-files using a command line tool 
! specified by ZIP_PROG. WAIT = .FALSE. simply means that the subroutine is allowed to continue regardless  of the compression
! being done. 

  use BASE_UTILS, only : TOSTR
  use params

  class(POPULATION), intent(in) :: population_birds
  integer, intent(in) :: gen_n

  !                 STEP_100000
  !                 12345678901
  character(len=11), dimension(TOTAL_TIME_STEP) :: var_name
  character(len=:), allocatable :: file_name_weight_history
  character(len=:), allocatable :: file_name_fear_history
  character(len=:), allocatable :: file_name_hunger_history

  real, dimension(POP_SIZE,TOTAL_TIME_STEP) :: weight_history_population
  real, dimension(POP_SIZE, TOTAL_TIME_STEP) :: fear_history_population, hunger_history_population
  
  integer :: i, j
 if (IS_DEBUG_DATA) then
     file_name_weight_history = WEIGHT_HISTORY_SUBFOLDER// "weight_history_gen_" // TOSTR(gen_n) // ".csv"
     file_name_fear_history = EMOTION_STATE_SUBFOLDER // "fear_history_gen_" // TOSTR(gen_n) // ".csv"
     file_name_hunger_history = EMOTION_STATE_SUBFOLDER // "hunger_history_gen_" // TOSTR(gen_n) // ".csv"


    do i = 1, TOTAL_TIME_STEP
       var_name(i) = "STEP_" // TOSTR(i, TOTAL_TIME_STEP) ! add second parameter that defines the number format
     end do

     do i = 1, POP_SIZE
       do j = 1, TOTAL_TIME_STEP
         weight_history_population(i,j) = population_birds%birds(i)%weight_history(j)
         fear_history_population(i,j) = population_birds%birds(i)%fear_history(j)
         hunger_history_population(i,j) = population_birds%birds(i)%hunger_history(j)

       end do
     end do

     call CSV_MATRIX_WRITE(matrix = weight_history_population,                     &
                           csv_file_name = file_name_weight_history,               &
                           colnames = var_name)

     if (IS_ZIP_DATA) call execute_command_line(ZIP_PROG//' '//file_name_weight_history, WAIT= .FALSE.)! compresses everything if IS_ZIP_DATA = .TRUE.

     call CSV_MATRIX_WRITE(matrix = fear_history_population,                        &
                           csv_file_name = file_name_fear_history,                  &
                           colnames = var_name)

     call CSV_MATRIX_WRITE(matrix = hunger_history_population,                        &
                           csv_file_name = file_name_hunger_history,                  &
                           colnames = var_name)


     if (IS_ZIP_DATA) call execute_command_line(ZIP_PROG//' '//file_name_fear_history, WAIT= .FALSE.)! compresses everything if IS_ZIP_DATA = .TRUE.
     if (IS_ZIP_DATA) call execute_command_line(ZIP_PROG//' '//file_name_hunger_history, WAIT= .FALSE.)! compresses everything if IS_ZIP_DATA = .TRUE.

  end if

end subroutine save_history_array

! collumns
! 1: generation nr (int)
! 2: Num alive (int)
! 3: average mass all (real)
! 4: average mass alive (real)
! 5: best mass(real)
! 6: average mass top 25% (real)
subroutine build_output_data(row)
 ! This subroutine compiles and organizes various statistics and metrics from a simulation
 ! of a population of birds for a given generation.

  character(len = 15) :: value_str  ! Temporary string for formatting real values
  real :: value_real                ! Temporary real variable for storing formatted values

! these values represent the values for each single generation.
  integer, intent(in) :: row        ! 1  
  integer :: n_alive                ! 2  
  real :: average_mass_all          ! 3  
  real :: average_mass_alive        ! 4  
  real :: best_mass                 ! 5  
  real :: average_mass_dominants    ! 6  
  real :: average_gene_fear_alive   ! 7     
  real :: average_gene_hunger_alive ! 8
  real :: average_hunger_alive      ! 9
  real :: average_fear_alive        ! 10
  real :: std_dev_mass_alive        ! 11 
!  real :: freq_predation            ! 12 
  real :: meet_predator             ! 12
  integer :: n_is_killed_by_pred    ! 13 
  real :: n_is_killed_average       ! 14 
  real :: MUTATION_PROBABILITY      ! 15 
  real :: average_gene_fear_dead    ! 16
  real :: average_gene_hunger_dead  ! 17    
  real :: average_fear_dead         ! 18
  real :: average_hunger_dead       ! 19
  real :: std_dev_gene_fear_alive   ! 20
  real :: std_dev_gene_hunger_alive ! 21
  real :: average_gene_fear_all     ! 22
  real :: median_mass_alive         ! 23
  real :: avg_repr_prob             ! 24
  real :: Full_Pop_Factor           ! 25





  ! temporary array keeping mass of birds that are alive, mass of dead birds
  ! is equal to MISSING
  real, dimension(POP_SIZE) :: mass_alive

  ! temporary array to keep population genes
  integer, dimension(POP_SIZE) :: gene_fear_alive, gene_hunger_alive, gene_fear_all

  !temporary array to keep the fear state (state_fear) of alive birds.
  real, dimension(POP_SIZE) :: fear_alive

  !temporary array to keep the hunger state (state_hunger) of alive birds.
  real, dimension(POP_SIZE) :: hunger_alive


  !for figuring out average gene of dead birds.
  real, dimension(POP_SIZE) :: gene_fear_dead, gene_hunger_dead

  !for figuring out average emotion of dead birds.
  real, dimension(POP_SIZE) :: fear_state_dead, hunger_state_dead

  !For figuring out how many birds have met predator

  integer, dimension(POP_SIZE) :: meet_predator_count

 ! Initialize temporary arrays with default values

 real, dimension(POP_SIZE) :: probability_repro

  mass_alive = MISSING

  gene_fear_dead = MISSING
  gene_hunger_dead = MISSING
  

  gene_fear_alive = -9999
  gene_hunger_alive = -9999
  

  fear_alive = -9999
  fear_state_dead = -9999

  hunger_alive = -9999

  meet_predator_count = UNDEFINED




! Calculate various stats and metrics. 
!NUMBER OF LIVING BIRDS AND THEIR AVERAGE MASS
  n_alive = count(parent_generation%birds%is_alive)
  average_mass_all = average(parent_generation%birds%weight)


  where(parent_generation%birds%is_alive)
    mass_alive = parent_generation%birds%weight
  end where


! AVERAGE MASS ALIVE
  average_mass_alive = average(mass_alive)
! BEST MASS
  best_mass = parent_generation%birds(1)%weight
!AVERAGE MASS DOMINANTS
  average_mass_dominants = average(parent_generation%birds(1:GA_REPRODUCE_N)%weight)

!MEDIAN OF MASS ALIVE BIRDS: 
! Calculate median of mass of alive birds

  median_mass_alive = median(mass_alive)

!AVG REPRODUCTION PROBABILITY (FITNESS)
  probability_repro = prob_repr(average_mass_alive, parent_generation%birds(n_alive)%weight, mass_alive(1))

  avg_repr_prob = average(probability_repro)



!AVERAGE GENE FEAR/HUNGER OF ALIVE BIRDS
  where(parent_generation%birds%is_alive)
    gene_fear_alive = parent_generation%birds%gene_fear
    gene_hunger_alive = parent_generation%birds%gene_hunger
  end where

  average_gene_fear_alive = average( real(gene_fear_alive) )
  average_gene_hunger_alive = average( real(gene_hunger_alive) )




!AVERAGE GENE FEAR OF DEAD BIRDS

  where(.not.parent_generation%birds%is_alive)
    gene_fear_dead = parent_generation%birds%gene_fear
    gene_hunger_dead = parent_generation%birds%gene_hunger
  end where

  average_gene_fear_dead = average( real(gene_fear_dead) )
  average_gene_hunger_dead = average( real(gene_hunger_dead) )

! AVERAGE GENE FEAR ALL
  gene_fear_all = real(parent_generation%birds%gene_fear)
  average_gene_fear_all = average( real(gene_fear_all) )


!AVERAGE EMOTION OF DEAD BIRDS

  where(.not.parent_generation%birds%is_alive)
    fear_state_dead = parent_generation%birds%state_fear
    hunger_state_dead = parent_generation%birds%state_hunger
    
  end where
  average_fear_dead = average( real( fear_state_dead ) )
  average_hunger_dead = average( real( hunger_state_dead ) )

!AVERAGE STATE FEAR /HUNGER OF ALIVE BIRDS
  where(parent_generation%birds%is_alive)
      fear_alive = parent_generation%birds%state_fear
      hunger_alive = parent_generation%birds%state_hunger
  end where

  average_fear_alive = average( real(fear_alive) )
  average_hunger_alive = average( real(hunger_alive) )


! STANDARD DEVIATION OF MASS
  std_dev_mass_alive = std_dev(mass_alive)


! STANDARD DEVIATION OF GENE

  std_dev_gene_fear_alive = std_dev(real(gene_fear_alive))
  std_dev_gene_hunger_alive = std_dev(real(gene_hunger_alive))

! MEET THE PREDATOR AVERAGE
! Calculate the average number of alive birds' number of predator encounters
   where(parent_generation%birds%is_alive)
     meet_predator_count = ( parent_generation%birds%bird_meets_predator_counter )

   elsewhere
     meet_predator_count = UNDEFINED !for all birds that are dead: assign UNDEFINED value

   end where

   meet_predator = average( real(meet_predator_count) ) 


! Update and display global variable
   Full_Pop_Factor = update_pop_regain_size_factor(n_alive, POP_SIZE)

  ! Update Last_Full_Pop_Factor if this is the last generation
  if (.not. Is_Evolutionary_Generations .and. row == SIGNAL_ONLY_GENERATIONS - (SIGNAL_ONLY_GENERATIONS + 1)) then
    Last_Full_Pop_Factor = update_pop_regain_size_factor(n_alive, size(parent_generation%birds))
    Full_Pop_Factor = Last_Full_Pop_Factor
  end if

    
!PREDATOR KILLS BIRD

 n_is_killed_by_pred = count(parent_generation%birds%is_killed_by_predator) !count works with arrays, and therefore we do not call birds(1)
 n_is_killed_average = (n_is_killed_by_pred / meet_predator) 

!
write(value_str, '(f6.3)') MUTATION_PROBABILITY
write(value_str, '(f6.3)') n_is_killed_average


 ! Add calculated values to the output array

  csv_generation_output(row, 1)  = row                       
  csv_generation_output(row, 2)  = n_alive                   
  csv_generation_output(row, 3)  = average_mass_all          
  csv_generation_output(row, 4)  = average_mass_alive        
  csv_generation_output(row, 5)  = best_mass                 
  csv_generation_output(row, 6)  = average_mass_dominants    
  csv_generation_output(row, 7)  = median_mass_alive
  csv_generation_output(row, 8)  = average_gene_fear_alive   
  csv_generation_output(row, 9)  = average_gene_hunger_alive 
  csv_generation_output(row, 10) = average_fear_alive        
  csv_generation_output(row, 11) = average_hunger_alive     
  csv_generation_output(row, 12) = std_dev_mass_alive               
  csv_generation_output(row, 13) = meet_predator            
  csv_generation_output(row, 14) = n_is_killed_by_pred      
  csv_generation_output(row, 15) = n_is_killed_average      
  csv_generation_output(row, 16) = MUTATION_PROBABILITY     
  csv_generation_output(row, 17) = average_gene_fear_dead   
  csv_generation_output(row, 18) = average_gene_hunger_dead 
  csv_generation_output(row, 19) = average_fear_dead        
  csv_generation_output(row, 20) = average_hunger_dead      
  csv_generation_output(row, 21) = std_dev_gene_fear_alive  
  csv_generation_output(row, 22) = std_dev_gene_hunger_alive
  csv_generation_output(row, 23) = average_gene_fear_all
  csv_generation_output(row, 24) = avg_repr_prob
  csv_generation_output(row, 25) = Full_Pop_Factor 


end subroutine build_output_data



subroutine select_and_reproduce()
  implicit none

  integer :: i, j, ii, alive_counter, new_popsize_fact, random_from_alive
  integer :: n_offspring_reproducing_parents
  type(BIRD), allocatable :: temp_array(:)
  
  ! Initialize new offspring
  j = 1
  n_offspring_reproducing_parents = 1 

  ! Calculate the number of alive birds
  alive_counter = parent_generation%alive_count()
  
  ! Update Full_Pop_Factor if in evolutionary generations
  if (Is_Evolutionary_Generations) then
    Full_Pop_Factor = update_pop_regain_size_factor(alive_counter, POP_SIZE)
  else
    Full_Pop_Factor = Last_Full_Pop_Factor
  end if

  ! Allocate temp_array if it's not already allocated and we are not in evolutionary generations
  if (.not. allocated(temp_array) .and. .not. Is_Evolutionary_Generations) then
    allocate(temp_array(alive_counter))
  end if

  ! Selection and reproduction process
  ! Loop through each bird in the parent generation
  do i = 1, parent_generation%alive_count()
    associate(m => parent_generation%birds(i)%weight, &
              m_0 => parent_generation%smallest_weight(), &
              m_max => parent_generation%biggest_weight())

      ! If in evolutionary generations, reproduce based on the probability of reproduction 
      if (Is_Evolutionary_Generations) then
        if (RAND() < prob_repr(m, m_0, m_max)) then
          offspring_generation%birds(j) = parent_generation%birds(i)
          j = j + 1
        end if
      else
        ! Otherwise, just copy alive birds to the offspring generation directly
        if (parent_generation%birds(i)%is_alive .and. &
            prob_repr(m, m_0, m_max) > GA_PROB_REPR_MIN) then
          offspring_generation%birds(i) = parent_generation%birds(i)
        end if
      end if
    end associate
  end do
  
  n_offspring_reproducing_parents = j - 1 



  ! After all alive birds have been allocated to the offspring array,
  ! fill the remaining slots in the offspring array with randomly chosen alive birds 
  ! based on the size determined by survivors * Full_Pop_Factor.  
  if (.not. Is_Evolutionary_Generations) then

    new_popsize_fact = int(alive_counter * Full_Pop_Factor)
    print*, "Full_Pop_Factor is: ", Full_Pop_Factor
    do i = alive_counter + 1, new_popsize_fact
      random_from_alive = RAND_I(1, alive_counter)
      offspring_generation%birds(i) = parent_generation%birds(random_from_alive)
    end do
  end if 

  ! In evolutionary generations, fill the remaining slots in the offspring array
  ! with randomly chosen birds from the parent generation up to POP_SIZE.
  if (Is_Evolutionary_Generations) then
    do i = n_offspring_reproducing_parents + 1, POP_SIZE
      ii = RAND_I(1, POP_SIZE/2)
      offspring_generation%birds(i) = parent_generation%birds(ii)
    end do
  end if

end subroutine select_and_reproduce






subroutine minimum_weight_exp_gens()

  integer :: i
  real :: m, m_0, m_max

  associate(  m =>  parent_generation%birds(i)%weight,          &
              m_0 => parent_generation%smallest_weight(),       &
              m_max => parent_generation%biggest_weight())


  ! need minimum threshold weight that is dynamically set by smallest surviving mother
  ! that produces eggs from generation "Genome_File_Name"
   if (.NOT.((parent_generation%birds(i)%is_alive) .and.   &
      (prob_repr(m, m_0, m_max) > GA_PROB_REPR_MIN ))) then  !make a parameter 

 
   end if

  end associate
end subroutine minimum_weight_exp_gens



subroutine generation_swap()

  if(associated(parent_generation, target=generation_one)) then
    parent_generation => generation_two
    offspring_generation => generation_one
  else
    parent_generation => generation_one
    offspring_generation => generation_two
  end if

end subroutine generation_swap



function reset_generation() result(popout)

 type(POPULATION) :: popout

 !This subroutine resets the parent generation to its initial state
 call popout%init_pop(is_gene_init = .TRUE.)
 call popout%load(Genome_File_Name) 

end function reset_generation

!Weighted average and it's "factor": 
!> Calculates the fitness factor for a population generation.
!!
!! The fitness factor is used to determine the probability of reproduction for each individual in the population.
!! It is calculated as the weighted average of the reproduction probabilities for all alive individuals in the generation.
!! The weight is determined by the difference between the current population size and the number of alive individuals.
!! If the full population size is provided, it is used to calculate the factor. Otherwise, the current population size is used.
!! The factor is ensured to be non-negative in case the current population exceeds the target population size.
function fitness_factor(p_gen, full_pop) result(factor)
  ! Local variables
  class(POPULATION), intent(in) :: p_gen
  integer :: i, alive_count 
  real :: total_probability
  real :: factor
  integer, optional, intent(in) :: full_pop


  !factor = Weighed_Avg_Factor
  
  associate(m => p_gen%birds(i)%weight, &
            m_0 => p_gen%smallest_weight(), &
            m_max => p_gen%biggest_weight())

    ! Check if the factor has already been calculated
      total_probability = 0.0

      ! Calculate the total weight and total probability
      alive_count = count(p_gen%birds%is_alive)
      do i = 1, alive_count
        total_probability = total_probability + prob_repr(m, m_0, m_max)
      end do

      ! Calculate the weighted average
      if (total_probability > 0.0) then
        if (present(full_pop)) then
        factor = real(full_pop - alive_count) / total_probability
         ! Ensure factor is not negative (in case current population exceeds POP_SIZE)
        else ! If full_pop is not provided, use the current population size
          factor = real(POP_SIZE - alive_count) / total_probability
        end if 
        factor = max(factor, 0.0)
      else
        factor = 0.0
      end if



  end associate
end function fitness_factor

end module GA

