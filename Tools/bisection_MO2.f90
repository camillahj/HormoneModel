!*******************************************************************************
! SVN Version ID: $Id: bisection_MO2.f90 6982 2018-04-04 07:04:57Z cje012 $
!*******************************************************************************
!!!!
!!!! Program that uses the bisection method to find wanted M_O2_size values per year
!!!! for plotting

!!!! It uses the same parameterfile and the same files as the main program so
!!!! it must be placed in the same folder as these files 
!!!!
!!!! NB! This program needs the resultfile called ValueOfE.bin to run
!!!!

program bisection_MO2
 
  !Module containing all parameter values
  use parameter_values
  
  !Module containing all functions as well as some parameters (like pi)
  ! used for calculations
  use functions

  !Module for printing and reading results to binary files
  use binary_IO

  implicit none
  
  !Variables used for loops
  integer :: i,E
  
  !Parameters
  real(kind=RP), parameter :: length            = 20._RP
  real(kind=RP), parameter :: reserve_fullness  = 0.2_RP 
  real(kind=RP), parameter :: THF               = (THF_max-THF_min)/2._RP
  
  real(kind=RP), parameter :: M_O2_size_min = 0.04_RP, M_O2_size_max = 0.22_RP, M_O2_size_step = 0.02_RP
  real(kind=RP), parameter :: E_real_min = 1._RP, E_real_max = E_categories, E_real_step = 0.1_RP
  
  !Making array with a sequence of M_O2_size per years values 
  real(kind=RP), dimension(*), parameter :: M_O2_size_array = [( M_O2_size_min + M_O2_size_step * (i-1), &
                                            i=1, floor((M_O2_size_max-M_O2_size_min)/M_O2_size_step + 1) )]
  !Making array with a sequence of E_real values
  real(kind=RP), dimension(*), parameter :: E_real_array = [( E_real_min + E_real_step * (i-1), &
                                            i=1, floor((E_real_max-E_real_min)/E_real_step + 1) )]

  !Variables
  real(kind=RP)                 :: M_O2_size_wanted, M_O2_size, M_O2_size_last
  real(kind=RP)                 :: M_size, M_O2, M_O2_size_tstep 
  real(kind=RP)                 :: OXF, OXF_llimit, OXF_ulimit
  real(kind=RP)                 :: SDA
  real(kind=RP)                 :: O2_used, O2max_std, O2max_THF
  real(kind=RP)                 :: SMR_coeff, SMR_exp_CJ, SMR_coeff_CJ, SMR_std, SMR_THF, SMR_somatic_std
  real(kind=RP)                 :: weight_somatic, weight
  real(kind=RP)                 :: E_real, ValueOfE_real
  real(kind=RP)                 :: reserves_max, reserves
  real(kind=RP)                 :: target_intake, intake, foraging_required, foraging_cost
  real(kind=RP)                 :: surplus_before_growth, growth_cost
  real(kind=RP)                 :: conversion_cost_via_reserves, conversion_cost_to_growth 
  real(kind=RP), dimension(1:3) :: phiVE !Array containing phi, ValueOfE and E_real, and is returned from the autocorrelatedE function
  real(kind=RP), dimension(1:34, 1:t_max+1, 1:n_max) :: ind !Individual matrix
  
  !Result array positions is as follows: environment, M_O2, surplus_before_growth
  real(kind=RP), dimension(1:size(E_real_array),1:size(M_O2_size_array)) ::  surplus_results
  
  !SMR
  !Clarke & Johnston 1999 - General teleost fish
  SMR_exp_CJ = 0.80_RP
  SMR_coeff_CJ = exp(-5.43_RP)                                         !mmol O2 g-1 h-1 Originally -5.43
  SMR_coeff_CJ = SMR_coeff_CJ * 434._RP * 24._RP                       !J g-1 d-1  Conversion from paper, p 895 %Units conversion
  SMR_coeff_CJ = SMR_coeff_CJ * (0.001_RP**(-SMR_exp_CJ)) * t_duration !J kg-1 timestep-1 %Units conversion
  SMR_coeff = SMR_coeff_CJ * (SMR_coeff_weight**(SMR_exp_CJ-SMR_exp))  !Correction according to new exponent, converted for a 3-kg fish

  !Calculations before running
  weight_somatic  = k_somatic*(length**3._RP)
  reserves_max = (k_max_reserves*(length**3._RP) - weight_somatic)*reserves_energy_density
  reserves = reserves_max * 0.2_RP
  weight = weight_somatic + reserves/reserves_energy_density  
  SMR_somatic_std = SMR_coeff * (weight_somatic**SMR_exp)
  SMR_std = SMR_coeff * (weight**SMR_exp)
  O2max_std = O2max_coeff * (weight_somatic**O2max_exp) 
  
  SMR_THF = SMR_std*(1._RP+((THF/THF_max)-0.5_RP)*THF_SMR_effect)
  O2max_THF = O2max_std * (1._RP + ((THF / THF_max) - 0.5_RP) * THF_O2max_effect)
  
  M_size = M_size_coeff*(length**M_size_exp) 
  
  
  !Looping for each environment
  do E=1, size(E_real_array)
  
    !Getting environment from array
    E_real = E_real_array(E)
  
    !Setting up the environment to get ValueOfE_categories
    phiVE = autocorrelatedE(E_real, 1._RP, StochScale, ValueOfE_min, ValueOfE_max, E_categories) !Since AutoCorr = 1, lastE=E_real and lastphi=phi
    ValueOfE_real = phiVE(2)  
    
    !Looping for each wanted M_O2 value
    do i = 1, size(M_O2_size_array)
    
      !Getting wanted value from array
      M_O2_size_wanted = M_O2_size_array(i)
      
      !Initialising M_O2 before loop
      M_O2      = -1000
      M_O2_size_last = -1000
      
      !Initialising OXF
      OXF_llimit = OXF_min !OXF lower limit
      OXF_ulimit = OXF_max !OXF upper limit
      OXF = (OXF_llimit + OXF_ulimit)/2._RP
      
      !Loop to find the suplus before growth for a given M_O2_size and environment
      do while (M_O2_size_wanted /= M_O2_size)
        
        !Calculating target_intake, intake, foraging_required and foraging_cost  
        target_intake = (OXF / OXF_max) * OXF_effect_on_intake
        intake = target_intake *  SMR_somatic_std
        foraging_required = target_intake / ValueOfE_real
        foraging_cost = foraging_cost_coeff * foraging_required * SMR_std
        
        !Calculating SDA
        SDA = SDA_coeff * intake
        
        !Calculating surplus before_growth
        surplus_before_growth = intake - SDA - SMR_THF - foraging_cost
        
        !Calculating conversion_cost_via_reserves
        if (surplus_before_growth > 0._RP) then !When surplus before growth is positive
          conversion_cost_via_reserves = surplus_before_growth * (1._RP-conversion_efficiency_reserves) 
        else !When surplus before growth is negative
          conversion_cost_via_reserves =  abs(surplus_before_growth) / conversion_efficiency_reserves * (1._RP-conversion_efficiency_reserves) 
        end if
        
        !Calculating growth cost when assuming no change to reserves
        growth_cost = conversion_efficiency_growth*(surplus_before_growth-conversion_cost_via_reserves)
        
        !Calculating conversion_cost_to_growth
        conversion_cost_to_growth = (growth_cost/conversion_efficiency_growth)*(1._RP-conversion_efficiency_growth)
        
        !Calculating O2_used
        O2_used = SMR_THF + SDA + foraging_cost + conversion_cost_via_reserves + conversion_cost_to_growth
        
        !Calculating M_O2
        M_O2 = M_O2_coeff * (O2_used / O2max_THF)**M_O2_exp
        
        !Calculating M_O2_size
        M_O2_size_tstep = M_O2 * M_size
        
        !Calculating M_O2_size per year
        M_O2_size = M_O2_size_tstep *(365._RP / t_duration) 
        
        !Print to see that the program is running
        print*, E_real

        if (M_O2_size < M_O2_size_wanted) then    !If M_O2_size is lower than M_O2_size_wanted
          OXF_llimit = OXF              !Set the lower limit to OXF
        else                            !If M_O2_size is lower than M_O2_size_wanted
          OXF_ulimit = OXF              !Set OXF as the new upper limit
        end if
        
        !Calculate new OXF 
        OXF = (OXF_llimit + OXF_ulimit)/2._RP 
        
        !If the last calculated M_O2_size is the same as the M_O2_size calculated this round
        ! exit the loop to avoid infinite loop
        if(M_O2_size_last==M_O2_size) exit
        
        !Saving last value to avoid that the loop gets stuck
        M_O2_size_last = M_O2_size

      end do
      
      
      
      !Saving the values that worked in the array
      surplus_results(E,i) = surplus_before_growth



    end do
      
     
    
  end do

  !Saving array to binary files
  call array_to_binary(E_real_array,"bisection_E_real.bin")
  call array_to_binary(M_O2_size_array,"bisection_M_O2_size.bin") 
  call array_to_binary(surplus_results,"bisection_surplus_before_growth.bin")
  
  print*, "Finished \(^o^)/"
  
end program bisection_MO2
