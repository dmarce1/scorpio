################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/hydro_FMM_grid/hydro_FMM_grid.cpp 

OBJS += \
./src/hydro_FMM_grid/hydro_FMM_grid.o 

CPP_DEPS += \
./src/hydro_FMM_grid/hydro_FMM_grid.d 


# Each subdirectory must supply rules for building sources it contributes
src/hydro_FMM_grid/%.o: ../src/hydro_FMM_grid/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -opt-prefetch=3 -ipo -inline-level=2 -I"/home/dmarce1/Scorpio/src" -DNDEBUG -fp-speculation=fast -fp-model fast=2 -xHost -openmp-stubs -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


