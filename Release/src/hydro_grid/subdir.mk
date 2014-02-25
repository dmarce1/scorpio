################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/hydro_grid/hydro_grid.cpp \
../src/hydro_grid/hydro_grid_static.cpp 

OBJS += \
./src/hydro_grid/hydro_grid.o \
./src/hydro_grid/hydro_grid_static.o 

CPP_DEPS += \
./src/hydro_grid/hydro_grid.d \
./src/hydro_grid/hydro_grid_static.d 


# Each subdirectory must supply rules for building sources it contributes
src/hydro_grid/%.o: ../src/hydro_grid/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -opt-prefetch=3 -ipo -inline-level=2 -I"/home/dmarce1/Scorpio/src" -DNDEBUG -diag-disable 161 -fp-speculation=fast -fp-model fast=2 -xHost -openmp-stubs -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


