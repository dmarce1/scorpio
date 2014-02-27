################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/hydro/hydro.cpp \
../src/hydro/hydro_static.cpp 

OBJS += \
./src/hydro/hydro.o \
./src/hydro/hydro_static.o 

CPP_DEPS += \
./src/hydro/hydro.d \
./src/hydro/hydro_static.d 


# Each subdirectory must supply rules for building sources it contributes
src/hydro/%.o: ../src/hydro/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -ipo -inline-level=2 -I"/home/dmarce1/Scorpio/src" -DNDEBUG -openmp-stubs -diag-disable 161 -fp-speculation=fast -fp-model fast=2 -xHost -unroll-aggressive -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


