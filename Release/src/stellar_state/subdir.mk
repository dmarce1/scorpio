################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/stellar_state/state.cpp 

OBJS += \
./src/stellar_state/state.o 

CPP_DEPS += \
./src/stellar_state/state.d 


# Each subdirectory must supply rules for building sources it contributes
src/stellar_state/%.o: ../src/stellar_state/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -opt-prefetch=2 -ipo -inline-level=2 -I"/home/dmarce1/Scorpio/src" -DNDEBUG -fp-model fast=2 -xHost -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


