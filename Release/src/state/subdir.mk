################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/state/state.cpp 

OBJS += \
./src/state/state.o 

CPP_DEPS += \
./src/state/state.d 


# Each subdirectory must supply rules for building sources it contributes
src/state/%.o: ../src/state/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	mpicc -O2 -DNDEBUG -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


