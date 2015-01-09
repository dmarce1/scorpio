################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/FMM/FMM.cpp 

OBJS += \
./src/FMM/FMM.o 

CPP_DEPS += \
./src/FMM/FMM.d 


# Each subdirectory must supply rules for building sources it contributes
src/FMM/%.o: ../src/FMM/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	mpicc -O2 -DNDEBUG -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


