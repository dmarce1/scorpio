################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/sedov/sedov.c 

OBJS += \
./src/sedov/sedov.o 

C_DEPS += \
./src/sedov/sedov.d 


# Each subdirectory must supply rules for building sources it contributes
src/sedov/%.o: ../src/sedov/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	mpicc -g -O0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


