################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/blast_wave/grid_blast_wave.cpp 

OBJS += \
./src/blast_wave/grid_blast_wave.o 

CPP_DEPS += \
./src/blast_wave/grid_blast_wave.d 


# Each subdirectory must supply rules for building sources it contributes
src/blast_wave/%.o: ../src/blast_wave/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	mpicc -g -O0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


