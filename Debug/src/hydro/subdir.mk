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
	@echo 'Invoking: Intel C++ Compiler'
	mpicc -g -O0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


