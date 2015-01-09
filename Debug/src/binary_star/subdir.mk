################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/binary_star/binary_star.cpp \
../src/binary_star/binary_star_read.cpp \
../src/binary_star/binary_star_static.cpp \
../src/binary_star/dwd.cpp \
../src/binary_star/lane_emden.cpp 

OBJS += \
./src/binary_star/binary_star.o \
./src/binary_star/binary_star_read.o \
./src/binary_star/binary_star_static.o \
./src/binary_star/dwd.o \
./src/binary_star/lane_emden.o 

CPP_DEPS += \
./src/binary_star/binary_star.d \
./src/binary_star/binary_star_read.d \
./src/binary_star/binary_star_static.d \
./src/binary_star/dwd.d \
./src/binary_star/lane_emden.d 


# Each subdirectory must supply rules for building sources it contributes
src/binary_star/%.o: ../src/binary_star/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	mpicc -g -O0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


