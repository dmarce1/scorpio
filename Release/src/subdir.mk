################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/assert.cpp \
../src/comm.cpp \
../src/legendre.cpp \
../src/main.cpp \
../src/physical_constants.cpp \
../src/program.cpp \
../src/reconstruct.cpp \
../src/tag.cpp 

C_SRCS += \
../src/backtrace_symbols.c 

OBJS += \
./src/assert.o \
./src/backtrace_symbols.o \
./src/comm.o \
./src/legendre.o \
./src/main.o \
./src/physical_constants.o \
./src/program.o \
./src/reconstruct.o \
./src/tag.o 

C_DEPS += \
./src/backtrace_symbols.d 

CPP_DEPS += \
./src/assert.d \
./src/comm.d \
./src/legendre.d \
./src/main.d \
./src/physical_constants.d \
./src/program.d \
./src/reconstruct.d \
./src/tag.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -ipo -inline-level=2 -I"/home/dmarce1/Scorpio/src" -DNDEBUG -openmp-stubs -diag-disable 161 -fp-speculation=fast -fp-model fast=2 -xHost -unroll-aggressive -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C Compiler'
	mpicc -DNDEBUG -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


