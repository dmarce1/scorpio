################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/assert.cpp \
../src/comm.cpp \
../src/legendre.cpp \
../src/main.cpp \
../src/parameter_reader.cpp \
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
./src/parameter_reader.o \
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
./src/parameter_reader.d \
./src/physical_constants.d \
./src/program.d \
./src/reconstruct.d \
./src/tag.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	mpicc -O2 -DNDEBUG -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	mpicc -O2 -DNDEBUG -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


