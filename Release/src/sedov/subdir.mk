################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/sedov/sedov.c 

F_SRCS += \
../src/sedov/sedov3.f 

OBJS += \
./src/sedov/sedov.o \
./src/sedov/sedov3.o 

C_DEPS += \
./src/sedov/sedov.d 


# Each subdirectory must supply rules for building sources it contributes
src/sedov/%.o: ../src/sedov/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C Compiler'
	mpicc -DNDEBUG -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/sedov/sedov3.o: ../src/sedov/sedov3.f
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) IA-64 Fortran Compiler'
	ifort -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/sedov/sedov3.o: ../src/sedov/sedov3.f


