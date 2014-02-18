################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/oct_node/child_index.cpp \
../src/oct_node/oct_connections.cpp \
../src/oct_node/oct_face.cpp \
../src/oct_node/oct_node.cpp \
../src/oct_node/oct_node_static.cpp 

OBJS += \
./src/oct_node/child_index.o \
./src/oct_node/oct_connections.o \
./src/oct_node/oct_face.o \
./src/oct_node/oct_node.o \
./src/oct_node/oct_node_static.o 

CPP_DEPS += \
./src/oct_node/child_index.d \
./src/oct_node/oct_connections.d \
./src/oct_node/oct_face.d \
./src/oct_node/oct_node.d \
./src/oct_node/oct_node_static.d 


# Each subdirectory must supply rules for building sources it contributes
src/oct_node/%.o: ../src/oct_node/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -opt-prefetch=3 -ipo -inline-level=2 -I"/home/dmarce1/Scorpio/src" -DNDEBUG -fp-speculation=fast -fp-model fast=2 -xHost -openmp-stubs -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


