################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/ba.cpp \
../src/baWithPrior.cpp \
../src/main.cpp 

OBJS += \
./src/ba.o \
./src/baWithPrior.o \
./src/main.o 

CPP_DEPS += \
./src/ba.d \
./src/baWithPrior.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -I/home/l_vis/Programas/Eigen_3.3.5 -I"/home/l_vis/git/toySlam/toyBA/inc" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


