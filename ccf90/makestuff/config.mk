# Compilers
CC := gcc
FC := gfortran
LD := $(FC)

# Flags
F90_MOD_FLAG := -J
CPPFLAGS := -DUSE_POPCNT #-DDISABLE_OPT_T3

FFLAGS := -O3 -g -pedantic -Wall -Wextra -std=f2018 -fdefault-real-8 -fdefault-integer-8

# Linking
LDFLAGS := -llapack -lblas
LIBS := -Wl, --start-group /opt/intel/lib/intel64_mac  
