# Brief Makefile for modern fortran projects
#
# Inspired (almost the same as) in James Spencer's
# http://github.com/jsspencer/generic_makefile
#
# 2019, J. Emiliano Deustua


# Directory structure and config
# ------------------------------
BUILD_DIR := build
BIN_DIR := bin
DEPEND_DIR := $(BUILD_DIR)/depend
VPATH := src:src/cc:src/utils

# Include architecture configuration file
include config.mk
ifeq ($(findstring DISABLE_OPT_T3, $(CPPFLAGS)),)
		VPATH := $(VPATH):src/cc_opt_ccsdt
endif


# Program naming
# --------------
PROGRAM := ccq
TEST_DIR := tests
ENTRY := main.f90

SHARED_LIBRARY := libccq.so
STATIC_LIBRARY := libccq.a


# Source handling
# ---------------
CEXTS := .c
FEXTS := .f .F .f90 .F90
EXTS := $(FEXTS) $(CEXTS)

# Load source dirs
SRCDIRS := $(subst :, , $(VPATH))

# Find all files
find_files = $(foreach ext, $(EXTS), $(wildcard $(dir)/*$(ext)))
SRCS := $(foreach dir, $(SRCDIRS), $(find_files))

# Add full path to an object
objects_path = $(addprefix $(BUILD_DIR)/, $(addsuffix .o,$(basename $(notdir $(1)))))
OBJS := $(call objects_path, $(SRCS))

# Dependency resolutions
F_FILES = $(filter $(addprefix %,$(FEXTS)), $(SRCS))
ifeq ($(F_FILES),)
F_DEPEND :=
else
F_DEPEND := $(DEPEND_DIR)/fortran.d
endif

# Library objects
ENTRY_OBJ := $(call objects_path, $(ENTRY))
LIB_OBJS := $(filter-out $(ENTRY_OBJ),$(OBJS))
AR := ar
ARFLAGS := -rcs

# Force good coding
WARNINGS := -Wall -std=f2018 -fall-intrinsics


# Enviroment data logging
# -----------------------
COMP_HOST := $(shell hostname --fqdn)
COMP_TIME := $(shell date)
GIT_SHA1 := $(shell git rev-parse HEAD 2> /dev/null || echo "not available.")
# Append -dirty to SHA1 if source code contains changes.
GIT_SHA1 := $(GIT_SHA1)$(shell test -z "$$(git status --porcelain 2>/dev/null)" || echo -dirty)

# Compilation macros
.SUFFIXES:
.SUFFIXES: $(EXTS)
.PHONY: all lib clean cleanall debug test

all: $(BIN_DIR)/$(PROGRAM)

lib: FFLAGS += -fpic
lib: CFLAGS += -fpic
lib: $(BIN_DIR)/$(SHARED_LIBRARY) $(BIN_DIR)/$(STATIC_LIBRARY)


# Fortran
# -------
$(BUILD_DIR)/printing.o: src/printing.F90
	$(FC) -c $(FFLAGS) -DCOMP_HOST="'$(COMP_HOST)'" -DCOMP_TIME="'$(COMP_TIME)'" -DVERSION="'$(GIT_SHA1)'" -DFLAGS="'$(FFLAGS)'" -o $@ $< $(F90_MOD_FLAG) $(BUILD_DIR)

$(BUILD_DIR)/%.o: %.F
	$(FC) $(CPPFLAGS) -c $(FFLAGS) $< -o $@ $(F90_MOD_FLAG) $(BUILD_DIR)

$(BUILD_DIR)/%.o: %.F90
	$(FC) $(CPPFLAGS) -c $(WARNINGS) $(FFLAGS) $< -o $@ $(F90_MOD_FLAG) $(BUILD_DIR)

$(BUILD_DIR)/%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@ $(F90_MOD_FLAG) $(BUILD_DIR)

$(BUILD_DIR)/%.o: %.f90
	$(FC) -c $(WARNINGS) $(FFLAGS) $< -o $@ $(F90_MOD_FLAG) $(BUILD_DIR)

%.mod: ;

# C
# -
$(BUILD_DIR)/%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<


# Goals
# -------

# Compile binary
$(BIN_DIR)/$(PROGRAM): $(OBJS) | $(BIN_DIR)
	$(LD) -o $@ $(LDFLAGS) -I $(BUILD_DIR) $(OBJS) $(LIBS)

$(BIN_DIR) $(BUILD_DIR) $(DEPEND_DIR):
	mkdir -p $@

$(F_DEPEND): $(F_FILES)
	@tools/sfmakedepend --file - --silent --objdir \$$\(BUILD_DIR\) --moddir \$$\(BUILD_DIR\) --depend=mod $^ > $@

$(BIN_DIR)/$(SHARED_LIBRARY): $(LIB_OBJS) | $(BIN_DIR)
	$(LD) -shared -o $@ $(LDFLAGS) -I $(BUILD_DIR) $(LIB_OBJS) $(LIBS)

$(BIN_DIR)/$(STATIC_LIBRARY): $(LIB_OBJS) | $(BIN_DIR)
	$(AR) $(ARFLAGS) $@ $(LIB_OBJS)

# Phonies
# -------
clean:
	rm -fr $(BUILD_DIR)/*
	rm -f $(BIN_DIR)/$(PROGRAM)

cleanall:
	rm -fr $(BUILD_DIR)
	rm -fr $(BIN_DIR)

debug: FFLAGS = $(FFLAGS_DEBUG)
debug: $(BIN_DIR)/$(PROGRAM)

test:
	cd $(TEST_DIR) && pytest -v --tb=line

# Dependencies
# -------------------
$(BUILD_DIR)/printing.o: $(SRCS)

$(F_DEPEND): | $(DEPEND_DIR)
ifneq ($(F_DEPEND),)
include $(F_DEPEND)
endif

$(OBJS): | $(BUILD_DIR)
