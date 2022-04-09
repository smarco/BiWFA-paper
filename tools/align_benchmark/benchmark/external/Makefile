###############################################################################
# Definitions
###############################################################################
FOLDER_WFA=../../../..
FOLDER_ROOT=../..
FOLDER_BUILD=$(FOLDER_ROOT)/build

###############################################################################
# Modules
###############################################################################
MODULES=benchmark_bitpal \
        benchmark_daligner \
        benchmark_diffutils \
        benchmark_edlib \
        benchmark_gaba \
        benchmark_ksw2 \
        benchmark_lh3 \
        benchmark_wfalm_bindings \
        benchmark_wfalm
        
SRCS=$(addsuffix .c, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))

LI_FLAGS=-I$(FOLDER_ROOT) -I$(FOLDER_WFA)

###############################################################################
# Rules
###############################################################################
all: $(OBJS)

$(FOLDER_BUILD)/benchmark_diffutils.o : benchmark_diffutils.c
	$(CC) -Wno-maybe-uninitialized $(CC_FLAGS) $(LI_FLAGS) -c $< -o $@

$(FOLDER_BUILD)/benchmark_gaba.o : benchmark_gaba.c
	$(CC) -march=native $(CC_FLAGS) $(LI_FLAGS) -c $< -o $@

$(FOLDER_BUILD)/benchmark_wfalm_bindings.o : benchmark_wfalm_bindings.cpp
	$(CPP) -std=c++14 $(CC_FLAGS) -Wno-sign-compare $(LI_FLAGS) -c $< -o $@

# General building rule
$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(CC_FLAGS) $(LI_FLAGS) -c $< -o $@
	
