FASTA := ../data/chm13v2.0.chr1.fa

# Arguments
KMER_SIZES   := $(shell seq 6 15)
WINDOW_SIZES := $(shell seq 4 10)

SKIP_KMER 	:= 3
SKIP_WINDOW := 2

BLEND_K_HIFI := 21
BLEND_W_HIFI := 11

# Programs
# fa
FAHMIN := minimizers
BLEND := fuzzy-seeds
SIMULATE := simulation

# Directories
CURRENT_DIR := $(shell pwd)
DATA_DIR := $(CURRENT_DIR)/../data
EXECUTABLE_DIR := $(CURRENT_DIR)/../bin
OUT_DIR := $(CURRENT_DIR)/../out

# Extension
CXX := .cpp

# Compiler
GXX := g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -Wpedantic
TIME := /usr/bin/time -v

.PHONY: all %

%:
	@:

all: mkdir_bin mkdir_out

mkdir_bin:
	@if [ ! -d "$(EXECUTABLE_DIR)" ]; then \
		echo "$(EXECUTABLE_DIR) does not exist. Creating..."; \
		mkdir $(EXECUTABLE_DIR); \
	fi

mkdir_out:
	@if [ ! -d "$(OUT_DIR)" ]; then \
		echo "$(OUT_DIR) does not exist. Creating..."; \
		mkdir $(OUT_DIR); \
	fi

check_fasta:
	@if [ -z "$(FASTA)" ]; then \
        echo "FASTA is not set or is empty. Exiting"; \
		exit 1; \
    fi
	@if [ ! -f $(FASTA) ]; then \
		echo "$(FASTA) does not exist. Exiting."; \
		exit 1; \
	fi

clean:
	@if [ -d "$(OUT_DIR)" ]; then \
		echo "Removing $(OUT_DIR)"; \
		rm -r $(OUT_DIR); \
	fi
	@if [ -d "$(EXECUTABLE_DIR)" ]; then \
		echo "Removing $(EXECUTABLE_DIR)"; \
		rm -r $(EXECUTABLE_DIR); \
	fi
	@echo "Cleaned."

######################################################################
# HIERARCHICAL MINIMIZER
######################################################################
faHmin: mkdir_bin mkdir_out
	$(GXX) $(CXXFLAGS) -c $(FAHMIN)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FAHMIN) $(FAHMIN).o
	@rm $(FAHMIN).o
	@mv $(FAHMIN) $(EXECUTABLE_DIR)
	@mkdir -p ../out/fa-hmin
	@for k in $(KMER_SIZES); do \
		if [ $$(( $$k % $(SKIP_KMER) )) -ne 0 ]; then continue; fi; \
			for w in $(WINDOW_SIZES); do \
				if [ $$(( $$w % $(SKIP_WINDOW) )) -ne 0 ]; then continue; fi; \
				echo "$(FASTA) k=$$k w=$$w -> ../out/fa-hmin/$(FAHMIN)-$$k-$$w-output.txt"; \
				date > $(OUT_DIR)/fa-hmin/$(FAHMIN)-$$k-$$w-output.txt; \
				$(TIME) $(EXECUTABLE_DIR)/$(FAHMIN) $(FASTA) $$k $$w >> $(OUT_DIR)/fa-hmin/$(FAHMIN)-$$k-$$w-output.txt 2>&1; \
		done; \
	done
	
######################################################################
# FUZZY SEEDS (BLEND)
######################################################################
faBlend: mkdir_bin mkdir_out
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -c $(BLEND)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(BLEND) $(BLEND).o
	@rm $(BLEND).o
	@mv $(BLEND) $(EXECUTABLE_DIR)
	@mkdir -p ../out/fa-blend
	echo "$(FASTA) -> ../out/fa-blend/$(BLEND)-k$(BLEND_K_HIFI)-w$(BLEND_W_HIFI)-output.txt"
	date > $(OUT_DIR)/fa-blend/$(BLEND)-k$(BLEND_K_HIFI)-w$(BLEND_W_HIFI)-output.txt
	$(TIME) $(EXECUTABLE_DIR)/$(BLEND) $(FASTA) $(BLEND_K_HIFI) $(BLEND_W_HIFI) >> $(OUT_DIR)/fa-blend/$(BLEND)-k$(BLEND_K_HIFI)-w$(BLEND_W_HIFI)-output.txt 2>&1
	
######################################################################
# ERROR TOLERANCE EXPERIMENT
######################################################################
errSim: mkdir_bin mkdir_out
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -c $(SIMULATE)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(SIMULATE) $(SIMULATE).o
	@rm $(SIMULATE).o
	@mv $(SIMULATE) $(EXECUTABLE_DIR)
	@rm -rf ../data/sim/*
	@mkdir -p ../out/err-sim
	echo "Output -> ../out/err-sim/$(SIMULATE)-short-blend-output.txt"
	date > $(OUT_DIR)/err-sim/$(SIMULATE)-short-blend-output.txt
 	$(TIME) $(EXECUTABLE_DIR)/$(SIMULATE) -g -d ../data/sim >> $(OUT_DIR)/err-sim/$(SIMULATE)-short-blend-output.txt 2>&1
	$(TIME) $(EXECUTABLE_DIR)/$(SIMULATE) -a --blend --short -d ../data/sim >> $(OUT_DIR)/err-sim/$(SIMULATE)-short-blend-output.txt 2>&1
	$(TIME) $(EXECUTABLE_DIR)/$(SIMULATE) -a --blend --hifi -d ../data/sim >> $(OUT_DIR)/err-sim/$(SIMULATE)-short-blend-output.txt 2>&1
	$(TIME) $(EXECUTABLE_DIR)/$(SIMULATE) -a --h-min -d ../data/sim >> $(OUT_DIR)/err-sim/$(SIMULATE)-short-blend-output.txt 2>&1
