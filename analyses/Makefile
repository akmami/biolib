FASTA := ../data/chm13v2.0.chr1.fa

# Arguments
KMER_SIZES   := $(shell seq 6 15)
WINDOW_SIZES := $(shell seq 4 10)

SKIP_KMER 	:= 3
SKIP_WINDOW := 2

BLEND_K_HIFI := 19
BLEND_W_HIFI := 50
BLEND_K_SHORT := 21
BLEND_W_SHORT := 11
BLEND_K_BEST := 15
BLEND_W_BEST := 10

K_VALUES := 15 17 19 21 23 25 27 29 31
W_VALUES := $(shell seq 10 50)

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
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -Wpedantic
BLENDLIB = -L ../lib -lblend 
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
faHmin: check_fasta mkdir_bin mkdir_out
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
				date > ../out/fa-hmin/$(FAHMIN)-$$k-$$w-output.txt; \
				$(TIME) ../bin/$(FAHMIN) $(FASTA) $$k $$w >> ../out/fa-hmin/$(FAHMIN)-$$k-$$w-output.txt 2>&1; \
		done; \
	done

######################################################################
# FUZZY SEEDS (BLEND)
######################################################################
faBlend: check_fasta mkdir_bin mkdir_out
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -c $(BLEND)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(BLEND) $(BLEND).o $(BLENDLIB)
	@rm $(BLEND).o
	@mv $(BLEND) $(EXECUTABLE_DIR)

	@mkdir -p ../out/fa-blend
	@echo "$(FASTA) -> ../out/fa-blend/$(BLEND)-output.txt"
	@date > ../out/fa-blend/$(BLEND)-output.txt
	@echo "" >> ../out/fa-blend/$(BLEND)-output.txt

	@$(foreach k,$(K_VALUES), \
		$(foreach w,$(W_VALUES), \
			echo "Running k=$(k), w=$(w)" >> ../out/fa-blend/$(BLEND)-output.txt; \
			$(TIME) ../bin/$(BLEND) $(FASTA) $(k) $(w) >> ../out/fa-blend/$(BLEND)-output.txt 2>&1; \
			echo "" >> ../out/fa-blend/$(BLEND)-output.txt; \
		) \
	)


######################################################################
# ERROR TOLERANCE EXPERIMENT
######################################################################
errSim: mkdir_bin mkdir_out
	$(GXX) $(CXXFLAGS) -mavx2 -msse4.1 -c $(SIMULATE)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(SIMULATE) $(SIMULATE).o $(BLENDLIB)
	@rm $(SIMULATE).o
	@mv $(SIMULATE) $(EXECUTABLE_DIR)
	@mkdir -p ../out/err-sim
	@mkdir -p ../data/sim
	@rm -rf ../data/sim/*
	@echo "Output -> ../out/err-sim/$(SIMULATE)-output.txt"
	@date > $(OUT_DIR)/err-sim/$(SIMULATE)-output.txt
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
	$(TIME) ../bin/$(SIMULATE) -g -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
	$(TIME) ../bin/$(SIMULATE) -a --blend --short -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
	$(TIME) ../bin/$(SIMULATE) -a --blend --hifi -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
	$(TIME) ../bin/$(SIMULATE) -a --blend -k $(BLEND_K_BEST) -w $(BLEND_W_BEST) -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
	$(TIME) ../bin/$(SIMULATE) -a --h-min -d ../data/sim >> ../out/err-sim/$(SIMULATE)-output.txt 2>&1
	@echo "" >> ../out/err-sim/$(SIMULATE)-output.txt
