OBJS := main.o memory.o utilities.o Ax.o cg.o \
        check_memory.o nodesused.o diag_time.o \
        outputs.o

EXECUTABLE := poisson

# choose compiler:
CC := mpiicx
#CC := mpicc
#CC := gcc

# choose flags:
# flags for Intel compiler icc on taki:
#CFLAGS := -std=c99 -Wall #-mkl -static# -qopenmp  ## NO optimization
#CFLAGS := -O3 -std=c99 -Wall #-mkl -static# -qopenmp 
#CFLAGS := -qopt-zmm-usage=high -xCORE-AVX512 -O3 -std=c99 -Wall -qmkl # -qopenmp 
CFLAGS := -xCORE-AVX512 -O3 -std=c99 -Wall -qmkl # -qopenmp 
#CFLAGS := -axSKYLAKE-AVX512 -O3 -std=c99 -Wall -mkl # -qopenmp 
#CFLAGS := -O3 -std=c99 -Wall -mkl -qopenmp
# flags for GNU compiler gcc anywhere:
#CFLAGS := -O3 -std=c99 -Wall # -Wno-unused-variable

DEFS := -DPARALLEL # -DBLAS
INCLUDES :=
LDFLAGS := -lm

# Used for executable tracking
BUILD_DATE  := $(shell date +%Y-%m-%d::T%H:%M:%S)
BUILD_DIR := $(shell pwd)
DEFS += -DBUILD_DATE=\"$(BUILD_DATE)\" -DBUILD_DIR=\"$(BUILD_DIR)\"

%.o: %.c %.h
	$(CC) $(CFLAGS) $(DEFS) $(INCLUDES) -c $< -o $@

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(DEFS) $(INCLUDES) $(OBJS) -o $@ $(LDFLAGS)

clean:
	-rm -f *.o $(EXECUTABLE)

