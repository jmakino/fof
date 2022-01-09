PS_PATH = -I../fdps/fdps-devel/src

CC = g++
#CC = mpicxx
CFLAGS = -std=c++17 -O3
CFLAGS += -Wall
CFLAGS += -ffast-math
CFLAGS += -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL


# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

SRC = fof.cpp
PROGRAM = $(SRC:%.cpp=%)

$(PROGRAM):$(SRC)
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -o $@ $<

clean:
	rm -f $(PROGRAM)

distclean: clean
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
