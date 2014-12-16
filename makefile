CC      = g++
#CFLAGS  = -g -pg #-fprofile-arcs
CFLAGS  = -O3 -fopenmp
SRC     = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
ALLOBJS = $(SRC:.cpp=.o)
COREOBJS= $(filter-out Test.o,$(ALLOBJS))
TESTOBJS= $(filter-out PrecipCal.o,$(ALLOBJS))
IFLAGS  = -I/usr/include/
LIBS    = -lnetcdf_c++
LFLAGS  = -L/usr/lib -L/usr/local/boost/lib/
INCS    = makefile $(HEADERS)

default: precipCal.exe

%.o : %.cpp $(HEADERS) makefile
	$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@

precipCal.exe: $(COREOBJS) makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(COREOBJS) $(LIBS) -o $@

test.exe: $(TESTOBJS) makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(TESTOBJS) $(LIBS) -o $@

clean: 
	rm *.o gmon.out *.gcda precipCal.exe test.exe
