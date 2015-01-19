CC      = g++
CFLAGS  = -g -pg -rdynamic #-fprofile-arcs
#CFLAGS  = -O3 -fopenmp
SRC     = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
CALSRC  = $(wildcard Calibrator/*.cpp)
FILESRC = $(wildcard File/*.cpp)
DOWNSRC = $(wildcard Downscaler/*.cpp)
ALLOBJS = $(SRC:.cpp=.o) $(CALSRC:.cpp=.o) $(FILESRC:.cpp=.o) $(DOWNSRC:.cpp=.o)
ALLHEADERS = $(ALLOBJS:.o=.h)
COREOBJS= $(filter-out Test.o,$(ALLOBJS))
TESTS   = $(wildcard Testing/*.cpp)
TESTEXE = $(TESTS:.cpp=.exe)
IFLAGS  = -I/usr/include/ -I/usr/local/boost/include/
LIBS    = -lnetcdf_c++ -lgtest
LFLAGS  = -L/usr/lib -L/usr/local/boost/lib/ -L/home/thomasn/local/lib/
INCS    = makefile $(HEADERS)

.PHONY: tags

default: precipCal.exe

%.o : %.cpp $(INCS)
	$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@

precipCal.exe: $(COREOBJS) Driver/PrecipCal.o makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(COREOBJS) Driver/PrecipCal.o $(LIBS) -o $@

statkraft.exe: $(COREOBJS) Driver/Statkraft.o makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(COREOBJS) Driver/Statkraft.o $(LIBS) -o $@

nn.exe: $(COREOBJS) Driver/TestNearestNeighbours.o makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(COREOBJS) Driver/TestNearestNeighbours.o $(LIBS) -o $@

test.exe: $(TESTOBJS) makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(TESTOBJS) $(LIBS) -o $@

test: $(TESTEXE)

Testing/%.exe: Testing/%.cpp $(INCS) $(COREOBJS)
	$(CC) $(CFLAGS) $(COREOBJS) $< $(LFLAGS) -lgtest $(LIBS) -o $@

clean: 
	rm *.o */*.o gmon.out *.gcda precipCal.exe test.exe Testing/*.exe Testing/*.o

tags:
	ctags -R --c++-kinds=+pl --fields=+iaS --extra=+q -f tags ./*.h ./*.cpp */*.h */*.cpp
