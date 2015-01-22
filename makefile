CC      = g++
#CFLAGS  = -g -pg -rdynamic #-fprofile-arcs
CFLAGS  = -O3 -fopenmp
SRC     = $(wildcard src/*.cpp)
HEADERS = $(wildcard src/*.h)
CALSRC  = $(wildcard src/Calibrator/*.cpp)
FILESRC = $(wildcard src/File/*.cpp)
DOWNSRC = $(wildcard src/Downscaler/*.cpp)
ALLOBJS = $(SRC:.cpp=.o) $(CALSRC:.cpp=.o) $(FILESRC:.cpp=.o) $(DOWNSRC:.cpp=.o)
ALLHEADERS = $(ALLOBJS:.o=.h)
COREOBJS= $(filter-out Test.o,$(ALLOBJS))
TESTS   = $(wildcard src/Testing/*.cpp)
TESTE   = $(patsubst src/Testing%,testing%,$(TESTS))
TESTEXE = $(TESTE:.cpp=.exe)
IFLAGS  = -I/usr/include/ -I/usr/local/boost/include/
LIBS    = -lnetcdf_c++ -lgtest
LFLAGS  = -L/usr/lib -L/usr/local/boost/lib/ -L/home/thomasn/local/lib/
INCS    = makefile $(HEADERS)

.PHONY: tags count

default: precipCal.exe

%.o : %.cpp $(INCS)
	$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@

postprocess.exe: $(COREOBJS) src/Driver/PostProcess.o makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(COREOBJS) src/Driver/PostProcess.o $(LIBS) -o $@

nn.exe: $(COREOBJS) src/Driver/TestNearestNeighbours.o makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(COREOBJS) src/Driver/TestNearestNeighbours.o $(LIBS) -o $@

test: $(TESTEXE)

testing/%.exe: src/Testing/%.cpp $(INCS) $(COREOBJS)
	$(CC) $(CFLAGS) $(COREOBJS) $< $(LFLAGS) -lgtest $(LIBS) -o $@

count:
	@wc src/*.h src/*.cpp src/*/*.h src/*/*.cpp -l | tail -1

clean: 
	rm *.o */*.o gmon.out *.gcda precipCal.exe test.exe testing/*.exe testing/*.o

tags:
	ctags -R --c++-kinds=+pl --fields=+iaS --extra=+q -f tags ./*.h ./*.cpp */*.h */*.cpp
