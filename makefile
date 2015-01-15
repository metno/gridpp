CC      = g++
#CFLAGS  = -g -pg #-fprofile-arcs
CFLAGS  = -O3 -fopenmp
SRC     = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
CALSRC  = $(wildcard Calibrator/*.cpp)
FILESRC = $(wildcard File/*.cpp)
ALLOBJS = $(SRC:.cpp=.o) $(CALSRC:.cpp=.o) $(FILESRC:.cpp=.o)
ALLHEADERS = $(ALLOBJS:.o=.h)
COREOBJS= $(filter-out Test.o,$(ALLOBJS))
TESTOBJS= $(filter-out PrecipCal.o,$(ALLOBJS))
IFLAGS  = -I/usr/include/ -I/usr/local/boost/include/
LIBS    = -lnetcdf_c++ -lgtest
LFLAGS  = -L/usr/lib -L/usr/local/boost/lib/ -L/home/thomasn/local/lib/
INCS    = makefile $(HEADERS)

default: precipCal.exe

%.o : %.cpp $(INCS)
	$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@

precipCal.exe: $(COREOBJS) makefile
	echo $(CALOBJS)
	$(CC) $(CFLAGS) $(LFLAGS) $(COREOBJS) $(LIBS) -o $@

test.exe: $(TESTOBJS) makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(TESTOBJS) $(LIBS) -o $@

clean: 
	rm *.o */*.o gmon.out *.gcda precipCal.exe test.exe

tags:
	cd src; ctags -R --c++-kinds=+pl --fields=+iaS --extra=+q -f tags ./ ./*
