CC      	= g++
CFLAGS  	= -g -pg -rdynamic -fprofile-arcs -ftest-coverage -coverage
#CFLAGS  = -O3 -fopenmp
SRCDIR   = src/
BUILDDIR = build
IFLAGS  	= -I/usr/include/ -I/usr/local/boost/include/
LIBS    	= -lnetcdf_c++ -lgtest
LFLAGS  	= -L/usr/lib -L/usr/local/boost/lib/ -L/home/thomasn/local/lib/
EXE     	= postprocess.exe

CORESRC 	= $(wildcard src/*.cpp)
CALSRC  	= $(wildcard src/Calibrator/*.cpp)
FILESRC 	= $(wildcard src/File/*.cpp)
DOWNSRC 	= $(wildcard src/Downscaler/*.cpp)
DRVSRC  	= src/Driver/PostProcess.cpp
DRVOBJ  	= $(BUILDDIR)/Driver/PostProcess.o
SRC     	= $(CORESRC) $(CALSRC) $(FILESRC) $(DOWNSRC)
HEADERS 	= $(SRC:.cpp=.h)
OBJ0    	= $(patsubst src/%,$(BUILDDIR)/%,$(SRC))
OBJ     	= $(OBJ0:.cpp=.o)
TESTSRC 	= $(wildcard src/Testing/*.cpp)
TESTS0  	= $(patsubst src/Testing%,testing%,$(TESTSRC))
TESTS   	= $(TESTS0:.cpp=.exe)
INCS    	= makefile $(HEADERS)

.PHONY: tags count coverage doxygen

default: $(EXE)

$(BUILDDIR):
	@mkdir build build/Calibrator build/Downscaler build/File build/Driver build/Testing

$(BUILDDIR)/%.o : src/%.cpp $(INCS)
	$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@

postprocess.exe: $(OBJ) $(DRVOBJ) makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(OBJ) $(DRVOBJ) $(LIBS) -o $@

test: $(TESTS)

testing/%.exe: $(BUILDDIR)/Testing/%.o $(INCS) $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $< $(LFLAGS) -lgtest $(LIBS) -o $@

count:
	@wc src/*.h src/*.cpp src/*/*.h src/*/*.cpp -l | tail -1

clean: 
	rm -rf build/*.o build/*/*.o gmon.out $(EXE) testing/*.exe\
		*.gcno build/*.gcda build/*.gcno build/*/*.gcda build/*/*.gcno\
		coverage/* coverage.* 

tags:
	ctags -R --c++-kinds=+pl --fields=+iaS --extra=+q -f tags ./*.h ./*.cpp */*.h */*.cpp

doxygen:
	doxygen doxygen/config

coverage:
	#rm -f build/*.gcno build/*.gcda build/*/*.gcno build/*/*.gcda
	lcov -b . -c -i -d . -o coverage.init
	csh runAllTests.csh
	lcov -b . -c -d . -o coverage.run
	lcov -a coverage.init -a coverage.run -o coverage.total
	lcov -e coverage.total "`pwd`/*" -o coverage.total.filtered
	lcov -r coverage.total.filtered "`pwd`/*Testing*" -o coverage.info
	genhtml -o ./coverage/ coverage.info

depend:
	makedepend -- $(CFLAGS) -- $(SRC)
