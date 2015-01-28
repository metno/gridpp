CC      	= g++
CFLAGS_D = -g -pg -rdynamic -fprofile-arcs -ftest-coverage -coverage -DDEBUG
CFLAGS_O = -O3 -fopenmp
IFLAGS  	= -I/usr/include/
LIBS_D   = -lnetcdf_c++ -lgtest
LIBS_O   = -lnetcdf_c++
LFLAGS   = -L/usr/lib
DEBUG    = 0  # Set to 1 to compile with debug flags, otherwise compile with optimization

ifeq ($(DEBUG), 1)
CFLAGS=$(CFLAGS_D) 
LIBS=$(LIBS_D)
else
CFLAGS=$(CFLAGS_O)
LIBS=$(LIBS_O)
endif
SRCDIR   = src/
BUILDDIR = build
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

$(BUILDDIR)/%.E : src/%.cpp $(INCS)
	$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@ -E

postprocess.exe: $(OBJ) $(DRVOBJ) makefile
	$(CC) $(CFLAGS) $(LFLAGS) $(OBJ) $(DRVOBJ) $(LIBS) -o $@

test: $(TESTS)

testing/%.exe: $(BUILDDIR)/Testing/%.o $(INCS) $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $< $(LFLAGS) -lgtest $(LIBS) -o $@

count:
	@wc src/*.h src/*.cpp src/*/*.h src/*/*.cpp -l | tail -1

clean: 
	rm -rf build/*.o build/*/*.o build/*.E build/*/*.E gmon.out $(EXE) testing/*.exe\
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
