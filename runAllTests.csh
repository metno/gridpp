#!/bin/csh
# Test suite

make test

set programs = `ls testing/*.exe`
foreach program ($programs)
   ./$program
end
