#!/bin/csh
# Test suite

make test

set programs = `ls Testing/*.exe`
foreach program ($programs)
   echo $program
   ./$program
end
