#include "Util.h"
#include <iostream>
#include <stdlib.h>

Util::Util() {

}
void Util::error(std::string iMessage) {
   std::cout << iMessage << std::endl;
   abort();
}
void Util::warning(std::string iMessage) {
   std::cout << iMessage << std::endl;
   abort();
}
