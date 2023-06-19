#define RUN_SEARCH      1

#include "Server.hpp"

int main(int argc, char **argv)
{
    Server Server(argv);
    
    Server.initialize_ORAM();

#ifndef RUN_SEARCH
    Server.test_update();
#else 
    Server.process();
#endif 
    
    return 0;
} 


