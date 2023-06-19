#define RUN_SEARCH      1

#include "Client.hpp"

// Note:
// 1. Change the number of ip address in cmpc_config.h when increasing the nP

int main(int argc, char **argv)
{	
    Client Client(argv);
    
    int bloom_filter_size = 16384;
    int num_files         = 65536;
    
    Client.initialize_ORAM(bloom_filter_size, num_files>>3, nP);
    
#ifndef RUN_SEARCH
    Client.test_update();
#else 
    Client.search_keyword("security");
#endif 
        
    return 0;
}
