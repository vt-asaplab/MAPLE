#include <iostream>
#include "ORAM.hpp"

using namespace std;

class Stash
{
    public:
        int             stash_size;
        uint8_t         **data;
        meta_data_t     *meta_data;

        Stash(int stash_size, int block_size)
        {
            this->stash_size = stash_size;
            this->meta_data  = new meta_data_t[stash_size];
            this->data       = new uint8_t*[stash_size];
            for(size_t i = 0; i < stash_size; i++)
            {
                this->meta_data[i].block_id = 0;
                this->meta_data[i].path_id  = 0;
                this->data[i] = new uint8_t[block_size];
                memset(this->data[i], 0, block_size);
            }
        }

        ~Stash()
        {
            if(this->meta_data != nullptr)
                delete [] meta_data;
            if(this->data != nullptr)
            {
                for(size_t i = 0; i < stash_size; i++)
                    delete [] data[i];
                delete [] data;
            }
        }
};
