#ifndef __ORAM_HPP_
#define __ORAM_HPP_
#include <iostream>
#include <cstring>
#include <cmath>

#define BUCKET_SIZE             (2)
#define BLOCK_ID_SIZE           (sizeof(int))
#define PATH_ID_SIZE            (sizeof(int))
#define META_DATA_SIZE          (BLOCK_ID_SIZE + PATH_ID_SIZE)

struct meta_data_t
{
    int block_id;
    int path_id;
};

class ORAM
{
    public:
        int          depth_word;
        int          depth_file;

        int          num_buckets_word;
        int          num_buckets_file;

        int          word_block_size;
        int          file_block_size;

        meta_data_t  **meta_data_word;
        meta_data_t  **meta_data_file;
        uint8_t      ***data;
        
        ORAM();
        ORAM(int file_block_size, int word_block_size);

        ~ORAM();
        
        void init_tree();
        bool insert_block(int block_id, uint8_t *data, int path_id, int path_position);
        void get_path_word(int path_id, int *full_path);
        void get_path_file(int path_id, int *full_path);
};

ORAM::ORAM()
{
    this->depth_word        = 0;
    this->depth_file        = 0;
    this->num_buckets_word  = 0;
    this->num_buckets_file  = 0;
    this->word_block_size   = 0;
    this->file_block_size   = 0;
    this->data              = nullptr;
    this->meta_data_word    = nullptr;
    this->meta_data_file    = nullptr;
}

ORAM::ORAM(int file_block_size, int word_block_size)
{
    this->depth_word        = ceil(log(file_block_size>>1)/log(2));
    this->depth_file        = ceil(log((word_block_size>>1)<<3)/log(2));
    this->num_buckets_word  = file_block_size - 1;
    this->num_buckets_file  = ((word_block_size>>1)<<3) - 1;
    this->word_block_size   = word_block_size;
    this->file_block_size   = file_block_size;
    this->data              = nullptr;
    this->meta_data_word    = nullptr;
    this->meta_data_file    = nullptr;
}

ORAM::~ORAM()
{
    if(this->meta_data_word != nullptr)
    {
        for(int i = 0; i < this->num_buckets_word; ++i)
        {
            delete [] meta_data_word[i];
        }
        delete [] meta_data_word;
    }
    
    if(this->meta_data_file != nullptr)
    {
        for(int i = 0; i < this->num_buckets_file; ++i)
        {
            delete [] meta_data_file[i];
        }
        delete [] meta_data_file;
    }

    if(this->data != nullptr)
    {
        for(int i = 0; i < this->num_buckets_word; ++i)
        {
            for(int j = 0; j < BUCKET_SIZE; ++j)
                delete [] data[i][j];
            delete [] data[i];
        }
        delete [] data;
    }
}

void ORAM::get_path_word(int path_id, int *full_path)
{
    int current_id = 0;
    full_path[0]   = 0;
    int count      = 0;

    while(count < this->depth_word)
    {
        full_path[count++] = current_id;

        path_id <<= 1;
        current_id <<= 1;
        if(path_id & (1 << (this->depth_word - 1)))
            current_id += 2;
        else 
            current_id += 1;
    }
}

void ORAM::get_path_file(int path_id, int *full_path)
{
    int current_id = 0;
    full_path[0]   = 0;
    int count      = 0;
    
    while(count < this->depth_file)
    {
        full_path[count++] = current_id;

        path_id <<= 1;
        current_id <<= 1;
        if(path_id & (1 << (this->depth_file - 1)))
            current_id += 2;
        else 
            current_id += 1;
    }
}

void ORAM::init_tree()
{
    // Compute other parameters
    this->depth_word        = ceil(log(file_block_size>>1)/log(2));
    this->depth_file        = ceil(log((word_block_size>>1)<<3)/log(2));
    this->num_buckets_word  = (file_block_size>>1) - 1;
    this->num_buckets_file  = ((word_block_size>>1)<<3) - 1;

    // Initialize meta-data word
    this->meta_data_word = new meta_data_t*[this->num_buckets_word];
    for(size_t i = 0; i < this->num_buckets_word; ++i)
    {
        this->meta_data_word[i] = new meta_data_t[BUCKET_SIZE];
        for(int j = 0; j < BUCKET_SIZE; j++)
        {
            this->meta_data_word[i][j].block_id = 0;
            this->meta_data_word[i][j].path_id  = 0;
        }
    }
    
    // Initialize meta-data file
    this->meta_data_file= new meta_data_t*[this->num_buckets_file];
    for(size_t i = 0; i < this->num_buckets_file; ++i)
    {
        this->meta_data_file[i] = new meta_data_t[BUCKET_SIZE];
        for(int j = 0; j < BUCKET_SIZE; j++)
        {
            this->meta_data_file[i][j].block_id = 0;
            this->meta_data_file[i][j].path_id  = 0;
        }
    }

    // Initialize a tree of data block  
    this->data = new uint8_t**[this->num_buckets_word];
    for(int i = 0; i < this->num_buckets_word; ++i)
    {
        data[i] = new uint8_t*[BUCKET_SIZE];
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            data[i][j] = new uint8_t[word_block_size];
            memset(data[i][j], 0, word_block_size);
        }
    }
}

bool ORAM::insert_block(int block_id, uint8_t *data, int path_id, int path_position)
{
    int current_id = 0;
    int current_position = 0;
    
    while(current_position + BUCKET_SIZE <= path_position)
    {
        path_id <<= 1;
        current_position += BUCKET_SIZE;
        current_id <<= 1;
        if(path_id & (1 << (this->depth_word - 1)))
            current_id += 2;
        else 
            current_id += 1;
    }

    int offset = path_position - current_position;
    
    if(this->meta_data_word[current_id][offset].block_id != 0)
        return false;

    this->meta_data_word[current_id][offset].block_id = block_id;
    this->meta_data_word[current_id][offset].path_id  = path_id;

    memcpy(this->data[current_id][offset], data, word_block_size);

    return true;
}

#endif 