#ifndef __CLIENT_HPP_
#define __CLIENT_HPP_
#include <iostream>
#include <cstring>
#include <string>
#include <random>
#include <ctime>
#include <zmq.hpp>
// Header file for Authenticated Garbling Circuit
#include "emp-tool/emp-tool.h"
#include "emp-agmpc/emp-agmpc.h"
// Other headers for ORAM 
#include "ORAM.hpp"
#include "Position_Map.hpp"
// The hash function used in Bloom Filter
#include "MurmurHash3.hpp"

using namespace std;

// Constants for ORAM procedures
const int nP             = 2;
const int SERVER_PORT    = 8888;

// Constants for Bloom Filter
const int NUM_HASH_FUNC  = 7;

// Constants for spliting database
const int DATABASE_SPLIT = 32;

random_device rd;
mt19937       mt(rd());

int           *result;

int           path_id[NUM_HASH_FUNC];
int           block_id[NUM_HASH_FUNC];

double        client_computation_time;

inline uint64_t nth_hash(uint8_t i, uint64_t hash_a, uint64_t hash_b, uint64_t filter_size)
{
    return (hash_a + i * hash_b) % filter_size;
}

array<uint64_t, 2> mm_hash(const uint8_t* data, size_t len)
{
    array<uint64_t, 2> hash_value;
    MurmurHash3_x64_128(data, len, 0, hash_value.data());
    return hash_value;
}

class Client
{
    public:
        Position_Map    *position_map_word;
        Position_Map    *position_map_file;

        int             word_block_size;
        int             file_block_size;

        zmq::context_t  **context;  
        zmq::socket_t   **socket;

        Client(char **client_info);
        ~Client();
        
        // For initialization
        void initialize_ORAM(int num_blocks, int word_block_size, int num_parties);
        void create_database_shares(uint8_t *real_database, uint8_t **database_shares, uint64_t total_database_size);
        
        // For retrieval procedure
        void request_search(int *path_id, int *block_id);
        void wait_response();

        // General processing
        void search_keyword(string keyword);

        // For testing ORAM
        void retrieve_data_block(int path_id, int block_id);
        void self_test();

        void test_update();
};

Client::Client(char **client_info)
{   
    cout << "Connecting to ORAM Server..." << endl;
    context = new zmq::context_t*[nP];
    socket  = new zmq::socket_t*[nP];
    
    string ip_addresses[nP] = {"127.0.0.1", "127.0.0.1"};
    
    for (int i = 0; i < nP; ++i)
    {
        context[i] = new zmq::context_t(1);
        socket[i]  = new zmq::socket_t(*context[i], ZMQ_REQ);
        string send_address = "tcp://" + ip_addresses[i] + ":" + to_string(SERVER_PORT+i*nP+i);
        cout << "Connecting to " << send_address << " for communication with Server " << (i+1) << " ..." << endl;
        socket[i]->connect(send_address);
    }
}

Client::~Client()
{
    
}

void Client::create_database_shares(uint8_t *real_database, uint8_t **database_shares, uint64_t total_database_size)
{
    uint64_t num_words = total_database_size >> 2;
    for(int n = 0; n < nP; ++n)
    {   
        int *share_ptr = (int*)(database_shares[n]);
        for(uint64_t i = 0; i < num_words; ++i)
            share_ptr[i] = mt();
    }
    
    int *share_ptr = (int*)(database_shares[nP-1]);
    for(uint64_t i = 0; i < num_words; ++i)
    {
        share_ptr[i] = *((int*)real_database + i);
        for(int j = 0; j < nP-1; ++j)
            share_ptr[i] ^= *((int*)(database_shares[j]) + i);
    }
}

void Client::initialize_ORAM(int file_block_size, int word_block_size, int num_parties)
{       
    // Initialize an empty ORAM tree
    this->word_block_size = word_block_size;
    this->file_block_size = file_block_size;
    
    ORAM *oram = new ORAM(file_block_size, word_block_size);
    oram->init_tree();
    
    // Initialize position map
    this->position_map_word = new Position_Map(file_block_size);
    this->position_map_file = new Position_Map(word_block_size<<3);

    int num_paths           = file_block_size>>1;
    int num_path_positions  = oram->depth_word * BUCKET_SIZE;
    
    int  path_id;
    int  path_position;
    bool is_success;

    uint8_t *buffer_data = new uint8_t[word_block_size + META_DATA_SIZE];

    // Insert some dummy data blocks
    FILE *in_omat = fopen("omat.bin", "rb");
    int rv;

    for(int i = 0; i < file_block_size; ++i)
    {
        is_success = false;
        fread(buffer_data, 1, word_block_size, in_omat);
        while(!is_success)
        {
            path_id       = mt() % num_paths;
            path_position = mt() % num_path_positions;
            is_success    = oram->insert_block(i + 1, buffer_data, path_id, path_position);
        }
        this->position_map_word->update_position_map(i, path_id);
        cout << "Insert block id: " << (i + 1) << " into path id: " << path_id << " at position: " << path_position << endl;
    }
    fclose(in_omat); 

    /* Following is the positions of file columns set as default to test updating feature */
    int count_bid          = 8;
    int count_pid          = 0;
    int n_buckets_on_level = 1;
    int n_paths            = (word_block_size<<3)>>1;
    int n_buckets          = 0;
    
    for(int i = 0; i < (word_block_size<<3)/BUCKET_SIZE; ++i)
    {
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            oram->meta_data_file[i][j].block_id = count_bid--; 
            oram->meta_data_file[i][j].path_id  = count_pid * (n_paths / n_buckets_on_level);
            position_map_file->update_position_map(oram->meta_data_file[i][j].block_id-1, oram->meta_data_file[i][j].path_id);
            cout << "Block ID: " << oram->meta_data_file[i][j].block_id << " is on the path: " << oram->meta_data_file[i][j].path_id << endl;
        }
        count_pid++;
        if(count_pid == n_buckets_on_level)
        {
            n_buckets_on_level <<= 1;
            count_pid            = 0;
        }
        n_buckets++;
        if(n_buckets == 4)
        {
            n_buckets  = 0;
            count_bid += 16;
        }
    }

#ifndef RUN_SEARCH 
    num_paths           = (word_block_size<<3)>>1;
    num_path_positions  = oram->depth_file * BUCKET_SIZE;
    int assigned_path_id;
    
    for(int i = 0; i < (word_block_size<<3)/BUCKET_SIZE; ++i)
    {
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            int updated_block_id = oram->meta_data_file[i][j].block_id;
            is_success = false;
            while(!is_success)
            {
                path_id       = mt() % num_paths;
                path_position = mt() % num_path_positions;
                is_success    = oram->swap_block(i, j, path_id, path_position);
            }
            position_map_file->update_position_map(updated_block_id-1, path_id);
        }
    } 
#endif 
    delete [] buffer_data;

    uint64_t total_bucket_size = ((META_DATA_SIZE + (word_block_size<<1)) * BUCKET_SIZE);
    uint64_t total_size_oram   = total_bucket_size * oram->num_buckets_word + META_DATA_SIZE * BUCKET_SIZE * oram->num_buckets_file;
    cout << "Num Buckets: " << oram->num_buckets_word << endl;
    cout << "Total OMAT Size: " << total_size_oram << endl;

    // Send initiliazed ORAM to server
    uint8_t *real_database = new uint8_t[total_size_oram];
    
    uint8_t **database_shares = new uint8_t*[num_parties];
    for(int i = 0; i < num_parties; ++i)
        database_shares[i] = new uint8_t[total_size_oram];

    size_t offset = 0;
    for(int i = 0; i < oram->num_buckets_word; ++i)
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            memcpy(real_database + offset, oram->meta_data_word[i]+j, META_DATA_SIZE);
            offset += META_DATA_SIZE;
            memcpy(real_database + offset, oram->data[i][j], word_block_size<<1);
            offset += (word_block_size<<1);
        }

    for(int i = 0; i < oram->num_buckets_file; ++i)
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            memcpy(real_database + offset, oram->meta_data_file[i]+j, META_DATA_SIZE);
            offset += META_DATA_SIZE;
        }

    cout << "Creating database shares..." << endl;

    create_database_shares(real_database, database_shares, total_size_oram);
    
    delete    oram;
    delete [] real_database;
    /*
    for(int n = 0; n < num_parties; ++n)
    {
        string db_file = "./DATABASE/DB" + to_string(word_block_size<<3) + "_" + to_string(n+1);
        FILE *fo = fopen(db_file.c_str(), "wb");
        fwrite(database_shares[n], total_size_oram, 1, fo);
        fclose(fo);
    } 
    */
    cout << "Sending OMAT configurations..." << endl;
    
    for(int n = 0; n < num_parties; ++n)
    {
        // Send ORAM configuration to server(s)
        int config_size = sizeof(oram->word_block_size) + sizeof(oram->file_block_size);
        zmq::message_t request(config_size);
        memcpy((void*)request.data(), &oram->word_block_size, sizeof(oram->word_block_size));
        memcpy((void*)((uint8_t*)request.data() + sizeof(oram->word_block_size)), &oram->file_block_size, sizeof(oram->file_block_size));
        socket[n]->send(request); 
        
        // Get the reply
        zmq::message_t reply;
        socket[n]->recv(&reply);
        cout << reply.to_string() << endl;
    }
    
    for(int n = 0; n < num_parties; ++n)
    {
        uint8_t *database_ptr = database_shares[n];
        cout << "Uploading ORAM to Server " << (n+1) << "..." << endl;
        for(int m = 0; m < DATABASE_SPLIT; ++m)
        {
            // Send ORAM database to server(s)
            zmq::message_t request(total_size_oram/DATABASE_SPLIT);
            memcpy((void*)request.data(), database_ptr, total_size_oram/DATABASE_SPLIT);
            database_ptr += total_size_oram/DATABASE_SPLIT;
            socket[n]->send(request); 
            
            // Get the reply
            zmq::message_t reply;
            socket[n]->recv(&reply);
            cout << reply.to_string() << endl;
        }
        delete [] database_shares[n];
    }

    delete [] database_shares;

    result = new int[word_block_size<<4];
}

void Client::request_search(int *path_id, int *block_id)
{
    static int shares_id[nP][NUM_HASH_FUNC];

    for(int n = 0; n < nP; ++n)
    {
        auto start = clock_start();

        string request = "";
        for(int i = 0; i < NUM_HASH_FUNC; ++i)
        {
            shares_id[n][i] = 0;
            if(n == nP - 1)
            {
                for(int j = 0; j < nP - 1; ++j)
                    shares_id[n][i] ^= shares_id[j][i];
                shares_id[n][i] ^= block_id[i];
            }
            else 
            {
                shares_id[n][i] = mt();
            }
            request += to_string(path_id[i]) + " " + to_string(shares_id[n][i]) + " ";
        }

        zmq::message_t search_request(request.length());
        memcpy((void*)search_request.data(), request.c_str(), request.length());

        client_computation_time += time_from(start);

        // Send request
        socket[n]->send(search_request);
    }
}

void Client::wait_response()
{
    memset(result, 0, word_block_size<<6);
    for(int n = 0; n < nP; ++n)
    {
        zmq::message_t reply;
        socket[n]->recv(&reply);
        auto start = clock_start();
        int *received_data = (int*)reply.data();
        for(int j = 0; j < word_block_size<<4; ++j)
            result[j] ^= received_data[j];
        client_computation_time += time_from(start);
    } 
}

void Client::search_keyword(string keyword)
{
    auto start_client_time = clock_start();

    array<uint64_t, 2> hash_value = mm_hash((uint8_t*)keyword.c_str(), keyword.length());
    vector<uint32_t> pos;

    int count_total = 0;

    for (int i = 0; i < NUM_HASH_FUNC; i++)
        pos.push_back(nth_hash(i, hash_value[0], hash_value[1], file_block_size));

    for(int i = 0; i < NUM_HASH_FUNC; i++)
    {
        block_id[i] = pos[i] + 1;
        path_id[i]  = position_map_word->path_id[pos[i]];
    }

    client_computation_time = time_from(start_client_time);

    request_search(path_id, block_id);
    wait_response();

    cout << "\n=====================================================\n";
    cout << "Keyword \"" << keyword << "\" appears in these file IDs: ";

    for(int i = 0; i < (word_block_size << 4); ++i)
    {
        if(result[i] != 0)
        {
            count_total++;
            cout << result[i] << " ";
        }
    }
    
    cout << "\nThere are " << count_total << " files containing \"" << keyword << "\" in database.\n";
    cout << "=====================================================\n\n";
    cout << "Client processing time: " << client_computation_time << endl;
}

void Client::retrieve_data_block(int path_id, int block_id)
{
    static int shares_id[nP];

    for(int n = 0; n < nP; ++n) 
    {
        string request = "";
        shares_id[n] = 0;
        
        if(n == nP - 1)
        {
            for(int i = 0; i < nP - 1; ++i)
                shares_id[n] ^= shares_id[i];
            shares_id[n] ^= block_id;
        }
        else 
        {
            shares_id[n] = mt();
        }
        
        request += to_string(path_id) + " " + to_string(shares_id[n]) + " ";
        zmq::message_t request_data_block(request.length());
        memcpy((void*)request_data_block.data(), request.c_str(), request.length());

        // Send request
        socket[n]->send(request_data_block);
    }

    int total_block_size = word_block_size + META_DATA_SIZE;
    memset(result, 0, total_block_size);

    for(int n = 0; n < nP; ++n)
    {
        zmq::message_t reply;
        socket[n]->recv(&reply);
        int *received_data = (int*)reply.data();
        for(int j = 0; j < total_block_size>>2; ++j)
            result[j] ^= received_data[j];
    } 

    int received_block_id = *(int*)result;
    cout << "Received block ID: " << received_block_id << endl;
    
    /* if(received_block_id != block_id)
    {
        cout << "An error occurred!!!" << endl;
        exit(0);
    } */
}

void Client::test_update()
{
    for(int r = 0; r < 10; ++r)
    {
        // for(int i = 1; i <= word_block_size<<3; ++i)
        for(int i = 1; i <= 100; ++i)
        {
            auto start = clock_start();

            int path_id = position_map_file->path_id[i-1];
            static int shares_id[nP];
            
            for(int n = 0; n < nP; ++n)
            {
                shares_id[n] = 0;
                if(n == nP - 1)
                {
                    for(int j = 0; j < nP - 1; ++j)
                        shares_id[n] ^= shares_id[j];
                    shares_id[n] ^= i;
                }
                else 
                {
                    shares_id[n] = mt();
                }
                string request = to_string(path_id) + " " + to_string(shares_id[n]) + " ";

                zmq::message_t update_request(request.length());
                memcpy((void*)update_request.data(), request.c_str(), request.length());

                // Send request
                socket[n]->send(update_request);
            }

            memset(result, 0, META_DATA_SIZE);

            double client_processing_time = time_from(start);

            for(int n = 0; n < nP; ++n)
            {
                zmq::message_t reply;
                socket[n]->recv(&reply);
                auto start = clock_start();
                int *received_data = (int*)reply.data();
                for(int j = 0; j < META_DATA_SIZE>>2; ++j)
                    result[j] ^= received_data[j];
                client_processing_time += time_from(start);
            } 

            int received_block_id = *(int*)result;

            cout << "Received block ID: " << received_block_id << endl;
            cout << "Client processing time: " << client_processing_time << endl;
            
            if(received_block_id != i)
            {
                cout << "An error occurred!!!" << endl;
                exit(0);
            }
        }
    }
}

void Client::self_test()
{
    int requested_block_id;
    for(int k = 0; k < 100; ++k)
    {
        cout << "Round test #" << (k + 1) << endl;
        for(int i = 0; i < file_block_size; ++i)
        {
            requested_block_id = i + 1;
            cout << "Request block ID: " << requested_block_id << endl;

            int path_id = position_map_word->path_id[requested_block_id - 1];
            retrieve_data_block(position_map_word->path_id[requested_block_id - 1], requested_block_id);
        }
    }
}

#endif 
