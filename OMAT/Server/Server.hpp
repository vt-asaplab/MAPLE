#ifndef __SERVER_HPP_
#define __SERVER_HPP_
#include <iostream>
#include <algorithm>
#include <zmq.hpp>
#include <string>
#include <random>
#include <vector>
#include <ctime>
#include <atomic>
#include <list>

// Header file for OMAT 
#include "OMAT.hpp"

using namespace std;

// Constants for spliting database
const int DATABASE_SPLIT        = 32;
const int SERVER_PORT           = 8888;

random_device rd;
mt19937       mt(rd());

uint8_t       **rin;
uint8_t       **rout;
uint8_t       **rtemp;

uint8_t       **ein;
uint8_t       **eout;
uint8_t       **etemp;

uint8_t       *in;
uint8_t       *out;
uint8_t       *temp;

atomic<int> counter;
mutex       mtx;
list<int>   paths;

class Server
{
    public:
        zmq::context_t        *context;
        zmq::socket_t         *socket;

        OMAT                  *keyword_file_table;
        int                   party;

        Server(char **server_info);
        ~Server();
        
        void initialize_ORAM();
#ifdef RUN_SEARCH
        void process();
        void self_test();
#else 
        // For testing ORAM
        void test_update();
#endif   
};

Server::Server(char **server_info)
{
    int party, port;
    parse_party_and_port(server_info, &party, &port);   
    this->party = party;

    keyword_file_table = new OMAT(party, port);

    context = new zmq::context_t(1);
    socket  = new zmq::socket_t(*context, ZMQ_REP);
    
    socket->bind("tcp://*:" + to_string(SERVER_PORT+(party-1)*nP+(party-1)));
}

Server::~Server()
{
    delete keyword_file_table;
}

void Server::initialize_ORAM()
{   
    // Receive ORAM from client
    cout << "Waiting for receiving database from client..." << endl;
    zmq::message_t oram_config;
    socket->recv(&oram_config);
    
    // Reply successfully received
    char *received_announcement = "SUCCESS";
    zmq::message_t reply(strlen(received_announcement));
    memcpy((void*)reply.data(), received_announcement, strlen(received_announcement));
    socket->send(reply);
    
    uint64_t word_block_size  = *(int*)oram_config.data();
    uint64_t file_block_size  = *((int*)oram_config.data() + 1);
    keyword_file_table->initialize(word_block_size<<1, file_block_size<<1);
    
    cout << "INITIALIZING OMAT HAS BEEN DONE!!!" << endl;
    
    uint64_t total_size_oram  = ((word_block_size<<1) + META_DATA_SIZE) * BUCKET_SIZE * keyword_file_table->oram->num_buckets_word + META_DATA_SIZE * BUCKET_SIZE * keyword_file_table->oram->num_buckets_file;
    uint8_t  *oram_info       = new uint8_t[total_size_oram];
    uint8_t  *oram_info_ptr   = oram_info;
    
    ORAM     *oram = keyword_file_table->oram;

    /*
    string db_file = "./DATABASE/DB" + to_string(word_block_size<<3) + "_" + to_string(party);
    FILE *fi = fopen(db_file.c_str(), "rb");

    for(int i = 0; i < oram->num_buckets_word; ++i)
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            fread(oram->meta_data_word[i]+j, META_DATA_SIZE, 1, fi);
            fread(oram->data[i][j], (word_block_size<<1), 1, fi);
        }
    
    for(int i = 0; i < oram->num_buckets_file; ++i)
        for(int j = 0; j < BUCKET_SIZE; ++j)
            fread(oram->meta_data_file[i]+j, META_DATA_SIZE, 1, fi);

    fclose(fi);
    */
    
    for(int i = 0; i < DATABASE_SPLIT; ++i)
    {
        // Receive ORAM database split
        zmq::message_t oram_data;   
        socket->recv(&oram_data);
        memcpy(oram_info_ptr, oram_data.data(), total_size_oram/DATABASE_SPLIT);
        oram_info_ptr += total_size_oram/DATABASE_SPLIT;
        
        // Reply successfully received
        char *received_announcement = "SUCCESS";
        zmq::message_t reply(strlen(received_announcement));
        memcpy((void*)reply.data(), received_announcement, strlen(received_announcement));
        socket->send(reply);
    }

    oram_info_ptr = oram_info;
    for(int i = 0; i < oram->num_buckets_word; ++i)
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            memcpy(oram->meta_data_word[i]+j, oram_info_ptr, META_DATA_SIZE);
            oram_info_ptr += META_DATA_SIZE;
            memcpy(oram->data[i][j], oram_info_ptr, (word_block_size<<1));
            oram_info_ptr += (word_block_size<<1);
        }
    
    for(int i = 0; i < oram->num_buckets_file; ++i)
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            memcpy(oram->meta_data_file[i]+j, oram_info_ptr, META_DATA_SIZE);
            oram_info_ptr += META_DATA_SIZE;
        }

    delete [] oram_info;
    
    cout << "Received " << oram->num_buckets_word << " buckets from client." << endl;
    
    // Buffer for MPC
#ifdef RUN_SEARCH
    uint64_t total_block_size = keyword_file_table->oram->word_block_size/MAX_NUM_CIRCUITS_RETRIEVAL + META_DATA_SIZE;

    rin   = new uint8_t*[MAX_NUM_CIRCUITS_RETRIEVAL];
    rout  = new uint8_t*[MAX_NUM_CIRCUITS_RETRIEVAL];
    rtemp = new uint8_t*[MAX_NUM_CIRCUITS_RETRIEVAL];

    for(int i = 0; i < MAX_NUM_CIRCUITS_RETRIEVAL; ++i)
    {
        rin[i]   = new uint8_t[LEVEL_SIZE+(LEVEL_SIZE*3+total_block_size)*(keyword_file_table->stash_word_size+oram->depth_word*BUCKET_SIZE)+total_block_size];
        rout[i]  = new uint8_t[total_block_size*((keyword_file_table->stash_word_size+oram->depth_word+1)*BUCKET_SIZE)];
        rtemp[i] = new uint8_t[(keyword_file_table->oram->word_block_size<<1)/MAX_NUM_CIRCUITS_RETRIEVAL]; 
    }
    
    ein   = new uint8_t*[MAX_NUM_CIRCUITS_EVICTION];
    eout  = new uint8_t*[MAX_NUM_CIRCUITS_EVICTION];
    etemp = new uint8_t*[MAX_NUM_CIRCUITS_EVICTION];

    for(int i = 0; i < MAX_NUM_CIRCUITS_EVICTION; ++i)
    {
        ein[i]   = new uint8_t[LEVEL_SIZE+(LEVEL_SIZE*3+total_block_size)*(keyword_file_table->stash_word_size+oram->depth_word*BUCKET_SIZE)+total_block_size];
        eout[i]  = new uint8_t[total_block_size*((keyword_file_table->stash_word_size+oram->depth_word+1)*BUCKET_SIZE)];
        etemp[i] = new uint8_t[(keyword_file_table->oram->word_block_size<<1)/MAX_NUM_CIRCUITS_EVICTION]; 
    }
#else   
    uint64_t total_block_size = (keyword_file_table->oram->file_block_size>>3) + META_DATA_SIZE;
    in    = new uint8_t[LEVEL_SIZE+(LEVEL_SIZE*3+total_block_size)*(keyword_file_table->stash_file_size+oram->depth_file*BUCKET_SIZE)+total_block_size];
    out   = new uint8_t[total_block_size*((keyword_file_table->stash_file_size+oram->depth_file+1)*BUCKET_SIZE)];
    temp  = new uint8_t[(keyword_file_table->oram->file_block_size>>3)<<1];
#endif 
}

#ifndef RUN_SEARCH
void Server::test_update()
{
    stringstream ss;
    int path_id, block_id;
    
    while(1)
    {
        for(int r = 0; r < 100; ++r)
        {
            zmq::message_t update_request;
            socket->recv(&update_request);
            ss << update_request.to_string();

            std::cout << "Request from client: " << update_request.to_string() << std::endl;
            
            int file_block_size_in_byte = keyword_file_table->oram->file_block_size >> 3;
            
            if(party == 1)
                memset(temp, 255, file_block_size_in_byte);
            else 
                memset(temp, 0, file_block_size_in_byte);
            
            auto start = clock_start();

            ss >> path_id >> block_id;

            // Retrieve and update stash
            keyword_file_table->retrieve_block_file(in, out, path_id, block_id);        

            zmq::message_t reply(META_DATA_SIZE);
            memcpy(reply.data(), out, META_DATA_SIZE);
            socket->send(reply);

            keyword_file_table->update_stash_file(in, out); 

            // Call eviction 
            int path_id = keyword_file_table->reverse_order_file((keyword_file_table->time_step_file << 1) % (1 << (keyword_file_table->oram->depth_file-1))); 
            keyword_file_table->get_deepest_blocks_file(in, out, path_id);
            keyword_file_table->prepare_deepest_file(in, out);
            keyword_file_table->prepare_target_file(in, out); 
            keyword_file_table->eviction_file(in, out);
            
            // Call eviction
            path_id = keyword_file_table->reverse_order_file(((keyword_file_table->time_step_file << 1) + 1) % (1 << (keyword_file_table->oram->depth_file-1))); 
            keyword_file_table->get_deepest_blocks_file(in, out, path_id);
            keyword_file_table->prepare_deepest_file(in, out);
            keyword_file_table->prepare_target_file(in, out); 
            keyword_file_table->eviction_file(in, out);
            
            // Increase time step for next eviction
            keyword_file_table->increase_time_step_file(); 
            
            for(int i = 0; i < keyword_file_table->stash_file_size; ++i)
                cout << "Stash File[" << i << "]: " << keyword_file_table->stash_file->meta_data[i].block_id << endl;;
            
            cout << "[TIME FOR REQUEST: " << time_from(start) << "]" << endl;
            
            // char *update_notice = "UPDATE FINISHED!";
            // zmq::message_t reply(strlen(update_notice));
            // memcpy((void*)reply.data(), update_notice, strlen(update_notice));
            // socket->send(reply);
        }
    }
}
#else 
void Server::process()
{
    stringstream ss;
    int path_id, block_id;
    
    vector<future<void>> res;
    ThreadPool pool(1);

    // Create a thread for listening and running background eviction
    res.push_back(pool.enqueue([this]() 
	{
        while(1) {
            if(counter > 0) {

                mtx.lock();

                counter--;
                int lpath = paths.front();
                paths.pop_front();
                int rpath = paths.front();
                paths.pop_front();

                mtx.unlock();

                // Perform left-eviction 
                cout << "Eviction on path " << lpath << endl;
                keyword_file_table->get_deepest_blocks_word(ein, eout, lpath);
                keyword_file_table->prepare_deepest_word(ein, eout);
                keyword_file_table->prepare_target_word(ein, eout); 
                keyword_file_table->eviction_word(ein, eout);

                // Perform right-eviction 
                cout << "Eviction on path " << rpath << endl; 
                keyword_file_table->get_deepest_blocks_word(ein, eout, rpath);
                keyword_file_table->prepare_deepest_word(ein, eout);
                keyword_file_table->prepare_target_word(ein, eout); 
                keyword_file_table->eviction_word(ein, eout);
            }
        }
    }));

    while(1)
    {
        zmq::message_t search_request;
        socket->recv(&search_request);
        ss << search_request.to_string();

        std::cout << "Request from client: " << search_request.to_string() << std::endl;

        for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
        {
            if(party == 1)
                memset(rtemp[t], 255, keyword_file_table->oram->word_block_size/MAX_NUM_CIRCUITS_RETRIEVAL);
            else 
                memset(rtemp[t], 0, keyword_file_table->oram->word_block_size/MAX_NUM_CIRCUITS_RETRIEVAL);
        }
        
        auto start = clock_start();

        keyword_file_table->reset_time();
        
        while(!ss.eof())
        {
            ss >> path_id >> block_id;
            if (ss.tellg() == -1) 
            {
                ss.clear();
                break;
            }
            
            // Retrieve and update stash
            keyword_file_table->retrieve_block_word(rin, rout, path_id, block_id); 
            // cout << "Block ID: " << *(int*)rout[0] << endl;
            keyword_file_table->and_2blocks(rin, rout, rtemp);
            keyword_file_table->update_stash_word(rin, rout); 
            
            // for(int i = 0; i < keyword_file_table->stash_word->stash_size; ++i)
            //     cout << "Stash_word[" << i << "]: " << keyword_file_table->stash_word->meta_data[i].block_id << endl;
            
            // Increment counter and add paths for background eviction
            mtx.lock();

            counter++;
            
            int evt_path = keyword_file_table->reverse_order_word((keyword_file_table->time_step_word << 1) % (1 << (keyword_file_table->oram->depth_word-1))); 
            paths.push_back(evt_path);

            evt_path = keyword_file_table->reverse_order_word(((keyword_file_table->time_step_word << 1) + 1) % (1 << (keyword_file_table->oram->depth_word-1))); 
            paths.push_back(evt_path);

            mtx.unlock();

            // Increase time step for next eviction
            keyword_file_table->increase_time_step_word(); 
        }

        keyword_file_table->and_meta_data(rin, rout, rtemp); 
        
        cout << "[PREPROCESSING TIME: " << keyword_file_table->offline_time << "]" << endl;
        
        cout << "[ONLINE PROCESSING TIME: " << keyword_file_table->online_time << "]" << endl;
        
        zmq::message_t reply(keyword_file_table->oram->word_block_size<<5);
        size_t result_size_per_thread = (keyword_file_table->oram->word_block_size<<5)/MAX_NUM_CIRCUITS_RETRIEVAL;

        for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
            memcpy((uint8_t*)reply.data() + t * result_size_per_thread, rout[t], result_size_per_thread);

        socket->send(reply);
    }
}
#endif 

#ifdef RUN_SEARCH
void Server::self_test()
{
    stringstream ss;
    int path_id, block_id;

    while(1)
    {
        zmq::message_t search_request;
        socket->recv(&search_request);
        ss << search_request.to_string();

        std::cout << "Request from client: " << search_request.to_string() << std::endl;

        for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
        {
            if(party == 1)
                memset(rtemp[t], 255, keyword_file_table->oram->word_block_size/MAX_NUM_CIRCUITS_RETRIEVAL);
            else 
                memset(rtemp[t], 0, keyword_file_table->oram->word_block_size/MAX_NUM_CIRCUITS_RETRIEVAL);
        }

        auto start = clock_start();

        ss >> path_id >> block_id;
        if (ss.tellg() == -1) 
        {
            ss.clear();
            break;
        }
        // Retrieve and update stash
        keyword_file_table->retrieve_block_word(rin, rout, path_id, block_id);

        size_t data_size_per_thread = keyword_file_table->oram->word_block_size/MAX_NUM_CIRCUITS_RETRIEVAL;

        zmq::message_t reply(keyword_file_table->oram->word_block_size);

        for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
            memcpy(reply.data() + t * data_size_per_thread, rout[t], data_size_per_thread);

        socket->send(reply);

        /* keyword_file_table->update_stash_word(in, out); 

        // Call eviction 
        int path_id = keyword_file_table->reverse_order_word((keyword_file_table->time_step_word << 1) % (1 << (keyword_file_table->oram->depth_word-1))); 
        
        keyword_file_table->get_deepest_blocks_word(in, out, path_id);
        keyword_file_table->prepare_deepest_word(in, out);
        keyword_file_table->prepare_target_word(in, out); 
        keyword_file_table->eviction_word(in, out);
        
        // Call eviction
        path_id = keyword_file_table->reverse_order_word(((keyword_file_table->time_step_word << 1) + 1) % (1 << (keyword_file_table->oram->depth_word-1)));
        
        keyword_file_table->get_deepest_blocks_word(in, out, path_id);
        keyword_file_table->prepare_deepest_word(in, out);
        keyword_file_table->prepare_target_word(in, out); 
        keyword_file_table->eviction_word(in, out);
        */
        // Increase time step for next eviction
        keyword_file_table->increase_time_step_word(); 
        
        for(int i = 0; i < keyword_file_table->stash_word->stash_size; ++i)
            cout << "Stash_word[" << i << "]: " << keyword_file_table->stash_word->meta_data[i].block_id << endl;
        
        cout << "[TIME FOR REQUEST: " << time_from(start) << "]" << endl;
    }
}
#endif 
#endif