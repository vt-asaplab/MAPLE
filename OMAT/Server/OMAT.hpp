#ifndef __OMAT_HPP_
#define __OMAT_HPP_

#include <iostream>
#include <algorithm>
#include <zmq.hpp>
#include <string>
#include <random>
#include <vector>
#include <ctime>
// Header file for Authenticated Garbling Circuit
#include <emp-tool/emp-tool.h>
#include <emp-agmpc/emp-agmpc.h>

// Header files for ORAM
#include "ORAM.hpp"
#include "Stash.hpp"
// Header file for building Bristol Fashion circuits
#include "Build_Circuit.hpp"

using namespace std;

const int nP                            = 2;

const int THREAD_PER_CIRCUIT            = 1;

const int MAX_NUM_CIRCUITS_RETRIEVAL    = 128;
const int CIRCUITS_PER_THREAD_RETRIEVAL = 32;

const int MAX_NUM_CIRCUITS_EVICTION     = 128;
const int CIRCUITS_PER_THREAD_EVICTION  = 32;

const int MAX_NUM_THREAD_UPDATE         = 8;

string ip_addresses[nP]                 = {"127.0.0.1", "127.0.0.1"};

// Input positions for swapping circuit
const int LEVEL_SIZE            = sizeof(int);
const int EMPTY_SIZE            = sizeof(uint8_t);

// Input positions for prepare deepest circuit
const int GOAL_OFFSET           = 0;
const int SOURCE_OFFSET         = GOAL_OFFSET + LEVEL_SIZE;
const int CURR_LEVEL_OFFSET     = SOURCE_OFFSET + LEVEL_SIZE;
const int DEEPEST_LEVEL_OFFSET  = CURR_LEVEL_OFFSET + LEVEL_SIZE;
const int DEEPEST_OFFSET        = DEEPEST_LEVEL_OFFSET + LEVEL_SIZE;
const int NEW_SOURCE_OFFSET     = DEEPEST_OFFSET + LEVEL_SIZE;

// Input positions for prepare target circuit
const int TARGET_SRC_OFFSET     = 0;
const int TARGET_DEST_OFFSET    = LEVEL_SIZE;
const int TARGET_I_OFFSET       = TARGET_DEST_OFFSET + LEVEL_SIZE;
const int TARGET_EMPTY_OFFSET   = TARGET_I_OFFSET + LEVEL_SIZE;
const int TARGET_DEEPEST_OFFSET = TARGET_EMPTY_OFFSET + EMPTY_SIZE;
const int TARGET_TARGET_OFFSET  = TARGET_DEEPEST_OFFSET + LEVEL_SIZE;

// Input positions for eviction circuit
const int DEST_LEVEL_OFFSET     = 0;
const int CURRENT_BLOCK_OFFSET  = BLOCK_ID_SIZE + LEVEL_SIZE + LEVEL_SIZE;

// Data for retrieval/update stash/combine
block (*rGTM[MAX_NUM_CIRCUITS_RETRIEVAL])[4][nP+1];
block (*rGTK[MAX_NUM_CIRCUITS_RETRIEVAL])[4][nP+1];
bool  (*rGTv[MAX_NUM_CIRCUITS_RETRIEVAL])[4];
block (*rGT[MAX_NUM_CIRCUITS_RETRIEVAL])[nP+1][4][nP+1];
block *rlabels[MAX_NUM_CIRCUITS_RETRIEVAL];
bool  *rvalue[MAX_NUM_CIRCUITS_RETRIEVAL];
block *(reval_labels[MAX_NUM_CIRCUITS_RETRIEVAL])[nP+1];
int   rvalue_temp;

// Data for eviction 
block (*eGTM[MAX_NUM_CIRCUITS_RETRIEVAL])[4][nP+1];
block (*eGTK[MAX_NUM_CIRCUITS_RETRIEVAL])[4][nP+1];
bool  (*eGTv[MAX_NUM_CIRCUITS_RETRIEVAL])[4];
block (*eGT[MAX_NUM_CIRCUITS_RETRIEVAL])[nP+1][4][nP+1];
block *elabels[MAX_NUM_CIRCUITS_RETRIEVAL];
bool  *evalue[MAX_NUM_CIRCUITS_RETRIEVAL];
block *(eeval_labels[MAX_NUM_CIRCUITS_RETRIEVAL])[nP+1];
int   evalue_temp;

mutex mtx_write_data;

class OMAT
{
    public:
        ORAM                  *oram;
        
        Stash                 *stash_word;
        Stash                 *stash_file;

        int                   stash_word_size;
        int                   stash_file_size;

        NetIOMP<nP>           *(rios[MAX_NUM_CIRCUITS_RETRIEVAL])[2];
        NetIOMP<nP>           *(eios[MAX_NUM_CIRCUITS_RETRIEVAL])[2];

        ThreadPool            *rpool[MAX_NUM_CIRCUITS_RETRIEVAL];
        ThreadPool            *epool[MAX_NUM_CIRCUITS_EVICTION];

        NetIOMP<nP>           *ios[2];
        ThreadPool            *pool;

        zmq::context_t        **context_req;
        zmq::socket_t         **socket_req;
        
        zmq::context_t        **context_rep;
        zmq::socket_t         **socket_rep;

        // These circuits are used for processing keywords (rows)
        BristolFormat         *cf_retrieval_word;
        BristolFormat         *cf_update_word;
        BristolFormat         *cf_eviction_word;
        BristolFormat         *cf_deepest_word;
        BristolFormat         *cf_prepare_deepest_word;
        BristolFormat         *cf_prepare_target_word;

        CMPC<nP>*             mpc_retrieval_word[MAX_NUM_CIRCUITS_RETRIEVAL]; 
        CMPC<nP>*             mpc_update_word[MAX_NUM_CIRCUITS_RETRIEVAL]; 
        CMPC<nP>*             mpc_eviction_word[MAX_NUM_CIRCUITS_EVICTION]; 
        CMPC<nP>*             mpc_prepare_deepest_word; 
        CMPC<nP>*             mpc_prepare_target_word; 
        CMPC<nP>*             mpc_deepest_word; 

        // These circuits are used for processing files (columns)
        BristolFormat         *cf_retrieval_file;
        BristolFormat         *cf_update_file;
        BristolFormat         *cf_eviction_file;
        BristolFormat         *cf_deepest_file;
        BristolFormat         *cf_prepare_deepest_file;
        BristolFormat         *cf_prepare_target_file;

        CMPC<nP>*             mpc_retrieval_file; 
        CMPC<nP>*             mpc_update_file; 
        CMPC<nP>*             mpc_eviction_file; 
        CMPC<nP>*             mpc_prepare_deepest_file; 
        CMPC<nP>*             mpc_prepare_target_file; 
        CMPC<nP>*             mpc_deepest_file; 

        // Additional circuits
        BristolFormat         *cf_and_2blocks;
        BristolFormat         *cf_and_meta_data;
        CMPC<nP>*             mpc_and_2blocks[MAX_NUM_CIRCUITS_RETRIEVAL];
        CMPC<nP>*             mpc_and_meta_data[MAX_NUM_CIRCUITS_RETRIEVAL];
        
        int                   party;
        int                   time_step_word;
        int                   time_step_file;

        int                   *deepest_blocks;
        int                   *deepest_levels;
        int                   *target;
        int                   *full_path;
        int                   *rfull_path;
        int                   *efull_path;
        int                   *deepest;
        char                  *available;

        int                   num_ands_retrieval_word;
        int                   num_ands_update_word;
        int                   num_ands_eviction_word;

        int                   num_ands_retrieval_file;
        int                   num_ands_update_file;
        int                   num_ands_eviction_file;
        
        int                   num_ands_2blocks;
        int                   num_ands_meta_data;

        double                online_time;
        double                offline_time;

        OMAT(int party, int port);
        ~OMAT();

        void initialize(int word_block_size, int file_block_size);
#ifdef RUN_SEARCH
        // These functions used for retrieving based on keywords (rows)
        void retrieve_block_word(uint8_t **in, uint8_t **out, int path_id, int block_id);
        void update_stash_word(uint8_t **in, uint8_t **out);
        void prepare_deepest_word(uint8_t **in, uint8_t **out);
        void prepare_target_word(uint8_t **in, uint8_t **out);
        void get_deepest_blocks_word(uint8_t **in, uint8_t **out, const int &eviction_path);
        void eviction_word(uint8_t **in, uint8_t **out);

        // Function support to combine results of Bloom Filter
        void and_2blocks(uint8_t **in, uint8_t **out, uint8_t **temp);
        void and_meta_data(uint8_t **in, uint8_t **out, uint8_t **temp);
#else 
        // These functions used for retrieving based on files (columns) 
        void retrieve_block_file(uint8_t *in, uint8_t *out, int path_id, int block_id);
        void update_stash_file(uint8_t *in, uint8_t *out);
        void prepare_deepest_file(uint8_t *in, uint8_t *out);
        void prepare_target_file(uint8_t *in, uint8_t *out);
        void get_deepest_blocks_file(uint8_t *in, uint8_t *out, const int &eviction_path);
        void eviction_file(uint8_t *in, uint8_t *out);
#endif 

        int  reverse_order_word(const int &t);
        int  reverse_order_file(const int &t);

        void increase_time_step_word();
        void increase_time_step_file();

        void reset_time();
};

OMAT::OMAT(int party, int port)
{
    this->party = party;
    oram        = new ORAM();
    
    // Initialize timestep for eviction
    time_step_word = 0;
    time_step_file = 0;

#ifndef RUN_SEARCH 
    // ios for circuits
    ios[0]     = new NetIOMP<nP>(party, port);
    ios[1]     = new NetIOMP<nP>(party, port+2*(nP+1)*(nP+1)+1);

    ios[0]->flush();
    ios[1]->flush();

    // Initialize MPC for swapping    
    pool = new ThreadPool(MAX_NUM_THREAD_UPDATE);

#else 
    // ios for circuits
    for(int i = 0; i < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++i)
    {
        rios[i][0] = new NetIOMP<nP>(party, port);
        rios[i][1] = new NetIOMP<nP>(party, port + 2*(nP+1)*(nP+1)+1);
        
        rios[i][0]->flush();
        rios[i][1]->flush();
        
        // Initialize thread pool 
        rpool[i] = new ThreadPool(nP+1);
    }

    int eport = port + 4096;

    for(int i = 0; i < MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_EVICTION; ++i)
    {
        eios[i][0] = new NetIOMP<nP>(party, eport);
        eios[i][1] = new NetIOMP<nP>(party, eport + 2*(nP+1)*(nP+1)+1);
        
        eios[i][0]->flush();    
        eios[i][1]->flush();

        // Initialize thread pool 
        epool[i] = new ThreadPool(nP+1);
    }
#endif 

}

void OMAT::initialize(int word_block_size, int file_block_size)
{
    // Init ORAM data     
    oram->word_block_size = word_block_size;
    oram->file_block_size = file_block_size;
    oram->init_tree();
    
    // Init stash
    this->stash_word_size = oram->depth_word;
    this->stash_file_size = oram->depth_file;
    stash_word = new Stash(stash_word_size, word_block_size);
    stash_file = new Stash(stash_file_size, file_block_size>>3);

 #ifdef RUN_SEARCH   
    // Circuits
    size_t data_size_per_circuit = (word_block_size<<3)/MAX_NUM_CIRCUITS_RETRIEVAL;

    cout << "Building prepare deepest word circuit..." << endl;
    build_prepare_deepest_circuit(oram->depth_word+1, "prepare_deepest_word.txt");
    cf_prepare_deepest_word  = new BristolFormat("prepare_deepest_word.txt");
    
    cout << "Building prepare target word circuit..." << endl;
    build_prepare_target_circuit(oram->depth_word+1, "prepare_target_word.txt");
    cf_prepare_target_word  = new BristolFormat("prepare_target_word.txt");
    
    cout << "Building deepest word circuit..." << endl;
    build_deepest_circuit(stash_word_size, oram->depth_word, BUCKET_SIZE, "deepest_word.txt");
    cf_deepest_word  = new BristolFormat("deepest_word.txt");
    
    cout << "Building combining circuit..." << endl;
    build_and_2blocks("and_2blocks.txt", data_size_per_circuit, THREAD_PER_CIRCUIT);
    cf_and_2blocks  = new BristolFormat("and_2blocks.txt");
    
    cout << "Building and meta data circuit..." << endl;
    build_and_meta_data("and_meta_data.txt", data_size_per_circuit, THREAD_PER_CIRCUIT);
    cf_and_meta_data  = new BristolFormat("and_meta_data.txt");

    cout << "Building retrieval word circuit..." << endl;
    uint64_t max_num_wires_retrieval = build_swapping_circuit(stash_word_size+oram->depth_word*BUCKET_SIZE, data_size_per_circuit, "retrieval_word.txt", THREAD_PER_CIRCUIT);
    cf_retrieval_word  = new BristolFormat("retrieval_word.txt");
    
    cout << "Building update stash word circuit..." << endl;
    build_swapping_circuit(stash_word_size, data_size_per_circuit, "update_word.txt", THREAD_PER_CIRCUIT);
    cf_update_word  = new BristolFormat("update_word.txt");
    
    cout << "Building eviction word circuit..." << endl;
    uint64_t max_num_wires_eviction = build_eviction_circuit(stash_word_size+oram->depth_word*BUCKET_SIZE, data_size_per_circuit, "eviction_word.txt", THREAD_PER_CIRCUIT);
    cf_eviction_word  = new BristolFormat("eviction_word.txt");
    
    uint64_t max_num_and_gates_retrieval = 0;
    for(int i = 0; i < cf_retrieval_word->num_gate; ++i) 
    {
        if (cf_retrieval_word->gates[4*i+3] == AND_GATE)
            ++max_num_and_gates_retrieval;
    }

    uint64_t max_num_and_gates_eviction = 0;
    for(int i = 0; i < cf_eviction_word->num_gate; ++i) 
    {
        if (cf_eviction_word->gates[4*i+3] == AND_GATE)
            ++max_num_and_gates_eviction;
    }
    
    // Retrieval
    if(party == 1) 
    {
        for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t) 
            rGTM[t] = new block[max_num_and_gates_retrieval][4][nP+1];

        for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t) 
            rGTK[t] = new block[max_num_and_gates_retrieval][4][nP+1];

        for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t) 
            rGTv[t] = new bool[max_num_and_gates_retrieval][4];

        for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t) 
            rGT[t] = new block[max_num_and_gates_retrieval][nP+1][4][nP+1];
    }

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        rlabels[t] = new block[max_num_wires_retrieval];
    
    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        rvalue[t] = new bool[max_num_wires_retrieval];

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        for(int i = 1; i <= nP; ++i)
            reval_labels[t][i] = new block[max_num_wires_retrieval];
    
    // Eviction
    if(party == 1) 
    {
        for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_RETRIEVAL; ++t) 
            eGTM[t] = new block[max_num_and_gates_eviction][4][nP+1];

        for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_RETRIEVAL; ++t) 
            eGTK[t] = new block[max_num_and_gates_eviction][4][nP+1];

        for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_RETRIEVAL; ++t) 
            eGTv[t] = new bool[max_num_and_gates_eviction][4];

        for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_RETRIEVAL; ++t) 
            eGT[t] = new block[max_num_and_gates_eviction][nP+1][4][nP+1];
    }

    for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        elabels[t] = new block[max_num_wires_eviction];
    
    for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        evalue[t] = new bool[max_num_wires_eviction];

    for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        for(int i = 1; i <= nP; ++i)
            eeval_labels[t][i] = new block[max_num_wires_eviction];
    
    deepest_blocks = new int[oram->depth_word + 1];
    deepest_levels = new int[oram->depth_word + 1];
    target         = new int[oram->depth_word + 1];
    deepest        = new int[oram->depth_word + 1];
    rfull_path     = new int[oram->depth_word];
    efull_path     = new int[oram->depth_word];
    available      = new char[oram->depth_word + 1]; 

#else 
    cout << "Preprocessing prepare deepest file circuit..." << endl;
    build_prepare_deepest_circuit(oram->depth_file+1, "prepare_deepest_file.txt");
    cf_prepare_deepest_file  = new BristolFormat("prepare_deepest_file.txt");
    mpc_prepare_deepest_file = new CMPC<nP>(ios, pool, party, cf_prepare_deepest_file);
    mpc_prepare_deepest_file->function_independent();
    ios[0]->flush();
    ios[1]->flush();
    mpc_prepare_deepest_file->function_dependent();
    ios[0]->flush();
    ios[1]->flush(); 

    cout << "Preprocessing prepare target file circuit..." << endl;
    build_prepare_target_circuit(oram->depth_file+1, "prepare_target_file.txt");
    cf_prepare_target_file  = new BristolFormat("prepare_target_file.txt");
    mpc_prepare_target_file = new CMPC<nP>(ios, pool, party, cf_prepare_target_file);
    mpc_prepare_target_file->function_independent();
    ios[0]->flush();
    ios[1]->flush();
    mpc_prepare_target_file->function_dependent();
    ios[0]->flush();
    ios[1]->flush(); 
    
    cout << "Preprocessing deepest file circuit..." << endl;
    build_deepest_circuit(stash_file_size, oram->depth_file, BUCKET_SIZE, "deepest_file.txt");
    cf_deepest_file  = new BristolFormat("deepest_file.txt");
    mpc_deepest_file = new CMPC<nP>(ios, pool, party, cf_deepest_file);
    mpc_deepest_file->function_independent();
    ios[0]->flush();
    ios[1]->flush();
    mpc_deepest_file->function_dependent();
    ios[0]->flush();
    ios[1]->flush(); 

    cout << "Preprocessing retrieval file circuit..." << endl;
    build_swapping_circuit(stash_file_size+oram->depth_file*BUCKET_SIZE, file_block_size, "retrieval_file.txt", MAX_NUM_THREAD_UPDATE);
    cf_retrieval_file  = new BristolFormat("retrieval_file.txt");
    mpc_retrieval_file = new CMPC<nP>(ios, pool, party, cf_retrieval_file, &num_ands_retrieval_file);
    mpc_retrieval_file->function_independent();
    ios[0]->flush();
    ios[1]->flush();
    mpc_retrieval_file->function_dependent();
    ios[0]->flush();
    ios[1]->flush(); 
    
    cout << "Preprocessing update stash file circuit..." << endl;
    build_swapping_circuit(stash_file_size, file_block_size, "update_file.txt", MAX_NUM_THREAD_UPDATE);
    cf_update_file  = new BristolFormat("update_file.txt");
    mpc_update_file = new CMPC<nP>(ios, pool, party, cf_update_file, &num_ands_update_file);
    mpc_update_file->function_independent();
    ios[0]->flush();
    ios[1]->flush();
    mpc_update_file->function_dependent();
    ios[0]->flush();
    ios[1]->flush(); 
    
    cout << "Preprocessing eviction file circuit..." << endl;
    build_eviction_circuit(stash_file_size+oram->depth_file*BUCKET_SIZE, file_block_size, "eviction_file.txt", MAX_NUM_THREAD_UPDATE);
    cf_eviction_file  = new BristolFormat("eviction_file.txt");
    mpc_eviction_file = new CMPC<nP>(ios, pool, party, cf_eviction_file, &num_ands_eviction_file);
    mpc_eviction_file->function_independent();
    ios[0]->flush();
    ios[1]->flush();
    mpc_eviction_file->function_dependent();
    ios[0]->flush();
    ios[1]->flush(); 
    
    deepest_blocks = new int[oram->depth_file + 1];
    deepest_levels = new int[oram->depth_file + 1];
    target         = new int[oram->depth_file + 1];
    deepest        = new int[oram->depth_file + 1];
    full_path      = new int[oram->depth_file];
    available      = new char[oram->depth_file + 1];
#endif 

}

OMAT::~OMAT()
{
    delete oram;
    delete stash_word;
    
    delete [] deepest_blocks;
    delete [] deepest_levels;
    delete [] target;
    delete [] rfull_path;
    delete [] efull_path;
    delete [] deepest;
    delete [] available;
    delete [] rpool;
    delete [] rios;
    delete [] epool;
    delete [] eios;

    delete cf_retrieval_word;
    delete cf_eviction_word;
    delete cf_deepest_word;
    delete cf_prepare_deepest_word;
    delete cf_prepare_target_word;

    delete mpc_retrieval_word; 
    delete mpc_update_word; 
    delete mpc_eviction_word; 
    delete mpc_prepare_deepest_word; 
    delete mpc_prepare_target_word; 
    delete mpc_deepest_word; 

    delete cf_retrieval_file;
    delete cf_eviction_file;
    delete cf_deepest_file;
    delete cf_prepare_deepest_file;
    delete cf_prepare_target_file;

    delete mpc_retrieval_file; 
    delete mpc_update_file; 
    delete mpc_eviction_file; 
    delete mpc_prepare_deepest_file; 
    delete mpc_prepare_target_file; 
    delete mpc_deepest_file; 

    delete cf_and_2blocks;
}

void OMAT::reset_time()
{
    offline_time = 0;
    online_time  = 0;    
}

#ifdef RUN_SEARCH
void OMAT::retrieve_block_word(uint8_t **in, uint8_t **out, int path_id, int block_id)
{
    auto start = clock_start();

    cout << "Preprocessing retrieval word circuit..." << endl;
    
    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
    {
        // int p = t/CIRCUITS_PER_THREAD_RETRIEVAL;
        mpc_retrieval_word[t] = new CMPC<nP>(rios[t], rpool[t], party, cf_retrieval_word, &num_ands_retrieval_word, 
                                             rGTM[t], rGTK[t], rGTv[t], rGT[t], rlabels[t], rvalue[t], reval_labels[t]); 
        mpc_retrieval_word[t]->function_independent();
        rios[t][0]->flush();
        rios[t][1]->flush();
        mpc_retrieval_word[t]->function_dependent();
        rios[t][0]->flush();
        rios[t][1]->flush(); 
    }

    offline_time += time_from(start);

    vector<future<void>> res;
    ThreadPool pool(MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL);

    cout << "Retrieving...";

    start = clock_start();

    oram->get_path_word(path_id, rfull_path);
    
    size_t data_size_per_circuit = oram->word_block_size/MAX_NUM_CIRCUITS_RETRIEVAL;

    mtx_write_data.lock();

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
    {
        memcpy(in[t], &block_id, BLOCK_ID_SIZE);
        memset(in[t]+BLOCK_ID_SIZE, 0, data_size_per_circuit + META_DATA_SIZE);
    }
    
    size_t offset;

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
    {
        offset = BLOCK_ID_SIZE + data_size_per_circuit + META_DATA_SIZE;
        for(int i = 0; i < stash_word_size; ++i) 
        {
            memcpy(in[t] + offset, stash_word->meta_data + i, META_DATA_SIZE); 
            offset += META_DATA_SIZE;
            memcpy(in[t] + offset, stash_word->data[i] + t * data_size_per_circuit, data_size_per_circuit); 
            offset += data_size_per_circuit;
        }
    }

    size_t prev_offset = offset;

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
    {
        offset = prev_offset;
        for(int i = 0; i < oram->depth_word; ++i)
            for(int j = 0; j < BUCKET_SIZE; ++j)
            {
                memcpy(in[t] + offset, oram->meta_data_word[rfull_path[i]] + j, META_DATA_SIZE); 
                offset += META_DATA_SIZE;
                memcpy(in[t] + offset, oram->data[rfull_path[i]][j] + t * data_size_per_circuit, data_size_per_circuit);
                offset += data_size_per_circuit;
            }
    }

    mtx_write_data.unlock();

	for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; t += CIRCUITS_PER_THREAD_RETRIEVAL) 
    {
        res.push_back(pool.enqueue([this, in, out, t]() 
		{
            auto start = clock_start();
            for(int v = t; v < t + CIRCUITS_PER_THREAD_RETRIEVAL; ++v) 
            {
                mpc_retrieval_word[v/CIRCUITS_PER_THREAD_RETRIEVAL]->online(in[v], out[v], num_ands_retrieval_word);
                rios[v/CIRCUITS_PER_THREAD_RETRIEVAL][0]->flush();
                rios[v/CIRCUITS_PER_THREAD_RETRIEVAL][1]->flush();
            }
            // cout << "Time thread #" << t/CIRCUITS_PER_THREAD_RETRIEVAL << ": " << time_from(start) << endl; 
        }));
    }
    
    joinNclean(res);

    mtx_write_data.lock();

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
    {
        offset = data_size_per_circuit + META_DATA_SIZE;
        for(int i = 0; i < stash_word_size; ++i)
        {
            memcpy(stash_word->meta_data + i, out[t]+offset, META_DATA_SIZE);
            memcpy(stash_word->data[i] + t * data_size_per_circuit, out[t]+offset+META_DATA_SIZE, data_size_per_circuit);
            offset += data_size_per_circuit + META_DATA_SIZE;
        }
    }

    prev_offset = offset;

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
    {   
        offset = prev_offset;
        for(int i = 0; i < oram->depth_word; ++i)
        {
            for(int j = 0; j < BUCKET_SIZE; ++j)
            {
                memcpy(oram->meta_data_word[rfull_path[i]] + j, out[t]+offset, META_DATA_SIZE);
                memcpy(oram->data[rfull_path[i]][j] + t * data_size_per_circuit, out[t]+offset+META_DATA_SIZE, data_size_per_circuit);
                offset += data_size_per_circuit + META_DATA_SIZE;
            }
        }
    }

    mtx_write_data.unlock();

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        delete mpc_retrieval_word[t];

    online_time += time_from(start);

    cout << "Fnished" << endl; 
}

void OMAT::update_stash_word(uint8_t **in, uint8_t **out)
{
    auto start = clock_start();

    cout << "Preprocessing update stash word circuit..." << endl;
    
    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
    {
        // int p = t/CIRCUITS_PER_THREAD_RETRIEVAL;
        mpc_update_word[t] = new CMPC<nP>(rios[t],rpool[t], party, cf_update_word, &num_ands_update_word, 
                                          rGTM[t], rGTK[t], rGTv[t], rGT[t], rlabels[t], rvalue[t], reval_labels[t]);    
    
        mpc_update_word[t]->function_independent();
        rios[t][0]->flush();
        rios[t][1]->flush();
        mpc_update_word[t]->function_dependent();
        rios[t][0]->flush();
        rios[t][1]->flush(); 
    }

    offline_time += time_from(start);

    vector<future<void>> res;
    ThreadPool pool(MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL);

    cout << "Updating stash...";

    start = clock_start();

    size_t data_size_per_circuit = oram->word_block_size/MAX_NUM_CIRCUITS_RETRIEVAL;

    mtx_write_data.lock();

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
    {
        memset(in[t], 0, BLOCK_ID_SIZE);
        memcpy(in[t] + BLOCK_ID_SIZE, out[t], data_size_per_circuit + META_DATA_SIZE);

        size_t offset = BLOCK_ID_SIZE + data_size_per_circuit + META_DATA_SIZE;

        for(int i = 0; i < stash_word_size; ++i)
        {
            memcpy(in[t] + offset, stash_word->meta_data + i, META_DATA_SIZE); 
            offset += META_DATA_SIZE;
            memcpy(in[t] + offset, stash_word->data[i] + t * data_size_per_circuit, data_size_per_circuit); 
            offset += data_size_per_circuit;
        }    
    }

    mtx_write_data.unlock();

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; t += CIRCUITS_PER_THREAD_RETRIEVAL) 
    {
        res.push_back(pool.enqueue([this, in, out, t]() 
		{
            for(int v = t; v < t + CIRCUITS_PER_THREAD_RETRIEVAL; ++v) 
            {
                mpc_update_word[v/CIRCUITS_PER_THREAD_RETRIEVAL]->online(in[v], out[v], num_ands_update_word);
                rios[v/CIRCUITS_PER_THREAD_RETRIEVAL][0]->flush();
                rios[v/CIRCUITS_PER_THREAD_RETRIEVAL][1]->flush();
            }
        }));
    }
    
    joinNclean(res);

    //cout << "ONLINE UPDATE STASH: \t" << party << "\t" << t2 << " \n" << flush;
    
    mtx_write_data.lock();

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
    {
        size_t offset = data_size_per_circuit + META_DATA_SIZE;
        for(int i = 0; i < stash_word_size; ++i)
        {
            memcpy(stash_word->meta_data+i, out[t]+offset, META_DATA_SIZE);
            memcpy(stash_word->data[i] + t * data_size_per_circuit, out[t] + offset + META_DATA_SIZE, data_size_per_circuit);
            offset += data_size_per_circuit + META_DATA_SIZE;
        }
    }

    mtx_write_data.unlock();

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        delete mpc_update_word[t];
    
    online_time += time_from(start);

    cout << "Fnished" << endl;
}

void OMAT::prepare_deepest_word(uint8_t **in, uint8_t **out)
{
    // cout << "Preprocessing prepare deepest word circuit..." << endl;
    mpc_prepare_deepest_word = new CMPC<nP>(eios[0], epool[0], party, cf_prepare_deepest_word, &evalue_temp, eGTM[0], 
                                            eGTK[0], eGTv[0], eGT[0], elabels[0], evalue[0], eeval_labels[0]);    
    mpc_prepare_deepest_word->function_independent();
    eios[0][0]->flush();
    eios[0][1]->flush();
    mpc_prepare_deepest_word->function_dependent();
    eios[0][0]->flush();
    eios[0][1]->flush(); 

    if(party == 1)
        memset(in[0]+SOURCE_OFFSET, 255, LEVEL_SIZE);
    else
        memset(in[0]+SOURCE_OFFSET, 0, LEVEL_SIZE);
        
    memset(in[0]+GOAL_OFFSET, 0, LEVEL_SIZE);
    
    uint8_t *in_ptr = in[0];
    int curr_level = (1<<(31-oram->depth_word))-1;
    
    for(int i = 0; i < oram->depth_word + 1; ++i)
    {
        curr_level = (curr_level << 1) | 1;
        if(party == 1)
        {
            memcpy(in_ptr+CURR_LEVEL_OFFSET, &curr_level, LEVEL_SIZE);
            memset(in_ptr+DEEPEST_OFFSET, 255, LEVEL_SIZE);
            memcpy(in_ptr+NEW_SOURCE_OFFSET, &i, LEVEL_SIZE);
        }
        else
        {   
            memset(in_ptr+CURR_LEVEL_OFFSET, 0, LEVEL_SIZE);
            memset(in_ptr+DEEPEST_OFFSET, 0, LEVEL_SIZE);
            memset(in_ptr+NEW_SOURCE_OFFSET, 0, LEVEL_SIZE);
        } 
        memcpy(in_ptr+DEEPEST_LEVEL_OFFSET, &deepest_levels[i], LEVEL_SIZE);
        in_ptr += (LEVEL_SIZE << 2);
    }
    
    // auto start = clock_start();

    mpc_prepare_deepest_word->online(in[0], out[0]);
    eios[0][0]->flush();
    eios[0][1]->flush();

    // online_time += time_from(start);

    int offset = LEVEL_SIZE;
    for(int i = 1; i < oram->depth_word + 1; ++i)
    {
        deepest[i] = *(int*)(out[0]+offset);
        offset    += LEVEL_SIZE;
    }
    
    delete mpc_prepare_deepest_word;

    // cout << "Deepest: ";
    // for(int i = 1; i <= oram->depth_word; ++i)
    //     cout << deepest[i] << " ";
    // cout << endl;
}

void OMAT::prepare_target_word(uint8_t **in, uint8_t **out)
{
    // cout << "Preprocessing prepare target word circuit..." << endl;
    mpc_prepare_target_word = new CMPC<nP>(eios[0], epool[0], party, cf_prepare_target_word, &evalue_temp, eGTM[0], 
                                           eGTK[0], eGTv[0], eGT[0], elabels[0], evalue[0], eeval_labels[0]);    
    mpc_prepare_target_word->function_independent();
    eios[0][0]->flush();
    eios[0][1]->flush();
    mpc_prepare_target_word->function_dependent();
    eios[0][0]->flush();
    eios[0][1]->flush(); 

    if(party == 1)
    {
        memset(in[0]+TARGET_SRC_OFFSET, 255, LEVEL_SIZE);
        memset(in[0]+TARGET_DEST_OFFSET, 255, LEVEL_SIZE);
    }
    else 
    {
        memset(in[0]+TARGET_SRC_OFFSET, 0, LEVEL_SIZE);
        memset(in[0]+TARGET_DEST_OFFSET, 0, LEVEL_SIZE);
    }
    
    uint8_t *in_ptr = in[0];

    for(int i = oram->depth_word; i >= 0; --i)
    {
        if(party == 1)
        {
            memcpy(in_ptr + TARGET_I_OFFSET, &i, LEVEL_SIZE);
            memset(in_ptr + TARGET_TARGET_OFFSET, 255, LEVEL_SIZE);
        }
        else
        {
            memset(in_ptr + TARGET_I_OFFSET, 0, LEVEL_SIZE);
            memset(in_ptr + TARGET_TARGET_OFFSET, 0, LEVEL_SIZE);
        }
        memcpy(in_ptr+TARGET_EMPTY_OFFSET, &available[i], EMPTY_SIZE);
        memcpy(in_ptr+TARGET_DEEPEST_OFFSET, &deepest[i], LEVEL_SIZE);
        in_ptr += TARGET_DEEPEST_OFFSET;
    }
    
    // auto start = clock_start();

    mpc_prepare_target_word->online(in[0], out[0]);
    eios[0][0]->flush();
    eios[0][1]->flush();

    // online_time += time_from(start);

    int offset = 0;
    for(int i = oram->depth_word; i >= 0; --i)
    {
        target[i] = *(int*)(out[0] + offset);
        offset      += LEVEL_SIZE;
    }

    delete mpc_prepare_target_word;

    //cout << "Target: ";
    //for(int i = 0; i <= oram->depth_word; ++i)
    //    cout << target[i] << " ";
    //cout << endl;
}

void OMAT::get_deepest_blocks_word(uint8_t **in, uint8_t **out, const int &eviction_path)
{
    // cout << "Preprocessing deepest word circuit..." << endl;
    mpc_deepest_word = new CMPC<nP>(eios[0], epool[0], party, cf_deepest_word, &evalue_temp, eGTM[0], 
                                    eGTK[0], eGTv[0], eGT[0], elabels[0], evalue[0], eeval_labels[0]);    
    mpc_deepest_word->function_independent();
    eios[0][0]->flush();
    eios[0][1]->flush();
    mpc_deepest_word->function_dependent();
    eios[0][0]->flush();
    eios[0][1]->flush(); 
    
    this->oram->get_path_word(eviction_path, efull_path);
    memset(in[0], 0, BLOCK_ID_SIZE+LEVEL_SIZE+EMPTY_SIZE);

    uint8_t *in_ptr = in[0]+BLOCK_ID_SIZE+LEVEL_SIZE+EMPTY_SIZE;
    for(int i = 0; i < stash_word_size; ++i) 
    {
        memcpy(in_ptr, &(stash_word->meta_data[i].path_id), PATH_ID_SIZE);
        memcpy(in_ptr+PATH_ID_SIZE, &(stash_word->meta_data[i].block_id), BLOCK_ID_SIZE);
        if(party == 1)
            memcpy(in_ptr+PATH_ID_SIZE+BLOCK_ID_SIZE, &eviction_path, PATH_ID_SIZE);
        else 
            memset(in_ptr+PATH_ID_SIZE+BLOCK_ID_SIZE, 0, PATH_ID_SIZE);
        in_ptr += PATH_ID_SIZE + BLOCK_ID_SIZE + PATH_ID_SIZE;
    }

    for(int i = 0; i < oram->depth_word; ++i)
    {
        memset(in_ptr, 0, BLOCK_ID_SIZE+LEVEL_SIZE+EMPTY_SIZE);
        in_ptr += BLOCK_ID_SIZE+LEVEL_SIZE+EMPTY_SIZE;
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            memcpy(in_ptr, &(oram->meta_data_word[efull_path[i]][j].path_id), PATH_ID_SIZE);
            memcpy(in_ptr+PATH_ID_SIZE, &(oram->meta_data_word[efull_path[i]][j].block_id), BLOCK_ID_SIZE);
            if(party == 1)
                memcpy(in_ptr+PATH_ID_SIZE+BLOCK_ID_SIZE, &eviction_path, PATH_ID_SIZE);
            else 
                memset(in_ptr+PATH_ID_SIZE+BLOCK_ID_SIZE, 0, PATH_ID_SIZE);
            in_ptr += PATH_ID_SIZE + BLOCK_ID_SIZE + PATH_ID_SIZE;
        }
    }

    // auto start = clock_start();

    mpc_deepest_word->online(in[0], out[0]);
    eios[0][0]->flush();
    eios[0][1]->flush();

    // online_time += time_from(start);

    int offset = 0;
    for(int i = 0; i <= oram->depth_word; ++i)
    {
        deepest_blocks[i] = *(int*)(out[0]+offset);
        deepest_levels[i] = *(int*)(out[0]+offset+BLOCK_ID_SIZE);
        available[i]      = *(char*)(out[0]+offset+BLOCK_ID_SIZE+LEVEL_SIZE);
        offset           += (LEVEL_SIZE+BLOCK_ID_SIZE+EMPTY_SIZE);
    }

    delete mpc_deepest_word;
}

void OMAT::eviction_word(uint8_t **in, uint8_t **out)
{
    // cout << "Preprocessing eviction word circuit..." << endl;

    for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_EVICTION; ++t)
    {
        // int p = t/CIRCUITS_PER_THREAD_EVICTION;
        mpc_eviction_word[t] = new CMPC<nP>(eios[t], epool[t], party, cf_eviction_word, &num_ands_eviction_word, 
                                            eGTM[t], eGTK[t], eGTv[t], eGT[t], elabels[t], evalue[t], eeval_labels[t]);
        mpc_eviction_word[t]->function_independent();
        eios[t][0]->flush();
        eios[t][1]->flush();
        mpc_eviction_word[t]->function_dependent();
        eios[t][0]->flush();
        eios[t][1]->flush(); 
    }
    
    //for(int i = 0; i <= oram->depth_word; ++i)
    //    cout << "Block ID: " << deepest_blocks[i] << " will be moved to the level: " << target[i] << endl;

    size_t data_size_per_circuit = oram->word_block_size/MAX_NUM_CIRCUITS_EVICTION;

    mtx_write_data.lock();

    for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION; ++t)
    {    
        memset(in[t]+DEST_LEVEL_OFFSET, 0, LEVEL_SIZE);
    
        uint8_t *in_ptr = in[t] + LEVEL_SIZE;
        for(int i = 0; i < stash_word_size; ++i) 
        {   
            memset(in_ptr, 0, LEVEL_SIZE);
            memcpy(in_ptr + LEVEL_SIZE, &deepest_blocks[0], BLOCK_ID_SIZE);
            memcpy(in_ptr + LEVEL_SIZE + BLOCK_ID_SIZE, &target[0], LEVEL_SIZE);
            memcpy(in_ptr + CURRENT_BLOCK_OFFSET, stash_word->meta_data + i, META_DATA_SIZE); 
            memcpy(in_ptr + CURRENT_BLOCK_OFFSET + META_DATA_SIZE, stash_word->data[i] + t * data_size_per_circuit, data_size_per_circuit); 
            in_ptr += (((LEVEL_SIZE<<1)+BLOCK_ID_SIZE)+data_size_per_circuit + META_DATA_SIZE);
        }

        for(int i = 1; i <= oram->depth_word; ++i)
        {
            for(int j = 0; j < BUCKET_SIZE; ++j)
            {
                if(party == 1)
                    memcpy(in_ptr, &i, LEVEL_SIZE);
                else 
                    memset(in_ptr, 0, LEVEL_SIZE);
                memcpy(in_ptr + LEVEL_SIZE, &deepest_blocks[i], BLOCK_ID_SIZE);
                memcpy(in_ptr + LEVEL_SIZE + BLOCK_ID_SIZE, &target[i], LEVEL_SIZE);
                memcpy(in_ptr + CURRENT_BLOCK_OFFSET, oram->meta_data_word[efull_path[i-1]] + j, META_DATA_SIZE); 
                memcpy(in_ptr + CURRENT_BLOCK_OFFSET + META_DATA_SIZE, oram->data[efull_path[i-1]][j] + t * data_size_per_circuit, data_size_per_circuit);
                in_ptr += (((LEVEL_SIZE<<1)+BLOCK_ID_SIZE)+data_size_per_circuit + META_DATA_SIZE);
            }
        }
        
        // Need to replace this by a random share of 0
        memset(in_ptr, 0, data_size_per_circuit + META_DATA_SIZE);
    }

    mtx_write_data.unlock();

    vector<future<void>> res;
    ThreadPool pool(MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_EVICTION);

    // auto start = clock_start();

    for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION; t += CIRCUITS_PER_THREAD_EVICTION) 
    {
        res.push_back(pool.enqueue([this, in, out, t]() 
		{
            for(int v = t; v < t + CIRCUITS_PER_THREAD_EVICTION; ++v) 
            {
                mpc_eviction_word[v/CIRCUITS_PER_THREAD_EVICTION]->online(in[v], out[v], num_ands_eviction_word);
                eios[v/CIRCUITS_PER_THREAD_EVICTION][0]->flush();
                eios[v/CIRCUITS_PER_THREAD_EVICTION][1]->flush();
            }
        }));
    }

    joinNclean(res);

    // online_time += time_from(start);

    mtx_write_data.lock();

    for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION; ++t)
    {
        size_t offset = 0;
        for(int i = 0; i < stash_word_size; ++i)
        {
            memcpy(stash_word->meta_data + i, out[t]+offset, META_DATA_SIZE);
            memcpy(stash_word->data[i] + t * data_size_per_circuit, out[t]+offset+META_DATA_SIZE, data_size_per_circuit);
            offset += data_size_per_circuit + META_DATA_SIZE;
        }

        for(int i = 0; i < oram->depth_word; ++i)
        {
            for(int j = 0; j < BUCKET_SIZE; ++j)
            {
                memcpy(oram->meta_data_word[efull_path[i]] + j, out[t]+offset, META_DATA_SIZE);
                memcpy(oram->data[efull_path[i]][j] + t * data_size_per_circuit, out[t]+offset+META_DATA_SIZE, data_size_per_circuit);
                offset += data_size_per_circuit + META_DATA_SIZE;
            }
        }
    }
    
    mtx_write_data.unlock();

    for(int t = 0; t < MAX_NUM_CIRCUITS_EVICTION/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        delete mpc_eviction_word[t];

    //cout << "ONLINE EVICTION: \t" << party << "\t" << t2 << " \n" << flush;
}
#endif 

void OMAT::increase_time_step_word()
{
    time_step_word = (time_step_word + 1) % (1 << (oram->depth_word-1));
}

void OMAT::increase_time_step_file()
{
    time_step_file = (time_step_file + 1) % (1 << (oram->depth_file-1));
}

int OMAT::reverse_order_word(const int &t)
{
    int no_of_bits = oram->depth_word - 1;
    int eviction_path = 0, temp;
    for(int i = 0; i < no_of_bits; ++i)
    {
        temp = t & (1 << i);
        if(temp)
            eviction_path |= (1 << ((no_of_bits - 1) - i));
    }
    return eviction_path;
}

int OMAT::reverse_order_file(const int &t)
{
    int no_of_bits = oram->depth_file - 1;
    int eviction_path = 0, temp;
    for(int i = 0; i < no_of_bits; ++i)
    {
        temp = t & (1 << i);
        if(temp)
            eviction_path |= (1 << ((no_of_bits - 1) - i));
    }
    return eviction_path;
}

#ifdef RUN_SEARCH
void OMAT::and_2blocks(uint8_t **in, uint8_t **out, uint8_t **temp)
{   
    auto start = clock_start();

    cout << "Preprocessing combining circuit..." << endl;

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
    {
        // int p = t/CIRCUITS_PER_THREAD_RETRIEVAL;
        mpc_and_2blocks[t] = new CMPC<nP>(rios[t], rpool[t], party, cf_and_2blocks, &num_ands_2blocks, 
                                          rGTM[t], rGTK[t], rGTv[t], rGT[t], rlabels[t], rvalue[t], reval_labels[t]);    
        mpc_and_2blocks[t]->function_independent();
        rios[t][0]->flush();
        rios[t][1]->flush();
        mpc_and_2blocks[t]->function_dependent();
        rios[t][0]->flush();
        rios[t][1]->flush();     
    }

    offline_time += time_from(start);

    vector<future<void>> res;
    ThreadPool pool(MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL);

    cout << "ANDing 2 blocks...";

    start = clock_start();

    size_t data_size_per_circuit = oram->word_block_size/MAX_NUM_CIRCUITS_RETRIEVAL;

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
        memcpy(temp[t] + data_size_per_circuit, out[t] + META_DATA_SIZE, data_size_per_circuit);

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; t += CIRCUITS_PER_THREAD_RETRIEVAL) 
    {
        res.push_back(pool.enqueue([this, temp, t]() 
		{
            for(int v = t; v < t + CIRCUITS_PER_THREAD_RETRIEVAL; ++v)
            {
                mpc_and_2blocks[v/CIRCUITS_PER_THREAD_RETRIEVAL]->online(temp[v], temp[v], num_ands_2blocks);
                rios[v/CIRCUITS_PER_THREAD_RETRIEVAL][0]->flush();
                rios[v/CIRCUITS_PER_THREAD_RETRIEVAL][1]->flush();
            }
        }));
    }
    
    joinNclean(res);

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        delete mpc_and_2blocks[t];

    online_time += time_from(start);

    cout << "Fnished" << endl;
}

void OMAT::and_meta_data(uint8_t **in, uint8_t **out, uint8_t **temp)
{   
    auto start = clock_start();

    cout << "Preprocessing and meta data circuit..." << endl; 

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
    {
        // int p = t/CIRCUITS_PER_THREAD_RETRIEVAL;
        mpc_and_meta_data[t] = new CMPC<nP>(rios[t], rpool[t], party, cf_and_meta_data, &num_ands_meta_data, 
                                            rGTM[t], rGTK[t], rGTv[t], rGT[t], rlabels[t], rvalue[t], reval_labels[t]);    
        mpc_and_meta_data[t]->function_independent();
        rios[t][0]->flush();
        rios[t][1]->flush();
        mpc_and_meta_data[t]->function_dependent();
        rios[t][0]->flush();
        rios[t][1]->flush(); 
    }
    
    offline_time += time_from(start);

    vector<future<void>> res;
    ThreadPool pool(MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL);

    cout << "ANDing meta-data...";
    
    start = clock_start();

    size_t data_size_per_circuit = oram->word_block_size/MAX_NUM_CIRCUITS_RETRIEVAL;
    size_t num_buckets_file_per_circuit = ((oram->word_block_size<<3)/BUCKET_SIZE)/MAX_NUM_CIRCUITS_RETRIEVAL;

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; ++t)
    {
        memcpy(in[t], temp[t], data_size_per_circuit);
        uint8_t *in_ptr = in[t] + data_size_per_circuit;

        size_t start_bucket = t * num_buckets_file_per_circuit;
        size_t end_bucket   = (t+1) * num_buckets_file_per_circuit;
        
        if(t == MAX_NUM_CIRCUITS_RETRIEVAL-1) 
            end_bucket -= 1;
        
        for(int i = start_bucket; i < end_bucket; ++i)
        {
            for(int j = 0; j < BUCKET_SIZE; ++j)
            {
                memcpy(in_ptr, &oram->meta_data_file[i][j].block_id, BLOCK_ID_SIZE);
                in_ptr += BLOCK_ID_SIZE;
            }
        }
    }
    
    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL; t += CIRCUITS_PER_THREAD_RETRIEVAL) 
    {
        res.push_back(pool.enqueue([this, in, out, t]() 
		{
            for(int v = t; v < t + CIRCUITS_PER_THREAD_RETRIEVAL; ++v)
            {
                mpc_and_meta_data[v/CIRCUITS_PER_THREAD_RETRIEVAL]->online(in[v], out[v], num_ands_meta_data);
                rios[v/CIRCUITS_PER_THREAD_RETRIEVAL][0]->flush();
                rios[v/CIRCUITS_PER_THREAD_RETRIEVAL][1]->flush();
            }
        }));
    }
    
    joinNclean(res);

    for(int t = 0; t < MAX_NUM_CIRCUITS_RETRIEVAL/CIRCUITS_PER_THREAD_RETRIEVAL; ++t)
        delete mpc_and_meta_data[t];

    online_time += time_from(start);

    cout << "Fnished" << endl;
}
#endif 

#ifndef RUN_SEARCH
void OMAT::retrieve_block_file(uint8_t *in, uint8_t *out, int path_id, int block_id)
{
    oram->get_path_file(path_id, full_path);
    int file_block_size_in_byte = oram->file_block_size>>3;
    
    memcpy(in, &block_id, BLOCK_ID_SIZE);
    memset(in+BLOCK_ID_SIZE, 0, file_block_size_in_byte + META_DATA_SIZE);

    size_t offset = BLOCK_ID_SIZE + file_block_size_in_byte + META_DATA_SIZE;

    for(int i = 0; i < stash_file_size; ++i) 
    {
        memcpy(in + offset, stash_file->meta_data + i, META_DATA_SIZE); 
        offset += META_DATA_SIZE;
        memcpy(in + offset, stash_file->data[i], file_block_size_in_byte); 
        offset += file_block_size_in_byte;
    }
    
    for(int i = 0; i < oram->depth_file; ++i)
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            memcpy(in + offset, oram->meta_data_file[full_path[i]] + j, META_DATA_SIZE); 
            offset           += META_DATA_SIZE;
            int position      = (full_path[i]<<1) + j;
            int byte_position = position >> 3;
            int bit_position  = (position % 8);
            uint8_t val       = 0;
            int count         = 0;
            // The following operation writes data row by row, which might cause memory access overhead
            for(int k = 0; k < oram->file_block_size/BUCKET_SIZE - 1; ++k)
            {
                for(int l = 0; l < BUCKET_SIZE; ++l)
                {
                    val <<= 1;
                    val |= (oram->data[k][l][byte_position] >> bit_position) & 0x1;
                }
                count += BUCKET_SIZE;
                if(count == 8) 
                {
                    memcpy(in + offset, &val, sizeof(uint8_t));
                    offset++;
                    count = 0;
                    val   = 0;
                }                
            } 
            offset++;
        }
    
    auto start = clock_start();
    
    mpc_retrieval_file->online(in, out, num_ands_retrieval_file);
    ios[0]->flush();
    ios[1]->flush();

    double t2 = time_from(start);

    cout << "ONLINE RETRIEVAL: \t" << party << "\t" << t2 << " \n" << flush;
    
    offset = file_block_size_in_byte + META_DATA_SIZE;
    for(int i = 0; i < stash_file_size; ++i)
    {
        memcpy(stash_file->meta_data + i, out+offset, META_DATA_SIZE);
        memcpy(stash_file->data[i], out+offset+META_DATA_SIZE, file_block_size_in_byte);
        offset += file_block_size_in_byte + META_DATA_SIZE;
    }
    
    for(int i = 0; i < oram->depth_file; ++i)
    {
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            memcpy(oram->meta_data_file[full_path[i]] + j, out+offset, META_DATA_SIZE);
            offset           += META_DATA_SIZE;
            int position      = (full_path[i]<<1) + j;
            int byte_position = position >> 3;
            int bit_position  = (position % 8);
            uint8_t val       = *(uint8_t*)(out+offset);
            int count         = 0;
            // The following operation writes data row by row, which might cause memory access overhead
            for(int k = 0; k < oram->file_block_size/BUCKET_SIZE - 1; ++k)
            {
                for(int l = 0; l < BUCKET_SIZE; ++l)
                {
                    oram->data[k][l][byte_position] &= ~(1<<bit_position);
                    oram->data[k][l][byte_position] |= ((val & 0x80) >> (7-bit_position));
                    val                            <<= 1;
                }
                count += BUCKET_SIZE;
                if(count == 8) 
                {
                    offset++;
                    count = 0;
                    val   = *(uint8_t*)(out+offset);
                }                
            }   
            offset++;
        }
    }
    cout << "File Block ID retrieved: " << *(int*)out << endl;
}

void OMAT::update_stash_file(uint8_t *in, uint8_t *out)
{
    memset(in, 0, BLOCK_ID_SIZE);
    int file_block_size_in_byte = oram->file_block_size>>3;

    memcpy(in+BLOCK_ID_SIZE, out, file_block_size_in_byte + META_DATA_SIZE);

    size_t offset = BLOCK_ID_SIZE + file_block_size_in_byte + META_DATA_SIZE;
    for(int i = 0; i < stash_file_size; ++i)
    {
        memcpy(in + offset, stash_file->meta_data + i, META_DATA_SIZE); 
        offset += META_DATA_SIZE;
        memcpy(in + offset, stash_file->data[i], file_block_size_in_byte); 
        offset += file_block_size_in_byte;
    }

    //auto start = clock_start();

    mpc_update_file->online(in, out, num_ands_update_file);
    ios[0]->flush();
    ios[1]->flush();

    //double t2 = time_from(start);
    //cout << "ONLINE UPDATE STASH: \t" << party << "\t" << t2 << " \n" << flush;
    
    offset = file_block_size_in_byte + META_DATA_SIZE;
    for(int i = 0; i < stash_file_size; ++i)
    {
        memcpy(stash_file->meta_data+i, out+offset, META_DATA_SIZE);
        memcpy(stash_file->data[i], out+offset+META_DATA_SIZE, file_block_size_in_byte);
        offset += file_block_size_in_byte + META_DATA_SIZE;
    }
}

void OMAT::prepare_deepest_file(uint8_t *in, uint8_t *out)
{
    if(party == 1)
        memset(in+SOURCE_OFFSET, 255, LEVEL_SIZE);
    else
        memset(in+SOURCE_OFFSET, 0, LEVEL_SIZE);
        
    memset(in+GOAL_OFFSET, 0, LEVEL_SIZE);
    
    uint8_t *in_ptr = in;
    int curr_level = (1<<(31-oram->depth_file))-1;
    
    for(int i = 0; i < oram->depth_file + 1; ++i)
    {
        curr_level = (curr_level << 1) | 1;
        if(party == 1)
        {
            memcpy(in_ptr+CURR_LEVEL_OFFSET, &curr_level, LEVEL_SIZE);
            memset(in_ptr+DEEPEST_OFFSET, 255, LEVEL_SIZE);
            memcpy(in_ptr+NEW_SOURCE_OFFSET, &i, LEVEL_SIZE);
        }
        else
        {   
            memset(in_ptr+CURR_LEVEL_OFFSET, 0, LEVEL_SIZE);
            memset(in_ptr+DEEPEST_OFFSET, 0, LEVEL_SIZE);
            memset(in_ptr+NEW_SOURCE_OFFSET, 0, LEVEL_SIZE);
        } 
        memcpy(in_ptr+DEEPEST_LEVEL_OFFSET, &deepest_levels[i], LEVEL_SIZE);
        in_ptr += (LEVEL_SIZE << 2);
    }
    
    mpc_prepare_deepest_file->online(in, out);
    ios[0]->flush();
    ios[1]->flush();

    int offset = LEVEL_SIZE;
    for(int i = 1; i < oram->depth_file + 1; ++i)
    {
        deepest[i] = *(int*)(out+offset);
        offset    += LEVEL_SIZE;
    }
    
    //cout << "Deepest: ";
    //for(int i = 1; i <= oram->depth_file; ++i)
    //    cout << deepest[i] << " ";
    //cout << endl;
}

void OMAT::prepare_target_file(uint8_t *in, uint8_t *out)
{
    if(party == 1)
    {
        memset(in+TARGET_SRC_OFFSET, 255, LEVEL_SIZE);
        memset(in+TARGET_DEST_OFFSET, 255, LEVEL_SIZE);
    }
    else 
    {
        memset(in+TARGET_SRC_OFFSET, 0, LEVEL_SIZE);
        memset(in+TARGET_DEST_OFFSET, 0, LEVEL_SIZE);
    }
    
    uint8_t *in_ptr = in;

    for(int i = oram->depth_file; i >= 0; --i)
    {
        if(party == 1)
        {
            memcpy(in_ptr + TARGET_I_OFFSET, &i, LEVEL_SIZE);
            memset(in_ptr + TARGET_TARGET_OFFSET, 255, LEVEL_SIZE);
        }
        else
        {
            memset(in_ptr + TARGET_I_OFFSET, 0, LEVEL_SIZE);
            memset(in_ptr + TARGET_TARGET_OFFSET, 0, LEVEL_SIZE);
        }
        memcpy(in_ptr+TARGET_EMPTY_OFFSET, &available[i], EMPTY_SIZE);
        memcpy(in_ptr+TARGET_DEEPEST_OFFSET, &deepest[i], LEVEL_SIZE);
        in_ptr += TARGET_DEEPEST_OFFSET;
    }
    
    mpc_prepare_target_file->online(in, out);
    ios[0]->flush();
    ios[1]->flush();

    int offset = 0;
    for(int i = oram->depth_file; i >= 0; --i)
    {
        target[i] = *(int*)(out + offset);
        offset   += LEVEL_SIZE;
    }

    //cout << "Target: ";
    //for(int i = 0; i <= oram->depth_file; ++i)
    //    cout << target[i] << " ";
    //cout << endl;
}

void OMAT::get_deepest_blocks_file(uint8_t *in, uint8_t *out, const int &eviction_path)
{
    this->oram->get_path_file(eviction_path, full_path);
    memset(in, 0, BLOCK_ID_SIZE+LEVEL_SIZE+EMPTY_SIZE);

    uint8_t *in_ptr = in+BLOCK_ID_SIZE+LEVEL_SIZE+EMPTY_SIZE;
    for(int i = 0; i < stash_file_size; ++i) 
    {
        memcpy(in_ptr, &(stash_file->meta_data[i].path_id), PATH_ID_SIZE);
        memcpy(in_ptr+PATH_ID_SIZE, &(stash_file->meta_data[i].block_id), BLOCK_ID_SIZE);
        if(party == 1)
            memcpy(in_ptr+PATH_ID_SIZE+BLOCK_ID_SIZE, &eviction_path, PATH_ID_SIZE);
        else 
            memset(in_ptr+PATH_ID_SIZE+BLOCK_ID_SIZE, 0, PATH_ID_SIZE);
        in_ptr += PATH_ID_SIZE + BLOCK_ID_SIZE + PATH_ID_SIZE;
    }

    for(int i = 0; i < oram->depth_file; ++i)
    {
        memset(in_ptr, 0, BLOCK_ID_SIZE+LEVEL_SIZE+EMPTY_SIZE);
        in_ptr += BLOCK_ID_SIZE+LEVEL_SIZE+EMPTY_SIZE;
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            memcpy(in_ptr, &(oram->meta_data_file[full_path[i]][j].path_id), PATH_ID_SIZE);
            memcpy(in_ptr+PATH_ID_SIZE, &(oram->meta_data_file[full_path[i]][j].block_id), BLOCK_ID_SIZE);
            if(party == 1)
                memcpy(in_ptr+PATH_ID_SIZE+BLOCK_ID_SIZE, &eviction_path, PATH_ID_SIZE);
            else 
                memset(in_ptr+PATH_ID_SIZE+BLOCK_ID_SIZE, 0, PATH_ID_SIZE);
            in_ptr += PATH_ID_SIZE + BLOCK_ID_SIZE + PATH_ID_SIZE;
        }
    }
    
    mpc_deepest_file->online(in, out);
    ios[0]->flush();
    ios[1]->flush();

    int offset = 0;
    for(int i = 0; i <= oram->depth_file; ++i)
    {
        deepest_blocks[i] = *(int*)(out+offset);
        deepest_levels[i] = *(int*)(out+offset+BLOCK_ID_SIZE);
        available[i]      = *(char*)(out+offset+BLOCK_ID_SIZE+LEVEL_SIZE);
        offset           += (LEVEL_SIZE+BLOCK_ID_SIZE+EMPTY_SIZE);
    }
}

void OMAT::eviction_file(uint8_t *in, uint8_t *out)
{
    int file_block_size_in_byte = oram->file_block_size>>3;
    memset(in+DEST_LEVEL_OFFSET, 0, LEVEL_SIZE);
    
    uint8_t *in_ptr = in + LEVEL_SIZE;
    for(int i = 0; i < stash_file_size; ++i) 
    {   
        memset(in_ptr, 0, LEVEL_SIZE);
        memcpy(in_ptr + LEVEL_SIZE, &deepest_blocks[0], BLOCK_ID_SIZE);
        memcpy(in_ptr + LEVEL_SIZE + BLOCK_ID_SIZE, &target[0], LEVEL_SIZE);
        memcpy(in_ptr + CURRENT_BLOCK_OFFSET, stash_file->meta_data + i, META_DATA_SIZE); 
        memcpy(in_ptr + CURRENT_BLOCK_OFFSET + META_DATA_SIZE, stash_file->data[i], file_block_size_in_byte); 
        in_ptr += (((LEVEL_SIZE<<1)+BLOCK_ID_SIZE)+file_block_size_in_byte + META_DATA_SIZE);
    }
    
    for(int i = 1; i <= oram->depth_file; ++i)
    {
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            if(party == 1)
                memcpy(in_ptr, &i, LEVEL_SIZE);
            else 
                memset(in_ptr, 0, LEVEL_SIZE);
            memcpy(in_ptr + LEVEL_SIZE, &deepest_blocks[i], BLOCK_ID_SIZE);
            memcpy(in_ptr + LEVEL_SIZE + BLOCK_ID_SIZE, &target[i], LEVEL_SIZE);
            memcpy(in_ptr + CURRENT_BLOCK_OFFSET, oram->meta_data_file[full_path[i-1]] + j, META_DATA_SIZE); 

            in_ptr           += CURRENT_BLOCK_OFFSET + META_DATA_SIZE;
            int position      = (full_path[i-1]<<1) + j;
            int byte_position = position >> 3;
            int bit_position  = (position % 8);
            uint8_t val       = 0;
            int count         = 0;
            // The following operation writes data row by row, which might cause memory access overhead
            for(int k = 0; k < oram->file_block_size/BUCKET_SIZE - 1; ++k)
            {
                for(int l = 0; l < BUCKET_SIZE; ++l)
                {
                    val <<= 1;
                    val |= (oram->data[k][l][byte_position] >> bit_position) & 0x1;
                }
                count += BUCKET_SIZE;
                if(count == 8) 
                {
                    memcpy(in_ptr, &val, sizeof(uint8_t));
                    in_ptr++;
                    count = 0;
                    val   = 0;
                }                
            } 
            in_ptr++;
        }
    }
    
    // Need to replace this by a random share of 0
    memset(in_ptr, 0, file_block_size_in_byte + META_DATA_SIZE);
    
    //auto start = clock_start();

    mpc_eviction_file->online(in, out, num_ands_eviction_file);
    ios[0]->flush();
    ios[1]->flush();    
    
    //double t2 = time_from(start);

    int offset = 0;
    for(int i = 0; i < stash_file_size; ++i)
    {
        memcpy(stash_file->meta_data + i, out+offset, META_DATA_SIZE);
        memcpy(stash_file->data[i], out+offset+META_DATA_SIZE, file_block_size_in_byte);
        offset += file_block_size_in_byte + META_DATA_SIZE;
    }

    for(int i = 0; i < oram->depth_file; ++i)
    {
        for(int j = 0; j < BUCKET_SIZE; ++j)
        {
            memcpy(oram->meta_data_file[full_path[i]] + j, out+offset, META_DATA_SIZE);
            offset           += META_DATA_SIZE;
            int position      = (full_path[i]<<1) + j;
            int byte_position = position >> 3;
            int bit_position  = (position % 8);
            uint8_t val       = *(uint8_t*)(out+offset);
            int count         = 0;
            // The following operation writes data row by row, which might cause memory access overhead
            for(int k = 0; k < oram->file_block_size/BUCKET_SIZE - 1; ++k)
            {
                for(int l = 0; l < BUCKET_SIZE; ++l)
                {
                    oram->data[k][l][byte_position] &= ~(1<<bit_position);
                    oram->data[k][l][byte_position] |= ((val & 0x80) >> (7-bit_position));
                    val                            <<= 1;
                }
                count += BUCKET_SIZE;
                if(count == 8) 
                {
                    offset++;
                    count = 0;
                    val   = *(uint8_t*)(out+offset);
                }                
            }        
            offset++;
        }
    }
}
#endif
#endif 