#include <fstream>

#define  ID_SIZE         (32)

using    namespace std; 

uint64_t build_swapping_circuit(int num_blocks, int block_size, string file_name, int num_threads)
{
    ofstream ofs;
    ofs.open(file_name);
    
    block_size += 64;

    int start_out              = ID_SIZE + block_size + block_size*num_blocks;
    int num_intermediate_wire  = (ID_SIZE*3 - 1)*num_threads + block_size*3;
    int output_replacing_block = start_out + num_intermediate_wire*num_blocks;
    int output_holding_block   = output_replacing_block - block_size;

    ofs << (num_intermediate_wire+block_size)*num_blocks << " " << output_replacing_block + block_size*num_blocks << endl;
    ofs << block_size + ID_SIZE << " " << block_size*num_blocks << " " << block_size*(num_blocks+1) << endl;

    int start_thread_pos = 0;

    for(int t = 0; t < num_threads; ++t)
    {
        int start_p2              = block_size + ID_SIZE;
        int current_holding_block = ID_SIZE   + start_thread_pos;
        for(int n = 0; n < num_blocks; ++n)
        {
            // Create comparison circuit
            int start_p1 = 0;
            int start_out_xor = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 2 << " " << 1 << " " << start_p1 + i << " " << start_p2 + i << " " << start_out++ << " XOR" << endl;
            
            int start_out_inv = start_out;
            for(int i = 0; i < ID_SIZE; ++i) 
                ofs << 1 << " " << 1 << " " << start_out_xor + i << " " << start_out++ << " INV" << endl;
            
            int end_index = ID_SIZE;
            while(end_index > 1)
            {
                for(int i = 0; i < end_index; i += 2)
                    ofs << 2 << " " << 1 << " " << start_out_inv + i << " " << start_out_inv + i + 1 << " " << start_out++ << " AND" << endl;
                end_index >>= 1;
                start_out_inv = start_out - end_index;
            }
            int sel_signal = start_out - 1;
            
            // Create swapping circuit
            start_p1      = current_holding_block;
            start_out_xor = start_out;

            for(int i = 0; i < block_size/num_threads; ++i)
                ofs << 2 << " " << 1 << " " << start_p1 + i << " " << (start_p2 + start_thread_pos + i) << " " << start_out++ << " XOR" << endl;
            
            int start_out_and = start_out;
            for(int i = 0; i < block_size/num_threads; ++i)
                ofs << 2 << " " << 1 << " " << sel_signal << " " << start_out_xor + i << " " << start_out++ << " AND" << endl;
            
            if(n == num_blocks-1)
            {
                for(int i = 0; i < block_size/num_threads; ++i)
                    ofs << 2 << " " << 1 << " " << start_p1 + i << " " << (start_out_and + i) << " " << output_holding_block++ << " XOR" << endl;
            }
            else 
            {
                current_holding_block = start_out;
                for(int i = 0; i < block_size/num_threads; ++i)
                    ofs << 2 << " " << 1 << " " << start_p1 + i << " " << (start_out_and + i) << " " << start_out++ << " XOR" << endl;
            }
    
            for(int i = 0; i < block_size/num_threads; ++i)
                ofs << 2 << " " << 1 << " " << (start_p2 + start_thread_pos + i) << " " << (start_out_and + i) << " " << (output_replacing_block+n*block_size+start_thread_pos+i) << " XOR" << endl;
            
            start_p2 += block_size;
        }
        start_thread_pos += block_size/num_threads;
    }

    ofs.close();

    return (output_replacing_block + block_size*num_blocks);
}

uint64_t build_eviction_circuit(int num_blocks, int block_size, string file_name, int num_threads)
{
    ofstream ofs;
    ofs.open(file_name);
    
    block_size += 64;

    int start_out              = 32+(96+block_size)*num_blocks+block_size;
    int num_intermediate_wire  = (ID_SIZE*12+5)*num_threads + block_size*3;
    int output_replacing_block = start_out + num_intermediate_wire*num_blocks;
    
    ofs << (num_intermediate_wire+block_size)*num_blocks << " " << output_replacing_block + block_size*num_blocks << endl;
    ofs << 32+(96+block_size)*num_blocks << " " << block_size << " " << block_size*num_blocks << endl;
    
    int start_thread_pos = 0;
    for(int t = 0; t < num_threads; ++t)
    {
        int start_current_level    = 32;
        int start_deepest_block_id = 64;
        int start_target           = 96;
        int start_dest_level       = 0;
        int start_current_block    = 128;
        int holding_block          = 32 + (96+block_size)*num_blocks + t*block_size/num_threads;
        for(int n = 0; n < num_blocks; ++n)
        {
            // Circuit compares target with -1
            int start_and = start_out;
            for(int i = 0; i < ID_SIZE; i+=2)
                ofs << 2 << " " << 1 << " " << start_target+i << " " << start_target+i+1 << " " << start_out++ << " AND" << endl;
            
            int end_index = 16;
            while(end_index > 1)
            {
                for(int i = 0; i < end_index; i+=2)
                    ofs << 2 << " " << 1 << " " << start_and+i << " " << start_and+i+1 << " " << start_out++ << " AND" << endl;
                end_index >>= 1;
                start_and = start_out - end_index;
            }
            
            ofs << 1 << " " << 1 << " " << start_out-1 << " " << start_out << " INV" << endl;
            int is_target_diff_minus_1 = start_out;
            start_out++;
            
            // Circuit compares deepest block ID and block ID
            int start_out_xor = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 2 << " " << 1 << " " << start_deepest_block_id+i << " " << start_current_block+i << " " << start_out++ << " XOR" << endl;
                
            int start_out_inv = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 1 << " " << 1 << " " << start_out_xor+i << " " << start_out++ << " INV" << endl;
            
            end_index = ID_SIZE;
            while(end_index > 1)
            {
                for(int i = 0; i < end_index; i+=2)
                    ofs << 2 << " " << 1 << " " << start_out_inv+i << " " << start_out_inv+i+1 << " " << start_out++ << " AND" << endl;
                end_index >>= 1;
                start_out_inv = start_out - end_index;
            }
            int is_deepest_block_id = start_out-1;
            int sel_target_and_block_id = start_out;
            
            ofs << 2 << " " << 1 << " " << is_target_diff_minus_1 << " " << is_deepest_block_id << " " << start_out++ << " AND" << endl;
            
            // Circuit compares current level and destination level
            start_out_xor = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 2 << " " << 1 << " " << start_current_level+i << " " << start_dest_level+i << " " << start_out++ << " XOR" << endl;
                
            start_and = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 1 << " " << 1 << " " << start_out_xor+i << " " << start_out++ << " INV" << endl;
            
            end_index = ID_SIZE;
            while(end_index > 1)
            {
                for(int i = 0; i < end_index; i+=2)
                    ofs << 2 << " " << 1 << " " << start_and+i << " " << start_and+i+1 << " " << start_out++ << " AND" << endl;
                end_index >>= 1;
                start_and = start_out - end_index;
            }
            
            int is_dest_level = start_out - 1;
            
            // Comparison current block id and zero
            start_and = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 1 << " " << 1 << " " << start_current_block+i << " " << start_out++ << " INV" << endl;
            
            end_index = ID_SIZE;
            while(end_index > 1)
            {
                for(int i = 0; i < end_index; i+=2)
                    ofs << 2 << " " << 1 << " " << start_and+i << " " << start_and+i+1 << " " << start_out++ << " AND" << endl;
                end_index >>= 1;
                start_and = start_out - end_index;
            }
            
            int is_block_id_zero = start_out-1;
            
            ofs << 2 << " " << 1 << " " << is_dest_level << " " << is_block_id_zero << " " << start_out++ << " AND" << endl;
            
            int sel_level_block_id = start_out-1;
            
            // Compute sel signal for swapping
            ofs << 2 << " " << 1 << " " << sel_target_and_block_id << " " << sel_target_and_block_id << " " << start_out++ << " AND" << endl;
            ofs << 2 << " " << 1 << " " << sel_level_block_id << " " << sel_level_block_id << " " << start_out++ << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out-2 << " " << start_out << " INV" << endl;
            start_out++;
            ofs << 1 << " " << 1 << " " << start_out-2 << " " << start_out << " INV" << endl;
            start_out++;
            ofs << 2 << " " << 1 << " " << start_out-2 << " " << start_out-1 << " " << start_out << " AND" << endl;
            start_out++;
            ofs << 1 << " " << 1 << " " << start_out-1 << " " << start_out << " INV" << endl;
            start_out++;
            
            int sel_swapping = start_out-1;
            
            // Update destination level
            start_out_xor = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 2 << " " << 1 << " " << start_target+i << " " << start_dest_level+i << " " << start_out++ << " XOR" << endl;
            
            int saved_start_out_dest = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 2 << " " << 1 << " " << start_out_xor+i << " " << sel_target_and_block_id << " " << start_out++ << " AND" << endl;
            
            // Out destination level
            int updated_dest_level = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 2 << " " << 1 << " " << saved_start_out_dest+i << " " << start_dest_level+i << " " << start_out++ << " XOR" << endl;
            
            // Swapping    
            start_and = start_out;
            for(int i = 0; i < block_size/num_threads; ++i)
                ofs << 2 << " " << 1 << " " << (holding_block+i) << " " << (start_current_block+start_thread_pos+i) << " " << start_out++ << " XOR" << endl;
            
            int saved_start_out_and = start_out;
            for(int i = 0; i < block_size/num_threads; ++i)
                ofs << 2 << " " << 1 << " " << sel_swapping << " " << start_and + i << " " << start_out++ << " AND" << endl;
                
            // Out holding block
            int updated_holding_block = start_out;
            for(int i = 0; i < block_size/num_threads; ++i)
                ofs << 2 << " " << 1 << " " << (holding_block+i) << " " << saved_start_out_and+i << " " << start_out++ << " XOR" << endl;
            
            // Out update block
            for(int i = 0; i < block_size/num_threads; ++i)
                ofs << 2 << " " << 1 << " " << (start_current_block+start_thread_pos+i) << " " << saved_start_out_and+i << " " << (output_replacing_block+n*block_size+start_thread_pos+i) << " XOR" << endl;
            
            start_dest_level = updated_dest_level;
            holding_block    = updated_holding_block;
            
            start_target           += (ID_SIZE*3 + block_size);
            start_current_level    += (ID_SIZE*3 + block_size);
            start_deepest_block_id += (ID_SIZE*3 + block_size);
            start_current_block    += (ID_SIZE*3 + block_size);
        }
        start_thread_pos += block_size/num_threads;
    }
    
    ofs.close();  

    return (output_replacing_block + block_size*num_blocks);
}

void build_prepare_deepest_circuit(int num, string file_name)
{
    ofstream ofs;
    ofs.open(file_name);

    int start_goal            = 0;
    int start_src             = 32;
    int start_i               = 64;
    int start_deepest_level   = 96;
    int start_deepest         = 128;
    int start_idx             = 160;
    int start_out             = 64 + 128*num;
    int num_intermediate_wire = (ID_SIZE * 26) - 7;
    int output_deepest_value  = start_out + num_intermediate_wire*num;
    
    ofs << (num_intermediate_wire+ID_SIZE)*num << " " << output_deepest_value + ID_SIZE*num << endl;
    ofs << 64 << " " << 128*num << " " << ID_SIZE*num << endl;
    
    for(int n = 0; n < num; ++n)
    {
        // Compare i and goal
        int start_out_xor = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << start_goal + i << " " << start_i + i << " " << start_out++ << " XOR" << endl;
        }
        
        int start_out_and = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << start_out_xor + i << " " << start_goal + i << " " << start_out++ << " AND" << endl;
        }
        
        int input_or     = start_out_and;
        int start_or     = input_or + 1;
        
        for(int i = 0; i < ID_SIZE-1; ++i)
        {
            ofs << 2 << " " << 1 << " " << input_or << " " << input_or << " " << start_out + i << " AND" << endl;
            ofs << 2 << " " << 1 << " " << start_or + i << " " << start_or + i << " " << start_out + i + 1 << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i << " " << start_out + i + 2 << " INV" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i + 1 << " " << start_out + i + 3 << " INV" << endl;
            ofs << 2 << " " << 1 << " " << start_out + i + 2 << " " << start_out + i + 3 << " " << start_out + i + 4 << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i + 4 << " " << start_out + i + 5 << " INV" << endl;
            
            start_out += 5;
            input_or  = start_out + i;
        }
        
        int sel_cmpr_less_than = start_out + ID_SIZE-2;
        start_out = sel_cmpr_less_than + 1;
        
        // Create comparison equal circuit: i == goal?
        int start_equal = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_i + i << " " << start_goal + i << " " << start_out++ << " XOR" << endl;
        
        int start_inv = start_out;
        for(int i = 0; i < ID_SIZE; ++i) 
            ofs << 1 << " " << 1 << " " << start_equal + i << " " << start_out++ << " INV" << endl;
        
        int end_index = ID_SIZE;
        while(end_index > 1)
        {
            for(int i = 0; i < end_index; i += 2)
                ofs << 2 << " " << 1 << " " << start_inv + i << " " << start_inv + i + 1 << " " << start_out++ << " AND" << endl;
            end_index >>= 1;
            start_inv = start_out - end_index;
        }
        int sel_cmpr_equal = start_out - 1;
        
        // OR CIRCUIT
        ofs << 2 << " " << 1 << " " << sel_cmpr_less_than << " " << sel_cmpr_less_than << " " << start_out++ << " AND" << endl;
        ofs << 2 << " " << 1 << " " << sel_cmpr_equal << " " << sel_cmpr_equal << " " << start_out++ << " AND" << endl;
        ofs << 1 << " " << 1 << " " << start_out-2 << " " << start_out << " INV" << endl;
        start_out++;
        ofs << 1 << " " << 1 << " " << start_out-2 << " " << start_out << " INV" << endl;
        start_out++;
        ofs << 2 << " " << 1 << " " << start_out-2 << " " << start_out-1 << " " << start_out << " AND" << endl;
        start_out++;
        ofs << 1 << " " << 1 << " " << start_out-1 << " " << start_out << " INV" << endl;
        start_out++;
        
        int saved_cmpr_i_goal = start_out - 1;
        
        start_out_xor = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_src + i << " " << start_deepest + i << " " << start_out++ << " XOR" << endl;
        
        int saved_start_out_deepest = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_xor + i << " " << saved_cmpr_i_goal << " " << start_out++ << " AND" << endl;
        
        // Comparison circuit between goal and deepest_levels
        start_out_xor = start_out;
        int saved_start_out_xor_goal_deepest_level = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << start_goal + i << " " << start_deepest_level + i << " " << start_out++ << " XOR" << endl;
        }
        
        start_out_and = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << start_out_xor + i << " " << start_deepest_level + i << " " << start_out++ << " AND" << endl;
        }
        
        input_or     = start_out_and;
        start_or     = input_or + 1;
        
        for(int i = 0; i < ID_SIZE-1; ++i)
        {
            ofs << 2 << " " << 1 << " " << input_or << " " << input_or << " " << start_out + i << " AND" << endl;
            ofs << 2 << " " << 1 << " " << start_or + i << " " << start_or + i << " " << start_out +i + 1 << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i << " " << start_out + i + 2 << " INV" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i + 1 << " " << start_out + i + 3 << " INV" << endl;
            ofs << 2 << " " << 1 << " " << start_out + i + 2 << " " << start_out + i + 3 << " " << start_out + i + 4 << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i + 4 << " " << start_out + i + 5 << " INV" << endl;
            
            start_out += 5;
            input_or  = start_out + i;
        }
        
        int sel_goal_lt_deepest_level = start_out + ID_SIZE-2;
        start_out = sel_goal_lt_deepest_level + 1;
        start_out_xor = start_out;
        
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << start_idx + i << " " << start_src + i << " " << start_out++ << " XOR" << endl;
        }
        
        int saved_start_out_goal = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << saved_start_out_xor_goal_deepest_level + i << " " << sel_goal_lt_deepest_level << " " << start_out++ << " AND" << endl;
        }
        
        int saved_start_out_src = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << start_out_xor + i << " " << sel_goal_lt_deepest_level << " " << start_out++ << " AND" << endl;
        }
        
        // Output goal    
        int updated_goal = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << saved_start_out_goal + i << " " << start_goal + i << " " << start_out++ << " XOR" << endl;
        }
        
        // Output src
        int updated_src = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << saved_start_out_src + i << " " << start_src + i << " " << start_out++ << " XOR" << endl;
        }
        
        // Output deepest
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << saved_start_out_deepest + i << " " << start_deepest + i << " " << output_deepest_value++ << " XOR" << endl;
        }
        
        start_goal = updated_goal;
        start_src  = updated_src;
        
        start_i             += 128;
        start_deepest_level += 128;
        start_deepest       += 128;
        start_idx           += 128;
    }

    ofs.close();
}

void build_prepare_target_circuit(int num, string file_name)
{
    int start_src             = 0;
    int start_dest            = 32;
    int start_i               = 64;
    int start_empty_slot      = 96;
    int start_deepest         = 104;
    int start_target          = 136;
    int start_out             = 64 + 104*num;
    int num_intermediate_wire = (ID_SIZE*20) + 6;
    int output_target_value   = start_out + num_intermediate_wire*num;
    
    ofstream ofs;
    ofs.open(file_name);

    ofs << (num_intermediate_wire+ID_SIZE)*num << " " << output_target_value + ID_SIZE*num << endl;
    ofs << 64 << " " << 104*num << " " << ID_SIZE*num << endl;
    
    for(int n = 0; n < num; ++n)
    {
        // Create comparison circuit
        int start_out_xor = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_i + i << " " << start_src + i << " " << start_out++ << " XOR" << endl;
            
        int start_inv       = start_out;
        int saved_inv_i_src = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 1 << " " << 1 << " " << start_out_xor + i << " " << start_out++ << " INV" << endl;
        
        int end_index = ID_SIZE;
        while(end_index > 1)
        {
            for(int i = 0; i < end_index; i+=2)
                ofs << 2 << " " << 1 << " " << start_inv + i << " " << start_inv + i + 1 << " " << start_out++ << " AND" << endl;
            end_index >>= 1;
            start_inv = start_out - end_index;
        }
        int sel_cmpr_i_src = start_out-1;
        
        // Start of swapping between target and dest
        start_out_xor = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_target + i << " " << start_dest + i << " " << start_out++ << " XOR" << endl;
        
        int start_out_and = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_xor + i << " " << sel_cmpr_i_src << " " << start_out++ << " AND" << endl;
        
        int final_out_target = output_target_value;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_and + i << " " << start_target + i << " " << output_target_value++ << " XOR" << endl;
            
        // Swapping dest
        start_out_xor = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_dest + i << " " << saved_inv_i_src + i << " " << start_out++ << " XOR" << endl;
        
        start_out_and = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_xor + i << " " << sel_cmpr_i_src << " " << start_out++ << " AND" << endl;
        
        int saved_out_dest_1 = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_and + i << " " << start_dest + i << " " << start_out++ << " XOR" << endl;
        
        // Swapping src
        start_out_xor = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_src + i << " " << saved_inv_i_src + i << " " << start_out++ << " XOR" << endl;
        
        start_out_and = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_xor + i << " " << sel_cmpr_i_src << " " << start_out++ << " AND" << endl;
        
        int saved_out_src_1 = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_and + i << " " << start_src + i << " " << start_out++ << " XOR" << endl;
        
        // Compare dest with -1
        start_dest = saved_out_dest_1;
        for(int i = 0; i < ID_SIZE; i += 2)
            ofs << 2 << " " << 1 << " " << start_dest + i << " " << start_dest + i + 1 << " " << start_out++ << " AND" << endl;
        
        int start_and = start_out - ID_SIZE/2;
        end_index = ID_SIZE/2;
        while(end_index > 1)
        {
            for(int i = 0; i < end_index; i+=2)
                ofs << 2 << " " << 1 << " " << start_and + i << " " << start_and + i + 1 << " " << start_out++ << " AND" << endl;
            end_index >>= 1;
            start_and = start_out - end_index;
        }
           
        int sel_dest_is_neg_1 = start_out-1;
            
        // Compare empty slot with 1
        ofs << 2 << " " << 1 << " " << sel_dest_is_neg_1 << " " << start_empty_slot << " " << start_out++ << " AND" << endl;
        
        int sel_dest_is_neg_1_and_has_empty_slot = start_out - 1;
        
        // Compare target with -1
        int start_target_new = final_out_target;
        for(int i = 0; i < ID_SIZE; i += 2)
            ofs << 2 << " " << 1 << " " << start_target_new + i << " " << start_target_new + i + 1 << " " << start_out++ << " AND" << endl;
        
        start_and = start_out - ID_SIZE/2;
        end_index = ID_SIZE/2;
        while(end_index > 1)
        {
            for(int i = 0; i < end_index; i += 2)
                ofs << 2 << " " << 1 << " " << start_and + i << " " << start_and + i + 1 << " " << start_out++ << " AND" << endl;
            end_index >>= 1;
            start_and = start_out - end_index;
        }
        
        ofs << 1 << " " << 1 << " " << start_out-1 << " " << start_out << " INV" << endl;
        int sel_target_diff_neg_1 = start_out;
        start_out++ ;
        
        // Compute or of two first conditions
        ofs << 2 << " " << 1 << " " << sel_target_diff_neg_1 << " " << sel_target_diff_neg_1 << " " << start_out++ << " AND" << endl;
        ofs << 2 << " " << 1 << " " << sel_dest_is_neg_1_and_has_empty_slot << " " << sel_dest_is_neg_1_and_has_empty_slot << " " << start_out++ << " AND" << endl;
        ofs << 1 << " " << 1 << " " << start_out-2 << " " << start_out << " INV" << endl;
        start_out++;
        ofs << 1 << " " << 1 << " " << start_out-2 << " " << start_out << " INV" << endl;
        start_out++;
        ofs << 2 << " " << 1 << " " << start_out-2 << " " << start_out-1 << " " << start_out << " AND" << endl;
        start_out++;
        ofs << 1 << " " << 1 << " " << start_out-1 << " " << start_out << " INV" << endl;
        
        int sel_or_2_conditions = start_out;
        start_out++;
        
        // Compare deepest with -1
        start_and = start_out;
        for(int i = 0; i < ID_SIZE; i+=2)
            ofs << 2 << " " << 1 << " " << start_deepest + i << " " << start_deepest + i + 1 << " " << start_out++ << " AND" << endl;
        
        end_index = ID_SIZE/2;
        while(end_index > 1)
        {
            for(int i = 0; i < end_index; i+=2)
                ofs << 2 << " " << 1 << " " << start_and + i << " " << start_and + i + 1 << " " << start_out++ << " AND" << endl;
            end_index >>= 1;
            start_and = start_out - end_index;
        }
        
        ofs << 1 << " " << 1 << " " << start_out-1 << " " << start_out << " INV" << endl;
        int sel_deepest_diff_minus_1 = start_out;
        start_out++;
        
        ofs << 2 << " " << 1 << " " << sel_or_2_conditions << " " << sel_deepest_diff_minus_1 << " " << start_out++ << " AND" << endl;
        
        int sel_swapping = start_out-1;
        
        // Swapping src
        start_out_xor = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << saved_out_src_1 + i << " " << start_deepest + i << " " << start_out++ << " XOR" << endl;
        
        start_out_and = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_xor + i << " " << sel_swapping << " " << start_out++ << " AND" << endl;
        
        int updated_src = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_and + i << " " << saved_out_src_1 + i << " " << start_out++ << " XOR" << endl;
            
        // Swapping dest
        start_out_xor = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << saved_out_dest_1 + i << " " << start_i + i << " " << start_out++ << " XOR" << endl;
        
        start_out_and = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_xor + i << " " << sel_swapping << " " << start_out++ << " AND" << endl;
            
        int updated_dest = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_out_and + i << " " << saved_out_dest_1 + i << " " << start_out++ << " XOR" << endl;
        
        start_src  = updated_src;
        start_dest = updated_dest;
        
        start_i          += 104;
        start_empty_slot += 104;
        start_deepest    += 104;
        start_target     += 104;
    }

    ofs.close();
}

void build_deepest_circuit(int stash_size, int depth, int bucket_size, string file_name)
{
    int start_block_id        = 0;
    int start_deepest_level   = 32;
    int start_empty           = 64;
    int start_path_id         = 72;
    int current_block_id      = 104;
    int eviction_path         = 136;
    int start_out             = 72*(depth+1) + 96*(stash_size+depth*bucket_size);
    int num_intermediate_wire = (26*ID_SIZE-10)*(stash_size+depth*bucket_size)+72*(stash_size-1+depth*(bucket_size-1));
    int output                = start_out + num_intermediate_wire;
    
    ofstream ofs;
    ofs.open(file_name);

    ofs << num_intermediate_wire+72*(depth+1) << " " << output + 72*(depth+1) << endl;
    ofs << 72*(depth+1) << " " << 96*(stash_size+depth*bucket_size) << " " << 72*(depth+1) << endl;
    
    for(int n = 1; n <= stash_size+depth*bucket_size; ++n)
    {
        int start_out_xor = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << start_path_id + ID_SIZE - 1 - i << " " << eviction_path + ID_SIZE - 1 - i << " " << start_out++ << " XOR" << endl;
        }
        
        int start_out_inv = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 1 << " " << 1 << " " << start_out_xor + i << " " << start_out++ << " INV" << endl;
        }
        
        int input_or = start_out_xor;
        int start_or = input_or + 1;
        
        int saved_start_out_or = start_out;    
        for(int i = 0; i < ID_SIZE-1; ++i)
        {
            ofs << 2 << " " << 1 << " " << input_or << " " << input_or << " " << start_out + i << " AND" << endl;
            ofs << 2 << " " << 1 << " " << start_or + i << " " << start_or + i << " " << start_out + i + 1 << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i << " " << start_out + i + 2 << " INV" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i + 1 << " " << start_out + i + 3 << " INV" << endl;
            ofs << 2 << " " << 1 << " " << start_out + i + 2 << " " << start_out + i + 3 << " " << start_out + i + 4 << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i + 4 << " " << start_out + i + 5 << " INV" << endl;
            
            start_out += 5;
            input_or   = start_out + i;
        }
        
        start_out += ID_SIZE-1;
        ofs << 2 << " " << 1 << " " << start_out_inv << " " << start_out_inv << " " << start_out++ << " AND" << endl;
        start_out_inv++;
        
        int start_out_and = start_out - 1;
        for(int i = 0; i < ID_SIZE-1; ++i)
        {
            ofs << 2 << " " << 1 << " " << start_out_inv + i << " " << saved_start_out_or + i + 4 << " " << start_out++ << " AND" << endl;
            saved_start_out_or += 5;
        }
        
        // Compare current deepest level with new deepest level
        start_out_xor = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << start_deepest_level + i << " " << start_out_and + i << " " << start_out++ << " XOR" << endl;
        }
        
        int start_out_and2 = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
        {
            ofs << 2 << " " << 1 << " " << start_out_xor+i << " " << start_out_and + i << " " << start_out++ << " AND" << endl;
        }
        
        input_or = start_out_and2;
        start_or = input_or + 1;
        
        for(int i = 0; i < ID_SIZE-1; ++i)
        {
            ofs << 2 << " " << 1 << " " << input_or << " " << input_or << " " << start_out + i << " AND" << endl;
            ofs << 2 << " " << 1 << " " << start_or + i << " " << start_or + i << " " << start_out + i + 1 << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i << " " << start_out + i + 2 << " INV" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i + 1 << " " << start_out + i + 3 << " INV" << endl;
            ofs << 2 << " " << 1 << " " << start_out + i + 2 << " " << start_out + i + 3 << " " << start_out + i + 4 << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i + 4 << " " << start_out + i + 5 << " INV" << endl;
            
            start_out += 5;
            input_or   = start_out + i;
        }
        
        int saved_sel_less_than = start_out + ID_SIZE-2;
        start_out += ID_SIZE-1;
        
        input_or = current_block_id;
        start_or = input_or + 1;
        for(int i = 0; i < ID_SIZE-1; ++i)
        {
            ofs << 2 << " " << 1 << " " << input_or << " " << input_or << " " << start_out + i << " AND" << endl;
            ofs << 2 << " " << 1 << " " << start_or + i << " " << start_or + i << " " << start_out + i + 1 << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i << " " << start_out + i + 2 << " INV" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i + 1 << " " << start_out + i + 3 << " INV" << endl;
            ofs << 2 << " " << 1 << " " << start_out + i + 2 << " " << start_out + i + 3 << " " << start_out + i + 4 << " AND" << endl;
            ofs << 1 << " " << 1 << " " << start_out + i + 4 << " " << start_out + i + 5 << " INV" << endl;
            
            start_out += 5;
            input_or   = start_out + i;
        }
        
        int saved_sel_bid_diff_0 = start_out + ID_SIZE-2;
        
        ofs << 2 << " " << 1 << " " << saved_sel_less_than << " " << saved_sel_bid_diff_0 << " " << saved_sel_bid_diff_0 + 1 << " AND" << endl;
        
        start_out = saved_sel_bid_diff_0 + 2;
        
        int start_out_xor2 = start_out;
        
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_block_id + i << " " << current_block_id + i << " " << start_out++ << " XOR" << endl;
        
        int start_out_bid  = start_out;
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << saved_sel_bid_diff_0 + 1 << " " << start_out_xor2 + i << " " << start_out++ << " AND" << endl;
        
        for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << saved_sel_bid_diff_0 + 1 << " " << start_out_xor + i << " " << start_out++ << " AND" << endl;
        
        int start_out_empty = start_out;
        
        ofs << 1 << " " << 1 << " " << saved_sel_bid_diff_0 << " " << start_out << " INV" << endl;
        ofs << 2 << " " << 1 << " " << start_empty << " " << start_empty << " " << start_out + 1 << " AND" << endl;
        ofs << 2 << " " << 1 << " " << start_out << " " << start_out << " " << start_out + 2 << " AND" << endl;
        ofs << 1 << " " << 1 << " " << start_out + 1 << " " << start_out + 3 << " INV" << endl;
        ofs << 1 << " " << 1 << " " << start_out + 2 << " " << start_out + 4 << " INV" << endl;
        ofs << 2 << " " << 1 << " " << start_out + 3 << " " << start_out + 4 << " " << start_out + 5 << " AND" << endl;
        ofs << 1 << " " << 1 << " " << start_out + 5 << " " << start_out + 6 << " INV" << endl;
        
        start_out += 7;
        
        if( (n == stash_size) || (n > stash_size && ((n-stash_size)%bucket_size == 0)) )
        {
            for(int i = 0; i < ID_SIZE; ++i)
            ofs << 2 << " " << 1 << " " << start_block_id + i << " " << start_out_bid++ << " " << output++ << " XOR" << endl;

            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 2 << " " << 1 << " " << start_deepest_level + i << " " << start_out_bid++ << " " << output++ << " XOR" << endl;
            
            ofs << 2 << " " << 1 << " " << start_out_empty + 6 << " " << start_out_empty + 6 << " " << output++ << " AND" << endl;
            for(int i = 1; i < 8; ++i)
                ofs << 2 << " " << 1 << " " << start_empty + i << " " << start_empty + i << " " << output++ << " AND" << endl;
            
            start_block_id      = eviction_path  + 32;
            start_deepest_level = start_block_id + 32;
            start_empty         = start_block_id + 64;
            start_path_id       = start_block_id + 72;
            current_block_id    = start_block_id + 104;
            eviction_path       = start_block_id + 136;
        }
        else 
        {
            int updated_block_id = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 2 << " " << 1 << " " << start_block_id + i << " " << start_out_bid++ << " " << start_out++ << " XOR" << endl;
            
            int updated_deepest_level = start_out;
            for(int i = 0; i < ID_SIZE; ++i)
                ofs << 2 << " " << 1 << " " << start_deepest_level + i << " " << start_out_bid++ << " " << start_out++ << " XOR" << endl;
            
            int updated_empty = start_out;
            ofs << 2 << " " << 1 << " " << start_out_empty + 6 << " " << start_out_empty + 6 << " " << start_out++ << " AND" << endl;
            for(int i = 1; i < 8; ++i)
                ofs << 2 << " " << 1 << " " << start_empty + i << " " << start_empty + i << " " << start_out++ << " AND" << endl;
            
            start_block_id       =  updated_block_id;
            start_deepest_level  =  updated_deepest_level;
            start_empty          =  updated_empty;
            start_path_id       +=  96;
            current_block_id    +=  96;
            eviction_path       +=  96;
        }
    }

    ofs.close();
}

void build_and_2blocks(string file_name, int block_size, int num_threads)
{
    ofstream ofs;
    ofs.open(file_name);

    int start_p1  = 0;
    int start_p2  = start_p1 + block_size;
    int start_out = start_p2 + block_size;
    int start_pos_thread = 0;

    ofs << block_size << " " << block_size*3 << endl;
    ofs << block_size << " " << block_size << " " << block_size << endl;

    for(int t = 0; t < num_threads; ++t)
    {
        for(int i = 0; i < block_size/num_threads; ++i)
            ofs << 2 << " " << 1 << " " << (start_p1+start_pos_thread+i) << " " << (start_p2+start_pos_thread+i) << " " << start_out++ << " AND" << endl;

        start_pos_thread += block_size/num_threads;
    }    

    ofs.close();
}

void build_and_meta_data(string file_name, int num_ids, int num_threads)
{
    ofstream ofs;
    ofs.open(file_name);

    int start_p1  = 0;
    int start_p2  = num_ids;
    int start_out = start_p2 + (num_ids<<5);
    int start_pos = 0;
    int end_pos   = num_ids/num_threads;
    
    ofs << (num_ids<<5) << " " << (num_ids<<6)+num_ids << endl;
    ofs << num_ids << " " << (num_ids<<5) << " " << (num_ids<<5) << endl;

    for(int t = 0; t < num_threads; ++t)
    {
        for(int i = start_pos; i < end_pos; ++i)
        {
            for(int j = 0; j < 32; ++j) 
                ofs << 2 << " " << 1 << " " << (start_p1 + i) << " " << (start_p2 + i * 32 + j) << " " << start_out++ << " AND" << endl;
        }
        start_pos = end_pos;
        end_pos  += num_ids/num_threads;
    }    

    ofs.close();
}