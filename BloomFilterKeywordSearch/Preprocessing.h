#pragma once
#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <ctime>
#include <random>
#include <string>

#include "MurmurHash3.h"

#define DICTIONARY_FILE		"words.txt"
#define DICTIONARY_SIZE		466550
#define FILTER_SIZE			16384
#define NUM_HASH_FUNC		7
#define NUM_KEY_WORDS		300
#define NUM_FILES			65536
#define MATRIX_FILE			"indices_matrix.txt"

using namespace std;

vector<string> read_all_words(string file_name);
vector<uint8_t> create_file_and_row(uint32_t file_id, uint32_t num_keywords, uint32_t num_func_hash, const vector<string>& dict, vector<bool>& chosen_keywords, int& count_keywords);
void create_database(uint32_t num_files, uint32_t num_keywords, uint32_t num_hash_func);
array<uint64_t, 2> mm_hash(const uint8_t* data, size_t len);
inline uint64_t nth_hash(uint8_t i, uint64_t hash_a, uint64_t hash_b, uint64_t filter_size);
vector<vector<uint8_t>> read_indices_matrix(string file_name);
void write_log(string file_name, string log);