// BloomFilterKeywordSearch.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Reference: findingprotopia.org/posts/how-to-write-a-bloom-filter-cpp/

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <ctime>
#include <random>
#include <string>
#include <cstring>

#include "Preprocessing.h"

random_device rd;
mt19937       mt    (rd());

bool is_false_positive(string keyword, string file_name)
{
    ifstream ifs;
    ifs.open(file_name);
    string word;
    while (!ifs.eof())
    {
        ifs >> word;
        std::transform(word.begin(), word.end(), word.begin(), ::tolower);
        if (word.compare(keyword) == 0)
        {
            ifs.close();
            return false;
        }
    }
    ifs.close();
    return true;
}

void search_keyword(string keyword, vector<vector<uint8_t>> &indices_matrix)
{
    array<uint64_t, 2> hash_value = mm_hash((uint8_t*)keyword.c_str(), keyword.length());
    vector<uint32_t> pos;
    int count_total = 0, count_false_positive = 0;

    for (int i = 0; i < NUM_HASH_FUNC; i++)
        pos.push_back(nth_hash(i, hash_value[0], hash_value[1], FILTER_SIZE));

    cout << "\n\n=====================================================\n";
    cout << "Keyword \"" << keyword << "\" appears in these file IDs: ";

    for (int i = 0; i < NUM_FILES; i++)
    {
        bool appear = true;
        for (int j = 0; j < NUM_HASH_FUNC; j++)
        {
            if (indices_matrix[i][pos[j]] == 0)
            {
                appear = false;
                break;
            }
        }
        if (appear)
        {
            count_total++;
            cout << i << " ";
            if (is_false_positive(keyword, "database/" + to_string(i) + ".txt"))
                count_false_positive++;
        }
    }
    cout << "\nThere are " << count_total << " files containing \"" << keyword << "\" in database.\n";
    cout << "There are " << count_false_positive << " false positives in results.\n";
    cout << "=====================================================\n\n";
}

void build_omat(const vector<vector<uint8_t>> &indices_matrix, uint8_t **omat)
{
    uint8_t row[NUM_FILES>>3];

    int cnt = 0;
    for(int i = 0; i < FILTER_SIZE; ++i)
    {
        memset(row, 0, NUM_FILES>>3);
        for(int j = 0; j < NUM_FILES; ++j)
            row[j>>3] = (row[j>>3] << 1) | indices_matrix[j][i];
        memcpy(omat[cnt++], row, NUM_FILES>>3);
    }
}

void search_omat(string keyword, uint8_t **omat)
{
    array<uint64_t, 2> hash_value = mm_hash((uint8_t*)keyword.c_str(), keyword.length());
    vector<uint32_t> pos;
    int count_total = 0, count_false_positive = 0;

    for (int i = 0; i < NUM_HASH_FUNC; i++)
        pos.push_back(nth_hash(i, hash_value[0], hash_value[1], FILTER_SIZE));

    cout << "\n\n=====================================================\n" << flush;
    cout << "Keyword \"" << keyword << "\" appears in these file IDs: " << flush;

    for (int i = 0; i < NUM_FILES; i++)
    {
        bool appear = true;
        int offset = i % 8;
        for (int j = 0; j < NUM_HASH_FUNC; j++)
        {
            if (((omat[pos[j]][i>>3] << offset) & 0x80) == 0)
            {
                appear = false;
                break;
            }
        }
        if (appear)
        {
            count_total++;
            cout << i << " ";
            if (is_false_positive(keyword, "database/" + to_string(i) + ".txt"))
               count_false_positive++;
        }
    }

    cout << "\nThere are " << count_total << " files containing \"" << keyword << "\" in database.\n" << flush;
    cout << "There are " << count_false_positive << " false positives in results.\n" << flush;
    cout << "=====================================================\n\n" << flush;
}

int main()
{
    uint8_t **omat = new uint8_t*[FILTER_SIZE];
    for(int i = 0; i < FILTER_SIZE; ++i)
        omat[i] = new uint8_t[NUM_FILES>>3];

    // Prepocessing: 
    // 1. Create database files and indices matrix
    create_database(NUM_FILES, NUM_KEY_WORDS, NUM_HASH_FUNC);

    // 2. Read indices matrix  
    vector<vector<uint8_t>> indices_matrix = read_indices_matrix(MATRIX_FILE);

    // 3. Build OMAT
    build_omat(indices_matrix, omat);

    FILE *out_omat = fopen("omat.bin", "wb");
    for(int i = 0; i < FILTER_SIZE; ++i)
        fwrite(omat[i], 1, NUM_FILES>>3, out_omat);
    fclose(out_omat); 

    /* search_omat("security", omat);
    search_omat("realize", omat);
    search_omat("design", omat);
    search_omat("beautiful", omat);
    search_omat("coffee", omat);
    search_omat("book", omat);
    search_omat("mobile", omat);
    search_omat("mountain", omat);  */
    
    // 3. Give it a try
    search_keyword("security", indices_matrix);
    search_keyword("realize", indices_matrix);
    search_keyword("design", indices_matrix);
    search_keyword("beautiful", indices_matrix);
    search_keyword("coffee", indices_matrix);
    search_keyword("mobile", indices_matrix);
    search_keyword("mountain", indices_matrix);
    search_keyword("realize", indices_matrix);
    search_keyword("standup", indices_matrix);
    search_keyword("hello", indices_matrix);
    search_keyword("holding", indices_matrix);
    search_keyword("design", indices_matrix);
    search_keyword("cryptography", indices_matrix);
    search_keyword("beautiful", indices_matrix);
    search_keyword("tea", indices_matrix);
    search_keyword("coffee", indices_matrix);
    search_keyword("culture", indices_matrix);
    search_keyword("book", indices_matrix);
    search_keyword("mobile", indices_matrix);
    search_keyword("laptop", indices_matrix);
    search_keyword("iphone", indices_matrix);
    search_keyword("monitor", indices_matrix);
    search_keyword("hacker", indices_matrix);
    search_keyword("security", indices_matrix);
    search_keyword("attendance", indices_matrix);
    search_keyword("mountain", indices_matrix);
    search_keyword("subject", indices_matrix);
    
    return 0;
}