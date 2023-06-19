#include "Preprocessing.h"

vector<string> read_all_words(string file_name)
{
    ifstream ifs;
    string word;
    vector<string> result;
    ifs.open(file_name);
    int cnt = 0;
    int i = 0;
    while (!ifs.eof())
    {
        ifs >> word;
        transform(word.begin(), word.end(), word.begin(), ::tolower);
        result.push_back(word);
    }
    ifs.close();
    return result;
}

vector<uint8_t> create_file_and_row(uint32_t file_id, 
                                    uint32_t num_keywords, 
                                    uint32_t num_func_hash, 
                                    const vector<string> &dict, 
                                    vector<bool> &chosen_keywords, 
                                    int &count_keywords)
{
    static mt19937 mt(time(nullptr));
    ofstream ofs;
    vector<uint8_t> result(FILTER_SIZE, '0');
    ofs.open("database/" + to_string(file_id) + ".txt", ios::out);
    for (int i = 0; i < num_keywords; i++)
    {
        int rand_num = mt() % DICTIONARY_SIZE;
        string kw = dict[rand_num];
        if (chosen_keywords[rand_num] == false)
        {
            count_keywords++;
            chosen_keywords[rand_num] = true;
        }
        ofs << kw << endl;
        array<uint64_t, 2> hash_value = mm_hash((uint8_t*)kw.c_str(), kw.length());
        for (int j = 0; j < num_func_hash; j++) {
            uint64_t pos = nth_hash(j, hash_value[0], hash_value[1], FILTER_SIZE);
            result[pos] = '1';
        }
    }
    ofs.close();
    return result;
}

void create_database(uint32_t num_files, uint32_t num_keywords, uint32_t num_hash_func)
{
    vector<string> dict = read_all_words(DICTIONARY_FILE);
    vector<bool> chosen_keywords(DICTIONARY_SIZE, false);
    int count_keywords = 0;

    ofstream ofs;
    ofs.open(MATRIX_FILE, ios::out);

    for (int i = 0; i < num_files; i++)
    {
        cout << "File #" << i << endl;
        vector<uint8_t> new_row = create_file_and_row(i, num_keywords, num_hash_func, dict, chosen_keywords, count_keywords);
        for (uint8_t v : new_row)
            ofs << v;
        ofs << endl;
    }
    ofs.close();
    write_log("statistics.txt",
              "Created " + to_string(num_files) + " in database.\nThere are " + to_string(count_keywords) + " different keywords in these files.\n");
}

array<uint64_t, 2> mm_hash(const uint8_t* data, size_t len)
{
    array<uint64_t, 2> hash_value;
    MurmurHash3_x64_128(data, len, 0, hash_value.data());
    return hash_value;
}

inline uint64_t nth_hash(uint8_t i, uint64_t hash_a, uint64_t hash_b, uint64_t filter_size)
{
    return (hash_a + i * hash_b) % filter_size;
}

vector<vector<uint8_t>> read_indices_matrix(string file_name)
{
    ifstream ifs;
    string word;
    vector<vector<uint8_t>> matrix;
    ifs.open(file_name);
    int cnt = 0;
    while (!ifs.eof())
    {
        ifs >> word;
        cout << "Reading row index #" << cnt++ << endl;
        vector<uint8_t> row(FILTER_SIZE, 0);
        for(int i = 0; i < row.size(); i++)
            row[i] = word[i] - '0';
        matrix.push_back(row);
    }
    ifs.close();
    return matrix;
}

void write_log(string file_name, string log)
{
    ofstream ofs;
    ofs.open(file_name, ios::out);
    ofs << log;
    ofs.close();
}