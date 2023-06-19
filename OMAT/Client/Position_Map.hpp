#include <iostream>

using namespace std;

class Position_Map
{
    public:
        int *path_id;
        
        Position_Map()
        {
            path_id = nullptr;
        }

        Position_Map(int size)
        {
            this->path_id = new int[size];
        }
        
        void update_position_map(int index, int path_id)
        {
            this->path_id[index] = path_id;
        }
};