#ifndef CMPC_H__
#define CMPC_H__
#include "fpremp.h"
#include "abitmp.h"
#include "netmp.h"
#include "flexible_input_output.h"
#include <emp-tool/emp-tool.h>

using namespace emp;

template<int nP>
class CMPC { public:
	const static int SSP = 5;//5*8 in fact...
	const block MASK = makeBlock(0x0ULL, 0xFFFFFULL);
	FpreMP<nP>* fpre = nullptr;
	block* mac[nP+1];
	block* key[nP+1];
	bool* value;

	block * preprocess_mac[nP+1];
	block * preprocess_key[nP+1];
	bool* preprocess_value;

	block * sigma_mac[nP+1];
	block * sigma_key[nP+1];
	bool * sigma_value;

	block * ANDS_mac[nP+1];
	block * ANDS_key[nP+1];
	bool * ANDS_value;

	block * labels;
	BristolFormat * cf;
	NetIOMP<nP> * io;
	int num_ands = 0, num_in;
	int party, total_pre, ssp;
	ThreadPool * pool;
	block Delta;
		
	block (*GTM)[4][nP+1];
	block (*GTK)[4][nP+1];
	bool (*GTv)[4];
	block (*GT)[nP+1][4][nP+1];
	block * eval_labels[nP+1];
	PRP prp;

	CMPC(NetIOMP<nP> * io[2], ThreadPool * pool, int party, BristolFormat * cf, 
		int * num_and_gates, bool * _delta = nullptr, int ssp = 40) {
		
		this->party = party;
		this->io = io[0];
		this->cf = cf;
		this->ssp = ssp;
		this->pool = pool;

		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == AND_GATE)
				++num_ands;
		}
		
		if(num_and_gates != nullptr)
			*num_and_gates = num_ands;
		
		num_in = cf->n1+cf->n2;
		total_pre = num_in + num_ands + 3*ssp;
		fpre = new FpreMP<nP>(io, pool, party, _delta, ssp);
		Delta = fpre->Delta;

		if(party == 1) {
			// this->GTM = GTM;
			// this->GTK = GTK;
			// this->GTv = GTv;
			// this-> GT = GT;
			GTM = new block[num_ands][4][nP+1];
			GTK = new block[num_ands][4][nP+1];
			GTv = new bool[num_ands][4];
			GT = new block[num_ands][nP+1][4][nP+1];
		}
		
		this->labels = new block[cf->num_wire];
		this->value = new bool[cf->num_wire];

		for(int i = 1; i <= nP; ++i)
			eval_labels[i] = new block[cf->num_wire];
		
		for(int i = 1; i <= nP; ++i) {
			key[i] = new block[cf->num_wire];
			mac[i] = new block[cf->num_wire];
			ANDS_key[i] = new block[num_ands*3];
			ANDS_mac[i] = new block[num_ands*3];
			preprocess_mac[i] = new block[total_pre];
			preprocess_key[i] = new block[total_pre];
			sigma_mac[i] = new block[num_ands];
			sigma_key[i] = new block[num_ands]; 
			// eval_labels[i] = new block[cf->num_wire];
		}
		// value = new bool[cf->num_wire];
		ANDS_value = new bool[num_ands*3];
		preprocess_value = new bool[total_pre];
		sigma_value = new bool[num_ands];
	}
	~CMPC() {	
		// delete fpre;
		if(party == 1) {
			delete[] GTM;
			delete[] GTK;
			delete[] GTv;
			delete[] GT;
		} 
		delete[] labels;
		for(int i = 1; i <= nP; ++i) {
			/* delete[] key[i];
			delete[] mac[i];
			delete[] ANDS_key[i];
			delete[] ANDS_mac[i];
			delete[] preprocess_mac[i];
			delete[] preprocess_key[i];
			delete[] sigma_mac[i];
			delete[] sigma_key[i]; */
			delete[] eval_labels[i];
		}
		delete[] value;
		// delete[] ANDS_value;
		// delete[] preprocess_value;
		// delete[] sigma_value;
	}
	PRG prg;

	void function_independent() {
		if(party != 1)
			prg.random_block(labels, cf->num_wire);
		cout << "Step 1" << endl;	
		fpre->compute(ANDS_mac, ANDS_key, ANDS_value, num_ands);
		cout << "Step 2" << endl;
		prg.random_bool(preprocess_value, total_pre);
		cout << "Step 3" << endl;
		fpre->abit->compute(preprocess_mac, preprocess_key, preprocess_value, total_pre);
		cout << "Step 4" << endl;
		auto ret = fpre->abit->check(preprocess_mac, preprocess_key, preprocess_value, total_pre);
		ret.get();
		cout << "Step 5" << endl;
		for(int i = 1; i <= nP; ++i) {
			memcpy(key[i], preprocess_key[i], num_in * sizeof(block));
			memcpy(mac[i], preprocess_mac[i], num_in * sizeof(block));
		}
		memcpy(value, preprocess_value, num_in * sizeof(bool));
#ifdef __debug
		check_MAC<nP>(io, ANDS_mac, ANDS_key, ANDS_value, Delta, num_ands*3, party);
		check_correctness<nP>(io, ANDS_value, num_ands, party);
#endif
//		ret.get();
	}

	void function_dependent() {
		int ands = num_in;
		bool * x[nP+1];
		bool * y[nP+1];
		for(int i = 1; i <= nP; ++i) {
			x[i] = new bool[num_ands];
			y[i] = new bool[num_ands];
		}

		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == AND_GATE) {
				for(int j = 1; j <= nP; ++j) {
					key[j][cf->gates[4*i+2]] = preprocess_key[j][ands];
					mac[j][cf->gates[4*i+2]] = preprocess_mac[j][ands];
				}
				value[cf->gates[4*i+2]] = preprocess_value[ands];
				++ands;
			}
		}
		
		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == XOR_GATE) {
				for(int j = 1; j <= nP; ++j) {
					key[j][cf->gates[4*i+2]] = key[j][cf->gates[4*i]] ^ key[j][cf->gates[4*i+1]];
					mac[j][cf->gates[4*i+2]] = mac[j][cf->gates[4*i]] ^ mac[j][cf->gates[4*i+1]];
				}
				value[cf->gates[4*i+2]] = value[cf->gates[4*i]] != value[cf->gates[4*i+1]];
				if(party != 1)
					labels[cf->gates[4*i+2]] = labels[cf->gates[4*i]] ^ labels[cf->gates[4*i+1]];
			} else if (cf->gates[4*i+3] == NOT_GATE) {
				for(int j = 1; j <= nP; ++j) {
					key[j][cf->gates[4*i+2]] = key[j][cf->gates[4*i]];
					mac[j][cf->gates[4*i+2]] = mac[j][cf->gates[4*i]];
				}
				value[cf->gates[4*i+2]] = value[cf->gates[4*i]];
				if(party != 1)
					labels[cf->gates[4*i+2]] = labels[cf->gates[4*i]] ^ Delta;
			}
		}

#ifdef __debug
		check_MAC<nP>(io, mac, key, value, Delta, cf->num_wire, party);
#endif
		ands = 0;
		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == AND_GATE) {
				x[party][ands] = value[cf->gates[4*i]] != ANDS_value[3*ands];
				y[party][ands] = value[cf->gates[4*i+1]] != ANDS_value[3*ands+1];	
				ands++;
			}
		}
		
		vector<future<void>>	 res;
		for(int i = 1; i <= nP; ++i) for(int j = 1; j <= nP; ++j) if( (i < j) and (i == party or j == party) ) {
			int party2 = i + j - party;
			res.push_back(pool->enqueue([this, x, y, party2]() {
				io->send_data(party2, x[party], num_ands);
				io->flush(party2);
				io->send_data(party2, y[party], num_ands);
				io->flush(party2);
			}));
			res.push_back(pool->enqueue([this, x, y, party2]() {
				io->recv_data(party2, x[party2], num_ands);
				io->recv_data(party2, y[party2], num_ands);
			}));
		}
		joinNclean(res);
		for(int i = 2; i <= nP; ++i) for(int j = 0; j < num_ands; ++j) {
			x[1][j] = x[1][j] != x[i][j];
			y[1][j] = y[1][j] != y[i][j];
		}

		ands = 0;
		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == AND_GATE) {
				for(int j = 1; j <= nP; ++j) {
					sigma_mac[j][ands] = ANDS_mac[j][3*ands+2];
					sigma_key[j][ands] = ANDS_key[j][3*ands+2];
				}
				sigma_value[ands] = ANDS_value[3*ands+2];

				if(x[1][ands]) {
					for(int j = 1; j <= nP; ++j) {
						sigma_mac[j][ands] = sigma_mac[j][ands] ^ ANDS_mac[j][3*ands+1];
						sigma_key[j][ands] = sigma_key[j][ands] ^ ANDS_key[j][3*ands+1];
					}
					sigma_value[ands] = sigma_value[ands] != ANDS_value[3*ands+1];
				}
				if(y[1][ands]) {
					for(int j = 1; j <= nP; ++j) {
						sigma_mac[j][ands] = sigma_mac[j][ands] ^ ANDS_mac[j][3*ands];
						sigma_key[j][ands] = sigma_key[j][ands] ^ ANDS_key[j][3*ands];
					}
					sigma_value[ands] = sigma_value[ands] != ANDS_value[3*ands];
				}
				if(x[1][ands] and y[1][ands]) {
					if(party != 1)
						sigma_key[1][ands] = sigma_key[1][ands] ^ Delta;
					else
						sigma_value[ands] = not sigma_value[ands];
				}
				ands++;
			}
		}//sigma_[] stores the and of input wires to each AND gates
#ifdef __debug_
		check_MAC<nP>(io, sigma_mac, sigma_key, sigma_value, Delta, num_ands, party);
		ands = 0;
		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == AND_GATE) {
				bool tmp[] = { value[cf->gates[4*i]], value[cf->gates[4*i+1]], sigma_value[ands]};
				check_correctness(io, tmp, 1, party);
				ands++;
			}
		}
#endif
		
		ands = 0;
		block H[4][nP+1];
		block K[4][nP+1], M[4][nP+1];
		bool r[4];
		if(party != 1) { 
			for(int i = 0; i < cf->num_gate; ++i) if(cf->gates[4*i+3] == AND_GATE) {
				r[0] = sigma_value[ands] != value[cf->gates[4*i+2]];
				r[1] = r[0] != value[cf->gates[4*i]];
				r[2] = r[0] != value[cf->gates[4*i+1]];
				r[3] = r[1] != value[cf->gates[4*i+1]];

				for(int j = 1; j <= nP; ++j) {
					M[0][j] = sigma_mac[j][ands] ^ mac[j][cf->gates[4*i+2]];
					M[1][j] = M[0][j] ^ mac[j][cf->gates[4*i]];
					M[2][j] = M[0][j] ^ mac[j][cf->gates[4*i+1]];
					M[3][j] = M[1][j] ^ mac[j][cf->gates[4*i+1]];

					K[0][j] = sigma_key[j][ands] ^ key[j][cf->gates[4*i+2]];
					K[1][j] = K[0][j] ^ key[j][cf->gates[4*i]];
					K[2][j] = K[0][j] ^ key[j][cf->gates[4*i+1]];
					K[3][j] = K[1][j] ^ key[j][cf->gates[4*i+1]];
				}
				K[3][1] = K[3][1] ^ Delta;

				Hash(H, labels[cf->gates[4*i]], labels[cf->gates[4*i+1]], ands);
				for(int j = 0; j < 4; ++j) {
					for(int k = 1; k <= nP; ++k) if(k != party) {
						H[j][k] = H[j][k] ^ M[j][k];
						H[j][party] = H[j][party] ^ K[j][k];
					}
					H[j][party] = H[j][party] ^ labels[cf->gates[4*i+2]];
					if(r[j]) 
						H[j][party] = H[j][party] ^ Delta;
				}
				for(int j = 0; j < 4; ++j)
					io->send_data(1, H[j]+1, sizeof(block)*(nP));
				++ands;
			}
			io->flush(1);
		} else {
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, party2]() {
					for(int i = 0; i < num_ands; ++i)
						for(int j = 0; j < 4; ++j)
							io->recv_data(party2, GT[i][party2][j]+1, sizeof(block)*(nP));
				}));
			}
			for(int i = 0; i < cf->num_gate; ++i) if(cf->gates[4*i+3] == AND_GATE) {
				r[0] = sigma_value[ands] != value[cf->gates[4*i+2]];
				r[1] = r[0] != value[cf->gates[4*i]];
				r[2] = r[0] != value[cf->gates[4*i+1]];
				r[3] = r[1] != value[cf->gates[4*i+1]];
				r[3] = r[3] != true;

				for(int j = 1; j <= nP; ++j) {
					M[0][j] = sigma_mac[j][ands] ^ mac[j][cf->gates[4*i+2]];
					M[1][j] = M[0][j] ^ mac[j][cf->gates[4*i]];
					M[2][j] = M[0][j] ^ mac[j][cf->gates[4*i+1]];
					M[3][j] = M[1][j] ^ mac[j][cf->gates[4*i+1]];

					K[0][j] = sigma_key[j][ands] ^ key[j][cf->gates[4*i+2]];
					K[1][j] = K[0][j] ^ key[j][cf->gates[4*i]];
					K[2][j] = K[0][j] ^ key[j][cf->gates[4*i+1]];
					K[3][j] = K[1][j] ^ key[j][cf->gates[4*i+1]];
				}
				memcpy(GTK[ands], K, sizeof(block)*4*(nP+1));
				memcpy(GTM[ands], M, sizeof(block)*4*(nP+1));
				memcpy(GTv[ands], r, sizeof(bool)*4);
				++ands;
			}
			joinNclean(res);
		}
		for(int i = 1; i <= nP; ++i) {
			delete[] x[i];
			delete[] y[i];
		}
		delete fpre;
		for(int i = 1; i <= nP; ++i) {
			delete[] key[i];
			delete[] mac[i];
			delete[] ANDS_key[i];
			delete[] ANDS_mac[i];
			delete[] preprocess_mac[i];
			delete[] preprocess_key[i];
			delete[] sigma_mac[i];
			delete[] sigma_key[i];
		}
		delete[] ANDS_value;
		delete[] preprocess_value;
		delete[] sigma_value;
	}
	
	void save_preprocessing(string folder_name)
	{
		// write Delta 
		ofstream wf_Delta(folder_name + "/Delta_" + to_string(party) + ".dat", ios::out | ios::binary);
		wf_Delta.write((char*)&Delta, sizeof(Delta));
   		wf_Delta.close();

		// write value
		ofstream wf_Value(folder_name + "/value_" + to_string(party) + ".dat", ios::out | ios::binary);
		wf_Value.write((char*)value, sizeof(bool)*cf->num_wire);
   		wf_Value.close();

		// write labels
		ofstream wf_Labels(folder_name + "/labels_" + to_string(party) + ".dat", ios::out | ios::binary);
		wf_Labels.write((char*)labels, sizeof(block)*cf->num_wire);
   		wf_Labels.close();

		// These following for P1 only
		
		if(party == 1)
		{
			// write cf->gates 
			ofstream wf_Gates(folder_name + "/gates.dat", ios::out | ios::binary);
			wf_Gates.write((char*)cf->gates.data(), sizeof(int)*(cf->num_gate << 2));
			wf_Gates.close();
		
			// write GTM
			ofstream wf_GTM(folder_name + "/GTM.dat", ios::out | ios::binary);
			for(int i = 0; i < num_ands; ++i)
			{
				for(int j = 0; j < 4; ++j)
					wf_GTM.write((char*)GTM[i][j], (nP+1)*sizeof(block));
			}
			wf_GTM.close();

			// write GTv
			ofstream wf_GTv(folder_name + "/GTv.dat", ios::out | ios::binary);
			for(int i = 0; i < num_ands; ++i)
			{
				wf_GTv.write((char*)GTv[i], 4*sizeof(bool));
			}
			wf_GTv.close();

			// write GT
			ofstream wf_GT(folder_name + "/GT.dat", ios::out | ios::binary);
			for(int i = 0; i < num_ands; ++i)
			{
				for(int j = 0; j < nP+1; ++j)
				{
					for(int k = 0; k < 4; ++k)
						wf_GT.write((char*)GT[i][j][k], (nP+1)*sizeof(block));
				}
			}
			wf_GT.close();

			// write GTK
			ofstream wf_GTK(folder_name + "/GTK.dat", ios::out | ios::binary);
			for(int i = 0; i < num_ands; ++i)
			{
				for(int j = 0; j < 4; ++j)
					wf_GTK.write((char*)GTK[i][j], (nP+1)*sizeof(block));
			}
			wf_GTK.close();
		} 
	}

	void fast_preprocessing(string folder_name)
    {
        ifstream wf_Delta(folder_name + "/Delta_" + to_string(party) + ".dat", ios::in | ios::binary);
		if(!wf_Delta.is_open())
			cout << "Could not open\n";
		
        wf_Delta.read((char*)&Delta, sizeof(Delta));
        wf_Delta.close();

        ifstream wf_Value(folder_name + "/value_" + to_string(party) + ".dat", ios::in | ios::binary);
		if(!wf_Value.is_open())
			cout << "Could not open\n";
        wf_Value.read((char*)value, sizeof(bool)*cf->num_wire);
        wf_Value.close();

        ifstream wf_Labels(folder_name + "/labels_" + to_string(party) + ".dat", ios::in | ios::binary);
		if(!wf_Labels.is_open())
			cout << "Could not open\n";
        wf_Labels.read((char*)labels, sizeof(block)*cf->num_wire);
        wf_Labels.close();

        if(party == 1)
        {
            ifstream wf_Gates(folder_name + "/gates.dat", ios::in | ios::binary);
			if(!wf_Gates.is_open())
				cout << "Could not open\n";
            wf_Gates.read((char*)cf->gates.data(), sizeof(int)*(cf->num_gate<<2));
            wf_Gates.close();

            ifstream wf_GTM(folder_name + "/GTM.dat", ios::in | ios::binary);
			if(!wf_GTM.is_open())
				cout << "Could not open\n";
            for(int i = 0; i < num_ands; ++i)
            {
                for(int j = 0; j < 4; ++j)
                    wf_GTM.read((char*)GTM[i][j], (nP+1)*sizeof(block));
            }
            wf_GTM.close();

            ifstream wf_GTv(folder_name + "/GTv.dat", ios::in | ios::binary);
			if(!wf_GTv.is_open())
				cout << "Could not open\n";
            for(int i = 0; i < num_ands; ++i)
            {
                wf_GTv.read((char*)GTv[i], 4*sizeof(bool));
            }
            wf_GTv.close();

            ifstream wf_GT(folder_name + "/GT.dat", ios::in | ios::binary);
			if(!wf_GT.is_open())
				cout << "Could not open\n";
            for(int i = 0; i < num_ands; ++i)
            {
                for(int j = 0; j < nP + 1; ++j)
                {
                    for(int k = 0; k < 4; ++k)
                        wf_GT.read((char*)GT[i][j][k], (nP+1)*sizeof(block));
                }
            }
            wf_GT.close();

            ifstream wf_GTK(folder_name + "/GTK.dat", ios::in | ios::binary);
			if(!wf_GTK.is_open())
				cout << "Could not open\n";
            for(int i = 0; i < num_ands; ++i)
            {
                for(int j = 0; j < 4; ++j)
                    wf_GTK.read((char*)GTK[i][j], (nP+1)*sizeof(block));
            }
            wf_GTK.close();
        }
    }

	void online (uint8_t *input, uint8_t *output) 
	{
		//auto start_eval = clock_start();
		vector<future<void>> res;
		bool * mask_input = new bool[cf->num_wire];
		int total_bytes = num_in >> 3;
		
		int cnt = 0;
		for(int i = 0; i < total_bytes; ++i)
		{
			for(int j = 0; j < 8; ++j)
			{
				mask_input[cnt] = ((input[i] >> j) & 0x01) != value[cnt];
				cnt++;
			}
		}

		if(party != 1) 
		{
			io->send_data(1, mask_input, num_in);
			io->flush(1);
			io->recv_data(1, mask_input, num_in);
		} 
		else 
		{
			// Step 6. P1 receives data from other parties
			bool *tmp[nP+1];
			for(int i = 2; i <= nP; ++i) tmp[i] = new bool[num_in];
			for(int i = 2; i <= nP; ++i) 
			{
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() 
				{
					io->recv_data(party2, tmp[party2], num_in);
				}));
			}
			joinNclean(res);

			// Step 6. P1 computes x (xor) lambda then sends to parties Pi
			for(int i = 0; i < num_in; ++i)
				for(int j = 2; j <= nP; ++j)
					mask_input[i] = tmp[j][i] != mask_input[i];
			for(int i = 2; i <= nP; ++i) 
			{
				int party2 = i;
				res.push_back(pool->enqueue([this, mask_input, party2]() 
				{
					io->send_data(party2, mask_input, num_in);
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 2; i <= nP; ++i) delete[] tmp[i];
		}
		
		if(party!= 1) 
		{
			// Step 6. Parties send labels to P1
			block *tmp = new block[num_in];
			for(int i = 0; i < num_in; ++i) 
			{
				tmp[i] = labels[i];
				if(mask_input[i]) tmp[i] ^= Delta;
			}
			io->send_data(1, tmp, num_in*sizeof(block));
			io->flush(1);
			delete [] tmp;
		} 
		else 
		{
			for(int i = 2; i <= nP; ++i) 
			{
				int party2 = i;
				res.push_back(pool->enqueue([this, party2]() {
					io->recv_data(party2, eval_labels[party2], num_in*sizeof(block));
				}));
			}
			joinNclean(res);
			
			int ands = 0;
			int index = 0;
			block H[nP+1];
			block t0;
			int num_gates = cf->num_gate << 2;
			auto start = clock_start();	
			for(int i = 0; i < num_gates; i+=4) 
			{
				switch(cf->gates[i+3])
				{
					case XOR_GATE:
						mask_input[cf->gates[i+2]] = mask_input[cf->gates[i]] != mask_input[cf->gates[i+1]];
						for(int j = 2; j <= nP; ++j)
							eval_labels[j][cf->gates[i+2]] = eval_labels[j][cf->gates[i]] ^ eval_labels[j][cf->gates[i+1]];		
						break;
					case AND_GATE:
						index = (mask_input[cf->gates[i]] << 1) + mask_input[cf->gates[i+1]];
						
						for(int j = 2; j <= nP; ++j)
							eval_labels[j][cf->gates[i+2]] = GTM[ands][index][j];
						mask_input[cf->gates[i+2]] = GTv[ands][index];

						for(int j = 2; j <= nP; ++j) 
						{
							Hash(H, eval_labels[j][cf->gates[i]], eval_labels[j][cf->gates[i+1]], ands, index);
							xorBlocks_arr(H, H, GT[ands][j][index], nP+1);

							for(int k = 2; k <= nP; ++k)
								eval_labels[k][cf->gates[i+2]] = H[k] ^ eval_labels[k][cf->gates[i+2]];

							t0 = GTK[ands][index][j] ^ Delta;
							
							if(cmpBlock(&H[1], &GTK[ands][index][j], 1))
								mask_input[cf->gates[i+2]] = mask_input[cf->gates[i+2]] != false;
							else if(cmpBlock(&H[1], &t0, 1))
								mask_input[cf->gates[i+2]] = mask_input[cf->gates[i+2]] != true;
							else { cout << ands << "no match GT!" << endl << flush; }
						}
						ands++;
						break;
					default:
						mask_input[cf->gates[i+2]] = not mask_input[cf->gates[i]];
						for(int j = 2; j <= nP; ++j)
							eval_labels[j][cf->gates[i+2]] = eval_labels[j][cf->gates[i]];
						break;
				}
			}
			cout << "Evaluation time: " << time_from(start) << endl;
		}

		if(party != 1) 
		{
			total_bytes = cf->n3 >> 3;
			bool *result = value + cf->num_wire - cf->n3;
			cnt = 0;
			for(int i = 0; i < total_bytes; ++i)
			{
				output[i] = 0;
				for(int j = 0; j < 8; ++j)
				{
					output[i] = output[i] | ((result[cnt]) << j);
					cnt++;
				}
			}
		} 
		else 
		{
			total_bytes = cf->n3 >> 3;
			bool *v = value + cf->num_wire - cf->n3;
			bool *m = mask_input + cf->num_wire - cf->n3;
			cnt = 0;
			for(int i = 0; i < total_bytes; ++i)
			{
				output[i] = 0;
				for(int j = 0; j < 8; ++j)
				{
					output[i] = output[i] | ((v[cnt] != m[cnt]) << j);
					cnt++;
				}
			}
		}
		delete [] mask_input;
		//cout << "Evaluation time: " << time_from(start_eval) << endl;
	}
	
	void online(uint8_t *input, uint8_t *output, int num_ands) 
	{
		//auto start_eval = clock_start();
		vector<future<void>> res;
		bool * mask_input = new bool[cf->num_wire];
		int total_bytes  = num_in >> 3;
		int num_threads  = 1;

		int cnt = 0;
		for(int i = 0; i < total_bytes; ++i)
		{
			for(int j = 0; j < 8; ++j)
			{
				mask_input[cnt] = ((input[i] >> j) & 0x01) != value[cnt];
				cnt++;
			}
		}

		if(party != 1) 
		{
			io->send_data(1, mask_input, num_in);
			io->flush(1);
			io->recv_data(1, mask_input, num_in);
		} 
		else 
		{
			// Step 6. P1 receives data from other parties
			bool *tmp[nP+1];
			for(int i = 2; i <= nP; ++i) tmp[i] = new bool[num_in];
			for(int i = 2; i <= nP; ++i) 
			{
				int party2 = i;
				// res.push_back(pool->enqueue([this, tmp, party2]() 
				// {
					io->recv_data(party2, tmp[party2], num_in);
				// }));
			}
			// joinNclean(res);

			// Step 6. P1 computes x (xor) lambda then sends to parties Pi
			int start = 0;
			int end   = num_in/num_threads;

			auto start_masking = clock_start();

			for(int t = 0; t < num_threads; ++t)
			{
				// res.push_back(pool->enqueue([this, mask_input, tmp, start, end]() 
				// {
					for(int i = start; i < end; ++i)
						for(int j = 2; j <= nP; ++j)
							mask_input[i] = tmp[j][i] != mask_input[i];
				// }));
				start = end;
				end  += num_in/num_threads;
			}
			// joinNclean(res);

			cout << "Input masking time: " << time_from(start_masking) << endl;

			for(int i = 2; i <= nP; ++i) 
			{
				int party2 = i;
				// res.push_back(pool->enqueue([this, mask_input, party2]() 
				// {
					io->send_data(party2, mask_input, num_in);
					io->flush(party2);
				// }));
			}
			// joinNclean(res);
			for(int i = 2; i <= nP; ++i) delete[] tmp[i];
		}
		
		if(party != 1) 
		{
			// Step 6. Parties send labels to P1
			block *tmp = new block[num_in];
			int start = 0;
			int end   = num_in/num_threads;
			for(int t = 0; t < num_threads; ++t)
			{
				// res.push_back(pool->enqueue([this, mask_input, tmp, start, end]() 
				// {
					for(int i = start; i < end; ++i) 
					{
						tmp[i] = labels[i];
						if(mask_input[i]) tmp[i] ^= Delta;
					
					}
				// }));
				start = end;
				end  += num_in/num_threads;
			}
			// joinNclean(res);

			io->send_data(1, tmp, num_in*sizeof(block));
			io->flush(1);
			delete [] tmp;
		} 
		else 
		{
			for(int i = 2; i <= nP; ++i) 
			{
				int party2 = i;
				// res.push_back(pool->enqueue([this, party2]() {
					io->recv_data(party2, eval_labels[party2], num_in*sizeof(block));
				// }));
			}
			// joinNclean(res);
			
			int num_gates  = cf->num_gate << 2;
			int start_gate = 0;
			int end_gate   = num_gates/num_threads;
			int start_ands = 0;
			// ThreadPool pool(num_threads);

			auto start_eval = clock_start();

			for(int t = 0; t < num_threads; ++t) 
			{
				// res.push_back(pool.enqueue([this, start_gate, end_gate, mask_input, start_ands]() 
				// {
					// printf("start_gate: %d, end_gate: %d, num_gates: %d\n", start_gate, end_gate, end_gate-start_gate);
					auto start = clock_start();

					int index = 0;
					int ands  = start_ands;
					block H[nP+1];
					block t0;
					// double hash_time = 0;
					for(int i = start_gate; i < end_gate; i+=4) 
					{
						switch(cf->gates[i+3])
						{
							case XOR_GATE:
								mask_input[cf->gates[i+2]] = mask_input[cf->gates[i]] != mask_input[cf->gates[i+1]];
								for(int j = 2; j <= nP; ++j)
									eval_labels[j][cf->gates[i+2]] = eval_labels[j][cf->gates[i]] ^ eval_labels[j][cf->gates[i+1]];		
								break;
							case AND_GATE:
								index = (mask_input[cf->gates[i]] << 1) + mask_input[cf->gates[i+1]];
								
								for(int j = 2; j <= nP; ++j)
									eval_labels[j][cf->gates[i+2]] = GTM[ands][index][j];
								mask_input[cf->gates[i+2]] = GTv[ands][index];

								for(int j = 2; j <= nP; ++j) 
								{
									// auto start_hash = clock_start();
									Hash(H, eval_labels[j][cf->gates[i]], eval_labels[j][cf->gates[i+1]], ands, index);
									// hash_time += time_from(start_hash);
									xorBlocks_arr(H, H, GT[ands][j][index], nP+1);

									for(int k = 2; k <= nP; ++k)
										eval_labels[k][cf->gates[i+2]] = H[k] ^ eval_labels[k][cf->gates[i+2]];

									t0 = GTK[ands][index][j] ^ Delta;
									
									if(cmpBlock(&H[1], &GTK[ands][index][j], 1))
										mask_input[cf->gates[i+2]] = mask_input[cf->gates[i+2]] != false;
									else if(cmpBlock(&H[1], &t0, 1))
										mask_input[cf->gates[i+2]] = mask_input[cf->gates[i+2]] != true;
									else { cout << ands << "no match GT!" << endl << flush; }
								}
								ands++;
								break;
							default:
								mask_input[cf->gates[i+2]] = not mask_input[cf->gates[i]];
								for(int j = 2; j <= nP; ++j)
									eval_labels[j][cf->gates[i+2]] = eval_labels[j][cf->gates[i]];
								break;
						}
					}
					
					cout << "Time of thread is: " << time_from(start) << endl;
					// cout << "Hash time: " << hash_time << endl;
				// }));
				start_gate  = end_gate;
				end_gate   += num_gates/num_threads;
				start_ands += num_ands/num_threads;
			}
			// joinNclean(res);
			
			cout << "Circuit evaluation time: " << time_from(start_eval) << endl;
		}

		if(party != 1) 
		{
			total_bytes = cf->n3 >> 3;
			bool *result = value + cf->num_wire - cf->n3;
			cnt = 0;
			for(int i = 0; i < total_bytes; ++i)
			{
				output[i] = 0;
				for(int j = 0; j < 8; ++j)
				{
					output[i] = output[i] | ((result[cnt]) << j);
					cnt++;
				}
			}
		} 
		else 
		{
			total_bytes = cf->n3 >> 3;
			bool *v = value + cf->num_wire - cf->n3;
			bool *m = mask_input + cf->num_wire - cf->n3;
			cnt = 0;
			for(int i = 0; i < total_bytes; ++i)
			{
				output[i] = 0;
				for(int j = 0; j < 8; ++j)
				{
					output[i] = output[i] | ((v[cnt] != m[cnt]) << j);
					cnt++;
				}
			}
		}
		delete [] mask_input;
		//cout << "Evaluation time: " << time_from(start_eval) << endl;
	}
	
	void Hash(block H[4][nP+1], const block & a, const block & b, uint64_t idx) {
		block T[4];
		T[0] = sigma(a);
		T[1] = sigma(a ^ Delta);
		T[2] = sigma(sigma(b));
		T[3] = sigma(sigma(b ^ Delta));
		
		H[0][0] = T[0] ^ T[2];  
		H[1][0] = T[0] ^ T[3];  
		H[2][0] = T[1] ^ T[2];  
		H[3][0] = T[1] ^ T[3];  
		for(int j = 0; j < 4; ++j) for(int i = 1; i <= nP; ++i) {
			H[j][i] = H[j][0] ^ makeBlock(4*idx+j, i);
		}
		for(int j = 0; j < 4; ++j) {
			prp.permute_block(H[j]+1, nP);
		}
	}

	void Hash(block H[nP+1], const block &a, const block& b, uint64_t idx, uint64_t row) {
		H[0] = sigma(a) ^ sigma(sigma(b));
		int pos = (idx << 2) + row;
		for(int i = 1; i <= nP; ++i) {
			H[i] = H[0] ^ makeBlock(pos, i);
		}
		prp.permute_block(H+1, nP);
	}

	string tostring(bool a) {
		if(a) return "T";
		else return "F";
	}

	void online (bool * input, bool * output, int* start, int* end) {
		bool * mask_input = new bool[cf->num_wire];
		bool * input_mask[nP+1];
		for(int i = 0; i <= nP; ++i) input_mask[i] = new bool[end[party] - start[party]];
		memcpy(input_mask[party], value+start[party], end[party] - start[party]);
		memcpy(input_mask[0], input+start[party], end[party] - start[party]);

		vector<future<bool>> res;
		for(int i = 1; i <= nP; ++i) for(int j = 1; j<= nP; ++j) if( (i < j) and (i == party or j == party) ) {
			int party2 = i + j - party;
			res.push_back(pool->enqueue([this, start, end, party2]() {
				char dig[Hash::DIGEST_SIZE];
				io->send_data(party2, value+start[party2], end[party2]-start[party2]);
				emp::Hash::hash_once(dig, mac[party2]+start[party2], (end[party2]-start[party2])*sizeof(block));
				io->send_data(party2, dig, Hash::DIGEST_SIZE);
				io->flush(party2);
				return false;
			}));
			res.push_back(pool->enqueue([this, start, end, input_mask, party2]() {
				char dig[Hash::DIGEST_SIZE];
				char dig2[Hash::DIGEST_SIZE];
				io->recv_data(party2, input_mask[party2], end[party]-start[party]);
				block * tmp = new block[end[party]-start[party]];
				for(int i =  0; i < end[party] - start[party]; ++i) {
					tmp[i] = key[party2][i+start[party]];
					if(input_mask[party2][i])tmp[i] = tmp[i] ^ Delta;
				}
				emp::Hash::hash_once(dig2, tmp, (end[party]-start[party])*sizeof(block));
				io->recv_data(party2, dig, Hash::DIGEST_SIZE);
				delete[] tmp;
				return strncmp(dig, dig2, Hash::DIGEST_SIZE) != 0;	
			}));
		}
		if(joinNcleanCheat(res)) error("cheat!");
		for(int i = 1; i <= nP; ++i)
			for(int j = 0; j < end[party] - start[party]; ++j)
				input_mask[0][j] = input_mask[0][j] != input_mask[i][j];


		if(party != 1) {
			io->send_data(1, input_mask[0], end[party] - start[party]);
			io->flush(1);
			io->recv_data(1, mask_input, num_in);
		} else {
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, mask_input, start, end , party2]() {
					io->recv_data(party2, mask_input+start[party2], end[party2] - start[party2]);
				}));
			}
			joinNclean(res);
			memcpy(mask_input, input_mask[0], end[1]-start[1]);
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, mask_input, party2]() {
					io->send_data(party2, mask_input, num_in);
					io->flush(party2);
				}));
			}
			joinNclean(res);
		}
	
		if(party!= 1) {
			for(int i = 0; i < num_in; ++i) {
				block tmp = labels[i];
				if(mask_input[i]) tmp = tmp ^ Delta;
				io->send_data(1, &tmp, sizeof(block));
			}
			io->flush(1);
		} else {
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, party2]() {
					io->recv_data(party2, eval_labels[party2], num_in*sizeof(block));
				}));
			}
			joinNclean(res);
	
			int ands = 0;	
			for(int i = 0; i < cf->num_gate; ++i) {
				if (cf->gates[4*i+3] == XOR_GATE) {
					for(int j = 2; j <= nP; ++j)
						eval_labels[j][cf->gates[4*i+2]] = eval_labels[j][cf->gates[4*i]] ^ eval_labels[j][cf->gates[4*i+1]];
					mask_input[cf->gates[4*i+2]] = mask_input[cf->gates[4*i]] != mask_input[cf->gates[4*i+1]];
				} else if (cf->gates[4*i+3] == AND_GATE) {
					int index = 2*mask_input[cf->gates[4*i]] + mask_input[cf->gates[4*i+1]];
					block H[nP+1];
					for(int j = 2; j <= nP; ++j)
						eval_labels[j][cf->gates[4*i+2]] = GTM[ands][index][j];
					mask_input[cf->gates[4*i+2]] = GTv[ands][index];
					for(int j = 2; j <= nP; ++j) {
						Hash(H, eval_labels[j][cf->gates[4*i]], eval_labels[j][cf->gates[4*i+1]], ands, index);
						xorBlocks_arr(H, H, GT[ands][j][index], nP+1);
						for(int k = 2; k <= nP; ++k)
							eval_labels[k][cf->gates[4*i+2]] = H[k] ^ eval_labels[k][cf->gates[4*i+2]];
					
						block t0 = GTK[ands][index][j] ^ Delta;

						if(cmpBlock(&H[1], &GTK[ands][index][j], 1))
							mask_input[cf->gates[4*i+2]] = mask_input[cf->gates[4*i+2]] != false;
						else if(cmpBlock(&H[1], &t0, 1))
							mask_input[cf->gates[4*i+2]] = mask_input[cf->gates[4*i+2]] != true;
						else 	{cout <<ands <<"no match GT!"<<endl<<flush;
						}
					}
					ands++;
				} else {
					mask_input[cf->gates[4*i+2]] = not mask_input[cf->gates[4*i]];	
					for(int j = 2; j <= nP; ++j)
						eval_labels[j][cf->gates[4*i+2]] = eval_labels[j][cf->gates[4*i]];
				}
			}
		}
		if(party != 1) {
			io->send_data(1, value+cf->num_wire - cf->n3, cf->n3);
			io->flush(1);
		} else {
			vector<future<void>> res;
			bool * tmp[nP+1];
			for(int i = 2; i <= nP; ++i) 
				tmp[i] = new bool[cf->n3];
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], cf->n3);
				}));
			}
			joinNclean(res);
			for(int i = 0; i < cf->n3; ++i)
				for(int j = 2; j <= nP; ++j)
					mask_input[cf->num_wire - cf->n3 + i] = tmp[j][i] != mask_input[cf->num_wire - cf->n3 + i];
			for(int i = 0; i < cf->n3; ++i)
					mask_input[cf->num_wire - cf->n3 + i] = value[cf->num_wire - cf->n3 + i] != mask_input[cf->num_wire - cf->n3 + i];

			for(int i = 2; i <= nP; ++i) delete[] tmp[i];
			memcpy(output, mask_input + cf->num_wire - cf->n3, cf->n3);
		}
		delete[] mask_input;
	}

	void online (FlexIn<nP> * input, FlexOut<nP> *output) {
		bool * mask_input = new bool[cf->num_wire];
		input->associate_cmpc(pool, value, mac, key, io, Delta);
		input->input(mask_input);

		if(party!= 1) {
			for(int i = 0; i < num_in; ++i) {
				block tmp = labels[i];
				if(mask_input[i]) tmp = tmp ^ Delta;
				io->send_data(1, &tmp, sizeof(block));
			}
			io->flush(1);
		} else {
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, party2]() {
					io->recv_data(party2, eval_labels[party2], num_in*sizeof(block));
				}));
			}
			joinNclean(res);

			int ands = 0;
			for(int i = 0; i < cf->num_gate; ++i) {
				if (cf->gates[4*i+3] == XOR_GATE) {
					for(int j = 2; j<= nP; ++j)
						eval_labels[j][cf->gates[4*i+2]] = eval_labels[j][cf->gates[4*i]] ^ eval_labels[j][cf->gates[4*i+1]];
					mask_input[cf->gates[4*i+2]] = mask_input[cf->gates[4*i]] != mask_input[cf->gates[4*i+1]];
				} else if (cf->gates[4*i+3] == AND_GATE) {
					int index = 2*mask_input[cf->gates[4*i]] + mask_input[cf->gates[4*i+1]];
					block H[nP+1];
					for(int j = 2; j <= nP; ++j)
						eval_labels[j][cf->gates[4*i+2]] = GTM[ands][index][j];
					mask_input[cf->gates[4*i+2]] = GTv[ands][index];
					for(int j = 2; j <= nP; ++j) {
						Hash(H, eval_labels[j][cf->gates[4*i]], eval_labels[j][cf->gates[4*i+1]], ands, index);
						xorBlocks_arr(H, H, GT[ands][j][index], nP+1);
						for(int k = 2; k <= nP; ++k)
							eval_labels[k][cf->gates[4*i+2]] = H[k] ^ eval_labels[k][cf->gates[4*i+2]];

						block t0 = GTK[ands][index][j] ^ Delta;

						if(cmpBlock(&H[1], &GTK[ands][index][j], 1))
							mask_input[cf->gates[4*i+2]] = mask_input[cf->gates[4*i+2]] != false;
						else if(cmpBlock(&H[1], &t0, 1))
							mask_input[cf->gates[4*i+2]] = mask_input[cf->gates[4*i+2]] != true;
						else 	{cout <<ands <<"no match GT!"<<endl<<flush;
						}
					}
					ands++;
				} else {
					mask_input[cf->gates[4*i+2]] = not mask_input[cf->gates[4*i]];
					for(int j = 2; j <= nP; ++j)
						eval_labels[j][cf->gates[4*i+2]] = eval_labels[j][cf->gates[4*i]];
				}
			}
		}

		output->associate_cmpc(pool, value, mac, key, eval_labels, labels, io, Delta);
		output->output(mask_input, cf->num_wire - cf->n3);

		delete[] mask_input;
	}
};
#endif// CMPC_H__
