#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <valarray>
#include <algorithm>
#include "./solvers/cpp/solver.h"
#define SIZE 100
const double epsilon = 1e-6;
using namespace std;
typedef struct _gate{
    int id;
    double x;
    double y;
    vector<int> connected_nets;
} gate;

typedef struct _net{
    // int id;
    int total_gates_pads;
    vector<int> connected_gates;
    vector<int> connected_pads;
    vector<int> connected_fake_pads;
} net;

typedef struct _pad{
    int connected_net_id;
    double y;
    double x;
} pad;

int xcompare(gate a, gate b) {
    if (a.x < b.x) return 1;
    if (a.x > b.x) return 0;
    if (a.y < b.y) return 1;
    else return 0;
}

int ycompare(gate a, gate b) {
    if (a.y < b.y) return 1;
    if (a.y > b.y) return 0;
    if (a.x < b.x) return 1;
    else return 0;
}

int main(int argc, char** argv) {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    if (argc < 2 || argc > 3) {
        cerr << "Error: Input Invalid!" << endl;
        return 1;
    }
    bool use_clique_weight = false;
    const char clique[] = "--clique";
    if (argc == 3) {
        if (!strcmp(argv[2], clique)) {
            use_clique_weight = true;
            cout << "Using Clique Weight" << endl;
        }
        else {
            cerr << "Error: Input Invalid!" << argv[2] << endl;
            return 1;
        }
    }
    else cout << "Using Default Weight" << endl;
    ifstream inputFile(argv[1]);
    if (!inputFile.is_open()) {
        cerr << "Error: Cannot Open Input File " << argv[1] << endl;
    }
    int num_gates;
    int num_nets;
    inputFile >> num_gates >> num_nets;
    vector<gate> gates(num_gates+1);
    vector<net> nets(num_nets+1);
    for (int i = 1; i <= num_gates; i++) {
        int gate_id, num_connected_nets;
        inputFile >> gate_id >> num_connected_nets;
        gates[gate_id].id = gate_id;
        for (int j = 1; j <= num_connected_nets; j++) {
            int connected_net_id;
            inputFile >> connected_net_id;
            gates[gate_id].connected_nets.push_back(connected_net_id);
            nets[connected_net_id].connected_gates.push_back(gate_id);
            nets[connected_net_id].total_gates_pads++;
        }
    }
    vector<int> inverse_gates(num_gates+1); // gate id를 넣으면 index와 매칭시킴. 

    int num_pads;
    inputFile >> num_pads;
    vector<pad> pads(num_pads+1);
    for (int i = 1; i <= num_pads; i++) {
        int pad_id, connected_net_id;
        inputFile >> pad_id >> connected_net_id;
        pads[pad_id].connected_net_id = connected_net_id;
        nets[connected_net_id].connected_pads.push_back(pad_id);
        nets[connected_net_id].total_gates_pads++;
        inputFile >> pads[pad_id].x >> pads[pad_id].y;
    }
    inputFile.close();
    
    int exponential_iter; // QP1까지는 1, QP3까지는 2, ...
    cout << "Enter number of exponential iterations: " << endl;
    cin >> exponential_iter;

    vector<int> QPstart(1<<exponential_iter);
    vector<int> QPend(1<<exponential_iter);
    vector<int> QPnum_gates(1<<exponential_iter);
    vector<double> min_x(1<<exponential_iter);
    vector<double> min_y(1<<exponential_iter);
    vector<double> max_x(1<<exponential_iter);
    vector<double> max_y(1<<exponential_iter);
    QPstart[1] = 1;
    QPend[1] = num_gates;
    QPnum_gates[1] = num_gates;
    min_x[1] = 0;
    max_x[1] = SIZE;
    min_y[1] = 0;
    max_y[1] = SIZE;

    for (int ei = 1; ei <= exponential_iter; ei++) { // exponential i
        int is_y = ei % 2;
        for (int iter = 1 << (ei - 1); iter < 1 << ei; iter++) {
            coo_matrix A;
            valarray<double> b_x(0.0, QPnum_gates[iter]);
            valarray<double> b_y(0.0, QPnum_gates[iter]);
            
            vector<int> A_rows, A_cols;
            vector<double> A_vals;

            // cout << "QP" << iter << " start" << endl;
            // debugging 
            // cout << "Gates in this QP: " << QPnum_gates[iter] << endl;
            // cout << "Boundary Gates: " << QPstart[iter] << " to " << QPend[iter] << endl;
            // cout << "Physical Boundaries: (" << min_x[iter] << ", " << min_y[iter] << ") to (" << max_x[iter] << ", " << max_y[iter] << ")" << endl;

            for (int i = 0; i < QPnum_gates[iter]; i++) {
                // 대각성분 계산
                int total_weight = 0;
                for (int j : gates[QPstart[iter] + i].connected_nets) {
                    if (use_clique_weight) total_weight += 1;
                    else total_weight += nets[j].total_gates_pads - 1;
                }
                A_rows.push_back(i);
                A_cols.push_back(i);
                A_vals.push_back(total_weight);
                
                vector<double> um(QPnum_gates[iter], 0.0); // 두 개의 gate가 여러 net으로 연결되어 있을 수 있으므로 weight를 누적시켜야 함. 
                for (int j : gates[QPstart[iter] + i].connected_nets) {
                    // nondiagonal 계산
                    for (int k : nets[j].connected_gates) { // k는 gate id, i와 다름
                        if (gates[QPstart[iter] + i].id != k) {
                            if (inverse_gates[k] >= QPstart[iter] && inverse_gates[k] <= QPend[iter]) { // index가 범위 안에 있음, k는 같은 QP gate임. 
                                /* // weight 누적시키지 않는 경우
                                A_rows.push_back(i);
                                A_cols.push_back(inverse_gates[k]-QPstart[iter]);
                                if (use_clique_weight) A_vals.push_back(-1/(nets[j].total_gates_pads - 1)); // weight = 1 / k-1 (using clique method)
                                else A_vals.push_back(-1);
                                */
                                um[inverse_gates[k]-QPstart[iter]] += (use_clique_weight ? -1.0/(nets[j].total_gates_pads - 1) : -1.0);
                            }
                            else { // index가 범위 밖에 있음. 즉 k는 다른 QP gate이므로 pad로 만들어 vector b에 추가. 
                                if (use_clique_weight) {
                                    if (gates[inverse_gates[k]].x > max_x[iter]) b_x[i] += max_x[iter] / (nets[j].total_gates_pads - 1);
                                    else if (gates[inverse_gates[k]].x < min_x[iter]) b_x[i] += min_x[iter] / (nets[j].total_gates_pads - 1);
                                    else b_x[i] += gates[inverse_gates[k]].x / (nets[j].total_gates_pads - 1);
                                    if (gates[inverse_gates[k]].y > max_y[iter]) b_y[i] += max_y[iter] / (nets[j].total_gates_pads - 1);
                                    else if (gates[inverse_gates[k]].y < min_y[iter]) b_y[i] += min_y[iter] / (nets[j].total_gates_pads - 1);
                                    else b_y[i] += gates[inverse_gates[k]].y / (nets[j].total_gates_pads - 1);
                                }
                                else {
                                    if (gates[inverse_gates[k]].x > max_x[iter]) b_x[i] += max_x[iter];
                                    else if (gates[inverse_gates[k]].x < min_x[iter]) b_x[i] += min_x[iter];
                                    else b_x[i] += gates[inverse_gates[k]].x;
                                    if (gates[inverse_gates[k]].y > max_y[iter]) b_y[i] += max_y[iter];
                                    else if (gates[inverse_gates[k]].y < min_y[iter]) b_y[i] += min_y[iter];
                                    else b_y[i] += gates[inverse_gates[k]].y;
                                }
                            }
                        }
                    }

                    // vector b 계산
                    for (int k : nets[j].connected_pads) { // pad는 앞과 변화 X
                        if (use_clique_weight) {
                            if (pads[k].x > max_x[iter]) b_x[i] += max_x[iter] / (nets[j].total_gates_pads - 1);
                            else if (pads[k].x < min_x[iter]) b_x[i] += min_x[iter] / (nets[j].total_gates_pads - 1);
                            else b_x[i] += pads[k].x / (nets[j].total_gates_pads - 1);
                            if (pads[k].y > max_y[iter]) b_y[i] += max_y[iter] / (nets[j].total_gates_pads - 1);
                            else if (pads[k].y < min_y[iter]) b_y[i] += min_y[iter] / (nets[j].total_gates_pads - 1);
                            else b_y[i] += pads[k].y / (nets[j].total_gates_pads - 1);
                        }
                        else {
                            if (pads[k].x > max_x[iter]) b_x[i] += max_x[iter];
                            else if (pads[k].x < min_x[iter]) b_x[i] += min_x[iter];
                            else b_x[i] += pads[k].x;
                            if (pads[k].y > max_y[iter]) b_y[i] += max_y[iter];
                            else if (pads[k].y < min_y[iter]) b_y[i] += min_y[iter];
                            else b_y[i] += pads[k].y;
                        }
                    }
                }
                for (int k = 0; k < QPnum_gates[iter]; k++) {
                    if (um[k] != 0.0) {
                        A_rows.push_back(i);
                        A_cols.push_back(k);
                        A_vals.push_back(um[k]);
                    }
                }
            }
            A.n = QPnum_gates[iter];
            A.nnz = A_vals.size();
            A.row.resize(A.nnz);
            A.col.resize(A.nnz);
            A.dat.resize(A.nnz);
            for (int i = 0; i < A.nnz; i++) {
                A.row[i] = A_rows[i];
                A.col[i] = A_cols[i];
                A.dat[i] = A_vals[i];
            }
            valarray<double> x_coords(0.0, QPnum_gates[iter]);
            valarray<double> y_coords(0.0, QPnum_gates[iter]);

            A.solve(b_x, x_coords);
            A.solve(b_y, y_coords);

            // min이나 max 넘어가는 경우 처리. epsilon은 10^-6
            for (int i = 0; i < QPnum_gates[iter]; i++) {
                gates[i + QPstart[iter]].x = x_coords[i];
                if (abs(x_coords[i] - min_x[iter]) < epsilon) gates[i + QPstart[iter]].x = min_x[iter];
                else if (abs(x_coords[i] - max_x[iter]) < epsilon) gates[i + QPstart[iter]].x = max_x[iter];
                gates[i + QPstart[iter]].y = y_coords[i];
                if (abs(y_coords[i] - min_y[iter]) < epsilon) gates[i + QPstart[iter]].y = min_y[iter];
                else if (abs(y_coords[i] - max_y[iter]) < epsilon) gates[i + QPstart[iter]].y = max_y[iter];
                // cout << gates[i + QPstart[iter]].id << " " << gates[i + QPstart[iter]].x << " " << gates[i + QPstart[iter]].y << endl;
            }

            if (ei == exponential_iter) continue; // 마지막 iteration이면 partition 나눌 필요 X

            // 다음 partition 계산
            QPstart[iter * 2] = QPstart[iter];
            QPend[iter * 2] = (QPstart[iter] + QPend[iter] - 1) / 2;
            QPstart[iter * 2 + 1] = QPend[iter * 2] + 1;
            QPend[iter * 2 + 1] = QPend[iter];
            QPnum_gates[iter * 2] = QPend[iter * 2] - QPstart[iter * 2] + 1;
            QPnum_gates[iter * 2 + 1] = QPend[iter * 2 + 1] - QPstart[iter * 2 + 1] + 1;

            // 다음 물리적 partition 계산
            if (is_y) {
                min_x[iter * 2] = min_x[iter];
                max_x[iter * 2] = (min_x[iter] + max_x[iter]) / 2;
                min_x[iter * 2 + 1] = (min_x[iter] + max_x[iter]) / 2;
                max_x[iter * 2 + 1] = max_x[iter];
                min_y[iter * 2] = min_y[iter];
                max_y[iter * 2] = max_y[iter];
                min_y[iter * 2 + 1] = min_y[iter];
                max_y[iter * 2 + 1] = max_y[iter];
            }
            else {
                min_x[iter * 2] = min_x[iter];
                max_x[iter * 2] = max_x[iter];
                min_x[iter * 2 + 1] = min_x[iter];
                max_x[iter * 2 + 1] = max_x[iter];
                min_y[iter * 2] = min_y[iter];
                max_y[iter * 2] = (min_y[iter] + max_y[iter]) / 2;
                min_y[iter * 2 + 1] = (min_y[iter] + max_y[iter]) / 2;
                max_y[iter * 2 + 1] = max_y[iter];
            }

            // dividing. 더 정확하게는 sort하고 index divide. 현재가 y로 분할이라면 그 다음번은 x로 분할, x로 값들을 정렬해야함. 
            if(is_y) sort(gates.begin() + QPstart[iter], gates.begin() + QPend[iter] + 1, xcompare);
            else sort(gates.begin() + QPstart[iter], gates.begin() + QPend[iter] + 1, ycompare);
        }
        
        // find inverse gate; index와 gate id가 달라짐!
        for (int i = 1; i <= num_gates; i++) {
            inverse_gates[gates[i].id] = i;
        }
    }
    ofstream outputFile("final_output.txt");
    outputFile << fixed;
    outputFile.precision(8);
    cout << "Final Gate Positions:" << endl;
    for (int i = 1; i <= num_gates; i++) {
        outputFile << i << " " << gates[inverse_gates[i]].x << " " << gates[inverse_gates[i]].y << endl;
    }
    outputFile.close();
    return 0;
}