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
    double y;
    double x;
    vector<int> connected_nets;
} gate;

typedef struct _net{
    int id;
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

/*
typedef struct _fake_pad{
    int og_id;
    vector<int> connected_nets;
    double y;
    double x;
} fake_pad;
*/

int xcompare(gate a, gate b) {
    if (a.x < b.x) return 1;
    if (a.x > b.x) return 0;
    if (a.y < b.y) return 1;
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
    
    // QP1
    {
        coo_matrix A;
        valarray<double> b_x(0.0, num_gates);
        valarray<double> b_y(0.0, num_gates);
        
        vector<int> A_rows, A_cols;
        vector<double> A_vals;

        for (int i = 1; i <= num_gates; i++) {
            // 대각성분 계산
            int total_weight = 0;
            for (int j : gates[i].connected_nets) {
                if (use_clique_weight) total_weight += 1;
                else total_weight += nets[j].total_gates_pads - 1;
            }
            A_rows.push_back(i-1);
            A_cols.push_back(i-1);
            A_vals.push_back(total_weight);
            
            
            for (int j : gates[i].connected_nets) {
                // nondiagonal 계산
                for (int k : nets[j].connected_gates) {
                    if (i != k) {
                        A_rows.push_back(i-1);
                        A_cols.push_back(k-1);
                        if (use_clique_weight) A_vals.push_back(-1/(nets[j].total_gates_pads - 1)); // weight = 1 / k-1 (using clique method)
                        else A_vals.push_back(-1);
                    }
                }

                // vector b 계산
                for (int k : nets[j].connected_pads) {
                    if (use_clique_weight) {
                        b_x[i - 1] += pads[k].x / (nets[j].total_gates_pads - 1);
                        b_y[i - 1] += pads[k].y / (nets[j].total_gates_pads - 1);
                    }
                    else {
                        b_x[i - 1] += pads[k].x;
                        b_y[i - 1] += pads[k].y;
                    }
                }
            }
        }
        A.n = num_gates;
        A.nnz = A_vals.size();
        A.row.resize(A.nnz);
        A.col.resize(A.nnz);
        A.dat.resize(A.nnz);
        for (int i = 0; i < A.nnz; i++) {
            A.row[i] = A_rows[i];
            A.col[i] = A_cols[i];
            A.dat[i] = A_vals[i];
        }
        valarray<double> x_coords(0.0, num_gates);
        valarray<double> y_coords(0.0, num_gates);
        double min_x = 0;
        double min_y = 0;
        double max_x = SIZE;
        double max_y = SIZE;

        A.solve(b_x, x_coords);
        A.solve(b_y, y_coords);
        // min이나 max 넘어가는 경우 처리. epsilon은 10^-6
        for (int i = 1; i <= num_gates; i++) {
            gates[i].x = x_coords[i-1];
            if (abs(x_coords[i-1] - min_x) < epsilon) gates[i].x = min_x;
            else if (abs(x_coords[i-1] - max_x) < epsilon) gates[i].x = max_x;
            gates[i].y = y_coords[i-1];
            if (abs(y_coords[i-1] - min_y) < epsilon) gates[i].y = min_y;
            else if (abs(y_coords[i-1] - max_y) < epsilon) gates[i].y = max_y;
            cout << i << " " << gates[i].x << " " << gates[i].y << endl;
        }
    }
    // dividing. 
    /*
    vector<int> left_gate(num_gates / 2 + 1); // hmmmmmmmmmmm index 붙여서 QP23, 4567, ... 이렇게 한번에 처리할수는 없음? size는 알 수 있지않나?
    vector<int> right_gate(num_gates / 2 + 1 + (num_gates % 2));
    vector<pair<int, gate>> arr(num_gates + 1);
    for (int i = 1; i <= num_gates; i++) {
        arr[i] = {i, gates[i]};
    }
    sort(arr.begin() + 1, arr.end(), xcompare);
    cout << "left" << endl;
    for (int i = 1; i <= num_gates / 2; i++) {
        cout << arr[i].first << " " << arr[i].second.x << " " << arr[i].second.y << endl;
        left_gate[i] = arr[i].first;
    }
    cout << "right" << endl;
    for (int i = num_gates / 2 + 1; i <= num_gates; i++) {
        cout << arr[i].first << " " << arr[i].second.x << " " << arr[i].second.y << endl;
        right_gate[i] = arr[i].first;
    }
    */
    // 이거 정렬까지 QP1 넣어야할듯?
    sort(gates.begin()+1, gates.end(), xcompare);
    cout << "left" << endl;
    for (int i = 1; i <= num_gates / 2; i++) {
        cout << gates[i].id << " " << gates[i].x << " " << gates[i].y << endl;
    }
    cout << "right" << endl;
    for (int i = num_gates / 2 + 1; i <= num_gates; i++) {
        cout << gates[i].id << " " << gates[i].x << " " << gates[i].y << endl;
    }
    for (int i = 1; i <= num_gates; i++) {
        for (int j = 1; j <= num_gates; j++) {
            if (gates[j].id == i) {
                inverse_gates[i] = j;
            }
        }
    }
    
    cout << "QP2 start" << endl;
    // QP2
    {
        coo_matrix A;
        int num_gates_QP2 = num_gates / 2;
        int QP2_start = 1;
        int QP2_end = num_gates / 2;
        valarray<double> b_x(0.0, num_gates_QP2);
        valarray<double> b_y(0.0, num_gates_QP2);
        
        vector<int> A_rows, A_cols;
        vector<double> A_vals;

        for (int i = 0; i < num_gates_QP2; i++) { // i는 index!!
            // 대각성분 계산
            int total_weight = 0;
            for (int j : gates[QP2_start + i].connected_nets) {
                if (use_clique_weight) total_weight += 1; // 연결된 gate나 pad의 수는 동일
                else total_weight += nets[j].total_gates_pads - 1;
            }
            A_rows.push_back(i);
            A_cols.push_back(i);
            A_vals.push_back(total_weight);
            
            for (int j : gates[QP2_start + i].connected_nets) {
                // nondiagonal 계산
                for (int k : nets[j].connected_gates) { // k는 gate id, i와 다름
                    if (gates[QP2_start + i].id != k) {
                        if (inverse_gates[k] >= QP2_start && inverse_gates[k] <= QP2_end) { // index가 범위 안에 있음, 즉 k는 QP2 gate임. 
                            A_rows.push_back(i);
                            A_cols.push_back(inverse_gates[k]-QP2_start);
                            if (use_clique_weight) A_vals.push_back(-1/(nets[j].total_gates_pads - 1)); // weight = 1 / k-1 (using clique method)
                            else A_vals.push_back(-1);
                        }
                        else { // index가 범위 밖에 있음. 즉 k는 QP3 gate이므로 pad로 만들어 vector b에 추가. fake pad가 필요한가? 단순히 값만 가져오면 되는거아님? 아님 말고
                            if (use_clique_weight) {
                                b_x[i] += ((gates[k].x > 50) ? 50 : gates[k].x) / (nets[j].total_gates_pads - 1);
                                b_y[i] += gates[k].y / (nets[j].total_gates_pads - 1);
                            }
                            else {
                                b_x[i] += (gates[k].x > 50) ? 50 : gates[k].x;
                                b_y[i] += gates[k].y;
                            }
                        }
                    }
                }

                // vector b 계산
                for (int k : nets[j].connected_pads) { // pad는 앞과 변화 X
                    if (use_clique_weight) {
                        b_x[i] += ((pads[k].x > 50) ? 50 : pads[k].x) / (nets[j].total_gates_pads - 1);
                        b_y[i] += pads[k].y / (nets[j].total_gates_pads - 1);
                    }
                    else {
                        b_x[i] += (pads[k].x > 50) ? 50 : pads[k].x;
                        b_y[i] += pads[k].y;
                    }
                }
            }
        }
        A.n = num_gates_QP2;
        A.nnz = A_vals.size();
        A.row.resize(A.nnz);
        A.col.resize(A.nnz);
        A.dat.resize(A.nnz);
        for (int i = 0; i < A.nnz; i++) {
            A.row[i] = A_rows[i];
            A.col[i] = A_cols[i];
            A.dat[i] = A_vals[i];
        }
        valarray<double> x_coords(0.0, num_gates_QP2);
        valarray<double> y_coords(0.0, num_gates_QP2);
        double min_x = 0;
        double min_y = 0;
        double max_x = SIZE / 2;
        double max_y = SIZE;

        A.solve(b_x, x_coords);
        A.solve(b_y, y_coords);
        for (int i = 0; i < num_gates_QP2; i++) {
            gates[i + QP2_start].x = x_coords[i];
            if (abs(x_coords[i] - min_x) < epsilon) gates[i + QP2_start].x = min_x;
            else if (abs(x_coords[i] - max_x) < epsilon) gates[i + QP2_start].x = max_x;
            gates[i + QP2_start].y = y_coords[i];
            if (abs(y_coords[i] - min_y) < epsilon) gates[i + QP2_start].y = min_y;
            else if (abs(y_coords[i] - max_y) < epsilon) gates[i + QP2_start].y = max_y;
            cout << gates[i + QP2_start].id << " " << gates[i + QP2_start].x << " " << gates[i + QP2_start].y << endl;
        }
    }

    // QP3
    cout << "QP3 start" << endl;
    {
        coo_matrix A;
        int num_gates_QP3 = num_gates / 2 + (num_gates % 2);
        int QP3_start = num_gates / 2 + 1;
        int QP3_end = num_gates;
        valarray<double> b_x(0.0, num_gates_QP3);
        valarray<double> b_y(0.0, num_gates_QP3);
        
        vector<int> A_rows, A_cols;
        vector<double> A_vals;

        for (int i = 0; i < num_gates_QP3; i++) {
            // 대각성분 계산
            int total_weight = 0;
            for (int j : gates[QP3_start + i].connected_nets) {
                if (use_clique_weight) total_weight += 1; // 연결된 gate나 pad의 수는 동일
                else total_weight += nets[j].total_gates_pads - 1;
            }
            A_rows.push_back(i);
            A_cols.push_back(i);
            A_vals.push_back(total_weight);
            
            for (int j : gates[QP3_start + i].connected_nets) {
                // nondiagonal 계산
                for (int k : nets[j].connected_gates) { // k는 gate id, i와 다름
                    if (gates[QP3_start + i].id != k) {
                        if (inverse_gates[k] >= QP3_start && inverse_gates[k] <= QP3_end) { // index가 범위 안에 있음, 즉 k는 QP2 gate임. 
                            A_rows.push_back(i);
                            A_cols.push_back(inverse_gates[k]-QP3_start);
                            if (use_clique_weight) A_vals.push_back(-1/(nets[j].total_gates_pads - 1)); // weight = 1 / k-1 (using clique method)
                            else A_vals.push_back(-1);
                        }
                        else { // index가 범위 밖에 있음. 즉 k는 QP3 gate이므로 pad로 만들어 vector b에 추가. fake pad가 필요한가? 단순히 값만 가져오면 되는거아님? 아님 말고
                            if (use_clique_weight) {
                                b_x[i] += ((gates[k].x > 50) ? gates[k].x : 50) / (nets[j].total_gates_pads - 1);
                                b_y[i] += gates[k].y / (nets[j].total_gates_pads - 1);
                            }
                            else {
                                b_x[i] += (gates[k].x > 50) ? gates[k].x : 50;
                                b_y[i] += gates[k].y;
                            }
                        }
                    }
                }

                // vector b 계산
                for (int k : nets[j].connected_pads) { // pad는 앞과 변화 X
                    if (use_clique_weight) {
                        b_x[i] += ((pads[k].x > 50) ? pads[k].x : 50) / (nets[j].total_gates_pads - 1);
                        b_y[i] += pads[k].y / (nets[j].total_gates_pads - 1);
                    }
                    else {
                        b_x[i] += (pads[k].x > 50) ? pads[k].x : 50;
                        b_y[i] += pads[k].y;
                    }
                }
            }
        }
        A.n = num_gates_QP3;
        A.nnz = A_vals.size();
        A.row.resize(A.nnz);
        A.col.resize(A.nnz);
        A.dat.resize(A.nnz);
        for (int i = 0; i < A.nnz; i++) {
            A.row[i] = A_rows[i];
            A.col[i] = A_cols[i];
            A.dat[i] = A_vals[i];
        }
        valarray<double> x_coords(0.0, num_gates_QP3);
        valarray<double> y_coords(0.0, num_gates_QP3);
        double min_x = SIZE / 2;
        double min_y = 0;
        double max_x = SIZE;
        double max_y = SIZE;

        A.solve(b_x, x_coords);
        A.solve(b_y, y_coords);
        for (int i = 0; i < num_gates_QP3; i++) {
            gates[i + QP3_start].x = x_coords[i];
            if (abs(x_coords[i] - min_x) < epsilon) gates[i + QP3_start].x = min_x;
            else if (abs(x_coords[i] - max_x) < epsilon) gates[i + QP3_start].x = max_x;
            gates[i + QP3_start].y = y_coords[i];
            if (abs(y_coords[i] - min_y) < epsilon) gates[i + QP3_start].y = min_y;
            else if (abs(y_coords[i] - max_y) < epsilon) gates[i + QP3_start].y = max_y;
            cout << gates[i + QP3_start].id << " " << gates[i + QP3_start].x << " " << gates[i + QP3_start].y << endl;
        }
    }
    return 0;
}