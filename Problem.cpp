//
// Created by drworms on 5/19/25.
//

#include "Problem.h"

#include <algorithm>
#include <climits>

Problem::Problem(const int n, const int m) {
    this->n = n;
    this->m = m;
    fill_test1();
}
void Problem::clear() {
    machines.clear();
    machines.resize(m, Machine(n));
    pi.clear();
    pi.resize(n);
}

void Problem::PZ() {
    int best_cMax = INT_MAX;
    std::vector<int> perm;
    perm.resize(n);
    for (int i = 0; i < n; i++) perm[i] = i;
    do {
        int curr_Cmax = CMax(perm);
        if (curr_Cmax < best_cMax) {
            best_cMax = curr_Cmax;
            pi = perm;
        }
    } while (std::next_permutation(perm.begin(),perm.end()));
}
void Problem::NEH() {
    std::vector<std::pair<int,int>> job_sum(n);
    for (int i=0;i<n;i++) {
        int total=0;
        for (int j=0;j<m;j++) {
            total += machines[j].operations[i];
        }
        job_sum[i] = {i,total};
    }
    std::sort(job_sum.begin(), job_sum.end(), [](auto &a, auto &b) {
        return a.second > b.second;
    });

    std::vector<int> my_seq;
    for (int k=0;k<n;k++) {
        int job = job_sum[k].first;
        int best_cMax = INT_MAX;
        std::vector<int> best_operations;
        for (int i=0;i<=my_seq.size();i++) {
            std::vector<int> temp = my_seq;
            temp.insert(temp.begin()+i,job);
            int new_cmax = CMax(temp);
            if (new_cmax < best_cMax) {
                best_cMax = new_cmax;
                best_operations = temp;
            }
        }
        my_seq = best_operations;
    }
    pi = my_seq;
}


int Problem::CMax(const std::vector<int> &perm) {
    int curr_n = perm.size();
    std::vector<std::vector<int>> C(curr_n, std::vector<int>(m, 0));
    for (int i = 0; i < curr_n; i++) {
        int job = perm[i];
        for (int j = 0; j < m; j++) {
            int time = machines[j].operations[job];
            if (i > 0 && j > 0)
                C[i][j] = std::max(C[i - 1][j], C[i][j - 1]) + time;
            else if (j > 0)
                C[i][j] = C[i][j - 1] + time;
            else if (i > 0)
                C[i][j] = C[i - 1][j] + time;
            else
                C[i][j] = time;
        }
    }
    return C[curr_n - 1][m - 1];
}
void Johnson() {
    
}
void Problem::fill_test1() {
    clear();
    machines[0].operations = {3,2,4};
    machines[1].operations = {2,1,3};

}