//
// Created by drworms on 5/19/25.
//

#include "Problem.h"

#include <algorithm>
#include <climits>

Problem::Problem(const int n, const int m) {
    this->n = n;
    this->m = m;
    clear();
    fill_test1();
}
void Problem::clear() {
    tasks.clear();
    tasks.resize(m, Machine(n));
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

int Problem::CMax(const std::vector<int> &perm) {
    std::vector<std::vector<int>> C(n, std::vector<int>(m, 0));
    for (int i = 0; i < n; i++) {
        int job = perm[i];
        for (int j = 0; j < m; j++) {
            int time = tasks[j].operations[job];
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

    return C[n - 1][m - 1];
}
void Problem::fill_test1() {
    tasks[0].operations = {3,2,4};
    tasks[1].operations = {2,1,3};
}