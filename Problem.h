//
// Created by drworms on 5/19/25.
//

#ifndef PROBLEM_H
#define PROBLEM_H

#include "Machine.h"

class Problem {
    int n; //liczba zadan
    int m; //procesor
    std::vector<Machine> tasks;
public:
    std::vector<int> pi;

    Problem(const int n, const int m);
    void fill_test1();
    void clear();
    void PZ();
    int CMax(const std::vector<int> &perm);
};



#endif //PROBLEM_H
