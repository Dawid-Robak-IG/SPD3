//
// Created by drworms on 5/19/25.
//

#ifndef PROBLEM_H
#define PROBLEM_H

#include "Machine.h"

#include <chrono>
#include <random>

class Problem {
    int n; //liczba zadan
    int m; //procesor
    std::vector<Machine> machines;
    std::vector<Machine> back_up;
public:
    double neh_time;
    double fneh_time;
    double pz_time;
    double john_time;
    std::vector<int> pi;

    Problem(const int n, const int m, int max_val,int min_val);
    void fill(int max_val,int min_val);
    void reload();
    void fill_test1();
    void clear();
    void PZ();
    void NEH();
    void Johnson();
    void FNEH();
    int CMax(const std::vector<int> &perm);
};

struct Instance {
    int machines;
    int tasks;
    int max_val;
    int min_val;

    Instance(int m, int t, int maxv, int minv)
        : machines(m), tasks(t), max_val(maxv), min_val(minv) {}
};


#endif //PROBLEM_H
