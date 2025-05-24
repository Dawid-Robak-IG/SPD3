//
// Created by drworms on 5/19/25.
//

#ifndef PROBLEM_H
#define PROBLEM_H

#include "Machine.h"

#include <chrono>
#include <random>
#include <queue>

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
    double bnb_time;
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

    void BnB();
    int calLB(const std::vector<int> &scheduled, const std::vector<int> &remaining);

};

struct Instance {
    int machines;
    int tasks;
    int max_val;
    int min_val;

    Instance(int m, int t, int maxv, int minv)
        : machines(m), tasks(t), max_val(maxv), min_val(minv) {}
};

struct Node {
    std::vector<int> scheduled;
    std::vector<int> remaining;
    int lb;
    int currentCmax;

    Node(std::vector<int> sched, std::vector<int> rem, int lower, int cmax)
            : scheduled(sched), remaining(rem), lb(lower), currentCmax(cmax) {}
};

struct CompareNode {
    bool operator()(const Node& a, const Node& b) {
        return a.lb > b.lb;
    }
};


#endif //PROBLEM_H
