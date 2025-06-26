//
// Created by drworms on 5/19/25.
//

#ifndef PROBLEM_H
#define PROBLEM_H

#include "Task.h"

#include <chrono>
#include <random>
#include <queue>
#include <fstream>

struct ProblemInstanceData {
    int n, m;
    std::vector<std::vector<int>> durations;
};

class Problem {
public:
    int n; //liczba zadan
    int m; //ilosc maszyn
    std::vector<Task> tasks;
    std::vector<Task> back_up;
    double neh_time;
    double fneh_time;
    double pz_time;
    double john_time;
    double bnb_time;
    double sa_time;
    double ta_time;
    std::vector<int> pi;

    Problem(const int n=-1, const int m=-1, int max_val=-1,int min_val=-1);
    Problem(const ProblemInstanceData& data);
    void fill(int max_val,int min_val);
    void reload();
    void fill_test1();
    void clear();
    void fill_by_file(int nr_instances=-1);


    void PZ();
    void NEH();
    void Johnson();
    void FNEH();
    int CMax(const std::vector<int> &perm);

    void BnB();
    int calLB(const std::vector<int> &scheduled, const std::vector<int> &remaining);

    void SimulatedAnnealing(double T0, double T_end, int maxIter);
    void ChangePerm(std::vector<int>& perm);

    void ThresholdAccepting(double T0, double T_end, int steps_per_threshold, int max_outer_iter);

    int computeCmaxQNEH(const std::vector<int>& sequence_prefix_end_times,
                              const std::vector<int>& sequence_suffix_start_times,
                              int job_to_insert,
                              int k_length);

    static std::vector<ProblemInstanceData> load_all_instances(std::string filename);
};

struct Instance {
    int machines;
    int tasks;
    int max_val;
    int min_val;

    Instance(int m, int n, int maxv, int minv)
        : machines(m), tasks(n), max_val(maxv), min_val(minv) {}
};

struct MetaInstance {
    double T0;
    double Tend;
    int max_iter;

    MetaInstance(double T_0, double T_end, int max_Iter)
        : T0(T_0), Tend(T_end), max_iter(max_Iter) {}
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
