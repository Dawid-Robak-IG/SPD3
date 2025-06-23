//
// Created by drworms on 5/19/25.
//

#include "Problem.h"

#include <algorithm>
#include <climits>
#include <iostream>

Problem::Problem(const int n, const int m, int max_val,int min_val) {
    if (n==-1) {
        fill_by_file(1);
    } else {
        this->n = 3;
        this->m = 2;
        fill(max_val,min_val);
    }
    back_up = tasks;
}
void Problem::clear() {
    tasks.clear();
    tasks.resize(n, Task(m));
    pi.clear();
    pi.resize(n);
}
void Problem::fill(int max_val,int min_val) {
    clear();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(min_val, max_val);

    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            tasks[j].tasks_durations[i] = dis(gen);
        }
    }
}
void Problem::reload() {
    clear();
    tasks = back_up;
}
void Problem::fill_test1() {
    clear();
    tasks[0].tasks_durations = {54, 83, 15, 71, 77, 36, 53, 38, 27, 87, 76, 91, 14, 29, 12, 77, 32, 87, 68, 94};
    tasks[1].tasks_durations = {79, 3, 11, 99, 56, 70, 99, 60, 5, 56, 3, 61, 73, 75, 47, 14, 21, 86, 5, 77};
    tasks[2].tasks_durations = {16, 89, 49, 15, 89, 45, 60, 23, 57, 64, 7, 1, 63, 41, 63, 47, 26, 75, 77, 40};
    tasks[3].tasks_durations = {66, 58, 31, 68, 78, 91, 13, 59, 49, 85, 85, 9, 39, 41, 56, 40, 54, 77, 51, 31};
    tasks[4].tasks_durations = {58, 56, 20, 85, 53, 35, 53, 41, 69, 13, 86, 72, 8, 49, 47, 87, 58, 18, 68, 28};

}
void Problem::fill_by_file(int nr_instances) {
    std::string filename = "../tail.dat";
    // std::string filename = "../test.dat";
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << "\n";
    }
    int instances;
    file >> instances;
    std::cout << instances << "\n";
    if (nr_instances>0)instances = nr_instances;
    for (int i = 0; i < instances; ++i) {
        file >> this->n;
        file >> this->m;
        clear();
        std::cout << this->n << "(n) | " << this->m << "(m)\n";
        for (int task_idx=0;task_idx<this->n;task_idx++) {
            for (int times=0;times<this->m;times++) {
                int time,idx;
                file >> idx;
                file >> time;
                tasks[task_idx].tasks_durations[idx] = time;
                std::cout << idx << " " << tasks[task_idx].tasks_durations[idx] << " | ";
            }
            std::cout << "\n";
        }
    }
    file.close();
}

void Problem::PZ() {
    std::chrono::time_point<std::chrono::steady_clock> start0 = std::chrono::steady_clock::now();
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

    std::chrono::time_point<std::chrono::steady_clock> end0 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end0 - start0;
    pz_time = elapsed_seconds.count();
}

void Problem::NEH() {
    std::chrono::time_point<std::chrono::steady_clock> start0 = std::chrono::steady_clock::now();

    std::vector<std::pair<int,int>> job_sum(n);
    for (int i=0;i<n;i++) {
        int total=0;
        for (int j=0;j<m;j++) {
            total += tasks[i].tasks_durations[j];
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

    std::chrono::time_point<std::chrono::steady_clock> end0 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end0 - start0;
    neh_time = elapsed_seconds.count();
}

int Problem::CMax(const std::vector<int> &perm) {
    int curr_n = perm.size();
    std::vector<std::vector<int>> C(curr_n, std::vector<int>(m, 0));
    for (int i = 0; i < curr_n; i++) {
        int job = perm[i];
        for (int j = 0; j < m; j++) {
            int time = tasks[job].tasks_durations[j];
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

void Problem::Johnson() {
    std::chrono::time_point<std::chrono::steady_clock> start0 = std::chrono::steady_clock::now();
    if (m!=2) {
        std::cout<<"Johnson cannot work for m != 2\n";
        return;
    }
    std::vector<int> left;
    std::vector<int> right;
    for (int i=0;i<n;i++) {
        if (tasks[i].tasks_durations[0] >= tasks[i].tasks_durations[1]) {
            right.emplace_back(i);
        } else {
            left.emplace_back(i);
        }
    }
    std::sort(left.begin(), left.end(), [&](auto &a, auto &b) {
        return tasks[a].tasks_durations[0] < tasks[b].tasks_durations[0];
    });
    std::sort(right.begin(), right.end(), [&](auto &a, auto &b) {
        return tasks[a].tasks_durations[1] > tasks[b].tasks_durations[1];
    });

    pi.clear();
    pi.insert(pi.end(), left.begin(), left.end());
    pi.insert(pi.end(), right.begin(), right.end());

    std::chrono::time_point<std::chrono::steady_clock> end0 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end0 - start0;
    john_time = elapsed_seconds.count();
}

void Problem::FNEH() {
    auto start = std::chrono::steady_clock::now();
    std::vector<std::pair<int, int>> job_sum(n);
    for (int i = 0; i < n; ++i) {
        int sum = 0;
        for (int j = 0; j < m; ++j) {
            sum += tasks[i].tasks_durations[j];
        }
        job_sum[i] = {i, sum};
    }
    std::sort(job_sum.begin(), job_sum.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    std::vector<int> sequence;
    std::vector<std::vector<int>> C;
    for (int k = 0; k < n; ++k) {
        int job_to_insert = job_sum[k].first;
        int best_cmax = INT_MAX;
        std::vector<int> best_seq;

        std::vector<std::vector<int>> F_matrix;
        if (!sequence.empty()) {
            F_matrix.resize(sequence.size(), std::vector<int>(m, 0));
            for (size_t i = 0; i < sequence.size(); ++i) {
                int job = sequence[i];
                for (int j = 0; j < m; ++j) {
                    if (i == 0 && j == 0)
                        F_matrix[i][j] = tasks[job].tasks_durations[j];
                    else if (i == 0)
                        F_matrix[i][j] = F_matrix[i][j - 1] + tasks[job].tasks_durations[j];
                    else if (j == 0)
                        F_matrix[i][j] = F_matrix[i - 1][j] + tasks[job].tasks_durations[j];
                    else
                        F_matrix[i][j] = std::max(F_matrix[i - 1][j], F_matrix[i][j - 1]) + tasks[job].tasks_durations[j];
                }
            }
        }

        std::vector<std::vector<int>> B_matrix;
        if (!sequence.empty()) {
            B_matrix.resize(sequence.size(), std::vector<int>(m, 0));
            for (int j = m - 1; j >= 0; --j) {
                if (j == m - 1)
                    B_matrix[sequence.size() - 1][j] = tasks[sequence.back()].tasks_durations[j];
                else
                    B_matrix[sequence.size() - 1][j] = B_matrix[sequence.size() - 1][j + 1] + tasks[sequence.back()].tasks_durations[j];
            }
            for (int i_rev = 1; i_rev < sequence.size(); ++i_rev) {
                int i = sequence.size() - 1 - i_rev;
                int job = sequence[i];
                for (int j = m - 1; j >= 0; --j) {
                    if (j == m - 1)
                        B_matrix[i][j] = B_matrix[i + 1][j] + tasks[job].tasks_durations[j];
                    else
                        B_matrix[i][j] = std::max(B_matrix[i + 1][j], B_matrix[i][j + 1]) + tasks[job].tasks_durations[j];
                }
            }
        }

        for (int pos = 0; pos <= sequence.size(); ++pos) {
            std::vector<int> temp_prefix_end_times(m, 0);
            std::vector<int> temp_suffix_start_times(m, 0);

            if (pos > 0) {
                temp_prefix_end_times = F_matrix[pos - 1];
            }
            if (pos < sequence.size()) {
                temp_suffix_start_times = B_matrix[pos];
            }
            int cmax = computeCmaxQNEH(temp_prefix_end_times, temp_suffix_start_times, job_to_insert, sequence.size() + 1);
            if (cmax < best_cmax) {
                best_cmax = cmax;
                best_seq = sequence;
                best_seq.insert(best_seq.begin() + pos, job_to_insert);
            }
        }
        sequence = best_seq;
    }

    pi = sequence;
    auto end = std::chrono::steady_clock::now();
    fneh_time = std::chrono::duration<double>(end - start).count();
}
int Problem::computeCmaxQNEH(const std::vector<int>& sequence_prefix_end_times,
                              const std::vector<int>& sequence_suffix_start_times,
                              int job_to_insert,
                              int k_length) {
    std::vector<int> current_job_start_times(m);
    std::vector<int> current_job_end_times(m);

    for (int j = 0; j < m; ++j) {
        int prev_task_end_on_this_machine,current_task_end_on_prev_machine;
        if (k_length>0) {
            prev_task_end_on_this_machine = sequence_prefix_end_times[j];
        } else {
            prev_task_end_on_this_machine = 0;
        }
        if (j>0) {
            current_task_end_on_prev_machine = current_job_end_times[j - 1];
        } else {
            current_task_end_on_prev_machine = 0;
        }

        current_job_start_times[j] = std::max(prev_task_end_on_this_machine, current_task_end_on_prev_machine);
        current_job_end_times[j] = current_job_start_times[j] + tasks[job_to_insert].tasks_durations[j];
    }
    int cmax = 0;
    for (int j = 0; j < m; ++j) {
        cmax = std::max(cmax, current_job_end_times[j] + sequence_suffix_start_times[j]);
    }
    return cmax;
}


void Problem::BnB() {
    std::chrono::time_point<std::chrono::steady_clock> start0 = std::chrono::steady_clock::now();
    NEH();
    int UB = CMax(pi);
    std::vector<int> allTasks(n);
    iota(allTasks.begin(), allTasks.end(), 0);

    std::priority_queue<Node, std::vector<Node>, CompareNode> Q;
    Q.emplace(std::vector<int>{}, allTasks, 0, 0);

    while (!Q.empty()) {
        Node node = Q.top(); Q.pop();

        if (node.lb >= UB)
            continue;

        if (node.remaining.empty()) {
            int cmax = CMax(node.scheduled);
            if (cmax < UB) {
                UB = cmax;
                pi = node.scheduled;
            }
            continue;
        }

        for (int i = 0; i < node.remaining.size(); ++i) {
            std::vector<int> nextScheduled = node.scheduled;
            std::vector<int> nextRemaining = node.remaining;
            int job = nextRemaining[i];

            nextScheduled.push_back(job);
            nextRemaining.erase(nextRemaining.begin() + i);

            int cmax = CMax(nextScheduled);
            int LB = calLB(nextScheduled, nextRemaining);

            if (LB < UB) {
                Q.emplace(nextScheduled, nextRemaining, LB, cmax);
            }
        }
    }

    std::chrono::time_point<std::chrono::steady_clock> end0 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end0 - start0;
    bnb_time = elapsed_seconds.count();
}
int Problem::calLB(const std::vector<int> &scheduled, const std::vector<int> &remaining) {
    int cmax = CMax(scheduled);
    int minSum = 0;

    for (int machine = 0; machine < m; ++machine) {
        int minTime = INT_MAX;
        for (int job : remaining) {
            minTime = std::min(minTime, tasks[job].tasks_durations[machine]);
        }
        if (!remaining.empty()) minSum += minTime;
    }

    return cmax + minSum;
}

void Problem::SimulatedAnnealing(double T0, double T_end, int maxIter) {
    std::chrono::time_point<std::chrono::steady_clock> start0 = std::chrono::steady_clock::now();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    NEH();
    //iota(pi.begin(), pi.end(), 0);

    int bestCost = CMax(pi);
    double T = T0;
    double lambda = std::pow(T_end / T0, 1.0 / maxIter);

    for (int iter = 0; iter < maxIter; ++iter) {
        std::vector<int> neighbor = pi;
        ChangePerm(neighbor); //mutacja inplace

        int neighborCost = CMax(neighbor);
        int delta = neighborCost - bestCost;

        if (delta < 0 || dis(gen) < std::exp(-delta / T)) {
            pi = std::move(neighbor);
            bestCost = neighborCost;
        }

        T *= lambda;
    }

    std::chrono::time_point<std::chrono::steady_clock> end0 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end0 - start0;
    sa_time = elapsed_seconds.count();
}
void Problem::ChangePerm(std::vector<int>& perm) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, m - 1);

    int i = dis(gen);
    int j = dis(gen);
    while (j == i) j = dis(gen);
    std::swap(perm[i], perm[j]);
}
void Problem::ThresholdAccepting(double T0, double T_end, int steps_per_threshold, int max_outer_iter) {
    std::chrono::time_point<std::chrono::steady_clock> start0 = std::chrono::steady_clock::now();

    NEH();
    int bestCost = CMax(pi);
    double T = T0;
    double threshold = T0;
    double decay = (T0 - T_end) / max_outer_iter;

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int outer = 0; outer < max_outer_iter; ++outer) {
        for (int inner = 0; inner < steps_per_threshold; ++inner) {
            std::vector<int> neighbor = pi;
            ChangePerm(neighbor);

            int neighborCost = CMax(neighbor);
            int delta = neighborCost - bestCost;

            if (delta <= threshold) {
                pi = std::move(neighbor);
                bestCost = neighborCost;
            }
        }
        threshold -= decay;
        if (threshold < 0) break;
    }

    std::chrono::time_point<std::chrono::steady_clock> end0 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end0 - start0;
    ta_time = elapsed_seconds.count();
}
