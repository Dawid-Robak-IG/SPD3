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
    // std::string filename = "../tail.dat";
    std::string filename = "../test.dat";
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
    std::chrono::time_point<std::chrono::steady_clock> start0 = std::chrono::steady_clock::now();
    std::vector<std::pair<int,int>> job_sum(n);
    for (int i=0;i<n;i++) {
        int total=0;
        for (int j=0;j<m;j++) {
            total += tasks[j].tasks_durations[i];
        }
        job_sum[i] = {i,total};
    }
    std::sort(job_sum.begin(), job_sum.end(), [](auto &a, auto &b) {
        return a.second > b.second;
    });

    std::vector<int> my_seq;
    std::vector<std::vector<int>> F(n + 1, std::vector<int>(m, 0));
    std::vector<std::vector<int>> B(n + 1, std::vector<int>(m, 0));

    for (int k = 0; k < n; k++) {
        int job = job_sum[k].first;
        int best_cmax = INT_MAX;
        std::vector<int> best_seq;

        int current_seq_size = my_seq.size();
        
        std::vector<std::vector<int>> current_F(current_seq_size, std::vector<int>(m, 0));
        std::vector<std::vector<int>> current_B(current_seq_size, std::vector<int>(m, 0));

        // Calculate F for the current my_seq
        for (int i = 0; i < current_seq_size; i++) {
            int task = my_seq[i];
            for (int j = 0; j < m; j++) {
                int t = tasks[j].tasks_durations[task];
                if (i == 0 && j == 0)
                    current_F[i][j] = t;
                else if (i == 0)
                    current_F[i][j] = current_F[i][j - 1] + t;
                else if (j == 0)
                    current_F[i][j] = current_F[i - 1][j] + t;
                else
                    current_F[i][j] = std::max(current_F[i - 1][j], current_F[i][j - 1]) + t;
            }
        }
        // Calculate B for the current my_seq
        for (int i = current_seq_size - 1; i >= 0; i--) {
            int task = my_seq[i];
            for (int j = m - 1; j >= 0; j--) {
                int t = tasks[j].tasks_durations[task];
                if (i == current_seq_size - 1 && j == m - 1)
                    current_B[i][j] = t;
                else if (j == m - 1)
                    current_B[i][j] = current_B[i + 1][j] + t;
                else if (i == current_seq_size - 1)
                    current_B[i][j] = current_B[i][j + 1] + t;
                else
                    current_B[i][j] = std::max(current_B[i + 1][j], current_B[i][j + 1]) + t;
            }
        }

        for (int pos = 0; pos <= current_seq_size; pos++) {
            std::vector<int> temp = my_seq;
            temp.insert(temp.begin()+pos, job);

            int b_task;
            std::vector<int> C(m, 0);

            for (int m_idx = 0; m_idx < m; m_idx++) {
                if (pos > 0) b_task = current_F[pos-1][m_idx]; // Use current_F
                else b_task = 0;

                int t = tasks[m_idx].tasks_durations[job];

                if (m_idx == 0)
                    C[m_idx] = b_task + t;
                else {
                    C[m_idx] = std::max(b_task, C[m_idx-1]) + t;
                }
            }
            int final_cmax = C[m - 1];
            if (pos < current_seq_size) final_cmax += current_B[pos][m - 1]; // Use current_B

            if (final_cmax < best_cmax) {
                best_cmax = final_cmax;
                best_seq = temp;
            }
        }
        my_seq = best_seq;
    }
    pi = my_seq;

    std::chrono::time_point<std::chrono::steady_clock> end0 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end0 - start0;
    fneh_time = elapsed_seconds.count();
}

void Problem::BnB() {
    std::chrono::time_point<std::chrono::steady_clock> start0 = std::chrono::steady_clock::now();

    NEH(); // Calc UB
    int UB = CMax(pi);
    //int UB = INT_MAX;

    std::vector<int> allTasks(n);
    iota(allTasks.begin(), allTasks.end(), 0);  // Fill with 0, 1, ..., n-1

    // Priority queue sorted by increasing lower bound
    std::priority_queue<Node, std::vector<Node>, CompareNode> Q;
    Q.emplace(std::vector<int>{}, allTasks, 0, 0);  // Start from empty schedule

    while (!Q.empty()) {
        Node node = Q.top(); Q.pop();

        // Skip branch
        if (node.lb >= UB)
            continue;

        // Check if all tasks are scheduled
        if (node.remaining.empty()) {
            int cmax = CMax(node.scheduled);
            if (cmax < UB) {
                UB = cmax;
                pi = node.scheduled;
            }
            continue;
        }

        // Branching
        for (int i = 0; i < node.remaining.size(); ++i) {
            std::vector<int> nextScheduled = node.scheduled;
            std::vector<int> nextRemaining = node.remaining;
            int job = nextRemaining[i];

            // Move one job from remaining to scheduled
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

    for (int j = 0; j < m; ++j) {
        int minTime = INT_MAX;

        for (int i : remaining) {
            minTime = std::min(minTime, tasks[j].tasks_durations[i]);
        }

        if (!remaining.empty()) minSum += minTime;
    }

    return cmax + minSum;  // Lower bound = actual time so far + optimistic estimate
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
    std::uniform_int_distribution<> dis(0, n - 1);

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
