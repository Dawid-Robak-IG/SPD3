//
// Created by drworms on 5/19/25.
//

#include "Problem.h"

#include <algorithm>
#include <climits>
#include <iostream>

Problem::Problem(const int n, const int m, int max_val,int min_val) {
    this->n = n;
    this->m = m;
    // fill_test1();
    fill(max_val,min_val);
    back_up = machines;
}
void Problem::clear() {
    machines.clear();
    machines.resize(m, Machine(n));
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
            machines[j].tasks_durations[i] = dis(gen);
        }
    }
}
void Problem::reload() {
    clear();
    machines = back_up;
}
void Problem::fill_test1() {
    clear();
    machines[0].tasks_durations = {3,2,4};
    machines[1].tasks_durations = {2,1,3};
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
            total += machines[j].tasks_durations[i];
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
            int time = machines[j].tasks_durations[job];
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
        if (machines[0].tasks_durations[i] >= machines[1].tasks_durations[i]) {
            right.emplace_back(i);
        } else {
            left.emplace_back(i);
        }
    }
    std::sort(left.begin(), left.end(), [&](auto &a, auto &b) {
        return machines[0].tasks_durations[a] < machines[1].tasks_durations[b];
    });
    std::sort(right.begin(), right.end(), [&](auto &a, auto &b) {
        return machines[0].tasks_durations[a] > machines[1].tasks_durations[b];
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
            total += machines[j].tasks_durations[i];
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

        int size = my_seq.size();

        for (int i = 0; i < size; i++) { //F
            int task = my_seq[i];
            for (int j = 0; j < m; j++) {
                int t = machines[j].tasks_durations[task];
                if (i == 0 && j == 0)
                    F[i][j] = t;
                else if (i == 0)
                    F[i][j] = F[i][j - 1] + t;
                else if (j == 0)
                    F[i][j] = F[i - 1][j] + t;
                else
                    F[i][j] = std::max(F[i - 1][j], F[i][j - 1]) + t;
            }
        }
        for (int i = size - 1; i >= 0; i--) { //B
            int task = my_seq[i];
            for (int j = m - 1; j >= 0; j--) {
                int t = machines[j].tasks_durations[task];
                if (i == size - 1 && j == m - 1)
                    B[i][j] = t;
                else if (i == size - 1)
                    B[i][j] = B[i][j + 1] + t;
                else if (j == m - 1)
                    B[i][j] = B[i + 1][j] + t;
                else
                    B[i][j] = std::max(B[i + 1][j], B[i][j + 1]) + t;
            }
        }

        for (int pos = 0; pos <= size; pos++) {
            std::vector<int> temp = my_seq;
            temp.insert(temp.begin()+pos, job);

            int b_task;
            std::vector<int> C(m, 0);

            for (int m_idx = 0; m_idx < m; m_idx++) {
                if (pos>0)b_task = F[pos-1][m_idx];
                else b_task = 0;

                int t = machines[m_idx].tasks_durations[job];

                if (m_idx == 0)
                    C[m_idx] = b_task + t;
                else {
                    C[m_idx] = std::max(b_task, C[m_idx-1]) + t;
                }
            }
            int final_cmax = C[m - 1];
            if (pos < size) final_cmax += B[pos][m - 1];

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
            minTime = std::min(minTime, machines[j].tasks_durations[i]);
        }

        if (!remaining.empty()) minSum += minTime;
    }

    return cmax + minSum;  // Lower bound = actual time so far + optimistic estimate
}