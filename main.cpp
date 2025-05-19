#include <iostream>
#include <vector>

#include "Problem.h"

void print_vec(const std::vector<int>& vec);

int main() {
    int n=3,m=2;
    Problem p(n,m);
    p.PZ();
    std::cout<<p.CMax(p.pi)<<std::endl;
    print_vec(p.pi);
    return 0;
}


void print_vec(const std::vector<int>& vec) {
    for (auto &i : vec) {
        std::cout<<i<<" ";
    }
    std::cout<<std::endl;
}