#include <iostream>
#include <vector>

#include "Problem.h"

void print_vec(const std::vector<int>& vec);

int main() {
    int n=3,m=2;
    Problem p(n,m);

    std::cout<<"PZ\n";
    p.PZ();
    std::cout<<p.CMax(p.pi)<<std::endl;
    print_vec(p.pi);
    std::cout<<"===============\n";

    std::cout<<"NEH\n";
    p.fill_test1();
    p.NEH();
    std::cout<<p.CMax(p.pi)<<std::endl;
    print_vec(p.pi);
    std::cout<<"===============\n";

    if (m==2) {
        std::cout<<"JOHNSON\n";
        p.fill_test1();
        p.Johnson();
        std::cout<<p.CMax(p.pi)<<std::endl;
        print_vec(p.pi);
        std::cout<<"===============\n";
    }
    std::cout<<"FNEH\n";
    p.fill_test1();
    p.FNEH();
    std::cout<<p.CMax(p.pi)<<std::endl;
    print_vec(p.pi);
    std::cout<<"===============\n";

    std::cout<<"PZ || NEH\n";
    return 0;
}


void print_vec(const std::vector<int>& vec) {
    for (auto &i : vec) {
        std::cout<<i<<" ";
    }
    std::cout<<std::endl;
}