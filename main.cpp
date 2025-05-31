#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "Problem.h"

std::vector<Instance> instances = {
    {2,10,10,1}
    //{2,10,20,10},
    //{2,20,20,10},
    //{2,20,100,50},
    // {2,30,10,1},
    // {2,30,100,50},
    // {2,50,10,1},
    // {2,50,100,50},
    // {2,100,100,50},
    // {3,10,10,1},
    // {5,10,10,1},
    // {10,10,10,1},
    // {3,40,70,30},
    // {5,40,70,30},
    // {10,40,70,30},
    // {30,40,70,30},
};

void print_vec(const std::vector<int>& vec);
void test1();
void test_machines(int machines,int tasks, int max_val, int min_val);
void test_machines_csv(const Instance& inst, const std::string& filename, int rep);


int main() {
    for (int i=0;i<instances.size();i++) {
        test_machines_csv(instances[i], "results.csv", 3);
    }
    return 0;
}

void test_machines(int machines,int tasks, int max_val, int min_val) {
    int common_seed = 1;
    Problem p(tasks,machines,max_val,min_val);
    std::cout<<"========================================\n";
    std::cout<<"Machines = "<<machines<<"\n";
    std::cout<<"Tasks = "<<tasks<<"\n";
    std::cout<<"Max Val = "<<max_val<<"\n";
    std::cout<<"Min Val = "<<min_val<<"\n";
    std::cout<<"kryt || time\n";
    if (machines==2)std::cout<<"PZ || NEH || JOHN || FNEH\n";
    else std::cout<<"PZ || NEH || FNEH\n";

    /*
    p.PZ();
    std::cout<<p.CMax(p.pi)<<";"<<p.pz_time<<";";

    p.reload();
    p.NEH();
    std::cout<<p.CMax(p.pi)<<";"<<p.neh_time<<";";

    if (machines==2) {
        p.reload();
        p.Johnson();
        std::cout<<p.CMax(p.pi)<<";"<<p.john_time<<";";
    }

    p.reload();
    p.FNEH();
    std::cout<<p.CMax(p.pi)<<";"<<p.fneh_time<<";";

    p.reload();
    p.BnB();
    std::cout<<p.CMax(p.pi)<<";BnB: "<<p.bnb_time<<";";
    */
    p.reload();
    p.SimulatedAnnealing(1000.0, 0.1, 10000);
    std::cout<<p.CMax(p.pi)<<";SA: "<<p.sa_time<<";";

    std::cout<<"\n========================================\n";
}

void test_machines_csv(const Instance& inst, const std::string& filename, int rep) {
    std::ofstream file(filename, std::ios::app); // append mode
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << "\n";
        return;
    }

    // CSV Header
    file << "algorithm,machines,tasks,max_val,min_val,cmax,time\n";

    for (int i = 0; i < rep; ++i) {
        Problem p(inst.tasks, inst.machines, inst.max_val, inst.min_val);

        // PZ
        p.reload();
        p.PZ();
        file << "PZ," << inst.machines << "," << inst.tasks << "," << inst.max_val << "," << inst.min_val << ","
             << p.CMax(p.pi) << "," << p.pz_time << "\n";

        // NEH
        p.reload();
        p.NEH();
        file << "NEH," << inst.machines << "," << inst.tasks << "," << inst.max_val << "," << inst.min_val << ","
             << p.CMax(p.pi) << "," << p.neh_time << "\n";

        // JOHNSON (only if 2 machines)
        if (inst.machines == 2) {
            p.reload();
            p.Johnson();
            file << "JOHN," << inst.machines << "," << inst.tasks << "," << inst.max_val << "," << inst.min_val << ","
                 << p.CMax(p.pi) << "," << p.john_time << "\n";
        }

        // FNEH
        p.reload();
        p.FNEH();
        file << "FNEH," << inst.machines << "," << inst.tasks << "," << inst.max_val << "," << inst.min_val << ","
             << p.CMax(p.pi) << "," << p.fneh_time << "\n";

        // BnB
        p.reload();
        p.BnB();
        file << "BnB," << inst.machines << "," << inst.tasks << "," << inst.max_val << "," << inst.min_val << ","
             << p.CMax(p.pi) << "," << p.bnb_time << "\n";

        /*
        // Simulated Annealing
        p.reload();
        p.SimulatedAnnealing(1000.0, 0.1, 10000);
        file << "SA," << inst.machines << "," << inst.tasks << "," << inst.max_val << "," << inst.min_val << ","
             << p.CMax(p.pi) << "," << p.sa_time << "\n";
             */
    }

    file.close();
}

void test1() {
    int n=3,m=2;
    Problem p(n,m,10,1);

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
}
void print_vec(const std::vector<int>& vec) {
    for (auto &i : vec) {
        std::cout<<i<<" ";
    }
    std::cout<<std::endl;
}