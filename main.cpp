#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "Problem.h"

std::vector<Instance> instances = {
    // {20,5,10,1}
    {2,10,10,1},
    {2,15,10,1},
    {2,50,10,5},
    {2,100,10,5},
    {2,200,10,5},
    {2,300,10,5},
{2,500,10,5}
    // {3,5,10,1},
    // {3,10,10,1},
    // {3,15,10,1}
};

std::vector<MetaInstance> meta_instances = {
    {1000.0,0.1,100},
    {1000.0,0.1,1000},
    {1000.0,0.1,10000},
    //======================================
    {100.0,0.1,1000},
    {1000.0,0.1,1000},
    {10000.0,0.1,1000},
    //======================================
    {1000.0,0.1,1000},
    {1000.0,1,1000},
    {1000.0,10,1000}
};

void print_vec(const std::vector<int>& vec);
void test1();
void test_machines(int machines,int tasks, int max_val, int min_val);
void test_machines_csv(const Instance& inst, const std::string& filename, int rep);
void test_meta(const MetaInstance& meta_inst, const std::string& filename, int rep);
void test_taylor();
void test_taylor_all();
void test_johnson(const Instance& inst, const std::string& filename, int rep);

int main() {
    std::string filename = "results.csv";
    std::ofstream file(filename, std::ios::app); // append mod
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << "\n";
    }
    file << "PZ, NEH, JOHN, FNEH, BnB, Ann, Thres\n";
    file.close();
    // test_taylor();
    // test_taylor_all();




    for (int i=0;i<instances.size();i++) {
        test_johnson(instances[i], "results.csv", 1);
    }
    // for (int i=0;i<meta_instances.size();i++) {
    //     std::cout<<"Going to do "<<i<<std::endl;
    //     test_meta(meta_instances[i], "results.csv", 3);
    // }
    return 0;
}



void test_johnson(const Instance& inst, const std::string& filename, int rep) {
    std::ofstream file(filename, std::ios::app); // append mode
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << "\n";
        return;
    }
    std::cout<<"Opened file\n";

    for (int i = 0; i < rep; ++i) {
        std::cout<<"GOING FOR REP\n";
        Problem p(inst.tasks, inst.machines, inst.max_val, inst.min_val);
        std::cout<<"created problem\n";
        file<<p.n<<","<<p.m<<",";
        if (inst.machines == 2) {
            std::cout<<"Going for Johnson \n";
            p.reload();
            std::cout<<"reloaded\n";
            p.Johnson();
            file << p.CMax(p.pi) << "," << p.john_time << ",";
        } else {
            file << "---,---,";
        }
    }
    file  << "\n";
    file.close();
    std::cout<<"Closed file\n";
}
void test_taylor_all() {
    auto instances = Problem::load_all_instances("../tail.dat");

    std::string filename = "results.csv";
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << "\n";
        return;
    }

    for (size_t i = 0; i < instances.size(); ++i) {
        Problem p(instances[i]);

        std::cout << "Processing instance " << i+1 << "/" << instances.size()
                  << " (n=" << p.n << ", m=" << p.m << ")\n";
        file<<p.n<<","<<p.m<<",";

        // if (p.n > 12)
        //     file << "---,---,";
        // else {
        //     p.reload();
        //     p.PZ();
        //     file << p.CMax(p.pi) << "," << p.pz_time << ",";
        // }

        p.reload();
        p.NEH();
        file << p.CMax(p.pi) << "," << p.neh_time << ",";

        // if (p.m == 2) {
        //     p.reload();
        //     p.Johnson();
        //     file << p.CMax(p.pi) << "," << p.john_time << ",";
        // } else
        //     file << "---,---,";

        p.reload();
        p.FNEH();
        file << p.CMax(p.pi) << "," << p.fneh_time << ",";

        // if (p.n > 12)
        //     file << "---,---,";
        // else {
        //     p.reload();
        //     p.BnB();
        //     file << p.CMax(p.pi) << "," << p.bnb_time << ",";
        // }

        p.reload();
        p.SimulatedAnnealing(1000.0, 0.1, 10000);
        file << p.CMax(p.pi) << "," << p.sa_time << ",";

        p.reload();
        p.ThresholdAccepting(1000.0, 0.1, 100, 100);
        file << p.CMax(p.pi) << "," << p.ta_time << "\n";
    }

    file.close();
}
void test_taylor() {
    Problem p;

    std::string filename = "results.csv";
    std::ofstream file(filename, std::ios::app); // append mode
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << "\n";
        return;
    }
    std::cout<<"Opened file\n";

    int rep = 1;

    for (int i = 0; i < rep; ++i) {
        // PZ
        if (p.n>12) {
            std::cout<<"NOT GOING FOR PZ\n";
            file << "---,---,";
        }
        else {
            std::cout<<"Going for PZ \n";
            p.reload();
            p.PZ();
            file << p.CMax(p.pi) << "," << p.pz_time << ",";
        }

        // NEH
        std::cout<<"Going for NEH \n";
        p.reload();
        p.NEH();
        file << p.CMax(p.pi) << "," << p.neh_time << ",";

        //johnson
        if (p.m == 2) {
            std::cout<<"Going for Johnson \n";
            p.reload();
            p.Johnson();
            file << p.CMax(p.pi) << "," << p.john_time << ",";
        } else {
            file << "---,---,";
        }

        // FNEH
        std::cout<<"Going for FNEH \n";
        p.reload();
        p.FNEH();
        file << p.CMax(p.pi) << "," << p.fneh_time << ",";

        // BnB
        if (p.n>12) {
            file << "---,---,";
        }
        else {
            std::cout<<"Going for BnB \n";
            p.reload();
            p.BnB();
            file <<p.CMax(p.pi) << "," << p.bnb_time << ",";
        }


        // Simulated Annealing
        std::cout<<"Going for sim Annealing \n";
        p.reload();
        p.SimulatedAnnealing(1000.0, 0.1, 10000);
        file << p.CMax(p.pi) << "," << p.sa_time << ",";

        std::cout<<"Going for threshold accepting \n";
        p.reload();
        p.ThresholdAccepting(1000.0,0.1,100,1000);
        file << p.CMax(p.pi) << "," << p.ta_time << "\n";
    }
    file.close();
    std::cout<<"Closed file\n";
}
void test_machines_csv(const Instance& inst, const std::string& filename, int rep) {
    std::ofstream file(filename, std::ios::app); // append mode
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << "\n";
        return;
    }
    std::cout<<"Opened file\n";
    // CSV Header
    // file << "algorithm,machines,tasks,max_val,min_val,cmax,time\n";

    for (int i = 0; i < rep; ++i) {
        std::cout<<"GOING FOR REP\n";
        Problem p(inst.tasks, inst.machines, inst.max_val, inst.min_val);
        file<<p.n<<","<<p.m<<",";

        // PZ
        // if (inst.tasks>1) {
        //     std::cout<<"NOT GOING FOR PZ\n";
        //     file << "---,---,";
        // }
        // else {
        //     std::cout<<"Going for PZ \n";
        //     p.reload();
        //     p.PZ();
        //     // file << "PZ," << inst.machines << "," << inst.tasks << "," << inst.max_val << "," << inst.min_val << ","
        //     //      << p.CMax(p.pi) << "," << p.pz_time << "\n";
        //     file << p.CMax(p.pi) << "," << p.pz_time << ",";
        // }

        // NEH
        // std::cout<<"Going for NEH \n";
        // p.reload();
        // p.NEH();
        // // file << "NEH," << inst.machines << "," << inst.tasks << "," << inst.max_val << "," << inst.min_val << ","
        // //      << p.CMax(p.pi) << "," << p.neh_time << "\n";
        // file << p.CMax(p.pi) << "," << p.neh_time << ",";
        //
        // // JOHNSON (only if 2 machines)
        // if (inst.machines == 2) {
        //     std::cout<<"Going for Johnson \n";
        //     p.reload();
        //     p.Johnson();
        //     file << p.CMax(p.pi) << "," << p.john_time << ",";
        // } else {
        //     file << "---,---,";
        // }

        // FNEH
        // std::cout<<"Going for FNEH \n";
        // p.reload();
        // p.FNEH();
        // file << p.CMax(p.pi) << "," << p.fneh_time << ",";

        // BnB
        // if (inst.tasks>1) {
        //     file << "---,---,";
        // }
        // else {
        //     std::cout<<"Going for BnB \n";
        //     p.reload();
        //     p.BnB();
        //     file <<p.CMax(p.pi) << "," << p.bnb_time << ",";
        // }


        // Simulated Annealing
        // std::cout<<"Going for sim Annealing \n";
        // p.reload();
        // p.SimulatedAnnealing(1000.0, 0.1, 10000);
        // file << p.CMax(p.pi) << "," << p.sa_time << ",";

        // std::cout<<"Going for threshold accepting \n";
        // p.reload();
        // p.ThresholdAccepting(1000.0,0.1,100,10000);
        // file << p.CMax(p.pi) << "," << p.ta_time << "\n";
    }
    file  << "\n";
    file.close();
    std::cout<<"Closed file\n";
}
void test_meta(const MetaInstance& meta_inst, const std::string& filename, int rep) {
    std::ofstream file(filename, std::ios::app); // append mode
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << "\n";
        return;
    }
    std::cout<<"Opened file\n";

    Instance inst = {3,8,10,1};

    int anneal_c_med = 0;
    double anneal_time_med = 0;
    int thres_c_med = 0;
    double thres_time_med = 0;

    for (int i = 0; i < rep; ++i) {
        Problem p(inst.tasks, inst.machines, inst.max_val, inst.min_val);

        // Simulated Annealing
        p.reload();
        p.SimulatedAnnealing(meta_inst.T0, meta_inst.Tend, meta_inst.max_iter);
        anneal_c_med += p.CMax(p.pi);
        anneal_time_med += p.sa_time;

        p.reload();
        p.ThresholdAccepting(meta_inst.T0, meta_inst.Tend,100,meta_inst.max_iter);
        thres_c_med += p.CMax(p.pi);
        thres_time_med += p.ta_time;
    }
    file << anneal_c_med/3 <<","<<anneal_time_med/3<<",";
    file << thres_c_med/3 <<","<<thres_time_med/3<<"\n";
    file.close();
    std::cout<<"Closed file\n";
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
