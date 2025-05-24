//
// Created by drworms on 5/19/25.
//

#include "Machine.h"


Machine::Machine(const int n) {
    this->n = n;
    tasks_durations.resize(n);
}
