//
// Created by drworms on 5/19/25.
//

#include "Task.h"


Task::Task(const int m) {
    this->m = m;
    tasks_durations.resize(m);
}
