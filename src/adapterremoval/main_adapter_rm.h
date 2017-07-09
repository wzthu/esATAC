#ifndef MAIN_ADAPTER_RM_H
#define MAIN_ADAPTER_RM_H

#include <cstring>
#include<string>
#include"scheduler.h"
#include "userconfig.h"
namespace ar
{
void add_write_step(const userconfig& config, scheduler& sch, size_t offset,
                    const std::string& name, analytical_step* step);
int remove_adapter_sequences(const userconfig& config);
 } // namespace ar
#endif
