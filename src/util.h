#ifndef UTIL_H
#define UTIL_H

#include <ctime>
#include <cstdlib>
#include <sys/resource.h>
#include <sys/time.h>
#include <vector>

double UTIL_get_time();
unsigned int UTIL_smallest_missing_number(std::vector<int> &arr);
void UTIL_randomize_array(std::vector<unsigned int> &arr);
void UTIL_print_array(const std::vector<unsigned int> &arr);
void UTIL_heapsort_values_by_keys(std::vector<unsigned int> &keys, std::vector<unsigned int> &values);
unsigned int UTIL_max_in_array(const std::vector<unsigned int> &arr);
void UTIL_stable_qsort_values_by_keys(std::vector<unsigned int> &degrees, std::vector<unsigned int> &indexes);

#endif
