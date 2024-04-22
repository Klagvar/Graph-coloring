#include "util.h"

#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <random>
#include <iostream>

static void swap_uint(unsigned int &a, unsigned int &b) {
  std::swap(a, b);
}

static void swap_int(int &a, int &b) {
  std::swap(a, b);
}

/* Utility function that puts all
non-positive (0 and negative) numbers on left
side of arr[] and return count of such numbers */
static int segregate(std::vector<int> &arr) {
  int j = 0;
  for (size_t i = 0; i < arr.size(); i++) {
    if (arr[i] <= 0) {
      swap_int(arr[i], arr[j]);
      j++;  // increment count of non-positive integers
    }
  }

  return j;
}

/* Find the smallest positive missing number
in an array that contains all positive integers */
static int findMissingPositive(std::vector<int>::iterator begin, std::vector<int>::iterator end) {
  int size = std::distance(begin, end);
  for (auto it = begin; it != end; ++it) {
    if (std::abs(*it) - 1 < size && *(begin + std::abs(*it) - 1) > 0)
      *(begin + std::abs(*it) - 1) = -*(begin + std::abs(*it) - 1);
  }

  // Return the first index value at which is positive
  for (int i = 0; i < size; i++)
    if (*(begin + i) > 0)
      // 1 is added because indexes start from 0
      return i + 1;

  return size + 1;
}



static void heapify(std::vector<unsigned int> &keys, std::vector<unsigned int> &values,
                    unsigned int n, int i) {
  int max = i;  // Initialize max as root
  unsigned int leftChild = 2 * i + 1;
  unsigned int rightChild = 2 * i + 2;

  // If left child is greater than root
  if (leftChild < n && keys[leftChild] > keys[max]) max = leftChild;

  // If right child is greater than max
  if (rightChild < n && keys[rightChild] > keys[max]) max = rightChild;

  // If max is not root
  if (max != i) {
    swap_uint(keys[i], keys[max]);
    swap_uint(values[i], values[max]);
    // heapify the affected sub-tree recursively
    heapify(keys, values, n, max);
  }
}


struct sort_item {
  unsigned int key;
  unsigned int value;
};

static int compare_keys_stable(const sort_item &a, const sort_item &b) {
  int diff = a.key - b.key;
  return diff != 0 ? diff : a.value - b.value;
}

/* EXPOSED FUNCTIONS */

unsigned int UTIL_smallest_missing_number(std::vector<int> &arr) {
  int shift = segregate(arr);
  return findMissingPositive(arr.begin() + shift, arr.end());
}


double UTIL_get_time() {
  auto now = std::chrono::system_clock::now();
  auto duration = now.time_since_epoch();
  return std::chrono::duration_cast<std::chrono::microseconds>(duration).count() * 1e-6;
}

void UTIL_randomize_array(std::vector<unsigned int> &arr) {
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(arr.begin(), arr.end(), g);
}

void UTIL_print_array(const std::vector<unsigned int> &arr) {
  for (const auto &i : arr) std::cout << i << " ";
  std::cout << "\n";
}

void UTIL_heapsort_values_by_keys(std::vector<unsigned int> &keys, std::vector<unsigned int> &values) {
  // сортировка кучей по возрастанию
  // Rearrange array (building heap)
  for (int i = keys.size() / 2 - 1; i >= 0; i--) heapify(keys, values, keys.size(), i);

  // Extract elements from heap one by one
  for (int i = keys.size() - 1; i >= 0; i--) {
    swap_uint(keys[0], keys[i]);  // Current root moved to the end
    swap_uint(values[0], values[i]);
    heapify(keys, values, i, 0);  // calling max heapify on the heap reduced
  }
}

unsigned int UTIL_max_in_array(const std::vector<unsigned int> &arr) {
  return *std::max_element(arr.begin(), arr.end());
}

void UTIL_stable_qsort_values_by_keys(std::vector<unsigned int> &keys, std::vector<unsigned int> &values) {
  std::vector<sort_item> a(keys.size());

  for (size_t i = 0; i < keys.size(); i++) {
    a[i].key = keys[i];
    a[i].value = values[i];
  }

  std::stable_sort(a.begin(), a.end(), compare_keys_stable);

  for (size_t i = 0; i < keys.size(); i++) {
    keys[i] = a[i].key;
    values[i] = a[i].value;
  }
}