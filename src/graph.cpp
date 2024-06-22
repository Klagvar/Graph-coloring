#include "graph.h"
#include <vector>
#include <string>
#include <memory>
#include <mutex>
#include <fstream>
#include <sstream>
#include <thread>
#include <iostream>
#include <climits>

#define N_COLORING_METHODS 6
enum ColoringMethod {
  seq_greedy,
  seq_ldf,
  rec_rlf,
  par_jp,
  par_opt_ver,
  par_opt_ed
};

const std::vector<std::string> coloring_methods = {
    "seq_greedy", "seq_ldf", "rec_rlf", "par_jp", "par_opt_ver", "par_opt_ed"};

struct Node {
  int index;
  std::shared_ptr<Node> next;
};

struct Graph {
  unsigned int V, E;
  std::vector<std::shared_ptr<Node>> ladj; 
  std::shared_ptr<Node> z; 
  std::vector<unsigned int> degree;  // количество соседей у вершины
  std::vector<unsigned int> color;
};

struct ParamStruct {
  std::shared_ptr<Graph> G;
  unsigned int index;
  std::vector<unsigned int> weights;
  std::vector<unsigned int> vertexes;
  std::mutex* lock;
  unsigned int n_threads;

  ParamStruct(std::mutex* lock) : lock(lock) {}

  ParamStruct(std::shared_ptr<Graph> G, unsigned int n_threads, std::vector<unsigned int> weights, std::mutex* lock)
    : G(G), index(0), weights(weights), vertexes(std::vector<unsigned int>(G->V)), lock(lock), n_threads(n_threads) {}
};

/* UTILITY FUNCTIONS */
static ColoringMethod method_str_to_enum(std::string coloring_method_str) {
  for (unsigned int i = 0; i < N_COLORING_METHODS; i++) {
    if (coloring_method_str == coloring_methods[i]) {
      return static_cast<ColoringMethod>(i);
    }
  }
  return static_cast<ColoringMethod>(-1);
}

static std::shared_ptr<Node> LINK_new(int index, std::shared_ptr<Node> next) {
  auto x = std::make_shared<Node>();
  x->index = index;
  x->next = next;
  return x;
}

static std::shared_ptr<Node> EDGE_insert(std::shared_ptr<Graph> G, int from, int to) {
  auto new_node = LINK_new(to, G->ladj[from]);
  G->ladj[from] = new_node;
  G->degree[from]++;
  G->E++;
  return G->ladj[from];
}

/* Алгоритмы расскраски */
/* Последовательный жадный */
std::vector<unsigned int> color_sequential_greedy(std::shared_ptr<Graph> G) {
  unsigned int n = G->V; 
  std::vector<unsigned int> random_order(n); 
  for (unsigned int i = 0; i < n; i++) { // Инициализируем вектора случайного порядка и вектор цветов.
    G->color[i] = 0;
    random_order[i] = i;
  }
  UTIL_randomize_array(random_order); // Случайно перемешаем вектор
  for (unsigned int i = 0; i < n; i++) {
    unsigned int ii = random_order[i];

    std::vector<int> neighbours_colors(G->degree[ii]); // Вектор цветов каждого соседа
    unsigned int j = 0;
    for (auto t = G->ladj[ii]; t != G->z; t = t->next) { // Инициализируем вектор цветов каждого соседа
      neighbours_colors[j++] = G->color[t->index];
    }

    G->color[ii] = UTIL_smallest_missing_number(neighbours_colors); // Находим минимальный доступный цвет
  }
  return G->color;
}

/* Последовательный LDF */
std::vector<unsigned int> color_sequential_ldf(std::shared_ptr<Graph> G) {
  unsigned int n = G->V; 
  std::vector<unsigned int> degree(n); 
  std::vector<unsigned int> vertex(n); 
  for (unsigned int i = 0; i < n; i++) { // Заполнение векторов степеней, номеров вершин и цветов 0 для всех вершин
    G->color[i] = 0;
    degree[i] = G->degree[i];
    vertex[i] = i;
  }
  UTIL_heapsort_values_by_keys(degree, vertex); // Сортировка кучей по возрастанию
  for (int i = n - 1; i >= 0; i--) {
    unsigned int ii = vertex[i];  // Сортировка кучей сортирует степени в порядке возрастания,
                                  // поэтому мы получаем к ним доступ в обратном порядке.
    std::vector<int> neighbours_colors(degree[i]); // Выделение памяти для вектора цветов соседей
    unsigned int j = 0;
    for (auto t = G->ladj[ii]; t != G->z; t = t->next) { // Заполнение вектора цветов соседей
      neighbours_colors[j++] = G->color[t->index];
    }

    G->color[ii] = UTIL_smallest_missing_number(neighbours_colors); // Заполнение цвета вершины c использованием минимального доступного цвета
  }
  return G->color;
}

/* Последовательный RLF */
std::vector<unsigned int> color_rlf(std::shared_ptr<Graph> G) {
  unsigned int n = G->V; 
  std::vector<unsigned int> order(n); 
  std::vector<unsigned int> degree_copy(n); 
  for (unsigned int i = 0; i < n; i++) { // Инициализируем вектор порядка, вектор цветов и копию степеней вершин.
    G->color[i] = 0;
    order[i] = i;
    degree_copy[i] = G->degree[i];
  }

  unsigned int color = 1; 
  while (1) {
    int max_degree_index = -1;
    for (unsigned int i = 0; i < n; i++) { // Ищем нераскрашенную вершину с максимальной степенью
      if (G->color[order[i]] == 0 && (max_degree_index == -1 || degree_copy[order[i]] > degree_copy[order[max_degree_index]])) {
        max_degree_index = i;
      }
    }
    if (max_degree_index == -1) { // Если все вершины раскрашены, выходим из цикла
      break;
    }

    std::vector<unsigned int> S(n); 
    unsigned int S_size = 0;
    S[S_size++] = order[max_degree_index]; // Добавляем вершину с максимальной степенью в S
    G->color[order[max_degree_index]] = color; 
    using link = std::shared_ptr<Node>;
    for (unsigned int i = 0; i < n; i++) {
      if (G->color[order[i]] == 0) { // Если вершина нераскрашена
        unsigned int j;
        for (j = 0; j < S_size; j++) { // Проверяем, смежна ли она с вершинами из S
          link t;
          for (t = G->ladj[order[i]]; t != G->z && static_cast<unsigned int>(t->index) != S[j]; t = t->next);
          if (t != G->z) {
            break;
          }
        }
        if (j == S_size) { // Если вершина не смежна с вершинами из S, добавляем её в S и раскрашиваем
          S[S_size++] = order[i];
          //printf("Расскрашиваем вершину %d в цвет %d\n", order[i] + 1, color);
          G->color[order[i]] = color;
        }
      }
    }

    for (unsigned int i = 0; i < n; i++) { // Обновляем степени вершин
      if (G->color[order[i]] == color) {
        for (auto t = G->ladj[order[i]]; t != G->z; t = t->next) {
          if (G->color[t->index] == 0) {
            degree_copy[t->index]--;
          }
        }
      }
    }

    color++; 
  }

  return G->color;
}

/* Параллельная расскраска Джонса-Плассмана */
void jp_color_vertex(std::shared_ptr<Graph> G, unsigned int index, unsigned int n_threads,
                     std::vector<unsigned int> &weights) {
  unsigned int n = G->V; 
  int uncolored = n / n_threads; // Количество вершин, которые должныы быть обработанны в этом потоке
  if (n % n_threads && (index < n % n_threads)) {
    uncolored++;  // это нужно для управления графами с количеством вершин
                  // не кратным количеству потоков
  }
  while (uncolored > 0) {
    for (unsigned int i = index; i < n; i += n_threads) {
      if (G->color[i] == 0) { // если вершина не раскрашена
        std::vector<int> neighbours_colors(G->degree[i]); 
        unsigned int has_highest_number = 1;
        unsigned int j = 0;
        for (auto t = G->ladj[i]; t != G->z; t = t->next) { // Заполнение вектора цветов соседей
          if (G->color[t->index] == 0 &&
              (weights[t->index] > weights[i] ||
               (weights[t->index] == weights[i] && (unsigned int)t->index > i))) {
            has_highest_number = 0;
            break;
          } else {
            neighbours_colors[j++] = G->color[t->index];
          }
        }

        if (has_highest_number) { // если вершина имеет максимальный вес среди соседей, то раскрашиваем ее
          G->color[i] =
              UTIL_smallest_missing_number(neighbours_colors);
          uncolored--;
        }
      }
    }
  }
}

void jp_color_vertex_wrapper(ParamStruct par) { // обёртка функции jp_color_vertex
  unsigned int index = par.index; // read the index
  par.lock->unlock(); // then unlock the mutex
  jp_color_vertex(par.G, index, par.n_threads, par.weights);
  
}

std::vector<unsigned int> color_parallel_jp(std::shared_ptr<Graph> G, unsigned int n_threads) {
  unsigned int n = G->V;
  std::vector<unsigned int> weights(n); // вектор весов

  for (unsigned int i = 0; i < n; i++) { // Инициализация цветов и заполнение вектора весов
    G->color[i] = 0;
    weights[i] = rand();
  }

  std::mutex mutex; // Позволяет синхронизировать доступ к вектору цветов
  ParamStruct par (G, n_threads, weights, &mutex); // структура для хранения необходимых параметров
  std::vector<std::thread> threads(n_threads); // выделение памяти для вектора потоков

  for (unsigned int i = 0; i < n_threads; i++) { // Инициализация вектора потоков
    par.lock->lock(); // lock the mutex or wait until it is unlocked
    par.index = i; // assign i
    threads[i] = std::thread(jp_color_vertex_wrapper, par); // запуск потока
  }

  for (auto &thread : threads) {
    thread.join();
  }
  return G->color;
}

/* Алгоритмы оптимизации расскраски */
// Раскраска графа для оптимизации расскраски
void color_for_optimize(std::shared_ptr<Graph> G, unsigned int current_colors, std::vector<unsigned int> &mistakes) {
  unsigned int n = G->V;
  for (unsigned int i = 0; i < n; i++) {
    std::vector<unsigned int> color_counts(current_colors, 0);
    for (auto t = G->ladj[i]; t != G->z; t = t->next) {
      if (G->color[t->index] != 0) {
        color_counts[G->color[t->index] - 1]++;
      }
    }
    unsigned int min_color = 0;
    for (unsigned int j = 1; j < current_colors; j++) {
      if (color_counts[j] < color_counts[min_color]) {
        min_color = j;
      }
    }
    G->color[i] = min_color + 1;

    // Обновляем количество ошибок
    mistakes[i] = 0;
    for (auto t = G->ladj[i]; t != G->z; t = t->next) {
      if (G->color[i] == G->color[t->index]) {
        mistakes[i]++;
      }
    }
  }
}

// Перекраска лучшей вершины и возврат количества ошибок после перекраски
unsigned int best_vertex_color(std::shared_ptr<Graph> G, unsigned int current_colors, int best_vertex, std::vector<unsigned int> &mistakes) {
  std::vector<unsigned int> color_counts(current_colors, 0);
  for (auto t = G->ladj[best_vertex]; t != G->z; t = t->next) {
    color_counts[G->color[t->index] - 1]++;
  }

  unsigned int min_color = 0;
  for (unsigned int i = 1; i < current_colors; i++) {
    if (color_counts[i] < color_counts[min_color]) {
      min_color = i;
    }
  }

  G->color[best_vertex] = min_color + 1;
  int new_errors = 0;
  for (auto t = G->ladj[best_vertex]; t != G->z; t = t->next) {
    if (G->color[t->index] == G->color[best_vertex]) {
      new_errors++;
    }
  }

  mistakes[best_vertex] = new_errors;
  for (auto t = G->ladj[best_vertex]; t != G->z; t = t->next) {
    int neighbor = t->index;
    int neighbor_errors = 0;
    for (auto t2 = G->ladj[neighbor]; t2 != G->z; t2 = t2->next) {
      if (G->color[t2->index] == G->color[neighbor]) {
        neighbor_errors++;
      }
    }
    mistakes[neighbor] = neighbor_errors;
  }

  return new_errors;
}

// Подсчёт количества ошибок (соседей с таким-же цветом)
unsigned int count_mist(std::shared_ptr<Graph> G, int v) {
  unsigned int mist = 0;
  for (auto t = G->ladj[v]; t != G->z; t = t->next) {
    if (G->color[v] == G->color[t->index]) {
      mist++;
    }
  }

  return mist;
}

// Функция для поиска лучшей вершины в каждом потоке
void find_best_vertex_thread(std::shared_ptr<Graph> G, unsigned int start, unsigned int end, int &best_vertex, unsigned int &max_mist, std::vector<int> &proc, std::vector<unsigned int> &mistakes, std::mutex &mutex) {
  int local_best_vertex = -1;
  unsigned int local_max_mist = 0;
  for (unsigned int i = start; i < end; i++) {
    if (proc[i]) continue;
    unsigned int mist = mistakes[i];
    if (mist > local_max_mist) {
      local_max_mist = mist;
      local_best_vertex = i;
    }
  }

  std::lock_guard<std::mutex> lock(mutex);
  if (local_max_mist > max_mist) {
    max_mist = local_max_mist;
    best_vertex = local_best_vertex;
  }
}

// Измененная функция color_ver_optimize
std::vector<unsigned int> color_ver_optimize_parallel(std::shared_ptr<Graph> G, unsigned int n_threads) {
  unsigned int n = G->V;
  for (unsigned int i = 0; i < n; i++) {
    G->color[i] = 0;
  }

  unsigned int current_colors = 2;
  std::mutex mutex;
  std::vector<unsigned int> mistakes(n, 0);
  do {
    color_for_optimize(G, current_colors, mistakes);
    unsigned int prev_mist = UINT_MAX;
    unsigned int cur_mist = UINT_MAX - 1;
    while(cur_mist < prev_mist) {
      std::vector<int> proc(n, 0);
      prev_mist = cur_mist;
      cur_mist = 0;
      
      while (1) {
        int best_vertex = -1;
        unsigned int max_mist = 0;
        std::vector<std::thread> threads(n_threads);
        unsigned int chunk_size = n / n_threads;
        for (unsigned int i = 0; i < n_threads; i++) {
          unsigned int start = i * chunk_size;
          unsigned int end = (i == n_threads - 1) ? n : start + chunk_size;
          threads[i] = std::thread(find_best_vertex_thread, G, start, end, std::ref(best_vertex), std::ref(max_mist), std::ref(proc), std::ref(mistakes), std::ref(mutex));
        }

        for (auto &thread : threads) {
          thread.join();
        }

        if (best_vertex == -1) break;
        proc[best_vertex] = 1;
        cur_mist += best_vertex_color(G, current_colors, best_vertex, mistakes);
      }
    }
    current_colors++;
  } while (!GRAPH_check_given_coloring_validity(G, G->color));

  return G->color;
}

/* Оптимизация расскраски для ребер */
// Функция для оптимизации раскраски рёбер
void color_for_optimize_edges(std::shared_ptr<Graph> G, unsigned int current_colors, std::vector<std::vector<unsigned int>> &edge_mistakes) {
  unsigned int n = G->V;
  for (unsigned int i = 0; i < n; i++) {
    std::vector<unsigned int> color_counts(current_colors, 0);
    for (auto t = G->ladj[i]; t != G->z; t = t->next) {
      if (G->color[t->index] != 0) {
        color_counts[G->color[t->index] - 1]++;
      }
    }

    unsigned int min_color = 0;
    for (unsigned int j = 1; j < current_colors; j++) {
      if (color_counts[j] < color_counts[min_color]) {
        min_color = j;
      }
    }

    G->color[i] = min_color + 1;

    // Обновляем количество ошибок рёбер
    for (auto t = G->ladj[i]; t != G->z; t = t->next) {
      unsigned int v1 = i;
      unsigned int v2 = t->index;
      edge_mistakes[v1][v2] = count_mist(G, v1) + count_mist(G, v2);
      edge_mistakes[v2][v1] = edge_mistakes[v1][v2];
    }
  }
}

unsigned int best_edge_color(std::shared_ptr<Graph> G, unsigned int current_colors, int best_ed_v1, int best_ed_v2, std::vector<std::vector<unsigned int>> &edge_mistakes) {
  std::vector<unsigned int> color_counts_v1(current_colors, 0);
  std::vector<unsigned int> color_counts_v2(current_colors, 0);
  for (auto t = G->ladj[best_ed_v1]; t != G->z; t = t->next) {
    color_counts_v1[G->color[t->index] - 1]++;
  }

  for (auto t = G->ladj[best_ed_v2]; t != G->z; t = t->next) {
    color_counts_v2[G->color[t->index] - 1]++;
  }

  unsigned int min_color_v1 = 0;
  unsigned int min_color_v2 = 0;
  for (unsigned int i = 1; i < current_colors; i++) {
    if (color_counts_v1[i] < color_counts_v1[min_color_v1]) {
      min_color_v1 = i;
    }
  }

  for (unsigned int i = 1; i < current_colors; i++) {
    if (color_counts_v2[i] < color_counts_v2[min_color_v2]) {
      min_color_v2 = i;
    }
  }

  G->color[best_ed_v1] = min_color_v1 + 1;
  G->color[best_ed_v2] = min_color_v2 + 1;

  int new_errors_v1 = count_mist(G, best_ed_v1);
  int new_errors_v2 = count_mist(G, best_ed_v2);

  edge_mistakes[best_ed_v1][best_ed_v2] = new_errors_v1 + new_errors_v2;
  edge_mistakes[best_ed_v2][best_ed_v1] = edge_mistakes[best_ed_v1][best_ed_v2];

  return edge_mistakes[best_ed_v1][best_ed_v2];
}

// Функция для поиска лучшего ребра в каждом потоке на основе уже подсчитанных ошибок
void find_best_edge_thread(std::shared_ptr<Graph> G, unsigned int start, unsigned int end, int &best_ed_v1, int &best_ed_v2, unsigned int &max_mist, std::vector<std::vector<int>> &proc_edges, std::vector<std::vector<unsigned int>> &edge_mistakes, std::mutex &mutex) {
  int local_best_ed_v1 = -1;
  int local_best_ed_v2 = -1;
  unsigned int local_max_mist = 0;
  for (unsigned int i = start; i < end; i++) {
    for (auto t = G->ladj[i]; t != G->z; t = t->next) {
      if(proc_edges[i][t->index]) continue;
      unsigned int mist = edge_mistakes[i][t->index];
      if (mist == UINT_MAX) continue; // Пропускаем отсутствие ребра
      if (mist > local_max_mist) {
        local_max_mist = mist;
        local_best_ed_v1 = i;
        local_best_ed_v2 = t->index;
      }
    }
  }

  std::lock_guard<std::mutex> lock(mutex);
  if (local_max_mist > max_mist) {
    max_mist = local_max_mist;
    best_ed_v1 = local_best_ed_v1;
    best_ed_v2 = local_best_ed_v2;
  }
}

// Основная функция для параллельной оптимизации раскраски рёбер
std::vector<unsigned int> color_ed_optimize_parallel(std::shared_ptr<Graph> G, unsigned int n_threads) {
  unsigned int n = G->V;
  for (unsigned int i = 0; i < n; i++) {
    G->color[i] = 0;
  }
  unsigned int current_colors = 2;
  std::mutex mutex;

  std::vector<std::vector<unsigned int>> edge_mistakes(n, std::vector<unsigned int>(n, UINT_MAX));

  do {
    color_for_optimize_edges(G, current_colors, edge_mistakes);

    unsigned int prev_mist = UINT_MAX;
    unsigned int cur_mist = UINT_MAX - 1;
    while(cur_mist < prev_mist) {
      std::vector<std::vector<int>> proc_edges(n, std::vector<int>(n, 0));
      prev_mist = cur_mist;
      cur_mist = 0;
      while (1) {
        int best_ed_v1 = -1;
        int best_ed_v2 = -1;
        unsigned int max_mist = 0;
        std::vector<std::thread> threads(n_threads);
        unsigned int chunk_size = n / n_threads;
        for (unsigned int i = 0; i < n_threads; i++) {
          unsigned int start = i * chunk_size;
          unsigned int end = (i == n_threads - 1) ? n : start + chunk_size;
          threads[i] = std::thread(find_best_edge_thread, G, start, end, std::ref(best_ed_v1), std::ref(best_ed_v2), std::ref(max_mist), std::ref(proc_edges), std::ref(edge_mistakes), std::ref(mutex));
        }
        for (auto &thread : threads) {
          thread.join();
        }
        if (best_ed_v1 == -1 || best_ed_v2 == -1) break;
        proc_edges[best_ed_v1][best_ed_v2] = 1;
        proc_edges[best_ed_v2][best_ed_v1] = 1;
        cur_mist += best_edge_color(G, current_colors, best_ed_v1, best_ed_v2, edge_mistakes);
      }
    }
    current_colors++;
  } while (!GRAPH_check_given_coloring_validity(G, G->color));
  return G->color;
}


/* EXPOSED FUNCTIONS */
bool GRAPH_check_given_coloring_validity(std::shared_ptr<Graph> G, std::vector<unsigned int> &colors) {
  for (unsigned int i = 0; i < G->V; i++) {
    for (auto t = G->ladj[i]; t != G->z; t = t->next) {
      if (colors[i] == colors[t->index] || colors[i] == 0) {
        return 0;
      }
    }
  }
  return 1;
}

bool GRAPH_check_current_coloring_validity(std::shared_ptr<Graph> G) {
  for (unsigned int i = 0; i < G->V; i++) {
    for (auto t = G->ladj[i]; t != G->z; t = t->next) {
      if (G->color[i] == G->color[t->index] || G->color[i] == 0) {
        return false;
      }
    }
  }
  return true;
}

unsigned int GRAPH_get_edge_count(std::shared_ptr<Graph> G) { return G->E; }

unsigned int GRAPH_get_vertex_count(std::shared_ptr<Graph> G) { return G->V; }

std::shared_ptr<Graph> GRAPH_load_from_file(std::string filename) {
  std::ifstream fp(filename);
  if (!fp.is_open()) {
    std::cout << "Error opening file " << filename << "\n";
    return nullptr;
  }
  std::string dot = filename.substr(filename.find_last_of("."));

  if (dot == ".graph") {
    unsigned int V, E, fmt = 0, ncon = 0;
    std::string line;
    std::getline(fp, line);
    std::istringstream iss(line);
    iss >> V >> E >> fmt >> ncon;

    auto G = GRAPH_init(V);

    if (G == nullptr) {
      fp.close();
      return nullptr;
    }

    unsigned int from = 1;

    while (std::getline(fp, line)) {
      if (line[0] == '%') {
        continue;
      }
      std::istringstream iss(line);
      unsigned int to;
      int alternate = 0;
      if (fmt == 100) {
        if (ncon != 0) {
          fmt = 10;
        } else {
          fmt = 0;
        }
      }
      switch (fmt) {
        case 10:
          for (unsigned int i = 0; i < ncon; i++) {
            iss >> to;
          }
        case 0:
          while (iss >> to) {
            if (from != to) {
              if (EDGE_insert(G, from - 1, to - 1) == nullptr) {
                std::cout << "Couldn't insert edge from " << from << " to " << to << "\n";
                fp.close();
                GRAPH_free(G);
                return nullptr;
              }
            }
          }
          break;
        case 11:
          iss >> to;
        case 1:
          while (iss >> to) {
            if (alternate % 2 == 0 && from != to) {
              if (EDGE_insert(G, from - 1, to - 1) == nullptr) {
                std::cout << "Couldn't insert edge from " << from << " to " << to << "\n";
                fp.close();
                GRAPH_free(G);
                return nullptr;
              }
            }
            alternate++;
          }
          break;
        default:
          std::cout << "Invalid fmt\n";
          fp.close();
          GRAPH_free(G);
          return nullptr;
      }
      from++;
    }
    fp.close();
    return G;
  } else if (dot == ".gra") {
    int V = -1;
    std::string buf;
    while (V <= 0) {
      fp >> buf;
      if (std::stoi(buf) > 0) {
        V = std::stoi(buf);
      }
    }

    auto G = GRAPH_init(V);
    if (G == nullptr) {
      fp.close();
      return nullptr;
    }
    for (int i = 0; i < V; i++) {
      fp.ignore(std::numeric_limits<std::streamsize>::max(), ':');
      fp >> buf;
      while (buf != "#") {
        int to = std::stoi(buf);
        if (EDGE_insert(G, i, to) == nullptr) {
          std::cout << "Couldn't insert edge from " << i << " to " << to << "\n";
          fp.close();
          GRAPH_free(G);
          return nullptr;
        }
        if (EDGE_insert(G, to, i) == nullptr) {
          std::cout << "Couldn't insert edge from " << to << " to " << i << "\n";
          fp.close();
          GRAPH_free(G);
          return nullptr;
        }
        fp >> buf;
      }
    }
    fp.close();
    return G;
  } else {
    std::cout << "Invalid extension " << dot << "\n";
    fp.close();
    return NULL;
  }
  return NULL;
}

std::shared_ptr<Graph> GRAPH_init(unsigned int V) {
  auto G = std::make_shared<Graph>();
  if (G == nullptr) {
    std::cerr << "Error while allocating the graph\n";
    return nullptr;
  }
  G->V = V;
  G->E = 0;
  G->z = LINK_new(-1, nullptr);
  G->ladj.resize(V); 
  G->degree.resize(V);
  G->color.resize(V);

  for (unsigned int i = 0; i < V; i++) {
    G->ladj[i] = G->z;
    G->degree[i] = 0;
    G->color[i] = 0;
  }
  return G;
}

void GRAPH_free(std::shared_ptr<Graph> G) {
  using link = std::shared_ptr<Node>;
  if (G == nullptr) {
    return;
  }
  link next;
  for (unsigned int i = 0; i < G->V; i++) {
    for (link t = G->ladj[i]; t != G->z; t = next) {
      next = t->next;
      t.reset();
    }
  }
  G->ladj.clear();
  G->degree.clear();
  G->color.clear();
  G->z.reset();
}

void GRAPH_ladj_print(std::shared_ptr<Graph> G) {
  using link = std::shared_ptr<Node>;
  for (unsigned int i = 0; i < G->V; i++) {
    std::cout << i + 1 << " -->";
    for (link t = G->ladj[i]; t != G->z; t = t->next) {
      std::cout << (t == G->ladj[i] ? ' ' : '-') << (t->index) + 1;
    }
    std::cout << "\n";
  }
}

void GRAPH_ladj_print_with_colors(std::shared_ptr<Graph> G, std::vector<unsigned int> &colors) {
  using link = std::shared_ptr<Node>;
  for (unsigned int i = 0; i < G->V; i++) {
    std::cout << i + 1 << "(" << colors[i] << ") -->";
    for (link t = G->ladj[i]; t != G->z; t = t->next) {
      std::cout << (t == G->ladj[i] ? ' ' : '-') << (t->index) + 1 << "(" << colors[t->index] << ")";
    }
    std::cout << "\n";
  }
}

std::vector<unsigned int> GRAPH_get_degrees(std::shared_ptr<Graph> G) { return G->degree; }

unsigned long GRAPH_compute_bytes(std::shared_ptr<Graph> G) {
  using link = std::shared_ptr<Node>;
  unsigned long long bytes = 0;
  bytes += sizeof(G);
  bytes += sizeof(*G);
  bytes += (G->V + G->E + 1) * sizeof(link);
  bytes += (G->V + G->E + 1) * sizeof(Node);
  bytes += (2 * G->V) * sizeof(unsigned int);
  return bytes;
}

std::vector<unsigned int> GRAPH_color(std::shared_ptr<Graph> G, std::string coloring_method_str,
                          unsigned int n_threads) {
  switch (method_str_to_enum(coloring_method_str)) {
    case seq_greedy:
      return color_sequential_greedy(G);
      break;
    case seq_ldf:
      return color_sequential_ldf(G);
      break;
    case rec_rlf:
      return color_rlf(G);
      break;
    case par_jp:
      return color_parallel_jp(G, n_threads);
      break;
    case par_opt_ver:
      return color_ver_optimize_parallel(G, n_threads);
      break;  
    case par_opt_ed:
      return color_ed_optimize_parallel(G, n_threads);
      break;
    default:
      std::cerr << "Passed coloring method '" << coloring_method_str << "' is not valid!\n";
      return std::vector<unsigned int>();
  }
}