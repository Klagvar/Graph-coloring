#include "graph.h"

#define N_COLORING_METHODS 6
typedef enum {
  seq_greedy,
  seq_ldf,
  rec_rlf,
  par_jp,
  color_op_ver,
  color_op_ed
} coloring_method;
const char *coloring_methods[N_COLORING_METHODS] = {
    "seq_greedy", "seq_ldf", "rec_rlf", "par_jp", "color_op_ver", "color_op_ed"};

typedef struct node *link;

struct node {
  int index;
  link next;
};

struct graph {
  unsigned int V, E;
  link *ladj; 
  link z; 
  unsigned int *degree;  // количество соседей у вершины
  unsigned int *color;
};

typedef struct param_struct {
  Graph G;
  unsigned int index;
  unsigned int *weights;
  unsigned int *vertexes;
  pthread_mutex_t *lock;
  unsigned int n_threads;
} param_t;

/* UTILITY FUNCTIONS */

static coloring_method method_str_to_enum(char *coloring_method_str) {
  for (unsigned int i = 0; i < N_COLORING_METHODS; i++) {
    if (!strcmp(coloring_method_str, coloring_methods[i])) {
      return i;
    }
  }
  return -1;
}

static link LINK_new(int index, link next) {
  link x = malloc(sizeof *x);
  if (x == NULL) {
    fprintf(stderr, "Error while allocating a link\n");
    return NULL;
  }
  x->index = index;
  x->next = next;
  return x;
}

static link EDGE_insert(Graph G, int from, int to) {
  link new = LINK_new(to, G->ladj[from]);
  if (new == NULL) {
    fprintf(stderr, "Error while inserting an edge\n");
    return NULL;
  }
  G->ladj[from] = new;
  G->degree[from]++;
  G->E++;
  return G->ladj[from];
}

/* Алгоритмы расскраски */

/* Последовательный жадный */
unsigned int *color_sequential_greedy(Graph G) {
  unsigned int n = G->V; 
  unsigned int *random_order = malloc(n * sizeof(unsigned int)); 
  if (random_order == NULL) {
    printf("Ошибка выделения массива случайного порядка!\n");
    return NULL;
  }
  for (unsigned int i = 0; i < n; i++) { // Инициализируем массив случайного порядка и массив цветов.
    G->color[i] = 0;
    random_order[i] = i;
  }
  UTIL_randomize_array(random_order, n); // Случайно перемешаем массив.
  for (int i = 0; i < n; i++) {
    unsigned int ii = random_order[i];

    int *neighbours_colors = malloc(G->degree[ii] * sizeof(int)); // Выделяем память для цветов каждого соседа
    unsigned int j = 0;
    for (link t = G->ladj[ii]; t != G->z; t = t->next) { // Инициализируем массив цветов каждого соседа
      neighbours_colors[j++] = G->color[t->index];
    }

    G->color[ii] = UTIL_smallest_missing_number(neighbours_colors, G->degree[ii]); // Находим минимальный доступный цвет
    free(neighbours_colors);
  }
  free(random_order);
  return G->color;
}

/* Последовательный LDF */
unsigned int *color_sequential_ldf(Graph G) {
  unsigned int n = G->V; 
  unsigned int *degree = malloc(n * sizeof(unsigned int)); 
  if (degree == NULL) {
    printf("Error allocating degrees array!\n");
    return NULL;
  }
  unsigned int *vertex = malloc(n * sizeof(unsigned int)); 
  if (vertex == NULL) {
    printf("Error allocating vertex array!\n");
    return NULL;
  }
  for (unsigned int i = 0; i < n; i++) { // Заполнение массивов степеней, номеров вершин и цветов 0 для всех вершин
    G->color[i] = 0;
    degree[i] = G->degree[i];
    vertex[i] = i;
  }
  UTIL_heapsort_values_by_keys(n, degree, vertex); // Сортировка кучей по возрастанию
  for (int i = n - 1; i >= 0; i--) {
    unsigned int ii = vertex[i];  // Сортировка кучей сортирует степени в порядке возрастания,
                                  // поэтому мы получаем к ним доступ в обратном порядке.
    int *neighbours_colors = malloc(degree[i] * sizeof(int)); // Выделение памяти для массива цветов соседей
    if (neighbours_colors == NULL) {
      printf("Error allocating neighbours_colors array!\n");
    }
    unsigned int j = 0;
    for (link t = G->ladj[ii]; t != G->z; t = t->next) { // Заполнение массива цветов соседей
      neighbours_colors[j++] = G->color[t->index];
    }

    G->color[ii] = UTIL_smallest_missing_number(neighbours_colors, degree[i]); // Заполнение цвета вершины c использованием минимального доступного цвета
    free(neighbours_colors);
  }
  free(vertex);
  free(degree);
  return G->color;
}



/* Последовательный RLF */
unsigned int *color_rlf(Graph G) {
  unsigned int n = G->V; 
  unsigned int *order = malloc(n * sizeof(unsigned int)); 
  unsigned int *degree_copy = malloc(n * sizeof(unsigned int)); 
  if (order == NULL || degree_copy == NULL) {
    printf("Ошибка выделения памяти!\n");
    return NULL;
  }
  for (unsigned int i = 0; i < n; i++) { // Инициализируем массив порядка, массив цветов и копию степеней вершин.
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

    unsigned int *S = malloc(n * sizeof(unsigned int)); 
    unsigned int S_size = 0;
    S[S_size++] = order[max_degree_index]; // Добавляем вершину с максимальной степенью в S
    //printf("\nРасскрашиваем макс вершину %d в цвет %d\n", order[max_degree_index] + 1, color);
    G->color[order[max_degree_index]] = color; 

    for (unsigned int i = 0; i < n; i++) {
      if (G->color[order[i]] == 0) { // Если вершина нераскрашена
        unsigned int j;
        for (j = 0; j < S_size; j++) { // Проверяем, смежна ли она с вершинами из S
          link t;
          for (t = G->ladj[order[i]]; t != G->z && t->index != S[j]; t = t->next);
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

    free(S); 

    for (unsigned int i = 0; i < n; i++) { // Обновляем степени вершин
      if (G->color[order[i]] == color) {
        for (link t = G->ladj[order[i]]; t != G->z; t = t->next) {
          if (G->color[t->index] == 0) {
            degree_copy[t->index]--;
          }
        }
      }
    }

    color++; 
  }

  free(order); 
  free(degree_copy); 
  return G->color;
}

/* Паралльная расскраска Джонса-Плассмана*/
void jp_color_vertex(Graph G, unsigned int index, unsigned int n_threads,
                     unsigned int *weights) {
  unsigned int n = G->V; 
  int uncolored = n / n_threads; // Количество вершин, которые должныы быть обработанны в этом потоке
  if (n % n_threads && (index < n % n_threads)) {
    uncolored++;  // это нужно для управления графами с количеством вершин
                  // не кратным количеству потоков
  }
  while (uncolored > 0) {
    for (int i = index; i < n; i += n_threads) {
      if (G->color[i] == 0) { // если вершина не раскрашена
        int *neighbours_colors = malloc(G->degree[i] * sizeof(int)); 
        unsigned int has_highest_number = 1;
        unsigned int j = 0;
        for (link t = G->ladj[i]; t != G->z; t = t->next) { // Заполнение массива цветов соседей
          if (G->color[t->index] == 0 &&
              (weights[t->index] > weights[i] ||
               (weights[t->index] == weights[i] && t->index > i))) {
            has_highest_number = 0;
            break;
          } else {
            neighbours_colors[j++] = G->color[t->index];
          }
        }

        if (has_highest_number) { // если вершина имеет максимальный вес среди соседей, то раскрашиваем ее
          G->color[i] =
              UTIL_smallest_missing_number(neighbours_colors, G->degree[i]);
          uncolored--;
        }
        free(neighbours_colors);
      }
    }
  }
}

void jp_color_vertex_wrapper(void *par) { // обёртка функции jp_color_vertex
  param_t *tD = (param_t *)par; // Приведение параметра к типу param_t для доступа к полям
  unsigned int index = tD->index; // read the index
  pthread_mutex_unlock(tD->lock); // then unlock the mutex
  jp_color_vertex(tD->G, index, tD->n_threads, tD->weights);
}

unsigned int *color_parallel_jp(Graph G, unsigned int n_threads) {
  unsigned int n = G->V;
  unsigned int *weights = malloc(n * sizeof(unsigned int)); // выделение памяти для массива весов
  if (weights == NULL) {
    printf("Error allocating weights array!\n");
    return NULL;
  }
  for (unsigned int i = 0; i < n; i++) { // Инициализация цветов и заполнение массива весов
    G->color[i] = 0;
    weights[i] = rand();
  }
  pthread_mutex_t mutex; // Позволяет синхронизировать доступ к массиву цветов

  param_t *par = malloc(sizeof(param_t)); // структура для хранения необходимых параметров
  par->G = G;
  par->n_threads = n_threads;
  par->weights = weights;
  par->lock = &mutex;
  pthread_t *threads = malloc(n_threads * sizeof(pthread_t)); // выделение памяти для массива потоков
  /* since we use a single param struct for all threads, we need a mutex to
   to make sure the value of par->index is NOT modified until the thread has
   read it
   поскольку мы используем одну структуру параметров для всех потоков, нам нужен мьютекс, 
   чтобы гарантировать, что значение par -> index НЕ изменяется, пока поток не прочитает его.
  */
  pthread_mutex_init(&mutex, NULL); 

  for (unsigned int i = 0; i < n_threads; i++) { // Инициализация массива потоков
    pthread_mutex_lock(&mutex); // lock the mutex or wait until it is unlocked
    par->index = i; // assign i
    pthread_create(&threads[i], NULL, (void *)&jp_color_vertex_wrapper, (void *)par); // запуск потока
  }

  for (int i = 0; i < n_threads; i++) {
    pthread_join(threads[i], NULL);
  }

  free(par);
  free(weights);
  free(threads);
  pthread_mutex_destroy(&mutex);
  return G->color;
}

/* Алгоритмы оптимизации расскраски */

// Раскраска графа для оптимизации расскраски
void color_for_optimize(Graph G, unsigned int current_colors) {
  unsigned int n = G->V;
  for (unsigned int i = 0; i < n; i++) {
    // Инициализируем массив подсчета цветов
    unsigned int *color_counts = malloc(current_colors * sizeof(unsigned int));
    for (unsigned int j = 0; j < current_colors; j++) {
      color_counts[j] = 0;
    }
    // Подсчитываем количество каждого цвета среди соседей
    for (link t = G->ladj[i]; t != G->z; t = t->next) {
      if (G->color[t->index] != 0) { // Если сосед раскрашен
        color_counts[G->color[t->index] - 1]++;
      }
    }
    // Выбор цвета с наименьшим количеством повторений
    unsigned int min_color = 0;
    for (unsigned int j = 1; j < current_colors; j++) {
      if (color_counts[j] < color_counts[min_color]) {
        min_color = j;
      }
    }
    // Раскрашиваем вершину в оптимальный цвет
    G->color[i] = min_color + 1;
    free(color_counts);
  }
}

// Перекраска лучшей вершины
void best_vertex_color(Graph G, unsigned int current_colors, int best_vertex) {
  // Инициализируем массив подсчета цветов
  unsigned int *color_counts = malloc(current_colors * sizeof(unsigned int));
  for (unsigned int i = 0; i < current_colors; i++) {
    color_counts[i] = 0;
  }
  // Подсчитываем количество каждого цвета среди соседей
  for (link t = G->ladj[best_vertex]; t != G->z; t = t->next) {
    color_counts[G->color[t->index] - 1]++;
  }
  // Выбор цвета с наименьшим количеством повторений
  unsigned int min_color = 0;
  for (unsigned int i = 1; i < current_colors; i++) {
    if (color_counts[i] < color_counts[min_color]) {
      min_color = i;
    }
  }
  // Раскрашиваем вершину в оптимальный цвет
  G->color[best_vertex] = min_color + 1;
  free(color_counts);
}

// Подсчёт количества ошибок (соседей с таким-же цветом)
unsigned int count_mist(Graph G, int v) {
  unsigned int mist = 0;
  for (link t = G->ladj[v]; t != G->z; t = t->next) {
    if (G->color[v] == G->color[t->index]) {
      mist++;
    }
  }
  return mist;
}

/* Оптимизация расскраски для вершин */
unsigned int *color_ver_optimize(Graph G) {
  unsigned int n = G->V;
  for (unsigned int i = 0; i < n; i++) {
    G->color[i] = 0;
  }
  unsigned int current_colors = 2;
  // Пока расскраска не будет валидна
  do
  {
    // Расскрашиваю граф в текущее количество цветов на подобии с жадной расскраской
    color_for_optimize(G, current_colors);
    
    int *proc = malloc(n * sizeof(int)); // Массив для отслеживания уже обработанных вершин
    for (unsigned int i = 0; i < n; i++) {
      proc[i] = 0;
    }

    while(1)
    {
      int best_vertex = -1;
      unsigned int max_mist = 0;
      // Выбираем вершину с максимальным количеством соседей с таким-же цветом
      for (unsigned int i = 0; i < n; i++) {
        // Если вершина уже была перекрашена пропускаем
        if (proc[i]) continue; 

        unsigned int mist = count_mist(G, i);
        if (mist > max_mist) {
          max_mist = mist;
          best_vertex = i;
        }
      }
      // Если не нашли лучшую вершину для перекраски выходим из цикла
      if (best_vertex == -1) break; 

      //Отмечаем вершиу как обработанную и перекрашиваем
      proc[best_vertex] = 1;
      best_vertex_color(G, current_colors, best_vertex);
    }
    free(proc);
    current_colors++;
  } while (!GRAPH_check_given_coloring_validity(G, G->color));

  return G->color;
}


/* Оптимизация расскраски для ребер */
// Работает ооооооооооооооооооооооооооооочень долго
unsigned int *color_ed_optimize(Graph G) {
  unsigned int n = G->V;
  for (unsigned int i = 0; i < n; i++) {
    G->color[i] = 0;
  }
  unsigned int current_colors = 2;
  // Пока расскраска не будет валидна
  do
  {
    // Расскрашиваю граф в текущее количество цветов на подобии с жадной расскраской
    color_for_optimize(G, current_colors);
    
    // Квадратный массив для поиска уже обработанных рёбер
    int **proc_edges = malloc(n * sizeof(int *));
    for (unsigned int i = 0; i < n; i++) {
      proc_edges[i] = malloc(n * sizeof(int));
      for (unsigned int j = 0; j < n; j++) {
        proc_edges[i][j] = 0;
      }
    }

    while(1)
    {
      int best_ed_v1 = -1;
      int best_ed_v2 = -1;
      unsigned int max_mist = 0;

      // Выбираем ребро с максимальным количеством соседей с таким-же цветом
      for (unsigned int i = 0; i < n; i++) {
        // Подсчёт ошибок для первой вершины
        unsigned int mist_v1 = count_mist(G, i);

        //Подсчёт количества соседей с таким же цветом для каждого соседа. И поиск лучшего ребра
        for (link t = G->ladj[i]; t != G->z; t = t->next) {
          // Если ребро уже было обработано пропускаем
          if(proc_edges[i][t->index]) continue;
          unsigned int mist_v2 = count_mist(G, t->index);
          if (mist_v1 + mist_v2 > max_mist) {
            max_mist = mist_v1 + mist_v2;
            best_ed_v1 = i;
            best_ed_v2 = t->index;
          }
        }
      }

      // Если не нашли лучшего ребра, выходим из цикла
      if (best_ed_v1 == -1 || best_ed_v2 == -1) break;
      // Перекрашиваем лучшее ребро
      proc_edges[best_ed_v1][best_ed_v2] = 1;
      proc_edges[best_ed_v2][best_ed_v1] = 1;
      best_vertex_color(G, current_colors, best_ed_v1);
      best_vertex_color(G, current_colors, best_ed_v2);
    }
    free(proc_edges);
    current_colors++;
  } while (!GRAPH_check_given_coloring_validity(G, G->color));

  return G->color;
}

/* EXPOSED FUNCTIONS */
unsigned int GRAPH_check_given_coloring_validity(Graph G, unsigned int *colors) {
  for (unsigned int i = 0; i < G->V; i++) {
    for (link t = G->ladj[i]; t != G->z; t = t->next) {
      if (colors[i] == colors[t->index] || colors[i] == 0) {
        return 0;
      }
    }
  }
  return 1;
}

unsigned int GRAPH_check_current_coloring_validity(Graph G) {
  for (unsigned int i = 0; i < G->V; i++) {
    for (link t = G->ladj[i]; t != G->z; t = t->next) {
      if (G->color[i] == G->color[t->index] || G->color[i] == 0) {
        return 0;
      }
    }
  }
  return 1;
}

unsigned int GRAPH_get_edge_count(Graph G) { return G->E; }

unsigned int GRAPH_get_vertex_count(Graph G) { return G->V; }

Graph GRAPH_load_from_file(char *filename) {
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    printf("Error opening file %s\n", filename);
    return NULL;
  }
  char *dot = strrchr(filename, '.');

  if (!strcmp(dot, ".graph")) {
    unsigned int V, E, fmt = 0, ncon = 0;
    char line[4096], *p, *e;
    fgets(line, sizeof(line), fp);
    sscanf(line, "%d %d %d %d", &V, &E, &fmt, &ncon);

    Graph G = GRAPH_init(V);

    if (G == NULL) {
      fclose(fp);
      return NULL;
    }

    unsigned int from = 1;

    while (fgets(line, sizeof(line), fp)) {
      if (line[0] == '%') {
        continue;
      }
      p = line;
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
          for (int i = 0; i < ncon; i++, p = e) {
            to = strtol(p, &e, 10);
          }
        case 0:
          for (;; p = e) {
            to = strtol(p, &e, 10);
            if (p == e) break;
            if (from != to) {
              if (EDGE_insert(G, from - 1, to - 1) == NULL) {
                printf("Couldn't insert edge from %d to %d\n", from, to);
                fclose(fp);
                GRAPH_free(G);
                return NULL;
              }
            }
          }
          break;
        case 11:
          to = strtol(p, &e, 10);
          p = e;
        case 1:
          for (;; p = e) {
            to = strtol(p, &e, 10);
            if (p == e) break;
            if (alternate % 2 == 0 && from != to) {
              if (EDGE_insert(G, from - 1, to - 1) == NULL) {
                printf("Couldn't insert edge from %d to %d\n", from, to);
                fclose(fp);
                GRAPH_free(G);
                return NULL;
              }
            }
            alternate++;
          }
          break;
        default:
          printf("Invalid fmt\n");
          fclose(fp);
          GRAPH_free(G);
          return NULL;
      }
      from++;
    }
    fclose(fp);
    return G;
  } else if (!strcmp(dot, ".gra")) {
    int V = -1;
    char buf[64];
    while (V <= 0) {
      fscanf(fp, "%s", buf);
      if (atoi(buf) > 0) {
        V = atoi(buf);
        // printf("V = %d\n", V);
      } else {
        // printf("Skipped line: %s\n", buf);
      }
    }

    Graph G = GRAPH_init(V);
    if (G == NULL) {
      fclose(fp);
      return NULL;
    }
    for (int i = 0; i < V; i++) {
      fscanf(fp, "%*d:");
      fscanf(fp, " %s ", buf);
      while (strcmp(buf, "#")) {
        int to = atoi(buf);
        if (EDGE_insert(G, i, to) == NULL) {
          printf("Couldn't insert edge from %d to %d\n", i, to);
          fclose(fp);
          GRAPH_free(G);
          return NULL;
        }
        if (EDGE_insert(G, to, i) == NULL) {
          printf("Couldn't insert edge from %d to %d\n", to, i);
          fclose(fp);
          GRAPH_free(G);
          return NULL;
        }
        fscanf(fp, "%s", buf);
      }
    }
    fclose(fp);
    return G;
  } else if (!strcmp(dot, ".col")) {
    size_t V = 0; // Количество вершин
    Graph G = NULL; // Инициализируем граф

    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
      if (line[0] == 'c') continue; // Пропускаем комментарии

      char ch;
      sscanf(line, " %c", &ch);

      switch(ch) {
        case 'p': 
          if (V != 0) {
            // Ошибка: если в файле несколько строк с "p", то проигнорируем все, кроме первой
            printf("Invalid DIMACS format: multiple 'p' lines\n");
            fclose(fp);
            return NULL;
          }

          size_t unused_edge_count;
          sscanf(line, " %*c %*s %zu %zu", &V, &unused_edge_count);
          G = GRAPH_init(V);
          printf("V = %zu\n", V);
          printf("unused_edge_count = %zu\n", unused_edge_count);
          if (G == NULL) {
            fclose(fp);
            return NULL;
          }
          break;
        case 'e':
          if (G == NULL) {
            printf("DIMACS format error: 'p' line missing\n");
            fclose(fp);
            return NULL;
          }
          size_t from, to;
          sscanf(line, " %*c %zu %zu", &from, &to);
          if (from == 0 || to == 0 || from > V || to > V) {
            printf("Invalid vertex ID in DIMACS format\n");
            fclose(fp);
            GRAPH_free(G);
            return NULL;
          }
          if (EDGE_insert(G, from - 1, to - 1) == NULL) {
            printf("Couldn't insert edge from %zu to %zu\n", from, to);
            fclose(fp);
            GRAPH_free(G);
            return NULL;
          }
          break;
        default:
          printf("Invalid DIMACS format\n");
          fclose(fp);
          GRAPH_free(G);
          return NULL;
      }
    }
    fclose(fp);
    return G;  
  } else {
    printf("Invalid extension %s\n", dot);
    fclose(fp);
    return NULL;
  }
  return NULL;
}


Graph GRAPH_init(unsigned int V) {
  Graph G = malloc(sizeof *G);
  if (G == NULL) {
    fprintf(stderr, "Error while allocating the graph\n");
    return NULL;
  }
  G->V = V;
  G->E = 0;
  G->z = LINK_new(-1, NULL);
  G->ladj = malloc(V * sizeof(link)); 
  G->degree = malloc(V * sizeof(unsigned int));
  G->color = malloc(V * sizeof(unsigned int));

  for (unsigned int i = 0; i < V; i++) {
    G->ladj[i] = G->z;
    G->degree[i] = 0;
    G->color[i] = 0;
  }
  return G;
}

void GRAPH_free(Graph G) {
  if (G == NULL) {
    return;
  }
  link next;
  for (int i = 0; i < G->V; i++) {
    for (link t = G->ladj[i]; t != G->z; t = next) {
      next = t->next;
      free(t);
    }
  }
  free(G->ladj);
  free(G->degree);
  free(G->color);
  free(G->z);
  free(G);
}

void GRAPH_ladj_print(Graph G) {
  for (int i = 0; i < G->V; i++) {
    printf("%d -->", i + 1);
    for (link t = G->ladj[i]; t != G->z; t = t->next) {
      printf("%c%d", t == G->ladj[i] ? ' ' : '-', (t->index) + 1);
    }
    putchar('\n');
  }
}
void GRAPH_ladj_print_with_colors(Graph G, unsigned int *colors) {
  for (int i = 0; i < G->V; i++) {
    printf("%d(%d) -->", i + 1, colors[i]);
    for (link t = G->ladj[i]; t != G->z; t = t->next) {
      printf("%c%d(%d)", t == G->ladj[i] ? ' ' : '-', (t->index) + 1, colors[t->index]);
    }
    putchar('\n');
  }
}

unsigned int *GRAPH_get_degrees(Graph G) { return G->degree; }

unsigned long GRAPH_compute_bytes(Graph G) {
  unsigned long long bytes = 0;
  bytes += sizeof(G);
  bytes += sizeof(*G);
  bytes += (G->V + G->E + 1) * sizeof(link);
  bytes += (G->V + G->E + 1) * sizeof(struct node);
  bytes += (2 * G->V) * sizeof(unsigned int);
  return bytes;
}

unsigned int *GRAPH_color(Graph G, char *coloring_method_str,
                          unsigned int n_threads) {
  switch (method_str_to_enum(coloring_method_str)) {
    case seq_greedy:
      //GRAPH_ladj_print(G);
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
    case color_op_ver:
      return color_ver_optimize(G);
      break;  
    case color_op_ed:
      return color_ed_optimize(G);
      break;
    default:
      fprintf(stderr, "Passed coloring method '%s' is not valid!\n", coloring_method_str);
      return NULL;
  }
}