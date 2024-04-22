#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <memory>
#include <string>

#include "util.h"

struct Edge {
  int from;
  int to;
};

struct Graph;

unsigned int GRAPH_get_edge_count(std::shared_ptr<Graph> G);
unsigned int GRAPH_get_vertex_count(std::shared_ptr<Graph> G);
std::shared_ptr<Graph> GRAPH_load_from_file(std::string filename);
std::shared_ptr<Graph> GRAPH_init(unsigned int V);
void GRAPH_free(std::shared_ptr<Graph> G);
void GRAPH_ladj_print(std::shared_ptr<Graph> G);
void GRAPH_ladj_print_with_colors(std::shared_ptr<Graph> G, std::vector<unsigned int> &colors);
std::vector<unsigned int> GRAPH_color(std::shared_ptr<Graph> G, std::string coloring_method_str, unsigned int n_threads);
std::vector<unsigned int> GRAPH_get_degrees(std::shared_ptr<Graph> G);
bool GRAPH_check_given_coloring_validity(std::shared_ptr<Graph> G, std::vector<unsigned int> &colors);
bool GRAPH_check_current_coloring_validity(std::shared_ptr<Graph> G);
unsigned long GRAPH_compute_bytes(std::shared_ptr<Graph> G);

#endif 
