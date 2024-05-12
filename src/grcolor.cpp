#include <dirent.h>
#include <iostream>
#include <string>
#include <sys/sysinfo.h>
#include <ctime>
#include <unistd.h>
#include <fstream>
#include <iomanip>

#include "graph.h"
#include "util.h"

#define N_COLORING_METHODS 6

struct Results {
  std::string graph_name;
  std::string coloring_method;
  unsigned int vertex_count;
  unsigned int colors_used;
  unsigned int n_threads;
  double coloring_time;
};

int main(int argc, char *argv[]) {
  srand(time(NULL));
  std::vector<std::string> graphs_filenames;
  int number_of_graphs = 0;
  int export_to_csv = 0;
  int par_only = 0;
  int n_threads = get_nprocs();
  int iterations = 1;
  if (argc > 1) {
    for (int i = 1; i < argc; i++) {
      std::string arg = argv[i];
      /* flag '-t' or '--threads' to specify how many threads to use in parallel
       * algorithms */
      if (arg == "--threads" || arg == "-t") {
        if (i + 1 != argc) {
          n_threads = std::stoi(argv[i + 1]);
          if (n_threads <= 0) {
            std::cout << "Error: '-t|--threads' flag is specified but the number of "
                      << "threads is invalid! (negative, zero or not numeric)\n";
            return 1;
          }
        } else {
          std::cout << "Error: '-t|--threads' flag is specified without the number of "
                    << "threads!\n";
          return 1;
        }

        i++;  // Move to the next flag
        continue;
      }

      /* flag '-n' to specify how many times the coloring should be repeated */
      if (arg == "-n") {
        if (i + 1 != argc) {
          iterations = std::stoi(argv[i + 1]);
          if (iterations <= 0) {
            std::cout << "Error: '-n' flag is specified but the number of "
                      << "iterations is invalid! (negative, zero or not numeric)\n";
            return 1;
          }
        } else {
          std::cout << "Error: '-n' flag is specified without the number of "
                    << "iterations!\n";
          return 1;
        }

        i++;  // Move to the next flag
        continue;
      }

      /* flag '--csv' to specify whether or not we want to export results to
       * csv */
      if (arg == "--csv") {
        export_to_csv = 1;
        continue;
      }

      /* flag '--par' to specify whether or not we want to use parallel
       * algorithms only */
      if (arg == "--par") {
        par_only = 1;
        continue;
      }

      /* if an argument is not a known flag, it's treated as a graph's filename
       */
      graphs_filenames.push_back(argv[i]);
    }
  }

/* if no graphs' filenames were passed as argument */
  if (graphs_filenames.empty()) {
    /* then use the graphs in the /graphs subfolder */
    DIR *d;
    struct dirent *dir;
    d = opendir("graphs");
    if (d) {
      while ((dir = readdir(d)) != NULL) {
        std::string filename(dir->d_name);
        size_t last_dot = filename.find_last_of('.');
        if (last_dot != std::string::npos) {
            std::string dot = filename.substr(last_dot);
            if (dot == ".graph" || dot == ".gra" || dot == ".col") {
                /* if the extension is .graph or .gra, add to graphs to be colored */
                graphs_filenames.push_back("graphs/" + filename);
            }
        }
    }
    closedir(d);
    } else {
      std::cout << "Error opening graphs/ folder\n";
      return 2;
    }
  }

  if (n_threads > get_nprocs()) {
    /* if the number of threads exceeds the number of available logic
     * processors, the coloring is slower and more error prone */
    std::cout << "Lowering the number of threads to " << get_nprocs() 
              << " (number of available logic processors in the system)\n";
    n_threads = get_nprocs();
  }

  std::vector<std::string> coloring_methods = { "seq_greedy", "seq_ldf", "rec_rlf", "par_jp", "par_opt_ver", "par_opt_ed"};

  Results res;  // struct to hold the results of a coloring (in terms of time
                // and colors used, not the coloring itself)
  std::ofstream csv_file;
  std::string csv_filename;
  if (export_to_csv) {
    /* create the csv file where to export results */
    time_t t = time(NULL);
    struct tm now = *localtime(&t);
    std::ostringstream oss;
    oss << "results/results_" 
        << std::to_string(now.tm_year + 1900) << "-" 
        << std::setw(2) << std::setfill('0') << std::to_string(now.tm_mon + 1) << "-" 
        << std::setw(2) << std::setfill('0') << std::to_string(now.tm_mday) << "_" 
        << std::to_string(now.tm_hour) << "-" 
        << std::to_string(now.tm_min) << "-" 
        << std::to_string(now.tm_sec) << ".csv";
    csv_filename = oss.str();
    csv_file.open(csv_filename);
    if (!csv_file.is_open()) {
      std::cout << "Error opening " << csv_filename << "\n";
      return 3;
    }
    csv_file << "graph_name,vertex_count,coloring_method,n_threads,coloring_time,"
             << "colors_used\n";  // add the header line
    csv_file.close();
  }

  res.n_threads = n_threads;

  if (number_of_graphs == 0) {
    printf("No graphs found in the /graphs subfolder!\n");
  }

 /* for each graph */
  for (std::vector<std::string>::size_type i = 0; i < graphs_filenames.size(); i++) {
    double start, finish;

    /* load the graph from file */
    start = UTIL_get_time();
    auto G = GRAPH_load_from_file(graphs_filenames[i]);
    finish = UTIL_get_time();

    /* take the portion of the filename after the last '/' slash */
    std::string s = graphs_filenames[i];
    std::string last = s.substr(s.find_last_of('/') + 1);
    res.graph_name = last;

    if (G != nullptr) {
      res.vertex_count = GRAPH_get_vertex_count(G);
      std::cout << std::left << std::setw(16) << "GRAPH NAME" << " | "
                << std::setw(11) << "LOADED IN" << " | "
                << std::setw(11) << "MAX DEGREE" << " | "
                << "ESTIMATED MEMORY FOOTPRINT\n";

      std::cout << std::setw(16) << last << " | "
                << std::setw(11) << finish - start << " | "
                << std::setw(11) << UTIL_max_in_array(GRAPH_get_degrees(G)) << " | "
                << ((double)GRAPH_compute_bytes(G)) / 1024 / 1024 << " MB\n";
      std::cout << std::string(60, '_') << "\n";

      for (int k = 0; k < iterations; k++) {
        /* for each iteration */
        if (iterations > 1) {
          std::cout << "Iteration " << k + 1 << " of " << iterations << "\n";
        }

        std::cout << std::setw(16) << "COLOR METHOD" << " | "
                  << std::setw(11) << "COLORED IN" << " | "
                  << std::setw(11) << "COLORS USED" << " | "
                  << "VALID?\n";

        for (int method_number = 0; method_number < N_COLORING_METHODS;
             method_number++) {
          /* for each coloring method */
          
          /* skip sequential methods if --par flag had been set */
          if (par_only &&
              (coloring_methods[method_number] == "seq_greedy" ||
               coloring_methods[method_number] == "seq_ldf")) {
            continue;
          }

          res.coloring_method = coloring_methods[method_number];
          /* color the graph */
          start = UTIL_get_time();
          std::vector<unsigned int> colors =
              GRAPH_color(G, res.coloring_method, n_threads);
          finish = UTIL_get_time();

          /* if the coloring succeeds (i.e: GRAPH_color() returns something !=
           * NULL)*/
          if (!colors.empty()) {
            res.colors_used = UTIL_max_in_array(colors);
            res.coloring_time = finish - start;
            std::cout << std::setw(16) << res.coloring_method << " | "
                      << std::fixed << std::setw(11) << res.coloring_time << " | "
                      << std::setw(11) << res.colors_used << " | ";

            /* check whether or not the produced coloring is valid */
            if (GRAPH_check_given_coloring_validity(G, colors)) {
              std::cout << "YES \n";

              if (export_to_csv) {
                // export to csv if flag had been set
                csv_file.open(csv_filename, std::ios::app);
                if (!csv_file.is_open()) {
                  std::cout << "Error opening " << csv_filename << " in append mode\n";
                } else {
                  csv_file << res.graph_name << "," << res.vertex_count << ","
                           << res.coloring_method << "," << res.n_threads << ","
                           << std::fixed << res.coloring_time << "," << res.colors_used << "\n";
                  csv_file.close();
                }
              }
            } else {
              std::cout << "NO \n";
            }
          }
        }
      }
      std::cout << "\n\n";
      /* after coloring the graph with each method, finally free it */
      GRAPH_free(G);
      sleep(1);
    }
  }

  return 0;
}
