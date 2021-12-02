#include "../GraphOS/src/graph.hpp"

using namespace std;;

int main() {
  string file_name = "P_STAR__N_NODES_50__N_ITERS_PER_LINK_100__t1_MF_13.151388__t2_MF_-75__t3_MF_25__INITIALIZE_RANDOMLY_0__TIME_1609962861.sif";
  graph gr;
  gr.load_from_sif(string("../data/p_star_model/graphs/") + file_name);
  
  // ofstream ofs_degree(string("../data/triad_model/degree_distributions/") + file_name.erase(117) + ".csv");
  // ofstream ofs_degree_distribution(string("../data/triad_model/degree_distributions/") + file_name.erase(117) + "_distribution.csv");
  // ofstream ofs_lcc(string("../data/triad_model/clustering_distributions/") + file_name.erase(117) + ".csv");

  ofstream ofs_degree(string("../data/triad_model/degree_distributions/") + file_name.erase(119) + ".csv");
  ofstream ofs_degree_distribution(string("../data/triad_model/degree_distributions/") + file_name.erase(119) + "_distribution.csv");
  ofstream ofs_lcc(string("../data/triad_model/clustering_distributions/") + file_name.erase(119) + ".csv");

  cout << gr;
  
  cout << "max_degree = " << gr.max_degree() << '\n'
       << "min_degree = " << gr.min_degree() << '\n';
  cout << "degree_assortativity = " << gr.degree_assortativity() << '\n';
  cout << "global_clustering_coefficient = " << gr.global_clustering_coefficient() << '\n';
  cout << "average_local_clustering_coefficient = " << gr.clustering_coefficient_sequence().avrg() << '\n';
  cout << "N_nodes_with_degree_4 = " << gr.nodes_with_degree_col_vec(4).size() << '\n'
       << "N_nodes_with_degree_5 = " << gr.nodes_with_degree_col_vec(5).size() << '\n'
       << "N_nodes_with_degree_6 = " << gr.nodes_with_degree_col_vec(6).size() << '\n';

  ofs_degree << gr.degree_sequence_col_vec();
  ofs_lcc << gr.clustering_coefficient_sequence();

  ofs_degree.close();
  ofs_lcc.close();
  
  return 0;
}
