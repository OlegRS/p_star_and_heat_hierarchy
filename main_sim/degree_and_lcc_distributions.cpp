#include "../GraphOS/src/matrices/A_matrix.hpp"
#include "../GraphOS/src/aux_math.hpp"
#include "../GraphOS/src/matrices/col_vector.tpp"
#include "../GraphOS/src/matrices/matrix.tpp"
#include "../GraphOS/src/matrices/symm_matrix.tpp"

#define N_NODES 1000
#define N_PAIRS (N_NODES*(N_NODES-1)/2.)
#define N_ITERS_PER_LINK 10
#define N_ITERS ((N_NODES*(N_NODES-1))/2*N_ITERS_PER_LINK)
#define N_AVRGING 100000

#define Ld .3
#define t2_MF -7500
#define t3_MF 4500
#define t1_MF ( 1/2.*log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld )

#define t1 t1_MF
#define t2 (2*N_NODES/(double)(N_NODES-2) * t2_MF) // Conversion from MF to PS couplings
#define t3 (6*N_NODES*N_NODES/(double)(N_NODES*N_NODES-5*N_NODES+6) * t3_MF) // Conversion from MF to PS couplings

#define INITIALIZE_RANDOMLY false

//#define MEAN_FIELD

using namespace std;;

#ifdef MEAN_FIELD
double t_rescaled[3] = {2*t1_MF, 2*t2_MF/N_PAIRS, 2*t3_MF/(N_PAIRS*N_PAIRS)};

double Hamiltonian(const unsigned int &L) {
  return -t_rescaled[0]*L - t_rescaled[1]*L*L - t_rescaled[2]*L*L*L;
}

#endif

int main() {
  prng rnd(RAND_MAX);
#ifdef MEAN_FIELD
  col_vector<double> params(3, (double[3]){t1_MF, t2_MF, t3_MF});
  cout << string("MF__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_MF_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv" << '\n';
  
  ofstream ofs_deg(string("../data/mean_field_model/degree_distributions/MF__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_MF_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv");
  ofstream ofs_cc(string("../data/mean_field_model/clustering_distributions/MF__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_MF_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv");
  
  A_matrix A(N_NODES);

  if(!INITIALIZE_RANDOMLY) // Initial thermalisation
    A.single_link_MF_GB_Metropolis_generator(Hamiltonian, N_ITERS*10, rnd);
  
  for(unsigned int i=0; i<N_AVRGING; ++i) {
    A.single_link_MF_GB_Metropolis_generator(Hamiltonian, N_ITERS, rnd, INITIALIZE_RANDOMLY);
    ofs_deg << A.degree_sequence_col_vec();
    ofs_deg.flush();
    ofs_cc << A.loc_clust_coeff_col_vec();
    ofs_cc.flush();
    cerr << i << '\t'; 
  }

#else
  col_vector<double> params(3, (double[3]){t1, t2, t3});
  cout << string("PS__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_MF_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv" << '\n';
    
  ofstream ofs_deg(string("../data/p_star_model/degree_distributions/PS__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_MF_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv");
  ofstream ofs_cc(string("../data/p_star_model/clustering_distributions/PS__N_NODES_") + to_string(N_NODES) + "__N_ITERS_PER_LINK_" + to_string(N_ITERS_PER_LINK) + "__t1_MF_" + to_string(t1_MF) + "__t2_MF_" + to_string(t2_MF) + "__t3_MF_" + to_string(t3_MF) + "__N_AVRGNG_" + to_string(N_AVRGING) + "__INITIALIZE_RANDOMLY_" + to_string(INITIALIZE_RANDOMLY) + ".csv");

  A_matrix A(N_NODES);

  if(!INITIALIZE_RANDOMLY) // Initial thermalisation
    A.sample_p_star_model_with_single_link_Metropolis(N_ITERS*10, rnd, params);
  
  for(unsigned int i=0; i<N_AVRGING; ++i) {
    A.sample_p_star_model_with_single_link_Metropolis(N_ITERS, rnd, params, INITIALIZE_RANDOMLY);
    ofs_deg << A.degree_sequence_col_vec();
    ofs_deg.flush();
    ofs_cc << A.loc_clust_coeff_col_vec();
    ofs_cc.flush();
    cerr << i << '\t'; 
  }
  
#endif
  ofs_deg.close();
  ofs_cc.close();
  return 0;
}
