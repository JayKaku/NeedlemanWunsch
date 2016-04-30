// Oðuz Kerem Tural
// 150130125

#ifndef NEEDLEMANWUNSCH_H_
#define NEEDLEMANWUNSCH_H_

#include <iostream>
#include <string>
#include <map>

enum AminoAcid
{
	A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, B, Z, X, MINSCORE
};

class NeedlemanWunsch
{
private:

	std::string _alignment_a;
	std::string _alignment_b;
	std::string _signs;

	std::string _dna_a, _dna_b;
	size_t _length_of_a, _length_of_b;

	int **_similarity_matrix;
	int _subsitition_matrix[24][24];
	int _similarity, _gaps, _identity;
	std::map<char, int> _dna_int_map;

public:
	// Constructor and destructor.
	NeedlemanWunsch(std::string dna_a, std::string dna_b);
	~NeedlemanWunsch();

	// Necessary methods.
	void calculate_similarity();
	void populate_subs_matrix(std::string subs_matrix_file);
	void dna_align();
	void print_results();
};

#endif
