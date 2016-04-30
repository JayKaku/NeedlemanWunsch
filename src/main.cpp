// Oðuz Kerem Tural
// 150130125

#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

#include "NeedlemanWunsch.h"

// Usage: ./program dna_file1 dna_file2
int main(int argc, char **argv) 
{
	// File Parsing.
	std::fstream compare_one, compare_two;
	if (argc != 3) 
	{
		std::cout << "Warning: Using default values! \n";
		compare_one.open("mouse-hemoglobin-sequence.fasta");
		compare_two.open("human-hemoglobin-sequence.fasta");
	}
	else
	{
		compare_one.open(argv[1]);
		compare_one.open(argv[2]);
	}

	
	if (!compare_one || !compare_two)
	{
		std::cerr << "Error (2): Necessary files not found! \n";
		exit(2);
	}
	
	std::string dna_a = "", dna_b = "", line;
	for (int i = 0; std::getline(compare_one, line); i++)
	{
		if (i != 0) 
		{
			dna_a += line;
		}
	}

	for (int i = 0; std::getline(compare_two, line); i++)
	{
		if (i != 0)
		{
			dna_b += line;
		}
	}

	// Application of algorithm.
	NeedlemanWunsch *nw = new NeedlemanWunsch(dna_a, dna_b);
	nw->populate_subs_matrix("BLOSUM62.txt");
	nw->calculate_similarity();
	nw->dna_align();
	nw->print_results();
	getchar();
	return EXIT_SUCCESS;
}