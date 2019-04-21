#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <unordered_set>

void getKmers( std::size_t& nKmers,
			   const std::size_t k,
			   std::vector< std::string >& reads,
			   std::vector< std::string >& kmers,
			   std::string kmer_filename )
{
	std::ifstream in( kmer_filename.c_str() );
	std::string sline;
	std::vector< std::string > vline;
	while ( std::getline( in, sline ) ) {
		vline.push_back( sline );
	}

	std::size_t pos = 0;
	std::string read;
	do {
		if (vline[pos][0] == '>') {
			//finish current read and start a new one
			if (!read.empty()) {
				reads.push_back(read);
				read.clear();
			}
		} else {
			read += vline[pos];
		}

		++pos;
	} while (pos != vline.size());

	if (!read.empty()) //handle the last read
		reads.push_back( read );

	if (k > 0) {
		for (std::size_t i = 0; i < reads.size(); ++i) {
			std::string sline = reads[i];
			std::size_t read_length = sline.size();

			std::size_t nMers = read_length - k + 1;
			for (std::size_t start = 0; start < nMers; ++start) {
				std::string kmer = sline.substr( start, k );
				kmers.push_back( kmer );
			}
		}
	}

	in.close();
	nKmers = kmers.size();
}


std::unordered_set< int > pickSet(const int& N, const int& k, std::mt19937& gen)
{
	std::unordered_set<int> elems;
	for (int r = N - k; r < N; ++r) {
		int v = std::uniform_int_distribution<>(1, r)(gen);
		if (!elems.insert(v).second) {
			elems.insert(r);
		}
	}
	return elems;
}


int main(int argc, char** argv) {
	if (argc < 2) {
		std::cerr << "USAGE: " << argv[0] << " <db> <Y> <Z> <gene> <contiguous>" << std::endl;
		return 1;
	}

	std::size_t nKmers;
	std::vector< std::string > reads;
	std::vector< std::string > kmers;

	getKmers( nKmers, 0, reads, kmers, argv[1] );

	std::cout << reads.size() << " genes in database." << std::endl;

	std::random_device rd;
	std::mt19937 gen( rd() );

	//uniform_int_distribution<size_t> dist1(0,reads.size());
	std::uniform_int_distribution< std::size_t > distChar(0,3);
	std::string sY = argv[2];
	double Y = stod( sY );
	std::string sZ = argv[3];
	std::size_t Z = stoi( sZ );
	std::string su = argv[4];
	std::size_t u = stoi( su );
	std::string contiguous = argv[5];

	std::string alphabet = "ACGT";
	std::string ofname = "fake-gene=" + su + "-Z=" + sZ + "-Y=" + sY + ".fa";

	std::ofstream ofile( ofname.c_str() );

	//for (size_t u = 0; u < nGenes; ++u) {
	std::size_t which = u;
	std::cout << "Selected gene: "
		 << reads[which] << std::endl;
	std::size_t len = Y * reads[which].size();
	std::uniform_int_distribution< std::size_t > dist2(0,reads[which].size());
	std::size_t pos = dist2( gen );
	std::string& read = reads[which];

	if(contiguous == "true") {
		for (std::size_t v = pos; (v < read.size()) && (v < pos + len); ++v) {
			read[v] = alphabet[ distChar( gen ) ];
		}
	}
	else if(contiguous == "false") {
		int n_values = ceil(Y * (double)read.length());
		std::unordered_set< int > random_idxs = pickSet(read.length(), n_values, gen);
		for(auto &idx : random_idxs) {
			read[idx] = alphabet[ distChar( gen ) ];
		}
	}
	else {
		std::cerr << "Argument <contiguous> [5] must be one of [true, false]" << std::endl;
	}

	std::cout << "modified gene: " << read << std::endl;

	std::string geneExpand;
	for (std::size_t i = 0; i < Z; ++i)
		geneExpand += read;

	ofile << "> " << u << std::endl;
	ofile << geneExpand << std::endl;
	//}
	ofile.close();
	return 0;
}
