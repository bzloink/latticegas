 /* Two-dimensional square lattice gas model
    by Andrew M. Launder
    
    Last updated 4.17.2018.
    
    See README for proper code usage. */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

class MySort {
public:
    MySort(int sortint) {this->sortint = sortint;}
    template <typename GenVec1, typename GenVec2>
    bool operator()(const GenVec1& a, const GenVec2& b) {
        return a[sortint] < b[sortint];
    }
    
    int sortint;
};

unsigned long int NLines(char* file) {
 // Determines number of lines in a file.
    std::ifstream data(file);
    
    unsigned long int nlines = 0;
    
    if (!data) {
        return nlines;
    }
    
    std::string line;
    
    while (std::getline(data, line)) {
        ++nlines;
    }
    
    return nlines;
}

std::vector<std::string> LineSplit(std::ifstream& data) {
 // Reads in a string and parses it into a vector split by space delimiters.
    std::string line, item;
    std::getline(data, line);
    std::istringstream linestream(line);
    std::vector<std::string> linesplitstring;
    while (linestream >> item) {
        linesplitstring.push_back(item);
    }
    
    return linesplitstring;
}

std::vector<std::vector<std::string>> ReadInput(char* inputfile, unsigned long int nlines) {
 // Reads in parameters from input.dat.
    std::ifstream data(inputfile);
    std::vector<std::vector<std::string>> params;
    std::vector<std::string>* linesplit = new std::vector<std::string>;
    
    unsigned long int i;
    for (i = 0; i < nlines; ++i) {
        *linesplit = LineSplit(data);
        params.push_back(*linesplit);
    }
    
    return params;
}

int ParamsCheck(std::vector<std::vector<std::string>> params) {
 /* Returns 1 if input.dat does not contain temp parameter.
    Returns 0 if input.dat has the expected format.         */
    unsigned long int i;
    for (i = 1; i < params.size(); ++i) {
        if (params[i][0] == "temp" && params[i].size() == 2) {
            return 0;
        }
    }
    
    return 1;
}

std::vector<std::vector<unsigned long int>> Lattice(unsigned long int nnodes, unsigned long int xdim, unsigned long int ydim, unsigned long int na, unsigned long int nb) {
 // Generates starting lattice.
    std::vector<std::vector<unsigned long int>> lattice(nnodes, std::vector<unsigned long int>(5));
    
    unsigned long int i, j;
    for (i = 0; i < nnodes; ++i) {
        if (i < na) {
            lattice[i][0] = 1;
        }
        else {
            lattice[i][0] = 2;
        }
    }
    for (i = 0; i < xdim; ++i) {
        for (j = 0; j < ydim; ++j) {
            if (j) {
                lattice[i * ydim + j][1] = i * ydim + j - 1;
            }
            else {
                lattice[i * ydim + j][1] = i * ydim + ydim - 1;
            }
            if (j != ydim - 1) {
                lattice[i * ydim + j][2] = i * ydim + j + 1;
            }
            else {
                lattice[i * ydim + j][2] = i * ydim;
            }
            if (i) {
                lattice[i * ydim + j][3] = (i - 1) * ydim + j;
            }
            else {
                lattice[i * ydim + j][3] = (xdim - 1) * ydim + j;
            }
            if (i != xdim - 1) {
                lattice[i * ydim + j][4] = (i + 1) * ydim + j;
            }
            else {
                lattice[i * ydim + j][4] = j;
            }
        }
    }
    
    return lattice;
}

unsigned long int NAB(std::vector<std::vector<unsigned long int>> lattice, unsigned long int nnodes, unsigned long int z) {
 // Determines the total number of A-B type interactions.
    unsigned long int moltype, adjid, nab;
    nab = 0;
    
    unsigned long int i, j;
    for (i = 0; i < nnodes; ++i) {
        moltype = lattice[i][0];
        for (j = 1; j < z + 1; ++j) {
            adjid = lattice[i][j];
            if (lattice[adjid][0] != moltype) {
                ++nab;
            }
        }
    }
    nab /= 2;
    
    return nab;
}

void PrintLattice(std::ofstream& outfile, std::vector<std::vector<unsigned long int>> lattice, unsigned long int nab, unsigned long int xdim, unsigned long int ydim, unsigned long int snap) {
 // Prints simulation parameters.
    if (!snap) {
        outfile << "Initial number of A-B intercell interactions ((n_ab)_i) = ";
    }
    else {
        outfile << "Final number of A-B intercell interactions ((n_ab)_f) in snapshot " << snap << " = ";
    }
    outfile << nab << std::endl;
    
    if (!snap) {
        outfile << "Initial lattice configuration:";
    }
    else {
        outfile << "Final lattice configuration in snapshot " << snap << ":";
    }
    outfile << std::endl;
    
    unsigned long int i, j;
    for (i = 0; i < xdim; ++i) {
        for (j = 0; j < ydim; ++j) {
            if (lattice[i * ydim + j][0] == 1) {
                outfile << " A";
            }
            else {
                outfile << " B";
            }
        }
        outfile << std::endl;
    }
}

std::vector<std::vector<unsigned long int>> EdgeSort(std::vector<std::vector<unsigned long int>> edges, unsigned long int nedges) {
 /* Sorts vectors of edge indices by i) 0th indices; then ii) 1th indices,
    then deletes doublecounted edges.                                      */
    unsigned long int nodelist, oldnodelist;
    oldnodelist = 0;
    
    unsigned long int i;
    for (i = 0; i < nedges; ++i) {
        std::sort(edges[i].begin(), edges[i].end());
    }
    
    std::sort(edges.begin(), edges.end(), MySort(0));
    
    for (i = 1; i < nedges; ++i) {
        if (edges[i][0] != edges[i - 1][0]) {
            nodelist = i;
            if (nodelist > oldnodelist + 1) {
                std::sort(edges.begin() + oldnodelist, edges.begin() + nodelist, MySort(1));
            }
            oldnodelist = nodelist;
        }
        if (i == nedges - 1) {
            std::sort(edges.begin() + oldnodelist, edges.end(), MySort(1));
        }
    }
    
    std::vector<std::vector<unsigned long int>> finaledges(nedges / 2, std::vector<unsigned long int>(2));
    for (i = 0; i < nedges / 2; ++i) {
        finaledges[i] = edges[2 * i];
    }
    
    return finaledges;
}

void PrintGraphs(std::ofstream& AAgraphfile, std::ofstream& BBgraphfile, std::ofstream& ABgraphfile, std::vector<std::vector<unsigned long int>> lattice, unsigned long int nnodes, int z) {
 // Writes .GraphGeod files.
    std::vector<std::vector<unsigned long int>>* AAedges = new std::vector<std::vector<unsigned long int>>;
    std::vector<std::vector<unsigned long int>>* BBedges = new std::vector<std::vector<unsigned long int>>;
    std::vector<std::vector<unsigned long int>>* ABedges = new std::vector<std::vector<unsigned long int>>;
    std::vector<unsigned long int>* edge = new std::vector<unsigned long int>(2);
    unsigned long int adjid;
    
    unsigned long int i;
    int a;
    for (i = 0; i < nnodes; ++i) {
        (*edge)[0] = i + 1;
        for (a = 1; a < z + 1; ++a) {
            adjid = lattice[i][a];
            (*edge)[1] = adjid + 1;
            if (lattice[adjid][0] == lattice[i][0]) {
                if (lattice[i][0] == 1) {
                    AAedges->push_back(*edge);
                }
                else {
                    BBedges->push_back(*edge);
                }
            }
            else {
                ABedges->push_back(*edge);
            }
        }
    }
    delete edge;
    
    *AAedges = EdgeSort(*AAedges, AAedges->size());
    for (i = 0; i < AAedges->size(); ++i) {
        AAgraphfile << (*AAedges)[i][0] << " " << (*AAedges)[i][1];
        for (a = 0; a < 7; ++a) {
            AAgraphfile << " " << 0;
        }
        AAgraphfile << std::endl;
    }
    delete AAedges;
    *BBedges = EdgeSort(*BBedges, BBedges->size());
    for (i = 0; i < BBedges->size(); ++i) {
        BBgraphfile << (*BBedges)[i][0] << " " << (*BBedges)[i][1];
        for (a = 0; a < 7; ++a) {
            BBgraphfile << " " << 0;
        }
        BBgraphfile << std::endl;
    }
    delete BBedges;
    *ABedges = EdgeSort(*ABedges, ABedges->size());
    for (i = 0; i < ABedges->size(); ++i) {
        ABgraphfile << (*ABedges)[i][0] << " " << (*ABedges)[i][1];
        for (a = 0; a < 7; ++a) {
            ABgraphfile << " " << 0;
        }
        ABgraphfile << std::endl;
    }
    delete ABedges;
}

int main(int argc, char** argv) {
 // Generates lattice and writes outputs.
    double temp, eaa, ebb, eab;
    eaa = ebb = 40;
    eab = 5;
    unsigned long int nsnaps, niters, xdim, ydim, na, nb;
    nsnaps = 1;
    niters = 50000;
    xdim = ydim = 10;
    na = nb = 50;
    
    std::string inputfilestr = "input.dat";
    char* inputfile = &inputfilestr[0u];
    unsigned long int nlines = NLines(inputfile);
    if (!nlines) {
        std::cerr << "Have you correctly specified your input file?" << std::endl;
        return EXIT_FAILURE;
    }
    std::vector<std::vector<std::string>>* params = new std::vector<std::vector<std::string>>;
    *params = ReadInput(inputfile, nlines);
    int paramscheck = ParamsCheck(*params);
    if (paramscheck) {
        std::cerr << "Please specify temp parameter." << std::endl;
        return EXIT_FAILURE;
    }
    
    unsigned long int i, j;
    int a;
    for (i = 1; i < params->size(); ++i) {
        if ((*params)[i].size() == 1) {
            std::cerr << "Warning: some keywords have unspecified values." << std::endl;
        }
        else if ((*params)[i][0] == "temp") {
            std::istringstream tempstr((*params)[i][1]);
            tempstr >> temp;
        }
        else if ((*params)[i][0] == "nsnaps") {
            std::istringstream nsnapsstr((*params)[i][1]);
            nsnapsstr >> nsnaps;
        }
        else if ((*params)[i][0] == "niters") {
            std::istringstream nitersstr((*params)[i][1]);
            nitersstr >> niters;
        }
        else if ((*params)[i][0] == "xdim") {
            std::istringstream xdimstr((*params)[i][1]);
            xdimstr >> xdim;
        }
        else if ((*params)[i][0] == "ydim") {
            std::istringstream ydimstr((*params)[i][1]);
            ydimstr >> ydim;
        }
        else if ((*params)[i][0] == "na") {
            std::istringstream nastr((*params)[i][1]);
            nastr >> na;
        }
        else if ((*params)[i][0] == "nb") {
            std::istringstream nbstr((*params)[i][1]);
            nbstr >> nb;
        }
        else if ((*params)[i][0] == "eaa") {
            std::istringstream eaastr((*params)[i][1]);
            eaastr >> eaa;
            eaa *= -1;
        }
        else if ((*params)[i][0] == "ebb") {
            std::istringstream ebbstr((*params)[i][1]);
            ebbstr >> ebb;
            ebb *= -1;
        }
        else if ((*params)[i][0] == "eab") {
            std::istringstream eabstr((*params)[i][1]);
            eabstr >> eab;
            eab *= -1;
        }
        else {
            std::cerr << "Warning: option \"" << (*params)[i][0] << "\" not yet implemented." << std::endl;
        }
    }
    delete params;
    if (xdim < 1 || ydim < 1) {
        std::cerr << "Invalid: lattice dimensions must be positive integers." << std::endl;
        return EXIT_FAILURE;
    }
    unsigned long int nnodes, nanb;
    nnodes = xdim * ydim;
    nanb = na + nb;
    if (nnodes != nanb) {
        std::cerr << "na and nb must sum to xdim * ydim!" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::random_device seed;
    std::mt19937 gen(seed());
    std::uniform_int_distribution<unsigned long int> distnodes(0, nnodes - 1);
    std::uniform_int_distribution<unsigned long int> distmax(0, RAND_MAX);
    
    int z = 4;
    double kb, w;
    kb = 1.38065;
    w = eab - 0.5 * eaa - 0.5 * ebb;
    std::vector<double>* probs = new std::vector<double>(z);
    for (a = 0; a < z; ++a) {
        (*probs)[a] = exp(-((2 * ((double) a) + 2) * w) / (kb * temp));
    }
    double redtemp = kb * temp / w;
    
    std::ofstream datafile("lattice.dat");
    datafile << "Interchange energy (w) = " << w << std::endl;
    datafile << "Reduced temperature (kT / w) = " << redtemp << std::endl;
    datafile << "Probability of adding [x] A-B intercell interactions:" << std::endl;
    for (i = 0; i < probs->size(); ++i) {
        datafile << " " << 2 * (i + 1) << ": " << (*probs)[i] << std::endl;
    }
    
    std::vector<std::vector<unsigned long int>>* initlattice = new std::vector<std::vector<unsigned long int>>(nnodes, std::vector<unsigned long int>(z + 1));
    *initlattice = Lattice(nnodes, xdim, ydim, na, nb);
    unsigned long int initnab = NAB(*initlattice, nnodes, z);
    PrintLattice(datafile, *initlattice, initnab, xdim, ydim, 0);
    
    std::vector<std::vector<unsigned long int>>* lattice = new std::vector<std::vector<unsigned long int>>(nnodes, std::vector<unsigned long int>(z + 1));
    double prob, totenergy;
    unsigned long int nab, molid1, moltype1, adjcount1, molid2, moltype2, adjcount2, adjid;
    int nabdiff, adjint, probint;
    adjint = probint = 0;
    for (i = 0; i < nsnaps; ++i) {
        *lattice = *initlattice;
        nab = initnab;
        
        std::ostringstream AAgraphfilestr;
        AAgraphfilestr << "AAlattice" << i + 1 << ".GraphGeod";
        std::ofstream AAgraphfile(AAgraphfilestr.str());
        std::ostringstream BBgraphfilestr;
        BBgraphfilestr << "BBlattice" << i + 1 << ".GraphGeod";
        std::ofstream BBgraphfile(BBgraphfilestr.str());
        std::ostringstream ABgraphfilestr;
        ABgraphfilestr << "ABlattice" << i + 1 << ".GraphGeod";
        std::ofstream ABgraphfile(ABgraphfilestr.str());
        std::ostringstream energyfilestr;
        energyfilestr << "lattice" << i + 1 << ".energy";
        std::ofstream energyfile(energyfilestr.str());
        
        for (j = 0; j < niters; ++j) {
            adjint = adjcount1 = adjcount2 = 0;
            molid1 = distnodes(gen);
            moltype1 = (*lattice)[molid1][0];
            for (a = 1; a < z + 1; ++a) {
                adjid = (*lattice)[molid1][a];
                if ((*lattice)[adjid][0] == moltype1) {
                    ++adjcount1;
                }
            }
            molid2 = distnodes(gen);
            moltype2 = (*lattice)[molid2][0];
            while (moltype1 == moltype2) {
                molid2 = distnodes(gen);
                moltype2 = (*lattice)[molid2][0];
            }
            for (a = 1; a < z + 1; ++a) {
                adjid = (*lattice)[molid2][a];
                if (adjid == molid1) {
                    adjint = 1;
                }
                if ((*lattice)[adjid][0] == moltype2) {
                    ++adjcount2;
                }
            }
            nabdiff = 2 * (adjcount1 + adjcount2 - z);
            if (adjint) {
                nabdiff += 2;
            }
            if (nabdiff > 0) {
                prob = (*probs)[nabdiff / 2 - 1];
                if ((distmax(gen) / ((double) RAND_MAX)) <= prob) {
                    probint = 1;
                }
            }
            if (nabdiff <= 0 || probint) {
                probint = 0;
                nab += nabdiff;
                if ((*lattice)[molid1][0] == 1) {
                    ++(*lattice)[molid1][0];
                    --(*lattice)[molid2][0];
                }
                else {
                    --(*lattice)[molid1][0];
                    ++(*lattice)[molid2][0];
                }
            }
            
            if (j == niters - 1) {
                PrintLattice(datafile, *lattice, nab, xdim, ydim, i + 1);
                PrintGraphs(AAgraphfile, BBgraphfile, ABgraphfile, *lattice, nnodes, z);
            }
            totenergy = w * nab + 0.5 * z * (eaa * na + ebb * nb);
            energyfile << j + 1 << " " << totenergy << std::endl;
        }
    }
    delete probs;
    delete initlattice;
    delete lattice;
    
    return 0;
}
