 /* Two-dimensional square lattice gas model
    by Andrew M. Launder
    
    Last updated 07.27.2018.
    
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
    Returns 0 if input.dat has the expected format. */
    unsigned long int i;
    for (i = 1; i < params.size(); ++i) {
        if (params[i][0] == "temp" && params[i].size() == 2) {
            return 0;
        }
    }
    
    return 1;
}

std::vector<std::vector<unsigned long int>> Lattice(unsigned long int nnodes, unsigned long int xdim, unsigned long int ydim, int order, int ncells, unsigned long int nsmall, int minortype, std::mt19937 gen, std::uniform_int_distribution<unsigned long int> distnodes, std::uniform_int_distribution<unsigned long int> distmax) {
 // Generates random starting lattice.
    std::vector<std::vector<unsigned long int>> lattice(nnodes, std::vector<unsigned long int>(ncells));
    unsigned long int molid, count;
    count = 0;
    while (count < nsmall) {
        molid = distnodes(gen);
        if (distmax(gen) / ((double) RAND_MAX) <= 0.5 && !lattice[molid][0]) {
            lattice[molid][0] = minortype;
            ++count;
        }
    }
    unsigned long int majortype;
    if (minortype == 1) {
        majortype = 2;
    }
    else {
        majortype = 1;
    }
    
    unsigned long int i, j;
    for (i = 0; i < nnodes; ++i) {
        if (!lattice[i][0]) {
            lattice[i][0] = majortype;
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
    if (order >= 2) {
        for (i = 0; i < xdim; ++i) {
            for (j = 0; j < ydim; ++j) {
                lattice[i * ydim + j][5] = lattice[lattice[i * ydim + j][1]][3];
                lattice[i * ydim + j][6] = lattice[lattice[i * ydim + j][2]][3];
                lattice[i * ydim + j][7] = lattice[lattice[i * ydim + j][1]][4];
                lattice[i * ydim + j][8] = lattice[lattice[i * ydim + j][2]][4];
                if (order == 3) {
                    lattice[i * ydim + j][9] = lattice[lattice[i * ydim + j][1]][1];
                    lattice[i * ydim + j][10] = lattice[lattice[i * ydim + j][2]][2];
                    lattice[i * ydim + j][11] = lattice[lattice[i * ydim + j][3]][3];
                    lattice[i * ydim + j][12] = lattice[lattice[i * ydim + j][4]][4];
                }
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

void PrintLattice(std::ostream& outfile, std::vector<std::vector<unsigned long int>> lattice, unsigned long int nab, unsigned long int xdim, unsigned long int ydim, unsigned long int snap, int color) {
 /* If color is 0, then lattice is written to file.
    If color is 1, then lattice is written to file
    and ANSI color-coded lattice is printed. */
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
                if (!color) {
                    outfile << " A";
                }
                else {
                    outfile << " \e[31;1mA\e[0m";
                }
            }
            else {
                if (!color) {
                    outfile << " B";
                }
                else {
                    outfile << " \e[30;1mB\e[0m";
                }
            }
        }
        outfile << std::endl;
    }
}

std::vector<int> ABNodes(std::vector<std::vector<unsigned long int>> lattice, unsigned long int nnodes, int z) {
 /* Returns list of ints indicating present (1) or not present (0) nodes
    in the given AB edge list.*/
    std::vector<std::vector<unsigned long int>>* initnodelist = new std::vector<std::vector<unsigned long int>>;
    std::vector<unsigned long int>* node = new std::vector<unsigned long int>(2);
    std::vector<unsigned long int>* adjnode = new std::vector<unsigned long int>(2);
    unsigned long int adjid;
    
    unsigned long int i, j;
    int a;
    for (i = 0; i < nnodes; ++i) {
        (*node)[0] = lattice[i][0];
        (*node)[1] = i + 1;
        for (a = 1; a < z + 1; ++a) {
            adjid = lattice[i][a];
            (*adjnode)[0] = lattice[adjid][0];
            (*adjnode)[1] = adjid + 1;
            if (lattice[adjid][0] == lattice[i][0]) {
                continue;
            }
            else {
                initnodelist->push_back(*node);
                initnodelist->push_back(*adjnode);
            }
        }
    }
    delete adjnode;
    
    std::vector<std::vector<unsigned long int>>* nodelist = new std::vector<std::vector<unsigned long int>>;
    (*node)[0] = (*initnodelist)[0][0];
    (*node)[1] = (*initnodelist)[0][1];
    nodelist->push_back(*node);
    for (i = 1; i < initnodelist->size(); ++i) {
        for (j = 0; j < nodelist->size(); ++j) {
            if ((*initnodelist)[i][1] == (*nodelist)[j][1]) {
                break;
            }
            if (j == nodelist->size() - 1) {
                (*node)[0] = (*initnodelist)[i][0];
                (*node)[1] = (*initnodelist)[i][1];
                nodelist->push_back(*node);
            }
        }
    }
    delete initnodelist;
    delete node;
    
    std::vector<int> nodes(nnodes);
    for (i = 0; i < nnodes; ++i) {
        nodes[i] = 0;
    }
    for (i = 0; i < nodelist->size(); ++i) {
        nodes[(*nodelist)[i][1] - 1] = (*nodelist)[i][0];
    }
    delete nodelist;
    
    return nodes;
}

void PrintXYZs(std::ofstream& xyzfile, std::ofstream& AAxyzfile, std::ofstream& BBxyzfile, std::ofstream& ABxyzfile, std::vector<std::vector<unsigned long int>> lattice, unsigned long int nnodes, unsigned long int na, unsigned long int nb, unsigned long int xdim, unsigned long int ydim, int z) {
 // Writes .xyz files.
    xyzfile << nnodes << std::endl;
    xyzfile << std::endl;
    AAxyzfile << na << std::endl;
    AAxyzfile << std::endl;
    BBxyzfile << nb << std::endl;
    BBxyzfile << std::endl;
    std::vector<int>* abnodes = new std::vector<int>(nnodes);
    *abnodes = ABNodes(lattice, nnodes, z);
    unsigned long int nabnodes = 0;
    
    unsigned long int i, j;
    for (i = 0; i < nnodes; ++i) {
        if ((*abnodes)[i]) {
            ++nabnodes;
        }
    }
    ABxyzfile << nabnodes << std::endl;
    ABxyzfile << std::endl;
    
    for (i = 0; i < xdim; ++i) {
        for (j = 0; j < ydim; ++j) {
            if (lattice[i * ydim + j][0] == 1) {
                xyzfile << "A " << i + 1 << " " << j + 1 << " 0" << std::endl;
                AAxyzfile << "A " << i + 1 << " " << j + 1 << " 0" << std::endl;
            }
            else {
                xyzfile << "B " << i + 1 << " " << j + 1 << " 0" << std::endl;
                BBxyzfile << "B " << i + 1 << " " << j + 1 << " 0" << std::endl;
            }
            if ((*abnodes)[i * ydim + j] == 1) {
                ABxyzfile << "A " << i + 1 << " " << j + 1 << " 0" << std::endl;
            }
            else if ((*abnodes)[i * ydim + j] == 2) {
                ABxyzfile << "B " << i + 1 << " " << j + 1 << " 0" << std::endl;
            }
        }
    }
    delete abnodes;
}

std::vector<std::vector<double>> EdgeSort(std::vector<std::vector<double>> edges, unsigned long int nedges) {
 /* Sorts vectors of edge indices by i) 0th indices; then ii) 1th indices,
    then deletes doublecounted edges. */
    unsigned long int nodelist, oldnodelist;
    oldnodelist = 0;
    
    unsigned long int i;
    for (i = 0; i < nedges; ++i) {
        std::sort(edges[i].begin(), edges[i].end() - 2);
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
    
    std::vector<std::vector<double>> finaledges(nedges / 2, std::vector<double>(4));
    for (i = 0; i < nedges / 2; ++i) {
        finaledges[i] = edges[2 * i];
    }
    
    return finaledges;
}

void PrintGraphs(std::ofstream& graphfile, std::ofstream& AAgraphfile, std::ofstream& BBgraphfile, std::ofstream& ABgraphfile, std::vector<std::vector<unsigned long int>> lattice, unsigned long int nnodes, int z, double eaa, double ebb, double eab) {
 // Writes .GraphGeod files.
    std::vector<double>* energies = new std::vector<double>(3);
    (*energies)[0] = -1 * eaa;
    (*energies)[1] = -1 * ebb;
    (*energies)[2] = -1 * eab;
    
    std::vector<std::vector<double>>* edges = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>>* AAedges = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>>* BBedges = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>>* ABedges = new std::vector<std::vector<double>>;
    std::vector<double>* edge = new std::vector<double>(4);
    unsigned long int adjid;
    
    unsigned long int i;
    int a;
    for (i = 0; i < nnodes; ++i) {
        (*edge)[0] = (double) (i + 1);
        for (a = 1; a < z + 1; ++a) {
            adjid = lattice[i][a];
            (*edge)[1] = (double) (adjid + 1);
            if (lattice[adjid][0] == lattice[i][0]) {
                if (lattice[i][0] == 1) {
                    (*edge)[2] = (*energies)[0];
                    (*edge)[3] = 1;
                    AAedges->push_back(*edge);
                }
                else {
                    (*edge)[2] = (*energies)[1];
                    (*edge)[3] = 2;
                    BBedges->push_back(*edge);
                }
            }
            else {
                (*edge)[2] = (*energies)[2];
                (*edge)[3] = 3;
                ABedges->push_back(*edge);
            }
            edges->push_back(*edge);
        }
    }
    delete energies;
    delete edge;
    
    *edges = EdgeSort(*edges, edges->size());
    for (i = 0; i < edges->size(); ++i) {
        graphfile << (*edges)[i][0] << " " << (*edges)[i][1];
        for (a = 0; a < 7; ++a) {
            graphfile << " 0";
        }
        graphfile << " " << (*edges)[i][2] << " " << (*edges)[i][3] << std::endl;
    }
    delete edges;
    *AAedges = EdgeSort(*AAedges, AAedges->size());
    for (i = 0; i < AAedges->size(); ++i) {
        AAgraphfile << (*AAedges)[i][0] << " " << (*AAedges)[i][1];
        for (a = 0; a < 7; ++a) {
            AAgraphfile << " 0";
        }
        AAgraphfile << " " << (*AAedges)[i][2] << " " << (*AAedges)[i][3] << std::endl;
    }
    delete AAedges;
    *BBedges = EdgeSort(*BBedges, BBedges->size());
    for (i = 0; i < BBedges->size(); ++i) {
        BBgraphfile << (*BBedges)[i][0] << " " << (*BBedges)[i][1];
        for (a = 0; a < 7; ++a) {
            BBgraphfile << " 0";
        }
        BBgraphfile << " " << (*BBedges)[i][2] << " " << (*BBedges)[i][3] << std::endl;
    }
    delete BBedges;
    *ABedges = EdgeSort(*ABedges, ABedges->size());
    for (i = 0; i < ABedges->size(); ++i) {
        ABgraphfile << (*ABedges)[i][0] << " " << (*ABedges)[i][1];
        for (a = 0; a < 7; ++a) {
            ABgraphfile << " 0";
        }
        ABgraphfile << " " << (*ABedges)[i][2] << " " << (*ABedges)[i][3] << std::endl;
    }
    delete ABedges;
}

int main(int argc, char** argv) {
 // Generates lattice and writes outputs.
    double temp, eaa, ebb, eab;
    eaa = ebb = 10000;
    eab = 9500;
    unsigned long int nsnaps, niters, xdim, ydim, na, nb;
    nsnaps = 1;
    niters = 100000;
    xdim = ydim = 10;
    na = nb = 50;
    int order, color;
    order = color = 0;
    
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
        else if ((*params)[i][0] == "entropy") {
            std::istringstream orderstr((*params)[i][1]);
            orderstr >> order;
            if (order < 0 || order > 3) {
                std::cerr << "If local entropy calculation is requested, order must be 1, 2, or 3. See README." << std::endl;
                return EXIT_FAILURE;
            }
        }
        else if ((*params)[i][0] == "color") {
            if ((*params)[i][1] == "y") {
                color = 1;
            }
            else if ((*params)[i][1] != "n") {
                std::cerr << "Invalid: please indicate whether or not you would like to print color-coded lattice(s) to std::cout (default value: no)." << std::endl;
            }
        }
        else {
            std::cerr << "Warning: option \"" << (*params)[i][0] << "\" not yet implemented." << std::endl;
        }
    }
    delete params;
    unsigned long int mindim;
    if (!order) {
        mindim = 1;
    }
    else if (order <= 2) {
        mindim = 3;
    }
    else {
        mindim = 5;
    }
    if (xdim < mindim || ydim < mindim) {
        if (mindim == 1) {
            std::cerr << "Invalid: lattice dimensions must be positive integers." << std::endl;
        }
        else {
            std::cerr << "Invalid: lattice dimensions must be at least " << mindim << " for local entropy order " << order << "." << std::endl;
        }
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
    
    int ncells, z;
    if (order <= 1) {
        ncells = 5;
    }
    else if (order == 2) {
        ncells = 9;
    }
    else {
        ncells = 13;
    }
    z = 4;
    double r, w;
    r = 8.31446;
    w = eab - 0.5 * eaa - 0.5 * ebb;
    std::vector<double>* probs = new std::vector<double>(z);
    
    int a;
    for (a = 0; a < z; ++a) {
        if (temp) {
            (*probs)[a] = exp(-((2 * ((double) a) + 2) * std::abs(w)) / (r * temp));
        }
        else {
            (*probs)[a] = 0;
        }
    }
    double redtemp;
    if (temp) {
        if (w) {
            redtemp = r * temp / w;
        }
        else {
            redtemp = std::numeric_limits<double>::infinity();
        }
    }
    else {
        redtemp = 0;
        nsnaps = niters = 1;
        std::cerr << "Warning: all possible moves have zero probability. Only initial conditions will be returned." << std::endl;
    }
    
    std::ofstream outfile("output.dat");
    outfile << "Interchange energy (w) = " << w << std::endl;
    outfile << "Reduced temperature (kT / w) = " << redtemp << std::endl;
    if (redtemp > 0) {
        outfile << "Probability of adding [x] A-B intercell interactions:" << std::endl;
    }
    else if (redtemp < 0) {
        outfile << "Probability of removing [x] A-B intercell interactions:" << std::endl;
    }
    for (a = 0; a < z; ++a) {
        outfile << " " << 2 * a + 2 << ": " << (*probs)[a] << std::endl;
    }
    
    unsigned long int nsmall;
    int minortype;
    if (na <= nb) {
        nsmall = na;
        minortype = 1;
    }
    else {
        nsmall = nb;
        minortype = 2;
    }
    std::vector<std::vector<unsigned long int>>* initlattice = new std::vector<std::vector<unsigned long int>>(nnodes, std::vector<unsigned long int>(z + 1));
    *initlattice = Lattice(nnodes, xdim, ydim, order, ncells, nsmall, minortype, gen, distnodes, distmax);
    unsigned long int initnab = NAB(*initlattice, nnodes, z);
    PrintLattice(outfile, *initlattice, initnab, xdim, ydim, 0, 0);
    if (color) {
        PrintLattice(std::cout, *initlattice, initnab, xdim, ydim, 0, color);
    }
    
    std::vector<std::vector<unsigned long int>>* lattice = new std::vector<std::vector<unsigned long int>>(nnodes, std::vector<unsigned long int>(ncells + 1));
    double prob, totenergy;
    unsigned long int nab, molid1, moltype1, adjcount1, molid2, moltype2, adjcount2, adjid;
    int nabdiff, adjint, probint;
    adjint = probint = 0;
    for (i = 0; i < nsnaps; ++i) {
        *lattice = *initlattice;
        nab = initnab;
        
        std::ostringstream energyfilestr;
        energyfilestr << "lattice" << i + 1 << ".energy";
        std::ofstream energyfile(energyfilestr.str());
        std::ostringstream graphfilestr;
        graphfilestr << "lattice" << i + 1 << ".GraphGeod";
        std::ofstream graphfile(graphfilestr.str());
        std::ostringstream AAgraphfilestr;
        AAgraphfilestr << "AAlattice" << i + 1 << ".GraphGeod";
        std::ofstream AAgraphfile(AAgraphfilestr.str());
        std::ostringstream BBgraphfilestr;
        BBgraphfilestr << "BBlattice" << i + 1 << ".GraphGeod";
        std::ofstream BBgraphfile(BBgraphfilestr.str());
        std::ostringstream ABgraphfilestr;
        ABgraphfilestr << "ABlattice" << i + 1 << ".GraphGeod";
        std::ofstream ABgraphfile(ABgraphfilestr.str());
        std::ostringstream xyzfilestr;
        xyzfilestr << "lattice" << i + 1 << ".xyz";
        std::ofstream xyzfile(xyzfilestr.str());
        std::ostringstream AAxyzfilestr;
        AAxyzfilestr << "AAlattice" << i + 1 << ".xyz";
        std::ofstream AAxyzfile(AAxyzfilestr.str());
        std::ostringstream BBxyzfilestr;
        BBxyzfilestr << "BBlattice" << i + 1 << ".xyz";
        std::ofstream BBxyzfile(BBxyzfilestr.str());
        std::ostringstream ABxyzfilestr;
        ABxyzfilestr << "ABlattice" << i + 1 << ".xyz";
        std::ofstream ABxyzfile(ABxyzfilestr.str());
        
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
            if ((nabdiff > 0 && w > 0) || (nabdiff < 0 && w < 0)) {
                prob = (*probs)[std::abs(nabdiff) / 2 - 1];
                if ((distmax(gen) / ((double) RAND_MAX)) <= prob) {
                    probint = 1;
                }
            }
            if ((nabdiff <= 0 && w > 0) || (nabdiff >= 0 && w < 0) || probint || !w) {
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
                if (temp) {
                    PrintLattice(outfile, *lattice, nab, xdim, ydim, i + 1, 0);
                    if (color) {
                        PrintLattice(std::cout, *lattice, nab, xdim, ydim, i + 1, color);
                    }
                    PrintXYZs(xyzfile, AAxyzfile, BBxyzfile, ABxyzfile, *lattice, nnodes, na, nb, xdim, ydim, z);
                }
                PrintGraphs(graphfile, AAgraphfile, BBgraphfile, ABgraphfile, *lattice, nnodes, z, eaa, ebb, eab);
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
