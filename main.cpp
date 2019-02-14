/*
 * main.cpp
 *
 *  Created on: 15 May 2018
 *      Author: shingwanchoi
 *      Biobank Simulation
 */

#include "genotype.h"
#include "misc.hpp"
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <unordered_set>
#include <vector>
#include <cmath>

void usage()
{

    std::cerr << "Biobank Simulation Tool\n";
    std::cerr << "Usage: BBS [options]\n";
    std::cerr << "Options:\n";
    std::cerr << "    --input   | -i    Input file prefix\n";
    std::cerr << "    --out     | -o    Output file prefix\n";
    std::cerr << "    --num     | -n    Number of Causal SNPs\n";
    std::cerr << "    --effect  | -e    Effect size distribution.\n";
    std::cerr << "                      0 for exponential\n";
    std::cerr << "                      1 for Chi-Square\n";
    std::cerr << "                      2 for Normal Distribution\n";
    std::cerr << "    --fix     | -f    Fixed effect\n";
    std::cerr << "    --extract | -E    List of SNPs to extract.\n";
    std::cerr << "    --keep    | -k    List of samples to keep\n";
    std::cerr << "    --thread  | -t    Number of thread used\n";
    std::cerr << "    --herit   | -H    List of heritability to construct\n";
    std::cerr << "    --seed    | -s    Seed for random generator\n";
    std::cerr << "    --help    | -h    Display this help message\n";
}


std::unordered_set<std::string> extract_ref(const std::string& extract_name,
                                            const size_t index)
{
    std::unordered_set<std::string> res;
    if (extract_name.empty()) {
        return res;
    }
    std::ifstream in;
    in.open(extract_name.c_str());
    if (!in.is_open()) {
        throw std::runtime_error(
            std::string("ERROR: Cannot open file: " + extract_name));
    }
    std::string line;
    while (std::getline(in, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if (token.size() <= index) {
            throw std::runtime_error(
                std::string("ERROR: Number of column is less than expected!"));
        }
        res.insert(token[index]);
    }
    in.close();
    return res;
}

int main(int argc, char* argv[])
{
    if (argc <= 1) {
        usage();
        fprintf(stderr, "Please provide the required parameters\n");
        exit(-1);
    }
    static const char* optString = "i:s:o:f:e:k:n:t:E:xH:h?";
    static const struct option longOpts[] = {
        {"input", required_argument, nullptr, 'i'},
        {"seed", required_argument, nullptr, 's'},
        {"out", required_argument, nullptr, 'o'},
        {"extract", required_argument, nullptr, 'E'},
        {"effect", required_argument, nullptr, 'e'},
        {"keep", required_argument, nullptr, 'k'},
        {"nsnp", required_argument, nullptr, 'n'},
        {"fix", required_argument, nullptr, 'f'},
        {"thread", required_argument, nullptr, 't'},
        {"standardize", no_argument, nullptr, 'x'},
        {"herit", no_argument, nullptr, 'H'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}};

    std::string prefix, out = "Out", extract, keep, herit;
    size_t seed = std::random_device()(), num_snp=0, thread = 1, effect = 0;
    double fixed_effect = 0.0;
    bool use_fixed = false, standardize = false;


    int longIndex = 0;
    int opt = 0;
    std::string command = "";
    opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    std::string error_message = "";
    // Start reading all the parameters and perform the qc at the same time
    while (opt != -1) {
        switch (opt)
        {
        case 'i': prefix = optarg; break;
        case 'o': out = optarg; break;
        case 'f':
            try{
                fixed_effect = misc::convert<double>(optarg);
                use_fixed = true;
            }catch(...){
                fprintf(stderr, "ERROR: Effect size must be numeric!\n");
            }
            break;
        case 's': seed = misc::convert<size_t>(optarg); break;
        case 'E': extract = optarg; break;
        case 'k': keep = optarg; break;
        case 'e': effect = misc::convert<size_t>(optarg); break;
        case 't': thread = misc::convert<size_t>(optarg); break;
        case 'x': standardize = true; break;
        case 'H': herit = optarg; break;
        case 'n':
        {
            try
            {
                num_snp = misc::convert<size_t>(optarg);
                if (num_snp == 0) {
                    fprintf(stderr, "0 SNPs required, programme terminate\n");
                    exit(0);
                }
            }
            catch (const std::runtime_error&)
            {
                fprintf(stderr, "ERROR: Number of SNPs must be numeric!\n");
            }
        }
        break;
        case 'h':
        case '?':
            usage();
            exit(0);
        default:
            throw "Undefined operator, please use --help for more information!";
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    std::vector<double> heritability;
    std::vector<std::string> token = misc::split(herit, ",");
    for(auto &&t :token){
            double h = misc::convert<double>(t);
            if(h  < 0.0 | h > 1.0){
                std::cerr << "Error: Heritability must be within 0 and 1" << std::endl;
                exit(-1);
            }
            heritability.push_back(h);
    }
    std::unordered_set<std::string> snp_list = extract_ref(extract, 0);
    std::unordered_set<std::string> sample_list = extract_ref(keep, 1);

    Genotype geno(prefix);
    geno.load_samples(sample_list);
    geno.load_snps(snp_list, num_snp, seed);
    std::vector<double> score(geno.sample_size(), 0.0);
    std::normal_distribution<double> norm_dist(0, 1);
    std::chi_squared_distribution<double> chi_dist(1);
    std::exponential_distribution<double> exp_dis(1);
    std::cerr << "Start generating X Beta" << std::endl;
    if (use_fixed) {
        // use fixed effect
        geno.get_xbeta(score, fixed_effect, standardize);
    }
    else
    {
        switch (effect)
        {
        case 0:
            geno.get_xbeta<std::exponential_distribution<double>>(
                score, exp_dis, standardize, seed, out);
            break;
        case 1:
            geno.get_xbeta<std::chi_squared_distribution<double>>(
                score, chi_dist, standardize, seed, out);
            break;
        case 2:
            geno.get_xbeta<std::normal_distribution<double>>(score, norm_dist,
                                                             standardize, seed, out);
            break;
        }
    }
// here we've got the XB stored in the score items we can then generate the desired phenotypes

    //std::vector<double> pheno(geno.sample_size(), 0.0);
    std::cerr << "Calculate Variane of X Beta" << std::endl;
    misc::RunningStat rs;
    for(auto &&s : score){
            rs.push(s);
    }
    double varXB = rs.var();
    std::mt19937 g(seed);

    std::ofstream output;
    output.open(out.c_str());
    if(!output.is_open()){
        std::cerr << "ERROR: Cannot open output file: " << out << std::endl;
        exit(-1);
    }

    std::vector<std::function< void() >> rand_dist;

    output << "FID\tIID";
    misc::vec2d<double> error_values(heritability.size(), geno.sample_size());
    size_t i_h = 0;
    std::cerr << "Generating error term for samples" << std::endl;
    for(auto &&h : heritability){
        std::cerr << "Heritability of " << h << std::endl;
        output << "\th_" <<h;
        if(h==0.0){
            auto rand = std::bind(std::normal_distribution<double>(0,1), g);
            for(size_t i_sample=0; i_sample < geno.sample_size(); ++i_sample){
                error_values(i_h, i_sample) =rand();
            }
        }else{
            double error_var = varXB*(1.0-h)/h;
            auto rand=std::bind(std::normal_distribution<double>(0,std::sqrt(error_var)), g);
            for(size_t i_sample=0; i_sample < geno.sample_size(); ++i_sample){
                error_values(i_h, i_sample) =rand();
            }
        }
        i_h++;
    }
    output << std::endl;
    std::cerr << "Generate the phenotype" << std::endl;
    for(size_t i_sample=0; i_sample < geno.sample_size(); ++i_sample){
        output << geno.name(i_sample);
        for(i_h=0; i_h < heritability.size(); ++i_h)
            if(misc::logically_equal(heritability[i_h],0.0)){
                output << "\t" << error_values(i_h, i_sample);
            }else{
                output << "\t" << error_values(i_h, i_sample)+score[i_sample];
            }
        output << std::endl;
    }
    output.close();

}