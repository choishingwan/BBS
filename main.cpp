
#include "genotype.h"
#include "misc.hpp"
#include <cmath>
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

void usage()
{
    fprintf(stderr, "Biobank Simulation Tool\n");
    fprintf(stderr, "Version 1.2 (21th July, 2020)\n");
    fprintf(stderr, "Usage: BBS [options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    --input   | -i    Input file prefix\n");
    fprintf(stderr, "    --out     | -o    Output file prefix\n");
    fprintf(stderr, "    --nsnp    | -n    Number of Causal SNPs\n");
    fprintf(stderr, "    --effect  | -e    Effect size distribution.\n");
    fprintf(stderr, "                      0 for exponential\n");
    fprintf(stderr, "                      1 for Chi-Square\n");
    fprintf(stderr, "                      2 for Normal Distribution\n");
    fprintf(stderr,
            "    --fix     | -f    Fixed all effect size to this number\n");
    fprintf(stderr, "    --file    | -F    Use effect size in this file\n");
    fprintf(stderr,
            "    --rand    | -r    Fixed all effect size to a random number\n");
    fprintf(stderr, "    --std     | -d    Standardize the genotype\n");
    fprintf(stderr, "    --extract | -E    List of SNPs to extract.\n");
    fprintf(stderr, "    --keep    | -k    List of samples to keep\n");
    fprintf(stderr, "    --x-var   | -x    List of samples to exclude from \n");
    fprintf(stderr, "                      MAF and Var(XB) calculation\n");
    fprintf(stderr, "    --thread  | -t    Number of thread used\n");
    fprintf(stderr,
            "    --herit   | -H    List of heritability to construct\n");
    fprintf(stderr, "    --seed    | -s    Seed for random generator\n");
    fprintf(stderr, "    --help    | -h    Display this help message\n");
}


std::unordered_set<std::string>
extract_ref(std::unique_ptr<std::istream> extract_file, const size_t index)
{
    std::unordered_set<std::string> res;
    std::vector<std::string> token;
    std::string line;
    while (std::getline(*extract_file, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() <= index)
        {
            throw std::runtime_error(
                "ERROR: Number of column is less than expected!");
        }
        res.insert(token[index]);
    }
    extract_file.reset();
    return res;
}

template <typename T>
std::vector<double> generate_data(size_t size, T rand, std::mt19937 g)
{
    auto effect = std::bind(rand, g);
    std::vector<double> data(size);
    std::generate(data.begin(), data.end(), [&effect]() { return effect(); });
    return data;
}

int main(int argc, char* argv[])
{
    if (argc <= 1)
    {
        usage();
        fprintf(stderr, "Please provide the required parameters\n");
        exit(-1);
    }
    static const char* optString = "i:s:o:f:e:k:n:t:E:F:x:rdH:h?";
    static const struct option longOpts[] = {
        {"input", required_argument, nullptr, 'i'},
        {"seed", required_argument, nullptr, 's'},
        {"out", required_argument, nullptr, 'o'},
        {"extract", required_argument, nullptr, 'E'},
        {"effect", required_argument, nullptr, 'e'},
        {"keep", required_argument, nullptr, 'k'},
        {"nsnp", required_argument, nullptr, 'n'},
        {"fix", required_argument, nullptr, 'f'},
        {"file", required_argument, nullptr, 'F'},
        {"rand-fix", no_argument, nullptr, 'r'},
        {"thread", required_argument, nullptr, 't'},
        {"std", no_argument, nullptr, 'd'},
        {"x-var", required_argument, nullptr, 'x'},
        {"herit", required_argument, nullptr, 'H'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}};

    std::string prefix, out = "Out", extract, keep, herit, xvar, fix_file;
    size_t seed = std::random_device()(), num_snp = 0, thread = 1, effect = 0;
    double fixed_effect = 0.0;
    bool use_fixed = false, standardize = false, use_rand_fixed = false;

    int longIndex = 0;
    int opt = 0;
    std::string command = "";
    opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    std::string error_message = "";
    // Start reading all the parameters and perform the qc at the same time
    while (opt != -1)
    {
        switch (opt)
        {
        case 'i': prefix = optarg; break;
        case 'o': out = optarg; break;
        case 'f':
            try
            {
                fixed_effect = misc::convert<double>(optarg);
                use_fixed = true;
            }
            catch (...)
            {
                fprintf(stderr, "ERROR: Effect size must be numeric!\n");
            }
            break;
        case 'F': fix_file = optarg; break;
        case 'r': use_rand_fixed = true; break;
        case 's': seed = misc::convert<size_t>(optarg); break;
        case 'E': extract = optarg; break;
        case 'k': keep = optarg; break;
        case 'e': effect = misc::convert<size_t>(optarg); break;
        case 't': thread = misc::convert<size_t>(optarg); break;
        case 'd': standardize = true; break;
        case 'x': xvar = optarg; break;
        case 'H': herit = optarg; break;
        case 'n':
        {
            try
            {
                num_snp = misc::convert<size_t>(optarg);
                if (num_snp == 0)
                {
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
        case '?': usage(); exit(0);
        default:
            throw "Undefined operator, please use --help for more information!";
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    std::cerr << "Seed: " << seed << std::endl;
    std::vector<double> heritability;
    std::vector<std::string> token = misc::split(herit, ",");
    for (auto&& t : token)
    {
        double h = misc::convert<double>(t);
        if (h<0.0 | h> 1.0)
        {
            std::cerr << "Error: Heritability must be within 0 and 1"
                      << std::endl;
            std::cerr << "       Observed: " << h << std::endl;
            exit(-1);
        }
        heritability.push_back(h);
    }
    if (heritability.empty())
    {
        std::cerr << "Error: Must provide at least one heritability information"
                  << std::endl;
        exit(-1);
    }
    std::unordered_set<std::string> snp_list;
    std::unordered_set<std::string> sample_list;
    std::unordered_set<std::string> no_varx_list;
    if (!extract.empty())
    {
        auto extract_file = misc::load_stream(extract);
        snp_list = extract_ref(std::move(extract_file), 0);
        std::cerr << "Keeping " << snp_list.size() << " SNPs" << std::endl;
    }
    if (!keep.empty())
    {
        auto keep_file = misc::load_stream(keep);
        sample_list = extract_ref(std::move(keep_file), 1);
    }
    if (!xvar.empty())
    {
        auto xvar_file = misc::load_stream(xvar);
        no_varx_list = extract_ref(std::move(xvar_file), 1);
    }
    // relatedness file should contain two column, ID1 and ID2, which
    // represents the related pair. Here, we will always exclude
    // sample in the second column from the variance calculation
    // Read in the PLINK file information
    std::mt19937 g(seed);
    std::normal_distribution<double> norm_dist(0, 1);
    std::chi_squared_distribution<double> chi_dist(1);
    std::exponential_distribution<double> exp_dis(1);
    if (use_rand_fixed) { fixed_effect = norm_dist(g); }

    Genotype geno(prefix);

    geno.load_samples(sample_list, no_varx_list);
    geno.load_snps(snp_list, num_snp, seed);
    std::vector<double> score(geno.sample_size(), 0.0);
    std::cerr << "Start generating X Beta" << std::endl;
    std::vector<double> effect_sizes;
    if (use_fixed)
    {
        // use fixed effect
        effect_sizes.resize(num_snp);
        std::fill(effect_sizes.begin(), effect_sizes.end(), fixed_effect);
    }
    else if (!fix_file.empty())
    {
        geno.order_effects(fix_file, effect_sizes);
    }
    else
    {
        switch (effect)
        {
        case 0:
            effect_sizes = generate_data<std::exponential_distribution<double>>(
                num_snp, exp_dis, g);
            break;
        case 1:
            effect_sizes = generate_data<std::chi_squared_distribution<double>>(
                num_snp, chi_dist, g);
            break;
        case 2:
            effect_sizes = generate_data<std::normal_distribution<double>>(
                num_snp, norm_dist, g);
            break;
        }
    }
    geno.get_xbeta(score, effect_sizes, standardize, out);
    std::cerr << score[0] << std::endl;
    // here we've got the XB stored in the score items we can then generate the
    // desired phenotypes

    // std::vector<double> pheno(geno.sample_size(), 0.0);
    std::cerr << "Calculate Variance of X Beta" << std::endl;
    misc::RunningStat rs;
    assert(geno.sample_size() == score.size());
    for (size_t i = 0; i < score.size(); ++i)
    {
        if (!geno.sample_xvar(i)) { rs.push(score[i]); }
    }
    double varXB = rs.var();

    std::ofstream output;
    std::ofstream out_xb;
    output.open(out.c_str());
    out_xb.open(std::string(out + ".xbeta").c_str());
    if (!output.is_open())
    {
        std::cerr << "ERROR: Cannot open output file: " << out << std::endl;
        exit(-1);
    }
    if (!out_xb.is_open())
    {
        std::cerr << "ERROR: Cannot open output file: " << out << ".xbeta"
                  << std::endl;
        exit(-1);
    }

    std::vector<std::function<void()>> rand_dist;

    output << "FID\tIID";
    out_xb << "FID\tIID";
    misc::vec2d<double> error_values(heritability.size(), geno.sample_size());
    size_t i_h = 0;
    std::cerr << "Generating error term for samples" << std::endl;
    for (auto&& h : heritability)
    {
        std::cerr << "Heritability of " << h << std::endl;
        output << "\th_" << h;
        out_xb << "\th_" << h;
        if (h == 0.0)
        {
            auto rand = std::bind(std::normal_distribution<double>(0, 1), g);
            for (size_t i_sample = 0; i_sample < geno.sample_size(); ++i_sample)
            { error_values(i_h, i_sample) = rand(); }
        }
        else
        {
            double error_var = varXB * (1.0 - h) / h;
            auto rand = std::bind(
                std::normal_distribution<double>(0, std::sqrt(error_var)), g);
            for (size_t i_sample = 0; i_sample < geno.sample_size(); ++i_sample)
            { error_values(i_h, i_sample) = rand(); }
        }
        i_h++;
    }
    output << std::endl;
    out_xb << std::endl;
    std::cerr << "Generate the phenotype" << std::endl;
    for (size_t i_sample = 0; i_sample < geno.sample_size(); ++i_sample)
    {
        output << geno.name(i_sample);
        out_xb << geno.name(i_sample);
        for (i_h = 0; i_h < heritability.size(); ++i_h)
            if (misc::logically_equal(heritability[i_h], 0.0))
            {
                output << "\t" << error_values(i_h, i_sample);
                out_xb << "\t0";
            }
            else
            {
                output << "\t" << error_values(i_h, i_sample) + score[i_sample];
                out_xb << "\t" << score[i_sample];
            }
        output << std::endl;
        out_xb << std::endl;
    }
    output.close();
    out_xb.close();
}
