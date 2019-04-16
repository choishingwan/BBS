/*
 * genotype.h
 *
 *  Created on: 15 May 2018
 *      Author: shingwanchoi
 */

#ifndef GENOTYPE_H_
#define GENOTYPE_H_
#include "misc.hpp"
#include "plink_common.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct Sample
{
    std::string FID;
    std::string IID;
    bool x_var = false;
};

class SNP
{
public:
    SNP() {}
    SNP(const std::string& f, const std::streampos& b) : file(f), byte_pos(b) {}
    SNP(const std::string& f, const std::streampos& b, const double& m)
        : file(f), byte_pos(b), maf(m)
    {
    }
    SNP(const SNP& s) : file(s.file), byte_pos(s.byte_pos), maf(s.maf) {}
    std::string file;
    std::streampos byte_pos;
    double maf = 0.0;
    void set_maf(const double& input_maf) { maf = input_maf; }
    bool operator==(const SNP& s)
    {
        return (file.compare(s.file) == 0) && (byte_pos == s.byte_pos);
    }

    bool operator!=(const SNP& s)
    {
        return !((file.compare(s.file) == 0) && (byte_pos == s.byte_pos));
    }
};

class Genotype
{
public:
    Genotype(const std::string& prefix);
    virtual ~Genotype();

    void load_samples(const std::unordered_set<std::string>& sample_list,
                      const std::unordered_set<std::string>& related_list);
    void load_snps(const std::unordered_set<std::string>& snp_list, const std::vector<double> &effect,
                   const size_t num_selected, std::vector<double> &scores, const size_t seed, const bool standardize);
    size_t sample_size() const { return m_sample_ct; }
    std::string name(size_t i_sample)
    {
        auto&& sample = m_sample_names.at(i_sample);
        return std::string(sample.FID + "\t" + sample.IID);
    }
    bool sample_xvar(size_t i_sample)
    {
        return m_sample_names.at(i_sample).x_var;
    }

private:
    std::vector<std::string> set_genotype_files(const std::string& prefix);
    std::vector<std::string> m_genotype_files;


    std::vector<Sample>
    gen_sample_vector(const std::unordered_set<std::string>& sample_list,
                      const std::unordered_set<std::string>& related_list);
    std::vector<Sample> m_sample_names;
    uintptr_t m_unfiltered_sample_ct = 0;
    uintptr_t m_sample_ct = 0;
    std::vector<uintptr_t> m_founder_info;
    std::vector<uintptr_t> m_sample_include;
    std::vector<uintptr_t> m_tmp_genotype;
    size_t m_num_male = 0;
    size_t m_num_female = 0;
    size_t m_num_ambig_sex = 0;
    size_t m_num_non_founder = 0;
    size_t m_num_unrelated = 0;
    uintptr_t m_founder_ct = 0;
    void get_xbeta(const std::vector<double>& effect,std::vector<double>& scores, const bool standardize);
    void check_bed(const std::string& bed_name, const size_t num_marker);
    std::vector<SNP>
    gen_snp_vector(const std::unordered_set<std::string>& snp_list);
    inline uintptr_t get_final_mask(uint32_t sample_ct)
    {
        uint32_t uii = sample_ct % BITCT2;
        if (uii) {
            return (ONELU << (2 * uii)) - ONELU;
        }
        else
        {
            return ~ZEROLU;
        }
    }
    uintptr_t m_bed_offset = 3;
    std::vector<SNP> m_existed_snps;


    uint32_t load_and_collapse_incl(uint32_t unfiltered_sample_ct,
                                    uint32_t sample_ct,
                                    const uintptr_t* __restrict sample_include,
                                    uintptr_t final_mask,
                                    std::ifstream& bedfile,
                                    uintptr_t* __restrict rawbuf,
                                    uintptr_t* __restrict mainbuf)
    {
        assert(unfiltered_sample_ct);
        uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
        // if we don't perform selection, we can directly perform the read on
        // the mainbuf
        if (unfiltered_sample_ct == sample_ct) {
            rawbuf = mainbuf;
        }
        // we try to read in the data and store it in rawbug
        if (!bedfile.read((char*) rawbuf, unfiltered_sample_ct4)) {
            return RET_READ_FAIL;
        }
        if (unfiltered_sample_ct != sample_ct) {
            // if we need to perform selection, we will remove all unwanted
            // sample and push the data forward
            copy_quaterarr_nonempty_subset(rawbuf, sample_include,
                                           unfiltered_sample_ct, sample_ct,
                                           mainbuf);
        }
        else
        {
            // if we dno't need filtering, then we simply mask out the unwanted
            // region (to avoid the leftover, if any)
            mainbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
        }
        return 0;
    }

};

#endif /* GENOTYPE_H_ */
