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
    SNP(const std::string& f, const std::string& in, const std::streampos& b)
        : file(f), name(in), byte_pos(b)
    {
    }
    SNP(const std::string& f, const std::streampos& b, const double& m)
        : file(f), byte_pos(b)
    {
    }
    SNP(const std::string& f, const std::string& in, const size_t i)
        : file(f), name(in), byte_pos(i)
    {
    }
    SNP(const SNP& s) : file(s.file), name(s.name), byte_pos(s.byte_pos) {}
    void set_pos(const std::streampos& b) { byte_pos = b; }
    bool operator==(const SNP& s)
    {
        return (file.compare(s.file) == 0) && (byte_pos == s.byte_pos);
    }

    bool operator!=(const SNP& s)
    {
        return !((file.compare(s.file) == 0) && (byte_pos == s.byte_pos));
    }
    void set_name(const std::string& in_name) { name = in_name; }
    std::string get_name() const { return name; }


    std::string file;
    std::string name;
    std::streampos byte_pos;
};

class Genotype
{
public:
    Genotype(const std::string& prefix);
    virtual ~Genotype();

    void load_samples(const std::unordered_set<std::string>& sample_list,
                      const std::unordered_set<std::string>& related_list);
    void load_snps(const std::unordered_set<std::string>& snp_list,
                   const size_t num_selected, const size_t seed);
    size_t sample_size() const { return m_sample_ct; }
    std::vector<SNP>
    gen_snp_vector(const std::unordered_set<std::string>& snp_list,
                   const size_t num_selected, const size_t seed);
    void get_xbeta(std::vector<double>& score, std::vector<double>& effect,
                   const bool standardize, const std::string& out)
    {
        assert(score.size() == m_sample_ct);
        const uintptr_t final_mask =
            get_final_mask(static_cast<uint32_t>(m_sample_ct));
        // this is use for initialize the array sizes
        const uintptr_t unfiltered_sample_ctl =
            BITCT_TO_WORDCT(m_unfiltered_sample_ct);
        const uintptr_t unfiltered_sample_ct4 =
            (m_unfiltered_sample_ct + 3) / 4;
        const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;

        std::ifstream bed_file;
        std::vector<uintptr_t> genotype_byte(unfiltered_sample_ctl * 2, 0);
        std::vector<uintptr_t> sample_include2(unfiltered_sample_ctv2);
        std::vector<uintptr_t> founder_include2(unfiltered_sample_ctv2);
        // fill it with the required mask (copy from PLINK2)
        init_quaterarr_from_bitarr(m_sample_include.data(),
                                   m_unfiltered_sample_ct,
                                   sample_include2.data());
        init_quaterarr_from_bitarr(m_founder_info.data(),
                                   m_unfiltered_sample_ct,
                                   founder_include2.data());
        std::string prev_file = "";
        std::streampos prev_loc = 0;
        std::ofstream output;
        uint32_t uii = 0;
        uintptr_t ulii = 0;
        uint32_t ujj;
        uint32_t ukk;
        uint32_t sample_idx = 0;
        uintptr_t* lbptr;
        double eff;
        double maf;
        double var = 1.0;
        double mean = 0.0;
        uint32_t ll_ct, lh_ct, hh_ct;
        uint32_t ref_founder_ct, het_founder_ct, alt_founder_ct;
        output.open(std::string(out + ".eff").c_str());
        if (!output.is_open())
        {
            throw std::runtime_error(
                std::string("ERROR: Cannot open file: " + out + ".eff"));
        }
        int num_completed = 0;
        double prev_completed = 0;
        double total_snp = static_cast<double>(m_existed_snps.size());

        fprintf(stderr, "\rProcessing %03.2f%%",
                num_completed / total_snp * 100);
        // loop through the SNPs. Use PLINK's founder MAF calculation script
        size_t idx = 0;
        for (auto&& snp : m_existed_snps)
        {
            if (prev_file != snp.file)
            {
                if (bed_file.is_open()) { bed_file.close(); }
                std::string cur_file = snp.file + ".bed";
                bed_file.open(cur_file.c_str(), std::ios::binary);
                if (!bed_file.is_open())
                {
                    throw std::runtime_error("Error: Cannot open bed file: "
                                             + cur_file);
                }
                prev_file = snp.file;
                prev_loc = 0;
            }

            if (prev_loc != snp.byte_pos
                && !bed_file.seekg(snp.byte_pos, std::ios_base::beg))
            { throw std::runtime_error("Error: Cannot read the bed file!"); }

            // loadbuf_raw is the temporary
            if (load_raw(unfiltered_sample_ct4, bed_file,
                         m_tmp_genotype.data()))
            {
                throw std::runtime_error(
                    "Error: Cannot read the bed file(read): " + prev_file);
            }
            single_marker_freqs_and_hwe(
                unfiltered_sample_ctv2, m_tmp_genotype.data(),
                sample_include2.data(), founder_include2.data(), m_sample_ct,
                &ll_ct, &lh_ct, &hh_ct, m_founder_ct, &ref_founder_ct,
                &het_founder_ct, &alt_founder_ct);
            // got the maf, now transform genotype to the require format for PRS
            if (m_unfiltered_sample_ct != m_sample_ct)
            {
                copy_quaterarr_nonempty_subset(
                    m_tmp_genotype.data(), m_sample_include.data(),
                    static_cast<uint32_t>(m_unfiltered_sample_ct),
                    static_cast<uint32_t>(m_sample_ct), genotype_byte.data());
            }
            else
            {
                genotype_byte = m_tmp_genotype;
                genotype_byte[(m_unfiltered_sample_ct - 1) / BITCT2] &=
                    final_mask;
            }

            prev_loc = static_cast<std::streampos>(unfiltered_sample_ct4)
                       + snp.byte_pos;

            eff = effect[idx];
            output << snp.name << "\t" << eff << std::endl;
            // get_score(score, genotype_byte, eff, standardize);
            maf = (het_founder_ct + alt_founder_ct * 2.0)
                  / (2.0 * (ref_founder_ct + het_founder_ct + alt_founder_ct));
            var = 1.0;
            mean = 0.0;
            double miss_dose = 2.0 * maf;
            if (standardize)
            {
                mean = miss_dose;
                var = (sqrt(2.0 * maf * (1.0 - maf)));
            }
            uii = 0;
            ujj = 0;
            // now start calculating the score
            lbptr = genotype_byte.data();
            eff /= var;
            do
            {
                // ulii contain the numeric representation of the current
                // genotype
                ulii = ~(*lbptr++);
                if (uii + BITCT2 > m_unfiltered_sample_ct)
                {
                    // this is PLINK, not sure exactly what this is about
                    ulii &=
                        (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
                }
                // ujj sample index of the current genotype block
                ujj = 0;
                while (ujj < BITCT)
                {
                    // go through the whole genotype block
                    // ukk is the current genotype
                    sample_idx = uii + (ujj / 2);
                    if (sample_idx >= m_sample_ct) { break; }
                    ukk = (ulii >> ujj) & 3;
                    // now we will get all genotypes (0, 1, 2, 3)
                    switch (ukk)
                    {
                    default: score[sample_idx] -= eff * mean; break;
                    case 1: score[sample_idx] += eff * (1.0 - mean); break;
                    case 3: score[sample_idx] += eff * (2.0 - mean); break;
                    case 2:
                        score[sample_idx] += eff * (miss_dose - mean);
                        break;
                    }

                    // ulii &= ~((3 * ONELU) << ujj);
                    // as each sample is represented by two byte, we will add 2
                    // to the index
                    ujj += 2;
                }
                // uii is the number of samples we have finished so far
                uii += BITCT2;
            } while (uii < m_sample_ct);

            if (num_completed / total_snp - prev_completed > 0.01)
            {
                fprintf(stderr, "\rProcessing %03.2f%%",
                        num_completed / total_snp * 100);
                prev_completed = num_completed / total_snp;
            }
            ++num_completed;
            idx++;
        }
        fprintf(stderr, "\rProcessing 100.0%%\n");
    }
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
    void get_maf();
    void check_bed(const std::string& bed_name, const size_t num_marker);
    inline uintptr_t get_final_mask(uint32_t sample_ct)
    {
        uint32_t uii = sample_ct % BITCT2;
        if (uii) { return (ONELU << (2 * uii)) - ONELU; }
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
                                    uintptr_t final_mask, uint32_t do_reverse,
                                    std::ifstream& bedfile,
                                    uintptr_t* __restrict rawbuf,
                                    uintptr_t* __restrict mainbuf)
    {
        assert(unfiltered_sample_ct);
        uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
        // if we don't perform selection, we can directly perform the read on
        // the mainbuf
        if (unfiltered_sample_ct == sample_ct) { rawbuf = mainbuf; }
        // we try to read in the data and store it in rawbug
        if (!bedfile.read((char*) rawbuf, unfiltered_sample_ct4))
        { return RET_READ_FAIL; }
        if (unfiltered_sample_ct != sample_ct)
        {
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
        if (do_reverse)
        {
            // this will never be callsed in PRSice
            reverse_loadbuf(sample_ct, (unsigned char*) mainbuf);
        }
        return 0;
    }


    inline uint32_t load_raw(uintptr_t unfiltered_sample_ct4,
                             std::ifstream& bedfile, uintptr_t* rawbuf)
    {
        // only use this if all accesses to the data involve
        // 1. some sort of mask, or
        // 2. explicit iteration from 0..(unfiltered_sample_ct-1).
        // otherwise improper trailing bits might cause a segfault, when we
        // should be ignoring them or just issuing a warning.
        if (!bedfile.read((char*) rawbuf, unfiltered_sample_ct4))
        { return RET_READ_FAIL; }
        return 0;
    }

    // modified version of the
    // single_marker_freqs_and_hwe function from PLINK (plink_filter.c)
    // we remove the HWE calculation (we don't want to include that yet,
    // as that'd require us to implement a version for bgen)
    inline void single_marker_freqs_and_hwe(
        uintptr_t unfiltered_sample_ctl2, uintptr_t* lptr,
        uintptr_t* sample_include2, uintptr_t* founder_include2,
        uintptr_t sample_ct, uint32_t* ll_ctp, uint32_t* lh_ctp,
        uint32_t* hh_ctp, uintptr_t sample_f_ct, uint32_t* ll_ctfp,
        uint32_t* lh_ctfp, uint32_t* hh_ctfp)
    {
        uint32_t tot_a = 0;
        uint32_t tot_b = 0;
        uint32_t tot_c = 0;
        uint32_t tot_a_f = 0;
        uint32_t tot_b_f = 0;
        uint32_t tot_c_f = 0;
        uintptr_t* lptr_end = &(lptr[unfiltered_sample_ctl2]);
        uintptr_t loader;
        uintptr_t loader2;
        uintptr_t loader3;
#ifdef __LP64__
        uintptr_t cur_decr = 120;
        uintptr_t* lptr_12x_end;
        unfiltered_sample_ctl2 -= unfiltered_sample_ctl2 % 12;
        while (unfiltered_sample_ctl2 >= 120)
        {
        single_marker_freqs_and_hwe_loop:
            lptr_12x_end = &(lptr[cur_decr]);
            count_3freq_1920b((__m128i*) lptr, (__m128i*) lptr_12x_end,
                              (__m128i*) sample_include2, &tot_a, &tot_b,
                              &tot_c);
            count_3freq_1920b((__m128i*) lptr, (__m128i*) lptr_12x_end,
                              (__m128i*) founder_include2, &tot_a_f, &tot_b_f,
                              &tot_c_f);
            lptr = lptr_12x_end;
            sample_include2 = &(sample_include2[cur_decr]);
            founder_include2 = &(founder_include2[cur_decr]);
            unfiltered_sample_ctl2 -= cur_decr;
        }
        if (unfiltered_sample_ctl2)
        {
            cur_decr = unfiltered_sample_ctl2;
            goto single_marker_freqs_and_hwe_loop;
        }
#else
        uintptr_t* lptr_twelve_end =
            &(lptr[unfiltered_sample_ctl2 - unfiltered_sample_ctl2 % 12]);
        while (lptr < lptr_twelve_end)
        {
            count_3freq_48b(lptr, sample_include2, &tot_a, &tot_b, &tot_c);
            count_3freq_48b(lptr, founder_include2, &tot_a_f, &tot_b_f,
                            &tot_c_f);
            lptr = &(lptr[12]);
            sample_include2 = &(sample_include2[12]);
            founder_include2 = &(founder_include2[12]);
        }
#endif
        while (lptr < lptr_end)
        {
            loader = *lptr++;
            loader2 = *sample_include2++;
            loader3 = (loader >> 1) & loader2;
            loader2 &= loader;
            // N.B. because of the construction of sample_include2, only
            // even-numbered bits can be present here.  So popcount2_long is
            // safe.
            tot_a += popcount2_long(loader2);
            tot_b += popcount2_long(loader3);
            tot_c += popcount2_long(loader & loader3);
            loader2 = *founder_include2++;
            loader3 = (loader >> 1) & loader2;
            loader2 &= loader;
            tot_a_f += popcount2_long(loader2);
            tot_b_f += popcount2_long(loader3);
            tot_c_f += popcount2_long(loader & loader3);
        }
        *hh_ctp = tot_c;
        *lh_ctp = tot_b - tot_c;
        *ll_ctp = sample_ct - tot_a - *lh_ctp;
        *hh_ctfp = tot_c_f;
        *lh_ctfp = tot_b_f - tot_c_f;
        *ll_ctfp = sample_f_ct - tot_a_f - *lh_ctfp;
    }
};

#endif /* GENOTYPE_H_ */
