/*
 * genotype.h
 *
 *  Created on: 15 May 2018
 *      Author: shingwanchoi
 */

#ifndef GENOTYPE_H_
#define GENOTYPE_H_
#include "misc.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifdef __LP64__
#include <emmintrin.h>
#define ZEROLU 0LLU
#define ONELU 1LLU
#define VEC_BYTES 16
#define BITCT 64
#else // not __LP64__
#define ZEROLU 0LU
#define ONELU 1LU
#define VEC_BYTES 8
#define BITCT 32
#endif // __LP64__

#define BITCT_TO_WORDCT(val) (((val) + (BITCT - 1)) / BITCT)
#define SET_BIT(idx, arr) ((arr)[(idx) / BITCT] |= ONELU << ((idx) % BITCT))
#define BITCT2 (BITCT / 2)
#define CTZLU __builtin_ctzl

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
    SNP(const SNP& s) : file(s.file), byte_pos(s.byte_pos) {}
    std::string file;
    std::streampos byte_pos;
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
    void load_snps(const std::unordered_set<std::string>& snp_list,
                   const size_t num_selected, const size_t seed);
    size_t sample_size() const { return m_sample_ct; }
    void get_xbeta(std::vector<double>& score, double fixed_effect,
                   bool standardize);

    template <typename T>
    void get_xbeta(std::vector<double>& score, T rand, const bool standardize,
                   const size_t seed, const std::string& out)
    {
        assert(score.size() == m_sample_ct);
        std::mt19937 g(seed);
        auto effect = std::bind(rand, g);

        const uintptr_t final_mask = get_final_mask(m_sample_ct);
        // for array size
        const uintptr_t unfiltered_sample_ctl =
            BITCT_TO_WORDCT(m_unfiltered_sample_ct);
        std::ifstream bed_file;
        std::vector<uintptr_t> genotype_byte(unfiltered_sample_ctl * 2, 0);
        std::string prev_file = "";
        std::streampos prev_loc = 0;
        std::ofstream output;
        uint32_t uii = 0;
        uintptr_t ulii = 0;
        uint32_t ujj;
        uint32_t ukk;
        uint32_t sample_idx = 0;
        size_t nmiss = 0;
        size_t num_not_xvar = 0;
        size_t total = 0;
        // do two pass. First pass get the MAF
        uintptr_t* lbptr;
        output.open(std::string(out + ".eff").c_str());
        if (!output.is_open()) {
            throw std::runtime_error(
                std::string("ERROR: Cannot open file: " + out + ".eff"));
        }
        int num_completed = 0;
        double prev_completed = 0;
        double total_snp = static_cast<double>(m_existed_snps.size());
        for (auto&& snp : m_existed_snps) {
            if (prev_file.compare(snp.file) != 0) {
                if (bed_file.is_open()) {
                    bed_file.close();
                }
                bed_file.open(snp.file.c_str(), std::ios::binary);
                if (!bed_file.is_open()) {
                    std::string error_message =
                        "Error: Cannot open bed file: " + snp.file;
                    throw std::runtime_error(error_message);
                }
                prev_file = snp.file;
                prev_loc = 0;
            }

            if (prev_loc != snp.byte_pos
                && !bed_file.seekg(snp.byte_pos, std::ios_base::beg))
            {
                throw std::runtime_error("Error: Cannot read the bed file!");
            }

            // loadbuf_raw is the temporary

            if (load_and_collapse_incl(m_unfiltered_sample_ct, m_sample_ct,
                                       m_sample_include.data(), final_mask,
                                       bed_file, m_tmp_genotype.data(),
                                       genotype_byte.data()))
            {
                throw std::runtime_error("Error: Cannot read the bed file!");
            }
            prev_loc = bed_file.tellg();
            double eff = effect();
            output << eff << std::endl;
            // get_score(score, genotype_byte, eff, standardize);
            lbptr = genotype_byte.data();
            do
            {
                ulii = ~(*lbptr++);
                if (uii + BITCT2 > m_unfiltered_sample_ct) {
                    ulii &=
                        (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
                }
                ujj = 0;
                while (ulii) {
                    if (uii + (ujj / 2) >= m_sample_ct) {
                        break;
                    }
                    ukk = (ulii >> ujj) & 3;
                    sample_idx = uii + (ujj / 2);
                    if (!m_sample_names[sample_idx].x_var) {
                        ++num_not_xvar;
                        switch (ukk)
                        {
                        default: break;
                        case 1: total += 1; break;
                        case 2: nmiss++; break;
                        case 3: total += 2; break;
                        }
                    }
                    ujj += 2;
                }
                uii += BITCT2;
            } while (uii < m_sample_ct);

            if (num_not_xvar - nmiss == 0) {
                throw std::runtime_error("ERROR: Genotype missingness of 1!");
            }
            double maf = (static_cast<double>(total)
                          / (static_cast<double>(num_not_xvar - nmiss)
                             * 2.0)); // MAF does not count missing
            double var = 1.0;
            double mean = 0.0;
            double miss_dose = maf * 2.0;
            if (standardize) {
                mean = maf * 2;
                var = (sqrt(2.0 * maf * (1.0 - maf)));
            }
            // now start calculating the score
            lbptr = genotype_byte.data();
            do
            {
                ulii = ~(*lbptr++);
                if (uii + BITCT2 > m_unfiltered_sample_ct) {
                    ulii &=
                        (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
                }
                ujj = 0;
                while (ulii) {
                    // ujj = CTZLU(ulii) & (BITCT - 2);
                    if (uii + (ujj / 2) >= m_sample_ct) {
                        break;
                    }
                    ukk = (ulii >> ujj) & 3;
                    sample_idx = uii + (ujj / 2);
                    switch (ukk)
                    {
                    default: break;
                    case 1:
                        score[sample_idx] += eff * (ukk - mean) / var;
                        break;
                    case 2:
                        score[sample_idx] += eff * (miss_dose - mean) / var;
                        break;
                    case 3: score[sample_idx] += eff * (2 - mean) / var; break;
                    }
                    ujj += 2;
                }
                uii += BITCT2;
            } while (uii < m_sample_ct);
            if (num_completed / total_snp - prev_completed > 0.01) {
                fprintf(stderr, "\rProcessing %03.2f%%",
                        num_completed / total_snp);
                prev_completed = num_completed / total_snp;
            }
            num_completed++;
        }
        fprintf(stderr, "\n");
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
    uintptr_t m_founder_ct = 0;

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
        if (unfiltered_sample_ct == sample_ct) {
            rawbuf = mainbuf;
        }
        if (!bedfile.read((char*) rawbuf, unfiltered_sample_ct4)) {
            return 1;
        }
        if (unfiltered_sample_ct != sample_ct) {
            copy_quaterarr_nonempty_subset(rawbuf, sample_include,
                                           unfiltered_sample_ct, sample_ct,
                                           mainbuf);
        }
        else
        {
            mainbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
        }
        // mainbuf should contains the information
        return 0;
    }


    void copy_quaterarr_nonempty_subset(
        const uintptr_t* __restrict raw_quaterarr,
        const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_size,
        uint32_t subset_size, uintptr_t* __restrict output_quaterarr)
    {
        // in plink 2.0, we probably want (0-based) bit raw_quaterarr_size
        // of subset_mask to be always allocated and unset.  This removes a
        // few special cases re: iterating past the end of arrays.
        assert(subset_size);
        assert(raw_quaterarr_size >= subset_size);
        uintptr_t cur_output_word = 0;
        uintptr_t* output_quaterarr_last =
            &(output_quaterarr[subset_size / BITCT2]);
        const uint32_t word_write_halfshift_end = subset_size % BITCT2;
        uint32_t word_write_halfshift = 0;
        // if < 2/3-filled, use sparse copy algorithm
        if (subset_size * (3 * ONELU) < raw_quaterarr_size * (2 * ONELU)) {
            uint32_t subset_mask_widx = 0;
            while (1) {
                const uintptr_t cur_include_word =
                    subset_mask[subset_mask_widx];
                if (cur_include_word) {
                    uint32_t wordhalf_idx = 0;
#ifdef __LP64__
                    uint32_t cur_include_halfword = (uint32_t) cur_include_word;
#else
                    uint32_t cur_include_halfword = (uint16_t) cur_include_word;
#endif
                    while (1) {
                        if (cur_include_halfword) {
                            uintptr_t raw_quaterarr_word =
                                raw_quaterarr[subset_mask_widx * 2
                                              + wordhalf_idx];
                            do
                            {
                                uint32_t rqa_idx_lowbits =
                                    __builtin_ctz(cur_include_halfword);
                                cur_output_word |=
                                    ((raw_quaterarr_word
                                      >> (rqa_idx_lowbits * 2))
                                     & 3)
                                    << (word_write_halfshift * 2);
                                if (++word_write_halfshift == BITCT2) {
                                    *output_quaterarr++ = cur_output_word;
                                    word_write_halfshift = 0;
                                    cur_output_word = 0;
                                }
                                cur_include_halfword &=
                                    cur_include_halfword - 1;
                            } while (cur_include_halfword);
                        }
                        if (wordhalf_idx) {
                            break;
                        }
                        wordhalf_idx++;
#ifdef __LP64__
                        cur_include_halfword = cur_include_word >> 32;
#else
                        cur_include_halfword = cur_include_word >> 16;
#endif
                    }
                    if (output_quaterarr == output_quaterarr_last) {
                        if (word_write_halfshift == word_write_halfshift_end) {
                            if (word_write_halfshift_end) {
                                *output_quaterarr_last = cur_output_word;
                            }
                            return;
                        }
                    }
                }
                subset_mask_widx++;
            }
        }
        // blocked copy
        while (1) {
            const uintptr_t cur_include_word = *subset_mask++;
            uint32_t wordhalf_idx = 0;
#ifdef __LP64__
            uintptr_t cur_include_halfword = (uint32_t) cur_include_word;
#else
            uint32_t cur_include_halfword = (uint16_t) cur_include_word;
#endif
            while (1) {
                uintptr_t raw_quaterarr_word = *raw_quaterarr++;
                while (cur_include_halfword) {
                    uint32_t rqa_idx_lowbits = CTZLU(cur_include_halfword);
                    uintptr_t halfword_invshifted =
                        (~cur_include_halfword) >> rqa_idx_lowbits;
                    uintptr_t raw_quaterarr_curblock_unmasked =
                        raw_quaterarr_word >> (rqa_idx_lowbits * 2);
                    uint32_t rqa_block_len = CTZLU(halfword_invshifted);
                    uint32_t block_len_limit = BITCT2 - word_write_halfshift;
                    cur_output_word |= raw_quaterarr_curblock_unmasked
                                       << (2 * word_write_halfshift);
                    if (rqa_block_len < block_len_limit) {
                        word_write_halfshift += rqa_block_len;
                        cur_output_word &=
                            (ONELU << (word_write_halfshift * 2)) - ONELU;
                    }
                    else
                    {
                        // no need to mask, extra bits vanish off the high
                        // end
                        *output_quaterarr++ = cur_output_word;
                        word_write_halfshift = rqa_block_len - block_len_limit;
                        if (word_write_halfshift) {
                            cur_output_word =
                                (raw_quaterarr_curblock_unmasked
                                 >> (2 * block_len_limit))
                                & ((ONELU << (2 * word_write_halfshift))
                                   - ONELU);
                        }
                        else
                        {
                            // avoid potential right-shift-64
                            cur_output_word = 0;
                        }
                    }
                    cur_include_halfword &=
                        (~(ONELU << (rqa_block_len + rqa_idx_lowbits))) + ONELU;
                }
                if (wordhalf_idx) {
                    break;
                }
                wordhalf_idx++;
#ifdef __LP64__
                cur_include_halfword = cur_include_word >> 32;
#else
                cur_include_halfword = cur_include_word >> 16;
#endif
            }
            if (output_quaterarr == output_quaterarr_last) {
                if (word_write_halfshift == word_write_halfshift_end) {
                    if (word_write_halfshift_end) {
                        *output_quaterarr_last = cur_output_word;
                    }
                    return;
                }
            }
        }
    }


    void get_score(std::vector<double>& score,
                   std::vector<uintptr_t>& geno_byte, const double effect,
                   const bool standardize)
    {
        if (score.size() == 0) {
            return;
        }
        uint32_t uii = 0;
        uintptr_t ulii = 0;
        uint32_t ujj;
        uint32_t ukk;
        uint32_t sample_idx = 0;
        size_t nmiss = 0;
        size_t num_not_xvar = 0;
        size_t total = 0;
        // do two pass. First pass get the MAF

        uintptr_t* lbptr = geno_byte.data();
        do
        {
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct) {
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
            }
            ujj = 0;
            while (ulii) {
                if (uii + (ujj / 2) >= m_sample_ct) {
                    break;
                }
                ukk = (ulii >> ujj) & 3;
                sample_idx = uii + (ujj / 2);
                if (!m_sample_names[sample_idx].x_var) {
                    ++num_not_xvar;
                    switch (ukk)
                    {
                    default: break;
                    case 1: total += 1; break;
                    case 2: nmiss++; break;
                    case 3: total += 2; break;
                    }
                }
                ujj += 2;
            }
            uii += BITCT2;
        } while (uii < m_sample_ct);

        if (num_not_xvar - nmiss == 0) {
            throw std::runtime_error("ERROR: Genotype missingness of 1!");
        }
        double maf = (static_cast<double>(total)
                      / (static_cast<double>(num_not_xvar - nmiss)
                         * 2.0)); // MAF does not count missing
        double var = 1.0;
        double mean = 0.0;
        double miss_dose = maf * 2.0;
        if (standardize) {
            mean = maf * 2;
            var = (sqrt(2.0 * maf * (1.0 - maf)));
        }
        // now start calculating the score
        lbptr = geno_byte.data();
        do
        {
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct) {
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
            }
            ujj = 0;
            while (ulii) {
                // ujj = CTZLU(ulii) & (BITCT - 2);
                if (uii + (ujj / 2) >= m_sample_ct) {
                    break;
                }
                ukk = (ulii >> ujj) & 3;
                sample_idx = uii + (ujj / 2);
                switch (ukk)
                {
                default: break;
                case 1: score[sample_idx] += effect * (ukk - mean) / var; break;
                case 2:
                    score[sample_idx] += effect * (miss_dose - mean) / var;
                    break;
                case 3: score[sample_idx] += effect * (2 - mean) / var; break;
                }
                ujj += 2;
            }
            uii += BITCT2;
        } while (uii < m_sample_ct);
    }
};

#endif /* GENOTYPE_H_ */
