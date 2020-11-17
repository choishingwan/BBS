/*
 * genotype.cpp
 *
 *  Created on: 15 May 2018
 *      Author: shingwanchoi
 */

#include "genotype.h"

Genotype::Genotype(const std::string& prefix)
{
    m_genotype_files = set_genotype_files(prefix);
}

Genotype::~Genotype() {}

std::vector<std::string> Genotype::set_genotype_files(const std::string& prefix)
{
    std::vector<std::string> genotype_files;
    if (prefix.find("#") != std::string::npos)
    {
        for (size_t chr = 1; chr <= 22; ++chr)
        {
            std::string name = prefix;
            misc::replace_substring(name, "#", std::to_string(chr));
            genotype_files.push_back(name);
        }
    }
    else
    {
        genotype_files.push_back(prefix);
    }
    return genotype_files;
}

void Genotype::load_samples(const std::unordered_set<std::string>& sample_list,
                            const std::unordered_set<std::string>& related_list)
{
    m_sample_names = gen_sample_vector(sample_list, related_list);
    std::string message = std::to_string(m_unfiltered_sample_ct) + " people ("
                          + std::to_string(m_num_male) + " male(s), "
                          + std::to_string(m_num_female)
                          + " female(s)) observed\n";
    message.append(std::to_string(m_founder_ct) + " founder(s) included\n");
}


std::vector<Sample>
Genotype::gen_sample_vector(const std::unordered_set<std::string>& sample_list,
                            const std::unordered_set<std::string>& related_list)
{
    assert(m_genotype_files.size() > 0);
    std::string fam_name = m_genotype_files.front() + ".fam";
    auto fam_file = misc::load_stream(fam_name);
    // number of unfiltered samples
    m_unfiltered_sample_ct = 0;
    std::string line;
    std::vector<std::string> token;
    while (std::getline(*fam_file, line))
    {
        misc::trim(line);
        if (!line.empty())
        {
            token = misc::split(line);
            if (token.size() < 6)
            {
                throw std::runtime_error(
                    "Error: Malformed fam file. Less than 6 column on "
                    "line: "
                    + std::to_string(m_unfiltered_sample_ct + 1) + "\n");
            }
            ++m_unfiltered_sample_ct;
        }
    }
    // now reset the fam file to the start
    fam_file->clear();
    fam_file->seekg(0, fam_file->beg);
    // the unfiltered_sampel_ct is used to define the size of all vector used
    // within the program
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    m_founder_info.resize(unfiltered_sample_ctl, 0);
    m_sample_include.resize(unfiltered_sample_ctl, 0);
    m_num_male = 0;
    m_num_female = 0;
    m_num_ambig_sex = 0;
    m_num_non_founder = 0;
    std::vector<Sample> sample_name;
    std::unordered_set<std::string> duplicated_samples;
    uintptr_t sample_index = 0; // this is just used for error message
    bool inclusion = false;
    const bool no_exclusion = sample_list.empty();
    while (std::getline(*fam_file, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() < 6)
        {
            throw std::runtime_error(
                "Error: Malformed fam file. Less than 6 column on line: "
                + std::to_string(sample_index + 1));
        }
        Sample cur_sample;
        cur_sample.FID = token[0];
        cur_sample.IID = token[1];
        inclusion = (no_exclusion
                     || sample_list.find(cur_sample.IID) != sample_list.end());
        if (inclusion)
        {
            if (related_list.find(cur_sample.IID) == related_list.end())
            {
                // only store samples not found in the related list
                // so that when we calculate MAF using the founder_info
                // vector, we will always be fine
                SET_BIT(sample_index, m_founder_info.data());
                ++m_num_unrelated;
                ++m_founder_ct;
            }
            else
            {
                cur_sample.x_var = true;
            }
            SET_BIT(sample_index, m_sample_include.data());
            ++m_sample_ct;
        }
        if (token[5].compare("1") == 0) { ++m_num_male; }
        else if (token[5].compare("2") == 0)
        {
            ++m_num_female;
        }
        else
        {
            ++m_num_ambig_sex;
        }
        ++sample_index;
        if (inclusion) { sample_name.push_back(cur_sample); }
    }
    fam_file.reset();
    // this is the temporary storage for reading in genotype
    // initializing here allow us to save time by not doing
    // initialization each time we read in a SNP
    m_tmp_genotype.resize(unfiltered_sample_ctl * 2, 0);
    return sample_name;
}

std::vector<SNP>
Genotype::gen_snp_vector(const std::unordered_set<std::string>& snp_list,
                         const size_t num_selected, const size_t seed)
{
    std::string line, bim_name, bed_name;
    std::vector<SNP> snp_in_bim;
    std::vector<std::string> bim_info;
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    for (auto prefix : m_genotype_files)
    {
        bim_name = prefix + ".bim";
        bed_name = prefix + ".bed";
        auto bim = misc::load_stream(prefix + ".bim");
        size_t num_snp_read = 0;
        size_t num_snp = 0;
        while (std::getline(*bim, line))
        {
            misc::trim(line);
            if (line.empty()) continue;
            ++num_snp;
        }
        check_bed(bed_name, num_snp);
        bim->clear();
        bim->seekg(0, bim->beg);
        std::string prev_chr = "";
        while (std::getline(*bim, line))
        {
            misc::trim(line);
            if (line.empty()) continue;
            bim_info = misc::split(line);
            if (bim_info.size() < 6)
            {
                throw std::runtime_error(
                    "Error: Malformed bim file. Less than 6 column on "
                    "line: "
                    + std::to_string(num_snp_read) + "\n");
            }
            if (snp_list.empty()
                || snp_list.find(bim_info[1]) != snp_list.end())
            {
                std::streampos byte_pos =
                    m_bed_offset
                    + (num_snp_read
                       * (static_cast<uint64_t>(unfiltered_sample_ct4)));
                snp_in_bim.emplace_back(SNP(prefix, bim_info[1], byte_pos));
            }
            ++num_snp_read;
        }
    }
    //   now we know what SNPs we have, we can do the SNP selection
    std::cerr << snp_in_bim.size() << " SNPs in bim file" << std::endl;
    std::mt19937 g(seed);
    size_t num_selected_snp = num_selected;
    if (num_selected_snp >= snp_in_bim.size())
    {
        // do nothing
    }
    else if (num_selected_snp < snp_in_bim.size())
    {
        // shuffle
        size_t index = 0;
        for (auto iter = snp_in_bim.begin();
             (iter != snp_in_bim.end()) && (num_selected_snp > 0);
             ++iter, index++, --num_selected_snp)
        {
            std::uniform_int_distribution<int> dist(index,
                                                    snp_in_bim.size() - 1);
            const size_t random_index = static_cast<size_t>(dist(g));
            if (*iter != snp_in_bim.at(random_index))
            { std::swap(snp_in_bim.at(random_index), *iter); }
        }
    }
    if (num_selected < snp_in_bim.size())
    {
        snp_in_bim.resize(num_selected);
        sort(snp_in_bim.begin(), snp_in_bim.end(),
             [](const SNP& i1, const SNP& i2) {
                 if (i1.file == i2.file)
                 {
                     if (i1.byte_pos == i2.byte_pos)
                         return i1.name < i2.name;
                     else
                         return i1.byte_pos < i2.byte_pos;
                 }
                 else
                     return i1.file < i2.file;
             });
    }
    return snp_in_bim;
}
void Genotype::load_snps(const std::unordered_set<std::string>& snp_list,
                         const size_t num_selected, const size_t seed)
{
    // do multi-pass analysis, read in the bim file and get the number of SNPs
    m_existed_snps = gen_snp_vector(snp_list, num_selected, seed);
    std::cerr << m_existed_snps.size() << " SNPs remaining" << std::endl;
    // get_maf();
}


void Genotype::check_bed(const std::string& bed_name, const size_t num_marker)
{
    uint32_t uii = 0;
    int64_t llxx = 0;
    int64_t llyy = 0;
    int64_t llzz = 0;
    uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    std::ifstream bed(bed_name.c_str(), std::ios::binary);
    if (!bed.is_open())
    {
        std::string error_message = "Cannot read bed file: " + bed_name;
        throw std::runtime_error(error_message);
    }
    bed.seekg(0, bed.end);
    llxx = bed.tellg();
    if (!llxx) { throw std::runtime_error("Error: Empty .bed file."); }
    bed.seekg(0, bed.beg);
    char version_check[3];
    bed.read(version_check, 3);
    uii = bed.gcount();
    llyy = (static_cast<uint64_t>(unfiltered_sample_ct4)) * num_marker;
    llzz = (static_cast<uint64_t>(m_unfiltered_sample_ct))
           * ((num_marker + 3) / 4);
    bool sample_major = false;
    // compare only the first 3 bytes
    if ((uii == 3) && (!memcmp(version_check, "l\x1b\x01", 3))) { llyy += 3; }
    else if ((uii == 3) && (!memcmp(version_check, "l\x1b", 3)))
    {
        // v1.00 sample-major
        sample_major = true;
        llyy = llzz + 3;
        m_bed_offset = 2;
    }
    else if (uii && (*version_check == '\x01'))
    {
        // v0.99 SNP-major
        llyy += 1;
        m_bed_offset = 1;
    }
    else if (uii && (!(*version_check)))
    {
        // v0.99 sample-major
        sample_major = true;
        llyy = llzz + 1;
        m_bed_offset = 2;
    }
    else
    {
        // pre-v0.99, sample-major, no header bytes
        sample_major = true;
        if (llxx != llzz)
        {
            // probably not PLINK-format at all, so give this error instead
            // of "invalid file size"
            throw std::runtime_error(
                "Error: Invalid header bytes in .bed file.");
        }
        llyy = llzz;
        m_bed_offset = 2;
    }
    if (llxx != llyy)
    {
        if ((*version_check == '#')
            || ((uii == 3) && (!memcmp(version_check, "chr", 3))))
        {
            throw std::runtime_error("Error: Invalid header bytes in PLINK "
                                     "1 .bed file.  (Is this a UCSC "
                                     "Genome\nBrowser BED file instead?)");
        }
        else
        {
            throw std::runtime_error("Error: Invalid .bed file size.");
        }
    }
    if (sample_major)
    {
        throw std::runtime_error(
            "Error: Currently do not support sample major format");
    }
    bed.close();
}

void Genotype::order_effects(const std::string& fix_file,
                             std::vector<double>& effect_sizes)
{
    auto fix = misc::load_stream(fix_file);
    std::unordered_map<std::string, double> snp_effects;
    std::string line;
    std::vector<std::string> token;
    while (std::getline(*fix, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() != 2)
        { std::runtime_error("Error: Expect 2 columns! Input: " + line); }
        auto num = misc::convert<double>(token[1]);
        if (snp_effects.find(token[0]) != snp_effects.end())
        {
            std::cerr << "Warning: Duplicated SNP " << token[0]
                      << " in effect size file" << std::endl;
            std::cerr << "will use the later effect size" << std::endl;
        }
        snp_effects[token[0]] = num;
    }
    fix.reset();
    std::size_t not_found = 0;
    for (auto&& snp : m_existed_snps)
    {
        auto found = snp_effects.find(snp.get_name());
        if (found != snp_effects.end())
        { effect_sizes.push_back(found->second); }
        else
        {
            ++not_found;
        }
    }
}
