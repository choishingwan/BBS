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

Genotype::~Genotype()
{
    // TODO Auto-generated destructor stub
}

std::vector<std::string> Genotype::set_genotype_files(const std::string& prefix)
{
    std::vector<std::string> genotype_files;
    if (prefix.find("#") != std::string::npos) {
        for (size_t chr = 1; chr <= 22; ++chr) {
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

void Genotype::load_samples(const std::unordered_set<std::string>& sample_list)
{
    m_sample_names = gen_sample_vector(sample_list);
    std::string message = std::to_string(m_unfiltered_sample_ct) + " people ("
                          + std::to_string(m_num_male) + " male(s), "
                          + std::to_string(m_num_female)
                          + " female(s)) observed\n";
    message.append(std::to_string(m_founder_ct) + " founder(s) included\n");
}


std::vector<Sample>
Genotype::gen_sample_vector(const std::unordered_set<std::string>& sample_list)
{
    assert(m_genotype_files.size() > 0);
    std::string fam_name = m_genotype_files.front() + ".fam";
    std::ifstream famfile;
    famfile.open(fam_name.c_str());
    bool no_exclusion = sample_list.empty();
    if (!famfile.is_open()) {
        std::string error_message = "Error: Cannot open fam file: " + fam_name;
        throw std::runtime_error(error_message);
    }
    // number of unfiltered samples
    m_unfiltered_sample_ct = 0;
    std::string line;
    while (std::getline(famfile, line)) {
        misc::trim(line);
        if (!line.empty()) {
            std::vector<std::string> token = misc::split(line);
            if (token.size() < 6) {
                std::string message =
                    "Error: Malformed fam file. Less than 6 column on "
                    "line: "
                    + std::to_string(m_unfiltered_sample_ct + 1) + "\n";
                throw std::runtime_error(message);
            }
            m_unfiltered_sample_ct++;
        }
    }
    // now reset the fam file to the start
    famfile.clear();
    famfile.seekg(0);
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

    uintptr_t sample_index = 0; // this is just for error message
    bool inclusion = false;
    while (std::getline(famfile, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if (token.size() < 6) {
            std::string error_message =
                "Error: Malformed fam file. Less than 6 column on line: "
                + std::to_string(sample_index + 1);
            throw std::runtime_error(error_message);
        }
        Sample cur_sample;
        cur_sample.FID = token[0];
        cur_sample.IID = token[1];
        inclusion = (no_exclusion
                     || sample_list.find(cur_sample.IID) != sample_list.end());
        if (inclusion) {
            SET_BIT(sample_index, m_founder_info.data());
            SET_BIT(sample_index, m_sample_include.data());
        }
        m_sample_ct += inclusion;
        if (token[5].compare("1") == 0) {
            m_num_male++;
        }
        else if (token[5].compare("2") == 0)
        {
            m_num_female++;
        }
        else
        {
            m_num_ambig_sex++;
        }
        sample_index++;
        if (inclusion) {
            sample_name.push_back(cur_sample);
        }
    }
    famfile.close();
    m_tmp_genotype.resize(unfiltered_sample_ctl * 2, 0);
    return sample_name;
}

void Genotype::load_snps(const std::unordered_set<std::string>& snp_list,
                         const size_t num_selected, const size_t seed)
{
    m_existed_snps = gen_snp_vector(snp_list);
    std::cerr << "Read in " << m_existed_snps.size() << " SNPs" << std::endl;
    // now randomly select SNPs
    std::mt19937 g(seed);
    size_t num_snp = num_selected;
    size_t index = 0;
    for (auto iter = m_existed_snps.begin();
         (iter != m_existed_snps.end()) && (num_snp > 0);
         ++iter, index++, --num_snp)
    {
        std::uniform_int_distribution<int> dist(index,
                                                m_existed_snps.size() - 1);
        const size_t random_index = static_cast<size_t>(dist(g));
        if (*iter != m_existed_snps.at(random_index)) {
            std::swap(m_existed_snps.at(random_index), *iter);
        }
    }
    m_existed_snps.resize(num_selected);
    std::sort(m_existed_snps.begin(), m_existed_snps.end(),
              [](const SNP i1, const SNP i2) {
                  if (i1.file.compare(i2.file) == 0) {
                      return i1.byte_pos < i2.byte_pos;
                  }
                  else
                      return i1.file.compare(i2.file) < 0;
              });

    std::cerr << m_existed_snps.size() << " SNPs remaining" << std::endl;
}

std::vector<SNP>
Genotype::gen_snp_vector(const std::unordered_set<std::string>& snp_list)
{
    std::vector<SNP> snp_info;
    std::string line;
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    for (auto prefix : m_genotype_files) {
        std::string bim_name = prefix + ".bim";
        std::string bed_name = prefix + ".bed";
        std::ifstream bim(bim_name.c_str());
        if (!bim.is_open()) {
            std::string error_message =
                "Error: Cannot open bim file: " + bim_name;
            throw std::runtime_error(error_message);
        }
        size_t num_snp_read = 0;
        std::string prev_chr = "";
        while (std::getline(bim, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            std::vector<std::string> bim_info = misc::split(line);
            if (bim_info.size() < 6) {
                std::string error_message =
                    "Error: Malformed bim file. Less than 6 column on "
                    "line: "
                    + std::to_string(num_snp_read) + "\n";
                throw std::runtime_error(error_message);
            }
            num_snp_read++;
        }
        bim.clear();
        bim.seekg(0, bim.beg);
        // check if the bed file is valid
        check_bed(bed_name, num_snp_read);

        std::ifstream bed(bed_name.c_str());
        if (!bed.is_open()) {
            std::string error_message =
                "Error: Cannot open bed file: " + bed_name;
            throw std::runtime_error(error_message);
        }
        bed.seekg(m_bed_offset, std::ios_base::beg);
        // now go through the bim & bed file and perform filtering
        num_snp_read = 0;
        int prev_snp_processed = 0;

        bool no_extraction = snp_list.empty();
        while (std::getline(bim, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            num_snp_read++;
            std::vector<std::string> bim_info = misc::split(line);
            if (!no_extraction && snp_list.find(bim_info[1]) == snp_list.end())
                continue;
            if (num_snp_read - prev_snp_processed > 1) {
                // skip unread lines
                if (!bed.seekg(m_bed_offset
                                   + ((num_snp_read - 1)
                                      * ((uint64_t) unfiltered_sample_ct4)),
                               std::ios_base::beg))
                {
                    std::string error_message =
                        "Error: Cannot read the bed file(seek): " + bed_name;
                    throw std::runtime_error(error_message);
                }
            }
            prev_snp_processed = (num_snp_read - 1);
            // get the location of the SNP in the binary file
            // this is used in clumping and PRS calculation which
            // allow us to jump directly to the SNP of interest
            std::streampos byte_pos = bed.tellg();
            snp_info.push_back(SNP(std::string(prefix + ".bed"), byte_pos));
        }
    }
    return snp_info;
}

void Genotype::check_bed(const std::string& bed_name, const size_t num_marker)
{
    uint32_t uii = 0;
    int64_t llxx = 0;
    int64_t llyy = 0;
    int64_t llzz = 0;
    uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    std::ifstream bed(bed_name.c_str(), std::ios::binary);
    if (!bed.is_open()) {
        std::string error_message = "Cannot read bed file: " + bed_name;
        throw std::runtime_error(error_message);
    }
    bed.seekg(0, bed.end);
    llxx = bed.tellg();
    if (!llxx) {
        throw std::runtime_error("Error: Empty .bed file.");
    }
    bed.seekg(0, bed.beg);
    char version_check[3];
    bed.read(version_check, 3);
    uii = bed.gcount();
    llyy = (static_cast<uint64_t>(unfiltered_sample_ct4)) * num_marker;
    llzz = (static_cast<uint64_t>(m_unfiltered_sample_ct))
           * ((num_marker + 3) / 4);
    bool sample_major = false;
    // compare only the first 3 bytes
    if ((uii == 3) && (!memcmp(version_check, "l\x1b\x01", 3))) {
        llyy += 3;
    }
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
        if (llxx != llzz) {
            // probably not PLINK-format at all, so give this error instead
            // of "invalid file size"
            throw std::runtime_error(
                "Error: Invalid header bytes in .bed file.");
        }
        llyy = llzz;
        m_bed_offset = 2;
    }
    if (llxx != llyy) {
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
    if (sample_major) {
        throw std::runtime_error(
            "Error: Currently do not support sample major format");
    }
    bed.close();
}


void Genotype::get_xbeta(std::vector<double>& score, double fixed_effect,
                         bool standardize)
{
    const uintptr_t final_mask = get_final_mask(m_sample_ct);
    // for array size
    const uintptr_t unfiltered_sample_ctl =
        BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    const uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
    std::ifstream bed_file;
    std::vector<uintptr_t> genotype_byte(unfiltered_sample_ctl * 2, 0);
    std::string prev_file = "";
    std::streampos prev_loc = 0;
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
        prev_loc = snp.byte_pos + (std::streampos) unfiltered_sample_ct4;
        // loadbuf_raw is the temporary

        if (load_and_collapse_incl(m_unfiltered_sample_ct, m_sample_ct,
                                   m_sample_include.data(), final_mask,
                                   bed_file, m_tmp_genotype.data(),
                                   genotype_byte.data()))
        {
            throw std::runtime_error("Error: Cannot read the bed file!");
        }
        get_score(score, genotype_byte, fixed_effect, standardize);
    }
}
