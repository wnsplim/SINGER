//
//  Data_reader.hpp
//  SINGER
//

#ifndef Data_reader_hpp
#define Data_reader_hpp

#include <map>
#include <string>
#include <utility>
#include <vector>

class Data_reader {

public:

    std::string input = "";
    std::string chrom_name = "";
    int ploidy = 2;
    double missing_thres = 0.5;
    std::vector<std::pair<double, double>> masked = {};
    std::map<std::string, std::vector<std::pair<double, double>>> sample_masks = {};
    std::vector<double> unassayed_sites = {};
    std::vector<std::vector<double>> derived_sites = {};
    std::vector<std::vector<double>> missing_sites = {};
    std::vector<std::vector<std::pair<double, double>>> masked_intervals = {};
    bool any_missing = false;
    int unphased_masked = 0;
    int multiallelic_skipped = 0;

    std::vector<std::string> sample_names();

    std::vector<double> read_tip_ages(const std::string &filename);

    static std::vector<double> read_numbers(const std::string &filename, const std::string &missing_message, bool *complete = nullptr);

    static std::vector<std::string> read_last_line(const std::string &filename);

    std::string selected_chrom();

    void read_mask(std::string filename, double start);

    void read(double start_pos, double end_pos, bool with_sites);

    bool in_mask(double x);

private:

    int parse_genotype(const std::string &field, int expected_ploidy, int *calls);

    void check_ploidy(const std::string &field, int n, int expected_ploidy, const int *calls, long long pos, int column);

    void drop_missing_sites(double start_pos, double end_pos);

    void naive_read_sites(double start_pos, double end_pos);

    void guide_read_sites(double start, double end);

    void scan_missing(double start_pos, double end_pos);
};

#endif
