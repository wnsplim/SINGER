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
#include <set>

struct Genotype_entry {
    double pos;
    int haplotypes[2];
    double L[4];
};

class Data_reader {

public:

    std::string input = "";
    std::string chrom_name = "";
    int ploidy = 2;
    double missing_thres = 0.5;
    bool genotype_lik = false;
    std::string calls = "";
    std::vector<std::pair<double, double>> masked = {};
    std::map<std::string, std::vector<std::pair<double, double>>> sample_masks = {};
    std::vector<double> unassayed_sites = {};
    std::vector<std::vector<double>> derived_sites = {};
    std::vector<std::vector<double>> missing_sites = {};
    std::vector<std::vector<std::pair<double, double>>> masked_intervals = {};
    std::vector<Genotype_entry> entries = {};
    std::vector<double> entry_sites = {};
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

    struct Call_record {
        std::string ref;
        std::vector<std::string> alts;
        bool gl;
        std::vector<std::string> values;
    };

    struct Record_likelihood {
        int field = -1;
        bool gl = false;
        bool swapped = false;
        bool call_gl = false;
        std::vector<std::string> call_values = {};
    };

    struct Genotype_read {
        int n = 0;
        int calls[2] = {-1, -1};
        bool likelihood = false;
        bool from_calls = false;
        bool entry = false;
        bool unphased_missing = false;
        int alleles = 0;
        double L[4] = {0, 0, 0, 0};
    };

    Genotype_read read_genotype(const std::string &field, int individual, double x, long long pos, const Record_likelihood *r);

    std::map<long long, std::vector<Call_record>> call_records = {};
    std::vector<int> call_column = {};
    std::set<long long> input_positions = {};
    std::vector<const std::vector<std::pair<double, double>> *> mask_of_sample = {};
    double dropped_bases_total = 0;
    int dropped_records = 0;
    double deleted_bases = 0;
    int likelihood_genotypes = 0;
    int call_likelihoods = 0;
    int swapped_records = 0;
    int calls_absent = 0;
    bool likelihood_seen = false;
    int sites_without_likelihood = 0;
    int call_sites = 0;

    int parse_genotype(const std::string &field, int expected_ploidy, int *calls, bool &unphased);

    void read_calls(double start_pos, double end_pos);

    Record_likelihood record_likelihood(long long pos, const std::string &ref, const std::string &alt, const std::string &format);

    bool sample_likelihood(const Record_likelihood &r, int sample, const std::string &field, double *G, bool &from_calls);

    bool state_likelihoods(int n, int *calls, bool unphased, const double *likelihoods, double *L, int &alleles);

    bool variant_site(int alt_count, int haplotypes, int alleles);

    bool sample_masked(int sample, double x);

    void add_call_sites(double start_pos, double end_pos, bool with_sites);

    void print_missing_summary();

    void check_ploidy(const std::string &field, int n, int expected_ploidy, const int *calls, long long pos, int column);

    void drop_missing_sites(double start_pos, double end_pos);

    void naive_read_sites(double start_pos, double end_pos);

    void guide_read_sites(double start, double end);

    void scan_missing(double start_pos, double end_pos);
};

#endif
