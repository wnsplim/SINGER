//
//  main.cpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#include <iostream>
#include "Sampler.hpp"

static bool parse_rate(const string &value, double &rate) {
    try {
        size_t end = 0;
        double x = stod(value, &end);
        if (end == value.size()) {
            rate = x;
            return true;
        }
    } catch (const exception &) {}
    return false;
}

int main(int argc, const char * argv[]) {
    bool resume = false;
    bool debug = false;
    bool no_data = false;
    string mode = "exact";
    string tip_ages_file = "";
    string mask_file = "";
    string chrom_name = "";
    double g = -1;
    double r = -1, m = -1, Ne = -1;
    int num_iters = 0;
    int spacing = 1;
    double start_pos = -1, end_pos = -1;
    string input_filename = "", output_prefix = "";
    string r_value = "", m_value = "";
    double penalty = 0.01;
    double polar = 0.5;
    int scaling_rep = -1;
    int scaling_bin = 100;
    double epsilon_hmm = -1;
    double epsilon_psmc = 0.05;
    int ploidy = 2;
    int seed = 42;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-resume") {
            if (i + 1 < argc && argv[i+1][0] != '-') {
                cerr << "Error: -resume flag doesn't take any value. " << endl;
                exit(1);
            }
            resume = true;
        }
        else if (arg == "-no_data") {
            if (i + 1 < argc && argv[i+1][0] != '-') {
                cerr << "Error: -no_data flag doesn't take any value. " << endl;
                exit(1);
            }
            no_data = true;
        }
        else if (arg == "-tip_ages") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -tip_ages flag cannot be empty. " << endl;
                exit(1);
            }
            tip_ages_file = argv[++i];
        }
        else if (arg == "-mask") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -mask flag cannot be empty. " << endl;
                exit(1);
            }
            mask_file = argv[++i];
        }
        else if (arg == "-chrom") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -chrom flag cannot be empty. " << endl;
                exit(1);
            }
            chrom_name = argv[++i];
        }
        else if (arg == "-g") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -g flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                g = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -g flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "--mode") {
            if (i + 1 >= argc || (string(argv[i+1]) != "original" && string(argv[i+1]) != "exact")) {
                cerr << "Error: --mode flag expects original or exact. " << endl;
                exit(1);
            }
            mode = argv[++i];
        }
        else if (arg == "-debug") {
            if (i + 1 < argc && argv[i+1][0] != '-') {
                cerr << "Error: -debug flag doesn't take any value. " << endl;
                exit(1);
            }
            debug = true;
        }
        else if (arg == "-Ne") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -Ne flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                Ne = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -Ne flag expects a number. " << endl;
                exit(1);
            }
            Ne = 2*Ne;
        }
        else if (arg == "-r") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -r flag cannot be empty. " << endl;
                exit(1);
            }
            r_value = argv[++i];
        }
        else if (arg == "-m") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -m flag cannot be empty. " << endl;
                exit(1);
            }
            m_value = argv[++i];
        }
        else if (arg == "-penalty") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -penalty flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                penalty = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -penalty flag expects a number. " << endl;
                exit(1);
            }
            if (penalty <= 0) {
                cerr << "Error: -penalty must be positive. " << endl;
                exit(1);
            }
        }
        else if (arg == "-polar") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -polar flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                polar = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -polar flag expects a number. " << endl;
                exit(1);
            }
            if (polar <= 0 or polar >= 1) {
                cerr << "Error: -polar flag expects a number between 0 and 1. " << endl;
                exit(1);
            }
        }
        else if (arg == "-scaling_rep") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -scaling_rep flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                scaling_rep = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -scaling_rep flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-scaling_bin") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -scaling_bin flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                scaling_bin = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -scaling_bin flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-hmm_epsilon") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -hmm_epsilon flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                epsilon_hmm = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -hmm_epsilon flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-psmc_bins") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -psmc_epsilon flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                epsilon_psmc = 1.0/stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -p flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-start") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -start flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                start_pos = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -start flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-end") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -end flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                end_pos = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -end flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-input") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -input flag cannot be empty. " << endl;
                exit(1);
            }
            input_filename = argv[++i];
        }
        else if (arg == "-output") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -output flag cannot be empty. " << endl;
                exit(1);
            }
            output_prefix = argv[++i];
        }
        else if (arg == "-n") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -n flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                num_iters = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -n flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-thin") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -thin flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                spacing = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -thin flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-seed") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -seed flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                seed = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -seed flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-ploidy") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -ploidy flag expects a value (1 or 2)." << endl;
                exit(1);
            }

            try {
                ploidy = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -ploidy expects 1 or 2." << endl;
                exit(1);
            }

            if (ploidy != 1 && ploidy != 2) {
                cerr << "Error: -ploidy must be 1 or 2." << endl;
                exit(1);
            }
        }
        else {
            cerr << "Error: Unknown flag. " << arg << endl;
            exit(1);
        }
    }
    if (m_value.size() == 0) {
        cerr << "-m flag missing. " << endl;
        exit(1);
    }
    bool m_is_rate = parse_rate(m_value, m);
    if (r_value.size() == 0) {
        if (!m_is_rate) {
            cerr << "Error: -r is required when -m is a rate map file. " << endl;
            exit(1);
        }
        r_value = m_value;
    }
    bool r_is_rate = parse_rate(r_value, r);
    if (Ne < 0) {
        cerr << "-Ne flag missing or invalid value. " << endl;
        exit(1);
    }
    if (input_filename.size() == 0) {
        cerr << "-input flag missing or invalid value. " << endl;
        exit(1);
    }
    if (output_prefix.size() == 0) {
        cerr << "-output flag missing or invalid value. " << endl;
        exit(1);
    }
    if (num_iters < 0) {
        cerr << "-num_iters flag is invalid. " << endl;
        exit(1);
    }
    if (spacing < 1) {
        cerr << "-thin flag is invalid. " << endl;
        exit(1);
    }
    bool exact = mode == "exact";
    if (tip_ages_file.size() > 0 and g <= 0) {
        cerr << "Error: -g (generation time in years) must be provided with -tip_ages. " << endl;
        exit(1);
    }
    if (scaling_rep < 0) {
        scaling_rep = exact ? 1 : 5;
    }
    if (epsilon_hmm < 0) {
        epsilon_hmm = exact ? 0.1 : 0.001;
    }
    Sampler sampler = Sampler(Ne, r, m);
    if (!r_is_rate) {
        sampler.recomb_map.load_map(r_value, start_pos, end_pos);
        sampler.recomb_rate = sampler.recomb_map.mean_rate(start_pos, end_pos)*Ne;
    }
    if (!m_is_rate) {
        sampler.mut_map.load_map(m_value, start_pos, end_pos);
        sampler.mut_rate = sampler.mut_map.mean_rate(start_pos, end_pos)*Ne;
    }
    if (sampler.mut_rate == 0) {
        cerr << "Error: the mutation rate is set to 0. " << endl;
        exit(1);
    }
    if (sampler.recomb_rate == 0) {
        cerr << "Warning: the recombination rate is set to 0; are you sure? " << endl;
    }
    sampler.penalty = penalty;
    sampler.polar = polar;
    sampler.scaling_rep = scaling_rep;
    Threader_smc::no_data = no_data;
    sampler.scaling_bin = scaling_bin;
    sampler.set_precision(epsilon_hmm, epsilon_psmc);
    sampler.set_input_file_prefix(input_filename);
    sampler.set_output_file_prefix(output_prefix);
    sampler.ploidy = ploidy;
    if (tip_ages_file.size() > 0) {
        sampler.read_tip_ages(tip_ages_file, g);
    }
    sampler.exact = exact;
    sampler.random_seed = seed;
    sampler.start = start_pos;
    sampler.end = end_pos;
    sampler.chrom_name = chrom_name;
    if (mask_file.size() > 0) {
        sampler.read_mask(mask_file);
    }
    if (resume) {
        sampler.sequence_length = end_pos - start_pos;
        sampler.resume_internal_sample(num_iters, spacing);
        return 0;
    } else if (debug) {
        sampler.sequence_length = end_pos - start_pos;
        sampler.debug_resume_internal_sample(num_iters, spacing);
        return 0;
    }
    if (ploidy == 1) {
        sampler.naive_read_vcf_haploid(input_filename, start_pos, end_pos);
    } else {
        sampler.load_vcf(input_filename, start_pos, end_pos);
    }
    sampler.iterative_start();
    sampler.internal_sample(num_iters, spacing);
    return 0;
}
