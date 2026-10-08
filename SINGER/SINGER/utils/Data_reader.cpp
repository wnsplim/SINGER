//
//  Data_reader.cpp
//  SINGER
//

#include "Data_reader.hpp"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <fstream>
#include <limits>
#include <numeric>
#include <sstream>

using namespace std;

static int iupac_mask(const string &allele) {
    const string seq_nt16_str = "=ACMGRSVTWYHKDBN";
    size_t mask = allele.size() == 1 ? seq_nt16_str.find(toupper(allele[0])) : 0;
    return mask == string::npos ? 0 : (int) mask;
}

static void merge_intervals(vector<pair<double, double>> &v) {
    sort(v.begin(), v.end());
    vector<pair<double, double>> merged;
    for (auto &x : v) {
        if (!merged.empty() and x.first <= merged.back().second) {
            merged.back().second = max(merged.back().second, x.second);
        } else {
            merged.push_back(x);
        }
    }
    v = merged;
}

static bool is_unambiguous(const string &allele) {
    int mask = iupac_mask(allele);
    return mask == 1 or mask == 2 or mask == 4 or mask == 8;
}

static bool padding_first(const string &ref, const string &alt) {
    istringstream as(alt);
    string a;
    while (getline(as, a, ',')) {
        if (a.empty() or a[0] != ref[0]) {
            return false;
        }
    }
    return true;
}

static bool padding_last(const string &ref, const string &alt) {
    istringstream as(alt);
    string a;
    while (getline(as, a, ',')) {
        if (a.empty() or a.back() != ref.back()) {
            return false;
        }
    }
    return true;
}

static bool next_at(Vcf_reader &file, const string &chrom, long long pos) {
    string next_line;
    if (!file.peek(next_line)) {
        return false;
    }
    istringstream next_iss(next_line);
    string next_chrom;
    long long next_pos;
    next_iss >> next_chrom >> next_pos;
    return next_pos == pos and next_chrom == chrom;
}

static string first_chrom(const string &prefix) {
    Vcf_reader file(prefix);
    string line;
    while (file.next(line)) {
        if (line[0] != '#') {
            return line.substr(0, line.find('\t'));
        }
    }
    return "";
}

vector<string> Data_reader::sample_names() {
    Vcf_reader file(input);
    string line;
    vector<string> names;
    while (file.next(line)) {
        if (line.substr(0, 6) == "#CHROM") {
            istringstream iss(line);
            string field;
            int i = 0;
            while (iss >> field) {
                if (i >= 9) {
                    names.push_back(field);
                }
                i += 1;
            }
            return names;
        }
    }
    return names;
}

vector<double> Data_reader::read_tip_ages(const string &filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "tip ages file not found: " << filename << endl;
        exit(1);
    }
    vector<string> names = sample_names();
    int n = ploidy*(int) names.size();
    vector<pair<string, double>> rows;
    int ncol = 0;
    string line;
    while (getline(fin, line)) {
        istringstream iss(line);
        vector<string> parts;
        string tok;
        while (iss >> tok) {
            parts.push_back(tok);
        }
        if (parts.empty() or parts[0][0] == '#') {
            continue;
        }
        if (ncol == 0) {
            ncol = (int) parts.size();
        }
        if ((int) parts.size() != ncol or ncol > 2) {
            cerr << "tip ages file: expected 1 column (ages in VCF sample order) or 2 columns (name age)" << endl;
            exit(1);
        }
        rows.push_back({ncol == 2 ? parts[0] : "", stod(parts.back())});
    }
    if (rows.size() != names.size()) {
        cerr << "tip ages file has " << rows.size() << " rows for " << names.size() << " samples in the VCF" << endl;
        exit(1);
    }
    vector<double> ages(n, 0.0);
    if (ncol == 1) {
        for (int i = 0; i < n; i++) {
            ages[i] = rows[i/ploidy].second;
        }
    } else {
        map<string, double> by_name;
        for (auto &r : rows) {
            if (by_name.count(r.first) > 0) {
                cerr << "tip ages file: duplicate name " << r.first << endl;
                exit(1);
            }
            by_name[r.first] = r.second;
        }
        for (int i = 0; i < n; i++) {
            auto it = by_name.find(names[i/ploidy]);
            if (it == by_name.end()) {
                cerr << "tip ages file has no row for " << names[i/ploidy] << endl;
                exit(1);
            }
            ages[i] = it->second;
        }
    }
    return ages;
}

vector<double> Data_reader::read_numbers(const string &filename, const string &missing_message, bool *complete) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << missing_message << endl;
        exit(1);
    }
    vector<double> v;
    double x;
    while (fin >> x) {
        v.push_back(x);
    }
    if (complete != nullptr) {
        *complete = fin.eof();
    }
    return v;
}

vector<string> Data_reader::read_last_line(const string &filename) {
    ifstream file(filename, ios::in);
    vector<string> words;
    if (!file.is_open()) {
        cerr << "Error opening the file: " << filename << endl;
        exit(1);
    }
    file.seekg(-1, ios_base::end);
    bool found = false;
    int newlines_encountered = 0;
    char ch;
    while (file.get(ch)) {
        if (ch == '\n') {
            newlines_encountered++;
            if (newlines_encountered == 2) {
                found = true;
                break;
            }
        }
        file.seekg(-2, ios_base::cur);
    }
    if (!found) {
        return words;
    }
    string last_line_with_content;
    getline(file, last_line_with_content);
    stringstream ss(last_line_with_content);
    string word;
    while (ss >> word) {
        words.push_back(word);
    }
    return words;
}

string Data_reader::selected_chrom() {
    if (chrom_name.empty()) {
        chrom_name = first_chrom(input);
    }
    return chrom_name;
}

int Data_reader::parse_genotype(const string &field, int expected_ploidy, int *calls) {
    size_t stop = field.find(':');
    if (stop == string::npos) {
        stop = field.size();
    }
    int n = 0;
    size_t i = 0;
    char separator = '|';
    while (i < stop and n < expected_ploidy) {
        size_t j = i;
        while (j < stop and field[j] != '|' and field[j] != '/') {
            j++;
        }
        if (n == 0 and j < stop) {
            separator = field[j];
        }
        if (j == i + 1 and field[i] == '0') {
            calls[n] = 0;
        } else if (j == i + 1 and field[i] == '1') {
            calls[n] = 1;
        } else {
            calls[n] = -1;
        }
        n++;
        i = j + 1;
    }
    if (n == 2 and separator == '/' and calls[0] != calls[1]) {
        calls[0] = calls[1] = -1;
        unphased_masked++;
    }
    return i < stop ? -1 : n;
}

void Data_reader::check_ploidy(const string &field, int n, int expected_ploidy, const int *calls, long long pos, int column) {
    if (n >= 0 and (n == expected_ploidy or calls[0] < 0)) {
        return;
    }
    vector<string> names = sample_names();
    cerr << "Error: genotype " << field.substr(0, field.find(':')) << " of sample " << (column < (int) names.size() ? names[column] : to_string(column))
         << " at position " << pos << " has " << (n < 0 ? "more" : "fewer") << " alleles than -ploidy " << expected_ploidy << ". " << endl;
    exit(1);
}

void Data_reader::read_mask(string filename, double start) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "mask file not found: " << filename << endl;
        exit(1);
    }
    string line, chrom, name;
    string vcf_chrom = selected_chrom();
    double lo, hi;
    while (getline(fin, line)) {
        istringstream iss(line);
        if (!(iss >> chrom >> lo >> hi) or chrom != vcf_chrom) {
            continue;
        }
        if (iss >> name) {
            sample_masks[name].push_back({lo + 1 - start, hi + 1 - start});
        } else {
            masked.push_back({lo + 1 - start, hi + 1 - start});
        }
    }
    merge_intervals(masked);
    for (auto &m : sample_masks) {
        merge_intervals(m.second);
    }
}

bool Data_reader::in_mask(double x) {
    auto it = upper_bound(masked.begin(), masked.end(), make_pair(x, numeric_limits<double>::infinity()));
    return it != masked.begin() and prev(it)->second > x;
}

void Data_reader::read(double start_pos, double end_pos, bool with_sites) {
    drop_missing_sites(start_pos, end_pos);
    if (with_sites) {
        ifstream index(input + ".index");
        if (ploidy == 2 and index.is_open() and !Vcf_reader::is_bcf(input)) {
            guide_read_sites(start_pos, end_pos);
        } else {
            naive_read_sites(start_pos, end_pos);
        }
    }
    scan_missing(start_pos, end_pos);
}

void Data_reader::drop_missing_sites(double start_pos, double end_pos) {
    Vcf_reader file(input, selected_chrom());
    string line;
    vector<const vector<pair<double, double>> *> mask_of;
    vector<pair<double, double>> dropped;
    long long prev_pos = -1;
    int calls[2], records = 0;
    double deleted = 0;
    while (file.next(line)) {
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                istringstream hs(line);
                string w;
                for (int k = 0; hs >> w; k++) {
                    if (k >= 9) {
                        auto m = sample_masks.find(w);
                        mask_of.push_back(m == sample_masks.end() ? nullptr : &m->second);
                    }
                }
            }
            continue;
        }
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        long long pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        if (pos < start_pos) {continue;}
        if (pos >= end_pos) {break;}
        double x = pos - start_pos;
        if (ref.size() > 1) {
            bool first = padding_first(ref, alt);
            double lo = first ? x + 1 : x;
            double hi = !first and pos == 1 and padding_last(ref, alt) ? x + ref.size() - 1 : x + ref.size();
            if (hi > lo) {
                dropped.push_back({lo, hi});
                deleted += hi - lo;
            }
        }
        if (pos == prev_pos or in_mask(x) or ref.size() > 1 or alt.size() > 1) {
            prev_pos = pos;
            continue;
        }
        prev_pos = pos;
        bool known[2] = {is_unambiguous(ref), is_unambiguous(alt)};
        int individual = 0, missing = 0;
        while (iss >> genotype) {
            int n = parse_genotype(genotype, ploidy, calls);
            const vector<pair<double, double>> *v = mask_of[individual];
            bool in_sample_mask = false;
            if (v != nullptr) {
                auto it = upper_bound(v->begin(), v->end(), make_pair(x, numeric_limits<double>::infinity()));
                in_sample_mask = it != v->begin() and prev(it)->second > x;
            }
            for (int k = 0; k < ploidy; k++) {
                int c = k < n ? calls[k] : -1;
                missing += in_sample_mask or c < 0 or !known[c];
            }
            individual += 1;
        }
        if (missing > missing_thres*ploidy*individual) {
            dropped.push_back({x, x + 1});
            records += 1;
        }
    }
    double haplotypes = ploidy*(double) mask_of.size(), bases = 0;
    vector<pair<double, int>> events;
    for (auto v : mask_of) {
        if (v != nullptr) {
            for (auto &m : *v) {
                events.push_back({m.first, ploidy});
                events.push_back({m.second, -ploidy});
            }
        }
    }
    sort(events.begin(), events.end());
    int count = 0;
    for (int i = 0; i < (int) events.size(); i++) {
        count += events[i].second;
        if (count > missing_thres*haplotypes and i + 1 < (int) events.size() and events[i + 1].first > events[i].first) {
            dropped.push_back({events[i].first, events[i + 1].first});
            bases += events[i + 1].first - events[i].first;
        }
    }
    masked.insert(masked.end(), dropped.begin(), dropped.end());
    merge_intervals(masked);
    cout << "dropped for missing data above " << missing_thres << ": " << records << " records, " << bases << " masked bases" << endl;
    cout << "masked bases removed by deletions: " << deleted << endl;
}

void Data_reader::naive_read_sites(double start_pos, double end_pos) {
    Vcf_reader file(input, selected_chrom());
    string line;
    int num_individuals = 0;
    long long prev_pos = -1;
    int valid_mutation = 0;
    int removed_mutation = 0;
    vector<double> genotypes = {};
    int calls[2];
    while (file.next(line)) {
        if (line.substr(0, 6) == "#CHROM") {
            istringstream iss(line);
            vector<string> fields;
            string field;
            while (iss >> field) {
                fields.push_back(field);
            }
            num_individuals = (int) fields.size() - 9;
            derived_sites.assign(ploidy*num_individuals, {});
            genotypes.resize(ploidy*num_individuals);
            continue;
        } else if (line[0] == '#') {
            continue;
        }
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        long long pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        if (pos < start_pos) {continue;}
        if (pos >= end_pos) {break;}
        if (pos == prev_pos) {continue;}
        if (in_mask(pos - start_pos)) {
            prev_pos = pos;
            continue;
        }
        if (!is_unambiguous(ref) or !is_unambiguous(alt)) {
            removed_mutation += 1;
            prev_pos = pos;
            continue;
        }
        if (next_at(file, chrom, pos)) {
            removed_mutation += 1;
            prev_pos = pos;
            continue;
        }
        int individual_index = 0;
        while (iss >> genotype) {
            int n = parse_genotype(genotype, ploidy, calls);
            check_ploidy(genotype, n, ploidy, calls, pos, individual_index);
            for (int k = 0; k < ploidy; k++) {
                genotypes[ploidy*individual_index + k] = (k < n and calls[k] == 1) ? 1 : 0;
            }
            individual_index += 1;
        }
        int genotype_sum = accumulate(genotypes.begin(), genotypes.end(), 0.0);
        if (genotype_sum >= 1 and genotype_sum < genotypes.size()) {
            valid_mutation += 1;
            for (int i = 0; i < genotypes.size(); i++) {
                if (genotypes[i] == 1) {
                    derived_sites[i].push_back(pos - start_pos);
                }
            }
        }
    }
    cout << "valid mutations: " << valid_mutation << endl;
    cout << "removed mutations: " << removed_mutation << endl;
}

void Data_reader::guide_read_sites(double start, double end) {
    string index_file = input + ".index";
    ifstream idx_stream(index_file);
    string line;
    long byte_offset = -1;
    while (getline(idx_stream, line)) {
        istringstream iss(line);
        vector<string> w;
        string field;
        while (iss >> field) {
            w.push_back(field);
        }
        if (w.size() < 2) {
            continue;
        }
        if (w.size() >= 3 and w[0] != selected_chrom()) {
            continue;
        }
        if (stod(w[w.size() - 2]) == start) {
            byte_offset = stol(w.back());
            break;
        }
    }
    if (byte_offset == -1) {
        cerr << "Start position not found in index file: " + index_file << endl;
        exit(1);
    }
    Vcf_reader vcf_stream(input, selected_chrom());
    vcf_stream.seek(byte_offset);
    long long prev_pos = -1;
    int valid_mutation = 0;
    int removed_mutation = 0;
    vector<double> genotypes = {};
    int calls[2];
    while (vcf_stream.next(line)) {
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        long long pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        if (pos == prev_pos) {continue;}
        if (pos >= end) {break;}
        if (in_mask(pos - start)) {
            prev_pos = pos;
            continue;
        }
        if (!is_unambiguous(ref) or !is_unambiguous(alt)) {
            removed_mutation += 1;
            prev_pos = pos;
            continue;
        }
        if (next_at(vcf_stream, chrom, pos)) {
            removed_mutation += 1;
            prev_pos = pos;
            continue;
        }
        int individual_index = 0;
        while (iss >> genotype) {
            if (genotypes.size() < 2*individual_index + 2) {
                genotypes.resize(2*individual_index + 2);
            }
            int n = parse_genotype(genotype, 2, calls);
            check_ploidy(genotype, n, 2, calls, pos, individual_index);
            for (int k = 0; k < 2; k++) {
                genotypes[2*individual_index + k] = (k < n and calls[k] == 1) ? 1 : 0;
            }
            individual_index += 1;
        }
        if (derived_sites.size() == 0) {
            derived_sites.assign(genotypes.size(), {});
        } else {
            assert(derived_sites.size() == genotypes.size());
        }
        int genotype_sum = accumulate(genotypes.begin(), genotypes.end(), 0.0);
        if (genotype_sum >= 1 and genotype_sum < genotypes.size()) {
            valid_mutation += 1;
            for (int i = 0; i < genotypes.size(); i++) {
                if (genotypes[i] == 1) {
                    derived_sites[i].push_back(pos - start);
                }
            }
        }
    }
    if (valid_mutation < 3) {
        cerr << "there are too few variants in this region, algorithm not run" << endl;
    }
    cout << "valid mutations: " << valid_mutation << endl;
    cout << "removed mutations: " << removed_mutation << endl;
}

void Data_reader::scan_missing(double start_pos, double end_pos) {
    Vcf_reader file(input, selected_chrom());
    string line;
    long long prev_pos = -1;
    int calls[2];
    unphased_masked = 0;
    multiallelic_skipped = 0;
    vector<string> names;
    while (file.next(line)) {
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                istringstream hs(line);
                string w;
                for (int k = 0; hs >> w; k++) {
                    if (k >= 9) {
                        names.push_back(w);
                    }
                }
                missing_sites.assign(ploidy*names.size(), {});
                masked_intervals.assign(ploidy*names.size(), {});
            }
            continue;
        }
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        long long pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        if (pos < start_pos) {continue;}
        if (pos >= end_pos) {break;}
        if (pos == prev_pos) {continue;}
        if (in_mask(pos - start_pos)) {
            prev_pos = pos;
            continue;
        }
        if (ref.size() > 1 or alt.size() > 1) {
            if (!padding_first(ref, alt) or next_at(file, chrom, pos)) {
                unassayed_sites.push_back(pos - start_pos);
            }
            multiallelic_skipped += alt.find(',') != string::npos;
            prev_pos = pos;
            continue;
        }
        if (next_at(file, chrom, pos)) {
            unassayed_sites.push_back(pos - start_pos);
            prev_pos = pos;
            continue;
        }
        int individual = 0;
        vector<int> missing_haplotypes;
        bool known[2] = {is_unambiguous(ref), is_unambiguous(alt)};
        while (iss >> genotype) {
            int n = parse_genotype(genotype, ploidy, calls);
            check_ploidy(genotype, n, ploidy, calls, pos, individual);
            for (int k = 0; k < ploidy; k++) {
                int c = k < n ? calls[k] : -1;
                if (c < 0 or !known[c]) {
                    missing_haplotypes.push_back(ploidy*individual + k);
                }
            }
            individual += 1;
        }
        if (missing_haplotypes.size() == missing_sites.size()) {
            unassayed_sites.push_back(pos - start_pos);
        } else if (missing_haplotypes.size() > 0) {
            any_missing = true;
            for (int h : missing_haplotypes) {
                missing_sites[h].push_back(pos - start_pos);
            }
        }
    }
    sort(unassayed_sites.begin(), unassayed_sites.end());
    for (auto &m : sample_masks) {
        auto it = find(names.begin(), names.end(), m.first);
        if (it == names.end()) {
            Vcf_reader::fail("the mask file names sample " + m.first + ", which is not in " + input);
        }
        any_missing = true;
        int i = (int) (it - names.begin());
        for (int k = 0; k < ploidy; k++) {
            int h = ploidy*i + k;
            masked_intervals[h] = m.second;
            vector<double> kept;
            for (double x : missing_sites[h]) {
                auto g = upper_bound(masked_intervals[h].begin(), masked_intervals[h].end(), make_pair(x, numeric_limits<double>::infinity()));
                if (g == masked_intervals[h].begin() or prev(g)->second <= x) {
                    kept.push_back(x);
                }
            }
            missing_sites[h] = kept;
        }
    }
    if (unphased_masked > 0) {
        cerr << "Warning: unphased heterozygous genotypes treated as missing: " << unphased_masked << ". " << endl;
    }
    if (multiallelic_skipped > 0) {
        cerr << "Warning: multiallelic sites skipped: " << multiallelic_skipped << ". " << endl;
    }
}
