//
//  Data_reader.cpp
//  SINGER
//

#include "Data_reader.hpp"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <sys/stat.h>
#include <zlib.h>
#include <htslib/vcf.h>

using namespace std;

class Vcf_reader {

public:

    static void fail(const string &message) {
        cerr << "Error: " << message << ". " << endl;
        exit(1);
    }

    static bool is_bcf(const string &p) {
        return p.size() > 4 and p.compare(p.size() - 4, 4, ".bcf") == 0;
    }

    Vcf_reader(const string &input, const string &chrom = "") : path(input), chrom(chrom) {
        struct stat st;
        if (stat(path.c_str(), &st) != 0 or !S_ISREG(st.st_mode)) {
            fail("file " + path + " not found");
        }
        if (is_bcf(path)) {
            hts_set_log_level(HTS_LOG_OFF);
            bcf = hts_open(path.c_str(), "r");
            header = bcf == nullptr ? nullptr : bcf_hdr_read(bcf);
            if (header != nullptr) {
                bcf_hdr_format(header, 0, &text);
                record = bcf_init();
            }
        } else {
            file = gzopen(path.c_str(), "rb");
        }
        if (file == nullptr and header == nullptr) {
            fail("cannot read " + path);
        }
    }

    ~Vcf_reader() {
        if (file != nullptr) {
            gzclose(file);
        }
        if (record != nullptr) {
            bcf_destroy(record);
        }
        if (header != nullptr) {
            bcf_hdr_destroy(header);
        }
        if (bcf != nullptr) {
            hts_close(bcf);
        }
        free(text.s);
        free(buffer.s);
    }

    bool next(string &line) {
        if (has_pending) {
            line = move(pending);
            has_pending = false;
            return true;
        }
        return read(line);
    }

    bool peek(string &line) {
        if (!has_pending) {
            has_pending = read(pending);
        }
        if (has_pending) {
            line = pending;
        }
        return has_pending;
    }

    void seek(long offset) {
        gzseek(file, offset, SEEK_SET);
        has_pending = false;
    }

    string path;

private:

    string chrom;
    bool seen = false;
    bool done = false;
    gzFile file = nullptr;
    htsFile *bcf = nullptr;
    bcf_hdr_t *header = nullptr;
    bcf1_t *record = nullptr;
    kstring_t text = {0, 0, nullptr};
    kstring_t buffer = {0, 0, nullptr};
    size_t at = 0;
    string pending;
    bool has_pending = false;

    bool read(string &line) {
        while (!done and read_line(line)) {
            if (chrom.empty() or line[0] == '#') {
                return true;
            }
            if (line.compare(0, chrom.size() + 1, chrom + '\t') == 0) {
                seen = true;
                return true;
            }
            done = seen;
        }
        if (!chrom.empty() and !seen) {
            fail("no record of chromosome " + chrom + " in " + path);
        }
        return false;
    }

    bool read_line(string &line) {
        line.clear();
        if (header != nullptr) {
            if (at < text.l) {
                const char *s = text.s + at;
                const char *e = (const char *) memchr(s, '\n', text.l - at);
                size_t n = e == nullptr ? text.l - at : e - s;
                line.assign(s, n);
                at += n + 1;
                return true;
            }
            int status = bcf_read(bcf, header, record);
            if (status == -1) {
                return false;
            }
            buffer.l = 0;
            if (status < -1 or vcf_format(header, record, &buffer) < 0) {
                fail("cannot read " + path);
            }
            line.assign(buffer.s, buffer.l);
            if (!line.empty() and line.back() == '\n') {
                line.pop_back();
            }
            return true;
        }
        char chunk[65536];
        while (gzgets(file, chunk, sizeof(chunk)) != nullptr) {
            line += chunk;
            if (!line.empty() and line.back() == '\n') {
                line.pop_back();
                if (line.empty()) {
                    continue;
                }
                return true;
            }
        }
        int status = Z_OK;
        gzerror(file, &status);
        if (status != Z_OK) {
            fail("cannot read " + path);
        }
        return !line.empty();
    }
};

static vector<string> split(const string &s, char sep) {
    vector<string> out;
    string part;
    istringstream ss(s);
    while (getline(ss, part, sep)) {
        out.push_back(part);
    }
    return out;
}

static bool parse_likelihood(const string &s, bool gl, int count, double *G) {
    vector<string> v = split(s, ',');
    if ((int) v.size() != count) {
        return false;
    }
    for (int k = 0; k < count; k++) {
        char *e = nullptr;
        double x = strtod(v[k].c_str(), &e);
        if (e == v[k].c_str() or *e != '\0' or !isfinite(x)) {
            return false;
        }
        G[k] = gl ? pow(10, x) : pow(10, -x/10);
    }
    return *max_element(G, G + count) > *min_element(G, G + count);
}

static string single_alt(const vector<string> &alts) {
    string found = "";
    for (const string &a : alts) {
        if (a == "<*>" or a == "<NON_REF>") {
            continue;
        }
        if (!found.empty()) {
            return "";
        }
        found = a;
    }
    return found;
}

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

int Data_reader::parse_genotype(const string &field, int expected_ploidy, int *calls, bool &unphased) {
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
    unphased = n == 2 and separator == '/' and calls[0] != calls[1];
    return i < stop ? -1 : n;
}

void Data_reader::read_calls(double start_pos, double end_pos) {
    call_records.clear();
    call_column.clear();
    if (!genotype_lik or calls.empty()) {
        return;
    }
    vector<string> names = sample_names();
    Vcf_reader file(calls, selected_chrom());
    string line;
    while (file.next(line)) {
        if (line.rfind("#CHROM", 0) == 0) {
            vector<string> w = split(line, '\t');
            map<string, int> column;
            for (int i = 9; i < (int) w.size(); i++) {
                column[w[i]] = i - 9;
            }
            if (column.empty()) {
                Vcf_reader::fail("-genotype_calls file " + calls + " has no sample columns");
            }
            vector<string> absent;
            for (const string &s : names) {
                auto c = column.find(s);
                call_column.push_back(c == column.end() ? -1 : c->second);
                if (c == column.end()) {
                    absent.push_back(s);
                }
            }
            calls_absent = (int) absent.size();
            if (calls_absent == (int) names.size()) {
                Vcf_reader::fail("no sample of the input is a sample of -genotype_calls " + calls);
            }
            if (calls_absent > 0) {
                cerr << "Warning: " << calls_absent << " input samples have no column in -genotype_calls (";
                for (int i = 0; i < calls_absent and i < 10; i++) {
                    cerr << (i > 0 ? ", " : "") << absent[i];
                }
                cerr << (calls_absent > 10 ? ", ..." : "") << "): their genotypes are read as written. " << endl;
            }
            continue;
        }
        if (line[0] == '#') {
            continue;
        }
        vector<string> w = split(line, '\t');
        if (w.size() < 10) {
            Vcf_reader::fail("-genotype_calls file " + calls + " has a record with fewer than 10 columns");
        }
        long long pos = stoll(w[1]);
        if (pos < start_pos) {
            continue;
        }
        if (pos >= end_pos) {
            break;
        }
        vector<string> keys = split(w[8], ':');
        int gi = (int) (find(keys.begin(), keys.end(), "GL") - keys.begin());
        int pi = (int) (find(keys.begin(), keys.end(), "PL") - keys.begin());
        int fi = gi < (int) keys.size() ? gi : pi;
        if (fi == (int) keys.size()) {
            continue;
        }
        Call_record r = {w[3], split(w[4], ','), gi < (int) keys.size(), {}};
        for (int i = 9; i < (int) w.size(); i++) {
            vector<string> parts = split(w[i], ':');
            r.values.push_back(fi < (int) parts.size() ? parts[fi] : "");
        }
        call_records[pos].push_back(r);
    }
}

Data_reader::Record_likelihood Data_reader::record_likelihood(long long pos, const string &ref, const string &alt, const string &format) {
    Record_likelihood r;
    vector<string> keys = split(format, ':');
    int gi = (int) (find(keys.begin(), keys.end(), "GL") - keys.begin());
    int pi = (int) (find(keys.begin(), keys.end(), "PL") - keys.begin());
    r.gl = gi < (int) keys.size();
    r.field = r.gl ? gi : pi < (int) keys.size() ? pi : -1;
    auto it = call_records.find(pos);
    if (it == call_records.end()) {
        return r;
    }
    for (const Call_record &c : it->second) {
        int j = (int) (find(c.alts.begin(), c.alts.end(), alt) - c.alts.begin()) + 1;
        bool direct = c.ref == ref and j <= (int) c.alts.size();
        bool swapped = !direct and single_alt(c.alts) == ref and c.ref == alt;
        if (!direct and !swapped) {
            continue;
        }
        if (swapped) {
            j = (int) (find(c.alts.begin(), c.alts.end(), ref) - c.alts.begin()) + 1;
        }
        int alleles = (int) c.alts.size() + 1;
        int b = j*(j + 1)/2;
        r.swapped = swapped;
        r.call_gl = c.gl;
        r.call_values.assign(call_column.size(), "");
        for (int i = 0; i < (int) call_column.size(); i++) {
            if (call_column[i] < 0 or call_column[i] >= (int) c.values.size()) {
                continue;
            }
            vector<string> v = split(c.values[call_column[i]], ',');
            if ((int) v.size() == alleles*(alleles + 1)/2) {
                r.call_values[i] = swapped ? v[b + j] + "," + v[b] + "," + v[0] : v[0] + "," + v[b] + "," + v[b + j];
            } else if ((int) v.size() == alleles) {
                r.call_values[i] = swapped ? v[j] + "," + v[0] : v[0] + "," + v[j];
            }
        }
        break;
    }
    return r;
}

bool Data_reader::sample_likelihood(const Record_likelihood &r, int sample, const string &field, double *G, bool &from_calls) {
    int count = ploidy == 1 ? 2 : 3;
    from_calls = false;
    if (r.field >= 0) {
        vector<string> parts = split(field, ':');
        if (r.field < (int) parts.size() and parse_likelihood(parts[r.field], r.gl, count, G)) {
            likelihood_seen = true;
            return true;
        }
    }
    if (sample < (int) r.call_values.size() and !r.call_values[sample].empty() and parse_likelihood(r.call_values[sample], r.call_gl, count, G)) {
        from_calls = true;
        likelihood_seen = true;
        return true;
    }
    return false;
}

bool Data_reader::state_likelihoods(int n, int *calls, bool unphased, const double *likelihoods, double *L, int &alleles) {
    int nc = ploidy == 1 ? 2 : 4;
    int ng = ploidy == 1 ? 2 : 3;
    double G[3] = {0, 0, 0};
    for (int g = 0; g < ng; g++) {
        G[g] = likelihoods[g];
    }
    int start = 0;
    if (ploidy == 1) {
        L[0] = G[0];
        L[1] = G[1];
        L[2] = L[3] = 0;
        start = calls[0] >= 0 ? calls[0] : G[1] > G[0] ? 1 : 0;
    } else if (n == 2 and calls[0] >= 0 and calls[1] >= 0 and !unphased) {
        int a = calls[0], b = calls[1];
        L[0] = G[0];
        L[1] = a == b or a == 1 ? G[1] : 0;
        L[2] = a == b or b == 1 ? G[1] : 0;
        L[3] = G[2];
        start = a + 2*b;
    } else {
        L[0] = G[0];
        L[1] = L[2] = G[1];
        L[3] = G[2];
        int best = G[1] > G[0] ? 1 : 0;
        best = G[2] > G[best] ? 2 : best;
        start = unphased and calls[0] >= 0 and calls[1] >= 0 ? 2 : best == 0 ? 0 : best == 1 ? 2 : 3;
    }
    if (L[start] == 0) {
        start = (int) (max_element(L, L + nc) - L);
    }
    double top = *max_element(L, L + nc);
    bool other = false;
    alleles = 0;
    for (int c = 0; c < nc; c++) {
        L[c] /= top;
        other = other or (c != start and L[c] > 0);
        if (L[c] > 0) {
            alleles |= ploidy == 1 ? 1 << c : (1 << (c & 1)) | (1 << (c >> 1));
        }
    }
    calls[0] = start & 1;
    if (ploidy == 2) {
        calls[1] = start >> 1;
    }
    return other;
}

bool Data_reader::variant_site(int alt_count, int haplotypes, int alleles) {
    if (alt_count >= 1 and alt_count < haplotypes) {
        return true;
    }
    return (alleles >> (alt_count == 0 ? 1 : 0)) & 1;
}

Data_reader::Genotype_read Data_reader::read_genotype(const string &field, int individual, double x, long long pos, const Record_likelihood *r) {
    Genotype_read g;
    bool unphased;
    g.n = parse_genotype(field, ploidy, g.calls, unphased);
    check_ploidy(field, g.n, ploidy, g.calls, pos, individual);
    if (sample_masked(individual, x)) {
        g.calls[0] = g.calls[1] = -1;
        return g;
    }
    double G[3];
    if (r != nullptr and sample_likelihood(*r, individual, field, G, g.from_calls)) {
        g.likelihood = true;
        g.entry = state_likelihoods(g.n, g.calls, unphased, G, g.L, g.alleles);
        g.n = ploidy;
    } else if (unphased) {
        g.calls[0] = g.calls[1] = -1;
        g.unphased_missing = true;
    }
    return g;
}

bool Data_reader::sample_masked(int sample, double x) {
    if (sample >= (int) mask_of_sample.size() or mask_of_sample[sample] == nullptr) {
        return false;
    }
    const vector<pair<double, double>> *v = mask_of_sample[sample];
    auto it = upper_bound(v->begin(), v->end(), make_pair(x, numeric_limits<double>::infinity()));
    return it != v->begin() and prev(it)->second > x;
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
    entries.clear();
    entry_sites.clear();
    likelihood_genotypes = call_likelihoods = swapped_records = calls_absent = sites_without_likelihood = call_sites = 0;
    mask_of_sample.clear();
    for (const string &s : sample_names()) {
        auto m = sample_masks.find(s);
        mask_of_sample.push_back(m == sample_masks.end() ? nullptr : &m->second);
    }
    read_calls(start_pos, end_pos);
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
    add_call_sites(start_pos, end_pos, with_sites);
    print_missing_summary();
    if (genotype_lik and !likelihood_seen) {
        cerr << "Warning: -genotype_lik is set but no usable GL or PL was found in the input"
             << (calls.empty() ? " and no -genotype_calls is given" : " or in -genotype_calls " + calls)
             << ": genotypes are read as written. " << endl;
    }
}

void Data_reader::add_call_sites(double start_pos, double end_pos, bool with_sites) {
    if (call_records.empty()) {
        return;
    }
    int haplotypes = ploidy*(int) call_column.size();
    if (with_sites and (int) derived_sites.size() < haplotypes) {
        derived_sites.resize(haplotypes);
    }
    bool added = false;
    for (auto &p : call_records) {
        long long pos = p.first;
        double x = pos - start_pos;
        if (pos < start_pos or pos >= end_pos or input_positions.count(pos) > 0 or in_mask(x)) {
            continue;
        }
        string ref = p.second.front().ref, alt = single_alt(p.second.front().alts);
        if (ref.size() != 1 or alt.size() != 1 or !is_unambiguous(ref) or !is_unambiguous(alt)) {
            continue;
        }
        Record_likelihood r = record_likelihood(pos, ref, alt, "");
        vector<int> state(haplotypes, -1);
        vector<Genotype_entry> found;
        int missing = 0, from_file = 0, alt_count = 0, alleles = 0;
        for (int i = 0; i < (int) call_column.size(); i++) {
            Genotype_read g = read_genotype("", i, x, pos, &r);
            if (!g.likelihood) {
                missing += ploidy;
                continue;
            }
            from_file += 1;
            alleles |= g.alleles;
            if (g.entry) {
                found.push_back({x, {ploidy*i, ploidy == 2 ? ploidy*i + 1 : -1}, {g.L[0], g.L[1], g.L[2], g.L[3]}});
            }
            for (int k = 0; k < ploidy; k++) {
                state[ploidy*i + k] = g.calls[k];
                alt_count += g.calls[k] == 1;
            }
        }
        if (!variant_site(alt_count, haplotypes, alleles)) {
            continue;
        }
        if (missing > missing_thres*haplotypes) {
            masked.push_back({x, x + 1});
            dropped_records += 1;
            continue;
        }
        added = true;
        call_sites += 1;
        call_likelihoods += from_file;
        likelihood_genotypes += (int) found.size();
        entries.insert(entries.end(), found.begin(), found.end());
        entry_sites.push_back(x);
        for (int h = 0; h < haplotypes; h++) {
            if (state[h] == 1 and with_sites) {
                derived_sites[h].push_back(x);
            } else if (state[h] < 0 and !sample_masked(h/ploidy, x)) {
                missing_sites[h].push_back(x);
                any_missing = true;
            }
        }
    }
    if (!added) {
        return;
    }
    merge_intervals(masked);
    for (auto &v : derived_sites) {
        sort(v.begin(), v.end());
    }
    for (auto &v : missing_sites) {
        sort(v.begin(), v.end());
    }
    stable_sort(entries.begin(), entries.end(), [](const Genotype_entry &a, const Genotype_entry &b) { return a.pos < b.pos; });
    sort(entry_sites.begin(), entry_sites.end());
}

void Data_reader::print_missing_summary() {
    cout << "missing data: " << dropped_records << " records dropped above -missing_thres " << missing_thres << ", "
         << dropped_bases_total << " bases masked by the sample masks above -missing_thres, " << deleted_bases << " bases removed by deletions, "
         << unphased_masked << " unphased heterozygous genotypes treated as missing, " << multiallelic_skipped << " multiallelic sites skipped";
    if (genotype_lik) {
        cout << "; genotype likelihoods: " << likelihood_genotypes << " genotypes with a likelihood, " << call_likelihoods << " taken from -genotype_calls, "
             << swapped_records << " records with REF and ALT swapped, " << calls_absent << " samples absent from -genotype_calls, "
             << sites_without_likelihood << " sites without genotype likelihoods taken as true, "
             << call_sites << " positions added from -genotype_calls";
    }
    cout << endl;
}

void Data_reader::drop_missing_sites(double start_pos, double end_pos) {
    Vcf_reader file(input, selected_chrom());
    string line;
    vector<pair<double, double>> dropped;
    long long prev_pos = -1;
    int records = 0;
    double deleted = 0;
    while (file.next(line)) {
        if (line[0] == '#') {
            continue;
        }
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        long long pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        if (pos < start_pos) {continue;}
        if (pos >= end_pos) {break;}
        input_positions.insert(pos);
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
        bool weigh = genotype_lik and known[0] and known[1];
        Record_likelihood r;
        if (weigh) {
            r = record_likelihood(pos, ref, alt, format);
        }
        int individual = 0, missing = 0;
        while (iss >> genotype) {
            Genotype_read g = read_genotype(genotype, individual, x, pos, weigh ? &r : nullptr);
            for (int k = 0; k < ploidy; k++) {
                int c = k < g.n ? g.calls[k] : -1;
                missing += c < 0 or !known[c];
            }
            individual += 1;
        }
        if (missing > missing_thres*ploidy*individual) {
            dropped.push_back({x, x + 1});
            records += 1;
        }
    }
    double haplotypes = ploidy*(double) mask_of_sample.size(), bases = 0;
    vector<pair<double, int>> events;
    for (auto v : mask_of_sample) {
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
    dropped_records = records;
    dropped_bases_total = bases;
    deleted_bases = deleted;
}

void Data_reader::naive_read_sites(double start_pos, double end_pos) {
    Vcf_reader file(input, selected_chrom());
    string line;
    int num_individuals = 0;
    long long prev_pos = -1;
    int valid_mutation = 0;
    int removed_mutation = 0;
    vector<double> genotypes = {};
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
        Record_likelihood r;
        if (genotype_lik) {
            r = record_likelihood(pos, ref, alt, format);
        }
        int individual_index = 0;
        while (iss >> genotype) {
            Genotype_read g = read_genotype(genotype, individual_index, pos - start_pos, pos, genotype_lik ? &r : nullptr);
            for (int k = 0; k < ploidy; k++) {
                genotypes[ploidy*individual_index + k] = (k < g.n and g.calls[k] == 1) ? 1 : 0;
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
        Record_likelihood r;
        if (genotype_lik) {
            r = record_likelihood(pos, ref, alt, format);
        }
        int individual_index = 0;
        while (iss >> genotype) {
            if (genotypes.size() < 2*individual_index + 2) {
                genotypes.resize(2*individual_index + 2);
            }
            Genotype_read g = read_genotype(genotype, individual_index, pos - start, pos, genotype_lik ? &r : nullptr);
            for (int k = 0; k < 2; k++) {
                genotypes[2*individual_index + k] = (k < g.n and g.calls[k] == 1) ? 1 : 0;
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
        bool use_likelihood = genotype_lik and known[0] and known[1];
        Record_likelihood r;
        if (use_likelihood) {
            r = record_likelihood(pos, ref, alt, format);
            swapped_records += r.swapped;
        }
        vector<Genotype_entry> found;
        double x = pos - start_pos;
        bool any_likelihood = false;
        int alt_count = 0, alleles = 0;
        while (iss >> genotype) {
            Genotype_read g = read_genotype(genotype, individual, x, pos, use_likelihood ? &r : nullptr);
            if (g.entry) {
                found.push_back({x, {ploidy*individual, ploidy == 2 ? ploidy*individual + 1 : -1}, {g.L[0], g.L[1], g.L[2], g.L[3]}});
            }
            call_likelihoods += g.from_calls;
            any_likelihood = any_likelihood or g.likelihood;
            alleles |= g.alleles;
            unphased_masked += g.unphased_missing;
            for (int k = 0; k < ploidy; k++) {
                int c = k < g.n ? g.calls[k] : -1;
                if (c < 0 or !known[c]) {
                    missing_haplotypes.push_back(ploidy*individual + k);
                }
                alt_count += c == 1;
            }
            individual += 1;
        }
        if (!found.empty() and variant_site(alt_count, ploidy*individual, alleles)) {
            entries.insert(entries.end(), found.begin(), found.end());
            entry_sites.push_back(x);
            likelihood_genotypes += (int) found.size();
        }
        sites_without_likelihood += use_likelihood and !any_likelihood and missing_haplotypes.size() < missing_sites.size();
        if (missing_haplotypes.size() == missing_sites.size()) {
            unassayed_sites.push_back(x);
        } else if (missing_haplotypes.size() > 0) {
            any_missing = true;
            for (int h : missing_haplotypes) {
                missing_sites[h].push_back(x);
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
}
