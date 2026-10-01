//
//  Vcf_reader.hpp
//  SINGER
//

#ifndef Vcf_reader_hpp
#define Vcf_reader_hpp

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <zlib.h>
#include <htslib/vcf.h>

class Vcf_reader {

public:

    static void fail(const std::string &message) {
        std::cerr << "Error: " << message << ". " << std::endl;
        exit(1);
    }

    static bool is_bcf(const std::string &p) {
        return p.size() > 4 and p.compare(p.size() - 4, 4, ".bcf") == 0;
    }

    Vcf_reader(const std::string &input, const std::string &chrom = "") : path(input), chrom(chrom) {
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

    bool next(std::string &line) {
        if (has_pending) {
            line = std::move(pending);
            has_pending = false;
            return true;
        }
        return read(line);
    }

    bool peek(std::string &line) {
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

    std::string path;

private:

    std::string chrom;
    bool seen = false;
    bool done = false;
    gzFile file = nullptr;
    htsFile *bcf = nullptr;
    bcf_hdr_t *header = nullptr;
    bcf1_t *record = nullptr;
    kstring_t text = {0, 0, nullptr};
    kstring_t buffer = {0, 0, nullptr};
    size_t at = 0;
    std::string pending;
    bool has_pending = false;

    bool read(std::string &line) {
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

    bool read_line(std::string &line) {
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

#endif
