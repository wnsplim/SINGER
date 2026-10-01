//
//  Vcf_reader.hpp
//  SINGER
//

#ifndef Vcf_reader_hpp
#define Vcf_reader_hpp

#include <string>
#include <sys/stat.h>
#include <zlib.h>

class Vcf_reader {

public:

    Vcf_reader(const std::string &prefix) {
        struct stat st;
        path = stat((prefix + ".vcf.gz").c_str(), &st) == 0 ? prefix + ".vcf.gz" : prefix + ".vcf";
        file = gzopen(path.c_str(), "rb");
    }

    ~Vcf_reader() {
        if (file != nullptr) {
            gzclose(file);
        }
    }

    bool is_open() const {
        return file != nullptr;
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

    gzFile file = nullptr;
    std::string pending;
    bool has_pending = false;

    bool read(std::string &line) {
        line.clear();
        if (file == nullptr) {
            return false;
        }
        char buffer[65536];
        while (gzgets(file, buffer, sizeof(buffer)) != nullptr) {
            line += buffer;
            if (!line.empty() and line.back() == '\n') {
                line.pop_back();
                if (line.empty()) {
                    continue;
                }
                return true;
            }
        }
        return !line.empty();
    }
};

#endif
