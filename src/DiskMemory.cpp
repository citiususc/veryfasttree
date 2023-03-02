
#include "DiskMemory.h"

#include "Utils.h"
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>

using namespace veryfasttree;

DiskMemory::DiskMemory(const std::string &name, size_t size) : name(name), size(size), file(-1), map(0) {
    file = open(name.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
    if (file == -1) {
        throw std::runtime_error("disk memory path is invalid");
    }
    lseek(file, size, SEEK_SET);
    if (write(file, "", 1) == -1) {
        throw std::runtime_error("disk memory truncation error");
    }
    void *men = mmap(0, size, PROT_READ | PROT_WRITE, MAP_PRIVATE, file, 0);
    if (men == MAP_FAILED) {
        throw std::runtime_error("memory mapping fails");
    }
    map = (uintptr_t) men;
}

uintptr_t DiskMemory::get() {
    return map;
}

DiskMemory::~DiskMemory() {
    if (map != 0) {
        munmap((void *) map, size);
    }
    if (file != -1) {
        close(file);
    }
    remove(name.c_str());
}

size_t DiskMemory::store(size_t offset, std::string &s) {
    if (!s.empty()) {
        size_t disp = 0;
        size_t len;
        const char *init;
        if (s[0] == 0) {
            disp = *((size_t *) (s.c_str() + 1 + sizeof(size_t))) - 1;
            len = s.size() - sizeof(size_t) * 2;
            init = s.c_str() + 1 + sizeof(size_t) * 2;
        } else {
            len = s.size() + 1;
            init = s.c_str();
        }
        if (offset + len >= this->size) {
            throw std::runtime_error("input file has a invalid size for disk computing. If the file is compressed and "
                                     "is not a header file format, it must be decompressed.");
        }
        if(len > 1){
            memcpy((void *) (map + offset + disp), init, len);
            s.resize(sizeof(size_t) * 2 + 1);
            s[0] = 0;
            *((size_t *) (s.c_str() + 1)) = offset;
            *((size_t *) (s.c_str() + 1 + sizeof(size_t))) = disp + len;
        }
        offset += disp + len;
    }
    return offset;
}

void DiskMemory::load(std::string &s, size_t &id) {
    if (!s.empty() && s[0] == 0) {
        id = *((size_t *) (s.c_str() + 1));
        size_t len = *((size_t *) (s.c_str() + 1 + sizeof(size_t))) - 1;
        s.resize(len);
        memcpy((void *) s.c_str(), (void *) (map + id), s.size());
    } else {
        id = SIZE_MAX;
    }
}

void DiskMemory::release(std::string &s, size_t &id, bool update) {
    if (id == SIZE_MAX) {
        return;
    }
    if (update) {
        memcpy((void *) (map + id), (void *) s.c_str(), s.size());
    }
    s[0] = 0;
    *((size_t *) (s.c_str() + 1)) = id;
    *((size_t *) (s.c_str() + 1 + sizeof(size_t))) = s.size() + 1;
    s.resize(sizeof(size_t) * 2 + 1);
    s.shrink_to_fit();
}
