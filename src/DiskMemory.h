

#ifndef VERYFASTTREE_DISKMEMORY_H
#define VERYFASTTREE_DISKMEMORY_H

#include <fcntl.h>
#include <sys/mman.h>
#include <string>

namespace veryfasttree {
    class DiskMemory {
    public:

        DiskMemory(const std::string &name, size_t size);

        uintptr_t get();

        ~DiskMemory();

        size_t store(size_t offset, std::string &s);

        void load(std::string &s, size_t &id);

        void release(std::string &s, size_t &id, bool update = false);

    private:
        std::string name;
        size_t size;
        int file;
        uintptr_t map;
    };
}


#endif //VERYFASTTREE_DISKMEMORY_H
