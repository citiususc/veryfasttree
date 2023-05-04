

#ifndef VERYFASTTREE_DISKMEMORY_H
#define VERYFASTTREE_DISKMEMORY_H

#include <fcntl.h>
#include <memory>
#include <sys/mman.h>
#include <string>

namespace veryfasttree {
    class DiskMemory {
    public:

        DiskMemory(const std::string &path, const std::string &name, size_t size);

        uintptr_t ptr();

        size_t getSize();

        ~DiskMemory();

        size_t store(size_t offset, std::string &s);

        void load(std::string &s, size_t &id);

        void release(std::string &s, size_t &id, bool update = false);

    private:
        std::string file;
        size_t size;
        uintptr_t map;
    };

    class DynDiskMemory {
    public:
        DynDiskMemory(const std::string &path, const std::string &name);

        template<typename T>
        T *allocate(T *src, size_t n) { return (T *) menAllocate(src, n * sizeof(T)); }

        template<typename T>
        T *alignAllocate(T *src, size_t n, size_t alg) {
            size_t sz = n * sizeof(T) + alg;
            void *mem = (void *) menAllocate(src, sz);
            return (T *) std::align(alg, sizeof(T), mem, sz);
        }

        void release();

        bool inAlloc();

        const std::string &getName();

    private:
        uintptr_t menAllocate(void *src, size_t size);

        const std::string path;
        const std::string name;
        bool allocated;
        std::unique_ptr<DiskMemory> disk;
    };
}


#endif //VERYFASTTREE_DISKMEMORY_H
