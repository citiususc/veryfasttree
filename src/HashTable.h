
#ifndef VERYFASTTREE_HASHTABLE_H
#define VERYFASTTREE_HASHTABLE_H
#define XXH_INLINE_ALL

#include <vector>
#include <tsl/robin_map.h>
#include "DiskMemory.h"
#include <xxhash.h>

namespace veryfasttree {

    class HashTable {
    public:

        HashTable(const std::vector<std::string> &keys, bool checkUnique = false) {
            table.reserve(keys.size());
            if (checkUnique) {
                for (int64_t i = 0; i < (int64_t) keys.size(); i++) {
                    auto it = table.insert({&keys[i], i});
                    if (!it.second) {
                        throw std::invalid_argument("Non-unique name '" + keys[i] + "' in the alignment");
                    }
                }
            } else {
                for (int64_t i = keys.size(); i > 0; i--) {
                    table.insert_or_assign(&keys[i - 1], i - 1);
                }
            }
        }

        HashTable(const std::vector<std::string> &seqs, std::unique_ptr<DiskMemory> &disk) :
                tableDisk(0, HashDisk{disk.get()}, EqualDisk{disk.get()}) {
            for (int64_t i = seqs.size(); i > 0; i--) {
                tableDisk.insert_or_assign(&seqs[i - 1], i - 1);
            }
        }


        inline int64_t size() {
            return table.size() + tableDisk.size();
        }

        inline int64_t operator[](const std::string &key) const {
            if (!table.empty()) {
                return table.find(&key)->second;
            } else {
                return tableDisk.find(&key)->second;
            }
        }

        inline const int64_t *find(const std::string &key) const {
            if (!table.empty()) {
                auto entry = table.find(&key);
                if (entry != table.end()) {
                    return &(entry->second);
                }
            } else {
                auto entry = tableDisk.find(&key);
                if (entry != tableDisk.end()) {
                    return &(entry->second);
                }
            }
            return nullptr;
        }

        inline void clear() {
            table = {};
            tableDisk = {};
        }

    private:

        union HashResult {
            XXH64_hash_t x;
            int64_t r;
        };

        struct Hash {
            int64_t operator()(const std::string *p) const noexcept {
                return HashResult{XXH64(p->data(), p->size(), 0)}.r;
            }
        };

        struct Equal {
            bool operator()(const std::string *x, const std::string *y) const { return *x == *y; }
        };

        struct HashDisk {
            DiskMemory *disk;

            int64_t operator()(const std::string *p) const noexcept {
                std::string tmp = *p;
                size_t id;
                disk->load(tmp, id);
                return HashResult{XXH64(tmp.data(), tmp.size(), 0)}.r;
            }
        };

        struct EqualDisk {
            DiskMemory *disk;

            bool operator()(const std::string *x, const std::string *y) const {
                if (*x == *y) {
                    return true;
                }
                if (x->empty() || y->empty()) {
                    return false;
                }
                std::string tmpX = *x, tmpY = *y;
                size_t idX, idY;
                disk->load(tmpX, idX);
                disk->load(tmpY, idY);
                return tmpX == tmpY;
            }
        };

        tsl::robin_map<const std::string *, int64_t, Hash, Equal> table;
        tsl::robin_map<const std::string *, int64_t, HashDisk, EqualDisk> tableDisk;
    };
}

#endif
