
#ifndef FASTTREE_HASHTABLE_H
#define FASTTREE_HASHTABLE_H

#include <vector>
#include "tsl/robin_map.h"

namespace fasttree {

    class HashTable {
    public:

        HashTable(std::vector<std::string> &seqs) {
            table.reserve(seqs.size());
            for (size_t i = 0; i < seqs.size(); i++) {
                auto it = table.insert({&seqs[i], i});
                if (!it.second) {
                    throw std::invalid_argument("Non-unique name '" + seqs[i] + "' in the alignment");
                }
            }
        }

    private:

        struct Hash {
            std::hash<std::string> imp;

            size_t operator()(std::string *p) const noexcept { return imp(*p); }
        };

        struct Equal {
            bool operator()(const std::string *x, const std::string *y) const { return *x == *y; }
        };

        tsl::robin_map<std::string *, size_t, Hash, Equal> table;
    };
}

#endif
