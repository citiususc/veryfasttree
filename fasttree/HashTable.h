
#ifndef FASTTREE_HASHTABLE_H
#define FASTTREE_HASHTABLE_H

#include <vector>
#include "tsl/robin_map.h"
#include "xxh64.hpp"

namespace fasttree {

    class HashTable {
    public:

        HashTable(const std::vector<std::string> &keys, bool checkUnique = false) {
            table.reserve(keys.size());
            if (checkUnique) {
                for (size_t i = 0; i < keys.size(); i++) {
                    auto it = table.insert({&keys[i], i});
                    if (!it.second) {
                        throw std::invalid_argument("Non-unique name '" + keys[i] + "' in the alignment");
                    }
                }
            } else {
                for (size_t i = keys.size(); i > 0; i--) {
                    table.insert_or_assign(&keys[i - 1], i - 1);
                }
            }
        }

        inline size_t size() {
            return table.size();
        }

        inline size_t operator[](const std::string &key) const {
            return table.find(&key)->second;
        }

        inline const size_t *find(const std::string &key) const {
            auto entry = table.find(&key);
            if (entry != table.end()) {
                return &(entry->second);
            }
            return nullptr;
        }

    private:

        struct Hash {
            size_t operator()(const std::string *p) const noexcept { return xxh64::hash(p->data(), p->size(), 0); }
        };

        struct Equal {
            bool operator()(const std::string *x, const std::string *y) const { return *x == *y; }
        };

        tsl::robin_map<const std::string *, size_t, Hash, Equal> table;
    };
}

#endif
