
#ifndef VERYFASTTREE_HASHTABLE_H
#define VERYFASTTREE_HASHTABLE_H

#include <vector>
#include <tsl/robin_map.h>
#include <xxh64.hpp>

namespace veryfasttree {

    class HashTable {
    public:

        HashTable(const std::vector<std::string> &keys, bool checkUnique = false) {
            table.reserve(keys.size());
            if (checkUnique) {
                for (int64_t i = 0; i < (int64_t)keys.size(); i++) {
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

        inline int64_t size() {
            return table.size();
        }

        inline int64_t operator[](const std::string &key) const {
            return table.find(&key)->second;
        }

        inline const int64_t *find(const std::string &key) const {
            auto entry = table.find(&key);
            if (entry != table.end()) {
                return &(entry->second);
            }
            return nullptr;
        }

    private:

        struct Hash {
            int64_t operator()(const std::string *p) const noexcept { return xxh64::hash(p->data(), p->size(), 0); }
        };

        struct Equal {
            bool operator()(const std::string *x, const std::string *y) const { return *x == *y; }
        };

        tsl::robin_map<const std::string *, int64_t, Hash, Equal> table;
    };
}

#endif
