

#ifndef FASTTREE_UTILS_H
#define FASTTREE_UTILS_H

typedef int8_t ibool;

#if (defined _WIN32 || defined WIN32 || defined WIN64 || defined _WIN64)

#include <io.h>
#include <ciso646>
#define FILE_SEP '\\'

namespace veryfasttree {
    inline bool isWindows(){return true;}

    inline bool isattyIn(){return ::_isatty( _fileno( stdin ) );}

    inline bool isattyOut(){return ::_isatty( _fileno( stdout ) );}

    inline bool isattyErr(){return ::_isatty( _fileno( stderr ) );}
}
#else

#include <unistd.h>
#define FILE_SEP '/'

namespace veryfasttree {
    inline bool isWindows() { return false; }

    inline bool isattyIn() { return ::isatty(STDIN_FILENO); }

    inline bool isattyOut() { return ::isatty(STDOUT_FILENO); }

    inline bool isattyErr() { return ::isatty(STDERR_FILENO); }
}
#endif

#include <iostream>
#include <chrono>
#include "cassert"
#include <omp.h>
#include <boost/sort/parallel_stable_sort/parallel_stable_sort.hpp>
#include <cmath>
#include <random>

namespace veryfasttree {
    class TeeStream : public std::streambuf {
    public:
        TeeStream(std::ostream &os1, std::ostream &os2) : os1(os1), os2(os2) {}

        virtual std::streamsize xsputn(const char_type *s, std::streamsize n) override {
            os1.write(s, n);
            os2.write(s, n);
            return n;
        }

        int sync() override {
            return os1.rdbuf()->pubsync() | os2.rdbuf()->pubsync();
        }

        int overflow(int c) {
            if (c != EOF) {
                os1 << (char) c;
                os2 << (char) c;
            }
            return c;
        }

    private:
        std::ostream &os1;
        std::ostream &os2;
    };

    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args &&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    template<typename T>
    inline T &asRef(T &value) {
        return value;
    }

    template<typename T>
    inline T &asRef(T *value) {
        return *value;
    }

    constexpr int alignsz(int sz, int aln) {
        return aln == 0 ? sz : ((sz / aln + ((sz % aln) != 0)) * aln);
    }

    template<typename T>
    inline void vecrelease(std::vector<T> &v) {
        v.resize(0);
        v.shrink_to_fit();
    }

    inline void strrelease(std::string &s) {
        if (s.size() > 0) {
            s.resize(0);
        }
        s.shrink_to_fit();
    }

    inline std::string randomString(std::size_t length) {
        const std::string chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
        std::string res(length, 0);

        std::random_device random_device;
        std::mt19937 generator(random_device());
        std::uniform_int_distribution<> distribution(0, (int) chars.size() - 1);

        for (std::size_t i = 0; i < length; ++i) {
            res[i] += chars[(int) distribution(generator)];
        }

        return res;
    }

    template<typename Iter, typename Compare>
    inline void psort(Iter first, Iter last, const Compare &comp, bool force_parallel = false) {
        uint32_t threads = 1;
        if (!omp_in_parallel() || force_parallel) {
            threads = omp_get_num_threads();
        }
        #ifndef NDEBUG
        static thread_local const Compare *_pcomp = nullptr;
        _pcomp = &comp;
        auto cmp = [](const void *a, const void *b) {
            if ((*_pcomp)(*(typename Iter::pointer) a, *(typename Iter::pointer) b)) {
                return -1;
            } else {
                return 1;
            }
        };
        qsort(&(*first), std::distance(first, last), sizeof(*first), cmp);
        return;
        #endif
        boost::sort::parallel_stable_sort(first, last, comp, threads);
    }

    template<typename ... Args>
    inline std::string strformat(const std::string &format, Args ... args) {
        int64_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1;//\0'
        std::string buf;
        buf.resize(size);
        snprintf(&buf.front(), size, format.c_str(), args ...);
        buf.resize(size - 1);
        return buf;
    }

    inline std::string strformat(const std::string &format) {
        return std::string(format);
    }

    inline std::ostream &quotes(std::ostream &s, const std::string str, bool bQuote) {
        if (bQuote) {
            s << '\'' << str << '\'';
        } else {
            s << str;
        }
        return s;
    }

    template<typename _CharT, typename _Traits, typename _Alloc>
    inline bool readline(std::basic_istream<_CharT, _Traits> &__is, std::basic_string<_CharT, _Traits, _Alloc> &__str) {
        if (!__is.eof()) {
            std::getline(__is, __str);
            if (!__str.empty() && __str.back() == '\r') {
                __str.resize(__str.size() - 1);
            }
            return true;
        }
        return false;
    }

    class ProgressReport {

    public:

        typedef std::chrono::high_resolution_clock Clock;

        ProgressReport(bool showProgress, int verbose) : clockStart(Clock::now()), timeLast(Clock::now()),
        showProgress(showProgress) , verbose(verbose){}

        inline double clockDiff() {
            auto timeNow = Clock::now();
            return std::chrono::duration_cast<std::chrono::milliseconds>(timeNow - clockStart).count() /
                   (double) (1000);
        }

        template<typename ... Args>
        inline void print(const std::string &format, Args ... args) {
            if (!showProgress) {
                return;
            }

            auto timeNow = Clock::now();
            int64_t mili = std::chrono::duration_cast<std::chrono::milliseconds>(timeNow - timeLast).count();

            if (mili > 100 || verbose > 1) {
                std::cerr << strformat("%7i.%2.2i seconds: ", (int) (mili / 1000), (int) ((mili % 1000) / 10));
                std::cerr << strformat(format, args...);
                if (verbose > 1 || !isattyErr()) {
                    std::cerr << std::endl;
                } else {
                    std::cerr << "   \r" << std::flush;
                }
                timeLast = timeNow;
            }
        }

    private:
        const Clock::time_point clockStart;
        Clock::time_point timeLast;
        bool showProgress;
        int verbose;

    };

}
#endif
