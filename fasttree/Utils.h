

#ifndef FASTTREE_UTILS_H
#define FASTTREE_UTILS_H

#include <iostream>
#include <chrono>

#if (defined _WIN32 || defined WIN32 || defined WIN64 || defined _WIN64)

#include <io.h>

namespace fasttree {
    inline bool isWindows(){return true;}

    inline bool isattyIn(){return ::_isatty( _fileno( stdin ) );}

    inline bool isattyOut(){return ::_isatty( _fileno( stdout ) );}

    inline bool isattyErr(){return ::_isatty( _fileno( stderr ) );}
}
#else

#include <unistd.h>

namespace fasttree {
    inline bool isWindows() { return false; }

    inline bool isattyIn() { return ::isatty(STDIN_FILENO); }

    inline bool isattyOut() { return ::isatty(STDOUT_FILENO); }

    inline bool isattyErr() { return ::isatty(STDERR_FILENO); }
}
#endif

namespace fasttree {
    class TeeStream : public std::streambuf {
    public:
        TeeStream(std::ostream &os1, std::ostream &os2) : os1(os1), os2(os2) {}

        virtual std::streamsize xsputn(const char_type *__s, std::streamsize __n) override{
            os1 << __s;
            os2 << __s;
            return __n;
        }

        int sync() override {
            return os1.rdbuf()->pubsync() | os2.rdbuf()->pubsync();
        }

        int overflow(int c ) {
            if (c != EOF) {
                os1 << (char)c;
                os2 << (char)c;
            }
            return c;
        }

    private:
        std::ostream &os1;
        std::ostream &os2;
    };

    template<typename ... Args>
    inline std::string strformat(const std::string &format, Args ... args) {
        size_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1;//\0'
        std::string buf;
        buf.resize(size);
        snprintf(&buf.front(), size, format.c_str(), args ...);
        buf.resize(size - 1);
        return std::move(buf);
    }

    inline std::string strformat(const std::string &format) {
        return std::string(format);
    }

    template<typename _CharT, typename _Traits, typename _Alloc>
    inline bool readline(std::basic_istream<_CharT, _Traits> &__is, std::basic_string<_CharT, _Traits, _Alloc> &__str) {
        std::getline(__is, __str);
        if (__str.back() == '\r') {
            __str.resize(__str.size() - 1);
        }
        return !__is.eof();
    }

    class ProgressReport {

    public:

        typedef std::chrono::high_resolution_clock Clock;
        const Clock::time_point clockStart;

        ProgressReport(const Options &options) : clockStart(Clock::now()), options(options) {}

        template<typename ... Args>
        inline void print(const std::string &format, Args ... args) {
            if (!options.showProgress) {
                return;
            }

            auto timeNow = Clock::now();

            if (std::chrono::duration_cast<std::chrono::seconds>(timeNow - timeLast).count() > 1) {
                size_t mili = std::chrono::duration_cast<std::chrono::milliseconds>(timeNow - clockStart).count();
                std::cerr << strformat("%7i.%2.2i seconds: ", (int) (mili / 1000), (int) ((mili % 1000) / 10));
                std::cerr << strformat(format, args...);
                if (options.verbose > 1 || !isattyErr()) {
                    std::cerr << std::endl;
                } else {
                    std::cerr << "   \r" << std::flush;
                }
                timeLast = timeNow;
            }
        }

    private:

        Clock::time_point timeLast;
        const Options &options;

    };

}
#endif
