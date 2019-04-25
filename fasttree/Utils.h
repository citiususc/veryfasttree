

#ifndef FASTTREE_UTILS_H
#define FASTTREE_UTILS_H

#include <iostream>


namespace fasttree {

    class TeeStream : public std::streambuf {
    public:
        TeeStream(std::ostream &os1, std::ostream &os2) : os1(os1), os2(os2) {}

        virtual std::streamsize xsputn(const char_type *__s, std::streamsize __n) {
            os1 << __s;
            os2 << __s;
            return __n;
        }


    private:
        std::ostream &os1;
        std::ostream &os2;
    };

    template<typename ... Args>
    std::string strformat(const std::string &format, Args ... args) {
        size_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1;//\0'
        char buf[size];
        snprintf(buf, size, format.c_str(), args ...);
        return std::string(buf, buf + size - 1);
    }

    template<typename _CharT, typename _Traits, typename _Alloc>
    inline bool readline(std::basic_istream<_CharT, _Traits> &__is, std::basic_string<_CharT, _Traits, _Alloc> &__str) {
        std::getline(__is, __str);
        if (__str.back() == '\r') {
            __str.resize(__str.size() - 1);
        }
        return !__is.eof();
    }

}
#endif
