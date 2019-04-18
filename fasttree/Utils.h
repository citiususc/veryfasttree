

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

}
#endif
