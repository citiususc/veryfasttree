
#ifndef FASTTREE_VALIDATORS_H
#define FASTTREE_VALIDATORS_H

#include "CLI11.hpp"
#include <iostream>

namespace CLI {

    std::string isNotEmpty(const std::string &input) {
        if (input.size() == 0) {
            return "Value is empty";
        }
        return std::string();
    }

    struct Min : Validator {

        template<typename T>
        Min(T min, bool equals = true) {
            func = [min,equals](std::string input) {
                T val;
                CLI::detail::lexical_cast(input, val);
                if (val < min)
                    return "Min value " + std::to_string(min);

                if (!equals && val == min)
                    return "Min value greater than " + std::to_string(min);

                return std::string();
            };
        }

    };

    void deprecated(Option* op){
        op->check([&op](const std::string& in){
            std::cerr << "Warning: "<< op->get_name() <<" is a deprecated option and it has no effect" << std::endl;
            return "";
        });
    }


}

#endif
