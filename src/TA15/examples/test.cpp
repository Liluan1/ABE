#include "abecontext.h"

using namespace lbcrypto;

int main(int argc, char const *argv[]) {
    typedef NativePoly::Integer Integer;
    Integer a = 10;
    for (usint i = 2; i < a; ++i) {
        if (GreatestCommonDivisor(Integer(i), a) == 1) {
            std::cout<<i<<std::endl;
            break;
        }
    }
    // std::cout<<a<<std::endl;
    return 0;
}



