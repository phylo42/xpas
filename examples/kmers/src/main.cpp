#include <iostream>
#include <core/seq.h>
#include <core/phylo_kmer.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>

/// Iterate over all the combinations with repetition
/// http://shoaib-ahmed.com/2018/for-each-combination-with-repetetion-c++/
template<typename V, typename Callable>
void for_each_combination(V &v, size_t gp_sz, Callable f) {
    V gp(gp_sz);
    auto total_n = std::pow(v.size(), gp.size());
    for (auto i = 0; i < total_n; ++i) {
        auto n = i;
        for (auto j = 0ul; j < gp.size(); ++j) {
            gp[gp.size() - j - 1] = v[n % v.size()];
            n /= v.size();
        }
        f(gp);
    }
}

void encode_decode(const std::string& kmer)
{
    const auto key = core::encode_kmer(kmer);
    std::cout << kmer << ": " << key << std::endl;
    assert(kmer == core::decode_kmer(key, kmer.size()));
}

int main()
{
    std::vector<char> alphabet = { 'A', 'C', 'G', 'T' };
    const size_t kmer_size = 3;

    for_each_combination(alphabet, kmer_size,
        [&](std::vector<char>& bases) {
            const auto kmer = std::string{ bases.begin(), bases.end() };
            encode_decode(kmer);
    });
}