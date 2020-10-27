#ifndef XPAS_PHYLO_KMER_H
#define XPAS_PHYLO_KMER_H

#include "seq.h"
#include <limits>
#include <string>
#include "optional.h"

namespace xpas
{
    /// \brief A phylo k-mer structure.
    /// \details A key-value pair for a phylo-kmer, where key is a key_type value of a k-mer, and value is
    /// posterior probability score of this k-mer. Branch node id, position etc. omitted here, because these values
    /// are shared among multiple phylo k-mers and can be stored more effectively.
    struct phylo_kmer
    {
        /// The same type used to store a k-mer value
        /// (essentially we do not distinguish between "k-mer value" and "phylo-kmer value")
        using key_type = seq_traits::key_type;

        /// The type of a "posterior probability" score of a phylokmer
        using score_type = float;

        /// The type of a branch node id, that a phylokmer is mapped to.
        using branch_type = uint32_t;

        /// The type of a phylokmer's position in the alignment
        using pos_type = int32_t;

        static constexpr key_type nan_key = std::numeric_limits<phylo_kmer::key_type>::max();
        static constexpr score_type nan_score = std::numeric_limits<phylo_kmer::score_type>::quiet_NaN();
        static constexpr branch_type nan_branch = std::numeric_limits<phylo_kmer::branch_type>::max();

        bool is_nan() const;

        key_type key;
        score_type score;
    };

    bool operator==(const phylo_kmer& lhs, const phylo_kmer& rhs) noexcept;

    /// Returns a minumum score
    phylo_kmer::score_type score_threshold(phylo_kmer::score_type omega, size_t kmer_size);

    /// \brief Returns one or more codes of the input k-mer (depending on the policy)
    /// \details Assumes that the size of input sequence equals k
    template<typename AmbiguityPolicy>
    optional<typename AmbiguityPolicy::value_type> encode_kmer(std::string_view kmer);

    template<typename AmbiguityPolicy>
    optional<phylo_kmer::key_type> encode_kmer(const std::string& kmer)
    {
        return encode_kmer<AmbiguityPolicy>(std::string_view{ kmer });
    }

    /// Creates a string of size kmer_size by given key
    std::string decode_kmer(phylo_kmer::key_type key, size_t kmer_size);
}

#endif