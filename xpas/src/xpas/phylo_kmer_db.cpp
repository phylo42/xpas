#include "xpas/phylo_kmer_db.h"

using namespace xpas;

phylo_kmer_db::phylo_kmer_db(size_t kmer_size, xpas::phylo_kmer::score_type omega, const std::string& tree)
    : _kmer_size{ kmer_size }, _omega{ omega }, _tree { tree }
{}

void phylo_kmer_db::insert(key_type key, const pkdb_value& value)
{
    _map[key].push_back(value);
}

void phylo_kmer_db::set_sum(phylo_kmer::branch_type branch_id, phylo_kmer::score_type value)
{
    if (branch_id >= _bernoulli_sums.size())
    {
        const auto total_kmers_num = std::pow(xpas::seq_traits::alphabet_size, _kmer_size);
        const auto threshold = xpas::score_threshold(_omega, _kmer_size);
        double bernoulli_score_sum = total_kmers_num * std::log(1 - threshold);

        for (size_t i = _bernoulli_sums.size(); i <= branch_id; ++i)
        {
            _bernoulli_sums.push_back(bernoulli_score_sum);
        }
    }
    _bernoulli_sums[branch_id] = value;
}

phylo_kmer::score_type phylo_kmer_db::get_sum(phylo_kmer::branch_type branch_id) const
{
    return _bernoulli_sums[branch_id];
}

phylo_kmer::branch_type phylo_kmer_db::num_branches() const
{
    return _bernoulli_sums.size();
}

std::optional<impl::search_result> phylo_kmer_db::search(key_type key) const noexcept
{
    if (auto it = _map.find(key); it != _map.end())
    {
        return impl::search_result{ it->second.begin(), it->second.end() };
    }
    else
    {
        return std::nullopt;
    }
}

phylo_kmer_db::const_iterator phylo_kmer_db::begin() const noexcept
{
    return std::begin(_map);
}

phylo_kmer_db::const_iterator phylo_kmer_db::end() const noexcept
{
    return std::end(_map);
}

size_t phylo_kmer_db::size() const noexcept
{
    return _map.size();
}

size_t phylo_kmer_db::kmer_size() const noexcept
{
    return _kmer_size;
}

xpas::phylo_kmer::score_type phylo_kmer_db::omega() const noexcept
{
    return _omega;
}

std::string_view phylo_kmer_db::tree() const noexcept
{
    return _tree;
}

phylo_kmer_db::storage::hasher phylo_kmer_db::hash_function() const noexcept
{
    return _map.hash_function();
}

impl::search_result::search_result(
    impl::search_result::const_iterator begin,
    impl::search_result::const_iterator end) noexcept
    : _begin{ begin }, _end{ end }
{}

impl::search_result::const_iterator impl::search_result::begin() const noexcept
{
    return _begin;
}

impl::search_result::const_iterator impl::search_result::end() const noexcept
{
    return _end;
}