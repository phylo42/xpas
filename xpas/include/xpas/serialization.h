#ifndef XPAS_SERIALIZATION_H
#define XPAS_SERIALIZATION_H

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <fstream>
#include "phylo_kmer_db.h"

namespace xpas
{
    static const unsigned int protocol_version = 2;

    xpas::phylo_kmer_db load(const std::string& filename)
    {
        std::ifstream ifs(filename);
        boost::archive::binary_iarchive ia(ifs);

        xpas::phylo_kmer_db db { 0, 0.0, 1.0, "" };
        ia & db;
        return db;
    }

    void save(const xpas::phylo_kmer_db& db, const std::string& filename)
    {
        std::ofstream ofs(filename);
        boost::archive::binary_oarchive oa(ofs);
        oa & db;
    }
}

namespace boost::serialization
{
    template<class Archive>
    inline void save(Archive& ar, const xpas::phylo_kmer_db& db, unsigned int /*version*/)
    {
        const auto original_tree_view = std::string{ db.tree() };
        ar & original_tree_view;

        size_t kmer_size = db.kmer_size();
        ar & kmer_size;

        xpas::phylo_kmer::score_type omega = db.omega();
        ar & omega;

        const auto f = db.f();
        ar & f;

        size_t table_size = db.size();
        ar & table_size;

        for (const auto& [key, entries] : db)
        {
            size_t entries_size = entries.size();
            ar & key & entries_size;
            for (const auto& [branch, score] : entries)
            {
                ar & branch & score;
            }
        }

        size_t bernoulli_size = db.num_branches();
        ar & bernoulli_size;

        for (size_t i = 0; i < bernoulli_size; ++i)
        {
            const auto score_sum = db.get_sum(i);
            ar & score_sum;
        }
    }

    template<class Archive>
    inline void load(Archive& ar, xpas::phylo_kmer_db& db, unsigned int version)
    {
        if (version < xpas::protocol_version)
        {
            throw std::runtime_error("Failed to load database: this database was built with older version of RAPPAS.");
        }
        std::string tree = "";
        ar & tree;
        db._tree = tree;

        size_t kmer_size = 0;
        ar & kmer_size;
        db._kmer_size = kmer_size;

        xpas::phylo_kmer::score_type omega = 0;
        ar & omega;
        db._omega = omega;

        double f = 1.0;
        ar & f;
        db._f = f;

        size_t table_size = 0;
        ar & table_size;
        for (size_t i = 0; i < table_size; ++i)
        {
            xpas::phylo_kmer::key_type key = xpas::phylo_kmer::nan_key;
            size_t entries_size = 0;
            ar & key;
            ar & entries_size;
            for (size_t j = 0; j < entries_size; ++j)
            {
                xpas::phylo_kmer::branch_type branch = xpas::phylo_kmer::nan_branch;
                xpas::phylo_kmer::score_type score = xpas::phylo_kmer::nan_score;
                ar & branch & score;
                db.insert(key, { branch, score });
            }
        }

        size_t bernoulli_sum_size = 0;
        ar & bernoulli_sum_size;
        db.set_sum(bernoulli_sum_size, 0);

        for (size_t i = 0; i < bernoulli_sum_size; ++i)
        {
            xpas::phylo_kmer::score_type sum = 0.0;
            ar & sum;
            db.set_sum(i, sum);
        }
    }

    // split non-intrusive serialization function member into separate
    // non intrusive save/load member functions
    template<class Archive>
    inline void serialize(Archive& ar, xpas::phylo_kmer_db& db, unsigned int file_version)
    {
        boost::serialization::split_free(ar, db, file_version);
    }
}

BOOST_CLASS_VERSION(xpas::phylo_kmer_db, xpas::protocol_version)

#endif //RAPPAS_CORE_SERIALIZATION_H
