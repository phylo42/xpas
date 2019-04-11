#ifndef RAPPAS_CORE_VERSION_H
#define RAPPAS_CORE_VERSION_H

namespace core
{
    /// \brief Core version class. Follows the semantic versioning 2.0 rules.
    /// \sa Semantic versioning 2.0: https://semver.org/
    struct version
    {
        static constexpr auto major = "0";
        static constexpr auto minor = "1";
        static constexpr auto revision = "0";
        static constexpr auto suffix = "alpha";

        /// \brief Returns the core version number as a string
        static constexpr std::string as_string() const
        {
            return major + "." + minor + "." + revision + (suffix ? "-" + suffix : "");
        }
    };
}

#endif //RAPPAS_CORE_VERSION_H
