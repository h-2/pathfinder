#pragma once

#include <ranges>
#include <string>

#include <bio/alphabet/fmt.hpp>
#include <bio/io/genomic_region.hpp>

inline constexpr auto ssize = std::ranges::ssize;

struct hash_string
{
    using is_transparent = void;

    std::size_t operator()(std::string const & v) const { return std::hash<std::string>{}(v); }
    std::size_t operator()(char const * v) const { return std::hash<std::string_view>{}(v); }
    std::size_t operator()(std::string_view const & v) const { return std::hash<std::string_view>{}(v); }
};

struct hash_regs
{
    std::size_t operator()(std::span<bio::io::genomic_region const> regs) const
    {
        size_t ret = 0;
        for (bio::io::genomic_region const reg : regs)
        {
            ret ^= std::hash<std::string>{}(reg.chrom);
            ret ^= std::hash<int64_t>{}(reg.beg);
            ret ^= std::hash<int64_t>{}(reg.end);
        }

        return ret;
    }
};

void assign_append(std::ranges::forward_range auto &&                                                  source,
                   std::ranges::output_range<std::ranges::range_reference_t<decltype(source)>> auto && sink)
{
    size_t const old_sink_size = sink.size();
    sink.resize(old_sink_size + source.size());
    std::ranges::copy(source, sink.begin() + old_sink_size);
}

template <>
struct fmt::formatter<bio::io::genomic_region> : fmt::formatter<std::string>
{
    constexpr auto format(bio::io::genomic_region const & reg, auto & ctx) const
    {
        std::string buf;
        fmt::format_to(std::back_insert_iterator{buf}, "{}:{}-{}", reg.chrom, reg.beg, reg.end);
        return fmt::formatter<std::string>::format(buf, ctx);
    }
};

inline void append_region(bio::meta::decays_to<bio::io::genomic_region> auto && reg,
                          std::vector<bio::io::genomic_region> &                dest)
{
    assert(reg.beg != reg.end);

    if (!dest.empty() &&                                                // there is previous region
        dest.back().chrom == reg.chrom && dest.back().end == reg.beg && // they touch
        ((dest.back().beg <= dest.back().end) == (reg.beg <= reg.end))) // same orientation
    {
        dest.back().end = reg.end;
    }
    else
    {
        dest.push_back(std::forward<decltype(reg)>(reg));
    }
}

inline void append_regions(bio::meta::decays_to<std::vector<bio::io::genomic_region>> auto && source,
                           std::vector<bio::io::genomic_region> &                             dest)
{
    for (bio::meta::decays_to<bio::io::genomic_region> auto && reg : source)
    {
        if constexpr (std::is_lvalue_reference_v<decltype(source)>)
            append_region(reg, dest);
        else
            append_region(std::move(reg), dest);
    }
}
