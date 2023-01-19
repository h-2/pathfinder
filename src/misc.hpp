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

void assign_append(std::ranges::forward_range auto &&                                                  source,
                   std::ranges::output_range<std::ranges::range_reference_t<decltype(source)>> auto && sink)
{
    size_t old_sink_size = sink.size();
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

inline void append_region(std::vector<bio::io::genomic_region> &                vec,
                          bio::meta::decays_to<bio::io::genomic_region> auto && reg)
{
    assert(reg.beg != reg.end);

    if (!vec.empty() &&                                               // there is previous region
        vec.back().chrom == reg.chrom && vec.back().end == reg.beg && // they touch
        ((vec.back().beg <= vec.back().end) == (reg.beg <= reg.end))) // same orientation
    {
        vec.back().end = reg.end;
    }
    else
    {
        vec.push_back(std::forward<decltype(reg)>(reg));
    }
}
