#include <string>
#include <filesystem>
#include <vector>

#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/fmt.hpp>
#include <bio/ranges/container/concatenated_sequences.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>

#include <bio/io/txt/misc.hpp>
#include <bio/io/txt/reader.hpp>

enum class orientation : bool
{
    plus,
    minus
};

struct hash_string
{
    using is_transparent = void;

    std::size_t operator()(const std::string& v) const
    {
        return std::hash<std::string>{}(v);
    }
    std::size_t operator()(const char*v) const
    {
        return std::hash<std::string_view>{}(v);
    }

    std::size_t operator()(const std::string_view &v) const
    {
        return std::hash<std::string_view>{}(v);
    }
};

using seq_map_t = std::unordered_map<std::string, size_t, hash_string, std::equal_to<void>>;

struct graph
{
    bio::ranges::concatenated_sequences<std::string> stable_sequence_names;

    /* segment lines; all have same size */
    bio::ranges::concatenated_sequences<std::vector<bio::alphabet::dna5>> seqs;
    bio::ranges::concatenated_sequences<std::string> names;
    std::vector<uint64_t> stable_name_i; // -> stable_sequence_names
    std::vector<uint64_t> stable_offsets;
    std::vector<uint64_t> ranks;

    /* link lines; all have same size */
    std::vector<std::pair<uint64_t, uint64_t>> from_to; // -> segements
    std::vector<std::pair<orientation, orientation>> orient;
};



inline void parse_segment_optional(std::string_view const field,
                                   seq_map_t & seq_map,
                                   auto & checks,
                                   graph & output)
{
    if (field.starts_with("SN:Z:"))
    {
        std::string_view sn = field.substr(5);
        if (sn.empty())
            throw std::runtime_error{"SN field may not be empty."};

        if (auto it = seq_map.find(sn); it != seq_map.end())
        {
            output.stable_name_i.push_back(it->second);
        }
        else
        {
            output.stable_sequence_names.push_back(sn);
            seq_map[static_cast<std::string>(sn)] = output.stable_sequence_names.size() - 1;
        }

        checks.have_SN = true;
    }
    else if (field.starts_with("SO:i:"))
    {
        std::string_view so = field.substr(5);
        if (so.empty())
            throw std::runtime_error{"SO field may not be empty."};

        uint64_t offset = 0;
        auto [ptr, ec] = std::from_chars(so.data(), so.data() + so.size(), offset);
        if (ec != std::errc{})
            throw std::runtime_error{fmt::format("Expected number, got: {}", so)};

        output.stable_offsets.push_back(offset);
        checks.have_SO = true;
    }
    else if (field.starts_with("SR:i:"))
    {
        std::string_view sr = field.substr(5);
        if (sr.empty())
            throw std::runtime_error{"SR field may not be empty."};

        uint64_t rank = 0;
        auto [ptr, ec] = std::from_chars(sr.data(), sr.data() + sr.size(), rank);
        if (ec != std::errc{})
            throw std::runtime_error{fmt::format("Expected number, got: {}", sr)};

        output.ranks.push_back(rank);
        checks.have_SR = true;
    }
    // other optionals are ignored
}

inline void read_gfa(std::filesystem::path const & path, graph & output)
{

    auto reader = path == "-" ?
                  bio::io::txt::reader{std::cin, '\t', bio::io::txt::header_kind::starts_with{'#'}} :
                  bio::io::txt::reader{path, '\t', bio::io::txt::header_kind::starts_with{'#'}};

    seq_map_t stable_sequence_names_map;

    for (bio::io::txt::record & r : reader)
    {
        if (r.fields.empty())
            throw std::runtime_error{"Don't know how to handle empty line in rGFA."};

        if (r.fields[0].size() != 1)
            throw std::runtime_error{fmt::format("Don't know how to handle {} as first field.", r.fields[0])};

        switch (r.fields[0][0])
        {
            case 'S': /* SEGMENT */
            {
                if (r.fields.size() < 6)
                    throw std::runtime_error{"At least 6 fields required for SEGMENT record."};

                output.names.push_back(r.fields[1]);
                output.seqs.push_back(r.fields[2] | bio::views::char_strictly_to<bio::alphabet::dna5>);

                struct
                {
                    bool have_SN = false;
                    bool have_SO = false;
                    bool have_SR = false;
                } checks;

                for (size_t i = 3; i < r.fields.size(); ++i)
                    parse_segment_optional(r.fields[i], stable_sequence_names_map, checks, output);

                if (!checks.have_SN)
                    throw std::runtime_error{"rGFA requires SN optional but none was present."};
                if (!checks.have_SO)
                    throw std::runtime_error{"rGFA requires SO optional but none was present."};
                if (!checks.have_SR)
                    throw std::runtime_error{"rGFA requires SR optional but none was present."};

                assert(output.seqs.size() == output.names.size());
                assert(output.seqs.size() == output.stable_name_i.size());
                assert(output.seqs.size() == output.stable_offsets.size());
                assert(output.seqs.size() == output.ranks.size());
                break;
            }
            case 'L': /* LINKE */

                break;

            case '#': /* COMMENT */
                break;
            default:
                throw std::runtime_error{fmt::format("Don't know how to handle {} as first field.", r.fields[0])};
                break;
        }
    }
}

