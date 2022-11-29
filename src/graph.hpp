#include <bits/ranges_algo.h>
#include <bits/ranges_base.h>
#include <ranges>
#include <string>
#include <filesystem>
#include <vector>

#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/fmt.hpp>
#include <bio/ranges/container/concatenated_sequences.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>
#include <bio/ranges/views/complement.hpp>

#include <bio/io/txt/misc.hpp>
#include <bio/io/txt/reader.hpp>

/* DETAIL*/

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

void assign_append(std::ranges::forward_range auto && source,
                   std::ranges::output_range<std::ranges::range_reference_t<decltype(source)>> auto && sink)
{

    size_t old_sink_size = sink.size();
    sink.resize(old_sink_size+ source.size());
    std::ranges::copy(source, sink.begin() + old_sink_size);
}

/* PROGRAM */

enum class orientation : bool
{
    plus,
    minus
};

struct arc
{
    uint64_t target_i         : 62 = 0ull;
    orientation orient_self   : 1  = orientation::plus;
    orientation orient_target : 1  = orientation::plus;
};



struct graph
{
    bio::ranges::concatenated_sequences<std::string> stable_sequence_names;

    /* segment lines; all have same size */
    bio::ranges::concatenated_sequences<std::vector<bio::alphabet::dna5>> seqs;
    bio::ranges::concatenated_sequences<std::string> names;
    std::vector<uint64_t> stable_name_i; // -> stable_sequence_names
    std::vector<uint64_t> stable_offsets;
    std::vector<uint64_t> ranks;

    // bio::ranges::concatenated_sequences<std::vector<arc>> arcs;
    std::vector<std::vector<arc>> arcs;
    std::vector<uint64_t> split_counter;
};


using seq_map_t = std::unordered_map<std::string, size_t, hash_string, std::equal_to<void>>;
using arc_map_t = std::unordered_map<std::string, std::vector<arc>, hash_string, std::equal_to<void>>;

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

    seq_map_t segment_names_map;
    seq_map_t stable_sequence_names_map;
    arc_map_t arc_map;

    for (bio::io::txt::record & r : reader)
    {
        auto error = [&r] <typename ... arg_ts> (fmt::format_string<arg_ts...> fmt_str, arg_ts ... args)
        {
            std::string s = "rGFA parse error.\n" + fmt::format(fmt_str, std::forward<arg_ts>(args)...);
            s += "\nSource line:\n";
            s += r.line;
            throw std::runtime_error{std::move(s)};
        };

        if (r.fields.empty())
            error("Don't know how to handle empty line in rGFA.");

        if (r.fields[0].size() != 1)
            error("Don't know how to handle {} as first field.", r.fields[0]);

        switch (r.fields[0][0])
        {
            case 'S': /* SEGMENT */
            {
                if (r.fields.size() < 6)
                    error("At least 6 fields required for SEGMENT record.");

                if (auto it = segment_names_map.find(r.fields[1]); it != segment_names_map.end())
                {
                    error("Segment names have to be unique, but this one isn't.");
                }
                else
                {
                    output.names.push_back(r.fields[1]);
                    segment_names_map[static_cast<std::string>(r.fields[1])] = output.names.size() - 1;
                }

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
                    error("rGFA requires SN optional but none was present.");
                if (!checks.have_SO)
                    error("rGFA requires SO optional but none was present.");
                if (!checks.have_SR)
                    error("rGFA requires SR optional but none was present.");

                assert(output.seqs.size() == output.names.size());
                assert(output.seqs.size() == output.stable_name_i.size());
                assert(output.seqs.size() == output.stable_offsets.size());
                assert(output.seqs.size() == output.ranks.size());
                break;
            }
            case 'L': /* LINKE */
            {
                if (r.fields.size() < 6)
                    error("At least 6 fields required for LINK record.");

                std::string_view source_seg_id     = r.fields[1];
                std::string_view source_seg_orient = r.fields[2];
                std::string_view target_seg_id     = r.fields[3];
                std::string_view target_seg_orient = r.fields[4];
                std::string_view cigar             = r.fields[5];

                arc a;

                if (auto it = segment_names_map.find(target_seg_id); it == segment_names_map.end())
                    error("Target segment ID '{}' unknown.", target_seg_id);
                else
                    a.target_i = it->second;

                if (source_seg_orient == "+")
                    a.orient_self = orientation::plus;
                else if (source_seg_orient == "-")
                    a.orient_self = orientation::minus;
                else
                    error("Source orientation '{}' invalid, must be '+' or '-'.", source_seg_orient);

                if (target_seg_orient == "+")
                    a.orient_target = orientation::plus;
                else if (target_seg_orient == "-")
                    a.orient_target = orientation::minus;
                else
                    error("Source orientation '{}' invalid, must be '+' or '-'.", target_seg_orient);

                if (!segment_names_map.contains(source_seg_id))
                    error("Source segment ID '{}' unknown.", source_seg_id);
                else
                    arc_map[static_cast<std::string>(source_seg_id)].push_back(a);

                if (cigar != "0M")
                    error("Overlaps are not supported.");

                break;
            }
            case '#': /* COMMENT */
                break;
            default:
                error("Don't know how to handle {} as first field.", r.fields[0]);
                break;
        }
    }

    for (std::string_view seg_name : std::as_const(output.names))
    {
        if (auto it = arc_map.find(seg_name); it == arc_map.end())
            output.arcs.push_back({});
        else
            output.arcs.push_back(it->second);
    }
}

enum class reg_placement
{
    ends_before,
    ends_in,
    lies_in,
    contains,
    begins_in,
    begins_after
};

inline reg_placement overlap(std::pair<size_t, size_t> region1, std::pair<size_t, size_t> region2)
{
    //TODO this can probably be optimised with arithmetics
    if (region1.second <= region2.first)
        return reg_placement::ends_before;
    else if (region1.first >= region2.second)
        return reg_placement::begins_after;
    else if (region1.first < region2.first)
        return (region1.second >= region2.second) ? reg_placement::contains : reg_placement::ends_in;
    else if (region1.second > region2.second)
        return reg_placement::begins_in;
    else
        return reg_placement::lies_in;
}

inline void print_all_paths(graph const & grph)
{
    auto recu = [&grph] (auto recu,
                         std::vector<bio::alphabet::dna5> seq,
                         std::string path,
                         size_t i,
                         orientation const o) -> void
    {

        path += "->";
        path += grph.names[i];

        if (o == orientation::plus)
            assign_append(grph.seqs[i], seq);
        else
            assign_append(grph.seqs[i] | std::views::reverse | bio::views::complement, seq);

        if (grph.arcs[i].empty()) // leaf node
        {

            fmt::print("Leaf node\n  path: {}\n  seq: {}\n", path, seq);
        }
        else
        {

            for (arc const a : grph.arcs[i])
            {
                recu(recu, seq, path, a.target_i, a.orient_target);
            }
        }
    };

    recu(recu, {}, {}, 0, orientation::plus);
}


inline void print_chunks(graph const & grph, size_t const w)
{
    auto recu = [&grph, &w] (auto recu,
                         std::vector<bio::alphabet::dna5> seq,
                         std::string path,
                         size_t i,
                         orientation const o,
                         size_t stable_offset,
                         size_t local_offset) -> void
    {
        path += "->";
        path += grph.names[i];

        int64_t needed = (int64_t)w - seq.size();
        int64_t have = (int64_t)seq.size() + grph.seqs[i].size() - local_offset;

        if (have > needed) // we are entering next chunk
        {
            std::span in  = grph.seqs[i] | std::views::drop(local_offset) | std::views::take(needed);
            std::span out = grph.seqs[i] | std::views::drop(local_offset + needed);

            if (o == orientation::plus)
                fmt::print("Chunk\n  path: {}\n  seq: {}|{}\n", path, seq, in);
            else
                fmt::print("Chunk\n  path: {}\n  seq: {}|{}\n", path, seq, in | std::views::reverse | bio::views::complement);

            seq.assign(out.begin(), out.end());
            recu(recu, std::move(seq), std::move(path), i, o, stable_offset + w, local_offset + w);
            return;
        }
        else if (grph.arcs[i].empty()) // leaf node
        {
            std::span in  = grph.seqs[i] | std::views::drop(local_offset);

            if (o == orientation::plus)
                fmt::print("Leaf Chunk\n  path: {}\n  seq: {}|{}\n", path, seq, in);
            else
                fmt::print("Leaf Chunk\n  path: {}\n  seq: {}|{}\n", path, seq, in | std::views::reverse | bio::views::complement);
        }
        else
        {
            std::span in  = grph.seqs[i] | std::views::drop(local_offset);

            if (o == orientation::plus)
                assign_append(in, seq);
            else
                assign_append(in | std::views::reverse | bio::views::complement, seq);

            for (arc const a : grph.arcs[i])
            {
                recu(recu, seq, path, a.target_i, a.orient_target, stable_offset, 0ull);
            }
        }
    };

    recu(recu, {}, {}, 0, orientation::plus, 0, 0);
}


inline void discretise(graph & grph, size_t const w)
{
    // size_t max_offset = *std::ranges::max_element(grph.stable_offsets);
    //
    // graph new_graph;

    /* 1st iteration
     *
     * All segments that span a window border are split into two.
     */

//     for (size_t i = 0; i < grph.names.size(); ++i)
//     {
//         size_t stable_b = grph.stable_offsets[i];
//         size_t stable_e = grph.stable_offsets[i] + grph.seqs[i].size();
//
//         if (stable_e / w > stable_b / w) // end is in different window than beginning
//         {
//             /* create new node */
//             // these are just copied
//             grph.names.push_back(grph.names[i]);
//             grph.stable_name_i.push_back(grph.stable_name_i[i]);
//             grph.ranks.push_back(grph.ranks[i]);
//
//             // these have modified values
//             size_t this_len = w - (stable_b % w);
//             grph.seqs.push_back(grph.seqs[i] | std::views::drop(this_len));
//             grph.stable_offsets.push_back(stable_b + this_len); //TODO can we do this?
//             grph.split_counter.push_back(grph.split_counter[i] + 1);
//
//             // rewire arcs
//             grph.arcs.push_back(std::move(grph.arcs[i])); // move current arcs to new node
//
//
//             /* update old node */
//
//
//
//
//             grph.arcs[i].clear();
//             grph.arcs[i].push_back({.target_i = grph.arcs.size()}); // add single arc from old node to new node
//
//
//         }
//     }



    /* 2nd iteration
     *
     * All segments that begin and end within a window are fully prepended to their successors and
     * then dissolved. This involves
     * and prepended to their successors in the graph. The arcs do not change.
     */


     /**
     *
     * As a result all segements except the left-most one begin on a window boundary

    */

    // for (size_t i = 0; i < max_offset; i += w)
    // {
    //     for (size_t j = 0; j < grph.names.size(); ++j)
    //     {
    //         size_t stable_b = grph.stable_offsets[j];
    //         size_t stable_e = grph.stable_offsets[j] + grph.seqs[j].size();
    //
    //         reg_placement place = overlap({stable_b, stable_e}, {i, i + w});
    //
    //         switch (place)
    //         {
    //             case reg_placement::ends_in:
    //             {
    //                 size_t in_reach = stable_e - (stable_e % w);
    //
    //
    //
    //                 // TODO
    //                 break;
    //             }
    //             case reg_placement::lies_in:
    //                 // TODO
    //                 /* pad the sequence */
    //                 new_graph.seqs.push_back();
    //
    //
    //
    //
    //                 break;
    //             case reg_placement::begins_in:
    //                 // TODO
    //                 break;
    //             case reg_placement::ends_before:
    //             case reg_placement::begins_after:
    //                 // NOTHING NEEDS TO BE DONE
    //                 break;
    //         }
    //
    //
    //
    //
    //     }
    //
    //
    //
    //
    // }
    //
    //


}

