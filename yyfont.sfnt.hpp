#ifndef __yy_hpp_std_font_sfnt__
#define __yy_hpp_std_font_sfnt__

#include "yyfont.hpp"

#include <bit>

namespace Font
{
    namespace SFNT
    {
        enum class ScalerType : uint32_t
        {
            None = 0,
            TrueType = 0x00010000,
            TrueType_OSX_iOS = Utils::Tag32("true"),
            OpenType = Utils::Tag32("OTTO"),
            OldPostScript = Utils::Tag32("typ1")
        };

        struct Header
        {
            ScalerType scaler_type;
            uint16_t num_tables;
            uint16_t search_range;
            uint16_t entry_selector;
            uint16_t range_shift;
        };

        struct DirectoryEntry
        {
            uint32_t tag;
            uint32_t checksum;
            uint32_t offset;
            uint32_t length;
        };

        enum class TableTag : uint32_t
        {
            GSUB = Utils::Tag32("GSUB"),
            glyf = Utils::Tag32("glyf"),
            cmap = Utils::Tag32("cmap"),
            head = Utils::Tag32("head"),
            hhea = Utils::Tag32("hhea"),
            hmtx = Utils::Tag32("hmtx"),
            maxp = Utils::Tag32("maxp"),
            loca = Utils::Tag32("loca")
        };

        struct Table
        {
            virtual ~Table() = default;
            DirectoryEntry directory;
        };

        struct HeadTable : public Table
        {
            int16_t version;
            int16_t revision;
            uint32_t checksum_adjustment;
            uint32_t magic_number;
            uint16_t flags;
            uint16_t units_per_em;
            uint64_t created;
            uint64_t modified;
            int16_t x_min;
            int16_t y_min;
            int16_t x_max;
            int16_t y_max;
            uint16_t mac_style;
            uint16_t lowest_rec_PPEM;
            int16_t font_direction_hint;
            int16_t index_to_loc_format;
            int16_t glyph_data_format;

            HeadTable(const void *mem, DirectoryEntry dir)
            {
                size_t index = dir.offset + 4;
                uint8_t *data = (uint8_t *)mem;

                version = Utils::read16BE(data, index, true);
                revision = Utils::read16BE(data, index, true);
                checksum_adjustment = Utils::read32BE(data, index, true);
                magic_number = Utils::read32BE(data, index, true); /// 1594834165 o 0x5F0F3CF5
                flags = Utils::read16BE(data, index, true);
                units_per_em = Utils::read16BE(data, index, true);
                created = Utils::read64BE(data, index, true);
                modified = Utils::read64BE(data, index, true);
                x_min = Utils::read16BE(data, index, true);
                y_min = Utils::read16BE(data, index, true);
                x_max = Utils::read16BE(data, index, true);
                y_max = Utils::read16BE(data, index, true);
                mac_style = Utils::read16BE(data, index, true);
                lowest_rec_PPEM = Utils::read16BE(data, index, true);
                font_direction_hint = Utils::read16BE(data, index, true);
                index_to_loc_format = Utils::read16BE(data, index, true);
                glyph_data_format = Utils::read16BE(data, index, true);
            }
        };

        struct HheaTable : public Table
        {
            int16_t version;
            int16_t ascent;
            int16_t descent;
            int16_t lineGap;
            uint16_t advance_width_max;
            int16_t min_left_side_bearing;
            int16_t min_right_side_bearing;
            int16_t x_max_extent;
            int16_t caret_slope_rise;
            int16_t caret_slope_run;
            int16_t caret_offset;
            int16_t reserved_0;
            int16_t reserved_1;
            int16_t reserved_2;
            int16_t reserved_3;
            int16_t metric_data_format;
            uint16_t num_of_long_hor_metrics;

            HheaTable(const void *mem, DirectoryEntry dir)
            {
                size_t index = dir.offset + 2;
                uint8_t *data = (uint8_t *)mem;

                version = Utils::read16BE(data, index, true);
                ascent = Utils::read16BE(data, index, true);
                descent = Utils::read16BE(data, index, true);
                lineGap = Utils::read16BE(data, index, true);
                advance_width_max = Utils::read16BE(data, index, true);
                min_left_side_bearing = Utils::read16BE(data, index, true);
                min_right_side_bearing = Utils::read16BE(data, index, true);
                x_max_extent = Utils::read16BE(data, index, true);
                caret_slope_rise = Utils::read16BE(data, index, true);
                caret_slope_run = Utils::read16BE(data, index, true);
                caret_offset = Utils::read16BE(data, index, true);
                reserved_0 = Utils::read16BE(data, index, true);
                reserved_1 = Utils::read16BE(data, index, true);
                reserved_2 = Utils::read16BE(data, index, true);
                reserved_3 = Utils::read16BE(data, index, true);
                metric_data_format = Utils::read16BE(data, index, true);
                num_of_long_hor_metrics = Utils::read16BE(data, index, true);
            }
        };

        struct MaxpTable : public Table
        {
            int16_t version;
            uint16_t num_glyphs;
            uint16_t max_points;
            uint16_t max_contours;
            uint16_t max_component_points;
            uint16_t max_component_contours;
            uint16_t max_zones;
            uint16_t max_twilight_points;
            uint16_t max_storage;
            uint16_t max_function_defs;
            uint16_t max_instruction_defs;
            uint16_t max_stack_elements;
            uint16_t max_size_of_instructions;
            uint16_t max_component_elements;
            uint16_t max_component_depth;

            MaxpTable(const void *mem, DirectoryEntry dir)
            {
                size_t index = dir.offset + 2;
                uint8_t *data = (uint8_t *)mem;

                version = Utils::read16BE(data, index, true);
                num_glyphs = Utils::read16BE(data, index, true);
                max_points = Utils::read16BE(data, index, true);
                max_contours = Utils::read16BE(data, index, true);
                max_component_points = Utils::read16BE(data, index, true);
                max_component_contours = Utils::read16BE(data, index, true);
                max_zones = Utils::read16BE(data, index, true);
                max_twilight_points = Utils::read16BE(data, index, true);
                max_storage = Utils::read16BE(data, index, true);
                max_function_defs = Utils::read16BE(data, index, true);
                max_instruction_defs = Utils::read16BE(data, index, true);
                max_stack_elements = Utils::read16BE(data, index, true);
                max_size_of_instructions = Utils::read16BE(data, index, true);
                max_component_elements = Utils::read16BE(data, index, true);
                max_component_depth = Utils::read16BE(data, index, true);
            }
        };

        struct HmtxTable : public Table
        {
            struct LongHorMetric
            {
                uint16_t advance_width;
                int16_t left_side_bearing;
            };

            std::vector<LongHorMetric> h_metrics;
            std::vector<int16_t> left_side_bearing;

            LongHorMetric getMetrics(GlyphIndex glyph_index) const
            {
                // Normal fonts
                if (glyph_index < h_metrics.size())
                {
                    return h_metrics.at(glyph_index);
                }
                else // ;omospaced fonts
                {
                    LongHorMetric metric;
                    metric.advance_width = h_metrics.at(h_metrics.size() - 1).advance_width;
                    metric.left_side_bearing = left_side_bearing.at(glyph_index - h_metrics.size());
                    return metric;
                }
            }

            HmtxTable(const void *mem, DirectoryEntry dir, uint16_t num_of_long_hor_metrics, uint16_t num_glyphs)
            {
                size_t index = dir.offset;
                uint8_t *data = (uint8_t *)mem;

                num_of_long_hor_metrics = std::max(uint16_t(1), num_of_long_hor_metrics);

                h_metrics.resize(num_of_long_hor_metrics);
                left_side_bearing.resize(num_glyphs - num_of_long_hor_metrics);

                for (size_t i = 0; i < num_of_long_hor_metrics; i++)
                {
                    h_metrics[i].advance_width = Utils::read16BE(data, index, true);
                    h_metrics[i].left_side_bearing = Utils::read16BE(data, index, true);
                }

                for (size_t i = 0; i < num_glyphs - num_of_long_hor_metrics; i++)
                {
                    left_side_bearing[i] = Utils::read16BE(data, index, true);
                }
            }
        };

        struct CmapTable : public Table
        {
            struct EncodingSubtable
            {
                uint16_t platform_id;
                uint16_t platform_specific_id;
                uint32_t offset;
            };

            uint16_t version;
            uint16_t number_of_subtables;
            std::vector<EncodingSubtable> encoding_subtables;
            std::unique_ptr<uint8_t[]> mapping_data;

            CmapTable(const void *mem, DirectoryEntry dir)
            {
                directory = dir;
                size_t index = directory.offset;
                const uint8_t *data = (uint8_t *)mem;

                version = Utils::read16BE(data, index, true);
                number_of_subtables = Utils::read16BE(data, index, true);

                encoding_subtables.clear();

                // EncodingSubtable is aligned and does not produce padding
                uint64_t data_offset = index + sizeof(EncodingSubtable) * number_of_subtables - directory.offset;

                for (size_t i = 0; i < number_of_subtables; i++)
                {
                    CmapTable::EncodingSubtable encs;

                    encs.platform_id = Utils::read16BE(data, index, true);
                    encs.platform_specific_id = Utils::read16BE(data, index, true);
                    encs.offset = Utils::read32BE(data, index, true) - data_offset;

                    encoding_subtables.push_back(encs);
                }

                uint64_t data_size = directory.length - data_offset;

                mapping_data = std::make_unique<uint8_t[]>(data_size);
                memcpy(mapping_data.get(), &((uint8_t *)(mem))[directory.offset + data_offset], data_size);
            }

            GlyphIndex map(CodePoint codepoint)
            {
                EncodingSubtable *chosen_table = nullptr;

                // select a subtable
                for (auto &encsub : encoding_subtables)
                {
                    /*
                    Unicode encoding subtables are used in preference to non-Unicode encoding subtables.
                    Unicode encoding subtables not restricted to the BMP are used in preference to subtables restricted to the BMP.
                    Unicode variation sequence subtables are always processed if another Unicode cmap of type 4 or 12 is present.
                    */
                    if (
                        !chosen_table ||                                                                                // Select any if present
                        (chosen_table->platform_id != 0 && encsub.platform_id == 0) ||                                  // Prefer unicode
                        (chosen_table->platform_id == 0 && encsub.platform_id == 0 && encsub.platform_specific_id == 4) // If don't BMP better
                    )
                    {
                        chosen_table = &encsub;
                    }
                }

                if (!chosen_table)
                {
                    return 0;
                }

                // I put my font's data_size as a uint8, seeing the format change to 45080 was so XD
                size_t index = chosen_table->offset;
                const uint8_t *data = mapping_data.get();
                uint16_t format = Utils::read16BE(data, index, true);

                // Ahh ugly switch sytax
                switch (format)
                {
                case 4:
                {
                    uint16_t length = Utils::read16BE(data, index, true);
                    uint16_t language = Utils::read16BE(data, index, true);
                    // Segments couTN
                    uint16_t twice_seg_count = Utils::read16BE(data, index, true);
                    uint16_t seg_count = twice_seg_count >> 1; // (/2)
                    // Values for binary search
                    uint16_t search_range = Utils::read16BE(data, index, true);
                    uint16_t entry_selector = Utils::read16BE(data, index, true);
                    uint16_t range_shift = Utils::read16BE(data, index, true);
                    // get end codes array
                    const uint16_t *end_codes = (uint16_t *)(data + index);
                    index += twice_seg_count; // (+= uint16 segc)
                    // this should be zero >:(
                    uint16_t reserved_pad = Utils::read16BE(data, index, true);
                    if (reserved_pad != 0) // I did this wrong? I didn't
                    {
                        return 0;
                    }
                    // get start codes array
                    const uint16_t *start_codes = (uint16_t *)(data + index);
                    index += twice_seg_count; // (+= uint16 segc)
                    // get id deltas
                    const uint16_t *id_deltas = (uint16_t *)(data + index);
                    index += twice_seg_count; // (+= uint16 segc)
                    // get id range offsets
                    const uint16_t *id_range_offsets = (uint16_t *)(data + index);
                    index += twice_seg_count; // (+= uint16 segc)

                    for (uint16_t seg = 0; seg < seg_count; seg++)
                    {
                        // Character falls into range
                        if (Utils::swap16(end_codes[seg]) >= uint16_t(codepoint))
                        {
                            if (Utils::swap16(start_codes[seg]) <= uint16_t(codepoint))
                            {
                                if (id_range_offsets[seg] != 0) // BE or LE, if it's zero it's zero
                                {
                                    uint16_t gi_index = id_range_offsets[seg + (uint16_t(codepoint) - Utils::swap16(start_codes[seg])) + Utils::swap16(id_range_offsets[seg]) / 2];
                                    if (gi_index == 0) // BE or LE, if it's zero it's zero
                                    {
                                        return 0;
                                    }
                                    else
                                    {
                                        return uint16_t(Utils::swap16(gi_index) + Utils::swap16(id_deltas[seg]));
                                    }
                                }
                                else
                                {
                                    return uint16_t(Utils::swap16(id_deltas[seg]) + uint16_t(codepoint));
                                }
                            }
                            else
                            {
                                return 0;
                            }
                        }
                    }
                }
                break;
                case 12:
                {
                    Utils::read16BE(data, index, true); // Read WTF reserved field
                    uint32_t length = Utils::read32BE(data, index, true);
                    uint32_t language = Utils::read32BE(data, index, true);
                    uint32_t n_groups = Utils::read32BE(data, index, true);

                    // There should be a catch, it's impossible for it to be that easy :O
                    // IDK if i should binary search.
                    for (uint32_t group = 0; group < n_groups; group++)
                    {
                        uint32_t start_ch_code = Utils::read32BE(data, index, true);
                        uint32_t end_ch_code = Utils::read32BE(data, index, true);
                        uint32_t start_gi_code = Utils::read32BE(data, index, true);
                        if (codepoint >= start_ch_code && codepoint <= end_ch_code)
                        {
                            return codepoint - start_ch_code + start_gi_code;
                        }
                    }
                }
                break;
                default:
                    break;
                }

                return 0;
            }
        };

        struct GlyfTable : public Table
        {
            union OutlineFlag
            {
                struct
                {
                    uint8_t on_curve : 1;

                    uint8_t x_short : 1;
                    uint8_t y_short : 1;

                    uint8_t repeat : 1;

                    uint8_t x_short_pos : 1;
                    uint8_t y_short_pos : 1;

                    uint8_t reserved1 : 1;
                    uint8_t reserved2 : 1;
                };
                uint8_t flag;
            };

            struct ComponentDescription
            {
                union Flags
                {
                    struct
                    {
                        uint8_t arg_1_and_2_are_words : 1; // 0x1
                        uint8_t args_are_xy_values : 1;    // 0x2
                        uint8_t round_xy_to_grid : 1;      // 0x4
                        uint8_t we_have_a_scale : 1;       // 0x8

                        uint8_t obsolete_0 : 1;               // 0x10
                        uint8_t more_components : 1;          // 0x20
                        uint8_t we_have_an_x_and_y_scale : 1; // 0x40
                        uint8_t we_have_a_two_by_two : 1;     // 0x80

                        uint8_t we_have_instructions : 1; // 0x100
                        uint8_t use_my_metrics : 1;       // 0x200
                        uint8_t overlap_compound : 1;     // 0x400
                        uint8_t reserved_0 : 1;           // 0x800

                        uint8_t reserved_1 : 1; // 0x1000
                        uint8_t reserved_2 : 1; // 0x2000
                        uint8_t reserved_3 : 1; // 0x4000
                        uint8_t reserved_4 : 1; // 0x8000
                    };
                    uint16_t flag;
                } flags;
                uint16_t glyph_index;
                int16_t argument_1;
                int16_t argument_2;
                struct TransformOption
                {
                    int16_t xy_scale;
                    int16_t x_scale;
                    int16_t scale_01;
                    int16_t scale_10;
                    int16_t y_scale;
                } transform_option;
            };

            std::vector<int64_t> index_to_loc;
            std::unique_ptr<uint8_t[]> glyph_data;

            GlyfTable(const void *mem, DirectoryEntry loca_dir, DirectoryEntry glyf_dir, uint16_t index_to_loc_format, uint16_t num_glyphs)
            {
                directory = glyf_dir;
                const uint8_t *data = (uint8_t *)mem;

                size_t index = loca_dir.offset;

                index_to_loc.resize(num_glyphs + 1); // There's an extra offset, why? IDK

                for (int i = 0; i < num_glyphs + 1; i++)
                {
                    if (index_to_loc_format == 1)
                    {
                        index_to_loc[i] = Utils::read32BE(data, index, true);
                    }
                    else if (index_to_loc_format == 0)
                    {
                        index_to_loc[i] = Utils::read16BE(data, index, true) * 2;
                    }
                }

                index = glyf_dir.offset;

                glyph_data = std::make_unique<uint8_t[]>(directory.length);
                memcpy(glyph_data.get(), &((uint8_t *)(mem))[directory.offset], directory.length);
            }
        };

        class SFNTFont : public Font
        {
        protected:
            bool loadHeader(const void *mem, size_t size)
            {
                if (mem == nullptr || size == 0)
                {
                    return false;
                }
                size_t index = 0;
                const uint8_t *data = (uint8_t *)mem;

                file_header.scaler_type = (ScalerType)Utils::read32BE(data, index, true);
                file_header.num_tables = Utils::read16BE(data, index, true);
                file_header.search_range = Utils::read16BE(data, index, true);
                file_header.entry_selector = Utils::read16BE(data, index, true);
                file_header.range_shift = Utils::read16BE(data, index, true);

                return true;
            }

            bool loadFontDirectory(const void *mem, size_t size)
            {
                size_t index = sizeof(Header);
                const uint8_t *data = (uint8_t *)mem;

                for (size_t i = 0; i < file_header.num_tables; i++)
                {
                    DirectoryEntry dir;
                    dir.tag = Utils::read32BE(data, index, true);
                    dir.checksum = Utils::read32BE(data, index, true);
                    dir.offset = Utils::read32BE(data, index, true);
                    dir.length = Utils::read32BE(data, index, true);

                    directories[TableTag(dir.tag)] = dir;
                }

                return true;
            }

            bool loadRequired(const void *mem, size_t size)
            {
                auto head_it = directories.find(TableTag::head);
                auto maxp_it = directories.find(TableTag::maxp);
                auto hhea_it = directories.find(TableTag::hhea);
                auto hmtx_it = directories.find(TableTag::hmtx);
                auto cmap_it = directories.find(TableTag::cmap);

                // Please account for bitmap fonts later
                if (
                    head_it != directories.end() &&
                    maxp_it != directories.end() &&
                    hhea_it != directories.end() &&
                    hmtx_it != directories.end() &&
                    cmap_it != directories.end())
                {
                    DirectoryEntry &head_entry = head_it->second;
                    DirectoryEntry &maxp_entry = maxp_it->second;
                    DirectoryEntry &hhea_entry = hhea_it->second;
                    DirectoryEntry &hmtx_entry = hmtx_it->second;
                    DirectoryEntry &cmap_entry = cmap_it->second;

                    tables[TableTag::head] = std::make_unique<HeadTable>(mem, head_entry);
                    tables[TableTag::cmap] = std::make_unique<CmapTable>(mem, cmap_entry);

                    auto maxp_table = (MaxpTable *)(tables[TableTag::maxp] = std::make_unique<MaxpTable>(mem, maxp_entry)).get();
                    auto hhea_tale = (HheaTable *)(tables[TableTag::hhea] = std::make_unique<HheaTable>(mem, hhea_entry)).get();

                    tables[TableTag::hmtx] = std::make_unique<HmtxTable>(mem, hmtx_entry, hhea_tale->num_of_long_hor_metrics, maxp_table->num_glyphs);

                    return true;
                }
                else
                {
                    return false;
                }
            }

        public:
            virtual ~SFNTFont() = default;

            Header file_header;
            std::map<TableTag, DirectoryEntry> directories;
            std::map<TableTag, std::unique_ptr<Table>> tables;

            GlyphIndex getId(CodePoint codepoint) override
            {
                const auto it = tables.find(TableTag::cmap);
                if (it != tables.end() && dynamic_cast<CmapTable *>(it->second.get()))
                {
                    return ((CmapTable *)(it->second.get()))->map(codepoint);
                }
                return 0;
            }
        };

        class TrueTypeFont : public SFNTFont
        {
        public:
            ~TrueTypeFont() override {}

            GlyphMetrics getMetrics(GlyphIndex g_index) override
            {
                GlyphMetrics metrics;

                /* At this point, any valid ttf ( NOT BITMAP ) has those tables */
                GlyfTable *glyf_table = dynamic_cast<GlyfTable *>(tables[TableTag::glyf].get());
                HmtxTable *hmtx_table = dynamic_cast<HmtxTable *>(tables[TableTag::hmtx].get());

                size_t index = glyf_table->index_to_loc[g_index];
                const uint8_t *data = glyf_table->glyph_data.get(); // Safety? fuck off

                Utils::read16BE(data, index, true);

                int16_t xMin = std::bit_cast<int16_t>(Utils::read16BE(data, index, true));
                int16_t yMin = std::bit_cast<int16_t>(Utils::read16BE(data, index, true));
                int16_t xMax = std::bit_cast<int16_t>(Utils::read16BE(data, index, true));
                int16_t yMax = std::bit_cast<int16_t>(Utils::read16BE(data, index, true));

                HmtxTable::LongHorMetric lhm = hmtx_table->getMetrics(g_index);

                metrics.bounds.x = xMin;
                metrics.bounds.y = yMin;
                metrics.bounds.w = xMax - xMin;
                metrics.bounds.h = yMax - yMin;

                metrics.advance_width = lhm.advance_width;
                metrics.left_side_bearing = lhm.left_side_bearing;
                metrics.right_side_bearing = lhm.advance_width - lhm.left_side_bearing - metrics.bounds.w;

                return metrics;
            }

            GlyphCurve getCurves(GlyphIndex g_index) override
            {
                /* At this point, any valid ttf ( NOT BITMAP ) has those tables */
                GlyfTable *glyf_table = dynamic_cast<GlyfTable *>(tables[TableTag::glyf].get());

                GlyphCurve glyph_curve = GlyphCurve();

                size_t index = glyf_table->index_to_loc[g_index];
                const uint8_t *data = glyf_table->glyph_data.get(); // Safety? fuck off

                int16_t numberOfContours = std::bit_cast<int16_t>(Utils::read16BE(data, index, true));

                int16_t xMin = std::bit_cast<int16_t>(Utils::read16BE(data, index, true));
                int16_t yMin = std::bit_cast<int16_t>(Utils::read16BE(data, index, true));
                int16_t xMax = std::bit_cast<int16_t>(Utils::read16BE(data, index, true));
                int16_t yMax = std::bit_cast<int16_t>(Utils::read16BE(data, index, true));

                if (numberOfContours == 0)
                {
                    return glyph_curve;
                }
                else if (numberOfContours > 0)
                {
                    std::vector<std::pair<uint64_t, uint64_t>> contour_limits;

                    //contour_limits
                    for (int i = 0; i < numberOfContours; i++)
                    {
                        uint16_t start_point = i == 0 ? 0 : contour_limits.back().second + 1;
                        uint16_t end_point = Utils::read16BE(data, index, true);

                        contour_limits.push_back({start_point, end_point});
                    }

                    size_t point_count = contour_limits.at(numberOfContours - 1).second + 1;

                    uint16_t instructionsLength = Utils::read16BE(data, index, true);
                    index += instructionsLength; // For now, skip instructions

                    std::vector<GlyfTable::OutlineFlag> flags;
                    flags.resize(point_count);

                    for (size_t i = 0; i < point_count; i++)
                    {
                        flags.at(i).flag = Utils::read8(data, index, true);
                        if (flags.at(i).repeat)
                        {
                            uint8_t repeat_count = Utils::read8(data, index, true);

                            while (repeat_count-- > 0)
                            {
                                i++;
                                flags.at(i).flag = flags.at(i - 1).flag;
                            }
                        }
                    }

                    std::vector<Math::Vector2i> points;
                    points.resize(point_count);

                    {
                        int64_t previous_coordinate = 0;
                        int64_t current_coordinate = 0;

                        for (size_t i = 0; i < point_count; i++)
                        {
                            int coordinate_flag = (flags.at(i).x_short << 1) | (flags.at(i).x_short_pos);
                            switch (coordinate_flag)
                            {
                            case 0:
                            {
                                uint16_t cc = Utils::read16BE(data, index, true);
                                current_coordinate = *reinterpret_cast<int16_t *>(&cc);
                            }
                            break;
                            case 1:
                            {
                                current_coordinate = 0;
                            }
                            break;
                            case 2:
                            {
                                current_coordinate = -1 * (int16_t)Utils::read8(data, index, true);
                            }
                            break;
                            case 3:
                            {
                                current_coordinate = Utils::read8(data, index, true);
                            }
                            break;
                            }

                            points.at(i).x = current_coordinate + previous_coordinate;
                            previous_coordinate = points.at(i).x;
                        }
                    }

                    {
                        int64_t previous_coordinate = 0;
                        int64_t current_coordinate = 0;

                        for (size_t i = 0; i < point_count; i++)
                        {
                            uint64_t coordinate_flag = (flags.at(i).y_short << 1) | (flags.at(i).y_short_pos);
                            switch (coordinate_flag)
                            {
                            case 0:
                            {
                                uint16_t cc = Utils::read16BE(data, index, true);
                                current_coordinate = *reinterpret_cast<int16_t *>(&cc);
                            }
                            break;
                            case 1:
                            {
                                current_coordinate = 0;
                            }
                            break;
                            case 2:
                            {
                                current_coordinate = -1 * (int16_t)Utils::read8(data, index, true);
                            }
                            break;
                            case 3:
                            {
                                current_coordinate = Utils::read8(data, index, true);
                            }
                            break;
                            }

                            points.at(i).y = current_coordinate + previous_coordinate;
                            previous_coordinate = points.at(i).y;
                        }
                    }

                    size_t first_point_on_curve_contour = 0;
                    size_t last_selected_contour = 0;

                    for (int point_index = 0; point_index < point_count; point_index++)
                    {
                        // If we reached a new contour, move the selection forward
                        if (point_index > contour_limits[last_selected_contour].second)
                        {
                            glyph_curve.contour_limits.push_back(
                                {first_point_on_curve_contour, glyph_curve.curves.size() - 1}
                            );
                            first_point_on_curve_contour = glyph_curve.curves.size();
                            last_selected_contour++;
                        }
                        // Select the appropriate contour limit indexes
                        auto limits = contour_limits[last_selected_contour];
                        // Get prev & next points indexes with contour wraparound
                        size_t prev_point_index = point_index == limits.first ? limits.second : point_index - 1;
                        size_t next_point_index = point_index == limits.second ? limits.first : point_index + 1;

                        Math::Vector2i prev_point = points[prev_point_index];
                        Math::Vector2i curr_point = points[point_index];
                        Math::Vector2i next_point = points[next_point_index];

                        if (!flags[point_index].on_curve) // Is this point bezier?
                        {
                            if (!flags[prev_point_index].on_curve)
                            {
                                prev_point = Math::Vector2i(
                                    (curr_point.x + prev_point.x) / 2,
                                    (curr_point.y + prev_point.y) / 2);
                            }

                            if (!flags[next_point_index].on_curve)
                            {
                                next_point = Math::Vector2i(
                                    (curr_point.x + next_point.x) / 2,
                                    (curr_point.y + next_point.y) / 2);
                            }

                            GlyphCurve::BezierCurve point;
                            point.center = curr_point;
                            point.tail = prev_point;
                            point.head = next_point;

                            glyph_curve.curves.push_back(point);
                        }
                        // This point is a true corner and thus garbage for this case.
                        else if (!flags[next_point_index].on_curve)
                        {
                            continue;
                        }
                        // This point is part of a collinear pair
                        else if (flags[next_point_index].on_curve)
                        {
                            GlyphCurve::BezierCurve point;

                            point.is_collinear = true; // Enable any further optimizations & give me less headaches

                            point.tail = curr_point;
                            point.head = next_point;
                            point.center = Math::Vector2i(
                                (curr_point.x + next_point.x) / 2,
                                (curr_point.y + next_point.y) / 2);

                            glyph_curve.curves.push_back(point); // And DONT skip points! :D, trust me
                        }
                    }

                    if (point_count > 0)
                    {
                        glyph_curve.contour_limits.push_back(
                                {first_point_on_curve_contour, glyph_curve.curves.size() - 1}
                        );
                    }
                }
                else
                {
                    numberOfContours = 0;

                    bool one_more = true;

                    do
                    {
                        double matrix[6] = {1, 0, 0, 1, 0, 0};

                        GlyfTable::ComponentDescription description;

                        description.flags.flag = Utils::read16BE(data, index, true);
                        description.glyph_index = Utils::read16BE(data, index, true);

                        one_more = description.flags.more_components;

                        if (description.flags.arg_1_and_2_are_words)
                        {
                            uint16_t argument_1 = Utils::read16BE(data, index, true);
                            uint16_t argument_2 = Utils::read16BE(data, index, true);
                            description.argument_1 = *(int16_t *)&argument_1;
                            description.argument_2 = *(int16_t *)&argument_2;
                        }
                        else
                        {
                            uint8_t argument_1 = Utils::read8(data, index, true);
                            uint8_t argument_2 = Utils::read8(data, index, true);
                            description.argument_1 = *(int8_t *)&argument_1;
                            description.argument_2 = *(int8_t *)&argument_2;
                        }

                        if (description.flags.args_are_xy_values)
                        {
                            matrix[4] = description.argument_1;
                            matrix[5] = description.argument_2;
                        }
                        else
                        {
                            std::printf("Warning from trutype's glyph loader: Please implement matching points compound glyph positioning.\n");
                        }

                        if (description.flags.we_have_a_scale)
                        {
                            description.transform_option.xy_scale = Utils::read16BE(data, index, true);

                            matrix[0] = matrix[3] = description.transform_option.xy_scale / 16384.0;
                            matrix[1] = matrix[2] = 0;
                        }
                        else if (description.flags.we_have_an_x_and_y_scale)
                        {
                            description.transform_option.x_scale = Utils::read16BE(data, index, true);
                            description.transform_option.y_scale = Utils::read16BE(data, index, true);

                            matrix[0] = description.transform_option.x_scale / 16384.0;
                            matrix[3] = description.transform_option.y_scale / 16384.0;
                            matrix[1] = matrix[2] = 0;
                        }
                        else if (description.flags.we_have_a_two_by_two)
                        {
                            description.transform_option.x_scale = Utils::read16BE(data, index, true);
                            description.transform_option.scale_01 = Utils::read16BE(data, index, true);
                            description.transform_option.scale_10 = Utils::read16BE(data, index, true);
                            description.transform_option.y_scale = Utils::read16BE(data, index, true);

                            matrix[0] = description.transform_option.x_scale / 16384.0;
                            matrix[1] = description.transform_option.scale_01 / 16384.0;
                            matrix[2] = description.transform_option.scale_10 / 16384.0;
                            matrix[3] = description.transform_option.y_scale / 16384.0;
                        }

                        float transformM = std::sqrt(matrix[0] * matrix[0] + matrix[1] * matrix[1]);
                        float transformN = std::sqrt(matrix[2] * matrix[2] + matrix[3] * matrix[3]);

                        GlyphCurve subglyph = getCurves(description.glyph_index);

                        for (GlyphCurve::BezierCurve &point : subglyph.curves)
                        {
                            point.tail.x = (transformM * (matrix[0] * point.tail.x + matrix[2] * point.tail.y + matrix[4]));
                            point.tail.y = (transformN * (matrix[1] * point.tail.x + matrix[3] * point.tail.y + matrix[5]));

                            point.center.x = (transformM * (matrix[0] * point.center.x + matrix[2] * point.center.y + matrix[4]));
                            point.center.y = (transformN * (matrix[1] * point.center.x + matrix[3] * point.center.y + matrix[5]));

                            point.head.x = (transformM * (matrix[0] * point.head.x + matrix[2] * point.head.y + matrix[4]));
                            point.head.y = (transformN * (matrix[1] * point.head.x + matrix[3] * point.head.y + matrix[5]));
                        }

                        glyph_curve.merge( subglyph );
                    } while (one_more);
                }

                return glyph_curve;
            }

            static std::shared_ptr<Font> Loader(const void *mem, size_t size)
            {
                auto font = std::make_unique<TrueTypeFont>();

                if (
                    !font->loadHeader(mem, size) ||
                    !(font->file_header.scaler_type == ScalerType::TrueType || font->file_header.scaler_type == ScalerType::TrueType_OSX_iOS) ||
                    !font->loadFontDirectory(mem, size) ||
                    !font->loadRequired(mem, size))
                {
                    return nullptr;
                }

                auto loca_it = font->directories.find(TableTag::loca);
                auto glyf_it = font->directories.find(TableTag::glyf);

                if (
                    loca_it != font->directories.end() &&
                    glyf_it != font->directories.end())
                {
                    DirectoryEntry &loca_entry = loca_it->second;
                    DirectoryEntry &glyf_entry = glyf_it->second;

                    HeadTable *head_table = dynamic_cast<HeadTable *>(font->tables[TableTag::head].get());
                    MaxpTable *maxp_table = dynamic_cast<MaxpTable *>(font->tables[TableTag::maxp].get());

                    font->tables[TableTag::glyf] = std::make_unique<GlyfTable>(mem, loca_entry, glyf_entry, head_table->index_to_loc_format, maxp_table->num_glyphs);
                }
                else
                {
                    return nullptr;
                }

                return font;
            }
        };

        class OpenTypeFont : public SFNTFont
        {
        public:
            ~OpenTypeFont() override {}

            // What an undefined reference ;)
            GlyphMetrics getMetrics(GlyphIndex index) override
            {
                return GlyphMetrics();
            }

            GlyphCurve getCurves(GlyphIndex index) override
            {
                return GlyphCurve();
            }

            static std::shared_ptr<Font> Loader(const void *mem, size_t size)
            {
                auto font = std::make_unique<OpenTypeFont>();

                if (
                    !font->loadHeader(mem, size) ||
                    font->file_header.scaler_type != ScalerType::OpenType ||
                    !font->loadFontDirectory(mem, size) ||
                    !font->loadRequired(mem, size))
                {
                    return nullptr;
                }

                return font;
            }
        };
    };
}

#endif