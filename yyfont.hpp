#ifndef __yy_hpp_std_font__
#define __yy_hpp_std_font__

#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <filesystem>
#include <optional>
#include <memory>
#include <array>
#include <map>
#include <utility>

namespace Font
{
    namespace Math
    {
        struct Vector2i
        {
            int64_t x = 0;
            int64_t y = 0;

            Vector2i() = default;
            Vector2i(const Vector2i &v) = default;
            Vector2i(int64_t _x, int64_t _y) : x(_x), y(_y) {};

            constexpr Vector2i operator-() const noexcept
            {
                return Vector2i(-x, -y);
            }

            constexpr Vector2i operator-(const Vector2i &rhs) const noexcept
            {
                return Vector2i(x - rhs.x, y - rhs.y);
            }

            constexpr Vector2i operator*(const Vector2i &rhs) const noexcept
            {
                return Vector2i(x * rhs.x, y * rhs.y);
            }

            constexpr Vector2i operator*(const double &rhs) const noexcept
            {
                return Vector2i(x * rhs, y * rhs);
            }

            constexpr Vector2i operator+(const Vector2i &rhs) const noexcept
            {
                return Vector2i(x + rhs.x, y + rhs.y);
            }
        };

        struct Vector2f
        {
            double x = 0;
            double y = 0;

            Vector2f() = default;
            Vector2f(const Vector2f &v) = default;
            Vector2f(const Vector2i &v) : x(v.x), y(v.y) {};
            Vector2f(double v) : x(v), y(v) {}
            Vector2f(double _x, double _y) : x(_x), y(_y) {};

            constexpr Vector2f operator-() const noexcept
            {
                return Vector2f(-x, -y);
            }

            constexpr Vector2f operator-(const Vector2f &rhs) const noexcept
            {
                return Vector2f(x - rhs.x, y - rhs.y);
            }

            constexpr Vector2f operator*(const Vector2f &rhs) const noexcept
            {
                return Vector2f(x * rhs.x, y * rhs.y);
            }

            constexpr Vector2f operator/(const double &rhs) const noexcept
            {
                return Vector2f(x / rhs, y / rhs);
            }

            constexpr Vector2f operator+(const Vector2f &rhs) const noexcept
            {
                return Vector2f(x + rhs.x, y + rhs.y);
            }
        };

        struct Vector3f
        {
            double x = 0;
            double y = 0;
            double z = 0;

            Vector3f() = default;
            Vector3f(const Vector3f &v) = default;
            Vector3f(const double v) : x(v), y(v), z(v) {};
            Vector3f(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};

            constexpr Vector3f operator*(const double &rhs) const noexcept
            {
                return Vector3f(x * rhs, y * rhs, z * rhs);
            }

            constexpr Vector3f operator+(const Vector3f &rhs) const noexcept
            {
                return Vector3f(x + rhs.x, y + rhs.y, z + rhs.z);
            }
        };

        struct Vector3i
        {
            int64_t x = 0;
            int64_t y = 0;
            int64_t z = 0;

            Vector3i() {};
            Vector3i(const Vector3i &v)
            {
                x = v.x;
                y = v.y;
                z = v.z;
            };
            Vector3i(int64_t _x, int64_t _y, int64_t _z) : x(_x), y(_y), z(_z) {};
        };

        struct Vector4f
        {
            double x = 0;
            double y = 0;
            double z = 0;
            double w = 0;

            Vector4f() {};
            Vector4f(const Vector4f &v)
            {
                x = v.x;
                y = v.y;
                z = v.z;
                w = v.w;
            };
            Vector4f(double _x, double _y, double _z, double _w) : x(_x), y(_y), z(_z), w(_w) {};
        };

        struct Vector4i
        {
            int64_t x = 0;
            int64_t y = 0;
            int64_t z = 0;
            int64_t w = 0;

            Vector4i() {};
            Vector4i(const Vector4i &v)
            {
                x = v.x;
                y = v.y;
                z = v.z;
                w = v.w;
            };
            Vector4i(int64_t _x, int64_t _y, int64_t _z, int64_t _w) : x(_x), y(_y), z(_z), w(_w) {};
        };

        struct Rect4i
        {
            int64_t x = 0;
            int64_t y = 0;
            int64_t w = 0;
            int64_t h = 0;

            Rect4i() {};
            Rect4i(int64_t _x, int64_t _y, int64_t _w, int64_t _h) : x(_x), y(_y), w(_w), h(_h) {};
        };

        struct Rect4f
        {
            double x = 0;
            double y = 0;
            double w = 0;
            double h = 0;

            Rect4f() {};
            Rect4f(double _x, double _y, double _w, double _h) : x(_x), y(_y), w(_w), h(_h) {};
        };

        Vector2f abs(const Vector2f &k)
        {
            return Vector2f(
                std::abs(k.x),
                std::abs(k.y));
        }

        Vector2f pow(const Vector2f &k, const double p)
        {
            return Vector2f(
                std::pow(k.x, p),
                std::pow(k.y, p));
        }

        double dot(const Vector2f &k, const Vector2f &p)
        {
            return k.x * p.x + k.y * p.y;
        }

        float sign(float v)
        {
            if (v == 0)
            {
                return 0;
            }
            else if (v < 0)
            {
                return -1;
            }
            return 1;
        }

        Vector2f sign(const Vector2f &k)
        {
            return Vector2f(
                sign(k.x),
                sign(k.y));
        }
    }

    namespace Utils
    {
        /**
        @brief Reads a byte from memory.
        @param data Memory address to read from.
        @param index Reading position.
        @return The read value.
        */
        uint8_t read8(const void *data, size_t index)
        {
            uint8_t val = *((uint8_t *)data + index);

            return val;
        }

        /**
        @brief Reads a Big Endiant uint16 from memory.
        @param data Memory address to read from.
        @param index Reading position.
        @return The read value.
        */
        uint16_t read16BE(const void *data, size_t index)
        {
            uint16_t val = 0x0000 |
                           uint16_t(*((uint8_t *)data + index + 1)) |
                           uint16_t(*((uint8_t *)data + index)) << 8;

            return val;
        }

        /**
        @brief Reads a Big Endiant uint32 from memory.
        @param data Memory address to read from.
        @param index Reading position.
        @return The read value.
        */
        uint32_t read32BE(const void *data, size_t index)
        {
            uint32_t val = uint32_t(0) |
                           uint32_t(*((uint8_t *)data + index + 3)) |
                           uint32_t(*((uint8_t *)data + index + 2)) << 8 |
                           uint32_t(*((uint8_t *)data + index + 1)) << 16 |
                           uint32_t(*((uint8_t *)data + index)) << 24;

            return val;
        }

        /**
        @brief Reads a Big Endiant uint64 from memory.
        @param data Memory address to read from.
        @param index Reading position.
        @return The read value.
        */
        uint64_t read64BE(const void *data, size_t index)
        {
            uint64_t val = uint64_t(0) |
                           uint64_t(*((uint8_t *)data + index + 7)) |
                           uint64_t(*((uint8_t *)data + index + 6)) << 8 |
                           uint64_t(*((uint8_t *)data + index + 5)) << 16 |
                           uint64_t(*((uint8_t *)data + index + 4)) << 24 |
                           uint64_t(*((uint8_t *)data + index + 3)) << 32 |
                           uint64_t(*((uint8_t *)data + index + 2)) << 40 |
                           uint64_t(*((uint8_t *)data + index + 1)) << 48 |
                           uint64_t(*((uint8_t *)data + index)) << 56;

            return val;
        }

        /**
        @brief Reads a byte from memory and move forward.
        @param data Memory address to read from.
        @param index Reading position.
        @return The read value.
        */
        uint8_t read8(const void *data, size_t &index, bool)
        {
            uint8_t val = read8(data, index);
            index++;
            return val;
        }

        /**
        @brief Reads a Big Endiant uint16 from memory and move forward.
        @param data Memory address to read from.
        @param index Reading position.
        @return The read value.
        */
        uint16_t read16BE(const void *data, size_t &index, bool)
        {
            uint32_t val = read16BE(data, index);
            index += 2;
            return val;
        }

        /**
        @brief Reads a Big Endiant uint32 from memory and move forward.
        @param data Memory address to read from.
        @param index Reading position.
        @return The read value.
        */
        uint32_t read32BE(const void *data, size_t &index, bool)
        {
            uint32_t val = read32BE(data, index);
            index += 4;
            return val;
        }

        /**
        @brief Reads a Big Endiant uint64 from memory and move forward.
        @param data Memory address to read from.
        @param index Reading position.
        @return The read value.
        */
        uint64_t read64BE(const void *data, size_t &index, bool)
        {
            uint64_t val = read64BE(data, index);
            index += 8;
            return val;
        }

        /**
        @brief Swaps a given uint16 byte by byte.
        @param value Value to swap.
        @return Swapped value.
        */
        uint16_t swap16(uint16_t value)
        {
            return read16BE(&value, 0);
        }

        /**
        @brief Swaps a given uint32 byte by byte.
        @param value Value to swap.
        @return Swapped value.
        */
        uint32_t swap32(uint32_t value)
        {
            return read32BE(&value, 0);
        }

        /**
        @brief Swaps a given uint64 byte by byte.
        @param value Value to swap.
        @return Swapped value.
        */
        uint64_t swap64(uint64_t value)
        {
            return read64BE(&value, 0);
        }

        /**
        @brief Tries to read the contents of a file into a string.
        @param path filesystem path.
        @return File contents string or none.
        */
        std::optional<std::string> readFile(std::filesystem::path path)
        {
            std::ifstream file(path, std::ios::binary | std::ios::ate);

            if (file.fail() || file.bad())
            {
                return {};
            }

            size_t file_size = file.tellg();
            file.seekg(0, std::ios::beg);

            std::string data;

            if (file_size == 0)
            {
                return "";
            }

            data.resize(file_size);
            file.read((char *)&data.at(0), file_size);
            file.close();

            return data;
        }

        /**
        @brief Creates a uint32 tag from 4 chars.
        */
        constexpr uint32_t Tag32(const char chars[4])
        {
            return uint32_t(chars[0]) << 24 |
                   uint32_t(chars[1]) << 16 |
                   uint32_t(chars[2]) << 8 |
                   uint32_t(chars[3]);
        }

        float Q_sin(float x)
        {
            return std::sin(x);
        }

        float Q_cos(float x)
        {
            return std::cos(x);
        }

        void solveQuadratic(double a, double b, double c, double &r_0, double &r_1)
        {
            if (a == 0)
            {
                r_0 = -c / b;
                r_1 = r_0;
                return;
            }

            double d = b * b - 4 * a * c;

            if (d == 0)
            {
                r_0 = -b / (2 * a);
                r_1 = r_0;
            }
            else if (d > 0)
            {
                double cd = std::sqrt(d);
                r_0 = (-b - cd) / (2 * a);
                r_1 = (-b + cd) / (2 * a);
            }
            else
            {
                r_0 = NAN;
                r_1 = NAN;
            }
        }

        void solveCubic(float c0, float c1, float c2, float c3, float &r_0, float &r_1, float &r_2)
        {
            float a = c1 / c0;
            float b = c2 / c0;
            float c = c3 / c0;

            float p = b - a * a / 3.0;
            float p3 = p * p * p;
            float q = a * (2.0 * a * a - 9.0 * b) / 27.0 + c;
            float d = q * q + 4.0 * p3 / 27.0;
            float offset = -a / 3.0;

            if (d >= 0.0)
            {
                float z = d == 0 ? 0 : std::sqrt(d);

                float x_0 = (z - q) / 2.0;
                float x_1 = (-z - q) / 2.0;

                float uvx = Math::sign(x_0) * std::cbrt(std::abs(x_0));
                float uvy = Math::sign(x_1) * std::cbrt(std::abs(x_1));

                r_0 = offset + uvx + uvy;
                r_1 = r_0;
                r_2 = r_0;
                return;
            }
            float v = std::acos(-std::sqrt(-27.0 / p3) * q / 2.0) / 3.0;
            float m = Q_cos(v);
            float n = Q_sin(v) * 1.732050808;
            float ttt = std::sqrt(-p / 3.0);
            r_0 = ttt * (m + m) + offset;
            r_1 = ttt * (-n - m) + offset;
            r_2 = ttt * (n - m) + offset;
        }

        double bezierLerp(double a, double b, double c, double t)
        {
            return (1 - t) * (1 - t) * a + 2 * (1 - t) * t * b + t * t * c;
        }

        double linearLerp(double a, double b, double t)
        {
            return b * t + a * (1.0 - t);
        }

        double bezierDerivative(double a, double b, double c, double t)
        {
            return 2 * (-a * (1 - t) - 2 * b * t + b + c * t);
        }

        int8_t signBit(int64_t v)
        {
            if (v == 0)
            {
                return 0;
            }
            else if (v < 0)
            {
                return -1;
            }
            return 1;
        }

        bool isCornerSharp(Math::Vector2i a_1, Math::Vector2i a_2, Math::Vector2i b_1, Math::Vector2i b_2)
        {
            int64_t a_dx = a_2.x - a_1.x;
            int64_t a_dy = a_2.y - a_1.y;
            int64_t b_dx = b_2.x - b_1.x;
            int64_t b_dy = b_2.y - b_1.y;

            if (
                signBit(a_dx) != signBit(b_dx) ||
                signBit(a_dy) != signBit(b_dy) ||
                (a_dx == 0 && b_dx == 0 && a_dy == 0 && b_dy == 0) // bdx & bdy impliead
            )
            {
                return true;
            }

            return std::abs(float(a_dy) / a_dx - float(b_dy) / b_dx) > 1.0;
        }
    }

    typedef int64_t GlyphIndex;
    typedef int64_t CodePoint;

    // Floats are a trascendental lie. I'm scared...
    struct GlyphMetrics
    {
        int64_t advance_width;
        int64_t left_side_bearing;
        int64_t right_side_bearing;
        Math::Rect4i bounds;
    };

    struct GlyphCurve
    {
        struct BezierCurve
        {
        private:
            const float squaredDistance(const float a, const float b) const
            {
                return a * a + b * b;
            }

        public:
            bool is_collinear = false;
            Math::Vector2i tail;
            Math::Vector2i center;
            Math::Vector2i head;

            Math::Vector2i lerp(double t) const
            {
                return Math::Vector2i(
                    std::round(Utils::bezierLerp(tail.x, center.x, head.x, t)),
                    std::round(Utils::bezierLerp(tail.y, center.y, head.y, t)));
            }

            float distance(Math::Vector2i P) const
            {
                if (is_collinear) // Line segment, easiest case
                {
                    int64_t p01_x = head.x - tail.x;
                    int64_t p01_y = head.y - tail.y;
                    int64_t p0p_x = P.x - tail.x;
                    int64_t p0p_y = P.y - tail.y;

                    float t = std::min(
                        1.0f,
                        std::max(
                            0.0f,
                            float(p0p_x * p01_x + p0p_y * p01_y) / float(p01_x * p01_x + p01_y * p01_y)));

                    return std::hypot(
                        Utils::linearLerp(tail.x, head.x, t) - P.x,
                        Utils::linearLerp(tail.y, head.y, t) - P.y);
                }
                else // The tangled
                {
                    // P - P0
                    int64_t p0_x = P.x - tail.x; //
                    int64_t p0_y = P.y - tail.y;
                    // P1 - P0
                    int64_t p1_x = center.x - tail.x; //
                    int64_t p1_y = center.y - tail.y;
                    // P2 - 2P1 + P0
                    int64_t p2_x = head.x - 2 * center.x + tail.x; //
                    int64_t p2_y = head.y - 2 * center.y + tail.y;

                    int64_t a = p2_x * p2_x + p2_y * p2_y;
                    int64_t b = 3 * (p1_x * p2_x + p1_y * p2_y);
                    int64_t c = 2 * (p1_x * p1_x + p1_y * p1_y) - (p2_x * p0_x + p2_y * p0_y);
                    int64_t d = -(p1_x * p0_x + p1_y * p0_y);

                    float t_0 = 0;
                    float t_1 = 0;
                    float t_2 = 0;

                    Utils::solveCubic(a, b, c, d, t_0, t_1, t_2);

                    t_0 = std::min(1.0f, std::max(0.0f, t_0));
                    t_1 = std::min(1.0f, std::max(0.0f, t_1));
                    t_2 = std::min(1.0f, std::max(0.0f, t_2));

                    float d_0 = squaredDistance(
                        P.x - Utils::bezierLerp(tail.x, center.x, head.x, t_0),
                        P.y - Utils::bezierLerp(tail.y, center.y, head.y, t_0));
                    float d_1 = squaredDistance(
                        P.x - Utils::bezierLerp(tail.x, center.x, head.x, t_1),
                        P.y - Utils::bezierLerp(tail.y, center.y, head.y, t_1));
                    float d_2 = squaredDistance(
                        P.x - Utils::bezierLerp(tail.x, center.x, head.x, t_2),
                        P.y - Utils::bezierLerp(tail.y, center.y, head.y, t_2));

                    return std::sqrt(
                        std::min(std::min(d_0, d_1), d_2));
                }
            }
        };

        std::vector<BezierCurve> curves;
        std::vector<std::pair<uint64_t, uint64_t>> contour_limits;

        void merge(const GlyphCurve &other)
        {
            uint32_t offset = contour_limits.empty() ? 0 : contour_limits.back().second + 1;

            curves.insert(curves.end(), other.curves.begin(), other.curves.end());

            for (const auto &ctl : other.contour_limits)
            {
                contour_limits.push_back({ctl.first + offset, ctl.second + offset});
            }
        }
    };

    class Font
    {
    public:
        virtual ~Font() = default;
        virtual GlyphIndex getId(CodePoint codepoint) const = 0;
        virtual GlyphMetrics getMetrics(GlyphIndex index) const = 0;
        virtual GlyphCurve getCurves(GlyphIndex index) const = 0;
        virtual double getScale(double font_size) const = 0;
    };

    /**
        @brief Font loader function type.
        @param mem memory to read.
        @param size memory size.
        @return An optional Font reference
    */
    typedef std::shared_ptr<Font> LoaderFn_t(const void *, size_t);

    class FontLoaderFactory
    {
        inline static std::vector<LoaderFn_t *> loaders{};

    public:
        /**
            @brief Registers a new loader.
            @param loader Loader function pointer.
        */
        static void registerLoader(LoaderFn_t *loader)
        {
            if (loader)
            {
                loaders.push_back(loader);
            }
        }

        friend const std::shared_ptr<Font> fromMemory(const void *mem, size_t size);
    };

    /**
        @brief Try to load a font from memory.
        @param mem memory to read.
        @param size memory size.
        @return An optional Font reference
    */
    const std::shared_ptr<Font> fromMemory(const void *mem, size_t size)
    {
        for (LoaderFn_t *ldr : FontLoaderFactory::loaders)
        {
            if (const auto fnt = ldr(mem, size))
            {
                return fnt;
            }
        }

        return nullptr;
    }

    /**
        @brief Try to load a font from a path.
        @param path filesystem path.
        @return An optional Font reference
    */
    const std::shared_ptr<Font> fromPath(std::filesystem::path path)
    {
        if (const auto &fco = Utils::readFile(path))
        {
            const auto &file_content = *fco;

            return fromMemory(file_content.data(), file_content.size());
        }

        return nullptr;
    }

    /* class GlyphAtlas
    {
        struct GlyphEntry
        {
            GlyphCurve curve;
            Math::Rect4f bounds;
            uint16_t page;
            double uv[4];
        };

        struct AtlasPage
        {
            std::unique_ptr<uint8_t[]> image;
            Math::Vector2i size;
        };

        std::vector<AtlasPage> pages;
        std::map<GlyphIndex, std::shared_ptr<GlyphEntry>> entries;

        std::shared_ptr<GlyphEntry> getEntry(GlyphIndex index)
        {
            auto k = entries.find(index);
            if (k != entries.end())
            {
                return k->second;
            }
            return nullptr;
        }

        GlyphAtlas() {}
    }; */

    namespace Packer
    {
        std::map<GlyphIndex, Math::Rect4i> getScaledRects(const Font *font, const std::vector<GlyphIndex> &glyph_indices, const double font_size)
        {
            std::map<GlyphIndex, Math::Rect4i> rects = {};
            // Get font scale
            double scale = font->getScale(font_size);
            // Get scaled rects
            for (const auto &gi : glyph_indices)
            {
                GlyphMetrics mt = font->getMetrics(gi);
                Math::Rect4i rect;

                rect.x = 0;
                rect.y = 0;
                rect.w = std::round(mt.bounds.w * scale);
                rect.h = std::round(mt.bounds.h * scale);
                rects[gi] = rect;
            }
            return rects;
        }

        void basicSortPack(std::map<GlyphIndex, Math::Rect4i> &rects, int64_t canvas_width, int64_t canvas_height)
        {
            if (rects.empty())
            {
                return;
            }
            // Get pair pointers
            typedef std::pair<const GlyphIndex, Math::Rect4i> pair_t;
            std::vector<pair_t *> pairs;
            for (auto &p : rects)
            {
                pairs.push_back(&p);
            }
            // Sort rects by width and then by reverse h
            std::sort(
                pairs.begin(),
                pairs.end(),
                [](pair_t *a, pair_t *b)
                {
                    return a->second.w > b->second.w || (a->second.w == b->second.w && a->second.h > b->second.h);
                });
            // Init packing vars
            int64_t advance_x = 0;
            int64_t advance_y = 0; // resets to 0
            int64_t bin_width = pairs[0]->second.w;
            int64_t padding = 1;
            for (auto &pair : pairs)
            {
                if (
                    pair->second.w > bin_width ||
                    advance_y + pair->second.h + padding > canvas_height)
                {
                    advance_x += bin_width + padding;
                    advance_y = 0;
                    bin_width = pair->second.w;
                    if (advance_x + bin_width > canvas_width)
                    {
                        break;
                    }
                }
                pair->second.x = advance_x;
                pair->second.y = advance_y;
                advance_y += pair->second.h + padding;
            }
        }
    }

    struct RasterParameters
    {
        GlyphCurve *curve;
        GlyphMetrics *metrics;
        Math::Rect4i region;
    };

    /**
        @brief Font loader function type.
        @param mem memory of the atlas image.
        @param width atlas width.
        @param height atlas height.
        @param params raster params
        @return An optional Font reference
    */
    typedef void RasterizerFn(const void *, size_t, size_t, RasterParameters);

    namespace Rasterizer
    {
        union Color32
        {
            struct
            {
                uint8_t r;
                uint8_t g;
                uint8_t b;
                uint8_t a;
            };
            uint32_t color;
        };

        void bulk(const Font *font, std::map<GlyphIndex, Math::Rect4i> &rects, const void *atlas, size_t atlas_width, size_t atlas_height, RasterizerFn raster_fn)
        {
            for (const auto &p : rects)
            {
                RasterParameters params;

                auto gm = font->getMetrics(p.first);
                auto gcur = font->getCurves(p.first);

                params.metrics = &gm;
                params.curve = &gcur;
                params.region = Math::Rect4i(
                    p.second.x,
                    p.second.y,
                    p.second.w,
                    p.second.h);
                
                raster_fn(atlas, atlas_width, atlas_height, params);
            }
        }

        void scanline(const void *data, size_t atlas_width, size_t atlas_height, RasterParameters params)
        {
            uint8_t *atlas = (uint8_t *)data;
            // Calculate paint region
            uint32_t region_min_x = std::max(0L, params.region.x);
            uint32_t region_min_y = std::max(0L, params.region.y);
            uint32_t region_max_x = std::min((long)atlas_width, region_min_x + std::abs(params.region.w));
            uint32_t region_max_y = std::min((long)atlas_height, region_min_y + std::abs(params.region.h));
            // Count of glyph space pixels per atlas space pixels.
            double scale = std::max(
                double(params.metrics->bounds.w) / (region_max_x - region_min_x),
                double(params.metrics->bounds.h) / (region_max_y - region_min_y));
            // Bake any transformed glyphs
            std::vector<GlyphCurve::BezierCurve> curves;
            for (const auto &curve : params.curve->curves)
            {
                GlyphCurve::BezierCurve cur = curve;
                // Move curve coordinates to first quadrant (No negatives)
                cur.tail.x = cur.tail.x - params.metrics->bounds.x;
                cur.tail.y = cur.tail.y - params.metrics->bounds.y;
                cur.center.x = cur.center.x - params.metrics->bounds.x;
                cur.center.y = cur.center.y - params.metrics->bounds.y;
                cur.head.x = cur.head.x - params.metrics->bounds.x;
                cur.head.y = cur.head.y - params.metrics->bounds.y;
                // Later use a better data struct
                curves.push_back(cur);
                // Good idea: Keep boundary iterators
            }

            for (uint32_t y = region_min_y; y < region_max_y; y++)
            {
                // With the offset, the scanline will never overlap with the curve's extremes (Unless rounding errors, floats are a trascendental lie)
                double gly_y = std::round((y - region_min_y) * scale) + 0.5;
                // Calculate any intersections in the half displaced grid
                std::map<float, std::pair<float, float>> intersections;

                for (const auto &cur : curves)
                {
                    // Skip any unneeded curves
                    if (
                        std::min(cur.tail.y, std::min(cur.center.y, cur.head.y)) > gly_y ||
                        std::max(cur.tail.y, std::max(cur.center.y, cur.head.y)) < gly_y)
                    {
                        continue;
                    }
                    // Handle collinear case
                    if (cur.is_collinear)
                    {
                        if (cur.head.y - cur.tail.y == 0)
                        {
                            continue;
                        }
                        double t = double(gly_y - cur.tail.y) / (cur.head.y - cur.tail.y);

                        if (t >= 0.0f && t <= 1.0f)
                        {
                            float x = Utils::bezierLerp(cur.tail.x, cur.center.x, cur.head.x, t) / scale;
                            intersections[x] = {
                                Utils::bezierDerivative(cur.tail.x, cur.center.x, cur.head.x, t),
                                Utils::bezierDerivative(cur.tail.y, cur.center.y, cur.head.y, t)};
                        }
                    }
                    // Handle bezier case
                    else
                    {
                        double t_0, t_1;
                        Utils::solveQuadratic(
                            cur.tail.y - 2.0f * cur.center.y + cur.head.y,
                            2.0f * (cur.center.y - cur.tail.y),
                            cur.tail.y - gly_y,
                            t_0, t_1);
                        if (!(std::isnan(t_0) || std::isnan(t_1)))
                        {
                            // Push intersections if not edges
                            if (t_0 >= 0 && t_0 <= 1)
                            {
                                float x = Utils::bezierLerp(cur.tail.x, cur.center.x, cur.head.x, t_0) / scale;
                                intersections[x] = {
                                    Utils::bezierDerivative(cur.tail.x, cur.center.x, cur.head.x, t_0),
                                    Utils::bezierDerivative(cur.tail.y, cur.center.y, cur.head.y, t_0)};
                            }
                            if (t_1 >= 0 && t_1 <= 1)
                            {
                                float x = Utils::bezierLerp(cur.tail.x, cur.center.x, cur.head.x, t_1) / scale;
                                intersections[x] = {
                                    Utils::bezierDerivative(cur.tail.x, cur.center.x, cur.head.x, t_1),
                                    Utils::bezierDerivative(cur.tail.y, cur.center.y, cur.head.y, t_1)};
                            }
                        }
                    }
                }
                // If anything goes wrong, skip the scanline
                if (intersections.size() % 2 != 0)
                {
                    for (const auto &i : intersections)
                    {
                        std::printf("III: %f %f %f\n", i.first, i.second.first, i.second.second);
                    }
                    std::printf("fucked up %f\n", gly_y);
                    exit(0);
                }
                // Get full x-segments with winding
                int winding_value = 0;
                for (auto it = intersections.begin(); it != intersections.end();)
                {
                    auto start = it++;
                    // Wind for first intersection
                    winding_value += std::signbit(start->second.second) ? -1 : 1;
                    // Go ahead until we exit the outline
                    auto end = it;
                    while (it != intersections.end() && winding_value != 0)
                    {
                        end = it++;
                        winding_value += std::signbit(end->second.second) ? -1 : 1;
                    }
                    // Get extrema
                    uint32_t x_A = std::max(0, std::min((int)atlas_width, (int)std::round(start->first)));
                    uint32_t x_B = std::max(0, std::min((int)atlas_width, (int)std::round(end->first)));
                    // Create routine fill x with derivative gradient
                    for (uint32_t x = x_A + region_min_x; x < x_B + region_min_x; x++)
                    {
                        Color32 &pix = ((Color32 *)atlas)[y * atlas_width + x];

                        pix.a = 255;
                        pix.r = 255;
                    }
                }
            }
        }

        void SDF(const void *data, size_t atlas_width, size_t atlas_height, RasterParameters params, int8_t thickness)
        {
            uint8_t *atlas = (uint8_t *)data;

            int8_t one_side_thickness = thickness / 2;

            uint32_t region_min_x = std::max(0L, params.region.x);
            uint32_t region_min_y = std::max(0L, params.region.y);
            uint32_t region_max_x = std::min((long)atlas_width, region_min_x + std::abs(params.region.w));
            uint32_t region_max_y = std::min((long)atlas_height, region_min_y + std::abs(params.region.h));

            double scale = std::min(
                (region_max_x - region_min_x - thickness) / double(params.metrics->bounds.w),
                (region_max_y - region_min_y - thickness) / double(params.metrics->bounds.h));

            float distance_scale = 255 / (thickness / scale);

            for (const auto &curve : params.curve->curves)
            {
                int64_t scaled_one_side_thickness = one_side_thickness / scale;
                GlyphCurve::BezierCurve cur = curve;
                // Move curve coordinates to first quadrant (No negatives)
                cur.tail.x = cur.tail.x - params.metrics->bounds.x + scaled_one_side_thickness;
                cur.tail.y = cur.tail.y - params.metrics->bounds.y + scaled_one_side_thickness;
                cur.center.x = cur.center.x - params.metrics->bounds.x + scaled_one_side_thickness;
                cur.center.y = cur.center.y - params.metrics->bounds.y + scaled_one_side_thickness;
                cur.head.x = cur.head.x - params.metrics->bounds.x + scaled_one_side_thickness;
                cur.head.y = cur.head.y - params.metrics->bounds.y + scaled_one_side_thickness;
                // Compute paint bounds
                uint32_t paint_min_x = std::max((double)region_min_x, std::min((double)region_max_x, std::min(cur.tail.x, std::min(cur.center.x, cur.head.x)) * scale + region_min_x - one_side_thickness));
                uint32_t paint_min_y = std::max((double)region_min_y, std::min((double)region_max_y, std::min(cur.tail.y, std::min(cur.center.y, cur.head.y)) * scale + region_min_y - one_side_thickness));
                uint32_t paint_max_x = std::max((double)region_min_x, std::min((double)region_max_x, std::max(cur.tail.x, std::max(cur.center.x, cur.head.x)) * scale + region_min_x + one_side_thickness));
                uint32_t paint_max_y = std::max((double)region_min_y, std::min((double)region_max_y, std::max(cur.tail.y, std::max(cur.center.y, cur.head.y)) * scale + region_min_y + one_side_thickness));

                for (uint32_t paint_y = paint_min_y; paint_y < paint_max_y; paint_y++)
                {
                    int32_t y = paint_y / scale;
                    int32_t x_start = paint_min_x, x_end = paint_max_x;

                    int32_t dst = 0;
                    do
                    {
                        dst = cur.distance(Math::Vector2i(x_start / scale, y)) * scale - scaled_one_side_thickness;
                        if (dst <= 0)
                        {
                            break;
                        }
                        x_start += dst;
                    } while (x_start < x_end);

                    do
                    {
                        dst = cur.distance(Math::Vector2i(x_end / scale, y)) * scale - scaled_one_side_thickness;
                        if (dst <= 0)
                        {
                            break;
                        }
                        x_end -= dst;
                    } while (x_start < x_end);

                    if (x_start > x_end)
                    {
                        continue;
                    }

                    for (uint32_t paint_x = x_start; paint_x < x_end; paint_x++)
                    {
                        Color32 &pix = ((Color32 *)atlas)[paint_y * atlas_width + paint_x];

                        uint8_t distance = 0x7f - (int8_t)std::min(127.0f, cur.distance(Math::Vector2i(paint_x / scale, y)) * distance_scale);

                        pix.a = 255;
                        pix.r = std::max(uint8_t(pix.r & 0xff), distance);
                    }
                }
            }
        }
    }
}

#endif