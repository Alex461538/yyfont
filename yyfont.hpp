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
        virtual GlyphIndex getId(CodePoint codepoint) = 0;
        virtual GlyphMetrics getMetrics(GlyphIndex index) = 0;
        virtual GlyphCurve getCurves(GlyphIndex index) = 0;
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
    typedef void Rasterizer2DFn(const void *, size_t, size_t, RasterParameters);

    namespace Rasterizer
    {
        void SDF(const void * data, size_t atlas_width, size_t atlas_height, RasterParameters params, uint8_t thickness)
        {
            uint8_t *atlas = (uint8_t*)data;

            uint8_t half_thickness = thickness / 2;

            double x_scale = params.region.w / double(params.metrics->bounds.w);
            double y_scale = params.region.h / double(params.metrics->bounds.h);

            int64_t metrics_dx = -params.metrics->bounds.x;
            int64_t metrics_dy = -params.metrics->bounds.y;

            for (const auto &curve : params.curve->curves)
            {
                GlyphCurve::BezierCurve cur = curve;
                // Convert outlines to pixel coordinates
                cur.tail.x = (cur.tail.x + metrics_dx) * x_scale + half_thickness + params.region.x;
                cur.tail.y = (cur.tail.y + metrics_dy) * y_scale + half_thickness + params.region.y;
                cur.center.x = (cur.center.x + metrics_dx) * x_scale + half_thickness + params.region.x;
                cur.center.y = (cur.center.y + metrics_dy) * y_scale + half_thickness + params.region.y;
                cur.head.x = (cur.head.x + metrics_dx) * x_scale + half_thickness + params.region.x;
                cur.head.y = (cur.head.y + metrics_dy) * y_scale + half_thickness + params.region.y;

                int64_t min_x = std::min(cur.tail.x, std::min(cur.center.x, cur.head.x)) - half_thickness;
                int64_t min_y = std::min(cur.tail.y, std::min(cur.center.y, cur.head.y)) - half_thickness;
                
                int64_t max_x = std::max(cur.tail.x, std::max(cur.center.x, cur.head.x)) + half_thickness;
                int64_t max_y = std::max(cur.tail.y, std::max(cur.center.y, cur.head.y)) + half_thickness;

                for (int y = min_y; y < max_y; y++)
                {
                    int64_t x_start = min_x;
                    int64_t x_end = max_x;

                    float dst = 0;
                    do
                    {
                        dst = cur.distance(Math::Vector2i(x_start, y)) - thickness / 2;
                        x_start += dst;
                    } while (dst > 1);

                    do
                    {
                        dst = cur.distance(Math::Vector2i(x_end, y)) - thickness / 2;
                        x_end -= dst;
                    } while (dst > 1);

                    for (int x = x_start; x < x_end; x++)
                    {
                        uint32_t *pix = (uint32_t *)&(atlas[(y * atlas_width + x) * 4]);

                        int8_t distance = 0x7f - (int8_t)std::min(
                                                     127.0f,
                                                     cur.distance(Math::Vector2i(x, y)) * 255 / thickness);

                        *pix = 0xff000000 | std::max(int8_t(*pix & 0xff), distance);
                    }
                }
            }
        }
    }
}

#endif