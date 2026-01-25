# ☔️ yyfont.hpp

A modular font parser written in moderately modern c++.

![](img/atlas.png)

## It's easy to use btw 🤙

```cpp
#include "yyfont.hpp" // Import the main module
#include "yyfont.sfnt.hpp" // And any format you'd like
...
Font::SFNT::init(); // Register loaders for those
```

```cpp
// Load your favorite font!
if (const auto &test_font = Font::fromPath("samples/nj.ttf"))
{
    // Get internal ID's for unicode code points.
    Font::GlyphIndex gi = test_font->getId(u'字');
    // Access any properties about them.
    Font::GlyphMetrics gm = test_font->getMetrics(gi);
    // Or get their curves too 😏.
    Font::GlyphCurve letter = test_font->getCurves(gi);
}
```

```cpp
// Get your preferred charset
std::vector<Font::GlyphIndex> indices = { ... };
// And create an atlas for it!
int64_t atlas_width = 512;
int64_t atlas_height = 512;

auto atlas = std::make_unique<uint8_t[]>(atlas_width * atlas_height * 4);

// Get scaled rects for the glyphs (font size = 50px)
auto rects = Font::Packer::getScaledRects(test_font.get(), indices, 50);
// Pack them with any algorithm 📦️
Font::Packer::basicSortPack(rects, atlas_width, atlas_height);
// And rasterize everything together 🖍️
Font::Rasterizer::bulk(test_font.get(), rects, atlas.get(), atlas_width, atlas_height, Font::Rasterizer::scanline);
```

## Register new formats

There are allways obscurer formats out there. If you're brave enough, go ahead and extend the system!

```cpp
class JairoFragmentsFont : Font
{
    public:
        // Tries to load a font of your type. If it fails, it should return nullptr.
        static std::shared_ptr<Font> Loader(const void *mem, size_t size) { ... }
};
```
**NOTE: You don't own the given memory, if you need something, copy it please 😊**

```cpp
// Register your new loader.
Font::FontLoaderFactory::registerLoader( JairoFragmentsFont::Loader );
// & without aditional format madness, that's it.
const auto &test_font = Font::fromPath("samples/font.jfp")
```

The factory will try to load a file with every loader until any claims it.

## License

Do whatever you want with this, I don't give a... [Mit License](LICENSE)

### Resources

- Online font inspector: [FontDrop](https://fontdrop.info/)
- TrueType spec (Apple): [TrueType Reference Manual](https://developer.apple.com/fonts/TrueType-Reference-Manual/)
- OpenType spec (Microsoft): [OpenType® Specification Version 1.9.1](https://www.microsoft.com/typography/otspec/)
- WTF opentype with a font inside a font: [The Compact Font Format Specification](https://adobe-type-tools.github.io/font-tech-notes/pdfs/5176.CFF.pdf)

#### Protect chocobos worldwide 🐥