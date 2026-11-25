module main;

import std.stdio;
import std.conv;

import dlib;
import bindbc.freeimage;

import wuquant;
import palette: sortPalette;
import atkinson: ditheringAtkinson;

import std.string: toStringz;

bool saveIndexedPNG(
    string filename,
    uint width,
    uint height,
    const(ubyte)[] rgba,
    const(Color)[] colors
)
{
    auto dib = FreeImage_Allocate(
        cast(int)width,
        cast(int)height,
        8,
        0, 0, 0
    );

    if (dib is null)
        return false;

    scope(exit) FreeImage_Unload(dib);

    // Build RGB palette and transparency table.
    auto palette = FreeImage_GetPalette(dib);
    ubyte[] alpha = new ubyte[colors.length];
    bool hasAlpha = false;

    foreach (i, c; colors)
    {
        palette[i].rgbRed = c.r;
        palette[i].rgbGreen = c.g;
        palette[i].rgbBlue = c.b;
        palette[i].rgbReserved = 0;

        alpha[i] = c.a;
        if (c.a != 255)
            hasAlpha = true;
    }

    if (hasAlpha)
    {
        FreeImage_SetTransparencyTable(
            dib,
            alpha.ptr,
            cast(int)alpha.length
        );
        FreeImage_SetTransparent(dib, TRUE);
    }

    // Map exact palette colors to their indices.
    uint[uint] indexMap;

    foreach (i, c; colors)
    {
        uint key =
            cast(uint)c.r |
            (cast(uint)c.g << 8) |
            (cast(uint)c.b << 16) |
            (cast(uint)c.a << 24);

        indexMap[key] = cast(uint)i;
    }

    // FreeImage scanlines are bottom-up.
    foreach (y; 0 .. height)
    {
        auto scan = FreeImage_GetScanLine(
            dib,
            cast(int)(height - 1 - y)
        );

        foreach (x; 0 .. width)
        {
            size_t p = (cast(size_t)y * width + x) * 4;

            uint key =
                cast(uint)rgba[p] |
                (cast(uint)rgba[p + 1] << 8) |
                (cast(uint)rgba[p + 2] << 16) |
                (cast(uint)rgba[p + 3] << 24);

            scan[x] = cast(ubyte)indexMap[key];
        }
    }

    return FreeImage_Save(
        FIF_PNG,
        dib,
        filename.toStringz
    ) != 0;
}

void main(string[] args)
{
    if (args.length < 2)
    {
        writeln("Usage: dpngquant <input.png> [output.png] [numcolors]");
        return;
    }

    string inputFilename = args[1];
    string outputFilename = args.length >= 3
        ? args[2]
        : "output.png";
    uint numColors = args.length >= 4
        ? args[3].to!uint
        : 230;
    if (numColors > 256)
        numColors = 256;

    FISupport support = loadFreeImage();
    if (support != FISupport.fi318)
    {
        writeln("Failed to load FreeImage: ", support);
        return;
    }

    auto img = loadImage(inputFilename);
    if (img is null)
    {
        writeln("Failed to load image: ", inputFilename);
        return;
    }

    uint width = cast(uint)img.width;
    uint height = cast(uint)img.height;

    if (width == 0 || height == 0)
    {
        writeln("Invalid image dimensions");
        return;
    }

    size_t pixelCount = cast(size_t)width * height;
    if (pixelCount > size_t.max / 4)
    {
        writeln("Image is too large");
        return;
    }

    writeln("Image: ", width, "x", height);
    writeln("Building WuQuant palette...");

    ubyte[] rgba = new ubyte[pixelCount * 4];

    foreach (y; 0 .. height)
    foreach (x; 0 .. width)
    {
        auto c = img[cast(int)x, cast(int)y];
        size_t i = (cast(size_t)y * width + x) * 4;

        rgba[i]     = cast(ubyte)(c.r * 255.0f + 0.5f);
        rgba[i + 1] = cast(ubyte)(c.g * 255.0f + 0.5f);
        rgba[i + 2] = cast(ubyte)(c.b * 255.0f + 0.5f);
        rgba[i + 3] = cast(ubyte)(c.a * 255.0f + 0.5f);
    }

    WuQuant wu;
    auto colors = wu.quantize(rgba, numColors);

    if (colors.length == 0)
    {
        writeln("WuQuant produced an empty palette");
        return;
    }

    sortPalette(colors);
    writeln("Palette colors: ", colors.length);

    // Optional palette preview, retaining the original workflow.
    auto paletteImage = image(cast(int)colors.length, 1);

    foreach (i, c; colors)
    {
        paletteImage[cast(int)i, 0] = Color4f(
            c.r / 255.0f,
            c.g / 255.0f,
            c.b / 255.0f,
            c.a / 255.0f
        );
    }

    paletteImage.saveImage("palette.png");

    writeln("Applying Atkinson dithering...");
    ditheringAtkinson(rgba, width, height, colors);
    
    if (!saveIndexedPNG(outputFilename, width, height, rgba, colors))
    {
        writeln("Failed to save indexed PNG");
        return;
    }
    
    writeln("Saved: ", outputFilename);
    writeln("Saved: palette.png");
}
