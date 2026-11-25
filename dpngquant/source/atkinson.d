module atkinson;

import wuquant: Color;

private ubyte clampRound(double value) @nogc nothrow
{
    if (!(value > 0.0))
        return 0;

    if (value >= 255.0)
        return 255;

    return cast(ubyte)(value + 0.5);
}

private uint distance(Color a, Color b) @nogc nothrow
{
    int dr = cast(int)a.r - cast(int)b.r;
    int dg = cast(int)a.g - cast(int)b.g;
    int db = cast(int)a.b - cast(int)b.b;
    int da = cast(int)a.a - cast(int)b.a;

    return cast(uint)(
        dr * dr +
        dg * dg +
        db * db +
        da * da
    );
}

private Color nearestColor(Color source, const(Color)[] colors)
    @nogc nothrow
{
    Color result = colors[0];
    uint bestDistance = uint.max;

    foreach (color; colors)
    {
        uint d = distance(source, color);

        if (d < bestDistance)
        {
            bestDistance = d;
            result = color;
        }
    }

    return result;
}

void ditheringAtkinson(
    ubyte[] rgba,
    uint width,
    uint height,
    const(Color)[] colors
)
{
    if (width == 0 || height == 0 || colors.length == 0)
        return;

    size_t expected = cast(size_t)width * height * 4;
    if (rgba.length < expected)
        return;

    const double divider = 8.0;

    foreach (y; 0 .. height)
    foreach (x; 0 .. width)
    {
        size_t pixel = (cast(size_t)y * width + x) * 4;

        Color original = Color(
            rgba[pixel],
            rgba[pixel + 1],
            rgba[pixel + 2],
            rgba[pixel + 3]
        );

        Color corrected = Color(
            clampRound(original.r),
            clampRound(original.g),
            clampRound(original.b),
            clampRound(original.a)
        );

        Color reduced = nearestColor(corrected, colors);

        rgba[pixel]     = reduced.r;
        rgba[pixel + 1] = reduced.g;
        rgba[pixel + 2] = reduced.b;
        rgba[pixel + 3] = reduced.a;

        // Error is calculated from the pixel value before quantization.
        double er = cast(int)original.r - cast(int)reduced.r;
        double eg = cast(int)original.g - cast(int)reduced.g;
        double eb = cast(int)original.b - cast(int)reduced.b;
        double ea = cast(int)original.a - cast(int)reduced.a;

        // Atkinson diffusion pattern.
        static immutable int[2][6] Kernel = [
            [ 1, 0],
            [ 2, 0],
            [-1, 1],
            [ 0, 1],
            [ 1, 1],
            [ 0, 2]
        ];

        foreach (k; 0 .. 6)
        {
            int targetX = cast(int)x + Kernel[k][0];
            int targetY = cast(int)y + Kernel[k][1];

            if (targetX < 0 || targetX >= cast(int)width ||
                targetY < 0 || targetY >= cast(int)height)
                continue;

            size_t target =
                (cast(size_t)targetY * width + cast(size_t)targetX) * 4;

            rgba[target] = clampRound(rgba[target] + er / divider);
            rgba[target + 1] = clampRound(rgba[target + 1] + eg / divider);
            rgba[target + 2] = clampRound(rgba[target + 2] + eb / divider);
            rgba[target + 3] = clampRound(rgba[target + 3] + ea / divider);
        }
    }
}
