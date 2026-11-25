module palette;

import wuquant : Color;

private struct HSL
{
    double h, s, l;
}

private double minD(double a, double b) @nogc nothrow
{
    return a < b ? a : b;
}

private double maxD(double a, double b) @nogc nothrow
{
    return a > b ? a : b;
}

private HSL rgb2hsl(Color c) @nogc nothrow
{
    double r = c.r;
    double g = c.g;
    double b = c.b;

    double mn = minD(r, minD(g, b));
    double mx = maxD(r, maxD(g, b));
    double delta = mx - mn;

    // Explicit initialization: floating .init is NaN in D.
    HSL result;
    result.h = 0.0;
    result.s = 0.0;
    result.l = (mn + mx) / 510.0;

    if (result.l > 0.0 && result.l < 1.0)
    {
        result.s = delta / (result.l < 0.5
            ? mx + mn
            : 510.0 - mx - mn);
    }

    if (delta > 0.0)
    {
        if (mx == r)
            result.h = (g - b) / delta;
        else if (mx == g)
            result.h = 2.0 + (b - r) / delta;
        else
            result.h = 4.0 + (r - g) / delta;

        result.h *= 60.0;

        if (result.h < 0.0)
            result.h += 360.0;
    }

    return result;
}

private int hueGroup(double hue) @nogc nothrow
{
    enum int Segments = 10;
    enum double Segment = 360.0 / Segments;
    enum double HalfSegment = Segment / 2.0;

    for (int i = 1; i < Segments; ++i)
    {
        double mid = Segment * i - HalfSegment;

        if (hue >= mid && hue < mid + Segment)
            return i;
    }

    return 0;
}

private double luminosity(Color c) @nogc nothrow
{
    // Composite the color over white, as image-q's luminosity(true) does.
    double alpha = c.a;
    double r = minD(255.0, 255.0 - alpha + alpha * c.r / 255.0);
    double g = minD(255.0, 255.0 - alpha + alpha * c.g / 255.0);
    double b = minD(255.0, 255.0 - alpha + alpha * c.b / 255.0);

    return r * 0.2126 + g * 0.7152 + b * 0.0722;
}

private int compare(Color a, Color b) @nogc nothrow
{
    HSL ha = rgb2hsl(a);
    HSL hb = rgb2hsl(b);

    int hueA = (a.r == a.g && a.g == a.b)
        ? 0 : 1 + hueGroup(ha.h);

    int hueB = (b.r == b.g && b.g == b.b)
        ? 0 : 1 + hueGroup(hb.h);

    if (hueA != hueB)
        return hueA < hueB ? -1 : 1;

    double la = luminosity(a);
    double lb = luminosity(b);

    if (la != lb)
        return la > lb ? -1 : 1;

    int sa = cast(int)(ha.s * 100.0);
    int sb = cast(int)(hb.s * 100.0);

    if (sa != sb)
        return sa < sb ? -1 : 1;

    return 0;
}

void sortPalette(ref Color[] colors) @nogc nothrow
{
    // Stable insertion sort; palettes are limited to at most 256 entries.
    foreach (i; 1 .. colors.length)
    {
        Color value = colors[i];
        size_t j = i;

        while (j > 0 && compare(value, colors[j - 1]) < 0)
        {
            colors[j] = colors[j - 1];
            --j;
        }

        colors[j] = value;
    }
}
