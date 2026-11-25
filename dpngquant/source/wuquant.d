module wuquant;

struct Color
{
    ubyte r, g, b, a;
}

private struct Cube
{
    int r0, r1;
    int g0, g1;
    int b0, b1;
    int a0, a1;
    int volume;
}

private struct MaxResult
{
    double value;
    int position;
}

private int toInt32(double x) @nogc nothrow
{
    // Values used by WuQuant are finite and within the useful ToInt32 range.
    return cast(int)cast(uint)cast(long)x;
}

private double absd(double x) @nogc nothrow
{
    return x < 0.0 ? -x : x;
}

struct WuQuant
{
    enum int Bits = 5;
    enum int MaxIndex = 1 << Bits;
    enum int Side = MaxIndex + 1;

    enum size_t CellCount =
        cast(size_t)Side * Side * Side * Side;

    enum size_t PlaneSize =
        cast(size_t)Side * Side * Side;

    private double[] weights;
    private double[] momentsR;
    private double[] momentsG;
    private double[] momentsB;
    private double[] momentsA;
    private double[] moments;

    private double[] xw, xr, xg, xb, xa, xm;
    private double[] area, areaR, areaG, areaB, areaA, areaM;

    private Cube[] cubes;
    private uint colors;

    private size_t index(int a, int r, int g, int b) const
        @nogc nothrow
    {
        return (((cast(size_t)a * Side + r) * Side + g) * Side + b);
    }

    private size_t index3(int r, int g, int b) const
        @nogc nothrow
    {
        return (cast(size_t)r * Side + g) * Side + b;
    }

    private void initialize(uint count)
    {
        colors = count;

        // D floating-point .init is NaN. Every accumulator is explicitly zeroed.
        weights = new double[CellCount];
        momentsR = new double[CellCount];
        momentsG = new double[CellCount];
        momentsB = new double[CellCount];
        momentsA = new double[CellCount];
        moments = new double[CellCount];

        weights[] = 0.0;
        momentsR[] = 0.0;
        momentsG[] = 0.0;
        momentsB[] = 0.0;
        momentsA[] = 0.0;
        moments[] = 0.0;

        xw = new double[PlaneSize];
        xr = new double[PlaneSize];
        xg = new double[PlaneSize];
        xb = new double[PlaneSize];
        xa = new double[PlaneSize];
        xm = new double[PlaneSize];

        xw[] = 0.0;
        xr[] = 0.0;
        xg[] = 0.0;
        xb[] = 0.0;
        xa[] = 0.0;
        xm[] = 0.0;

        area = new double[Side];
        areaR = new double[Side];
        areaG = new double[Side];
        areaB = new double[Side];
        areaA = new double[Side];
        areaM = new double[Side];

        area[] = 0.0;
        areaR[] = 0.0;
        areaG[] = 0.0;
        areaB[] = 0.0;
        areaA[] = 0.0;
        areaM[] = 0.0;

        cubes = new Cube[count];

        cubes[0].r1 = MaxIndex;
        cubes[0].g1 = MaxIndex;
        cubes[0].b1 = MaxIndex;
        cubes[0].a1 = MaxIndex;
        cubes[0].volume = MaxIndex * MaxIndex * MaxIndex * MaxIndex;
    }

    private void addColor(Color c) @nogc nothrow
    {
        enum int Shift = 8 - Bits;

        int r = (c.r >> Shift) + 1;
        int g = (c.g >> Shift) + 1;
        int b = (c.b >> Shift) + 1;
        int a = (c.a >> Shift) + 1;

        size_t i = index(a, r, g, b);

        weights[i] += 1.0;
        momentsR[i] += c.r;
        momentsG[i] += c.g;
        momentsB[i] += c.b;
        momentsA[i] += c.a;

        moments[i] +=
            cast(double)c.r * c.r +
            cast(double)c.g * c.g +
            cast(double)c.b * c.b +
            cast(double)c.a * c.a;
    }

    // Wu's four-dimensional cumulative moments.
    private void calculateMoments()
    {
        foreach (a; 1 .. Side)
        {
            xw[] = 0.0;
            xr[] = 0.0;
            xg[] = 0.0;
            xb[] = 0.0;
            xa[] = 0.0;
            xm[] = 0.0;

            foreach (r; 1 .. Side)
            {
                area[] = 0.0;
                areaR[] = 0.0;
                areaG[] = 0.0;
                areaB[] = 0.0;
                areaA[] = 0.0;
                areaM[] = 0.0;

                foreach (g; 1 .. Side)
                {
                    double line = 0.0;
                    double lineR = 0.0;
                    double lineG = 0.0;
                    double lineB = 0.0;
                    double lineA = 0.0;
                    double lineM = 0.0;

                    foreach (b; 1 .. Side)
                    {
                        size_t i = index(a, r, g, b);
                        size_t j = index3(r, g, b);
                        size_t prevR = index3(r - 1, g, b);
                        size_t prevA = index(a - 1, r, g, b);

                        line += weights[i];
                        lineR += momentsR[i];
                        lineG += momentsG[i];
                        lineB += momentsB[i];
                        lineA += momentsA[i];
                        lineM += moments[i];

                        area[b] += line;
                        areaR[b] += lineR;
                        areaG[b] += lineG;
                        areaB[b] += lineB;
                        areaA[b] += lineA;
                        areaM[b] += lineM;

                        xw[j] = xw[prevR] + area[b];
                        xr[j] = xr[prevR] + areaR[b];
                        xg[j] = xg[prevR] + areaG[b];
                        xb[j] = xb[prevR] + areaB[b];
                        xa[j] = xa[prevR] + areaA[b];
                        xm[j] = xm[prevR] + areaM[b];

                        weights[i] = weights[prevA] + xw[j];
                        momentsR[i] = momentsR[prevA] + xr[j];
                        momentsG[i] = momentsG[prevA] + xg[j];
                        momentsB[i] = momentsB[prevA] + xb[j];
                        momentsA[i] = momentsA[prevA] + xa[j];
                        moments[i] = moments[prevA] + xm[j];
                    }
                }
            }
        }
    }

    private double volumeFloat(ref Cube c, const(double)[] m)
        const @nogc nothrow
    {
        int[4] lo = [c.a0, c.r0, c.g0, c.b0];
        int[4] hi = [c.a1, c.r1, c.g1, c.b1];

        double result = 0.0;

        foreach (mask; 0 .. 16)
        {
            int[4] p;
            int bits = 0;

            foreach (d; 0 .. 4)
            {
                if (mask & (1 << d))
                {
                    p[d] = lo[d];
                    ++bits;
                }
                else
                {
                    p[d] = hi[d];
                }
            }

            double v = m[index(p[0], p[1], p[2], p[3])];

            if ((bits & 1) == 0)
                result += v;
            else
                result -= v;
        }

        return result;
    }

    private int volume(ref Cube c, const(double)[] m)
        const @nogc nothrow
    {
        return toInt32(volumeFloat(c, m));
    }

    private double faceSum(
        ref Cube c,
        int direction,
        int position,
        const(double)[] m
    ) const @nogc nothrow
    {
        int[4] lo = [c.a0, c.r0, c.g0, c.b0];
        int[4] hi = [c.a1, c.r1, c.g1, c.b1];

        double result = 0.0;

        foreach (mask; 0 .. 8)
        {
            int[4] p;
            int bits = 0;
            int k = 0;

            foreach (d; 0 .. 4)
            {
                if (d == direction)
                {
                    p[d] = position;
                    continue;
                }

                if (mask & (1 << k))
                {
                    p[d] = lo[d];
                    ++bits;
                }
                else
                {
                    p[d] = hi[d];
                }

                ++k;
            }

            double v = m[index(p[0], p[1], p[2], p[3])];

            if ((bits & 1) == 0)
                result += v;
            else
                result -= v;
        }

        return result;
    }

    private int top(
        ref Cube c,
        int direction,
        int position,
        const(double)[] m
    ) const @nogc nothrow
    {
        return toInt32(faceSum(c, direction, position, m));
    }

    private int bottom(
        ref Cube c,
        int direction,
        const(double)[] m
    ) const @nogc nothrow
    {
        int[4] lo = [c.a0, c.r0, c.g0, c.b0];
        return toInt32(-faceSum(c, direction, lo[direction], m));
    }

    private double variance(ref Cube c) @nogc nothrow
    {
        int w = volume(c, weights);
        if (w <= 0)
            return 0.0;

        double dr = volume(c, momentsR);
        double dg = volume(c, momentsG);
        double db = volume(c, momentsB);
        double da = volume(c, momentsA);
        double m2 = volumeFloat(c, moments);

        double v = m2 -
            (dr * dr + dg * dg + db * db + da * da) / w;

        return v > 0.0 ? v : 0.0;
    }

    private MaxResult maximize(
        ref Cube c,
        int direction,
        int first,
        int last,
        int wholeR,
        int wholeG,
        int wholeB,
        int wholeA,
        int wholeW
    ) @nogc nothrow
    {
        int br = bottom(c, direction, momentsR);
        int bg = bottom(c, direction, momentsG);
        int bb = bottom(c, direction, momentsB);
        int ba = bottom(c, direction, momentsA);
        int bw = bottom(c, direction, weights);

        MaxResult result;
        result.value = 0.0;
        result.position = -1;

        foreach (position; first .. last)
        {
            double hr = br + top(c, direction, position, momentsR);
            double hg = bg + top(c, direction, position, momentsG);
            double hb = bb + top(c, direction, position, momentsB);
            double ha = ba + top(c, direction, position, momentsA);
            double hw = bw + top(c, direction, position, weights);

            if (hw <= 0.0 || hw >= wholeW)
                continue;

            double distance = hr * hr + hg * hg + hb * hb + ha * ha;
            double temp = distance / hw;

            hr = wholeR - hr;
            hg = wholeG - hg;
            hb = wholeB - hb;
            ha = wholeA - ha;
            hw = wholeW - hw;

            if (hw > 0.0)
            {
                distance = hr * hr + hg * hg + hb * hb + ha * ha;
                temp += distance / hw;

                if (temp > result.value)
                {
                    result.value = temp;
                    result.position = position;
                }
            }
        }

        return result;
    }

    private bool cut(ref Cube first, ref Cube second) @nogc nothrow
    {
        int wr = volume(first, weights);
        if (wr <= 1)
            return false;

        int rr = volume(first, momentsR);
        int gg = volume(first, momentsG);
        int bb = volume(first, momentsB);
        int aa = volume(first, momentsA);

        auto red = maximize(
            first, 1, first.r0 + 1, first.r1,
            rr, gg, bb, aa, wr
        );

        auto green = maximize(
            first, 2, first.g0 + 1, first.g1,
            rr, gg, bb, aa, wr
        );

        auto blue = maximize(
            first, 3, first.b0 + 1, first.b1,
            rr, gg, bb, aa, wr
        );

        auto alpha = maximize(
            first, 0, first.a0 + 1, first.a1,
            rr, gg, bb, aa, wr
        );

        int direction;
        int position;

        // image-q tie-breaking: alpha, red, green, blue.
        if (alpha.value >= red.value &&
            alpha.value >= green.value &&
            alpha.value >= blue.value)
        {
            direction = 0;
            position = alpha.position;
        }
        else if (red.value >= green.value &&
                 red.value >= blue.value)
        {
            direction = 1;
            position = red.position;
        }
        else if (green.value >= blue.value)
        {
            direction = 2;
            position = green.position;
        }
        else
        {
            direction = 3;
            position = blue.position;
        }

        if (position < 0)
            return false;

        second = first;

        switch (direction)
        {
            case 0:
                second.a0 = first.a1 = position;
                break;
            case 1:
                second.r0 = first.r1 = position;
                break;
            case 2:
                second.g0 = first.g1 = position;
                break;
            default:
                second.b0 = first.b1 = position;
                break;
        }

        first.volume =
            (first.r1 - first.r0) *
            (first.g1 - first.g0) *
            (first.b1 - first.b0) *
            (first.a1 - first.a0);

        second.volume =
            (second.r1 - second.r0) *
            (second.g1 - second.g0) *
            (second.b1 - second.b0) *
            (second.a1 - second.a0);

        return true;
    }

    private void preparePalette()
    {
        calculateMoments();

        double[] volumeVariance = new double[colors];
        volumeVariance[] = 0.0;
        volumeVariance[0] = variance(cubes[0]);

        uint next = 0;

        for (uint cubeIndex = 1; cubeIndex < colors; ++cubeIndex)
        {
            if (cut(cubes[next], cubes[cubeIndex]))
            {
                volumeVariance[next] =
                    cubes[next].volume > 1
                    ? variance(cubes[next])
                    : 0.0;

                volumeVariance[cubeIndex] =
                    cubes[cubeIndex].volume > 1
                    ? variance(cubes[cubeIndex])
                    : 0.0;
            }
            else
            {
                volumeVariance[next] = 0.0;
                --cubeIndex;
            }

            next = 0;
            double best = volumeVariance[0];

            for (uint i = 1; i <= cubeIndex; ++i)
            {
                if (volumeVariance[i] > best)
                {
                    best = volumeVariance[i];
                    next = i;
                }
            }

            if (best <= 0.0)
            {
                colors = cubeIndex + 1;
                break;
            }
        }
    }

    Color[] quantize(const(ubyte)[] rgba, uint requestedColors)
    {
        if (requestedColors < 1 || requestedColors > 256)
            return null;

        if (rgba.length == 0 || (rgba.length & 3) != 0)
            return null;

        initialize(requestedColors);

        // Histogram.
        foreach (i; 0 .. rgba.length / 4)
        {
            size_t p = i * 4;
            addColor(Color(
                rgba[p],
                rgba[p + 1],
                rgba[p + 2],
                rgba[p + 3]
            ));
        }

        preparePalette();

        int[] lookupR = new int[colors];
        int[] lookupG = new int[colors];
        int[] lookupB = new int[colors];
        int[] lookupA = new int[colors];

        lookupR[] = 0;
        lookupG[] = 0;
        lookupB[] = 0;
        lookupA[] = 0;

        foreach (k; 0 .. colors)
        {
            double w = volume(cubes[k], weights);

            if (w > 0.0)
            {
                lookupR[k] = toInt32(volume(cubes[k], momentsR) / w);
                lookupG[k] = toInt32(volume(cubes[k], momentsG) / w);
                lookupB[k] = toInt32(volume(cubes[k], momentsB) / w);
                lookupA[k] = toInt32(volume(cubes[k], momentsA) / w);
            }
        }

        // Reassign all source pixels to the nearest cube centroid.
        double[] sumsR = new double[colors];
        double[] sumsG = new double[colors];
        double[] sumsB = new double[colors];
        double[] sumsA = new double[colors];
        double[] sumsW = new double[colors];

        sumsR[] = 0.0;
        sumsG[] = 0.0;
        sumsB[] = 0.0;
        sumsA[] = 0.0;
        sumsW[] = 0.0;

        foreach (i; 0 .. rgba.length / 4)
        {
            size_t p = i * 4;

            int r = rgba[p];
            int g = rgba[p + 1];
            int b = rgba[p + 2];
            int a = rgba[p + 3];

            int bestDistance = int.max;
            size_t bestMatch = 0;

            foreach (k; 0 .. colors)
            {
                int d =
                    cast(int)absd(lookupR[k] - r) +
                    cast(int)absd(lookupG[k] - g) +
                    cast(int)absd(lookupB[k] - b) +
                    cast(int)absd(lookupA[k] - a);

                if (d < bestDistance)
                {
                    bestDistance = d;
                    bestMatch = k;
                }
            }

            sumsR[bestMatch] += r;
            sumsG[bestMatch] += g;
            sumsB[bestMatch] += b;
            sumsA[bestMatch] += a;
            sumsW[bestMatch] += 1.0;
        }

        Color[] result;
        foreach (k; 0 .. colors)
        {
            if (sumsW[k] > 0.0)
            {
                result ~= Color(
                    cast(ubyte)toInt32(sumsR[k] / sumsW[k]),
                    cast(ubyte)toInt32(sumsG[k] / sumsW[k]),
                    cast(ubyte)toInt32(sumsB[k] / sumsW[k]),
                    cast(ubyte)toInt32(sumsA[k] / sumsW[k])
                );
            }
        }

        return result;
    }
}
