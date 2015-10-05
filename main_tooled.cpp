#include <cmath>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <random>
#include <set>
#include <vector>
#include <unordered_map>


namespace {

inline size_t GenerateMaskWithRightOnes(size_t onesCount) {
    return (1ull << onesCount) - 1ull;
}

class TCacheSet {
public:
    TCacheSet(size_t waysCount)
        : TicksTimer(0)
        , CacheWays(waysCount)
        , CacheHits(0)
    {}

    bool Touch(size_t cacheLineNum) {
        ++TicksTimer;
        TCacheLine* minTickLine = nullptr;
        size_t minTickTime = TicksTimer;
        for (auto& activeLine : CacheWays) {
            if (activeLine.LineNum != cacheLineNum) {
                if (activeLine.LineNum == -1) {
                    activeLine.TouchTick = TicksTimer;
                    activeLine.LineNum = cacheLineNum;
                    return false;
                }
                if (activeLine.TouchTick < minTickTime) {
                    minTickTime = activeLine.TouchTick;
                    minTickLine = &activeLine;
                }
            } else {
                activeLine.TouchTick = TicksTimer;
                ++CacheHits;
                return true;
            }
        }
        minTickLine->LineNum = cacheLineNum;
        minTickLine->TouchTick = TicksTimer;
        return false;
    }

    size_t GetCacheHits() const {
        return CacheHits;
    }

    size_t GetTicks() const {
        return TicksTimer;
    }

    size_t GetCacheMisses() const {
        return TicksTimer - CacheHits;
    }

    double GetCacheMissRate() const {
        return GetCacheMisses() / (TicksTimer * 1.0);
    }

private:
    struct TCacheLine {
        size_t LineNum = -1;
        size_t TouchTick = 0;
    };

private:
    size_t TicksTimer;
    std::vector<TCacheLine> CacheWays; // Contains LineNum == -1 if free

    size_t CacheHits;
};

class TCache {
public:
    TCache(size_t cacheSize, size_t cacheLineSize, size_t waysCount)
        : CacheSize(cacheSize)
        , CacheLineSize(cacheLineSize)
        , WaysCount(waysCount)
        , BlocksCount(CacheSize / cacheLineSize)
        , SetsCount(BlocksCount / WaysCount)
        , BlockOffsetSize(static_cast<size_t>(log2(CacheLineSize)))
        , ShiftedIndexMask(GenerateMaskWithRightOnes(static_cast<size_t>(log2(SetsCount))))
        , InvBlockMask(~GenerateMaskWithRightOnes(BlockOffsetSize))
        , CacheSets(SetsCount, {WaysCount})
    {
        std::cout << "Params: " << std::endl
                  << "\tCacheSize:\t"        << CacheSize        << ' ' << std::endl
                  << "\tCacheLineSize:\t"    << CacheLineSize    << ' ' << std::endl
                  << "\tWaysCount:\t"        << WaysCount        << ' ' << std::endl
                  << "\tBlocksCount:\t"      << BlocksCount      << ' ' << std::endl
                  << "\tSetsCount:\t"        << SetsCount        << ' ' << std::endl
                  << "\tBlockOffsetSize:\t"  << BlockOffsetSize  << ' ' << std::endl
                  << "\tShiftedIndexMask:\t" << ShiftedIndexMask << ' ' << std::endl
                  << "\tInvBlockMask:\t"     << InvBlockMask     << ' ' << std::endl;
    }

    /// @return true if needed cache line is in cache
    bool Touch(const void* ptr) {
        return CacheSets[GetIndexFromPtr(ptr)].Touch(GetCacheLineNumFromPtr(ptr));
    }

    void PrintStats() const {
        size_t cacheSetIdx = 0;
        double mean = 0.0;
        double meanSq = 0.0;
        size_t totalTicks = 0;
        for (const auto& cacheSet : CacheSets) {
            std::cout << cacheSetIdx << ")\t"
                      << cacheSet.GetCacheMisses() << "\t"
                      << cacheSet.GetTicks() << "\t= "
                      << std::setprecision(3) << cacheSet.GetCacheMissRate() * 100.0 << '%' << std::endl;
            ++cacheSetIdx;
            totalTicks += cacheSet.GetTicks();
            double misses = cacheSet.GetCacheMisses() * 1.0;
            mean += misses / CacheSets.size();
            meanSq += misses * misses / CacheSets.size();
        }
        std::cout << std::fixed << std::setprecision(3)
                  << "Misses:" << std::endl
                  << "\tmean=" << mean << std::endl
                  << "\tvar=" << (meanSq - mean * mean) << std::endl
                  << "\trate=" << (mean * CacheSets.size() / totalTicks) * 100 << '%'
                  << std::endl;

    }

private:
    size_t GetIndexFromPtr(const void* ptr) {
        return ((size_t)ptr >> BlockOffsetSize) & ShiftedIndexMask;
    }

    size_t GetCacheLineNumFromPtr(const void* ptr) {
        return (size_t)ptr & InvBlockMask;
    }

private:
    const size_t CacheSize;
    const size_t CacheLineSize;
    const size_t WaysCount;
    const size_t BlocksCount;
    const size_t SetsCount;
    const size_t BlockOffsetSize;

    const size_t ShiftedIndexMask;
    const size_t InvBlockMask;

    std::vector<TCacheSet> CacheSets;
};

void MultSimpleTooled(TCache& cache,
                      const float* __restrict a,
                      const float* __restrict b,
                      float* __restrict c,
                      int n)
{
    /// We suppose that i j k, as iteration vars, would be placed to registers
    /// So cache will work only with n, a, b, c and corresponding arrays
    cache.Touch(&n); // for next loop start
    for (int i = 0; i < n; ++i) {
        cache.Touch(&n); // for next loop start
        for (int j = 0; j < n; ++j) {
            //>-- c[i * n + j]
            // cache.Touch(&i);
            // cache.Touch(&j);
            cache.Touch(&n);
            cache.Touch(&c);
            cache.Touch(&c[i * n + j]);
            //<-- c[i * n + j]

            c[i * n + j] = 0.f;

            cache.Touch(&n); // for next loop start
            for (int k = 0; k < n; ++k) {
                // cache.Touch(&i);
                // cache.Touch(&j);
                // cache.Touch(&k);
                cache.Touch(&n);
                cache.Touch(&a);
                cache.Touch(&b);
                cache.Touch(&c);
                cache.Touch(&a[i * n + k]);
                cache.Touch(&b[k * n + j]);
                cache.Touch(&c[i * n + j]);

                c[i * n + j] += a[i * n + k] * b[k * n + j];

                //>-- k < n
                // cache.Touch(&k);
                cache.Touch(&n);
                //<-- k < n
            }

            //>-- j < n
            // cache.Touch(&j);
            cache.Touch(&n);
            //<-- j < n
        }

        //>-- i < n
        // cache.Touch(&i);
        cache.Touch(&n);
        //<-- i < n
    }
}

void MultSimple(const float* __restrict a, const float* __restrict b, float* __restrict c, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            c[i * n + j] = 0.f;
            for (int k = 0; k < n; ++k) {
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
}

void FillRandom(float* a, int n)
{
    std::default_random_engine eng;
    std::uniform_real_distribution<float> dist;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i * n + j] = dist(eng);
        }
    }
}

}


int main(int argc, char* argv[])
{
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " MATRIX_SIZE CACHE_SIZE_KB CACHE_LINE_SIZE_B WAYS_COUNT" << std::endl;
        return 1;
    }
    const int n = atoi(argv[1]);
    const size_t cacheSizeKB    = atoi(argv[2]);
    const size_t cacheLineSizeB = atoi(argv[3]);
    const size_t waysCount      = atoi(argv[4]);
    std::cout << "n = " << n << std::endl;
    std::cout << "Cache: " << cacheSizeKB << "KB, "
              << waysCount << " ways with "
              << cacheLineSizeB << "b cache line." << std::endl;
    float* a = new float[n * n];
    float* b = new float[n * n];
    float* c = new float[n * n];

    FillRandom(a, n);
    FillRandom(b, n);
    TCache cache(cacheSizeKB * 1024 /* kb */, cacheLineSizeB /* b */, waysCount /* ways */);
    {
        const auto startTime = std::clock();
        MultSimpleTooled(cache, a, b, c, n);
        const auto endTime = std::clock();

        std::cout << "timeSimple: " << double(endTime - startTime) / CLOCKS_PER_SEC << std::endl;
    }
    cache.PrintStats();

    delete[] a;
    delete[] b;
    delete[] c;
}

