#include <cmath>
#include <ctime>
#include <iostream>
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
    {}

    bool Touch(size_t cacheLineNum) {
        ++TicksTimer;
        // CacheLine2LastTouch[cacheLineNum] = TicksTimer;
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

    size_t GetCacheMisses() const {
        return TicksTimer - CacheHits;
    }

private:
    struct TCacheLine {
        size_t LineNum = -1;
        size_t TouchTick = 0;
    };

private:
    size_t TicksTimer;
    std::vector<TCacheLine> CacheWays; // Contains LineNum == -1 if free
    // std::unordered_map<size_t, size_t> CacheLine2LastTouch;

    size_t CacheHits;
};

class TCache {
public:
    TCache(size_t cacheSize, size_t cacheLineSize, size_t waysCount)
        : CacheSize(cacheSize)
        , CacheLineSize(cacheLineSize)
        , WaysCount(waysCount)
        , BlocksCount(CacheSize / CacheLineSize)
        , SetsCount(BlocksCount / WaysCount)
        , BlockOffsetSize(static_cast<size_t>(log2(CacheLineSize)))
        , ShiftedIndexMask(GenerateMaskWithRightOnes(static_cast<size_t>(log2(SetsCount))))
        , InvBlockMask(~GenerateMaskWithRightOnes(BlockOffsetSize))
        , CacheSets(SetsCount, {WaysCount})
    {}

    /// @return true if needed cache line is in cache
    bool Touch(const void* ptr) {
        // std::cerr << std::oct << GetIndexFromPtr(ptr) << "\t->\t" << GetCacheLineNumFromPtr(ptr) << std::endl;
        return CacheSets[GetIndexFromPtr(ptr)].Touch(GetCacheLineNumFromPtr(ptr));
    }

    std::vector<size_t> GetCacheMisses() const {
        std::vector<size_t> result;
        result.reserve(SetsCount);
        for (const auto& cacheSet : CacheSets) {
            result.push_back(cacheSet.GetCacheMisses());
        }
        return result;
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

void PrintStats(const TCache& cache) {
    auto misses = cache.GetCacheMisses();
    // auto hits = cache.GetCacheHits();
    for (size_t i = 0; i < misses.size(); ++i) {
        std::cerr << i << ")\t" << misses[i] << std::endl;
    }
}

void MultSimpleTooled(TCache& cache,
                      const float* __restrict a,
                      const float* __restrict b,
                      float* __restrict c,
                      int n)
{
    for (int i = 0; i < n; ++i) {
        if (i == 0) {
            //>-- i = 0
            cache.Touch(&i);
            //<-- i = 0
            //>-- i < n
            cache.Touch(&i);
            cache.Touch(&n);
            //<-- i < n
        }

        for (int j = 0; j < n; ++j) {
            if (j == 0) {
                //>-- j = 0
                cache.Touch(&j);
                //<-- j = 0
                //>-- j < n
                cache.Touch(&j);
                cache.Touch(&n);
                //<-- j < n
            }

            //>-- c[i * n + j]
            cache.Touch(&i);
            cache.Touch(&j);
            cache.Touch(&n);
            cache.Touch(&c);
            cache.Touch(&c[i * n + j]);
            //<-- c[i * n + j]

            c[i * n + j] = 0.f;
            for (int k = 0; k < n; ++k) {
                if (k == 0) {
                    //>-- k = 0
                    cache.Touch(&k);
                    //<-- k = 0
                    //>-- k < n
                    cache.Touch(&k);
                    cache.Touch(&n);
                    //<-- k < n
                }

                cache.Touch(&i);
                cache.Touch(&j);
                cache.Touch(&k);
                cache.Touch(&n);
                cache.Touch(&a);
                cache.Touch(&b);
                cache.Touch(&c);
                cache.Touch(&a[i * n + k]);
                cache.Touch(&b[k * n + j]);
                cache.Touch(&c[i * n + j]);

                c[i * n + j] += a[i * n + k] * b[k * n + j];

                //>-- ++k
                cache.Touch(&k);
                //<-- ++k
                //>-- k < n
                cache.Touch(&k);
                cache.Touch(&n);
                //<-- k < n
            }

            //>-- ++j
            cache.Touch(&j);
            //<-- ++j
            //>-- j < n
            cache.Touch(&j);
            cache.Touch(&n);
            //<-- j < n
        }

        //>-- ++i
        cache.Touch(&i);
        //<-- ++i
        //>-- i < n
        cache.Touch(&i);
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
    const int n = atoi(argv[1]);
    std::cerr << "n = " << n << std::endl;

    float* a = new float[n * n];
    float* b = new float[n * n];
    float* c = new float[n * n];

    FillRandom(a, n);
    FillRandom(b, n);

    TCache cache(3*1024*1024 /* 32kb */, 64 /* 64b */, 12 /* ways */);
    {
        const auto startTime = std::clock();
        MultSimpleTooled(cache, a, b, c, n);
        const auto endTime = std::clock();

        std::cerr << "timeSimple: " << double(endTime - startTime) / CLOCKS_PER_SEC << '\n';
    }
    PrintStats(cache);

    delete[] a;
    delete[] b;
    delete[] c;
}

