#include <iostream>
#include <random>
#include <ctime>


//__declspec(align(16))
// __restrict
namespace {
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

    {
        const auto startTime = std::clock();
        MultSimple(a, b, c, n);
        const auto endTime = std::clock();

        std::cerr << "timeSimple: " << double(endTime - startTime) / CLOCKS_PER_SEC << '\n';
    }

    delete[] a;
    delete[] b;
    delete[] c;
}

