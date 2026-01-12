#ifndef __UTILS_SKETCHES_H__
#define __UTİLS_SKETCHES_H__

#include <iomanip>
#include <sstream>
#include <vector>
#include <math.h>
#include <stdint.h>

#define DISTANCE_LENGTH 65536

#define min(x, y) (x < y ? x : y)
#define max(x, y) (x > y ? x : y)

typedef struct {
	uint64_t minimum;
	uint64_t maximum;
} gap;

std::string format_int(int value) {
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << value;
    return ss.str();
};

std::string format_double(double value, size_t precision = 2) {
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
};

double mean(int (&numbers)[DISTANCE_LENGTH], std::vector<uint32_t> numbersXL = {}) {
    double sum = 0;
    double count = 0;
    for (size_t i = 0; i < DISTANCE_LENGTH; i++) {
        sum += (i * numbers[i]);
        count += numbers[i];
    }
    for (size_t i = 0; i < numbersXL.size(); i++) {
        sum += numbersXL[i];
    }
    count += numbersXL.size();
    return sum / count;
};

double stdev(int (&numbers)[DISTANCE_LENGTH], std::vector<uint32_t> numbersXL = {}) {
    double mean_value = mean(numbers, numbersXL);
    double count = 0;
    for (size_t i = 0; i < DISTANCE_LENGTH; i++) {
        count += numbers[i];
    }
    count += numbersXL.size();
    double variance = 0;
    for (size_t i = 0; i < DISTANCE_LENGTH; i++) {
        variance += ((mean_value - i) * (mean_value - i) * numbers[i]);
    }
    for (size_t i = 0; i < numbersXL.size(); i++) {
        variance += ((mean_value - numbersXL[i]) * (mean_value - numbersXL[i]));
    }
    return sqrt(variance / count);
};


#endif