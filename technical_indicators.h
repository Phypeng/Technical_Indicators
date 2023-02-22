#ifndef TECHNICAL_INDICATORS_H
#define TECHNICAL_INDICATORS_H

#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>

using namespace std;


template<typename T>
inline vector<T> RangeCopy(vector<T> v, int pos, int length = -1) {
    if (v.empty() || pos > int(v.size() - 1)) {
        return vector<T>();
    }
    if (pos < 0) {
        pos = 0;
    }
    if (length < 0 || pos + length > int(v.size())) {
        length = int(v.size()) - pos;
    }
    return vector<T>(v.begin() + pos, v.begin() + pos + length);
}

//四则运算方法

template <class T>
inline auto ti_IF(const vector<bool> &cond, const vector<T> &a, const vector<T> &b) {
    std::size_t cSize = cond.size();
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(min(aSize, bSize), cSize);
    vector<T> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = cond[i + cSize - minSize] ? a[i + aSize - minSize] : b[i + bSize - minSize];
    }
    return res;
}

template <class T>
inline auto ti_IF(const vector<bool> &cond, const vector<T> &a, const T &b) {
    std::size_t cSize = cond.size();
    std::size_t aSize = a.size();
    std::size_t minSize = min(aSize, cSize);
    vector<T> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = cond[i + cSize - minSize] ? a[i + aSize - minSize] : b;
    }
    return res;
}

template <class T>
inline auto ti_IF(const vector<bool> &cond, const T &a, const vector<T> &b) {
    std::size_t cSize = cond.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(bSize, cSize);
    vector<T> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = cond[i + cSize - minSize] ? a : b[i + bSize - minSize];
    }
    return res;
}

template <class T>
inline auto ti_IF(const vector<bool> &cond, const T &a, const T &b) {
    std::size_t cSize = cond.size();
    vector<T> res(cSize);

    for (std::size_t i = 0; i < cSize; i++) {
        res[i] = cond[i] ? a : b;
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Sub(const vector<T1> &a, const vector<T2> &b) {
    T1 t1; T2 t2;
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<decltype(t1 - t2)> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = a[i + aSize - minSize] - b[i + bSize - minSize];
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Sub(const vector<T1> &a, const T2 &b) {
    T1 t1; T2 t2;
    std::size_t aSize = a.size();
    vector<decltype(t1 - t2)> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = a[i] - b;
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Sub(const T2 &a, const vector<T1> &b) {
    T1 t1; T2 t2;
    std::size_t Size = b.size();
    vector<decltype(t1 - t2)> res(Size);

    for (std::size_t i = 0; i < Size; i++) {
        res[i] = a - b[i];
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Add(const vector<T1> &a, const vector<T2> &b) {
    T1 t1; T2 t2;
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<decltype(t1 + t2)> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = a[i + aSize - minSize] + b[i + bSize - minSize];
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Add(const vector<T1> &a, const T2 &b) {
    T1 t1; T2 t2;
    std::size_t aSize = a.size();
    vector<decltype(t1 - t2)> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = a[i] + b;
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Add(const T2 &a, const vector<T1> &b) {
    return ti_Add(b, a);
}


template <class T1, class T2>
inline auto ti_Mul(const vector<T1> &a, const vector<T2> &b) {
    T1 t1; T2 t2;
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<decltype(t1 * t2)> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = a[i + aSize - minSize] * b[i + bSize - minSize];
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Mul(const vector<T1> &a, const T2 &b) {
    T1 t1; T2 t2;
    std::size_t aSize = a.size();
    vector<decltype(t1 * t2)> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = a[i] * b;
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Mul(const T2 &a, const vector<T1> &b) {
    return ti_Mul(b, a);
}

template <class T1, class T2>
inline auto ti_Div(const vector<T1> &a, const vector<T2> &b) {
    T1 t1 = 1; T2 t2 = 1;
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<decltype(t1 / t2)> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        if (std::abs(double(b[i + bSize - minSize])) > 0.000001) {
            res[i] = a[i + aSize - minSize] / b[i + bSize - minSize];
        }
        else {
            res[i] = std::numeric_limits<decltype(t1 / t2)>::quiet_NaN();
        }
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Div(const vector<T1> &a, const T2 &b) {
    T1 t1 = 1; T2 t2 = 1;
    std::size_t aSize = a.size();
    if (abs(b) > 0.000001) {
        vector<decltype(t1 / t2)> res(aSize);

        for (std::size_t i = 0; i < aSize; i++) {
            res[i] = a[i] / b;
        }
        return res;
    }
    else {
        return vector<decltype(t1 / t2)>(aSize, std::numeric_limits<decltype(t1 / t2)>::quiet_NaN());
    }
}

template <class T1, class T2>
inline auto ti_Div(const T2 &a, const vector<T1> &b) {
    T1 t1 = 1; T2 t2 = 1;
    std::size_t Size = b.size();
    vector<decltype(t1 / t2)> res(Size);

    for (std::size_t i = 0; i < Size; i++) {
        if (abs(b[i]) > 0.000001) {
            res[i] = a / b[i];
        }
        else {
            res[i] = std::numeric_limits<decltype(t1 / t2)>::quiet_NaN();
        }
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_AndAnd(const vector<T1> &a, const vector<T2> &b) {
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<bool> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = a[i + aSize - minSize] && b[i + bSize - minSize];
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_AndAnd(const vector<T1> &a, const T2 &b) {
    std::size_t aSize = a.size();
    vector<bool> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = a[i] && b;
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_AndAnd(const T2 &a, const vector<T1> &b) {
    return ti_AndAnd(b, a);
}


template <class T1, class T2>
inline vector<bool> ti_OrOr(const vector<T1> &a, const vector<T2> &b) {
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<bool> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = a[i + aSize - minSize] || b[i + bSize - minSize];
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_OrOr(const vector<T1> &a, const T2 &b) {
    std::size_t aSize = a.size();
    vector<bool> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = a[i] || b;
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_OrOr(const T2 &a, const vector<T1> &b) {
    return ti_OrOr(b, a);
}


inline vector<bool> ti_Not(const vector<bool> &a) {
    std::size_t aSize = a.size();
    vector<bool> res(aSize);
    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = !a[i];
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_EqualEqual(const vector<T1> &a, const vector<T2> &b, double deviation = 0.00001) {
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<bool> res(minSize);

    if (std::is_floating_point<T1>::value || std::is_floating_point<T2>::value) {
        for (std::size_t i = 0; i < minSize; i++) {
            res[i] = abs(a[i + aSize - minSize] - b[i + bSize - minSize]) <= deviation;
        }
    }
    else {
        for (std::size_t i = 0; i < minSize; i++) {
            res[i] = a[i + aSize - minSize] == b[i + bSize - minSize];
        }
    }

    return res;
}

template <class T1, class T2>
inline vector<bool> ti_EqualEqual(const vector<T1> &a, const T2 &b, double deviation = 0.00001) {
    std::size_t aSize = a.size();
    vector<bool> res(aSize);

    if (std::is_floating_point<T1>::value || std::is_floating_point<T2>::value) {
        for (std::size_t i = 0; i < aSize; i++) {
            res[i] = abs(a[i] - b) <= deviation;
        }
    }
    else {
        for (std::size_t i = 0; i < aSize; i++) {
            res[i] = a[i] == b;
        }
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_EqualEqual(const T2 &a, const vector<T1> &b, double deviation = 0.00001) {
    return ti_EqualEqual(b, a, deviation);
}


template <class T1, class T2>
inline vector<bool> ti_NotEqual(const vector<T1> &a, const vector<T2> &b, double deviation = 0.00001) {
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<bool> res(minSize);

    if (std::is_floating_point<T1>::value || std::is_floating_point<T2>::value) {
        for (std::size_t i = 0; i < minSize; i++) {
            res[i] = abs(a[i + aSize - minSize] - b[i + bSize - minSize]) > deviation;
        }
    }
    else {
        for (std::size_t i = 0; i < minSize; i++) {
            res[i] = a[i + aSize - minSize] != b[i + bSize - minSize];
        }
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_NotEqual(const vector<T1> &a, const T2 &b, double deviation = 0.00001) {
    std::size_t aSize = a.size();
    vector<bool> res(aSize);

    if (std::is_floating_point<T1>::value || std::is_floating_point<T2>::value) {
        for (std::size_t i = 0; i < aSize; i++) {
            res[i] = abs(double(a[i] - b)) > deviation;
        }
    }
    else {
        for (std::size_t i = 0; i < aSize; i++) {
            res[i] = a[i] != b;
        }
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_NotEqual(const T2 &a, const vector<T1> &b, double deviation = 0.00001) {
    return ti_NotEqual(b, a, deviation);
}


template <class T1, class T2>
inline vector<bool> ti_LessThan(const vector<T1> &a, const vector<T2> &b) {
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<bool> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = a[i + aSize - minSize] < b[i + bSize - minSize];
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_LessThan(const vector<T1> &a, const T2 &b) {
    std::size_t aSize = a.size();
    vector<bool> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = a[i] < b;
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_LessThan(const T2 &a, const vector<T1> &b) {
    std::size_t Size = b.size();
    vector<bool> res(Size);

    for (std::size_t i = 0; i < Size; i++) {
        res[i] = a < b[i];
    }
    return res;
}


template <class T1, class T2>
inline vector<bool> ti_MoreThan(const vector<T1> &a, const vector<T2> &b) {
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<bool> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = a[i + aSize - minSize] > b[i + bSize - minSize];
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_MoreThan(const vector<T1> &a, const T2 &b) {
    std::size_t aSize = a.size();
    vector<bool> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = a[i] > b;
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_MoreThan(const T2 &a, const vector<T1> &b) {
    std::size_t Size = b.size();
    vector<bool> res(Size);

    for (std::size_t i = 0; i < Size; i++) {
        res[i] = a > b[i];
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_LessEqual(const vector<T1> &a, const vector<T2> &b) {
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<bool> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = a[i + aSize - minSize] <= b[i + bSize - minSize];
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_LessEqual(const vector<T1> &a, const T2 &b) {
    std::size_t aSize = a.size();
    vector<bool> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = a[i] <= b;
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_LessEqual(const T2 &a, const vector<T1> &b) {
    std::size_t Size = b.size();
    vector<bool> res(Size);

    for (std::size_t i = 0; i < Size; i++) {
        res[i] = a <= b[i];
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_MoreEqual(const vector<T1> &a, const vector<T2> &b) {
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<bool> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = a[i + aSize - minSize] >= b[i + bSize - minSize];
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_MoreEqual(const vector<T1> &a, const T2 &b) {
    std::size_t aSize = a.size();
    vector<bool> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = a[i] >= b;
    }
    return res;
}

template <class T1, class T2>
inline vector<bool> ti_MoreEqual(const T2 &a, const vector<T1> &b) {
    std::size_t Size = b.size();
    vector<bool> res(Size);

    for (std::size_t i = 0; i < Size; i++) {
        res[i] = a >= b[i];
    }
    return res;
}

template <class T>
inline vector<T> ti_Abs(const vector<T> &a) {
    std::size_t aSize = a.size();
    vector<T> res(aSize);
    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = abs(a[i]);
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Max(const vector<T1> &a, const vector<T2> &b) {
    T1 t1; T2 t2;
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<decltype(t1 + t2)> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = max(a[i + aSize - minSize], b[i + bSize - minSize]);
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Max(const vector<T1> &a, const T2 &b) {
    T1 t1; T2 t2;
    std::size_t aSize = a.size();
    auto bt = decltype(t1 + t2)(b);
    vector<decltype(t1 + t2)> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = max(a[i], bt);
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Max(const T2 &a, const vector<T1> &b) {
    return ti_Max(b, a);
}


template <class T1, class T2>
inline auto ti_Min(const vector<T1> &a, const vector<T2> &b) {
    T1 t1; T2 t2;
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t minSize = min(aSize, bSize);
    vector<decltype(t1 + t2)> res(minSize);

    for (std::size_t i = 0; i < minSize; i++) {
        res[i] = min(a[i + aSize - minSize], b[i + bSize - minSize]);
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Min(const vector<T1> &a, const T2 &b) {
    T1 t1; T2 t2;
    std::size_t aSize = a.size();
    vector<decltype(t1 + t2)> res(aSize);

    for (std::size_t i = 0; i < aSize; i++) {
        res[i] = min(a[i], b);
    }
    return res;
}

template <class T1, class T2>
inline auto ti_Min(const T2 &a, const vector<T1> &b) {
    return ti_Min(b, a);
}


struct ti_Point {
    double x, y;
};

inline double ti_Direction(ti_Point p1, ti_Point p2, ti_Point p0) {
    return (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y);
}

inline bool ti_OnSegment(ti_Point p1, ti_Point p2, ti_Point p0) {
    return (min(p1.x, p2.x) <= p0.x && p0.x <= max(p1.x, p2.x) &&
        min(p1.y, p2.y) <= p0.y && p0.y <= max(p1.y, p2.y));
}

inline bool ti_SegmentsIntersert(ti_Point p1, ti_Point p2, ti_Point p3, ti_Point p4) {
    double d1 = ti_Direction(p3, p4, p1);
    double d2 = ti_Direction(p3, p4, p2);
    double d3 = ti_Direction(p1, p2, p3);
    double d4 = ti_Direction(p1, p2, p4);

    if (d1 * d2 < 0 && d3 * d4 < 0) {
        return true;
    }
    else if (d1 == 0.0 && ti_OnSegment(p3, p4, p1)) {
        return true;
    }
    else if (d2 == 0.0 && ti_OnSegment(p3, p4, p2)) {
        return true;
    }
    else if (d3 == 0.0 && ti_OnSegment(p1, p2, p3)) {
        return true;
    }
    else if (d4 == 0.0 && ti_OnSegment(p1, p2, p4)) {
        return false;
    }
    return false;
}

inline vector<double> ti_DRAWLINE(const vector<bool> &cond1, const vector<double> &data1, const vector<bool> &cond2, const vector<double> &data2) {
    std::size_t cond1Size = cond1.size();
    std::size_t data1Size = data1.size();
    std::size_t cond2Size = cond2.size();
    std::size_t data2Size = data2.size();
    std::size_t minSize = min(min(min(cond1Size, cond2Size), data1Size), data2Size);

    double nan = numeric_limits<double>::quiet_NaN();
    bool startFlag = false;
    std::size_t lastStartIndex = 0, lastEndIndex = 0;
    vector<double> res(minSize, nan);
    vector<int> flag(minSize, -1);//-1:空位置 0:既是起点又是终点 1:起点  2:终点

    for (std::size_t i = 0; i < minSize; i++) {
        if (cond1[i + cond1Size - minSize]) {
            if (!startFlag) { startFlag = true; }
            if (lastStartIndex > lastEndIndex) {
                res[lastStartIndex] = nan;
                flag[lastStartIndex] = -1;
            }
            lastStartIndex = i;
            res[i] = data1[i + data1Size - minSize];
            flag[i] = 1;
        }

        if (cond2[i + cond2Size - minSize] && startFlag) {
            if (lastStartIndex < lastEndIndex) {
                res[lastEndIndex] = nan;
                flag[lastEndIndex] = -1;
            }
            if (flag[lastStartIndex] == 0) {
                res[lastStartIndex] = data1[lastStartIndex + data1Size - minSize];
                flag[lastStartIndex] = 1;
            }
            lastEndIndex = i;
            res[i] = data2[i + data2Size - minSize];

            if (flag[i] == 1) {
                flag[i] = 0;
            }
            else {
                flag[i] = 2;
            }
        }
    }

    //填充
    double k, b;
    std::size_t startIndex = 0, endIndex = 0;
    for (std::size_t i = 0; i < minSize; i++) {
        if (flag[i] == 1) {
            startIndex = i;
            for (std::size_t n = i + 1; n < minSize; n++) {
                if (flag[n] == 2) {
                    endIndex = n;
                    if (endIndex - startIndex > 1) {
                        k = (res[endIndex] - res[startIndex]) / (endIndex - startIndex);
                        b = res[startIndex] - k * startIndex;
                        for (std::size_t t = startIndex + 1; t < endIndex; t++) {
                            res[t] = k * t + b;
                        }
                    }
                    i = n;
                    break;
                }
            }
        }
    }

    return res;
}

inline vector<double> ti_MA(const vector<double> &data, std::size_t period) {
    std::size_t dataCount = data.size();
    if (dataCount < period) {
        return vector<double>();
    }

    vector<double> res(dataCount - period + 1);

    double temp = 0;
    for (std::size_t i = period - 1; i < dataCount; i++) {
        temp = 0;
        for (std::size_t j = 0; j < period; j++) {
            temp += data[i - j];
        }
        res[i - period + 1] = temp / period;
    }
    return res;
}

inline vector<double> ti_SUM(const vector<double> &data, std::size_t period) {
    std::size_t size = data.size();
    vector<double> res(data);

    for (std::size_t i = 1; i < size; i++) {
        for (std::size_t j = 1; j < (period == 0 ? i + 1 : period); j++) {
            if (i >= j) {
                res[i] += data[i - j];
            }
        }
    }

    return res;
}

inline vector<bool> ti_CROSS(const vector<double> &a, const vector<double> &b) {
    std::size_t aSize = a.size();
    std::size_t bSize = b.size();
    std::size_t size = min(aSize, bSize);

    vector<bool> res(size, false);
    for (std::size_t i = 1; i < size; i++) {
        res[i] = a[i + aSize - size] > b[i + bSize - size] && a[i + aSize - size - 1] < b[i + bSize - size - 1];
    }
    return res;
}

inline vector<bool> ti_CROSS(const vector<double> &a, double b) {
    std::size_t size = a.size();
    vector<bool> res(size, false);
    for (std::size_t i = 1; i < size; i++) {
        res[i] = a[i] > b && a[i - 1] < b;
    }
    return res;
}

inline vector<std::size_t> ti_COUNT(const vector<bool> &a, const vector<std::size_t> &period) {
    std::size_t aSize = a.size();
    std::size_t periodSize = period.size();
    std::size_t size = min(aSize, periodSize);

    vector<std::size_t> res(size, false);
    for (std::size_t i = 1; i < size; i++) {
        for (std::size_t j = 0; j < period[i]; j++) {
            if (i >= j && a[i - j]) {
                res[i]++;
            }
        }
    }
    return res;
}

inline vector<std::size_t> ti_COUNT(const vector<bool> &a, std::size_t period) {
    std::size_t size = a.size();
    vector<std::size_t> res(size, 0);
    for (std::size_t i = 1; i < size; i++) {
        for (std::size_t j = 0; j < period; j++) {
            if (i >= j && a[i - j]) {
                res[i]++;
            }
        }
    }
    return res;
}



template <class T>
inline vector<T> ti_REF(const vector<T> &data, std::size_t period) {
    if (period > data.size()) {
        return vector<T>();
    }
    return RangeCopy(data, 0, int(data.size() - period));
}

template <class T>
inline vector<T> ti_REF(const vector<T> &data, const vector<std::size_t> &period) {
    std::size_t size_d = data.size();
    std::size_t size_p = period.size();
    std::size_t size_min = min(size_d,size_p);
    vector<T> res(size_min);
    for (std::size_t i = 0; i < size_min; i++) {
        if (i >= period[i + size_p - size_min]) {
            res[i] = data[i + size_d - size_min - period[i]];
        }
    }
    return res;
}

inline vector<double> ti_LLV(const vector<double> &data, const vector<std::size_t> &period) {
    std::size_t dsize = data.size();
    std::size_t psize = period.size();
    std::size_t minSize = min(dsize, psize);
    vector<double> res = RangeCopy(data, int(dsize - minSize));

    for (std::size_t i = 0; i < minSize; i++) {
        for (std::size_t j = 0; j < period[i + psize - minSize]; j++) {
            if (i >= j) {
                res[i] = min(res[i], data[i - j + dsize - minSize]);
            }
        }
    }
    return res;
}

inline vector<double> ti_HHV(const vector<double> &data, const vector<std::size_t> &period) {
    std::size_t dsize = data.size();
    std::size_t psize = period.size();
    std::size_t minSize = min(dsize, psize);
    vector<double> res = RangeCopy(data, int(dsize - minSize));

    for (std::size_t i = 0; i < minSize; i++) {
        for (std::size_t j = 0; j < period[i + psize - minSize]; j++) {
            if (i >= j) {
                res[i] = max(res[i], data[i - j + dsize - minSize]);
            }
        }
    }
    return res;
}

inline vector<double> ti_LLV(const vector<double> &data, std::size_t period) {
    std::size_t size = data.size();
    vector<double> res(data);

    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = 0; j < period; j++) {
            if (i >= j) {
                if (res[i] > data[i - j]) {
                    res[i] = data[i - j];
                }
            }
        }
    }
    return res;
}

inline vector<double> ti_HHV(const vector<double> &data, std::size_t period)
{
    std::size_t size = data.size();
    vector<double> res(data);

    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = 0; j < period; j++) {
            if (i >= j) {
                if (res[i] < data[i - j]) {
                    res[i] = data[i - j];
                }
            }
        }
    }
    return res;
}

inline vector<double> ti_SMA(const vector<double> &data, std::size_t period, double weight) {
    if (data.empty()) {
        return vector<double>();
    }
    std::size_t count = data.size();

    std::size_t startIndex = 0;
    for (std::size_t i = 0; i < count; i++) {
        if (!std::isnan(data[i])) {
            startIndex = i;
            break;
        }
    }
    vector<double> res = vector<double>(count - startIndex);
    if (res.empty()) {
        return res;
    }
    double alpha = weight / period;
    double temp = data.at(startIndex);
    double perValue = temp;
    res[0] = temp;
    for (std::size_t i = startIndex + 1; i < count - startIndex; i++) {
        if (!std::isnan(data[i + startIndex])) {
            perValue = data[i + startIndex];
        }
        temp = temp * (1 - alpha) + perValue * alpha;
        res[i] = temp;
    }
    return res;
}

inline vector<double> ti_EMA(const vector<double> &data, std::size_t period) {
    return ti_SMA(data, period + 1, 2);
}

inline vector<double> ti_AVEDEV(const vector<double> &data, std::size_t period) {
    vector<double> mean = ti_MA(data, period);
    std::size_t meanSize = mean.size();
    vector<double> res(meanSize);
    double sum = 0;

    for (std::size_t i = 0; i < meanSize; i++) {
        sum = 0;
        for (std::size_t j = 0; j < period; j++) {
            sum += abs(data[i + j] - mean[i]);
        }
        res[i] = sum / period;
    }
    return res;
}

inline vector<double> ti_STDDEV(const vector<double> &close, std::size_t period) {
    vector<double> mean = ti_MA(close, period);
    std::size_t meanCount = mean.size();
    vector<double> res = vector<double>(meanCount);
    double sum = 0;

    for (std::size_t i = 0; i < meanCount; i++) {
        sum = 0;
        for (std::size_t j = 0; j < period; j++) {
            sum += pow(close.at(i + j) - mean.at(i), 2);
        }
        res[i] = sqrt(sum / (period - 1));
    }
    return res;
}

inline vector<bool> ti_BACKSET(vector<bool> &data, std::size_t period) {
    vector<bool> res = data;
    std::size_t count = res.size();
    for (std::size_t i = 0; i < count; i++) {
        if (res.at(i)) {
            std::size_t dif = i >= period - 1 ? i - period + 1 : i;
            for (std::size_t c = dif; c < i; c++) {
                res[c] = true;
            }
        }
    }
    return res;
}

inline vector<bool> ti_FILTER(const vector<bool> &data, std::size_t period) {
    vector<bool> res = data;
    std::size_t count = res.size();
    for (std::size_t i = 0; i < count; i++) {
        if (res.at(i)) {
            std::size_t dif = i + period < count ? period : count - 1 - i;
            for (std::size_t c = i + 1; c < dif + i + 1; c++) {
                res[c] = false;
            }
            i += dif;
        }
    }
    return res;
}

inline vector<std::size_t> ti_BARSLAST(const vector<bool> &data) {
    std::size_t count = data.size();
    vector<std::size_t> BARSLAST(count, 0);
    std::size_t index = 0;
    bool flag = false;
    for (std::size_t i = 0; i < count; i++) {
        if (data[i]) {
            index = i;
            flag = true;
        }
        else if (flag) {
            BARSLAST[i] = i - index;
        }
    }
    return BARSLAST;
}

inline vector<std::size_t> ti_BARSCOUNT(const vector<double> &data){
    auto size = data.size();
    vector<std::size_t> res(size);
    for(std::size_t i = 0;i < size;i++){
        res[i] = i;
    }
    return res;
}

inline vector<double> ti_DMA(const vector<double> &data,const vector<double> &period){
    auto dsize = data.size();
    auto psize = period.size();
    auto minSize = min(dsize, psize);
    vector<double> res = RangeCopy(data,int(dsize - minSize));

    for (std::size_t i = 1; i < minSize; i++) {
        double p = std::min(std::max(period[i + psize - minSize],.0),1.);
        res[i] = res[i] * p + (1 - p) * res[i - 1];
    }
    return res;
}

inline vector<double> ti_DMA(const vector<double> &data,double period){
    vector<double> res = data;
    double p = std::min(std::max(period,.0),1.);
    auto size = data.size();
    for (std::size_t i = 1; i < size; i++) {
        res[i] = res[i] * p + (1 - p) * res[i - 1];
    }
    return res;
}

template <class T>
vector<T> ti_SUM(const vector<T> &data,unsigned int period){
    T tmp = 0;
    auto size = data.size();
    vector<T> res(size);


    for(unsigned int i = 0;i < size;i++){
        tmp += data[i];
        if(period != 0 && i >= period){
            tmp -= data[i - period];
        }
        res[i] = tmp;
    }

    return res;
}

inline vector<double> ti_SAR(const vector<double> &close, const vector<double> &low, const vector<double> &high, std::size_t period, double af = 0.02, double incre = 0.02, double critical = 0.2) {

    std::size_t size = close.size();

    if (size < 2) {
        return vector<double>();
    }

    vector<double> res(size - 1);

    vector<double> llvs = ti_LLV(low, period);
    vector<double> hhvs = ti_HHV(high, period);

    bool isUp = close[1] > close[0];
    double m_af = af;
    double sarTemp = isUp ? high[0] : low[0];
    res[0] = sarTemp;

    for (std::size_t i = 2; i < size; i++) {

        if (sarTemp < close[i - 1] == isUp) {
            if (isUp) {
                if (high[i] > high[i - 1]) {
                    m_af = (m_af + incre > critical) ? af : (m_af + incre);
                }
            }else {
                if (low[i] < low[i - 1]) {
                    m_af = (m_af + incre > critical) ? af : (m_af + incre);
                }
            }
        }else {
            m_af = af;
            isUp = !isUp;
        }

        sarTemp = sarTemp + m_af * ((isUp ? high[i - 1] : low[i - 1]) - sarTemp);

        if (isUp) {
            if (sarTemp >= close[i]) {
                sarTemp = hhvs[i];
            }
        }else {
            if (sarTemp <= close[i]) {
                sarTemp = llvs[i];
            }
        }
        res[i - 1] = sarTemp;
    }
    return res;
}


#endif // TECHNICAL_INDICATORS_H
