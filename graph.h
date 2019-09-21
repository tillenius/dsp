#pragma once
#include <Windows.h>
#include <gdiplus.h>
#include <vector>
#include <complex>

class Graph {
public:
	virtual void Paint(Gdiplus::Graphics & graphics) = 0;

	Gdiplus::RectF m_rect;
	Gdiplus::Color m_colBackground;

};

class RealGraph : public Graph {
public:
	void Paint(Gdiplus::Graphics & graphics) final override;

	std::vector<double> m_data;
};

class ComplexGraph : public Graph {
public:

	void Paint(Gdiplus::Graphics & graphics) final override;

	std::vector<std::complex<double>> m_data;
};
