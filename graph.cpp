#include "graph.h"
#define _USE_MATH_DEFINES
#include <math.h>

void GetGraph(std::vector<double> buf, size_t width, std::vector<std::pair<double, double> > & out) {
	size_t bs = buf.size();

	size_t old_index = 0;
	for (size_t i = 0; i < width; ++i) {
		size_t idx = i / (double)width * bs + 0.5;

		if (old_index >= buf.size()) {
			old_index = buf.size() - 1;
		}
		double minval = buf[old_index];
		double maxval = buf[old_index];
		for (size_t j = old_index; j < idx; ++j) {
			if (buf[j] > maxval)
				maxval = buf[j];
			else if (buf[j] < minval)
				minval = buf[j];
		}
		old_index = idx;
		out.push_back(std::make_pair(maxval, minval));
	}
}

struct polar_t {
	double r;
	double theta;
	polar_t(double r, double theta) :r(r), theta(theta) {}
};

void GetGraph(std::vector<std::complex<double>> & buf, size_t width, std::vector< std::pair<polar_t, polar_t> > & out) {
	size_t bs = buf.size();

	size_t old_index = 0;
	double ampl = 1.0;
	for (size_t i = 0; i < width; ++i) {
		size_t idx = i / (double)width * bs + 0.5;

		if (old_index >= buf.size()) {
			old_index = buf.size() - 1;
		}

		polar_t minval{ std::abs(buf[old_index]), std::arg(buf[old_index]) };
		polar_t maxval{ std::abs(buf[old_index]), std::arg(buf[old_index]) };
		if (maxval.r > ampl)
			ampl = maxval.r;
		for (size_t j = old_index; j < idx; ++j) {
			polar_t v(std::abs(buf[j]), std::arg(buf[j]));

			if (v.r > maxval.r) {
				maxval.r = v.r;
				if ( maxval.r > ampl )
					ampl = maxval.r;
			}
			else if (v.r < minval.r)
				minval.r = v.r;
			if (v.theta > maxval.theta)
				maxval.theta = v.theta;
			else if (v.theta < minval.theta)
				minval.theta = v.theta;
		}
		old_index = idx;
		out.push_back(std::make_pair(maxval, minval));
	}

	for (size_t i = 0; i < width; ++i) {
		out[i].first.r /= ampl;
		out[i].first.theta /= 2.0 * M_PI;
		out[i].second.r /= ampl;
		out[i].second.theta /= 2.0 * M_PI;
	}
}

void RealGraph::Paint(Gdiplus::Graphics & graphics) {
	graphics.FillRectangle(&Gdiplus::SolidBrush(m_colBackground), m_rect);

	double mid = m_rect.GetTop() + 0.5 * m_rect.Height;

	Gdiplus::Point p0(m_rect.GetLeft(), mid);
	Gdiplus::Point p1(m_rect.GetRight(), mid);

	graphics.DrawLine(&Gdiplus::Pen(Gdiplus::Color(255, 0, 0, 0)), p0, p1);

	if (m_data.empty()) {
		return;
	}
	std::vector<std::pair<double, double>> graph;
	GetGraph(m_data, m_rect.Width, graph);

	double scale = m_rect.Height / 2.0;

	for (int x = 0; x < graph.size(); ++x) {
		graphics.DrawLine(&Gdiplus::Pen(Gdiplus::Color(128, 80, 80, 80)), Gdiplus::Point(x, mid), Gdiplus::Point(x, mid - graph[x].first * scale));
		graphics.DrawLine(&Gdiplus::Pen(Gdiplus::Color(128, 80, 80, 80)), Gdiplus::Point(x, mid), Gdiplus::Point(x, mid - graph[x].second * scale));
	}
}

void ComplexGraph::Paint(Gdiplus::Graphics & graphics) {
	graphics.FillRectangle(&Gdiplus::SolidBrush(m_colBackground), m_rect);

	double mid = m_rect.GetTop() + 0.5 * m_rect.Height;

	Gdiplus::Point p0(m_rect.GetLeft(), mid);
	Gdiplus::Point p1(m_rect.GetRight(), mid);

	graphics.DrawLine(&Gdiplus::Pen(Gdiplus::Color(255, 0, 0, 0)), p0, p1);

	if (m_data.empty()) {
		return;
	}
	std::vector<std::pair<polar_t, polar_t>> graph;
	GetGraph(m_data, m_rect.Width, graph);

	double scale = m_rect.Height / 2.0;

	for (int x = 0; x < graph.size(); ++x) {
		graphics.DrawLine(&Gdiplus::Pen(Gdiplus::Color(128, 0, 80, 0)), Gdiplus::Point(x, mid), Gdiplus::Point(x, mid - graph[x].first.r * scale));
		graphics.DrawLine(&Gdiplus::Pen(Gdiplus::Color(128, 0, 80, 0)), Gdiplus::Point(x, mid), Gdiplus::Point(x, mid - graph[x].second.r * scale));
		graphics.DrawLine(&Gdiplus::Pen(Gdiplus::Color(128, 80, 0, 0)), Gdiplus::Point(x, mid), Gdiplus::Point(x, mid - graph[x].first.theta * scale));
		graphics.DrawLine(&Gdiplus::Pen(Gdiplus::Color(128, 80, 0, 0)), Gdiplus::Point(x, mid), Gdiplus::Point(x, mid - graph[x].second.theta * scale));
	}
}
