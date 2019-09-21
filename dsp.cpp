#include "dsp.h"

#include <Windowsx.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include <complex>
#include <chrono>
#include <iostream>
#include "fftw3.h"

using namespace Gdiplus;
#pragma comment(lib, "Gdiplus.lib")
#pragma comment(lib, "libfftw3-3.lib")

DSP g_app;

static constexpr double g_sampleRate = 48000.0;


class LinkwitzRileyBPF { // https://github.com/dimtass/DSP-Cpp-filters/blob/master/src/so_linkwitz_riley_lpf.h
public:
	LinkwitzRileyBPF(int fc, int fs = 48000);
	double filter(double sample);

private:
	struct tp_coeffs {
		double a0, a1, a2, b1, b2;
	};
	double m_xnz1 = 0.0, m_xnz2 = 0.0, m_ynz1 = 0.0, m_ynz2 = 0.0;
	tp_coeffs m_coeffs;
};

LinkwitzRileyBPF::LinkwitzRileyBPF(int fc, int fs) {
	double th = M_PI * fc / fs;
	double Wc = M_PI * fc;
	double k = Wc / tan(th);

	double d = pow(k, 2.0) + pow(Wc, 2.0) + 2.0 * k * Wc;
	m_coeffs.a0 = pow(Wc, 2.0) / d;
	m_coeffs.a1 = 2.0 * pow(Wc, 2.0) / d;
	m_coeffs.a2 = m_coeffs.a0;
	m_coeffs.b1 = (-2.0 * pow(k, 2.0) + 2.0 * pow(Wc, 2.0)) / d;
	m_coeffs.b2 = (-2.0 * k * Wc + pow(k, 2.0) + pow(Wc, 2.0)) / d;
}

double LinkwitzRileyBPF::filter(double sample) {
	double xn = sample;
	double yn = m_coeffs.a0 * xn + m_coeffs.a1 * m_xnz1 + m_coeffs.a2 * m_xnz2
		- m_coeffs.b1 * m_ynz1 - m_coeffs.b2 * m_xnz2;

	m_xnz2 = m_xnz1;
	m_xnz1 = xn;
	m_xnz2 = m_ynz1;
	m_ynz1 = yn;

	return(yn);
}

void Zero(std::vector<double> & data) {
	const size_t num = data.size();
	for (int i = 0; i < num; ++i) {
		data[i] = 0.0;
	}
}

void Sine(std::vector<double> & data, double sampleRate, double freq) {
	const size_t num = data.size();
	for (int i = 0; i < num; ++i) {
		data[i] += sin((i / sampleRate) * freq * 2.0 * M_PI);
	}
}

void Hanning(std::vector<double> & data) {
	const size_t num = data.size();
	for (int i = 0; i < num; ++i) {
		double multiplier = 0.5 * (1.0 - cos(2.0 * M_PI * i / (num - 1.0)));
		data[i] = multiplier * data[i];
	}
}

void Normalize(std::vector<double> & data) {
	double ampl = 1.0;
	for (size_t i = 0; i < data.size(); ++i) {
		if (fabs(data[i]) > ampl)
			ampl = fabs(data[i]);
	}
	if (ampl > 1.0) {
		for (size_t i = 0; i < data.size(); ++i) {
			data[i] /= ampl;
		}
	}
}

void Scale(std::vector<double> & data, double s) {
	for (size_t i = 0; i < data.size(); ++i) {
		data[i] *= s;
	}
}

void DSP::OnCreate(HWND hwnd) {
	m_mainWindow = hwnd;
	RealGraph * input = new RealGraph();
	m_graphs.emplace_back(input);
	input->m_colBackground = Color(255, 255, 245, 240);
	ComplexGraph * transformed = new ComplexGraph();
	m_graphs.emplace_back(transformed);
	transformed->m_colBackground = Color(255, 235, 255, 250);
	RealGraph * back = new RealGraph();
	m_graphs.emplace_back(back);
	back->m_colBackground = Color(255, 255, 245, 240);

	RealGraph * expected = new RealGraph();
	m_graphs.emplace_back(expected);
	expected->m_colBackground = Color(255, 255, 245, 240);


	int N = 512;
	input->m_data.resize(N);
	Zero(input->m_data);
	Sine(input->m_data, g_sampleRate, 440.0);
	Sine(input->m_data, g_sampleRate, 1000.0);
	Sine(input->m_data, g_sampleRate, 2000.0);
	Sine(input->m_data, g_sampleRate, 4000.0);
	Normalize(input->m_data);
	Hanning(input->m_data);

	expected->m_data.resize(N);
	Zero(expected->m_data);
	Sine(expected->m_data, g_sampleRate, 440.0);
	Sine(expected->m_data, g_sampleRate, 1000.0);
	Sine(expected->m_data, g_sampleRate, 2000.0);
	Sine(expected->m_data, g_sampleRate, 4000.0);
	Normalize(expected->m_data);
	Hanning(expected->m_data);

	transformed->m_data.resize(N / 2 + 1);
	back->m_data.resize(N);
	fftw_plan plan_forward = fftw_plan_dft_r2c_1d(N, input->m_data.data(), (fftw_complex *) transformed->m_data.data(), FFTW_ESTIMATE);
	fftw_plan plan_backward = fftw_plan_dft_c2r_1d(N, (fftw_complex *) transformed->m_data.data(), back->m_data.data(), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

	auto t0 = std::chrono::high_resolution_clock::now();
	unsigned long long t0r = __rdtsc();
	fftw_execute(plan_forward);
	unsigned long long t1r = __rdtsc();
	auto t1 = std::chrono::high_resolution_clock::now();
	for (int i = 33; i < transformed->m_data.size(); ++i) {
		transformed->m_data[i] = 0;
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	unsigned long long t2r = __rdtsc();
	fftw_execute(plan_backward);
	unsigned long long t3r = __rdtsc();
	auto t3 = std::chrono::high_resolution_clock::now();

	typedef std::chrono::duration<double> fs;
	typedef std::chrono::duration<double, std::micro> fus;
	fus d0 = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
	fus d1 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
	char dbg[200];
	sprintf_s(dbg, "time: %.2f us fft: %.2f us ifft: %.2f us\n", std::chrono::duration_cast<fus>(fs(N/g_sampleRate)).count(), d0.count(), d1.count() );
	OutputDebugStringA(dbg);

	LinkwitzRileyBPF lr_filter(1000);
	auto t4 = std::chrono::high_resolution_clock::now();
	unsigned long long t4r = __rdtsc();
	for (int i = 0; i < expected->m_data.size(); ++i) {
		expected->m_data[i] = lr_filter.filter(expected->m_data[i]);
	}
	unsigned long long t5r = __rdtsc();
	auto t5 = std::chrono::high_resolution_clock::now();
	fus d2 = std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4);

	sprintf_s(dbg, "time: %.2f us LinkwitzRileyBPF: %.2f us\n", std::chrono::duration_cast<fus>(fs(N / g_sampleRate)).count(), d2.count());
	OutputDebugStringA(dbg);

	sprintf_s(dbg, "rdtsc: %lld %lld %lld\n", t1r - t0r, t3r - t2r, t5r - t4r);
	OutputDebugStringA(dbg);

	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_backward);

	Scale(back->m_data, 1.0 / N);

	for (int i = 0; i < 100 && i < N; ++i) {
		double freq = i * g_sampleRate / (double) N;
		double ampl = std::abs(transformed->m_data[i]) / N;
		sprintf_s(dbg, "%d freq %.1f ampl %f\n", i, freq, ampl);
		OutputDebugStringA(dbg);
	}

	//for (int i = 0; i < N/2 + 1; ++i) {
	//	transformed->m_data[i] = std::complex<double>(out[i][0], out[i][1]);
	//}

	//back->m_data.resize(N);
	//for (int i = 0; i < N; ++i) {
	//	back->m_data[i] = inverse[i];
	//}
}

void DSP::OnSize(HWND hWnd) {
	RECT rect;
	GetClientRect(hWnd, &rect);
	Layout(RectF((Gdiplus::REAL) rect.left, (Gdiplus::REAL) rect.top, (Gdiplus::REAL) (rect.right - rect.left), (Gdiplus::REAL) (rect.bottom - rect.top)));
	RedrawWindow(m_mainWindow, NULL, NULL, RDW_INVALIDATE | RDW_UPDATENOW);
}

void DSP::OnPaint(HWND hWnd) {
	RECT rect;
	GetClientRect(hWnd, &rect);

	UINT width = rect.right - rect.left;
	UINT height = rect.bottom - rect.top;
	if (m_bitmap == NULL || m_bitmap->GetWidth() != width || m_bitmap->GetHeight() != height) {
		delete m_bitmap;
		m_bitmap = new Gdiplus::Bitmap(width, height);
	}

	PAINTSTRUCT ps;
	HDC hdc = BeginPaint(hWnd, &ps);
	Graphics graphics(hdc);
	Gdiplus::Graphics * g1 = graphics.FromImage(m_bitmap);
	g1->Clear(Color(255, 255, 255, 255));
	Paint(*g1);
	graphics.DrawImage(m_bitmap, 0, 0);
	EndPaint(hWnd, &ps);
	delete g1;
}

void DSP::Layout(const RectF & rect) {
	const float h = rect.Height / m_graphs.size();
	float y = rect.GetTop();
	for (int i = 0; i < m_graphs.size(); ++i) {
		m_graphs[i]->m_rect = rect;
		m_graphs[i]->m_rect.Y = y;
		m_graphs[i]->m_rect.Height = h;
		y += h;
	}
}

void DSP::Paint(Graphics & graphics) {
	for (int i = 0; i < m_graphs.size(); ++i) {
		m_graphs[i]->Paint(graphics);
	}
}

LRESULT CALLBACK DSP::WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam) {
	switch (message) {
	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	case WM_CLOSE:
		DestroyWindow(hWnd);
		return 0;
	case WM_CREATE:
		g_app.OnCreate(hWnd);
		return 0;
	case WM_SIZE:
		g_app.OnSize(hWnd);
		return 0;
	case WM_ERASEBKGND:
		return 1;
	case WM_PAINT:
		g_app.OnPaint(hWnd);
		return 0;
	//case WM_MOUSEMOVE:
	//	g_app.OnMouseMove(hWnd, wParam, GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
	//	return 0;
	//case WM_LBUTTONDOWN:
	//	g_app.OnLButtonDown(hWnd, wParam, GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
	//	return 0;
	//case WM_LBUTTONUP:
	//	g_app.OnLButtonUp(hWnd, wParam, GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
	//	return 0;
	//case WM_RBUTTONUP:
	//	g_app.OnRButtonUp(hWnd, wParam, GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
	//	return 0;
	//case WM_MOUSEWHEEL:
	//	g_app.OnMouseWheel(hWnd, GET_KEYSTATE_WPARAM(wParam), GET_WHEEL_DELTA_WPARAM(wParam));
	//	return 0;
	//case WM_COMMAND:
	//	if (HIWORD(wParam) == 0) {
	//		g_app.OnMenuCommand(hWnd, LOWORD(wParam));
	//		return 0;
	//	}
	//	return DefWindowProc(hWnd, message, wParam, lParam);
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
}

int WINAPI wWinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ PWSTR pCmdLine, _In_ int nCmdShow) {

	g_app.m_hInstance = hInstance;

	if (OleInitialize(NULL) != S_OK) {
		return 0;
	}

	WNDCLASS wClass{};
	wClass.lpszClassName = L"DSP";
	wClass.hInstance = hInstance;
	wClass.lpfnWndProc = DSP::WndProc;
	wClass.hCursor = LoadCursor(NULL, IDC_ARROW);
	if (RegisterClass(&wClass) == 0) {
		MessageBox(0, L"RegisterClass() failed", 0, MB_OK);
		return 0;
	}

	GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR gdiplusToken;
	if (GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL) != Gdiplus::Ok) {
		MessageBox(0, L"GdiplusStartup() failed", 0, MB_OK);
		return 0;
	}

	HWND hWnd = CreateWindowEx(0, L"DSP", L"DSP", WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, NULL, NULL, hInstance, NULL);
	if (hWnd == NULL) {
		MessageBox(0, L"CreateWindow failed", 0, MB_OK);
		return 0;
	}

	ShowWindow(hWnd, nCmdShow);
	UpdateWindow(hWnd);

	MSG msg;
	while (BOOL bRet = GetMessage(&msg, NULL, 0, 0) != 0) {
		if (bRet == -1) {
			return 0;
		}
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}

	OleUninitialize();

	return (int)msg.wParam;
}


// https://github.com/dimtass/DSP-Cpp-filters
// https://github.com/tap/TapTools/blob/master/Core/DSP/extensions/FilterLib/source/TTHighpassLinkwitzRiley4.cpp
// https://github.com/mbrucher/AudioTK
