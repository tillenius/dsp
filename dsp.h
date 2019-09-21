#pragma once
#pragma once
#include <Windows.h>
#include <gdiplus.h>
#include "graph.h"

class DSP {
public:
	HINSTANCE m_hInstance;
	HWND m_mainWindow;

	void OnCreate(HWND hwnd);
	void OnPaint(HWND hWnd);
	void OnSize(HWND hWnd);

	void Layout(const Gdiplus::RectF & rect);
	void Paint(Gdiplus::Graphics & graphics);
	static LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);

	Gdiplus::Bitmap * m_bitmap = NULL;

	std::vector< std::unique_ptr< Graph > > m_graphs;
};

extern DSP g_app;
