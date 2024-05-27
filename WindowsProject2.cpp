#include <windows.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <unordered_map>

#define M_PI 3.14159

struct vec2 {
    float x;
    float y;
    float dx;
    float dy;
    float density;
    float pressure;
};

// Global variables
HBITMAP hBitmap = NULL;
const int bitmapWidth = 800; // Width of the buffer
const int bitmapHeight = 600; // Height of the buffer

const int particleBoxSizeX = 20;
const int particleBoxSizeY = 20;
const int numberOfBoxes = (bitmapHeight / particleBoxSizeY) * (bitmapWidth / particleBoxSizeX);

std::vector<vec2> particles[numberOfBoxes];

void generate_particles() {
    for (int i = 0; i < 10000; i++) {
        int r1 = rand() % 800;
        int r2 = rand() % 600;
        vec2 vec = { static_cast<float>(r1), static_cast<float>(r2) };
        int R1 = r1 / particleBoxSizeX;
        int R2 = r2 / particleBoxSizeY;
        int index = R1 + R2 * (bitmapHeight / particleBoxSizeX);
        particles[index].push_back(vec);
    }
}

float r(vec2 p1, vec2 p2) {
    return pow(p1.x - p2.x, 2)+pow(p1.y-p2.y,2);
}

float W(float r, float h) {
    if (r < h / 2) {
        return ((8 / (M_PI * pow(h, 2))) * (1 - 6 * pow(r / h, 2) + 6 * pow(r / h, 3)));
    }
    else if (r < h) {
        return (8 / (M_PI * pow(h, 2))) * 2 * pow(1 - r / h, 3);
    }
    else {
        return 0;
    }
}

float dW_dr(float r, float h) {
    if (r < h / 2) {
        return (8 / (M_PI * pow(h, 3))) * (-12 * r / h + 18 * pow(r, 2) / pow(h, 2));
    }
    else if (r < h) {
        return (8 / (M_PI * pow(h, 3))) * (-6 * pow(1 - r / h, 2) / h);
    }
    else {
        return 0;
    }
}

vec2 gradW(vec2 p1, vec2 p2, float h) {
    vec2 r = { p1.x - p2.x, p1.y - p2.y };
    float r_magnitude = sqrt(r.x * r.x + r.y * r.y);
    if (r_magnitude == 0) return { 0, 0 }; // Avoid division by zero
    float dWdr = dW_dr(r_magnitude, h);
    return { dWdr * (r.x / r_magnitude), dWdr * (r.y / r_magnitude) };
}

void update_density() {
    
    for (int j = 0; j < numberOfBoxes; j++) {
        for (size_t i = 0; i < particles[j].size(); ++i) {
            vec2& p = particles[j][i];
            p.density = 0;
            for (int offsetY = -1; offsetY <= 1; ++offsetY) {
                for (int offsetX = -1; offsetX <= 1; ++offsetX) {
                    // Calculate neighboring box indices
                    int neighborX = (j % (bitmapWidth / particleBoxSizeX)) + offsetX;
                    int neighborY = (j / (bitmapWidth / particleBoxSizeX)) + offsetY;

                    // Check if the neighboring box indices are within bounds
                    if (neighborX >= 0 && neighborX < (bitmapWidth / particleBoxSizeX) &&
                        neighborY >= 0 && neighborY < (bitmapHeight / particleBoxSizeY)) {
                        int neighborIndex = neighborX + neighborY * (bitmapWidth / particleBoxSizeX);

                        // Iterate over particles in the neighboring box
                        for (size_t k = 0; k < particles[neighborIndex].size(); ++k) {
                            // Reference to neighboring particle
                            const vec2& p2 = particles[neighborIndex][k];
                            //std::cout << W(1, 5) << std::endl;
                            p.density += W(r(p, p2), 20);
                        }
                    }
                }
            }
            //std::cout << p.density << std::endl;
        }
    }
}

void update_pressure() {
    for (int j = 0; j < numberOfBoxes; j++) {
        for (size_t i = 0; i < particles[j].size(); ++i) {
            vec2& p_i = particles[j][i];
            vec2 force = { 0, 0 };
            for (int offsetY = -1; offsetY <= 1; ++offsetY) {
                for (int offsetX = -1; offsetX <= 1; ++offsetX) {
                    int neighborX = (j % (bitmapWidth / particleBoxSizeX)) + offsetX;
                    int neighborY = (j / (bitmapWidth / particleBoxSizeX)) + offsetY;

                    if (neighborX >= 0 && neighborX < (bitmapWidth / particleBoxSizeX) &&
                        neighborY >= 0 && neighborY < (bitmapHeight / particleBoxSizeY)) {
                        int neighborIndex = neighborX + neighborY * (bitmapWidth / particleBoxSizeX);

                        for (size_t k = 0; k < particles[neighborIndex].size(); ++k) {
                            const vec2& p_j = particles[neighborIndex][k];
                            float pressure_i = p_i.density-1;
                            float pressure_j = p_j.density-1;
                            float density_j = p_j.density;

                            vec2 gradW_ij = gradW(p_i, p_j, 20);
                            force.x += (pressure_i + pressure_j) / (2.0f * density_j) * gradW_ij.x;
                            force.y += (pressure_i + pressure_j) / (2.0f * density_j) * gradW_ij.y;
                        }
                    }
                }
            }
            p_i.dx = force.x;
            p_i.dy = force.y;
            //std::cout << p_i.dx << " " << p_i.dy << std::endl;
        }
    }
}

void update_particles() {
    update_density();
    update_pressure();
    // Create a temporary vector to store updated particles
    //std::vector<std::vector<vec2>> updatedParticles(numberOfBoxes);

    // Iterate over particles in each box

    std::vector<std::vector<vec2>> updatedParticles(numberOfBoxes);

    for (int j = 0; j < numberOfBoxes; j++) {
        for (size_t i = 0; i < particles[j].size(); ++i) {
            vec2& p1 = particles[j][i];
            p1.x += p1.dx * 1;
            p1.y += p1.dy * 1;
            
            if (p1.x < 0) p1.x = 0;
            if (p1.x >= bitmapWidth) p1.x = bitmapWidth - 1;
            if (p1.y < 0) p1.y = 0;
            if (p1.y >= bitmapHeight) p1.y = bitmapHeight - 1;

            int R1 = p1.x / particleBoxSizeX;
            int R2 = p1.y / particleBoxSizeY;
            int index = R1 + R2 * (bitmapWidth / particleBoxSizeX);
            //std::cout << R2-p1.y/20 << std::endl;
            updatedParticles[index].push_back(p1);
        }
    }
    for (int j = 0; j < numberOfBoxes; j++) {
        // Copy the vector of particles from updatedParticles to particles
        particles[j] = updatedParticles[j];
    }

    for (int j = 0; j < numberOfBoxes; j++) {
        for (size_t i = 0; i < particles[j].size(); ++i) {
            vec2 p = particles[j][i];
            if (abs(p.x / 20 - j % 40) >= 1 || abs(p.y / 20-j/40) >=1) {
                std::cout << ":( "<<j << std::endl;
            }
        }
    }
}


// Function to draw onto the buffer
void draw_buffer(HDC hdcBuffer) {
    // Create a BITMAPINFO structure to describe the bitmap format
    BITMAPINFO bmi = {};
    bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bmi.bmiHeader.biWidth = bitmapWidth;
    bmi.bmiHeader.biHeight = -bitmapHeight; // Negative height to indicate a top-down bitmap
    bmi.bmiHeader.biPlanes = 1;
    bmi.bmiHeader.biBitCount = 24; // 24-bit color depth (RGB)
    bmi.bmiHeader.biCompression = BI_RGB;

    // Create a buffer to hold pixel data (24 bits per pixel)
    BYTE* buffer = new BYTE[bitmapWidth * bitmapHeight * 3]; // 3 bytes per pixel (RGB)

    // Clear the buffer to white
    memset(buffer, 255, bitmapWidth * bitmapHeight * 3); // White color (RGB: 255, 255, 255)

    // Draw particles onto the buffer
    for (int i = 0; i < numberOfBoxes; i++) {
        std::vector<vec2>& boxParticles = particles[i];
        if (boxParticles.size() != 0){
            //std::cout << boxParticles[0].density << std::endl;
        }
        for (const vec2& p : boxParticles) {
            // Calculate the offset of the current pixel in the buffer
            int x = static_cast<int>(p.x);
            int y = static_cast<int>(p.y);
            //std::cout << p.density<<std::endl;

            // Check if particle is within buffer bounds
            if (x >= 0 && x < bitmapWidth && y >= 0 && y < bitmapHeight) {
                int offset = (y * bitmapWidth + x) * 3; // 3 bytes per pixel

                // Set pixel color to black at particle coordinates
                buffer[offset] = 0;  // Blue component
                buffer[offset + 1] = 0; // Green component
                vec2 R = { p.dx * pow(10,5), p.dy * pow(10,5) };
                vec2 O = { 0,0 };
                //std::cout << r(R, O) << std::endl;
                buffer[offset + 2] = r(R, O); ; // Red component
            }
        }
    }

    // Set the pixel data in the buffer to the device context, covering the entire client area
    SetDIBitsToDevice(hdcBuffer, 0, 0, bitmapWidth, bitmapHeight, 0, 0, 0, bitmapHeight, buffer, &bmi, DIB_RGB_COLORS);

    // Clean up
    delete[] buffer;
}

// Window procedure
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    switch (uMsg) {
    case WM_CREATE:
    {
        // Create a compatible bitmap for drawing
        HDC hdc = GetDC(hwnd);
        hBitmap = CreateCompatibleBitmap(hdc, bitmapWidth, bitmapHeight);
        ReleaseDC(hwnd, hdc);
    }
    break;
    case WM_PAINT:
    {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);

        // Draw the buffer onto the window
        HDC hdcBuffer = CreateCompatibleDC(hdc);
        SelectObject(hdcBuffer, hBitmap);
        draw_buffer(hdcBuffer); // Draw onto the buffer
        BitBlt(hdc, 0, 0, bitmapWidth, bitmapHeight, hdcBuffer, 0, 0, SRCCOPY); // Transfer buffer to window
        DeleteDC(hdcBuffer);

        EndPaint(hwnd, &ps);
    }
    break;
    case WM_DESTROY:
        if (hBitmap != NULL)
            DeleteObject(hBitmap);
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hwnd, uMsg, wParam, lParam);
    }
    return 0;
}

// WinMain function
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PSTR lpCmdLine, int nCmdShow) {
    AllocConsole();
    // Redirect standard input and output to the console window
    FILE* pConsoleIn, * pConsoleOut;
    freopen_s(&pConsoleIn, "CONIN$", "r", stdin);
    freopen_s(&pConsoleOut, "CONOUT$", "w", stdout);
    freopen_s(&pConsoleOut, "CONOUT$", "w", stderr);
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    generate_particles();

    // Register window class
    const wchar_t CLASS_NAME[] = L"SampleWindowClass";

    WNDCLASS wc = {};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;

    RegisterClass(&wc);

    // Create window
    HWND hwnd = CreateWindowEx(
        0,
        CLASS_NAME,
        L"My Window",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 800, 600,
        NULL, NULL, hInstance, NULL
    );

    if (hwnd == NULL)
        return 0;

    ShowWindow(hwnd, nCmdShow);

    // Message loop
    MSG msg = {};
    while (true) {
        while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
            if (msg.message == WM_QUIT)
                return 0;
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
        std::cout << "hi" << std::endl;
        update_particles();
        InvalidateRect(hwnd, NULL, TRUE); // Invalidate the entire window
        UpdateWindow(hwnd); // Force the window to redraw
    }

    return 0;
}
